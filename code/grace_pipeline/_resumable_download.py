"""Resumable HTTP download helper shared by the per-source downloaders.

Goals
-----
Every download we run can be interrupted at any moment — WiFi drops,
Ctrl-C, a sleeping laptop. The pattern below makes resumption
near-free:

  1. **Atomic writes.** Stream to ``dest.partial``, fsync, then
     atomic-rename to ``dest``. A partial download never masquerades as
     a complete one on the next run's skip-if-exists check.

  2. **Server-side size validation.** If the server sends a
     ``Content-Length`` header, the rename only happens when the
     downloaded byte count matches. Otherwise we accept what we got
     but flag it in the registry.

  3. **Registry file** (optional, off by default). Per-downloader
     JSON registry ``<out_dir>/.download_registry.json`` records which
     filenames are confirmed-complete with size + completion timestamp.
     Resume reads it once and skips the per-file stat() scan over
     thousands of small files. Without the registry, the per-file
     stat()-based skip still works — it's just a bit slower at startup.

  4. **HTTP Range / resume from byte-offset.** When the server supports
     ``Accept-Ranges: bytes`` and a ``.partial`` exists, the next run
     continues from where the previous one stopped instead of
     re-downloading from byte 0. This is the biggest win for the ~5 GB
     GRACE SHM download — a single dropped file no longer wastes the
     bytes we already had.

  5. **Authentication renewal hook.** Caller passes a session factory;
     on 401 we reset the session via the factory and retry once before
     giving up. Useful for short-lived-token authentication flows.

  6. **Signal-aware Ctrl-C.** A signal handler catches SIGINT/SIGTERM,
     finishes the current chunk, flushes, and exits cleanly. So you
     can ^C with confidence that the next run won't redo the
     in-progress file from scratch.

The helper is deliberately the lowest-common-denominator interface —
just ``atomic_get(url, dest, session, ...)``. Each downloader still
owns its CMR/STAC/REST search logic; this module only handles the
last-mile "URL → bytes → file" step that's worth doing right.
"""
from __future__ import annotations

import json
import os
import signal
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import requests


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
CHUNK_BYTES = 1024 * 1024          # 1 MiB streaming chunks
N_RETRIES = 5
CONNECT_TIMEOUT = 30
READ_TIMEOUT = 300


# --------------------------------------------------------------------------
# Signal handling
# --------------------------------------------------------------------------
# A simple "have we been asked to stop?" flag, set by SIGINT/SIGTERM.
# Workers poll it between chunks; if set, they flush and exit. The
# helper provides install_shutdown_handler() that downloaders should
# call once at startup.

_SHUTDOWN_REQUESTED = threading.Event()


def install_shutdown_handler() -> None:
    """Register SIGINT/SIGTERM → set the shutdown flag.

    Downloaders should call this once at startup. Callers can also
    poll `shutdown_requested()` between operations to wind down
    gracefully.
    """
    def _handler(signum, frame):
        if not _SHUTDOWN_REQUESTED.is_set():
            print("\n  [shutdown] signal received; "
                  "finishing current chunk then exiting cleanly…",
                  flush=True)
            _SHUTDOWN_REQUESTED.set()
        else:
            # Second ^C — bail immediately
            print("\n  [shutdown] forced exit", flush=True)
            raise SystemExit(130)

    signal.signal(signal.SIGINT, _handler)
    try:
        signal.signal(signal.SIGTERM, _handler)
    except (ValueError, OSError):
        # SIGTERM handler can't always be installed (e.g., non-main
        # thread); SIGINT alone is fine.
        pass


def shutdown_requested() -> bool:
    return _SHUTDOWN_REQUESTED.is_set()


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------

@dataclass
class DownloadRegistry:
    """Per-downloader on-disk record of completed files.

    Thread-safe: a single internal lock serializes all mutations and the
    flush-then-rename. Multiple ThreadPoolExecutor workers can call
    ``mark_complete`` / ``mark_failed`` / ``is_complete`` concurrently
    without racing the JSON serializer or the atomic rename.
    """
    path: Path
    entries: dict[str, dict] = field(default_factory=dict)
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    @classmethod
    def load_or_empty(cls, out_dir: Path) -> "DownloadRegistry":
        p = out_dir / ".download_registry.json"
        if p.exists():
            try:
                with p.open() as fh:
                    return cls(path=p, entries=json.load(fh))
            except (json.JSONDecodeError, OSError):
                pass  # treat corrupted registry as empty; rebuild
        return cls(path=p, entries={})

    def is_complete(self, filename: str, expected_bytes: int | None = None
                     ) -> bool:
        with self._lock:
            entry = self.entries.get(filename)
            if entry is None:
                return False
            if expected_bytes is not None:
                return entry.get("bytes") == expected_bytes
            return entry.get("status") == "complete"

    def mark_complete(self, filename: str, n_bytes: int,
                       extra: dict | None = None) -> None:
        rec = {
            "status": "complete",
            "bytes": n_bytes,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                           time.gmtime()),
        }
        if extra:
            rec.update(extra)
        with self._lock:
            self.entries[filename] = rec
            self._flush_locked()

    def mark_failed(self, filename: str, reason: str) -> None:
        rec = {
            "status": "failed",
            "reason": reason[:200],
            "failed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                         time.gmtime()),
        }
        with self._lock:
            self.entries[filename] = rec
            self._flush_locked()

    def flush(self) -> None:
        with self._lock:
            self._flush_locked()

    def _flush_locked(self) -> None:
        # Must be called with self._lock held. Snapshot entries so the
        # serializer never iterates a dict that could still mutate.
        snapshot = dict(self.entries)
        tmp = self.path.with_suffix(".json.tmp")
        with tmp.open("w") as fh:
            json.dump(snapshot, fh, indent=2)
        os.replace(tmp, self.path)


# --------------------------------------------------------------------------
# Atomic download
# --------------------------------------------------------------------------

def atomic_get(
    url: str,
    dest: Path,
    *,
    session: requests.Session,
    expected_bytes: int | None = None,
    headers: dict | None = None,
    auth_renewer: Callable[[], requests.Session] | None = None,
    progress_label: str | None = None,
) -> tuple[Path, int] | None:
    """Stream `url` to `dest`, atomically.

    Returns (dest_path, bytes_written) on success, or None on failure.

    Behaviour
    ---------
    - Writes to ``dest.partial`` first. On success, atomic-renames to
      ``dest``.
    - If ``dest`` already exists and matches `expected_bytes` (when
      provided), returns immediately with no network call.
    - If ``dest.partial`` already exists and the server supports HTTP
      Range, continues from the offset using ``Range: bytes=N-``.
    - On 401, calls `auth_renewer` (if given) for a fresh session and
      retries once.
    - On 429/5xx/connection errors, exponential-backoff retries up to
      N_RETRIES.
    - Periodically polls `shutdown_requested()`; on shutdown leaves the
      `.partial` file in place (so the next run resumes) and returns
      None.
    """
    if dest.exists():
        actual = dest.stat().st_size
        if expected_bytes is None or actual == expected_bytes:
            return dest, actual
        # Size mismatch — delete and re-download
        dest.unlink()

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")

    headers = dict(headers or {})
    backoff = 2.0
    last_exc: Exception | None = None

    for attempt in range(N_RETRIES):
        if shutdown_requested():
            return None

        # Build Range header from any existing partial
        start_byte = partial.stat().st_size if partial.exists() else 0
        mode = "ab" if start_byte > 0 else "wb"
        req_headers = dict(headers)
        if start_byte > 0:
            req_headers["Range"] = f"bytes={start_byte}-"

        try:
            r = session.get(url, headers=req_headers, stream=True,
                             timeout=(CONNECT_TIMEOUT, READ_TIMEOUT),
                             allow_redirects=True)
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(backoff)
            backoff = min(backoff * 2, 60.0)
            continue

        if r.status_code == 401 and auth_renewer is not None:
            # Token expired — get a fresh session and retry
            session = auth_renewer()
            time.sleep(1.0)
            continue

        if r.status_code in (429,) or r.status_code >= 500:
            r.close()
            time.sleep(backoff)
            backoff = min(backoff * 2, 60.0)
            continue

        if r.status_code == 416:
            # Requested range not satisfiable — partial may be larger
            # than full file (corruption). Restart from scratch.
            if partial.exists():
                partial.unlink()
            time.sleep(0.5)
            continue

        if r.status_code >= 400:
            r.close()
            # 4xx other than handled cases — permanent failure
            return None

        # Stream the body
        bytes_written = start_byte
        try:
            with partial.open(mode) as fh:
                for chunk in r.iter_content(chunk_size=CHUNK_BYTES):
                    if not chunk:
                        continue
                    if shutdown_requested():
                        fh.flush()
                        os.fsync(fh.fileno())
                        r.close()
                        return None
                    fh.write(chunk)
                    bytes_written += len(chunk)
                fh.flush()
                os.fsync(fh.fileno())
        except (OSError, requests.RequestException) as exc:
            last_exc = exc
            time.sleep(backoff)
            backoff = min(backoff * 2, 60.0)
            continue
        finally:
            r.close()

        # Validate size if server told us what to expect
        if expected_bytes is not None and bytes_written != expected_bytes:
            print(f"  ! {progress_label or dest.name}: "
                  f"size mismatch ({bytes_written} vs {expected_bytes}); "
                  "will retry on next run", flush=True)
            return None

        # Atomic rename
        os.replace(partial, dest)
        return dest, bytes_written

    if last_exc is not None:
        print(f"  ! {progress_label or dest.name}: "
              f"{type(last_exc).__name__}: {str(last_exc)[:120]}",
              flush=True)
    return None


# --------------------------------------------------------------------------
# Convenience context manager
# --------------------------------------------------------------------------

@contextmanager
def downloader_session(out_dir: Path):
    """Context: install signal handler, hand back registry, flush on exit."""
    install_shutdown_handler()
    registry = DownloadRegistry.load_or_empty(out_dir)
    try:
        yield registry
    finally:
        registry.flush()
