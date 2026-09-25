"""Import every module of grace_pipeline and report failures.

Run from code/:  python check_imports.py
Exit status 1 if any module fails to import.
"""
import importlib
import pkgutil
import sys

sys.path.insert(0, ".")
import grace_pipeline  # noqa: E402

failed = []
mods = list(pkgutil.iter_modules(grace_pipeline.__path__))
for m in mods:
    name = f"grace_pipeline.{m.name}"
    try:
        importlib.import_module(name)
    except Exception as e:  # noqa: BLE001
        failed.append((name, f"{type(e).__name__}: {e}"))
print(f"{len(mods) - len(failed)}/{len(mods)} modules import")
for name, err in failed:
    print(f"FAIL {name}: {err}")
sys.exit(1 if failed else 0)
