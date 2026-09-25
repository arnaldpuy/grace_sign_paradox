theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.key.width = unit(0.4, "cm"),
          legend.key.height = unit(0.5, "lines"),
          legend.key.spacing.y = unit(0, "lines"),
          legend.box.spacing = unit(0, "pt"),
          legend.spacing.y  = unit(0.1, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7.3),
          axis.title.y = element_text(size = 7.3),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 8),
          strip.text.x = element_text(size = 7.4),
          strip.text.y = element_text(size = 7.4),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(margin = margin(t = 1.5, b = 1.5)))
}
