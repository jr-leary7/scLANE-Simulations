theme_analysis <- function(base.size = 12, 
                           base.lwd = 0.75, 
                           base.font = "sans") {
  ggplot2::theme_classic(base_size = base.size, 
                         base_family = base.font, 
                         base_line_size = base.lwd, 
                         base_rect_size = base.lwd) + 
  ggplot2::theme(strip.clip = "off", 
                 strip.background = ggplot2::element_rect(linewidth = base.lwd))
}
