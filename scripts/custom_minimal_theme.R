#' Custom Minimal Theme with Grid
#'
#' Creates a clean, minimal ggplot2 theme with a white background and subtle grid lines.
#' This theme is designed to provide a professional and readable appearance for data visualizations.
#'
#' @return A ggplot2 theme object that can be added to any ggplot
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   custom_minimal_theme_with_grid()
custom_minimal_theme_with_grid <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 14, face = "plain"),
      axis.title.y = element_text(size = 14, face = "plain"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    #   legend.title = element_blank() # Remove legend title
    )
}
