#' Custom Minimal Theme with Grid
#'
#' Creates a clean, minimal ggplot2 theme with a white background and subtle grid lines.
#' This theme is designed to provide a professional and readable appearance for data visualizations.
#'
#' @param base_size Base font size for the theme (default: 12)
#' @param base_family Base font family for the theme (default: "")
#'
#' @return A ggplot2 theme object that can be added to any ggplot
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   custom_minimal_theme_with_grid()
custom_minimal_theme_with_grid <- function(base_size = 12, base_family = "") {
  # Ensure ggplot2 is loaded
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Background elements
      panel.background = ggplot2::element_rect(fill = "white", color = NA), 
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Grid lines
      panel.grid.major = ggplot2::element_line(color = "grey90", size = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Axis elements
      axis.line = ggplot2::element_line(color = "black", size = 0.3),
      axis.ticks = ggplot2::element_line(color = "black", size = 0.3),
      
      # Text elements
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      
      # Plot title
      plot.title = ggplot2::element_text(
        hjust = 0.5, 
        margin = ggplot2::margin(b = 10)
      ),
      
      # Overall plot margins
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}
