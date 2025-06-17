library(ggplot2)

visualize <- function(data, entry_row, entry_col, true_mat = matrix(1/45*c(120, 50, 50, 21), nrow=2)) {
  param <- names(data)[4]
  ggplot(subset(data, row == entry_row & column == entry_col), 
         aes(x = factor(.data[[param]]), y = value, fill = factor(.data[[param]]))) +
    geom_violin(scale = "width", trim = FALSE) +
    geom_hline(yintercept = true_mat[entry_row, entry_col],
               col = "red", linetype = "dashed", linewidth = 1) +
    labs(x = param, 
         y = "Entry value", 
         title = sprintf("Distribution of [%s,%s] Matrix Entry Across Iterations ([%s])", entry_row, entry_col, deparse(substitute(data))),
         subtitle = sprintf("True value: %.2f (red dashed line)", true_mat[entry_row, entry_col])) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_viridis_d()
}
