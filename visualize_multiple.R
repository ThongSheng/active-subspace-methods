library(ggplot2)
library(dplyr)

visualize_multiple <- function(data_list, entry_row, entry_col, true_mat = matrix(1/45*c(120, 50, 50, 21), nrow=2)) {
  
  combined_data <- bind_rows(data_list, .id = "dataset_source")
  combined_data$dataset_source <- factor(combined_data$dataset_source, levels = names(data_list))
  
  param <- names(data_list[[1]])[4]
  
  ggplot(subset(combined_data, row == entry_row & column == entry_col),
         aes(x = factor(.data[[param]]), y = value, fill = dataset_source)) +
    geom_violin(scale = "width", trim = FALSE) +
    geom_hline(yintercept = true_mat[entry_row, entry_col],
               col = "red", linetype = "dashed", linewidth = 1) +
    labs(x = param,
         y = "Entry value",
         title = sprintf("Distribution of [%s,%s] Matrix Entry Across Iterations", entry_row, entry_col),
         subtitle = sprintf("True value: %.2f (red dashed line)", true_mat[entry_row, entry_col]),
         fill = "Dataset") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_viridis_d()
}
