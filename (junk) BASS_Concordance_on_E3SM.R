library(ncdf4) # For working with nc files
library(BASS) # For fitting BASS model
library(concordance) # For computing C matrix and concordance
library(ggplot2) # For plotting
library(tictoc) # For timing
library(reshape2) # For reshaping data to plot heatmap

# Import data file
df <- nc_open("AutoCalibration-main/data/lat_lon_10yr_24x48_DJF.nc")

# Extract variables
for (var in names(df$var)) {
  assign(var, ncvar_get(df, var))
}

# Take log of zmconv_tau to stablize its variance
lhs[4,] <- log(lhs[4,])

# Assign X variables
X <- t(lhs)

# Min-max scale
min_max_scale <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
X_scaled <- apply(X, 2, min_max_scale)

# Assign Y variable (put in quotation to allow dynamic label in biplot)
selected_y <- c("LWCF", "SWCF", "PRECT", "PRECL", "PRECC", "PSL", "U200", "U850", "Z500", "TREFHT", "T500", "RH500", "FLNT", "FSNT")

# Initiate lists
mod <- C <- C_eigen <- Z <- list()

# Perform PCA for each Y variable and create biplot visualization
for (i in seq_along(selected_y)) {
  
  # Compute area weighted mean and scale it
  Y <- apply(get(selected_y[[i]]) * area, 3, sum) / apply(area, 3, sum) 
  Y_scaled <- min_max_scale(Y)
  
  # Fit BASS model and record time
  set.seed(42)
  tic()
  if (selected_y[[i]] == "U850") {
    mod[[i]] <- bass(X_scaled, Y_scaled, h2 = 0.1) # define small value of h2 to favor more basis functions
  } else {
    mod[[i]] <- bass(X_scaled, Y_scaled)
  }
  toc()
  
  #plot(mod[[i]])
  
  # Estimate C matrix and perform eigenvalue decomposition
  C[[i]] <- C_bass(mod[[i]])
  C_eigen[[i]] <- eigen(C[[i]])
  
  # Transform the input variables
  Z[[i]] <- (X_scaled - 0.5) %*% C_eigen[[i]]$vectors
  
  # Draw the biplot
  # Create a data frame for points
  df_points <- data.frame(
    Z1 = Z[[i]][, 1],
    Z2 = Z[[i]][, 2],
    Y = Y
  )
  
  # Create a data frame for arrows (loadings)
  df_arrows <- data.frame(
    PC1 = C_eigen[[i]]$vectors[, 1],
    PC2 = C_eigen[[i]]$vectors[, 2],
    variable = ncvar_get(df, "x")
  )
  
  # Plot the biplot
  print(ggplot() +
          geom_point(data = df_points, aes(x = Z1, y = Z2, color = Y)) +
          scale_color_gradient(low = "blue", high = "red") +
          geom_segment(data = df_arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                       arrow = arrow(length = unit(0.2, "cm")), color = "grey30") +
          geom_text(data = df_arrows, aes(x = PC1, y = PC2, label = variable), hjust = 1.2, color = "grey30") +
          labs(x = "First Principal Component", y = "Second Principal Component", title = paste("Biplot of", selected_y[[i]]), color = selected_y[[i]]) +
          theme_minimal())
}

# Initiate concordance matrix
concordance_results <- matrix(NA, length(selected_y), length(selected_y))

# Compute concordance
for (i in 1:(length(selected_y)-1)) {
  for (j in (i+1):length(selected_y)) {
    C_ij <- Cfg_bass(mod[[i]], mod[[j]])
    concordance <- sum(diag(C_ij)) / sqrt(sum(diag(C[[i]])) * sum(diag(C[[j]])))
    concordance_results[i, j] <- concordance
    concordance_results[j, i] <- concordance
  }
}

# Print results
rownames(concordance_results) <- selected_y; colnames(concordance_results) <- selected_y
print(concordance_results)

# Convert data to long format for plotting
concordance_df <- melt(concordance_results, varnames = c("var1", "var2"), value.name = "concordance")
concordance_df["abs_conc"] <- abs(concordance_df$concordance)

# Plot the standard concordance heatmap
ggplot(concordance_df, aes(x = var1, y = var2, fill = concordance)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, na.value = "gray80") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Concordance Heatmap", fill = "Concordance")

# Plot the absolute concordance heatmap
ggplot(concordance_df, aes(x = var1, y = var2, fill = abs_conc)) +
  geom_tile() +
  scale_fill_gradient2(high = "red", na.value = "gray80") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Absolute Concordance Heatmap", fill = "Abs Concordance")

# Compute inverse of C matrix
C_inv <- list()
for (i in seq_along(selected_y)){
  C_inv[[i]] <- solve(C[[i]])
}
print(C_inv)

# Compute inverse of concordance matrix M
M <- concordance_results
M[is.na(M)] <- 0
M_inv <- solve(M)
print(M_inv)

# Close data file
nc_close(df)
