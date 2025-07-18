library(testthat)
library(omixVizR)
library(ggplot2)

test_that("plot_gwas_power_bt executes without errors and returns correct structure", {
  # Run the function without saving the plot
  suppressMessages({power_results <- plot_gwas_power_bt(
        n_cases = 4324,
        n_controls = 93945,
        maf_levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
        or_range = seq(1.01, 2.00, 0.01),
        save_plot = TRUE
    )
  })
  
  # 1. Check the returned object structure
  expect_type(power_results, "list")
  expect_named(power_results, c("plot", "power_data"))
  
  # 2. Check the plot object
  expect_s3_class(power_results$plot, "ggplot")
  
  # 3. Check the data object
  expect_s3_class(power_results$power_data, "data.table")
})
