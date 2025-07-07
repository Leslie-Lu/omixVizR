test_that("manhattan returns ggplot", {
  skip_if_not_installed("qqman")
  data <- qqman::gwasResults
  data.table::fwrite(
    data,
    file = file.path(tempdir(), "gwasResults.txt"),
    sep = "\t",
    row.names = FALSE
  )
  txt_path = file.path(tempdir(), "gwasResults.txt")
  plot_qqman(
    txt_path,
    pheno_name = "Test",
    output_graphics = "png"
  )
  out_file1 <- file.path(tempdir(), "Test_Manhattan_plot.png")
  out_file2 <- file.path(tempdir(), "Test_QQ_plot.png")
  expect_true(file.exists(out_file1))
  expect_true(file.exists(out_file2))
  unlink(out_file1)
  unlink(out_file2)
})
