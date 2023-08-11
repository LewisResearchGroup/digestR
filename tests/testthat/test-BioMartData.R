library(testthat)
library(withr)

test_that("test in a temp directory", {
  # Create a unique temp directory
  tmp_dir <- tempfile(pattern = "test_dir")
  dir.create(tmp_dir)

  # Run code within the temp directory and auto-cleanup afterwards
  with_dir(tmp_dir, {
    
    biomart <- BioMartData$new(biomart = "ensembl", dataset = "btaurus_gene_ensembl")

    # Retrieve and process the data
    biomart$get_data(filepath = tmp_dir, chromosomes = c("1",))

    # check file exists
    expect_true(file.exists(file.path(tmp_dir, "btaurus_gene_ensembl.txt")))

    # After the test code, cleanup the directory
    on.exit(recursive_delete(tmp_dir), add = TRUE)
  })
  
  # Here, you're outside the temp directory and can add further tests or checks
  # expect_false(dir.exists(tmp_dir))  # for example, confirming the directory is deleted
})