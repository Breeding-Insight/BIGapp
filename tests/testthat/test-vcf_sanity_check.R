# Test for vcf_sanity_check function

test_that("vcf_sanity_check works on a valid VCF", {
  # Use a small example VCF file from inst/ or create a minimal one
  vcf_path <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
  
  res <- vcf_sanity_check(vcf_path, n_data_lines = 100, max_markers = 10000, verbose = FALSE)
  expect_s3_class(res, "vcf_sanity_check")
  expect_true(res$checks["VCF_header"])
  expect_true(res$checks["VCF_columns"])
  expect_true(res$checks["max_markers"])
  expect_true(res$checks["samples"])
  expect_true(res$checks["chrom_info"])
  expect_true(res$checks["pos_info"])
  expect_false(res$checks["ref_alt"])
  
  expect_error(vcf_sanity_check("nonexistent.vcf"), "File does not exist")
})


