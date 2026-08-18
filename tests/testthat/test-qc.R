test_that("link QC reports coverage and distance", {
  links <- data.frame(id_st = c("s1", "s2"), dist = c(0.1, 0.3))
  unmatched <- data.frame(id_st = "s3")
  qc <- spatialatac:::.link_qc(
    c("s1", "s2", "s3"), c("c1", "c2"), links, unmatched, 50L, 123L
  )
  expect_equal(qc$n_links, 2L)
  expect_equal(qc$n_unmatched_spots, 1L)
  expect_equal(qc$link_fraction, 2 / 3)
  expect_equal(qc$median_distance, 0.2)
})

