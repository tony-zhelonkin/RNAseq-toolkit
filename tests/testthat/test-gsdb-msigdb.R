test_that("gsdb_msigdb validates arguments before touching the network", {
  expect_error(gsdb_msigdb(collection = character()), "`collection` must be")
  expect_error(gsdb_msigdb(species = NA), "`species` must be")
  expect_error(gsdb_msigdb(db_species = "rat"), "'arg' should be one of")
})

test_that("gsdb_msigdb returns a gs_db for the hallmarks", {
  skip_if_not_installed("msigdbr")
  db <- tryCatch(
    gsdb_msigdb("Mus musculus", collection = "H"),
    error = function(e) skip(paste("MSigDB unavailable:", conditionMessage(e)))
  )
  expect_s3_class(db, "gs_db")
  expect_identical(attr(db, "database"), "msigdb_H")  # join key
  expect_identical(attr(db, "database_label"), "MSigDB H")
  expect_identical(attr(db, "species"), "Mus musculus")
  expect_identical(attr(db, "gene_id_type"), "symbol")
  expect_equal(length(db), 50)
  expect_true(all(grepl("^HALLMARK_", names(db))))
  expect_identical(names(attr(db, "pathway_names")), names(db))
  expect_true(all(vapply(db, is.character, logical(1))))

  sub <- gsdb_msigdb("Mus musculus", collection = "H", min_size = 200)
  expect_true(all(vapply(sub, length, integer(1)) >= 200))
  expect_lt(length(sub), length(db))
})

test_that("an unknown collection errors with guidance", {
  skip_if_not_installed("msigdbr")
  err <- tryCatch(gsdb_msigdb(collection = "ZZZ"), error = conditionMessage)
  skip_if(grepl("download|connect|internet|cache", err, ignore.case = TRUE),
          "MSigDB unavailable")
  expect_match(err, "msigdbr_collections|no gene sets|collection")
})
