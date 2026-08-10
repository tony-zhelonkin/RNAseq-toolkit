test_that("gs_db carries the frozen contract", {
  db <- bulkiRNA:::gs_db(
    list(A = c("Actb", "Actb", "Gapdh"), B = c("Sdha")),
    database = "testdb", species = "Mus_musculus",
    pathway_names = c(A = "Set A")
  )
  expect_s3_class(db, "gs_db")
  expect_true(is.list(db))
  expect_identical(names(db), c("A", "B"))
  expect_identical(db$A, c("Actb", "Gapdh"))  # de-duplicated
  expect_identical(attr(db, "database"), "testdb")
  expect_identical(attr(db, "species"), "Mus musculus")  # normalised
  expect_identical(attr(db, "gene_id_type"), "symbol")
  expect_identical(attr(db, "pathway_names"),
                   c(A = "Set A", B = "B"))
  expect_identical(names(attr(db, "pathway_names")), names(db))
})

test_that("gs_db drops empty sets and rejects bad input", {
  db <- bulkiRNA:::gs_db(list(A = "Actb", B = character(), C = c(NA, "")),
                         database = "d", species = "Mus musculus")
  expect_identical(names(db), "A")

  expect_error(bulkiRNA:::gs_db(list(1, 2), database = "d",
                                species = "Mus musculus"),
               "named list")
  expect_error(bulkiRNA:::gs_db(list(A = "x", A = "y"), database = "d",
                                species = "Mus musculus"),
               "unique")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "",
                                species = "Mus musculus"),
               "non-empty string")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "d", species = NA),
               "`species` must be")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "d",
                                species = "Mus musculus",
                                gene_id_type = "entrez"),
               "gene_id_type")
})

test_that("subsetting and summary keep the contract", {
  db <- bulkiRNA:::gs_db(list(A = "a", B = c("b", "c")),
                         database = "d", species = "Mus musculus")
  sub <- db["B"]
  expect_s3_class(sub, "gs_db")
  expect_identical(names(sub), "B")
  expect_identical(names(attr(sub, "pathway_names")), "B")
  expect_identical(attr(sub, "database"), "d")

  s <- summary(db)
  expect_identical(names(s), c("pathway_id", "pathway_name", "n_genes"))
  expect_identical(s$n_genes, c(1L, 2L))
  expect_output(print(db), "<gs_db>")
})

test_that("filter_by_size honours open bounds", {
  db <- bulkiRNA:::gs_db(
    list(A = "a", B = c("b", "c"), C = c("d", "e", "f")),
    database = "d", species = "Mus musculus"
  )
  expect_identical(names(bulkiRNA:::filter_by_size(db)), c("A", "B", "C"))
  expect_identical(names(bulkiRNA:::filter_by_size(db, min_size = 2)),
                   c("B", "C"))
  expect_identical(names(bulkiRNA:::filter_by_size(db, max_size = 2)),
                   c("A", "B"))
  expect_identical(names(bulkiRNA:::filter_by_size(db, 2, 2)), "B")
  expect_error(bulkiRNA:::filter_by_size(db, 5, 2), "must not exceed")
  expect_error(bulkiRNA:::filter_by_size(list(A = "a"), 1), "must be a `gs_db`")
})

test_that(".gsdb_as_t2g round-trips the legacy shape", {
  db <- bulkiRNA:::gs_db(list(A = c("a", "b"), B = "c"),
                         database = "d", species = "Mus musculus",
                         pathway_names = c(A = "Set A", B = "Set B"))
  legacy <- bulkiRNA:::.gsdb_as_t2g(db)
  expect_identical(names(legacy), c("T2G", "T2N"))
  expect_identical(names(legacy$T2G), c("gs_name", "gene_symbol"))
  expect_identical(names(legacy$T2N), c("gs_name", "description"))
  expect_identical(nrow(legacy$T2G), 3L)
  expect_identical(legacy$T2N$description, c("Set A", "Set B"))

  back <- bulkiRNA:::.gsdb_from_t2g(legacy, database = "d",
                                    species = "Mus musculus")
  expect_identical(unclass(back)[names(db)], unclass(db)[names(db)])
  expect_identical(attr(back, "pathway_names"), attr(db, "pathway_names"))
})

test_that(".gsdb_from_t2g validates its input", {
  expect_error(bulkiRNA:::.gsdb_from_t2g(list(), "d", "Mus musculus"), "T2G")
  expect_error(
    bulkiRNA:::.gsdb_from_t2g(list(T2G = data.frame(a = 1)), "d",
                              "Mus musculus"),
    "gene_symbol"
  )
})
