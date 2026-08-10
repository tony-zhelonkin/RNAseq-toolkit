# Argument validation only — nothing here touches the network.

test_that("download_gatom_references keeps its frozen formals", {
  f <- formals(download_gatom_references)
  expect_identical(names(f),
                   c("dest_dir", "species", "networks", "overwrite"))
  expect_identical(f$dest_dir, "00_data/references/gatom")
  expect_identical(f$species, "Mus_musculus")
  expect_identical(eval(f$networks), c("kegg", "combined"))
  expect_false(eval(f$overwrite))
})

test_that("an unsupported species fails before anything is written", {
  dest <- file.path(tempdir(), "gatom-unsupported")
  expect_error(
    download_gatom_references(dest_dir = dest, species = "Xenopus"),
    "`species` must be \"Mus_musculus\" or \"Homo_sapiens\"; got \"Xenopus\"",
    fixed = TRUE
  )
  expect_false(dir.exists(dest))
})

test_that("networks must be a non-empty character vector", {
  dest <- file.path(tempdir(), "gatom-badnetworks")
  expect_error(
    download_gatom_references(dest_dir = dest, networks = character()),
    "`networks` must be a non-empty character vector"
  )
  expect_error(
    download_gatom_references(dest_dir = dest, networks = NA_character_),
    "`networks` must be a non-empty character vector"
  )
  expect_false(dir.exists(dest))
})

test_that("gsdb_load points GATOM users at the downloader", {
  expect_error(gsdb_load("gatom"),
               "download_gatom_references", fixed = TRUE)
  expect_error(gsdb_load("GATOM"), "not bundled")
})
