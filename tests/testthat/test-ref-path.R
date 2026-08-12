test_that(".ref_path resolves current and records its snapshot", {
  .clear_ref_resolutions()
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot <- file.path(source_dir, "coresh_snapshot_test")
  dir.create(snapshot, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  withr::local_envvar(c(REFCACHE_ROOT = root))

  out <- .ref_path("coresh")

  expect_identical(as.vector(out), file.path(source_dir, "current"))
  expect_identical(attr(out, "snapshot"), basename(snapshot))
  expect_identical(attr(out, "source"), "coresh")
  expect_identical(attr(out, "root"), root)
  expect_false(attr(out, "caller_supplied"))
})

test_that(".ref_path appends path components below current", {
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot <- file.path(source_dir, "coresh_snapshot_test")
  leaf <- file.path(snapshot, "chunks", "part-01")
  dir.create(leaf, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  withr::local_envvar(c(REFCACHE_ROOT = root))

  out <- .ref_path("coresh", "chunks", "part-01")

  expect_identical(
    as.vector(out),
    file.path(source_dir, "current", "chunks", "part-01")
  )
  expect_identical(attr(out, "snapshot"), basename(snapshot))
})

test_that("an explicit path wins and is marked caller-supplied", {
  root <- withr::local_tempdir()
  explicit <- file.path(root, "pinned-copy")
  dir.create(explicit)
  withr::local_envvar(c(REFCACHE_ROOT = file.path(root, "unused")))

  out <- .ref_path("coresh", "ignored-component", path = explicit)

  expect_identical(as.vector(out), explicit)
  expect_identical(attr(out, "snapshot"), NA_character_)
  expect_identical(attr(out, "root"), NA_character_)
  expect_true(attr(out, "caller_supplied"))
})

test_that("an unset REFCACHE_ROOT names both remedies", {
  withr::local_envvar(c(REFCACHE_ROOT = NA_character_))

  expect_error(
    .ref_path("coresh"),
    "REFCACHE_ROOT.*path="
  )
})

test_that("a missing source and missing current have different errors", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(REFCACHE_ROOT = root))

  expect_error(.ref_path("coresh"), "source directory is missing")

  source_dir <- file.path(root, "coresh")
  dir.create(source_dir)
  expect_error(.ref_path("coresh"), "current.*missing or dangling")
})

test_that("a dangling current identifies a failed refresh", {
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  dir.create(source_dir)
  withr::local_envvar(c(REFCACHE_ROOT = root))

  linked <- suppressWarnings(
    file.symlink(
      file.path(source_dir, "absent-snapshot"),
      file.path(source_dir, "current")
    )
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  expect_error(.ref_path("coresh"), "refcache refresh may have failed")
})

test_that("current may be a real directory", {
  root <- withr::local_tempdir()
  current <- file.path(root, "coresh", "current")
  dir.create(current, recursive = TRUE)
  withr::local_envvar(c(REFCACHE_ROOT = root))

  out <- .ref_path("coresh")

  expect_identical(as.vector(out), current)
  expect_identical(attr(out, "snapshot"), basename(current))
})

test_that("source must be a single safe path segment", {
  expect_error(.ref_path("../etc", path = tempdir()), "single path segment")
})

test_that("the session record keeps one row per source and clears", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  first <- file.path(root, "first")
  second <- file.path(root, "second")
  other <- file.path(root, "other")
  dir.create(first)
  dir.create(second)
  dir.create(other)

  .ref_path("coresh", path = first)
  .ref_path("coresh", path = second)
  .ref_path("cistarget", path = other)
  record <- .ref_resolutions()

  expect_identical(record$source, c("cistarget", "coresh"))
  expect_equal(nrow(record), 2L)
  expect_identical(record$path[record$source == "coresh"], second)

  .clear_ref_resolutions()
  expect_equal(nrow(.ref_resolutions()), 0L)
})
