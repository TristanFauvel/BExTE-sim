## NAMESPACE and man/ are generated from the roxygen blocks in R/, while the
## reference index in _pkgdown.yml is maintained by hand. Nothing in an R
## session reads any of them, so when they drift apart the only thing that
## notices is a pkgdown build - which publishes the site, and fails long after
## the commit that broke it. Both directions had already drifted on main:
## vs_tie_legend_text carried an @export that NAMESPACE had never been
## regenerated to include, and NAMESPACE still exported
## operating_characteristic_vs_tie although no roxygen block was left to
## document it - so the next devtools::document() deleted its help topic and
## took the docs build down with it, because _pkgdown.yml still pointed there.

project_dir <- testthat::test_path("..", "..")

## These checks read the package sources, which an installed package (as
## under R CMD check) does not have.
skip_if_not(
  file.exists(file.path(project_dir, "_pkgdown.yml")),
  "needs the package source tree"
)

namespace_exports <- function() {
  lines <- readLines(file.path(project_dir, "NAMESPACE"), warn = FALSE)
  sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", lines, value = TRUE))
}

## The S3 methods NAMESPACE registers, named generic.class as in R/.
namespace_s3_methods <- function() {
  lines <- readLines(file.path(project_dir, "NAMESPACE"), warn = FALSE)
  methods <- grep("^S3method\\(", lines, value = TRUE)
  sub("^S3method\\(([^,]+),([^)]+)\\)$", "\\1.\\2", methods)
}

## The object an `@export` tag applies to: the first thing assigned after the
## roxygen block it sits in.
roxygen_exports <- function() {
  r_files <- list.files(file.path(project_dir, "R"), pattern = "\\.R$",
                        full.names = TRUE)
  names <- unlist(lapply(r_files, function(path) {
    lines <- readLines(path, warn = FALSE)
    tags <- grep("^#+' *@export *$", lines)
    vapply(tags, function(i) {
      j <- i + 1
      while (j <= length(lines) && grepl("^#", lines[j])) j <- j + 1
      if (j > length(lines)) return(NA_character_)
      matched <- regmatches(
        lines[j], regexec("^([A-Za-z._][A-Za-z0-9._]*)\\s*<-", lines[j])
      )[[1]]
      if (length(matched) < 2) NA_character_ else matched[2]
    }, character(1))
  }))
  names
}

pkgdown_index <- function() {
  config <- yaml::read_yaml(file.path(project_dir, "_pkgdown.yml"))
  unlist(lapply(config$reference, function(section) section$contents))
}

help_topics <- function() {
  rd_files <- list.files(file.path(project_dir, "man"), pattern = "\\.Rd$",
                         full.names = TRUE)
  unique(unlist(lapply(rd_files, function(path) {
    lines <- readLines(path, warn = FALSE)
    sub("^\\\\alias\\{(.*)\\}$", "\\1",
        grep("^\\\\alias\\{", lines, value = TRUE))
  })))
}

test_that("the files these checks read are all there", {
  # A silent empty read would make every assertion below vacuously true.
  expect_gt(length(namespace_exports()), 15)
  expect_gt(length(roxygen_exports()), 15)
  expect_gt(length(pkgdown_index()), 80)
  expect_gt(length(help_topics()), 100)
})

test_that("every @export tag resolves to the object it is meant to export", {
  # Guards the parsing the two checks below rely on, rather than the package.
  expect_equal(sum(is.na(roxygen_exports())), 0)
})

test_that("NAMESPACE exports exactly what the roxygen blocks ask it to", {
  # Either direction means NAMESPACE is stale: run devtools::document().
  tagged <- roxygen_exports()
  tagged <- tagged[!is.na(tagged)]
  # An @export on an S3 method registers it rather than exporting it.
  tagged <- setdiff(tagged, namespace_s3_methods())
  exported <- namespace_exports()

  expect_equal(
    setdiff(tagged, exported), character(0),
    info = "tagged @export in R/ but missing from NAMESPACE"
  )
  expect_equal(
    setdiff(exported, tagged), character(0),
    info = "exported by NAMESPACE with no @export left in R/"
  )
})

test_that("every exported function is listed in the pkgdown reference index", {
  missing <- setdiff(namespace_exports(), pkgdown_index())

  expect_equal(
    missing, character(0),
    info = paste("add to _pkgdown.yml:", paste(missing, collapse = ", "))
  )
})

test_that("every name in the pkgdown reference index has a help topic", {
  orphans <- setdiff(pkgdown_index(), help_topics())

  expect_equal(
    orphans, character(0),
    info = paste("no help topic for:", paste(orphans, collapse = ", "))
  )
})

test_that("the reference index lists nothing twice", {
  index <- pkgdown_index()
  repeated <- unique(index[duplicated(index)])

  expect_equal(
    repeated, character(0),
    info = paste("listed more than once:", paste(repeated, collapse = ", "))
  )
})
