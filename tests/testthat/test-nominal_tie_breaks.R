## The nominal TIE tick is the reference the vs-TIE plots are read against,
## so it has to survive break selection and label drawing.

axis_label_grob <- function(plot) {
  gtable <- ggplot2::ggplotGrob(plot)
  axis <- gtable$grobs[[which(gtable$layout$name == "axis-b")]]
  collect <- function(grob) {
    if (inherits(grob, "text")) return(list(grob))
    children <- c(grob$grobs, as.list(grob$children))
    unlist(lapply(children, collect), recursive = FALSE)
  }
  collect(axis)[[1]]
}


test_that("nominal_tie_breaks keeps the nominal value among the breaks", {
  breaks <- nominal_tie_breaks(0.025)(c(0, 0.08))

  expect_true(0.025 %in% breaks)
})


test_that("nominal_tie_breaks drops pretty breaks that would collide with the nominal value", {
  breaks <- nominal_tie_breaks(0.025)(c(0, 0.08))

  expect_false(0.02 %in% breaks)
  expect_gte(min(abs(setdiff(breaks, 0.025) - 0.025)), 0.5 * min(diff(pretty(c(0, 0.08)))))
})


test_that("nominal_tie_breaks keeps pretty breaks that are clear of the nominal value", {
  breaks <- nominal_tie_breaks(0.025)(c(0, 0.08))

  expect_true(all(c(0.04, 0.06, 0.08) %in% breaks))
})


test_that("nominal_tie_breaks returns sorted breaks", {
  breaks <- nominal_tie_breaks(0.025)(c(0, 0.08))

  expect_identical(breaks, sort(breaks))
})


test_that("nominal_tie_breaks falls back to pretty breaks when no nominal value is given", {
  limits <- c(0, 0.08)

  expect_identical(nominal_tie_breaks(NULL)(limits), pretty(limits))
  expect_identical(nominal_tie_breaks(NA_real_)(limits), pretty(limits))
})


test_that("an axis using nominal_tie_breaks draws the nominal TIE label", {
  df <- data.frame(tie = c(0.005, 0.02, 0.03, 0.08), power = c(0.2, 0.5, 0.6, 0.8))
  plot <- ggplot2::ggplot(df, ggplot2::aes(tie, power)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = nominal_tie_breaks(0.025))

  labels <- axis_label_grob(plot)

  expect_true("0.025" %in% as.character(labels$label))
  ## check.overlap lets grid discard a label at draw time, and the nominal
  ## tick is the one it discards: it sits next to a pretty break and loses
  ## the priority contest against it.
  expect_false(isTRUE(labels$check.overlap))
})
