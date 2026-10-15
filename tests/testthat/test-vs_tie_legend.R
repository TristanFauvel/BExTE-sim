make_vs_tie_legend_plot <- function() {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  methods <- c('p-PP','Separate','RMP','Pooling','NPP','EBPP',
    'Conditional PP','Com. PP','Com. prior','Empirical RMP','NPP (KL)')
  p <- ggplot2::ggplot()
  for (method in methods) {
   n <- switch(method, RMP=9, 'Conditional PP'=3, 'Com. PP'=3,
     'Com. prior'=3, 'NPP (KL)'=2, 1)
   params <- if (method %in% c('Pooling','Separate','EBPP','Empirical RMP'))
     list(expression('')) else lapply(seq_len(n),function(i) bquote(w == .(i/10)))
   d <- data.frame(x=seq_len(n)/20,y=seq_len(n)/10,k=factor(seq_len(n)))
   p <- p + ggplot2::geom_point(data=d,ggplot2::aes(x,y,colour=k),shape=match(method,methods)) + ggplot2::scale_colour_discrete(name=if(n==1) NULL else method,labels=vs_tie_key_labels(method,params),guide=ggplot2::guide_legend(ncol=min(n,3),byrow=TRUE)) + ggnewscale::new_scale_color()
  }
  p <- p + ggplot2::theme_bw(base_size=12) + ggplot2::theme(legend.key.size=grid::unit(.1,'cm'),legend.margin=ggplot2::margin(2,2,2,2))
  list(plot = p, methods = methods)
}

test_that("legend placement follows method names, including inline math labels", {
  fixture <- make_vs_tie_legend_plot()
  legend <- vs_tie_legend_grid(fixture$plot, rev(fixture$methods))
  expected <- PAPER_VS_TIE_LEGEND_LAYOUT
  for (i in seq_along(legend$grobs)) {
    texts <- vs_tie_legend_text(legend$grobs[[i]])
    cell <- legend$layout[i, ]
    method <- expected[cell$t, cell$l]
    expect_true(any(texts == method | startsWith(texts, paste0('"', method, '" ~'))))
    expect_true(all(expected[cell$t:cell$b, cell$l:cell$r] == method))
  }
  expect_length(legend$grobs, 11)
  expect_equal(nrow(expected), 6)
  expect_identical(expected[4, 3], "Empirical RMP")
  expect_identical(expected[5, 3], "NPP (KL)")
  ## Test-then-pool gets the last row: its six equivalence settings span two
  ## columns, beside the four difference settings.
  expect_identical(expected[6, ], c("TtP (diff.)", "TtP (eq.)", "TtP (eq.)"))
})

test_that("spanning legends fit without counting their size repeatedly", {
  skip_on_cran()
  fixture <- make_vs_tie_legend_plot()
  legend <- vs_tie_legend_grid(fixture$plot, fixture$methods)
  for (i in seq_along(legend$grobs)) {
    cell <- legend$layout[i, ]
    g <- legend$grobs[[i]]$children[[1]]
    expect_gte(sum(as.numeric(legend$widths[cell$l:cell$r])) + 1e-8,
      grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE))
    expect_gte(sum(as.numeric(legend$heights[cell$t:cell$b])) + 1e-8,
      grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE))
  }
  composed <- vs_tie_compose(fixture$plot, fixture$methods, 1, 2)
  expect_gte(composed$width, sum(as.numeric(legend$widths)))
})

test_that("missing methods leave no empty rows and extra methods are retained", {
  p <- ggplot2::ggplot(data.frame(x = 1, k = "a"), ggplot2::aes(x, x, colour = k)) +
    ggplot2::geom_point() + ggplot2::scale_colour_discrete(name = "Extra")
  legend <- vs_tie_legend_grid(p, c("Pooling", "Extra"))
  expect_length(legend$grobs, 1)
  expect_length(legend$heights, 1)
  expect_null(vs_tie_legend_grid(p + ggplot2::guides(colour = "none"), "Extra"))
})
