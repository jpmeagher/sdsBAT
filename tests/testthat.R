library(testthat)
library(sdsBAT)

test_that("pou_covariance produces a valid covariance function", {
  h <-  log(runif(3, 0, 1))

  expect_equal(pou_covariance(h, phylogeny), t(pou_covariance(h, phylogeny)))
  expect_equal(dim(pou_covariance(h, phylogeny)), rep(length(phylogeny$tip.label), 2))
  expect_equal(all(eigen(pou_covariance(h, phylogeny), symmetric = TRUE)$values > 0), TRUE)

  remove(h)
})

test_that("pou_simulated_data produces data of the correct dimension", {
  h <-  log(runif(3, 0, 1))

  expect_equal(dim(pou_simulated_data(h, phylogeny, 10)), c(length(phylogeny$tip.label), 10))

  remove(h)
})

test_that("pou_logl_slow and pou_logl_fast behave themselves", {
  h <-  log(runif(3, 0, 1))
  s <- pou_simulated_data(h, phylogeny, 1)

  expect_equal(pou_logl_slow(h, phylogeny, s) > 0 , TRUE)
  expect_equal(pou_logl_fast(h, phylogeny, s) > 0 , TRUE)
  expect_equal(round(pou_logl_slow(h, phylogeny, s), 2), round(pou_logl_fast(h, phylogeny, s), 2))

  S <- pou_simulated_data(h, phylogeny, 10)

  expect_equal(all(apply(S, 2, pou_logl_slow, ln_hyperparameters = h, phylogenetic_tree = phylogeny) > 0), TRUE)
  expect_equal(all(apply(S, 2, pou_logl_fast, ln_hyperparameters = h, phylogenetic_tree = phylogeny) > 0), TRUE)
  expect_equal(round(apply(S, 2, pou_logl_slow, ln_hyperparameters = h, phylogenetic_tree = phylogeny), digits = 2) , round(apply(S, 2, pou_logl_fast, ln_hyperparameters = h, phylogenetic_tree = phylogeny), digits = 2))

  remove(h, s, S)
})

test_check("sdsBAT")
