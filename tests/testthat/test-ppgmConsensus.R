test_that("ppgmConsensus works with original paleoclimate", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  testtree <- testtree[[1]]
  bounds <- list(a = c(min = -1, max = 5), delta = c(min = 0, max = 100))
  #test models
  tBM  <- ppgmConsensus(testocc, trees=testtree, model="BM", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim)
  #test est
  tEs <- ppgmConsensus(testocc, trees=testtree, model="estimate", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, bounds=bounds)
  #test fossil
  tfos <- ppgmConsensus(testocc, trees=testtree, fossils=testfos, which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, plot.GeoRates=TRUE)
  #checking outputs
  expect_equal(length(tBM$node_est[[1]]),4)
  expect_equal(length(tfos$node_est[[1]]),4)
})


test_that("ppgmConsensus works with new paleoclimate", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  testtree <- testtree[[1]]
  newclim <- list(testclim[[2]],testclim[[5]],testclim[[11]])
  layerAge <- c(1,4,10)
  #test BM
  tBM  <- ppgmConsensus(testocc, trees=testtree, model="BM", which.biovars=1, use.paleoclimate=F, paleoclimateUser=newclim, layerAge=layerAge)
  #test fossil
  tfos <- ppgmConsensus(testocc, trees=testtree, fossils=testfos, which.biovars=1, use.paleoclimate=F, paleoclimateUser=newclim, layerAge=layerAge)
  #checking outputs
  expect_equal(length(tBM$node_est[[1]]),4)
  expect_equal(length(tfos$node_est[[1]]),4)
})