test_that("animatedppgmMultiPhylo works", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testclim.rda"))
  testclim <- list(testclim[[1]],testclim[[2]],testclim[[3]])
  #test models
  tBM  <- ppgm(testocc, trees=testtree, model="BM", which.biovars=c(1,4,15), use.paleoclimate=F, paleoclimateUser=testclim, layerAge=c(0:2))
  expect_snapshot(plotAnimatedPPGMMultiPhylo(tBM$envelope, tree=testtree, which.biovars=c(1,4,15), use.paleoclimate=F, paleoclimateUser=testclim, layerAge=c(0:2)))
})

