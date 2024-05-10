test_that("ppgm works with original paleoclimate", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  bounds <- list(a = c(min = -1, max = 5), delta = c(min = 0, max = 1000))
  #test models
  tBM  <- ppgm(testocc, trees=testtree, model="BM", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim)
  tOU  <- ppgm(testocc, trees=testtree, model="OU", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim)
  tEB  <- ppgm(testocc, trees=testtree, model="EB", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim)
  tLa  <- ppgm(testocc, trees=testtree, model="lambda", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, plot.GeoRates = TRUE)
  tKa  <- ppgm(testocc, trees=testtree, model="kappa", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, plot.TraitGram = TRUE)
  tDe  <- ppgm(testocc, trees=testtree, model="delta", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, verbose=FALSE)
  #test est
  tEs <- ppgm(testocc, trees=testtree, model="estimate", which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim, bounds=bounds)
  #test fossil
  tfos <- ppgm(testocc, trees=testtree, fossils=testfos, which.biovars=1, use.paleoclimate=F, paleoclimateUser=testclim)
  #checking outputs
  expect_equal(length(tBM$node_est),length(testtree))
  expect_equal(length(tfos$node_est),length(testtree))
  expect_equal(length(tEs$node_est),length(testtree))
})


test_that("ppgm works with new paleoclimate", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  newclim <- list(testclim[[2]],testclim[[5]],testclim[[11]])
  layerAge <- c(1,4,10)
  bounds <- list(a = c(min = -1, max = 5), delta = c(min = 0, max = 1000))
  #test BM
  tBM  <- ppgm(testocc, trees=testtree, model="BM" ,which.biovars=1, use.paleoclimate=F, paleoclimateUser=newclim, layerAge=layerAge)
  #test fossil
  tfos <- ppgm(testocc, trees=testtree, fossils=testfos, which.biovars=1, use.paleoclimate=F, paleoclimateUser=newclim, layerAge=layerAge)
  #checking outputs
  expect_equal(length(tBM$node_est),length(testtree))
  expect_equal(length(tfos$node_est),length(testtree))
})