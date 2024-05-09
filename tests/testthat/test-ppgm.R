test_that("ppgm works with original paleoclimate", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  #test BM
  tBM  <- ppgm(testocc ,fossils=FALSE, trees=testtree, model="BM", which.biovars=1, paleoclimateUser=testclim)
  #test fossil
  tfos <- ppgm(testocc, trees=testtree, fossils=testfos, which.biovars=1, paleoclimateUser=testclim)
  #checking outputs
  expect_equal(length(tBM$node_est),length(testtree))
  expect_equal(length(tfos$node_est),length(testtree))
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
  #test BM
  tBM  <- ppgm(testocc,fossils=FALSE,trees=testtree,model="BM",which.biovars=1,paleoclimateUser=newclim,layerAge=layerAge)
  #test fossil
  tfos <- ppgm(testocc,trees=testtree,fossils=testfos,which.biovars=1,paleoclimateUser=newclim,layerAge=layerAge)
  #checking outputs
  expect_equal(length(tBM$node_est),length(testtree))
  expect_equal(length(tfos$node_est),length(testtree))
})