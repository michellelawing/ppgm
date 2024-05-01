test_that("ppgm works", {
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  #test BM
  tBM  <- ppgm(testocc,fossils=FALSE,trees=testtree,model="BM",which.biovars=1,paleoclimateUser=testclim)
  #test estimate
  tE   <- ppgm(testocc,fossils=FALSE,trees=testtree,model="estimate",which.biovars=1,paleoclimateUser=testclim)
  #test fossil
  tfos <- ppgm(testocc,trees=testtree,fossils=testfos,which.biovars=1,paleoclimateUser=testclim)
  #checking outputs
  expect_equal(length(tBM$node_est),length(testtree))
  expect_equal(length(tE$node_est),length(testtree))
  expect_equal(length(tfos$node_est),length(testtree))
  
  
})

