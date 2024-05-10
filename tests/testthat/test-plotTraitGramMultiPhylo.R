test_that("plotTraitGram-MultiPhylo works", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  #test BM
  tBM  <- ppgm(testocc, trees=testtree, fossils=testfos, model="BM", which.biovars=1, use.paleoclimate=FALSE, paleoclimateUser=testclim)
  expect_snapshot(plotTraitGramMultiPhylo(tBM$treedata_min,tBM$treedata_max,tBM$node_est,fossils=testfos,which.biovars=1,use.paleoclimate=FALSE,paleoclimateUser = testclim))
})
