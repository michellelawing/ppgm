test_that("ppgmMESS works", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  #test BM
  tBM  <- ppgm(testocc, trees=testtree, model="BM", which.biovars=c(1,2), use.paleoclimate = FALSE, paleoclimateUser=testclim)
  c_mi <- cbind(tBM$cem[,1],tBM$cem[,2])
  c_ma <- cbind(tBM$cem[,5],tBM$cem[,6])
  rownames(c_mi) <- rownames(c_ma) <- rownames(tBM$cem)
  sco <- ppgmMESS(c_mi,c_ma,tBM$node_est,tree=testtree,timeslice=10,which.biovars=c(1,2),which.plot="all",use.paleoclimate=FALSE,paleoclimateUser=testclim)
  #test BMfos
  tfos  <- ppgm(testocc, trees=testtree, model="BM", which.biovars=2, fossils=testfos)
  c_mi <- cbind(tfos$cem[,1])
  c_ma <- cbind(tfos$cem[,3])
  rownames(c_mi) <- rownames(c_ma) <- rownames(tfos$cem)
  sco2 <- ppgmMESS(c_mi,c_ma,tfos$node_est,tree=testtree,timeslice=10,which.biovars=1,fossils=testfos,which.plot="all")
  #checking outputs
  expect_length(sco[[1]],210)
  expect_length(sco2[[1]],9756)
})