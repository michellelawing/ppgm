test_that("getBioclimVars works modern occ", {
  load(test_path("testdata","testocc.rda"))
  bio <- getBioclimVars(testocc,which.biovars=c(2,4))
  expect_length(bio,5)
  occ2 <- testocc[,1:3]
  bio2 <- getBioclimVars(occ2,which.biovars=c(2,4))
  expect_length(bio2,5)
})

test_that("getBioclimVars works fossilocc", {
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  testfos <- testfos[,-1,drop=F]
  bio <- getBioclimVars(testfos,which.biovars=c(2,4))
  expect_length(bio,5)
})

