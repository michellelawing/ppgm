test_that("getLineageClimate works", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testfos.rda"))
  load(test_path("testdata","testclim.rda"))
  bio <- getBioclimVars(testocc,which.biovars=2)
  min<- tapply(bio[,4],bio$Species,min)
  max<- tapply(bio[,4],bio$Species,max)
  tdmin <- treedata(testtree[[1]],min,sort=TRUE,warnings=FALSE)
  tdmax <- treedata(testtree[[1]],max,sort=TRUE,warnings=FALSE)
  est <- nodeEstimateEnvelopes(tdmin,tdmax)
  est2 <- est$est
  env <- getEnvelopes(tdmin,tdmax,est2)
  linCl <- getLineageClimate(env,testtree[[1]],which.biovars=2)
  expect_length(linCl,2)
  })
