test_that("nodeEstimate works", {
  skip_on_cran()
  #load test data
  load(test_path("testdata","testocc.rda"))
  load(test_path("testdata","testtree.rda"))
  load(test_path("testdata","testclim.rda"))
  #set up test
  bounds <- list(a = c(min = -1, max = 5), delta = c(min = 0, max = 100))
  bio <- getBioclimVars(testocc,which.biovars=1)
  sp <- tapply(bio[,4],bio$Species,min)
  tre <- geiger::treedata(testtree[[1]],sp)
  ne <- nodeEstimate(tre,1,model="BM")
  ne1 <- nodeEstimate(tre,1,model="estimate",bounds=bounds)
  #runtests
  expect_equal(ne$model,"BM")
  expect_length(ne$fitted,7)
  expect_length(ne1$model,1)
  })


