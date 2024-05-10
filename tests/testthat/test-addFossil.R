test_that("addFossil works", {
  mytree <- phytools::pbtree(n=3)
  newtree <- addFossil(mytree, mintime = max(mytree$edge.length)/2, maxtime= max(mytree$edge.length))
  newtree2 <- addFossil(mytree, mintime = max(mytree$edge.length)/2, maxtime= max(mytree$edge.length),edge = 1)
  expect_equal(length(newtree$tip.label), 4)
  expect_equal(length(newtree2$tip.label), 4)
})
