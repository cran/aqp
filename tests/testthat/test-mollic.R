context("estimate mollic epipedon bounds")
# construct a fake profile
spc <- data.frame(
  id = 1,
  taxsubgrp = "Lithic Haploxeralfs",
  hzname   = c("A", "AB", "Bt", "BCt", "R"),
  hzdept   = c(0,  20, 32, 42,  49),
  hzdepb   = c(20, 32, 42, 49, 200),
  prop     = c(18, 22, 28, 24,  NA),
  texcl    = c("l", "l", "cl", "l", "br"),
  d_value  = c(5,   5,  5,  6,  NA),
  m_value  = c(2.5, 3,  3,  4,  NA),
  m_chroma = c(2,   3,  4,  4,  NA)
)

spc2 <- data.frame(
  id = c(1, 1, 1, 1),
  hzname = c("A", "Bt1", "Bt2", "Bk"),
  hzdept = c(0L, 5L, 16L, 36L),
  hzdepb = c(5L, 16L, 36L, 66L),
  claytotest = c(11L, 30L, 32L, 11L),
  texcl = c("sl", "scl", "scl", "sl"),
  d_value = c(6L, 4L, 4L, 4L),
  m_chroma = c(2L, 2L,  3L, 3L)
)

# promote to SoilProfileCollection
depths(spc) <- id ~ hzdept + hzdepb
hzdesgnname(spc) <- 'hzname'
hztexclname(spc) <- 'texcl'
hzmetaname(spc, "clay") <- 'prop'

depths(spc2) <- id ~ hzdept + hzdepb
hzdesgnname(spc2) <- 'hzname'
hztexclname(spc2) <- 'texcl'
hzmetaname(spc2, "clay") <- 'claytotest'

spc3 <- data.frame(pedon_key = c("10016", "10016", "10016", "10047",  "10047",
                                   "10047", "10047", "10047", "10047", "10047"), 
                   hzn_top = c(0, 18, 51, 0, 26, 45, 55, 90, 102, 130), 
                   hzn_bot = c(18, 51, 107, 26, 45, 55, 90, 102, 130, 185), 
                   hzn_desgn = c("A", "C1",  "C2", "Ap", "E", "Btn", "Bt", 
                                 "Eb", "2Btk1", "2Btk2"), 
                   texture_lab = c("sl", "sil", "sil", "sil", "sil", 
                                   "sil", "sil", "sil", "sil", "sil"), 
                   clay_total = c(4.4, 7.7, 11.9, 16.7, 9.9, 
                                  19.5, 14.2, 7.3, 11.3, 7.6))
depths(spc3) <-pedon_key ~ hzn_top + hzn_bot
hzdesgnname(spc3) <- "hzn_desgn"
hztexclname(spc3) <- "texture_lab"
hzmetaname(spc3, "clay") <- "clay_total"

test_that("mollic.thickness.requirement", {
  expect_equal(mollic.thickness.requirement(spc, clay.attr = 'prop'), 18)
  expect_equal(mollic.thickness.requirement(spc, clay.attr = 'prop', truncate = FALSE), 49 / 3)
  expect_equal(mollic.thickness.requirement(trunc(spc, 0, 9), clay.attr = 'prop'), 10)

  # this is a controversial/undefined case:
  # 
  #  the Bk is identified as a cambic below an argillic 
  #  
  #  so most limiting crit is the bottom depth of cambic: 66/3 = 22 
  #  rather than 36/3 = 12 truncated to 18 which comes from argillic bottom depth/2nd carbonate upperbound
  expect_equal(mollic.thickness.requirement(spc2, clay.attr = 'claytotest'), 22)
  
  expect_equal(mollic.thickness.requirement(spc3), c(18, 25))
})

