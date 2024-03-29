context("soil profile simulation")

## sample data
data(sp3)
depths(sp3) <- id ~ top + bottom

# select a profile to use as the basis for simulation
s <- sp3[3, ]

# reset horizon names
s$name <- paste('H', seq_along(s$name), sep = '')
test_that("perturb() (by thickness) works as expected (single profile)", {

  # simulate 25 new profiles
  expect_warning({sim.1 <- sim(s, n = 25)}) 
  #deprecated version creates sd column for perturb using hz.sd=2

  # manually create hz.sd for perturb()
  s$hz.sd = c(1, 2, 5, 5, 5, 10, 3)
  expect_silent({sim.2 <- perturb(s, n = 25, thickness.attr = "hz.sd")})

  # result is an SPC
  expect_true(inherits(sim.1, 'SoilProfileCollection'))
  expect_true(inherits(sim.2, 'SoilProfileCollection'))

  # expected lengths
  expect_true(length(sim.1) == 25)
  expect_true(length(sim.2) == 25)
})

test_that("expected errors (single profile)", {

  # only 1 seed can be used
  expect_error(perturb(sp3[1:2, ], n = 25))

  # sd must recycle evenly over number of original horizons

  # NOTE: test removed; this recycling is now handled by horizons()<-
  # in the deprecated version of sim()
})

test_that("perturb (by boundaries) works as expected (single profile)", {
  # simulate 25 new profiles with a sd boundary thickness of 0.5 - 2.5cm
  s$bdy <- round(runif(nrow(s), 1, 5)) / 2
  diagnostic_hz(s) <- data.frame(id = profile_id(s),
                                 featkind = "foo",
                                 featdept = 0, featdepb = 10)
  restrictions(s) <- data.frame(id = profile_id(s),
                                 restrkind = "bar",
                                 restrdept = 0, restrdepb = 10)
  perp <- perturb(s, n = 25, boundary.attr = "bdy")

  # result is an SPC
  expect_true(inherits(perp, 'SoilProfileCollection'))

  # expected lengths
  expect_true(length(perp) == 25)

  perp2 <- perturb(s, id = 26:50, boundary.attr = "bdy", new.idname = "foo")

  # result is an SPC
  expect_true(inherits(perp2, 'SoilProfileCollection'))

  # expected lengths
  expect_true(length(perp2) == 25)

  # custom idname
  expect_equal(idname(perp2), "foo")

  # custom IDs
  expect_equal(profile_id(perp2), as.character(26:50))
})

test_that("perturb() (by thickness) works as expected (multiple profiles)", {
  
  horizons(sp3)$hz.sd <- 2
  expect_silent({ sim.1 <- perturb(sp3, n = 25, thickness.attr = "hz.sd") })
  expect_equal(length(sim.1), 250)
  
  # test min.thickness arg (no horizon thinner than min.thickness)
  thk <- sim.1$bottom - sim.1$top
  expect_true(min(thk) >= 1)
  
  # all results are logical
  expect_true(all(checkHzDepthLogic(sim.1)$valid))
})

test_that("perturb() (by boundary) works as expected (multiple profiles)", {
  
  horizons(sp3)$hz.sd <- 2
  expect_silent({ sim.1 <- perturb(sp3[9:10,], n = 250, boundary.attr = "hz.sd", min.thickness = 3, max.depth = 50) })
  expect_equal(length(sim.1), 500)
  
  # test min.thickness arg (no horizon thinner than min.thickness)
  thk <- sim.1$bottom - sim.1$top
  expect_true(min(thk) >= 3)
  
  # all results are logical
  expect_true(all(checkHzDepthLogic(sim.1)$valid))
  
  # test max.depth (thicknesses of 4th horizon below 50cm are unaltered)
  expect_equal(thk[seq_along(thk) %% 4 == 0], c(rep(34, 250), rep(22, 250)))
})

