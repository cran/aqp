context("color conversion")

## sample data
x <- '10YR 3/4'
x.p <- parseMunsell(x, convertColors = FALSE)

# Munsell --> sRGB conversion
m <- munsell2rgb(x.p$hue, x.p$value, x.p$chroma)
m.rgb <- munsell2rgb(x.p$hue, x.p$value, x.p$chroma, return_triplets = TRUE)

# sRGB --> Munsell
x.back <- col2Munsell(col = m.rgb, nClosest = 1)
# using truncated sRGB values
x.back.trunc <- col2Munsell(col = data.frame(r = 0.36, g = 0.26, b = 0.13))

# neutral colors map to shades of gray
x.neutral <- parseMunsell('N 2/', return_triplets = TRUE)

## tests

test_that("col2Munsell works as expected", {
  
  # 10YR 3/4 as hex notation sRGB 
  .res <- col2Munsell(m)
  expect_equal(.res$hue, '10YR')
  expect_equal(.res$value, 3L)
  expect_equal(.res$chroma, 4L)
  
  
  # N 2/ as sRGB [0,1]
  .res <- col2Munsell(x.neutral, space = 'sRGB')
  expect_equal(.res$hue, 'N')
  expect_equal(.res$value, 2L)
  expect_equal(.res$chroma, 0L)
  
})



test_that("parseMunsell()", {
  
  # parsing bogus notation generates NA
  # will also generate a warning from munsell2rgb()
  expect_equal(suppressWarnings(parseMunsell('10YZ 4/5')), NA_character_)
  expect_equal(suppressWarnings(parseMunsell('10YR /5')), NA_character_)
  expect_equal(suppressWarnings(parseMunsell('10YR ')), NA_character_)
  expect_equal(suppressWarnings(parseMunsell('10YR 4/')), NA_character_)
  expect_equal(suppressWarnings(parseMunsell('G1 6/N')), NA_character_)
  
  # parsing bogus notation without conversion
  bogus <- parseMunsell('G1 3/X', convertColors = FALSE)
  expect_equal(bogus$hue, NA_character_)
  expect_equal(bogus$value, NA_real_)
  expect_equal(bogus$chroma, NA_real_)
  
  # test NA
  some.NA <- parseMunsell(c(NA, '10YR 3/3'))
  expect_true(inherits(some.NA, 'character'))
  expect_true(length(some.NA) == 2)
  
  some.NA <- parseMunsell(c(NA, '10YR 3/3'), convertColors = FALSE)
  expect_true(inherits(some.NA, 'data.frame'))
  expect_true(nrow(some.NA) == 2)
  
  # splitting of text into columns within data.frame
  expect_identical(x.p, data.frame(hue = "10YR", value = 3, chroma = 4, stringsAsFactors = FALSE))
  
  # Test not using spaces
  expect_equal(suppressWarnings(parseMunsell('2.5YR 3/4')), suppressWarnings(parseMunsell('2.5YR3/4')))
  
  # Test different delimiters
  expect_equal(suppressWarnings(parseMunsell('2.5YR 3/4')), suppressWarnings(parseMunsell('2.5YR 3_4', delim = "_")))
  expect_equal(suppressWarnings(parseMunsell('2.5YR 3/4')), suppressWarnings(parseMunsell('2.5YR 3_4')))
})


# addresses #66 (https://github.com/ncss-tech/aqp/issues/66)
test_that("Munsell hue parsing", {
  
  # normal operation
  res <- .parseMunsellHue('10YR')
  expect_true(inherits(res, 'data.frame'))
  expect_equal(res$hue.numeric, 10L)
  expect_equal(res$hue.character, 'YR')
  expect_equal(nrow(res), 1)
  
  # white space is trimmed
  res <- .parseMunsellHue('10 YR')
  expect_true(inherits(res, 'data.frame'))
  expect_equal(res$hue.numeric, 10L)
  expect_equal(res$hue.character, 'YR')
  expect_equal(nrow(res), 1)
  
  # decimal, won't convert correctly, but should be split
  res <- .parseMunsellHue('10.1YR')
  expect_true(inherits(res, 'data.frame'))
  expect_equal(res$hue.numeric, 10.1)
  expect_equal(res$hue.character, 'YR')
  expect_equal(nrow(res), 1)
  
  # bogus hue
  res <- .parseMunsellHue('G1 ')
  expect_true(inherits(res, 'data.frame'))
  expect_true(is.na(res$hue.numeric))
  expect_true(is.na(res$hue.character))
  expect_equal(nrow(res), 1)
})


test_that("non-integer value and chroma are snapped to valid possibilities", {
  
  # rounding of value, throws warning
  expect_warning(res <- parseMunsell('10YR 3.3/4'), regexp = 'non-standard notation')
  
  # this will not throw a warning
  res <- parseMunsell('10YR 3.3/4', convertColors = FALSE)
  
  # results should be the same
  expect_equal(
    suppressWarnings(parseMunsell('10YR 3.3/4')),
    parseMunsell('10YR 3/4')
  )
  
  # 2.5,  8.5, 9.5 are valid Munsell value
  # snap nearby to these
  expect_equal(
    suppressWarnings(parseMunsell('10YR 2.6/4')),
    parseMunsell('10YR 2.5/4')
  )
  
  expect_equal(
    suppressWarnings(parseMunsell('10YR 8.4/4')),
    parseMunsell('10YR 8.5/4')
  )
  
  expect_equal(
    suppressWarnings(parseMunsell('10YR 9.55/4')),
    parseMunsell('10YR 9.5/4')
  )
  
  # rounding of chroma, throws warning
  expect_warning(res <- parseMunsell('10YR 3/4.6'), regexp = 'non-standard notation')
  # this will not throw a warning
  res <- parseMunsell('10YR 3/4.6', convertColors = FALSE)
  
  # results should be the same
  expect_equal(
    suppressWarnings(parseMunsell('10YR 3/4.6')),
    parseMunsell('10YR 3/5')
  )
  
  # double-check no snapping of 2.5 values
  res <- parseMunsell('10YR 2.5/2')
  res.test <- col2Munsell(col = '#493A2BFF')
  
  expect_true(res.test$value == 2.5)
  
  
  # test .snapValid() directly
  .valid <- c(2, 2.5, 3, 4, 5, 6, 7, 8, 8.5, 9, 9.5, 10)
  v <- c(2, 4, 6, 8, NA, 2.5, 2.6, 8.1, 9.6, 12)
  res <- aqp:::.snapValid(v, .valid)
  
  # ensure snapping worked as expected
  expect_equal(
    res, 
    c(2, 4, 6, 8, NA, 2.5, 2.5, 8, 9.5, 10)
  )
  
})


test_that("Munsell <--> sRGB and back again", {
  
  # sRGB in hex notation
  expect_equal(m, '#5C4222FF')
  expect_equal(parseMunsell(x), m)
  
  # sRGB triplets
  expect_equal(m.rgb$r, 0.3618738, tolerance=0.0001)
  expect_equal(m.rgb$g, 0.2598939, tolerance=0.0001)
  expect_equal(m.rgb$b, 0.1337521, tolerance=0.0001)
  
  # neutral colors
  expect_equal(x.neutral$r, 0.03278, tolerance=0.001)
  expect_equal(x.neutral$g, 0.03305, tolerance=0.001)
  expect_equal(x.neutral$b, 0.03305, tolerance=0.001)
  
  # sRGB --> Munsell
  expect_equal(x.back$hue, '10YR')
  expect_equal(x.back$value, 3)
  expect_equal(x.back$chroma, 4)
  expect_equal(x.back$sigma, 0)
  
  expect_equal(x.back.trunc$hue, '10YR')
  expect_equal(x.back.trunc$value, 3)
  expect_equal(x.back.trunc$chroma, 4)
})


test_that("missing data", {
  
  # data with missing sRGB coordinates
  color <- rbind(
    cbind(NA, NA, NA),
    cbind(0.5, 0.2, 0.2),
    cbind(1, 1, 1),
    cbind(NA, NA, NA)
  )
  
  # conversion should work without error
  res <- col2Munsell(color, space = 'sRGB')
  
  # same number of rows in / out
  expect_true(nrow(res) == nrow(color))
  
  # row order preserved
  expect_true(is.na(res$hue[1]) & is.na(res$hue[4]))
  
  # data with missing character representation
  color <- c('red', NA, NA, 'blue')
  res <- col2Munsell(color, space = 'sRGB')
  expect_true(nrow(res) == length(color))
  
  
  # all NA
  color <- c(NA, NA, NA, NA)
  res <- col2Munsell(color, space = 'sRGB')
  expect_true(nrow(res) == length(color))
  
  color <- rbind(
    cbind(NA, NA, NA),
    cbind(NA, NA, NA)
  )
  
  res <- col2Munsell(color, space = 'sRGB')
  expect_true(nrow(res) == nrow(color))
  
})

test_that("neutral hues", {
  
  # neutral colors
  res <- parseMunsell('N 2/', convertColors = FALSE)
  expect_true(inherits(res, 'data.frame'))
  
  # chroma is arbitrarily set to 0
  expect_true(res$chroma == 0)
  
  N2 <- parseMunsell('N 2/')
  N6 <- parseMunsell('N 6/')
  expect_equal(N2, '#080808FF')
  expect_equal(N6, '#464848FF')
  
  # alternative notation
  N2 <- parseMunsell('N 2')
  N6 <- parseMunsell('N 6/0')
  expect_equal(N2, '#080808FF')
  expect_equal(N6, '#464848FF')
  
  # missing, 0, or bogus chroma are all ignored
  N2 <- munsell2rgb('N', 2, NA)
  N4 <- munsell2rgb('N', 4, 0)
  N6 <- munsell2rgb('N', 6, 4)
  expect_equal(N2, '#080808FF')
  expect_equal(N4, '#1B1C1CFF')
  expect_equal(N6, '#464848FF')
  
  # 0 chroma with any hue -> neutral
  N <- munsell2rgb('10YR', 4, 0)
  expect_equal(N, '#1B1C1CFF')
})


test_that("closest Munsell chip based on sRGB coordinates", {
  
  # closest chip in aqp LUT
  expect_equal(getClosestMunsellChip('10YR 3.3/5', convertColors = FALSE), '10YR 3/5')
  expect_equal(getClosestMunsellChip('9YR 3.8/3', convertColors = FALSE), '10YR 4/3')
  
  # 2022-03-30: now we get correct selection of 2.5 value
  expect_equal(getClosestMunsellChip('7.9YR 2.7/2.0', convertColors = FALSE), '7.5YR 2.5/2')
})


# https://github.com/ncss-tech/aqp/issues/69
test_that("Munsell --> LAB + sRGB coordinates", {
  
  # sRGB
  test.1 <- parseMunsell("10YR 3/5", return_triplets=TRUE)
  expect_equal(names(test.1), c('r', 'g', 'b'))
  
  
  # sRGB and LAB
  test.2 <- parseMunsell("10YR 3/5", return_triplets=TRUE, returnLAB=TRUE)
  expect_equal(names(test.2), c('r', 'g', 'b', 'L', 'A', 'B'))
  
  # LAB
  test.3 <- parseMunsell("10YR 3/5", return_triplets=FALSE, returnLAB=TRUE)
  expect_equal(names(test.3), c('L', 'A', 'B'))
  
  # test the LAB ---> sRGB is close
  test.4 <- grDevices::convertColor(test.3, from = 'Lab', to='sRGB')
  
  # sRGB (r)
  expect_equal(test.1[, 1], test.4[, 1], tolerance=0.1)
  # sRGB (g)
  expect_equal(test.1[, 2], test.4[, 2], tolerance=0.1)
  # sRGB (b)
  expect_equal(test.1[, 3], test.4[, 3], tolerance=0.1)
})


test_that("similar colors result in same, closest chip", {
  
  cols <- t(col2rgb(c('#5F5335', '#5F5236'))) / 255
  res <-  col2Munsell(cols)
  
  expect_equal(res$hue[1], res$hue[2])
  expect_equal(res$value[1], res$value[2])
  expect_equal(res$chroma[1], res$chroma[2])
})


test_that("munsell2spc wrapper method works as expected", {
  
  data(sp3)
  depths(sp3) <- id ~ top + bottom
  
  # inspect input data
  # horizons(sp3)[,c("hue","value","chroma")]
  
  # do color conversions to sRGB and LAB, join into horizon data
  expect_silent( {sp3 <- munsell2spc(sp3)})
  expect_true(inherits(sp3, 'SoilProfileCollection'))
  
  # # plot rgb "R" coordinate by horizon
  # plot(sp3, color = "rgb_R")
  #
  # # plot lab "A" coordinate by horizon
  # plot(sp3, color = "lab_A")
  
  # test returning profile+horizon ID data.frame with results
  expect_silent( { dftest <- munsell2spc(sp3, as.spc = FALSE) } )
  expect_true(inherits(dftest, 'data.frame'))
  
  # foo is not a column in horizons()
  expect_error( { err1 <- munsell2spc(sp3, hue = "foo") } )
  
  # chip is not a column in horizons
  expect_error( { d1 <- horizons(munsell2spc(sp3, .data = "chip")) } )
  
  # create chip as a combination of hue value/chroma
  sp3$chip <- with(horizons(sp3), sprintf("%s %s/%s", hue, value, chroma))
  
  # calculate from: column name, vector, data.frame
  expect_silent( { d1 <- horizons(munsell2spc(sp3, .data = "chip")) } )
  expect_silent( { d2 <- horizons(munsell2spc(sp3, .data = sp3$chip)) } )
  expect_silent( { d3 <- horizons(munsell2spc(sp3, .data = parseMunsell(sp3$chip, convertColors = FALSE))) } )
  expect_silent( { d4 <- horizons(munsell2spc(sp3, .data = data.frame(foo = sp3$chip))) } )
})

