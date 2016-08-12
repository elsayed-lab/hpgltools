library(testthat)
library(hpgltools)

context("Can I create a organismdbi from a problematic species (tcruzi-nonesmer)?")

difficult <- make_organismdbi(id="tcruzi_nonesmer")
