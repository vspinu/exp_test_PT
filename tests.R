library(testthat)
source("../utils.R")
source("../models.R")
source("../data_init.R")
attach(data_ind)

test_that("v_PT is computed correctly",{
  w_Prelec <- function(p, theta = c(a = 1, b = 1), constr = expression(a > 0.01, b > .01)){
    exp(-(-log(p))^theta[["a"]])^theta[["b"]]
  }
  
  ix <- get_ix_PT()
  expect_equal(
    v_PT(Rf[1:10, ], theta = c(p1 = .1, p2 = .3, u = .2, l = .5, a = 1.1), ix = ix$Rf[1:10, ], w = w_Prelec_a),
    c(-0.0600373186740005, -0.442704831487583, -0.377224638901965, 
      -0.152277753868822, 0.104804400454769, -0.294684082047867, 0.143518866269238, 
      -0.442704831487583, 0.00742259856540682, -0.217757946454441), info = "(With utility)")
  expect_equal(
    v_PT(Lf[1:10, ], theta = c(p1 = .1, p2 = .9, u = .5, l = 2.1, a = .9), ix = ix$Lf[1:10, ], w = w_Prelec_a),
    c(-2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1), info = "(With utility 2)")
  expect_equal(
    v_PT(Rf[sub, ], theta = c(p1 = .1, p2 = .3, u = .2, l = .5, a = 1.1), ix = ix$Rf[sub, ], w = w_Prelec_a, use_u = F), 
    c(-0.0600373186740005, -0.450889855560786, -0.377224638901965, 
      -0.184413023159271, 0.104804400454769, -0.3240149274696, 0.0515875996593468, 
      -0.450889855560786, -0.0650663440867942, -0.258078239818092), info = "(Linear utility.)")
  expect_equal(
    v_PT(Rf[sub, ], theta = c(w1 = .1, w2 = .3, w12 = .5, w3 = .2, w23 = .4, w13 = .6, u = .14, l = .5, a = .7), ix = ix$Rf[sub, ], use_u = T), 
    c(-0.05, -0.436, -0.35, -0.094, 0.25, -0.308, -0.244, -0.436, 
      -0.116, -0.18), info =  "(General PT, with utility)")
  expect_equal(
    v_CU(Rf[sub, ], theta = c(w1 = .1, w2 = .3, w12 = .5, w3 = .2, w23 = .4, w13 = .6, u = .14, l = .5, a = .7), ix = ix$Rf[sub, ], use_u = T), 
    c(0.23, -0.076, 0.01, 0.106, 0.45, -0.028, -0.00399999999999999,  -0.076, 0.044, 0.02),
    info = "(General CEU, with utilty)")
  
})
  
  
  
