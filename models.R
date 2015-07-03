
#################### EV #####################################
v_EV <- function(x, theta, ...){ # Value function Identity
  p <- c(theta[["p1"]], theta[["p2"]], 1 - theta[["p1"]] - theta[["p2"]])
  x %*% p
}

EV <- function(subjects = 1:3, sig = 1, p1 = 0.5, p2 = .02){
  theta <- c(sig = sig, p1 = p1, p2 = p2)
  con <- constr(expression(sig > .01, p1 > 0, p2 > 0, p1 + p2 < 1))

  get_estim(subjects, theta, con$ui, con$ci, v_EV)
}


#################### EU #####################################
v_EU <- function(x, theta, use_u = T, ...){ 
  if(use_u) x[x == middle_val] <- theta[["u"]]
  p <- c(theta[["p1"]], theta[["p2"]], 1 - theta[["p1"]] - theta[["p2"]])
  x %*% p
}


EU <- function(subjects = 1:3, sig = 1, p1 = 0.5, p2 = .2, u = .7){
  theta <- c(sig = sig, p1 = p1, p2 = p2, u = u)
  con <- constr(expression(sig > .01, p1 > 0, p2 > 0, u > -.1, p1 + p2 < 1))
  
  get_estim(subjects, theta, con$ui, con$ci, v_EU)
}


#################### CEU && PT utils #########################

## orders (events in decreasing order). For example "231" means 2nd color is the
## best, 3rd is second best,  1st is worst outcome
w_pos <- list("123" = c("w1", "d12_1", "dA_12"),
              "132" = c("w1", "d13_1", "dA_13"),
              "312" = c("w3", "d13_3", "dA_13"),
              "321" = c("w3", "d23_3", "dA_23"),
              "231" = c("w2", "d23_2", "dA_23"),
              "213" = c("w2", "d12_2", "dA_12"))

w_neg <- list("123" = c("n_p3", "n_d23_3", "n_dA_23"),
              "132" = c("n_p2", "n_d23_2", "n_dA_23"),
              "312" = c("n_p2", "n_d12_2", "n_dA_12"),
              "321" = c("n_p1", "n_d12_1", "n_dA_12"),
              "231" = c("n_p1", "n_d13_1", "n_dA_13"),
              "213" = c("n_p3", "n_d13_3", "n_dA_13"))


d_from_w <- function(w)
  c(w[c("w1", "w2", "w3")],
    structure(
      c(w[c("w12", "w12", "w13", "w13", "w23", "w23")], 1    , 1    , 1    ) -
      c(w[c("w1" , "w2" , "w1" , "w3" , "w2" , "w3"   , "w12", "w13", "w23")]), names =
      c("d12_1","d12_2","d13_1","d13_3","d23_2","d23_3","dA_12","dA_13","dA_23")))

get_ix_PT <- function(CPT = FALSE){ 
  .local <- function(x){
    out <- t(apply(x, 1, function(r){
      nN <- sum(r == 0)
      order <- order(r, decreasing = T)
      inv_order <- order(order)
      order_name <- paste(order, collapse = "")
      wp <- w_pos[[order_name]]
      wn <- w_neg[[order_name]]
      w <-
        if( CPT || nN == 0L ) wp
        else c(wp[1:(3-nN)], rev(wn[1:nN]))
      w[inv_order] # put back in proper order
    }))
    dim <- dim(out)
    out <- as.factor(out)
    dim(out) <- dim
    out
  }

  list(Rp = .local(Rp), 
       Lp = .local(Lp),
       Rf = .local(Rf),
       Lf = .local(Lf))
}

get_ix_CEU <- function()
  get_ix_PT(CPT = TRUE)


### WEIGHTING FUNCTIONS
## THETA - a numeric vectors, DEFAULTS ARE STARTINTG VALUES!!
## CONSTR - an expression of constraints
## fixme: should be vectorised!!!
w_pow <- function(p, theta = c(a = 1), constr = expression(a > .01))
  p^theta[["a"]]

w_TK <- function(p,  theta = c(a = 1), constr = expression(a > 0.2, a < 4)){ 
  a <- theta[["a"]]
  p^a/(p^a + (1 - p)^a)^(1/a)
}

w_NA <- function(p,  theta = c(a = .1, b = .1), constr = expression(a > -.5, a < 1, b >  -.5, b < 1)){
  ## neo additive
  a <- theta[["a"]]
  pmin(pmax(a + p*(1-theta[["b"]] -a), 0), 1)
}

w_NA_b <- function(p,  theta = c(b = .1), constr = expression(b < 1), d = 0){
  pmin(pmax(d + p*(1-theta[["b"]] - d), 0), 1)
}

w_NA_a <- function(p,  theta = c(a = .1), constr = expression(a < 1), d = 0){
  pmin(pmax(theta[["a"]] + p*(1 - d - theta[["a"]]), 0), 1)
}

w_NA1 <- function(p,  theta = c(a = .1), constr = expression(a < 1)){
  ## neo additive
  a <- theta[["a"]]
  pmin(pmax(a + p*(1-2*a), 0), 1)
}

w_Beta <- function(p, theta = c(a = .5, b = 1),
                   constr = expression(a > .1, a < 10, b > .1, b<10)){
  pbeta(p, theta[["a"]], theta[["b"]])
}

w_Beta_a <- function(p, theta = c(a = .5),
                    constr = expression(a > .1, a < 10), d = 1.5){
  pbeta(p, theta[["a"]], d)
}

w_Beta_b <- function(p, theta = c(b = .5),
                     constr = expression(b > .1, b < 7), d = 1){
  pbeta(p, d, theta[["b"]])
}

w_Prelec <- function(p, theta = c(a = 1, b = 1), constr = expression(a > 0.01, b > .01)){
  exp(-(-log(p))^theta[["a"]])^theta[["b"]]
}

w_Prelec_b <- function(p, theta = c(b = 1), constr = expression(b > .01, b < 3), d = 1.18){
  exp(-(-log(p))^d)^theta[["b"]]
}

w_Prelec_a <- function(p, theta = c(a = 1), constr = expression(a > 0), d = 1){
  exp(-(-log(p))^theta[["a"]])^d
}

w_Prelec_a1 <- function(p, theta = c(a1 = 1), constr = expression(a1 > 0), d = 1){
  exp(-(-log(p))^theta[["a1"]])^d
}

w_GE <- function(p,  theta = c(a = 1, b = 1), constr = expression(a > 0.2, a < 3, b > 0)){ 
  a <- theta[["a"]]
  b <- theta[["b"]]
  b*p^a/(b*p^a + (1 - p)^a)
}

w_GE_a <- function(p,  theta = c(a = 1), constr = expression(a > 0.2, a < 4), d = 1){ 
    a <- theta[["a"]]
    d*p^a/(d*p^a + (1 - p)^a)
  }


#################### CEU #####################################

v_CU <- function(x, theta, ix, w = NULL, kadd = F, use_u = T, use_l = F, ...){
  if(!is.null(w)){ # sophistication
    p <- 
      if( kadd ) theta[c("p1", "p2", "p3")]
      else c(theta[c("p1", "p2")], p3 = 1 - sum(theta[c("p1", "p2")]))
    theta <- c(theta, w(c(w1 = p[[1]], w2 = p[[2]], w3 = p[[3]],
                          w12 = p[[1]] + p[[2]], w13 = p[[1]] + p[[3]], w23 = p[[2]] + p[[3]]),
                        theta))
  }
  if( use_u ) x[x == middle_val] <- theta[["u"]]
  if( use_l )   x[x == -.1] <- -theta[["l"]]
  d <- d_from_w(theta)[levels(ix)] # reorder 

  rowSums( x * array(d[ix], dim = dim(x)) ) # rowsum is faster but trickier
}

CEU <- function(subjects = 1:3, use_u = T, use_l = F, 
                sig = 1, w1 = 0.2, w2 = .2, w3 = .2,  w12 = .4, w13 = .4, w23 = .5, u = .7, l = .1){
  ix <- get_ix_CEU()
  theta <- c(sig = sig, w1 = w1, w2 = w2, w3 = w3, w12 = w12, w13 = w13, w23 = w23,
             if(use_u) c(u = u), if(use_l) c(l = l))
  ## theta = c(sig = 2, w1 = 0.2, w2 = 0.5, w3 = 0.3, w12 = 0.7, w13 = 0.5, w23 = 0.8, u = 0.4)
  con <- constr(expression(sig > .01, w1 > 0, w2 > 0, w3 > 0,
                           w12 > 0, w13 > 0, w23 > 0, 
                           w12 < 1, w13 < 1, w23 < 1, 
                           w12 - w1 > 0, w23 - w2 > 0, w13 - w1 > 0,
                           w12 - w2 > 0, w23 - w3 > 0, w13 - w3 > 0),
                if(use_u) expression(u < 1, u > -.1),
                if(use_l) expression(l > -.1, l < 4))
  
  get_estim(subjects, theta, V = v_CU, con$ui, con$ci, ix = ix, use_u = use_u)
}

CEU_w <- function(subjects = 1:3, w = w_Prelec_a, kadd = FALSE, 
                  use_u = T, use_l = F, sig = 1, p1 = 0.2, p2 = .2, u = .7, l = .1){
  ix <- get_ix_CEU()
  if( kadd ) w <- function(x, ...) x
  theta <- c(sig = sig, p1 = p1, p2 = p2, if(use_u) c(u = u), if(use_l) c(l = l), 
             if(kadd) c(p3 = .5), eval(formals(w)[["theta"]]))
  con <- constr(expression(sig > .01, p1 > 0, p2 > 0, p1 + p2 < 1),
                if(kadd) expression(p3 > 0, p1 + p3 < 1, p2 + p3 < 1),
                if(use_u) expression(u < 1, u > -.1),
                if(use_l) expression(l > -.1, l < 4), 
                eval(formals(w)[["constr"]]))
    
  get_estim(subjects, theta, V = v_CU, con$ui, con$ci,
            w = w, ix = ix, kadd = kadd, use_u = use_u, use_l = use_l)
}



#################### PT #####################################

v_PT <- function(x, theta, ix, w = NULL, use_u = TRUE, use_l = T, kadd = F, wneg = NULL, ...){
  if( !is.null(w) ){ # sophistication or kadd
    p <-
      if( kadd ) theta[c("p1", "p2", "p3")]
      else c(theta[c("p1", "p2")], p3 = 1 - sum(theta[c("p1", "p2")]))
    P <- c(w1 = p[[1]], w2 = p[[2]], w3 = p[[3]], w12 = p[[1]] + p[[2]],
           w13 = p[[1]] + p[[3]], w23 = p[[2]] + p[[3]])
    d <- d_from_w(w(P, theta))
    if( !is.null( wneg) ){
      d_n <- d_from_w(wneg(P, theta))
      names(d) <- paste("n_", names(d), sep = "")
      d <- c(d, d_n)
    }
  }else{
    ## weights are part of theta
    d <- d_from_w(theta)
  }
    
  ## theta[["w3"]] <- 1-theta[["w1"]]-theta[["w2"]]
  if( use_u ) x[x == middle_val] <- theta[["u"]]
  if( use_l ) x[x == -.1] <- -theta[["l"]]
  d <- d[levels(ix)] # reorder
  
  rowSums( x * array(d[ix], dim = dim(x)))
}

PT_w <- function(subjects = 1:3, use_u = T, use_l = T, w = w_Prelec_a, wneg = NULL, kadd = F, 
                 sig = .1, p1 = .2, p2 = .5, u = .5, l = 1){
  ix <- get_ix_PT()
  if( kadd ) w <- function(x, ...) x
  theta <- c(sig = sig, p1 = p1, p2 = p2, if (kadd) c(p3 = .5),
             if(use_u) c(u = u), if(use_l) c(l = l),
             if(!kadd) c(eval(formals(w)[[2]]), eval(formals(wneg)[[2]])))
  con <- constr(expression(sig > 0.01, p1 > 0, p2 > 0, p1 + p2 < 1),
                if(kadd) expression(p3 > 0, p1 + p3 < 1, p2 + p3 < 1, p1 + p2 + p3 > .4), 
                if(use_u) expression(u < 1, u > -.1),
                if(use_l) expression(l > 0, l < 4), 
                eval(formals(w)[["constr"]]), eval(formals(wneg)[["constr"]]))
  get_estim(subjects, theta, V = v_PT, con$ui, con$ci, ix = ix, w = w, wneg = wneg, 
            use_u = use_u, use_l = use_l, kadd = kadd)
}
  
  
PT <- function(subjects = 1:3, use_u = TRUE,
               sig = 1, w1 = 0.2, w2 = .2, w3 = .2,
               w12 = .6, w13 = .6, w23 = .6, u = .7, l = 1){
  ix <- get_ix_PT()
  theta <- c(sig = sig, w1 = w1, w2 = w2, w3 = w3,
             w12 = w12, w13 = w13, w23 = w23, if(use_u) c(u = u), l = l )
  con <- constr(expression(sig > .01, w1 > 0, w2 > 0,  l > 0, l < 4, w3 > 0,
                           w12 > 0, w13 > 0, w23 > 0, 
                           w12 < 1, w13 < 1, w23 < 1, 
                           w12 - w1 > 0, w23 - w2 > 0, w13 - w1 > 0,
                           w12 - w2 > 0, w23 - w3 > 0, w13 - w3 > 0), 
                if(use_u) expression(u < 1, u > -.1))
  
  get_estim(subjects, theta, V = v_PT, con$ui, con$ci, ix = ix, use_u = use_u)
}



###################### MULTIPLE PRIORS ############################

v_vectEU <- function(x, theta, use_u = T, ...){
  ## theta[["p1"]] and theta[["p2"]] are vectors
  ## return a list of computed EU (one per element of p1); easy to pass to pmin, pmax
  eu_vals <- list()
  u <- theta[["u"]]
  for (i in seq_along(theta[["p1"]]))
    eu_vals[[i]] <- v_EU(x, c(u = u, p1 = theta[["p1"]][[i]], p2 = theta[["p2"]][[i]]),
                         use_u = use_u)
  eu_vals
}

v_MP <- function(x, theta, mtype, set_type, use_u = T, ...){ ## simple
  new_theta <-
    switch(set_type[[1]], 
           hex = list(
             u = if(use_u) theta[["u"]],
             p1 = c(theta[["B1"]], 1 - theta[["B2"]] - theta[["b3"]], theta[c("b1", "b1")],
               1 - theta[["B3"]] - theta[["b2"]], theta[["B1"]]), 
             p2 = c(1 - theta[["B1"]] - theta[["b3"]], theta[c("B2", "B2")],
               1 - theta[["B3"]] - theta[["b1"]], theta[c("b2", "b2")])),
           delta = list(
             ## naming is syclical, 
             u = if(use_u) theta[["u"]],
             p1 = c(theta[["B1"]], theta[["b1"]], 1- theta[["b2"]] -theta[["B3"]]),
             p2 = c(1 - theta[["B1"]] - theta[["b3"]], theta[["B2"]], theta[["b2"]])),
           tri =  list(
             u = if(use_u) theta[["u"]],
             p1 = c(1 - theta[["b2"]] - theta[["b3"]], theta[c("b1", "b1")]), 
             p2 = c(theta[["b2"]], 1 - theta[["b1"]] - theta[["b3"]], theta[["b2"]])), 
           stop("unknown set_type: ", set_type))

  switch(mtype[[1]],
         alpha = {
           min <- do.call(pmin, v_vectEU(x, new_theta, use_u = use_u))
           max <- do.call(pmax, v_vectEU(x, new_theta, use_u = use_u))
           min*theta[["alpha"]] +  max*(1 - theta[["alpha"]])
         },
         min = do.call(pmin, v_vectEU(x, new_theta, use_u = use_u)),
         max = do.call(pmax, v_vectEU(x, new_theta, use_u = use_u)), 
         stop("uncnown mtype: ", mtype))
}


MP <- function(subjects = 1:3, set_type = c("tri", "hex", "delta"),
                   mtype = c("min", "max", "alpha"), use_u = T, sig = 1, u = .7,
                   b1 = .2, b2 = .2, b3 = .2, B1 = .5, B2 = .5, B3 = .5, alpha = .5){
  ## bN, BN stand for lower/upper *b*oundaries
  alpha_restr <- if(mtype[[1]] == "alpha") expression(alpha < 1, alpha > 0)
  alpha_init <- if (mtype[[1]] == "alpha") c(alpha = alpha)
  switch(set_type[[1]],
         tri = {
           theta <- c(sig = sig, b1 = b1, b2 = b2, b3 = b3, if(use_u) c(u = u), alpha_init)
           con <- constr(expression(sig > .01, b1 > 0, b2 > 0, b3 > 0, 
                                    b1 + b2 + b3 < 1),
                         if(use_u) expression(u > -.1, u < 1), 
                         alpha_restr)
         },
         delta = {
           theta <- c(sig = sig, b1 = b1, B1 = B1, b2 = b2,
                      B2 = B2, b3 = b3, B3 = B3, if(use_u) c(u = u), alpha_init)
           con <- constr(expression(sig > .01, b1 > 0, b2 > 0, b3 > 0, 
                                    B1 > 0, B2 > 0, B3 > 0, b3 + B1 < 1, b1 + B2 < 1,
                                    b2 + B3 < 1, B1 - b1 >  0, B1 + b2 + B3 > 1),
                         if(use_u) expression(u > -.1, u < 1), 
                         alpha_restr)
         }, 
         hex = {
           theta <- c(sig = sig, b1 = b1, B1 = B1, b2 = b2,
                      B2 = B2, b3 = b3, B3 = B3, if(use_u) c(u = u), alpha_init) 
           con <- constr(expression(sig > .01, b1 > 0, b2 > 0, b3 > 0, 
                                    B1 > 0, B2 > 0, B3 > 0, B1 < 1, B2 < 1, B3 < 1,
                                    b1 + b2 + b3 < 1, B1 + B2 + B3 > 1),
                         if(use_u) expression(u > -.1, u < 1), 
                         alpha_restr)
         })
  
  get_estim(subjects, theta, con$ui, con$ci, V = v_MP,
            mtype = mtype, set_type = set_type, use_u = use_u)
}

## MP_tri(1:3, "min")
## MP_tri(1:3, "max")
## MP(mtype = "max", set_type= "tri")
