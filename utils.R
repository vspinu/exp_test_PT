library(utils)
library(plyr)
library(ggplot2)

### Log Likelihoods
log_like <- function(theta,  sub,  estimate = TRUE,  V = v_U, ix = NULL, ...){
    ## V - vectorized value function of at least 3 arguments: matrix of values,
    ## matrix of theta parameters and index
    ## 
    ## theta - parameter matrix passed to V
    ##
    ## estimate - T if to operate on fit data. F to operate on predict data
    ##
    ## ix is a list. It stores additional preprocesed data depending on the
    ## model. Most commonly usually indexes.
    if( estimate ){
        D <- Df[sub, ]
        dV <- (V(Lf, theta = theta, ix = ix$Lf, ...) - V(Rf, theta = theta, ix = ix$Rf, ...))
    } else {
        D <- Dp[sub, ]
        dV <- (V(Lp, theta = theta, ix = ix$Lp, ...) - V(Rp, theta = theta, ix = ix$Rp, ...))
    }
    sum( log(c(pnorm(dV[D], sd = theta[["sig"]]), pnorm(-dV[!D], sd = theta[["sig"]]))) )
}

get_estim <- function(subjects, theta, ui, ci,  V, ix, ...){ ## theta = init pars, V = value function
    est <- cbind(t(theta), lfit = 0, lpred = 0, conver = 0, iter = 0)[-1, ]
    st <- system.time(
        for(sub in subjects){
            fit <- constrOptim(theta = theta, 
                               f = function(theta) -log_like(theta, sub, V = V, ix = ix, ...),
                               grad = NULL, ui = ui, ci = ci,
                               control = list(trace = F, maxit = 2000))
            lpredict <- log_like(fit$par, sub, estimate = FALSE, V = V, ix = ix, ...)
            est <- rbind(est, c(fit$par, fit$value, lpredict,
                                conver = as.logical(fit$convergence), iter = fit$counts[[1]]))
        })
    print(st)
    ##   fit <- maxNM(function(theta) log_like(theta, sub, V = V, ix = ix, ...),
    ##                start = theta, iterlim = 5000, finalHessian = F, 
    ##                constraints = list(ineqA = ui, ineqB = -ci))
    ##   lpredict <- log_like(fit$estimate, sub, estimate = FALSE, V = V, ix = ix, ...)
    ##   est <- rbind(est, c(fit$estimate, fit$maximum, lpredict,
    ##                       conver = as.logical(fit$code), iter = fit$iterations))
    est
}


### CONSTRAINTS EXPRESSION HANDLING

e <- function(i, n){ #return a vector of length n with one at i'th position and 0 elsewhere
  if(i<1 || i > n)
    return("error: i should be between 1 to n")
  return(c(rep(0, i-1), 1, rep(0, n-i)))
}

write_results <- function(res, model, sufix = "")
  write.table(round(res, 2), paste("estim/", model, sufix, ".csv", sep = ""),
              sep = "\t", row = F)

constr <- function(..., theta = get("theta", parent.frame())){
  Z <- diag(length(theta))
  E <- list()
  for( i in seq_along(theta) )
    E[[names(theta)[i]]] <- Z[i, ]
  
  expr <- c(...)
  ui <- list()
  ci <- list()
  for( i in seq_along(expr) ){
    if( !is.null(expr[[i]]) ){
      ineq <- substituteDirect(expr[[i]], E)
      sign <- 
        if ( ineq[[1]] == ">" || ineq[[1]] == ">=" ) 1L
        else if ( ineq[[1]] == "<" || ineq[[1]] == "<=" ) -1L
        else stop("expressions should be iniequalities, found ", ineq[[1]])

      ls <- sign*eval(ineq[[2]])
      rs <- sign*eval(ineq[[3]])

      if( length(rs) == 1 ){
        ui[[i]] <- ls
        ci[[i]] <- rs
      }else if ( length(rs) == 1 ){      
        ui[[i]] <- rs
        ci[[i]] <- ls
      }else stop(" at least one side of equation must evaluate to a constant")
    }
  }
  list(ui = do.call(rbind, ui), ci = unlist(ci))
}  


### OTHER

get_central_df <- function(tplike){
  predl <- ldply(tplike[-35, ], function(x) data.frame("mean$_{.1}$" = mean(x, trim = .1), 
                                                       "mean$_{.05}$" = mean(x, trim = .05), 
                                                       mean = mean(x), 
                                                       median = median(x), check.names = F))
  rownames(predl) <- predl$.id
  predl$.id <- NULL
  predl[order(predl$'mean$_{.1}$', decreasing = T),, drop = F]
}

get_max_df <- function(df, tr = treat){
  treat <- tr[as.numeric(row.names(df))]
  paraorder <- c(MinReg = 1, MaxMax = 1, MaxMin = 3, 
                 Alpha = 6, "$alpha$MM" = 6, 
                 "MnEU" = 5, "MxEU" = 5, DFT = 5, 
                 EV = 3, CEU = 8, EU = 4,
                 CEU_w = 5, PT_w = 5, 
                 SCEU = 5, SPT = 5)
  .local <- function(x){
    l <- structure(x == max(x), names = names(x))
    if(sum(l) > 1){
      ## order ties 
      nm <- names(x)[l]
      l <- names(x) == nm[which.min(paraorder[nm])]
    }
    structure(l, names = names(x))
  }
  ldf <- adply(df, 1,  .local, .expand = F)[-1]
  out <- rbind(All = colSums(ldf), 
               do.call(rbind, by(ldf, tr, colSums)))
  out[, order(out["All", ], decreasing = T)]
}

get_rank_df <- function(df){
  treat <- treat[as.numeric(row.names(df))]
  .local <- function(x) structure(rank(-x), names = names(x))
  ldf <- adply(df, 1, .local , .expand = F)[-1]
  out <- rbind(All = colMeans(ldf), 
               do.call(rbind, by(ldf, treat, colMeans)))
  out[, order(out["All", ])]
}

wilk <- function(tplike){
  wilk_df <- melt(tplike, id.vars = NULL, 
                  variable.name = "Model", value.name = "tplike")
  options(scipen = 0)
  wilk <- with(wilk_df, 
               pairwise.wilcox.test(tplike, Model, p.adjust.method = "none", 
                                    paired = TRUE))
  pp <- format.pval(wilk$p.value, 1, na.form = "-", eps = 1e-3)
  attributes(pp) <- attributes(wilk$p.value)
  pp
}

count <- function(tplike){
  nms <- names(tplike)
  counts <- outer(nms, nms, 
                  Vectorize(function(nm1, nm2){ 
                    sum(tplike[, nm1] <  tplike[, nm2])
                  }))
  dimnames(counts) <- list(nms, nms)
  counts
}

count_sig <- function(tplike){
  nms <- names(tplike)
  counts <- outer(nms, nms, 
                  Vectorize(function(nm1, nm2){ 
                    if(nm1 == nm2)
                      "-"
                    else{
                      wilk <- wilcox.test(tplike[, nm1], tplike[, nm2], paired = T)$p.value
                      wilk <- if(wilk < .001) "***" else if (wilk < .01) "**" else if (wilk < .05) "*" else ""
                      sprintf("$%d^{\\stackrel{\\textstyle%d\\hfill}{\\scriptscriptstyle%s\\hfill}}$", 
                              sum(tplike[, nm1] > tplike[, nm2]), 
                              sum(tplike[, nm1] < tplike[, nm2]), 
                              wilk)
                    }
                  }))
  counts[upper.tri(counts)] <- "-"
  dimnames(counts) <- list(nms, nms)
  counts[-1, -ncol(counts)]
}

count_table <- function(cplike, rplike){
    sapply(names(rplike),
           function(rn){
             sapply(names(cplike), function(cn){
               wilk <- wilcox.test(rplike[, rn], cplike[, cn], paired = T)$p.value
               wilk <- if(wilk < .01) "***" else if (wilk < .05) "**" else if (wilk < .1) "*" else ""
               sprintf("$%d^{\\stackrel{\\textstyle%d\\hfill}{\\scriptscriptstyle%s\\hfill}}$", 
                       sum(rplike[, rn] < cplike[, cn]), 
                       sum(rplike[, rn] > cplike[, cn]), 
                       wilk)
             })})}


print_tables <- function(df, file_name = NULL, caption, sub_ = F, ...){
    tb <- xtable(df, caption = caption, ..., label = paste("table", file_name, sep = ""))
    sanit <- function(x){ if(sub_) x <- gsub("_", "\\_", x, fixed = T)
        x <- gsub(">", "$>$", x, fixed = TRUE)
        x <- gsub("<", "$<$", x, fixed = TRUE)
    }

    dir.create("tables/", showWarnings = F)
    
    if(!is.null(file_name))
        print(tb, file = sprintf("tables/table_%s.tex", file_name),
              sanitize.text.function = sanit, floating.environment = "table")
    print(tb, ..., sanitize.text.function = sanit)
    df
}

sane_names <- function(x){
  x <- rename(x, c(ChoquetEU = "CEU", G.SMN = "MN-EU", G.SMX = "MX-EU", 
                   G.S.MM = "MnEU", G.S.MX = "MxEU", Min.Reg = "MinReg", 
                   Alpha = "$\\alpha$MM"), warn_missing = F)
  names(x) <- gsub("_no_u", "",
                   gsub("CEU_no_u", "CEV",
                        gsub("CEU_w_u|CEU_w", "SCEU",
                             gsub("PT_w", "SPT",
                                  gsub("CEU_w_u", "SCEU",
                                       gsub("CEU_w_no_u", "SCEV", names(x)))))))
  names(x) <- gsub("MP_max", "MxEU", gsub("MP_min", "MnEU", gsub("MP_alpha", "$\\alpha$MM", names(x), fixed = T)))
  names(x) <- gsub("_no_l", "$_{\\\\lambda = 1}$", names(x))
  names(x) <- gsub("_Prel1", "", names(x))
  names(x) <- gsub("_Prel2", "$_{2}$", names(x))
  names(x) <- gsub("_TK", "$_{TK}$", names(x))
  names(x) <- gsub("_GE2", "$_{GE2}$", names(x))
  names(x) <- gsub("_GE1", "$_{GE1}$", names(x))
  names(x) <- gsub("_GE", "$_{GE}$", names(x))
  names(x) <- gsub("_NA", "$_{NA}$", names(x))
  names(x) <- gsub("_s", "$_{\\\\pm}$", names(x))
  names(x) <- gsub("_tri", "", names(x))
  ## should be last
  names(x) <- gsub("_u", "$_{u}$", names(x))
  x
}

get_var_df <- function(var, res_ind){
  out <- sapply(res_ind, 
                function(el){
                  if(var %in% colnames(el)) el[, var]
                  else NULL
                }, simplify=F)
  as.data.frame(do.call(cbind, out))
}

P <- function(var, res_ind){
  treatments <- data_ind[["treatments"]]
  vals <- cbind(id = factor(1:length(treatments)), Treatment = treatments, abs(get_var_df(var, res_ind)))
  df <- melt(vals, id = c("id", "Treatment"), variable.name = "Model", 
             value.name = var)
  df$Treatment = factor(df$Treatment, levels = rev(levels(df$Treatment)))
  ## df$Treatment <-
  ##   reorder(df$Treatment), order = function(x) -mean(x))
  ggplot(df, aes_string(x = "Model", y = var)) +
    scale_fill_grey(start = .5, end = .95) + theme_bw() +
      scale_x_discrete(labels = prepare_labs(names(res_ind)))
}

prepare_labs <- function(names){
  labs <- label_parsed(
    , gsub("\\\\pm", "'+-'",
           gsub("\\}\\$", "]",
                gsub("\\$_\\{","[", names))))
  names(labs) <- names
  labs <- unlist(labs)
}

## p <- ggplot(mtcars, aes(factor(cyl), mpg))
## p + geom_boxplot(aes(fill = factor(am))) + coord_flip()
## savePlot("ggplot.png")

## (p <- qplot(1:5, 1, colour = letters[1:5]))
## (p <- qplot(1:5, 1, colour = letters[1:5]) + coord_flip() + guides(col = guide_legend( reverse = T)))
## savePlot("ggplot2.png")

## savePlot("ggplot3.png")



## +  guides(col = guide_legend( reverse = T))
