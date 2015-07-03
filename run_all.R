## run from this folder
## binary sore is in ./save subdir
## estimation output is in ./estim

source("./utils.R")
source("./models.R")
## source("./data_init.R")
load("./data/data.RData")
attach(data_ind)
## load("./data_init.RData")

run_model <- function(model, save = FALSE, sufix = "",  ...){
    cat("\nModel:  ", model, sufix, "\n", sep="")
    out <- do.call(model, list( 1:nrow(Df), ...))
    if( save ) write_results(out, model, sufix)
    as.data.frame(out)
}

run_PT <- function(){
    res_PT <- list()

    res_PT[["EU"]] <- run_model("EU")
    res_PT[["CEU_w_u"]] <- run_model("CEU_w", w = w_Prelec_a,
                                     use_u = T, use_l = F)
    res_PT[["CEU_w_no_u"]] <- run_model("CEU_w", w = w_Prelec_a,
                                        use_u = F, use_l = F)

    ## unqualified weighting function is Prelec1

    res_PT[["PT_w_u"]] <- run_model("PT_w", w = w_Prelec_a)
    res_PT[["PT_w_u_TK"]] <- run_model("PT_w", w = w_TK)
    res_PT[["PT_w_u_NA"]] <- run_model("PT_w", w = w_NA)
    res_PT[["PT_w_u_Prel2"]] <- run_model("PT_w", w = w_Prelec)

    res_PT[["PT_w_no_u"]] <- run_model("PT_w", w = w_Prelec_a, use_u = F)
    res_PT[["PT_w_no_l_no_u"]] <- run_model("PT_w", w = w_Prelec_a, use_l = F, use_u = F)
    res_PT[["PT_w_no_l_s_no_u"]] <- run_model("PT_w", w = w_Prelec_a, wneg = w_Prelec_a1, use_l = F, use_u = F)
    res_PT[["PT_w_s_no_u"]] <- run_model("PT_w", w = w_Prelec_a, wneg = w_Prelec_a1, use_u = F)

    res_PT[["PT_w_TK_no_u"]] <- run_model("PT_w", w = w_TK, use_u = F)
    res_PT[["PT_w_NA_no_u"]] <- run_model("PT_w", w = w_NA, use_u = F)
    res_PT[["PT_w_NA_a_no_u"]] <- run_model("PT_w", w = w_NA_a, use_u = F)
    res_PT[["PT_w_NA_b_no_u"]] <- run_model("PT_w", w = w_NA_b, use_u = F)
    res_PT[["PT_w_GE_no_u"]] <- run_model("PT_w", w = w_GE, use_u = F)
    res_PT[["PT_w_Prel2_no_u"]] <- run_model("PT_w", w = w_Prelec, use_u = F)

    wt <- w_Prelec_a
    formals(wt)[["d"]] <- .983 ## best "common-b-model"
    res_PT[["PT_w_B_no_u"]] <- run_model("PT_w", w = wt, use_u = F)

    ## general PT
    res_PT[["PT_u"]] <- run_model("PT")
    res_PT[["PT_no_u"]] <- run_model("PT", use_u = F)
    res_PT
}

run_all_models <- function(save = FALSE, sufix = ""){
    res <- list() 
    st <- system.time({
        for( M in c("EV", "EU") )
            Respr[[M]] <- run_model(M, save, sufix = sufix)

        for( M in c("CEU", "CEU_w", "PT", "PT_w") )
            res[[paste(M, "_u", sep = "")]] <-
                run_model(M, save, sufix = paste(sufix, "_u", sep = ""))


        for( M in c("PT", "PT_w", "CEU", "CEU_w"))
            res[[paste(M, "_no_u", sep = "")]] <-
                run_model(M, save, sufix = paste(sufix, "_no_u", sep = ""), use_u = FALSE)

        ## for( M in c("CEU_w", "PT_w") )
        ##   res[[paste(M, "_kadd", sep = "")]] <-
        ##     run_model(M, save, sufix = paste(sufix, "_kadd", sep = ""), kadd = TRUE)

        ## mm_refinement <- c("tri", "hex", "delta")
        mm_types <- c("tri")
        for( mtype in c("min", "max", "alpha"))
            for( stype in mm_types){
                suf <- sprintf("_%s_%s%s", mtype, stype, sufix)
                mname <- paste("MP",  mtype, stype, sep="_") 
                res[[mname]] <- run_model("MP", save, suf,
                                          mtype = mtype, set_type = stype)
            }

        ## for( mtype in c("min", "max", "alpha"))
        ##   for( stype in c("tri", "hex", "delta" )){
        ##     suf <- sprintf("_%s_%s%s_no_u", mtype, stype, sufix)
        ##     mname <- paste("MP",  mtype, stype, "no_u", sep="_")
        ##     res[[mname]] <- run_model("MP", save, suf,
        ##                               mtype = mtype, set_type = stype, use_u = FALSE)
        ##   }
        
    })
    cat("\nAll models time:\n")
    print(st)
    res
}

dir.create("./estim/tables", showWarnings = F, recursive = T)

## ALL PT MODELS
res_PT <- run_PT()
save(res_PT, file = "./estim/res_PT.RData")
for(nm in names(res_PT))
    write.table(round(res_PT[[nm]], 4), paste("./estim/tables/", nm, ".txt", sep = ""), sep = "\t")

## ALL MAIN MODELS
res_ind <- run_all_models(sufix="")
save(res_ind, file = "./estim/res_ind.RData")
for(nm in names(res_ind))
    write.table(round(res_ind[[nm]], 4), paste("./estim/tables/", nm, ".txt", sep = ""), sep = "\t")

### CROSS VALIDATION
res_CV <- list()
for (i in seq_along(data_CV)) {
    cat("#####  Running batch", i, "######\n")
    if("data_ind" %in% search()) detach("data_ind")
    attach(data_CV[[i]], name = "data_ind")
    res_CV[[i]] <- run_all_models(sufix = "_cv")
}
save(res_CV, data_CV, file = "./estim/res_CV.RData")
