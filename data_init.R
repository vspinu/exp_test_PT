## The original maximum likelihood code operates on old version of the data
## represented as matrices. This builds this matrices from the new publicly
## available long form of the data.

source("./utils.R")

idQ <- 1:162

make_data_env <- function(qPredict = NULL){
    ## qPredict is a vector of questions used for prediction. All others are
    ## used for fitting. If NULL, the original HLM set is used.
    
    env <- new.env()
    with(env, {
        
        ## Reconstruct old version of the data:
        hlm <- read.csv("./data/HLM.csv")
        ## We were not aware that there were 15 indifferences 
        hlm <- within(hlm, Choice <- Choice != "L")
        library(reshape2)
        hlmWide <- dcast(hlm, Subject + Treatment ~ Question, value.var = "Choice")
        choiceData <- as.matrix(hlmWide[, -c(1, 2)])
        Lnames <- c("Lp", "Ly", "Lb")
        Rnames <- c("Rp", "Ry", "Rb")
        questions <- unique(hlm[, c("Question", "Test", Lnames, Rnames)])
        ## normalize
        for(nm in c(Lnames, Rnames))
            questions[[nm]] <- questions[[nm]]/100
        middle_val <- .1
        if(is.null(qPredict))
            qPredict <- which(questions$Test)
        qFit <- setdiff(idQ, qPredict)
        questions <- as.matrix(questions[, -c(1, 2)])
        treatments <- as.factor(hlmWide$Treatment)

        ## Data in matrix form:

        ## fit (train) data
        Df <- choiceData[, qFit]
        Lf <- questions[qFit, Lnames]
        Rf <- questions[qFit, Rnames]

        ## prediction (test) data
        Dp <- choiceData[, qPredict]
        Lp <- questions[qPredict, Lnames]
        Rp <- questions[qPredict, Rnames]
        
    })
    env
}

## HLM set
data_ind <- make_data_env()

## REPRESENTATIVE AGENT 
data_rep <- new.env(parent = data_ind)
with(data_rep, {
    ## stack individuals in one long row
    Df <- matrix(t(Df), nrow = 1)
    Dp <- matrix(t(Dp), nrow = 1)
    Lf <- Lf[ rep(1:nrow(Lf), 48), ]
    Rf <- Rf[ rep(1:nrow(Rf), 48), ]
    Lp <- Lp[ rep(1:nrow(Lp), 48), ]
    Rp <- Rp[ rep(1:nrow(Rp), 48), ]
})

## CROSS VALIDATION
CV_cuts <- cut(idQ, 10, labels = FALSE)
CV_groups <- sample(CV_cuts, length(CV_cuts))
data_CV <- tapply(idQ, CV_groups, make_data_env)

## unqualified _w models use Prelec 1
treat <- factor(c(rep("Tr 1", 15), rep("Tr 2", 17), rep("Tr 3", 16)))
HLM_pll <- read.table("./data/predictedLogLikelihoods.csv", T, sep = ",")
HLM_pll <- sane_names(HLM_pll[-grep("PT|Prospect", names(HLM_pll))])

## used model names
nms_main <- c("PT_w_u",  "PT_w_Prel2_no_u", 
              "PT_w_no_u", "PT_w_TK_no_u", "PT_w_NA_no_u", "PT_no_u", 
              "PT_w_s_no_u", "PT_w_GE_no_u", "CEU_w_u")
nms_other <- c("EU", "CEU_w_no_u", "PT_w_no_l_no_u")
nms_all <- c(nms_main, nms_other)

save(HLM_pll, data_ind, data_rep, data_CV, nms_all, nms_main, nms_other, file = "./data/data.RData")
