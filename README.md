
Code and data for [An experimental test of prospect theory for predicting choice under ambiguity (2014)](http://link.springer.com/article/10.1007/s11166-014-9185-0)

## Data

Original data in csv format is in `./data/` subdirectory. Sourcing `data_init.R`
will transform the original data and produce the binary file `./data/data.RData`
containing enviroments with matrix representation of the data suitable for the
analysis.


## Running Models

File `utils.R` contains various utility functions. File `models.R` contains the
functions to compute all the models in the paper and some other, more exotic,
models. Please see the code in `run_all.R` for how to run the models.

Sourcing `run_all.R` runs all the models and saves estimation results into human
readable tables in `./estim/tables/` folder. It will also save the results of
the estimation into three binary files: `./estim/res_ind.RData` (individual
level estimates for all models), `./estim/res_PT.RData` (individual level
estimates for all PT models), `./estim/res_CV.RData` (crossvalidation results).

## Reproducing Tables and Plots

Sub-directory `./paper/` contains two `.Rnw` files `tables.Rnw` and
`wappendix.Rnw`. Weaving these files with `knitr` will produce `tables.tex` with
all the tables used in the main paper and `wappendix.tex` which is the entire
web appendix. `tables.Rnw` will also save all the latex tables in subdirectory
`./paper/tables/`.

Actual results might be slightly different from the original paper because of
the inherent imprecision of the maximum likelihood iterative algorithms, changes
in R and different seeds for cross validation.
