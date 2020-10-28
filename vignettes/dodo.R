dodo = function(sx=1, fixpars = TRUE, ICAR=TRUE, n_cores=20, re_order =c(2,2), 
    test=FALSE, sub_set=0) {
# remotes::install_github("kklot/dbtmb")
library(ktools)
library(data.table)
library(tidyverse)
library(magrittr)
set.seed(2020)

meta <- list() # track used metadata for post-processing
meta$seed <- 2020

# Read and prep data

pth = "data/pooled_AFS_lite.csv.bz2"
dhs = data.table::fread(pth)
dhs[, age := svy - yob]

# excluding old
dhs = dhs[age > 14 & age <= 50 & yob >=1950 & svy >=1990]
dhs = dhs[yob!=2004] # the one and only

# SSA only
ISO_SSA = name2iso(ktools:::.UN_SSA)
dt = dhs[sex==sx & ISO_A3 %in% ISO_SSA]

# remove fhs and BAIS IV under zero
dt = dt[!ver %in% c("BAISIV", "FHS")]

if (sub_set > 0)
    dt = dt[, .SD[sample(.N, ifelse(.N > sub_set,  sub_set, .N))], ISO_A3]

# Model indices

age_id = range2seq(dt[, age])
n_age_id = length(age_id)

meta$age <- age_id

yob_id = range2seq(dt[, yob])
n_yob_id = length(yob_id)

meta$yob <- yob_id

# Rs
yob_order = re_order[1]
age_order = re_order[2]
R_age = genR(n_age_id, order = age_order)
R_yob = genR(n_yob_id, order = yob_order)

# countries/ now remove those without data
if (sx==1)
    nbmat = readRDS('data/nbmat_SSA_men.rds')
else if (sx==2)
    nbmat = readRDS('data/nbmat_SSA_women.rds')
n_cc_id = nrow(nbmat)
cc_id = data.table(ISO_A3 = rownames(nbmat), cc_id = 1:n_cc_id)

# find singletons to constraint/different in each sex
svs = dt[, .(
    svy  = sort(unique(svy)), 
    svy2 = remove_consecutive(sort(unique(svy)))), ISO_A3]
svs = svs[, length(unique(svy2))==1, ISO_A3]
cc_id[svs, on="ISO_A3", singleton := i.V1]
cc_id[is.na(singleton), singleton := TRUE]
# cc_id[singleton==TRUE]
singlesvy <- as.numeric(cc_id$singleton)

meta$cc_id <- cc_id

# ICAR
# R_cc = INLA::inla.scale.model(R_cc, constr=list(A=matrix(1, 1, n_cc_id), e=0))
# R_cc = readRDS('data/R_cc_scaled_constr.rds') %>% as.matrix # 'm reading all matrix as not sparse

if (ICAR) {
    if (sx==1) {
        R_cc = readRDS('data/Q_men_ssa.rds') %>% as.matrix 
    } 
    else {
        R_cc = readRDS('data/Q_women_ssa.rds') %>% as.matrix 
    }
} 
else {
    R_cc = diag(n_cc_id)
}
R_cc_rank = qr(R_cc)$rank # not really needed

# Interaction with birth cohort
R_ccxyob      = kronecker(R_cc, R_yob)
R_ccxyob_rank = qr(R_ccxyob)$rank # not really needed

# Interaction with age
R_ccxage      = kronecker(R_cc, R_age)
R_ccxage_rank = qr(R_ccxage)$rank # not really needed

# interaction id
# - expand.grid run the first argument first which is not what we want
# - for the kronecker product, need to switch position of expand.grid

ccxyob_id = expand.grid(yob_id, 1:n_cc_id) %>% 
  set_colnames(c('yob', 'cc_id')) %>% 
  mutate(ccxyob_id = 1:n()) %>%
  data.table

ccxage_id = expand.grid(age_id, 1:n_cc_id) %>% 
  set_colnames(c('age', 'cc_id')) %>% 
  mutate(ccxage_id = 1:n()) %>%
  data.table

# add indices to data
dt[cc_id, on="ISO_A3", cc_id := i.cc_id]
dt[ccxyob_id, on=c('cc_id', 'yob'), ccxyob_id := i.ccxyob_id]
dt[ccxage_id, on=c('cc_id', 'age'), ccxage_id := i.ccxage_id]

# TMB metadata and data

data = list(
    # data
    afs           = dt[, afs], # right censor 
    afs_l         = dt[, afs - 0], #! interval
    afs_u         = dt[, afs + 1], 
    age           = dt[, age - min(age)],
    event         = dt[, event],
    yob           = dt[, yob - min(yob)], 
    cc_id         = dt[, cc_id-1],
    ccxyob_id     = dt[, ccxyob_id-1],
    ccxage_id     = dt[, ccxage_id-1],
    svw           = dt[, scaled_w],
    singlesvy     = singlesvy, 
    # priors
    sd_beta       = c(-2, 0.5), # mostly non positive intercept
    # RW order
    yob_order     = yob_order,
    age_order     = age_order,
    # RE precision
    sd_yob        = c(1e-2, 0.1), # precision not sd but specify in sd
    sd_age        = c(1e-2, 0.1),
    sd_cc         = c(1e-2, 0.1),
    sd_ccxyob     = c(1e-2, 0.1),
    sd_ccxage     = c(1e-2, 0.1),
    # Hyper Skewed LL
    palpha        = c(2.2, 0.25), # mean alpha
    p_a           = c(2.8, 0.35), # skewness * mean
    # RE structure
    R_age         = R_age,
    R_yob         = R_yob,
    R_cc          = R_cc,
    R_cc_rank     = R_cc_rank,
    R_ccxyob      = R_ccxyob,
    R_ccxyob_rank = R_ccxyob_rank,
    R_ccxage      = R_ccxage,
    R_ccxage_rank = R_ccxage_rank
)

init = list(
    intercept        = -2.0,
    log_alpha_vec    = rep(log(10), n_cc_id),
    a_vec_star       = rep(3, n_cc_id),

    yob_rw2          = rep(0, n_yob_id),
    age_rw2          = rep(0, n_age_id),
    cc_vec           = rep(0, n_cc_id),
    ccxyob           = rep(0, n_cc_id * n_yob_id),
    ccxage           = rep(0, n_cc_id * n_age_id),

    log_yob_rw2_e    = log(sd2prec(1e-2)),
    log_age_rw2_e    = log(sd2prec(1e-2)),
    log_cc_e         = log(sd2prec(1e-2)),
    log_ccxyob_e     = log(sd2prec(1e-2)),
    log_ccxage_e     = log(sd2prec(1e-2))
)

if (fixpars)
    fixpars = tmb_fixit(init, c(
        'log_ccxyob_e', 
        'log_ccxage_e', 
        'log_age_rw2_e', 
        'log_yob_rw2_e', 
        'log_cc_e'))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = c('intercept', 'age_rw2', 'yob_rw2', 'cc_vec', 'ccxyob', 'ccxage', 
        'log_alpha_vec', 'a_vec_star'),
    silent     = 0,
    DLL        = 'dbtmb', 
    map        = fixpars
)

if (test) return(0)

# Fit
library(TMB)
openmp(n_cores)
library(dbtmb)
config(tape.parallel=0, DLL="dbtmb")
obj = do.call('MakeADFunSafe', opts)
obj$env$inner.control$tol10 = 0
fit = nlminb(obj$par, obj$fn, obj$gr)
# Report
o = list(obj=obj, fit=fit, meta=meta)
class(o) <- 'dbtmb'
o
}