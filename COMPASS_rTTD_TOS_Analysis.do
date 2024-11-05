***************************************************************************
*** Stata codes to illustrate analysis using TOS and rTTD as time-scale 
*** in analysis of COMPASS decedent cohort data                                  

*******************************************
** data in multiple-record format:
** person with no event has one record of censored outcome
** person with X event(s) has X record(s) of observed outcomes plus 
**  one record of censored outcome
** entry is the calendar date of enrollment to study
** dod is the calendar date of death
** start and end are calendar dates of start and end of each record
** records are further splitted according to pre- and post-exposure to PC
** records are further splitted whenever time-varying covariate values are updated

*******************************************
** specification of covariates
 
** time-costant covariates

global TCC age gender medifund ib2.cancertype i.education

** time-varying variables

global TVC FACT_G_PWB FACT_G_FWB 
 
***********************************************
** AG model 1: Time-on-study (TOS) as timescale

gen start2=start-entry
gen end2=end-entry

stset end2, failure(event) exit(t .) id(pid) enter(start2) 

* with no covariates adjustment
stcox PC,  cluster(pid)
* with TCC adjustment only
stcox PC $TCC,  cluster(pid)
* with TCC + TVC adjustment
stcox PC $TCC $TVC,  cluster(pid)

********************************************************
** AG model 2: reverse time-to-death (rTTD) as timescale

stset, clear

gen deathtime=dod-entry
gen surv_start=deathtime-start2
gen surv_end=deathtime-end2

sum deathtime
local d=r(max)

gen death_start=`d'-surv_start
gen death_end=`d'-surv_end

stset death_end, failure(event) exit(t .) id(pid) enter(death_start) 

* with no covariates adjustment
stcox PC,  cluster(pid)
* with TCC adjustment only
stcox PC $TCC,  cluster(pid)
* with TCC + TVC adjustment
stcox PC $TCC $TVC,  cluster(pid)

