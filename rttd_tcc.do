*** rttd_tcc.do
*** program rttd_tcc
*** simulate time-constant confounder only 

 
cd "C:\workig"
   
capture log close
log using rttd_tcc,text replace

version 15.1 
clear
set more off
set seed 123

** program rttd

capture program drop rttd_tcc
program define rttd_tcc,rclass
 syntax [, obs(real 600) r(real 0.5) hr(real 0.5) corr(real 0.5)]
 
drop _all
set obs `obs'
gen id=_n

matrix C=(1,`corr',`corr' \ `corr',1,.5 \ `corr',.5,1)
drawnorm x1 x2 x3,cov(C)

gen hazDeath=exp(-18+x1)
gen dod=ceil((-ln(runiform())/hazDeath)^(1/3))

sum dod,det
global TTDmax=ceil(r(max))

expand $TTDmax
sort id
qui by id: gen tosStart=_n-1
qui by id: gen tosEnd=_n
drop if tosEnd>dod
gen ttd=dod-tosStart
gen rttd=$TTDmax-ttd

* confoudner zt correlated with ttd

sum ttd
local mu_ttd=r(mean)
local sd_ttd=r(sd)
gen z_ttd=(ttd-`mu_ttd')/`sd_ttd'
gen zt=`r'*z_ttd+sqrt(1-`r'^2)*rnormal(0,1)

gen trtLogOdds=-8-2*zt+x2
gen trtRisk=exp(trtLogOdds)/(1+exp(trtLogOdds))
gen trt=rbinomial(1,trtRisk)
qui forvalues i=2(1)$TTDmax{
replace trt=1 if trt[_n-1]==1 & tosEnd==`i'
}
sort id trt tosStart
qui by id trt: gen n=_n
sum tosStart if trt==1 & n==1,det

gen eventLogOdds=-6+ln(`hr')*trt-0*zt+x3
gen eventRisk=exp(eventLogOdds)/(1+exp(eventLogOdds) )
gen event=rbinomial(1,eventRisk)
tab event

* tos
stset tosEnd,failure(event) exit(time .) id(id) enter(tosStart)

stcox trt, cluster(id) nolog
return scalar tos_crude_est=_b[trt]
return scalar tos_crude_CP=cond((_b[trt]-1.96*_se[trt])<ln(`hr') & (_b[trt]+1.96*_se[trt])>ln(`hr'),1,0)

stcox trt x3, cluster(id) nolog
return scalar tos_x_est=_b[trt]
return scalar tos_x_CP=cond((_b[trt]-1.96*_se[trt])<ln(`hr') & (_b[trt]+1.96*_se[trt])>ln(`hr'),1,0)

* rTTD
gen rttdEnd=rttd
gen rttdStart=rttd-1
stset rttdEnd,failure(event) exit(time .) id(id) enter(rttdStart)

stcox trt, cluster(id) nolog
return scalar rttd_crude_est=_b[trt]
return scalar rttd_crude_CP=cond((_b[trt]-1.96*_se[trt])<ln(`hr') & (_b[trt]+1.96*_se[trt])>ln(`hr'),1,0)

stcox trt x3, cluster(id) nolog
return scalar rttd_x_est=_b[trt]
return scalar rttd_x_CP=cond((_b[trt]-1.96*_se[trt])<ln(`hr') & (_b[trt]+1.96*_se[trt])>ln(`hr'),1,0)

ereturn clear
  
end



*** set HR as a macro
*** HR = 0.5

global hr=0.5
 
display "$S_TIME  $S_DATE"

simulate, reps(1000) seed(123) dots(10): rttd_tcc, hr($hr) corr(0.5)
qui foreach m in crude x {
gen rttd_`m'_rmse=(exp(rttd_`m'_est)-$hr)^2
sum rttd_`m'_rmse
replace rttd_`m'_rmse=sqrt(r(mean))
gen rttd_`m'_bias=100*(exp(rttd_`m'_est)-$hr)/$hr
gen tos_`m'_rmse=(exp(tos_`m'_est)-$hr)^2
sum tos_`m'_rmse
replace tos_`m'_rmse=sqrt(r(mean))
gen tos_`m'_bias=100*(exp(tos_`m'_est)-$hr)/$hr
}
sum *est, sep(3)
sum *bias, sep(3)
sum *CP, sep(3) 
sum *rmse, sep(3)  
 
display "$S_TIME  $S_DATE"

simulate, reps(1000) seed(123) dots(10): rttd_tcc, hr($hr) corr(0.25)
qui foreach m in crude x {
gen rttd_`m'_rmse=(exp(rttd_`m'_est)-$hr)^2
sum rttd_`m'_rmse
replace rttd_`m'_rmse=sqrt(r(mean))
gen rttd_`m'_bias=100*(exp(rttd_`m'_est)-$hr)/$hr
gen tos_`m'_rmse=(exp(tos_`m'_est)-$hr)^2
sum tos_`m'_rmse
replace tos_`m'_rmse=sqrt(r(mean))
gen tos_`m'_bias=100*(exp(tos_`m'_est)-$hr)/$hr
}
sum *est, sep(3)
sum *bias, sep(3)
sum *CP, sep(3) 
sum *rmse, sep(3)  

display "$S_TIME  $S_DATE"

simulate, reps(1000) seed(123) dots(10): rttd_tcc, hr($hr) corr(0)
qui foreach m in crude x {
gen rttd_`m'_rmse=(exp(rttd_`m'_est)-$hr)^2
sum rttd_`m'_rmse
replace rttd_`m'_rmse=sqrt(r(mean))
gen rttd_`m'_bias=100*(exp(rttd_`m'_est)-$hr)/$hr
gen tos_`m'_rmse=(exp(tos_`m'_est)-$hr)^2
sum tos_`m'_rmse
replace tos_`m'_rmse=sqrt(r(mean))
gen tos_`m'_bias=100*(exp(tos_`m'_est)-$hr)/$hr
}
sum *est, sep(3)
sum *bias, sep(3)
sum *CP, sep(3) 
sum *rmse, sep(3)   
 
log close  
   
  