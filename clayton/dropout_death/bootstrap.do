* Program to recalculate table 3 for bootstrap samples and post estimates
* to file. Note that we  need to generate new id's (bsid) in these samples
* or the weights will be incorrectly calculated

cap program drop bstab3
program define bstab3
if "`1'"=="?" {
  global S_1 "cohort1 cohort2 cohort3 male agem75 malexage"
  exit
}
else {
  quietly {
   sort bsid wave

   gen seen = (outcome == 1)
   gen wave2 = (wave==2)
   gen wave3 = (wave==3)
   if !$MORTAL {
     logit seen wave2 wave3 age_p mmse_p if wave>0, nocons
   }
   else {
     logit seen wave2 wave3  age_p mmse_p  if wave>0 & outcome!=3
   }
   predict pseen
   replace pseen = 1 if wave==0
   sort bsid wave
   by bsid: gen surv = exp(sum(log(pseen)))
   gen weight = 1/surv
   reg mmse cht* male agem75 malexage [pw=weight], nocons
  }
  post `1' (_b[cht75]) (_b[cht80]) (_b[cht85])  (_b[male]) (_b[agem75]) (_b[malexage])
 drop seen pseen surv weight
}
end

* Read data and derive variables

use cc75css, clear
sort id wave
by id: gen mmse_p = mmse[_n-1]
by id: gen age_p = age[_n-1]
by id: gen age_entry = age[1]
replace age_entry = 75 if age_entry < 75
egen cohort = cut(age_entry), at(75,80,85,120)
gen agem75 = age-75
gen malexage = male*agem75
gen cht75 = (cohort==75)
gen cht80 = (cohort==80)
gen cht85 = (cohort==85)

bstrap bstab3, dots reps(1000 ) cl(id) idcl(bsid) saving(bstab3) replace 

