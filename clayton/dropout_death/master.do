* Decide whether analysis is "mortal" or "immortal"

* global MORTAL 1
* global MORTAL 0

if $MORTAL {di "Mortal analysis"} else {di "Immortal analysis"}

* Data file: subset of monotone missing data from CC75C

use cc75css, clear

*
* Generate variables for MMSE and age at _previous_ visit
*

sort id wave
by id: gen mmse_p = mmse[_n-1]
by id: gen age_p = age[_n-1]

*
* Calculate weights. For immortal analysis model probability of being seen
* conditional upon having been seen at the previous visit. For mortal analysis
* model this probability after also conditioning upon the subject being alive
* at the date the visit was due. 
*

gen seen = (outcome == 1)
gen wave2 = (wave==2)
gen wave3 = (wave==3)
if !$MORTAL {
  logit seen wave2 wave3 age_p mmse_p if wave>0, nocons
}
else {
  logit seen wave2 wave3  age_p mmse_p  if wave>0 & outcome!=3
}

* Cumulative "survival" probabilities, and inverse probability weights

predict pseen
replace pseen = 1 if wave==0
sort id wave
by id: gen surv = exp(sum(log(pseen)))
gen weight = 1/surv

* Drop temporary variables used in calculating weights

drop wave2 wave3 mmse_p age_p seen pseen surv

* Now continue teh analysis only on persons seen

keep if outcome==1

* Some variable derivations

by id: gen age_entry = age[1]
* One person was (just) <75 at entry --- make him/her 75!
replace age_entry = 75 if age_entry < 75
egen cohort = cut(age_entry), at(75,80,85,120)
gen agem75 = age-75
gen malexage = male*agem75
gen cht75 = (cohort==75)
gen cht80 = (cohort==80)
gen cht85 = (cohort==85)

* Weighted means of age and mmse by cohort and wave

table cohort wave  [pw=weight] if male, contents(freq mean age mean mmse)
table cohort wave  [pw=weight] if !male, contents(freq mean age mean mmse)

* Save for plotting

preserve
collapse (mean) mmse age [pw=weight], by(male cohort wave)
list
outsheet using cohort_means.dat, replace
restore

* Regression analysis (Table 3 of paper 1)

reg mmse cht75 cht80 cht85 male agem75 malexage [pw=weight], nocons cluster(id)

