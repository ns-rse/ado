use pngsubj, clear
merge id using pngevent
drop _merge
stset timeout, f(diag) or(dob) en(timein) sc(365.25) id(id) exit(time .)
stsplit ageband, at(0(2)8)
gen dalli=_d if _st==1
gen y=_t - _t0 if _st==1
sort id ageband
save tmp,replace

use pngsubj, clear
merge id using pngevent
drop _merge
stset timeout, f(diag=80,81) or(dob) en(timein) sc(365.25) id(id) exit(time .)
stsplit ageband, at(0(2)8)
gen dalri=_d if _st==1
gen y=_t - _t0 if _st==1
sort id ageband

merge id ageband using tmp

collapse (sum) dalli dalri y, by(id ageband vacc)
