use diet,clear
stset dox, f(chd) or(dob) en(doe) sc(365.25) id(id)
stsplit ageband, at(40(5)70) trim
keep if _st==1
replace chd=_d 
replace y=_t-_t0 
drop _*
save diet_age, replace

