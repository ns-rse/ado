use pngsubj, clear
merge id using pngevent
stset timeout, f(diag) or(timein) sc(365.25) exit(time .)
gen dalli=_d if _st==1
gen y=_t - _t0 if _st==1
stset timeout, f(diag==80,81) or(timein) sc(365.25) exit(time .)
gen dalri=_d if _st==1
collapse (sum) dalli dalri y, by(id vacc)
