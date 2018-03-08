use diet,clear
gen rate=chd/y*1000
table hieng [iw=y], c(mean rate)
gen u=uniform()
keep if u<=0.25 | chd==1
table hieng [iw=y], c(mean rate)
gen invpw=cond(chd==1,1,4)
gen w=invpw*y
table hieng [iw=w], c(mean rate)
poisson chd hieng [iw=invpw], e(y) robust
poisson, irr