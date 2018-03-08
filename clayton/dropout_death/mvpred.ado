*! Version 1.0

program define mvpred
syntax varlist [if] [in] [,PRed(string) TOL(real 0.0001) MAXit(integer 20)]
if "`pred'"!="" {
	confirm new var `pred'
}
tokenize `varlist'
local mvar `1'
mac shift
local indv
local cons = 1
while "`1'"!="" {
	local indv "`indv' `1'"
	mac shift
	local cons = `cons' + 1

}
di "Logistic regression to predict missing values of " _continue 
di "`mvar' " _continue
di "(Using `mvar' as offset)"
tempvar touse miss offs exb ome xb w hw
mark `touse' `if' `in' 
markout `touse' `indv'
gen `miss' = (`mvar' == .)
quietly {

	/* First approximation 
	   Gaussian imputation of mvar, followed by logit  

	regr `mvar' `indv' if (`touse' & !`miss')
	local rvar = e(rmse) /* Residual variance */
	predict `w' if `touse'
	gen `offs' = cond(`miss', `w'+`rvar', `mvar')
	logit `miss' `indv' if `touse', offset(`offs')
	drop `w'
	predict `w' if `touse', p
	matrix b = e(b)

	*/

	/* First approximation --- MAR */

	logit `miss' `indv' if `touse'
	predict `w' if `touse', p
	matrix b = e(b)
	replace `w' = 1 - `w'

	/* Iterative refinement */

	local nocnv = 1
	local iter = 1
	while `nocnv' {
		if `iter'>1 {
			drop `exb' `hw' `ome'
		}
		matrix score `exb' = b 
		replace `exb' = `exb' + `mvar' 
		replace `exb' = exp(`exb') 
		gen `hw' = `w'*`exb' 
		/* Adjust estimate of _cons */
		summ `hw' if !`miss'		
		local s1 = r(sum)
		summ `w' if `miss'
		local s2 = r(sum)
		local my = `s2'/`s1'
		matrix b[1,`cons'] = b[1,`cons'] + ln(`my')
		replace `hw' = `hw'*`my'
		/* Newton's method to adjust remaining coefficients */
		gen `ome' = cond(`miss', -`w', `hw')
		matrix vecaccum u = `ome' `indv' if `touse'
		matrix accum H = `indv' [aw = `hw'] if `touse' & !`miss'
 		matrix db = u*syminv(H)
		matrix b = b - db
		matrix qf = db*u'
		local mqf = qf[1,1]
		local nocnv = `mqf'>`tol'
		if `iter'==`maxit' & `nocnv' {
			di in re "No convergence after `maxit' iterations"
			exit
		}
		noi di "Iteration `iter': convergence test = `mqf'"
		local iter = `iter' + 1
	}
	/* Save predicted probabilities */
	if "`pred'"!="" {
		gen `pred' = 1/(1 + `exb') if `touse'
	}
	/* Robust variance estimate */
	matrix H = syminv(H)
	_robust `ome' [aw=`w'], v(H)
	estimates post b H, esample(`touse')
}
estimates display
end



