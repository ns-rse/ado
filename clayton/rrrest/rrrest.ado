*! Name    : rrrest.ado (estimate relative recurrence risk ratio)
*! Date    : 11 Oct 2002
*! Version : 0.1
*! Author  : Chris Wallace
*! Email   : chris.wallace@lshtm.ac.uk

program define rrrest, eclass
version 7
	set trace off
	if replay() {
		if "`e(cmd)'" ~= "rrres_lf" {
			error 301 /* last-estimates-not-found error */
			}
		Replay `0'
		}
	else  Estimate `0'
end

program define Replay
	/* display estimation results respecting the display options */
	syntax [, Level(int $S_level) ]
	ml display, level(`level')
end

program define Estimate, eclass
	/* parse equations */
	forvalues i=1/3 {
		gettoken (global)RRRv`i' 0:0, match(parens) /* lists within brackets */
		}
	forvalues i=1/2 {
		gettoken (global)RRRa`i' (global)RRRc`i': (global)RRRv`i', parse("=") /* aff & cov for pair */
		global RRRc`i' :  subinstr global RRRc`i' "=" ""
		}

	/* process remaining options */
	syntax [if] [in] [fweight pweight] [, Robust CLuster(varname num) SCore(passthru) /*
	*/ Level(passthru) CONStraint(string) DIFFicult /*
	*/ CI NSim(integer 1000) ]

	if "`cluster'"~="" {
		local clopt "cluster(`cluster')"
		}
	if "`constraint'"~="" {
		local consopt "constraint(`constraint')"
		}
	if "`difficult'"~="" {
		local difficult "difficult"
		}

	/* identify estimation subsample */
	qui marksample touse
	qui markout `touse' `cluster', strok

	/* perform estimation using ml */
	di _n "Fitting constant-only model:"
	ml model lf rrres_lf ($RRRa1=) ($RRRa2=) () [`weight'`exp'] if `touse', /*
	*/ maximize search(quietly) nopreserve

	di _n "Fitting full model:"
	ml model lf rrres_lf ($RRRa1=$RRRc1) ($RRRa2=$RRRc2) ($RRRv3) [`weight'`exp'] if `touse', /*
	*/ continue maximize `difficult' search(off) nopreserve `consopt'/*
	*/ `robust' `clopt' `score'

	/* display estimation results */
	estimate local cmd "rrrest"
	Replay

	/* calculate and display parameters of interest */
	rrres_sum $RRRa1 $RRRa2 [`weight'`exp'] /*u1(`u1') u2(`u2') th(`th') del(`delta') l1(`l1') l2(`l2') */

	/* do we need to calculate CI's ? */
	if "`ci'"~="" {
		Rci, ns(`nsim')
		}
end

program define Rci, rclass
	di _n "Calculating CIs for lambda, delta"
	syntax [, NSim(integer 1000) ]

	local u1 "u1"
	local u2 "u2"
	local delta "delta"
	local l1 "lambda1"
	local l2 "lambda2"
	local th "th"

	tempvar gp tg yp1 yp2 yp3 ii con
	tempname resdata   	/* results dataset */
	tempfile resfile 	   /* stored in  resfile  */

	mat V=e(V)  /* var-covar matrix */
	postfile `resdata' `ii' `delta'_lci `delta'_uci `l1'_lci `l1'_uci `l2'_lci `l2'_uci using `resfile', double

	forvalues i=1/3 {
		capture drop `yp`i''
		predict double `yp`i'', xb eq(eq`i')
		}
	
	/* only want pairs with at least one affected for CI estimates - drop the other obs, and add them back at end */
	egen `ii' = group($RRRc1 $RRRc2 $RRRc3) /* will merge according to this var */
	loc nwords : word count $RRRc1 $RRRc2 $RRRc3
	if "`nwords'" == "" {
		Rstderr `yp1' `yp2' `yp3' if ($RRRa1!=0 & $RRRa2!=0) ,  ns(`nsim') filename(`resdata')
		}
	else {
		bys $RRRc1 $RRRc2 $RRRc3: Rstderr `yp1' `yp2' `yp3' if ($RRRa1!=0 & $RRRa2!=0) ,/*
		*/ ns(`nsim') filename(`resdata')
		}


	/* JOIN NEW DATA */

	qui {
		cap restore, not /* prevent further restores wiping out new data */
		postclose `resdata'
		preserve
		sort `ii'
		use `resfile', clear
		sort `ii'
		save `resfile', replace
		restore
		sort `ii'
		merge `ii' using `resfile'
		drop _merge
		}
end

program define Rstderr, byable(recall, noheader)

	syntax varlist(min=3 max=3) [if] [in] , Filename(string) [NSim(integer 1000)] 
	qui marksample touse, novarlist
	tempvar u1 u2 th delta l1 l2 rt  y1 y2 y3 cons
	set more off

	/* restrict to single observation of `touse' */
	preserve
	qui keep if `touse'
	if (_N>0) {		
		scalar qwerty = `ii'[1]
		qui gen `cons'=1
		mkmat `cons' $RRRc1 in 1, matrix(X1)
		mkmat `cons' $RRRc2 in 1, matrix(X2)
		mkmat `cons' $RRRc3 in 1, matrix(X3)
 		mkmat `varlist' in 1, matrix(y_pred)
		
		drop _all
		qui set obs `nsim'

		/* get the new var-covar matrix */
		mat RRR_VC=J(3,3,0)
		forvalues i=1/3 {
			forvalues j=1/3 {
				mat RRR_VC[`i',`j'] = X`i' * V["eq`i':","eq`j':"] * X`j''
			if `i' != `j' {
				mat RRR_VC[`j',`i'] = RRR_VC[`i',`j']
				}
			}
		}

	/* simulate variables with required var-covar structure, but zero mean */
	mkcorr RRR_VC
	
	/* generate variables, adding y_pred to set correct mean */
	
	gen double `u1'=exp(y1+y_pred[1,1])/(1+exp(y1+y_pred[1,1]))
	gen double `u2'=exp(y2+y_pred[1,2])/(1+exp(y2+y_pred[1,2]))
	gen double `th'=exp(y3+y_pred[1,3])
	gen double `rt'=sqrt(((`u1'+`u2')*(`th'-1)+1)^2-4*`u1'*`u2'*`th'*(`th'-1))
	gen double `delta' = ((`th'-1)*(`u1'+`u2')-`rt'+1)/(2*(`th'-1))
	gen double `l1' = `delta'/(`u1'*`u2')
	gen double `l2' = `delta'/(`u1'*`u2')
	/* vars with no CIs */
	qui for var `l1' `l2' `delta': qui centile X, centile(2.5 97.5) \ scalar X_lci = r(c_1) \ scalar X_uci = r(c_2)
		
	/* STORE RESULTS */

	post `filename' (qwerty) (`delta'_lci) (`delta'_uci) (`l1'_lci) (`l1'_uci) (`l2'_lci) (`l2'_uci)

	/* restore original data */
	restore
	}
end

/*

Want to simulate 3 normals from the variance-covariance matrix.  These can
then be used to generate simulates of delta, lambda and analysed to get the
standard deviations out.

This is edited from an archive message on the statalist.

*/

* capture program drop mkcorr
program define mkcorr
args mymat
tempname A

*qui {
	local k = rowsof(`mymat')
	matrix `A' = cholesky(`mymat')
	local i 1
	capture drop c1 c2 c3
	while `i'<=`k' {
		gen c`i' = invnorm(uniform())
		local i=`i'+1
		}
	local i 1
	while `i'<=`k' {
		matrix row = `A'[`i',1...]
		matrix score y`i' = row
		local i=`i'+1
		}
	local i 1
	while `i' <= `k' {
		drop c`i'
		local i=`i'+1
		}
*	}
end
