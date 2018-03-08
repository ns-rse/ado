*! Name    : rrres_sum.ado (summarize results from rrrest.ado)
*! Date    : 11 Oct 2002
*! Version : 0.1
*! Author  : Chris Wallace
*! Email   : chris.wallace@lshtm.ac.uk

program define rrres_sum, byable(onecall)
	if _by() { 
		local by "bysort `_byvars': "
		}
	syntax varlist(min=2 max=2) [pweight fweight] [, u1(string) u2(string) DELta(string) TH(string) l1(string) l2(string) ]

	tokenize `varlist'
	local aff1 `1'
	local aff2 `2'
	
	/* use fweights for pweights */
	if "`weight'"~="" {
		di _n "Weighted summary statistics:"
		local weight="fweight" }
	
	/* names of output variables - change these if you want to change their names, but make sure rrrest.ado agrees */
	if "`u1'"=="" {
		local u1 "u1"	}
	if "`u2'"=="" {
		local u2 "u2"	}
	if "`delta'"=="" {
		local delta "delta"	}
	if "`th'"=="" {
		local th "theta"	}
	if "`l1'"=="" {
		local l1 "lambda1"	}
	if "`l2'"=="" {
		local l2 "lambda2"	}

	tempvar uu1 uu2
	tempvar rt eta se_uu1 se_uu2 se_eta
	tempvar lci_uu1 lci_uu2 uci_uu1 uci_uu2
	set more off
	/* generate results for whole dataset */
	quietly {
		for any `th' `u1' `u2' `th'*ci `u1'*ci `u2'*ci `delta' `l1' `l2', noheader: capture drop X
		predict double `uu1', equation(#1)
		predict double `se_uu1', equation(#1) stdp
		predict double `uu2', equation(#2)
		predict double `se_uu2', equation(#2) stdp
		predict double `eta', equation(#3)
		predict double `se_eta', equation(#3) stdp
		gen double `lci_uu1' = `uu1' - 1.96*`se_uu1'
		gen double `lci_uu2' = `uu2' - 1.96*`se_uu2'
		gen double `uci_uu1' = `uu1' + 1.96*`se_uu1'
		gen double `uci_uu2' = `uu2' + 1.96*`se_uu2'
		gen double `u1' = exp(`uu1')/(1+exp(`uu1'))
		gen double `u2' = exp(`uu2')/(1+exp(`uu2'))
		gen double `th'=exp(`eta')
  
		eta2u `lci_uu1' `u1'_lci
		eta2u `uci_uu1' `u1'_uci
		eta2u `lci_uu2' `u2'_lci
		eta2u `uci_uu2' `u2'_uci
		gen double `th'_lci = exp(`eta' - 1.96*`se_eta')
		gen double `th'_uci = exp(`eta' + 1.96*`se_eta')
  
		gen double `rt'=sqrt(((`u1'+`u2')*(`th'-1)+1)^2-4*`u1'*`u2'*`th'*(`th'-1))
		gen double `delta' = ((`th'-1)*(`u1'+`u2')-`rt'+1)/(2*(`th'-1))
		gen double `l1' = `delta'/(`u1'*`u2') if `aff1'==1
		gen double `l2' = `delta'/(`u1'*`u2') if `aff2'==1
		}
	/* summarize, taking account of by if necessary */
	`by' sum `th' `u1' `u2' `delta' `l1' `l2' [`weight'`exp']
end

program define eta2u
	args eta u
	capture drop `u'
	gen double `u' = exp(`eta') / (1+exp(`eta'))
end
