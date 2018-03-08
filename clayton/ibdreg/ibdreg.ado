!* Version 0.2 David Clayton Aug 1  2002
program define ibdreg, eclass
  version 7.0
  syntax [varlist (default=none)] [if] [in], [POsterior(varlist min=2 max=3) /*
     */ PRior(string) noCONStant OFFset(varname) Robust CLuster(varname) noLOG]
  marksample use                     
  if "`posterior'"!="" {
    global ibdreg_post "`posterior'"
  }
  else if "$ibdreg_post"=="" {
    di in red "No posterior IBD probabilities are defined"
    exit
  }
  if "`prior'"!="" {
    global ibdreg_prior "`prior'"
  }
  else if "$ibdreg_prior'"=="" {
    di "Sibling pairs are assumed"
    global ibdreg_prior "0.25 0.5 0.25"
  }  
  local npo : word count $ibdreg_post
  local npr : word count $ibdreg_prior
  if `npo'!=`npr' {
    di in red "Unequal number of prior and posterior IBD probabilities defined"
  }
  cap confirm var $ibdreg_post
  if _rc!=0 {
    di in red "Unknown variable(s) for posterior probabilities: " */
    */ "$ibdreg_post"
    exit
  }
  markout `use' $ibdreg_post
  local wprior : word 1 of $ibdreg_prior
  cap confirm number `wprior'
  if _rc==0 {
    forvalues i = 2/`npo' {
      local wprior : word `i' of $ibdreg_prior
      cap confirm number `wprior'
      if _rc!=0 {
        di in red "Syntax error defining prior probabilities: " /* 
        */ "$ibdreg_prior"
        exit
      }
    }
  }
  else {
    cap confirm var $ibdreg_prior
    if _rc!=0 {
      di in red "Unknown variable(s) defining posterior probabilities: " /*
      */ "$ibdreg_post"
      exit
    }
    markout `use' $ibdreg_prior
  }
  forvalues i = 1/`npo' {
    local wpost : word `i' of $ibdreg_post
    local wprior : word `i' of $ibdreg_prior
    quietly count if `use' & `wpost'!=. & (`wpost'<0 | `wpost'>1)
    if r(N) > 0 {
      di in red "Element(s) of {bf:`wpost'} not a probability"
      exit
    }
    quietly count if `use' & `wprior'!=. & (`wprior'<0 |  `wprior'>1)
    if r(N) > 0 {
      di in red "Element(s) of {bf:`wprior'} not a probability"
      exit
    }
    tempvar rat
    quietly gen `rat' = `wpost'/`wprior' if `use'
    local rats "`rats' `rat'"
  }

  global rat0 : word 1 of `rats'
  global rat1 : word 2 of `rats'
  global rat2 : word 3 of `rats'
  if "$rat2"=="" {
    tempvar rat
    quietly gen `rat' = .
    global rat2 "`rat'"
  }
  if "`offset'"!="" {
    local offset "offset(`offset')"
  }
  if "`cluster'"!="" {
    local robust "robust"
    local cluster "cluster(`cluster')"
  }

  /* ML estimation */

  ml model d2 _d2_ibd (ibd: `varlist', `constant' `offset') if `use',/*
    */  nopreserve maximize `robust' `cluster' `log'

  /* Print results */

  ml display
  local ll = e(ll)
  local df = e(k)
  if "`constant'"=="" {
    di
    di "Tests of all parameters, including constant (`df' df)"
    di "Likelihood ratio chi-squared =" %8.3f = 2*`ll' _continue
    di ", p =" %9.2g = chi2tail(`df', 2*`ll')
    quietly test _cons `varlist'
    di "Wald chi-squared             =" %8.3f = r(chi2) _continue
    di ", p =" %9.2g = chi2tail(`df', r(chi2))
  }
  else if "`robust'"=="" {
    di
    di "Likelihood ratio chi-squared (`df' df) =" %8.3f = 2*`ll' _continue
    di ", p =" %9.2g = chi2tail(`df', 2*`ll')
  }
  mac drop rat*
  estimates local cmd "ibdreg"
  end

