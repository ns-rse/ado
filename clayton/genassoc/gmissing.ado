*! Version 1.0, DGC Dec 2003
program def gmissing, sortpreserve
syntax [varlist (default=none)] [if] [in] [, POSTfix(string)]
marksample use, novarlist
if "`postfix'"=="" {
  local postfix _1 _2
}
else {
  local wc : word count `postfix'
  if `wc'!=2 {
    di in red "Two postfix strings expected"
    exit
  }
}
local pf1 : word 1 of `postfix'
local pf2 : word 2 of `postfix'
local nv: word count `varlist'
if mod(`nv',2) {
  di in red "Odd number of variables specified"
  exit
}
quietly gloci `varlist', postfix(`postfix')
local loci `r(loci)'
quietly count if `use'
local N = r(N)
local mvars
tempvar any
quietly gen byte `any'=0
foreach locus of local loci {
  tempvar miss
  gen byte `miss' = ((`locus'`pf1'==.)|(`locus'`pf2'==.))
  local mvars "`mvars' `miss'"
  quietly replace `any'=1 if `miss'
}
quietly count if `any' & `use'
local M = r(N)
sort `use' `mvars'
tempvar count
quietly {
  by `use' `mvars': gen `count' = cond(_n==_N & `use', _N, 0)
  gsort -`count'
}
di
di "Missing data for at least one locus in `M'/`N' (" _continue
di %5.1f 100*`M'/`N' "%) subjects:"
di
local nl : word count `mvars'
local i `nl'
while `i'>0 {
  local locus : word `i' of `loci'
  local miss : word `i' of `mvars'
  quietly count if `use' & `miss'
  local M = r(N)
  local space = 20 + `i'
  di %`space's "`locus'" _continue
  local j = `i'
  while `j'<`nl' {
    di "|" _continue
    local j = `j'+1
  }
  local i = `i'-1
  di %10.0f `M' " missing values (" %5.1f 100*`M'/`N' "%)"
}
di _skip(20) _continue
forvalues i=1/`nl' {
  di "|" _continue
}
di %10s "Count"
di
local i 1
while `count'[`i']!= 0 {
  di _skip(20) _continue
  foreach var of local mvars {
    if `var'[`i'] {
      di "-" _continue
    }
    else {
      di "+" _continue
    }
  }
  di %10.0f `count'[`i']
  local i = `i'+1
}
end
