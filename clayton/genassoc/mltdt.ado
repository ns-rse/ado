*! Version 2.2 DGC Apr 27, 2004
program define mltdt, rclass sortpreserve
syntax [varlist (default=none)] [using /] [fw pw iw] [, GType POSTfix(string)/*
          */ Ped(varname) Id(varname) MOther(varname) FAther(varname) /*
          */ SEX(varname) AFfect(varname) ACode(numlist) noWARN  /*
          */ RObust CLuster(string) noIMPute NOIsily /*
          */ MOntecarlo(integer 0) DOTS(integer -1) /*
          */ SAVing(string) REPlace APPEND TRanspose]

/* Variables, codes etc. */

unab vars : _all
tokenize `vars'
if "`ped'"=="" {
  local ped "`1'"
}
if "`ped'"=="" {
  local ped "`1'"
}
if "`id'"=="" {
  local id "`2'"
}
if "`mother'"=="" {
  local mother "`4'"
}
if "`father'"=="" {
  local father "`3'"
}
if "`sex'"=="" {
  local sex "`5'"
}
if "`affect'"=="" {
  local affect "`6'"
}
if "`acode'"=="" {
  local acode 2
}
if "`cluster'"!="" {
  local robust "robust"
}
if "`robust'"!="" {
  if "`cluster'"=="" {
    local cl "`ped' `mother' `father'"
  }
  else {
    local cluster = substr("`cluster'", 1, 1)
    if "`cluster'"=="n" {
      local cl "`ped' `mother' `father'"
    }
    else if "`cluster'"=="p" {
      local cl "`ped'"
    }
    else {
      di in red "Invalid cluster() option"
      exit
    }
  }
}
if `dots'>0 & `montecarlo'==0 {
  di in red "dots option only used in conjunction with montecarlo option"
  exit
}

/* Weights */

if "`weight'"!="" {
  tempvar wt
  quietly gen `wt' `exp'
}

/* Saving options */
  
if "`replace'"!="" & "`append'"!="" {
  di in red "replace and append options cannot be used together"
  exit
}
if "`transpose'"!="" & "`saving'"=="" {
  di in red "transpose option requires saving option to be set"
  exit
}
if "`append'"!="" & "`saving'"=="" {
  di in red "append option requires saving option to be set"
  exit
}                                                         

/* New variable names */
  
tempname lstub mstub fstub tustub
local cstub _t`tustub'
local ustub _u`tustub'
if "`postfix'"=="" {
  local pf1 _1
  local pf2 _2
}
else {
  local npf: word count `postfix'
  if `npf'!=2 {
    di in red "Postfix option must contain two strings"
    exit
  }
  local pf1 : word 1 of `postfix'
  local pf2 : word 2 of `postfix'
}
if "`saving'"!="" {
  tempname post
  tempfile save
  postfile `post' htSNP impute score variance t p using "`save'", /*
                  */ double `replace'
}
local nv: word count  `varlist'
tokenize `varlist'

/* Type of pedigree variable */

local pedtype : type `ped'

/* Find which snps to analyses */

local maxln 6
local htsnps
local topred
local loc 0
preserve
if "`using'"=="" {
  if `nv'==1 {
    di in red "You must indicate SNPs to be tested"
    exit
  }
  if "`gtype'"!="" {
    while "`1'"!="" {
      local loc = `loc'+1
      local newv `lstub'_`loc'
      quietly count if `1'!=0 & `1'!=1 & `1'!=2 & `1'!=.
      if r(N)!=0 {
        di "Variable {bf:`1'} not coded 0, 1, or 2 --- omitted"
      }
      local len = length("`1'")
      if `len'>`maxln' {
        local maxln `len'
      }
      local htsnps "`htsnps' `1'"
      mac shift
    }
  }
  else {
    if mod(`nv',2) {
      di in red "Genotypes expected as {bf:pairs} of variables"
      exit
    }
    while "`1'"!="" {
      quietly count if (`1'!=1 & `1'!=2 & `1'!=.)|(`2'!=1 & `2'!=2 & `2'!=.)
      if r(N)!=0 {
        di "Variable {bf:`1'} and/or {bf:`2'} not coded 1 or 2 --- omitted"
      }
      else {
        local loc = `loc'+1
        local newv `lstub'_`loc'
        quietly gen `newv' = `1'+`2'-2
        local len 0
        while substr("`1'", `len'+1, 1)==substr("`2'", `len'+1, 1) {
          local len = `len' + 1
        }
        if `len'==0 {
          local htsnps "`htsnps' `1'"
        }
        else {
          local last = substr("`1'", `len', 1)
          if "`last'"=="_" {
            local len = `len' - 1
          } 
          local stub  = substr("`1'", 1, `len')
          local htsnps "`htsnps' `stub'"
        }
      }
      mac shift
      mac shift
    }
  }
  local nht : word count `htsnps'
  local ntp  0
}
else {
  if `nv'!=0 {
    di in red "Too many arguments"
    exit
  }
  quietly {
    use `using'
    drop _hap_* _hcl_*
    unab topred : _all_*
    unab htsnps : _all
    foreach var of local topred {
      local htsnps : subinstr local htsnps " `var'" ""
    }
    local topred : subinstr local topred "_all_" "", all
    local ntp : word count `topred'
    local nht : word count `htsnps'
    local nb = 1 + `nht'
    tempname b beta
    matrix `beta' = J(`ntp', `nb', 0)
    local row 0
    foreach var of local topred {
      local len = length("`var'")
      if `len'>`maxln' {
        local maxln `len'
      }
      reg _all_`var' `htsnps'
      local row = `row'+1
      matrix `b' = e(b)
      local sb 0
      forvalues col = 1/`nht' {
        local sb = `sb' + `b'[1,`col']
        matrix `beta'[`row',`col'] = `b'[1, `col']
      }
      matrix `beta'[`row',`nb'] = 2*(`sb' + `b'[1,`nb'])
    }
    restore
    preserve
    local loc 0
    foreach var of local htsnps {
      local loc = `loc'+1
      local newv `lstub'_`loc'
      if "`gtype'"!="" {
        count if (`var'!=.) &  (`var'!=0) &  (`var'!=1) &  (`var'!=2)  
        if r(N)!=0 {
          noi di "htSNP genotype {bf:`1'} not coded 0, 1 or 2"
          exit
        }
        ren `var' `newv'
      }
      else {
        local v1 `var'`pf1'
        local v2 `var'`pf2'
        count if (`v1'!=1 & `v1'!=2 & `v1'!=.)|(`v2'!=1 & `v2'!=2 & `v2'!=.)
        if r(N)!=0 {
          noi di "htSNP allele {bf:`1'} and/or {bf:`2'} not coded 1 or 2"
          exit
        }
        else {
          quietly gen `newv' = `v1'+`v2'-2
        }
      }
    }
  }
}
di
if "`topred'" != "" {
  di "Haplotype tagging SNPs: `htsnps'"
  di "SNPs to be predicted: `topred'"
}
else {
  di "SNPs to be tested: `htsnps'"
}

/* Organize trios as one line */

tempvar use mid fid cid
quietly recode `affect' `acode' = 2 * = 1 
keep `ped' `id' `mother' `father' `affect' `lstub'_*
tempfile idord
sort `ped' `id'
quietly save "`idord'"

gen byte `use' = cond(`affect'==2, 1, 0)
markout `use'  `mother' `father'  
quietly keep if `use'
quietly drop `affect' `use'
di
di "Total number of affected offspring = " _N
rename `mother' `mid'
rename `father' `fid'
rename `id' `cid'
renpfix `lstub' `cstub'

rename `mid' `id'
sort `ped' `id'
foreach var of local mvars {
  quietly replace `var'=.
}
quietly merge `ped' `id' using "`idord'", update
quietly drop if _merge==2
quietly drop _merge `mother' `father'
renpfix `lstub' `mstub'
rename `id'  `mid'

rename `fid' `id'
sort `ped' `id'
foreach var of local fvars {
  quietly replace `var'=.
}
quietly merge `ped' `id' using "`idord'", update
quietly drop if _merge==2
quietly drop _merge `mother' `father'
renpfix `lstub' `fstub'
rename `id'  `fid'

rename `cid' `id'
rename `mid' `mother'
rename `fid' `father'
quietly drop `affect'

tempvar cmp mih mall sum
gen byte `cmp' = 1
markout `cmp' `mstub'* `fstub'*

/* Check for unknown genotypes and misinheritances
   Generate transmitted and untransmitted counts for htsnps */

gen byte `mih'=0
gen byte `mall'=1
forvalues i = 1/`nht' {
  local m `mstub'_`i'
  local f `fstub'_`i'
  local c `cstub'_`i'
  local u `ustub'_`i'
  quietly {
    gen `sum' = `m' + `f'
    replace `c' = `sum'/2 if (`c'==.)&(`m'!=.)&(`f'!=.)&(`m'!=1)&(`f'!=1)
    replace `c' = . if `sum'==.
    gen `u' = `sum' - `c'
    replace `cmp' = 0 if `c'==.
    replace `mih' = 1 if (`c'!=.) & (`m'!=.) & (`f'!=.) & ( /*
                   */ ((`m'==0)&(`f'==2)&(`c'!=1))| /*
                   */ ((`m'==2)&(`f'==0)&(`c'!=1))| /*
                   */ (`c'>`sum') | (`c'<(`sum'-2)) )
    replace `mall' = 0 if (`c'!=.) & (`m'!=.) & (`f'!=.)
    drop `m' `f' `sum'
  }
}
quietly count if `mall'
if r(N)>0 {
  di "There are " r(N) " trios which could not be scored at at least one locus"
}
quietly count if `mih'
if r(N)>0 {
  di "There are " r(N) " complete trios with misinheritance at at least one locus"
  if "`warn'"=="" {
    li `ped' `id' `mother' `father' if `mih'
    di "{sf}" 
  }
}
quietly drop if `mih'|`mall'

di "Number of trios analysed = " _N
quietly {
  tempvar incmp
  gen byte `incmp' = !`cmp'
  count if `incmp'
  local Ninc = r(N)
}
if "`robust'"!="" {
  tempvar last lastcmp
  quietly {
    sort `cl' `incmp'
    by `cl' `incmp': gen byte `lastcmp' = (_n==_N) & `cmp'
    by `cl': gen byte `last' = (_n==_N)
    count if `last'
    noi di "Number of clusters = " r(N)
    count if `lastcmp'
    noi di "Number with complete data for at least one trio = " r(N)
  }
}

/* Tests */
  
di
di "Single locus tests " _continue
if "`topred'"!="" {
  di "for tagging SNPs" _continue
}
else {
  di "for typed loci" _continue
}
di " using trios with complete data at the locus"
di "(Transmissions of allele 2 at each locus)"
di
di %20s "SNP" %8s "N" %13s "Transmitted" %15s "Untransmitted" _continue
di %12s "z-value" %12s "p-value"
local minp 1.0
local locmin
local nloc =`nht'+`ntp'
local allloc "`htsnps' `topred'"
tempvar tw uw d d2
local cvars
local uvars
forvalues i = 1/`nht' {
  quietly {
    local c `cstub'_`i'
    local u `ustub'_`i'
    local cvars "`cvars' `c'"
    local uvars "`uvars' `u'"
    local locus : word `i' of `allloc'
    if "`weight'"!="" {
      gen `tw' = `c'*`wt'
      gen `uw' = `u'*`wt'
    }
    else {
      gen `tw'=`c'
      gen `uw'=`u'
    }
    count if `c'!=. & `u'!= .
    local en = r(N)
    summ `tw', meanonly
    local sc = r(sum)
    summ `uw', meanonly
    local su = r(sum)
    gen `d' = `tw' - `uw'
    replace `d' = 0 if `c'==.
    if "`robust'"!="" {
      by `cl': gen `d2' = sum(`d')
      replace `d' = `d2'
      replace `d2' = `d'^2
      summ `d' if `last', meanonly
      local score = r(sum)
      summ `d2' if `last', meanonly
      local var = r(sum)
    }
    else {
      gen `d2' = `d'^2
      summ `d', meanonly
      local score = r(sum)
      summ `d2', meanonly
      local var = r(sum)
    }
    drop `tw' `uw' `d' `d2' 
    local zv = `score'/sqrt(`var')
    local pv = 2*(1-norm(abs(`zv')))
  }
  di %20s "`locus'" %8.0f `en' _skip(3) %10.3f `sc' _skip(5) _continue
  di  %10.3f `su' _skip(2) %10.3f `zv' _skip(2) %10.4g `pv'
  if `pv'<`minp' {
    local minp `pv'
    local locmin `locus'
  }
  if "`saving'"!="" {
    local ifht = (`i'<=`nht')
    post `post' (1) (1) (`score') (`var') (`zv') (`pv')
  }
  if `pv'<`minp' {
    local minp `pv'
    local locmin "`locus'"
  }
}

/* ????? */
  
if "`impute'"=="" | "`using'"!="" {
  di
  di "Smallest p-value is " `minp' " (`locmin')"
  return scalar p_min_ty = `minp'
  return local locus_p_min_ty = "`locmin'"
}

/* Imputation */

if "`impute'"=="" & `Ninc'>0 {
  di
  di "Number of trios with incomplete data at any locus = " `Ninc'
  di 
  di "Imputing missing data using Stata {bf:impute} command: " _continue
  if "`noisily'"!="" {
    di
  }
  local alist
  local ilist
  forvalues i = 1/`nht' {
    local alist "`alist'  @`tustub'_`i'"
    local ilist "`ilist' `tustub'_`i'"
  }
  tempvar tr 
  tokenize `allloc'
  quietly {
    reshape long `alist', i(`ped' `id') j(`tr') string
    local done
    foreach i of local ilist {
      local others : subinstr local ilist "`i'" "", word
      noi di "[`1']" _continue 
      mac shift
      tempvar im
      `noisily' impute `i' `others', gen(`im')
      local done "`done' `im'"
    }
    tokenize `done'
    foreach i of local ilist {
      drop `i'
      rename `1' `i'
      mac shift
    }
    reshape wide `alist', i(`ped' `id') j(`tr') string
  }
  di
  local impt 3
  quietly {
    if "`robust'"!="" {
      replace `lastcmp' = `last'
    }
    replace `cmp' = 1
  }
}
else {
  di
  di "Number of trios with incomplete data at any locus = " `Ninc'
  if `Ninc'>0 {
    di "(these are omitted for the remaining analyses)"
    local impt 2
  }
}
if "`robust'"!="" {
   sort `cl' `incmp' `lastcmp' `last'
}


/* Calculate predicted transmitted and untransmitted counts for other loci */

if `ntp'>0 {
  matrix colnames `beta' = `cvars' _cons
  forvalues i = 1/`ntp' {
    local loc = `nht' + `i'
    matrix `b' = `beta'[`i',1...]
    matrix score `cstub'_`loc' = `b'
  }
  matrix colnames `beta' = `uvars' _cons
  forvalues i = 1/`ntp' {
    local loc = `nht' + `i'
    matrix `b' = `beta'[`i',1...]
    matrix score `ustub'_`loc' = `b'
  }
}

di
di %20s "SNP" %8s "N" %13s "Transmitted" %15s "Untransmitted" _continue
di %12s "z-value" %12s "p-value"
local minp 1
forvalues i = 1/`nloc' {
  quietly {
    local c `cstub'_`i'
    local u `ustub'_`i'
    local locus : word `i' of `allloc'
    if "`weight'"!="" {
      replace `c' = `c'*`wt'
      replace `u' = `u'*`wt'
    }
    summ `c' if `cmp', meanonly
    local sc = r(sum)
    summ `u' if `cmp', meanonly
    local su = r(sum)
    replace `c' = (`c' - `u')/2 
    count if `cmp' &  `c' != .
    local en = r(N)
    replace `c' = 0 if (!`cmp') | (`c'==.)
    if "`robust'"!="" {
      by `cl': replace `u' = sum(`c')
      replace `c' = `u'
      replace `u' = `c'^2
      summ `c' if `last', meanonly
      local score = r(sum)
      summ `u' if `last', meanonly
      local var = r(sum)
    }
    else {
      replace `u' = `c'^2
      summ `c', meanonly
      local score = r(sum)
      summ `u', meanonly
      local var = r(sum)
    }
    local zv = `score'/sqrt(`var')
    local pv = 2*(1-norm(abs(`zv')))
  }
  if `pv'<`minp' {
    local minp `pv'
    local locmin `locus'
  }
  di %20s "`locus'" %8.0f `en' _skip(3) %10.3f `sc' _skip(5) _continue
  di  %10.3f `su' _skip(2) %10.3f `zv' _skip(2) %10.4g `pv'
  if "`saving'"!="" {
    local ifht = (`i'<=`nht')
    post `post' (`ifht') (`impt') (`score') (`var') (`zv') (`pv')
  }
}
di
di "Smallest p-value (all loci) is " `minp' " (`locmin')"
di

quietly {
  tempvar one
  gen `one'=1
  tempname U V A Q
  if "`robust'"!="" {
    matrix accum `V' = `cstub'_1 - `cstub'_`nht' if `lastcmp', noconst
    matrix vecaccum `U' = `one' `cstub'_1 - `cstub'_`nht' /*
         */ if `lastcmp' , noconst
  }
  else {
    matrix accum `V' = `cstub'_1 - `cstub'_`nht' if `cmp', noconst
    matrix vecaccum `U' = `one' `cstub'_1 - `cstub'_`nht' if `cmp', noconst
  }
  matrix colnames `U' = `htsnps'
  matrix colnames `V' = `htsnps'
  matrix rownames `V' = `htsnps'
  matrix `A' = syminv(`V')
  matrix `Q' = `U'*`A'*`U' '
  local chi2 = `Q'[1,1]
  local df = `nht' - diag0cnt(`A')
  local Cp = chi2tail(`df', `chi2')
}
di "Global test for all `nht' htSNP loci: Chi-squared (" `df' " df) " _continue
di `chi2' ", p = " `Cp'
if `Ninc'>0 {
  di "Note that there were incomplete data for " `Ninc' " trios"
  if "`complete'"!="" {
    di "These have been excluded"
  }
  else {
    if "`impute'" == "" {
      di "Missing scores were imputed"
    }
    else {
      di "Missing score contributions have been set to zero"
    }
  }
}
return scalar p_min_all = `minp'
return scalar p_global = `Cp'
return scalar df = `df'
return scalar chi2 = `chi2'
return local locus_p_min_all = "`locmin'"
if "`using'"!="" {
  tempname I T
  matrix `I' = I(`nht')
  matrix rownames `I' = `htsnps'
  matrix colnames `I' = `htsnps'
  matrix rownames `beta' = `topred' 
  matrix colnames `beta' = `htsnps' _cons
  matrix `T' = `I' \ `beta'[1...,1..`nht']
  matrix drop `b' `beta' `I'
}
if `montecarlo'>0 {
  tempname VU UT X2 P RU 
  tempvar rsign
  if "`using'"!="" {
    matrix `VU' = vecdiag(`T'*`V'*`T'')
    matrix `UT' = `U'*`T''
  }
  else {
    matrix `UT' = `U'
    matrix `VU' = vecdiag(`V')
  }
  matrix `X2' = J(1,`nloc',0)
  matrix `P' = J(1, `nloc',0)
  local maxX2 0
  forvalues i = 1/`nloc' {
    local X2i =  `UT'[1,`i']^2/`VU'[1,`i']
    matrix `X2'[1,`i'] = `X2i'
   di `i' " = " `X2i' 
    if `X2i'>`maxX2' {
      local maxX2 `X2i'
      local locX2 : word `i' of `allloc'
      di "`locX2'"
    }
  }
  di "Most significant single locus: `locX2', chi-squared on 1 df = " `maxX2' ")"
  di
  di "Monte Carlo test, using `montecarlo' random sign changes ... please wait"
  local mc_global 0
  local mc_most 0
  forvalues sim = 1/`montecarlo' {
    quietly {
      gen `rsign' = cond(uniform()>0.5, 1, -1)
      if "`robust'"!="" {
        matrix vecaccum `RU' = `rsign' `cstub'_1 - `cstub'_`nht' /*
              */  if `lastcmp' , noconst
      }
      else {
        matrix vecaccum `RU' = `rsign' `cstub'_1 - `cstub'_`nht' if `cmp' /*
              */ , noconst
      }
      drop `rsign'
      matrix `Q' = `RU'*`A'*`RU''
      if `Q'[1,1]>=`chi2' {
        local mc_global = `mc_global'+1
      }
      if "`using'"!="" {
        matrix `RU' = `RU'*`T''
      }
      local notyet 1
      forvalues i = 1/`nloc' {
        local X2i =  `RU'[1,`i']^2/`VU'[1,`i']
        matrix `X2'[1,`i'] = `X2i'
        if `notyet' & `X2i'>=`maxX2' {
          local mc_most = `mc_most' + 1
          local notyet 0
        }
      }
    }
    if mod(`sim',`dots')==0 {
      di "." _continue
    }
  }
  di
  di "Monte Carlo estimates of permutation p-values"
  di "Global test: " `mc_global'/`montecarlo'
  di "Most significant single locus: " `mc_most'/`montecarlo'
  di "(Locus `locX2', chi-squared on 1 df = " `maxX2' ")"
}

if "`using'"!="" {
  return matrix expand `T'
}
return matrix score_variance `V'
return matrix score `U'

if "`saving'"!="" {
  quietly {
    post `post' (.) (.) (.) (.) (.) (`Cp')
    postclose `post'
    use "`save'", clear
    gen str`maxln' locus = ""
    if "`impute'"=="" & `Ninc'>0{
      forvalues i = 1/`nloc' {
        local loci : word `i' of `allloc'
        replace locus = "`loci'" if _n==`i'
      }
    }
    else {
      forvalues i = 1/`nloc' {
        local loci : word `i' of `allloc'
        replace locus = "`loci'" if (_n==(`nht'+`i')) | /*
                                 */ (`i'<=`nht'  & _n==`i')
      }
    }
    replace locus = "_Global" if _n==_N
    sort locus
    if "`transpose'"!="" {
      local lsrt 
      local nl = _N
      forvalues i=1/`nl' {
        local li = locus[`i']
        if impute[`i'] {
          local li "`li'_I"
        }
        local lsrt "`lsrt' `li'"
      }
      drop locus
      xpose , clear varname
      ren _varname statistic
      tokenize "`lsrt'"
      local loc 0
      while "`1'"!="" {
        local loc = `loc' + 1
        ren v`loc' `1'
        mac shift
      }
    }
  }
  di
  if "`append'"!="" {
    local fxt = substr("`saving'", -4, 4)
    if "`fxt'"!=".dta" {
      local saving "`saving'.dta"
    }
    cap confirm file `saving'
    if _rc==0 {
      di "Appending results to file {bf:`saving'}"
      append using `saving'
      local replace replace
    }
    else {
      di "File {bf:`saving'} does not yet exist"
    }
  }
  save `saving', `replace'
  di "(File coutains " _N " records)"
}
restore
end 
