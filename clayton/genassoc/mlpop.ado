*! Version 2.0, DGC Sept 2004
program define mlpop, rclass sortpreserve
syntax varlist (min=1) [using /], [SAVing(string)] [REPlace] [TRanspose] /*
  */ [APPEND] [GType] [POSTfix(string)] [noIMPute] [WIthin(varlist)] /*
  */ [AOV] [NOIsily]
marksample use
if "`aov'"!="" & "`within'"=="" {
  di "{bf:aov} option only operates in conjunction with {bf:within} option"
  exit
}
if "`noisily'"!="" & "`impute'"!="" {
  di "{bf:noisily} option only effective if {bf:impute} is in force"
}
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
if "`saving'"!="" {
  tempname post
  tempfile save
  postfile `post' htSNP use t_F N p using `save', double `replace'
}
local nv: word count  `varlist'
tokenize `varlist'
tempvar cc
local aff `1'
quietly gen `cc' = `aff'
mac shift

/* Find which snps to analyse */
  
local htsnps
local topred
local xvars
if "`using'"=="" {
  if `nv'==1 {
    di in red "You must indicate SNPs to be tested"
    exit
  }
  if "`gtype'"!="" {
    while "`1'"!="" {
      quietly count if `1'!=0 & `1'!=1 & `1'!=2 & `1'!=.
      if r(N)!=0 {
        di "Variable {bf:`1'} not coded 0, 1, or 2 --- omitted"
      }
      local htsnps "`htsnps' `1'"
      local xvars "`xvars' `1'"
      mac shift
    }
  }
  else {
    if !mod(`nv',2) {
      di in red "Genotypes expected as {bf:pairs} of variables"
      exit
    }
    while "`1'"!="" {
      quietly count if (`1'!=1 & `1'!=2 & `1'!=.)|(`2'!=1 & `2'!=2 & `2'!=.)
      if r(N)!=0 {
        di "Variable {bf:`1'} and/or {bf:`2'} not coded 1 or 2 --- omitted"
      }
      else {
        tempvar newv
        quietly gen double `newv' = `1'+`2'-2
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
          local stub = substr("`1'", 1, `len')
          local htsnps "`htsnps' `stub'"
        }
        local xvars "`xvars' `newv'"
      }
      mac shift
      mac shift
    }
  }
  markout `use' `cc' `within' `xvars'
  local nht : word count `htsnps'
  local ntp  0
  local xvht `xvars'
}
else {
  if `nv'!=1 {
    di in red "Too many arguments"
    exit
  }
  quietly {
    preserve
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
    matrix beta = J(`ntp', `nb', 0)
    local row 0
    foreach var of local topred {
      reg _all_`var' `htsnps' 
      local row = `row'+1
      matrix b = e(b)
      forvalues col = 1/`nb' {
        matrix beta[`row',`col'] = b[1, `col']
      }
    }
    restore
    foreach var of local htsnps {
      if "`gtype'"!="" {
        count if (`var'!=.) &  (`var'!=0) &  (`var'!=1) &  (`var'!=2)  
        if r(N)!=0 {
          noi di "htSNP genotype {bf:`1'} not coded 0, 1 or 2"
          exit
        }
        else {
          local xvars "`xvars' `var'"
        }
      }
      else {
        tempvar newv
        local v1 `var'`pf1'
        local v2 `var'`pf2'
        count if (`v1'!=1 & `v1'!=2 & `v1'!=.)|(`v2'!=1 & `v2'!=2 & `v2'!=.)
        if r(N)!=0 {
          noi di "htSNP allele {bf:`1'} and/or {bf:`2'} not coded 1 or 2"
          exit
        }
        else {
          tempvar newv
          quietly gen double `newv' = `v1'+`v2'-2
          local xvars "`xvars' `newv'"
        }
      }
    }
    matrix colnames beta = `xvars' _cons
    local xvht `xvars'
    markout `use' `cc' `within' `xvars'
    local row 0
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
if "`within'"!="" {
  di
  di "The analysis will be stratified by: `within'"
  di
  quietly {
    gsort `within' -`use'
    tempvar nw gm 
    by `within': gen int `nw' = _n
    count if `nw'==1
    local nstrata = r(N)
    noi di "Number of strata =" `nstrata'
    by `within': egen `gm' = mean(`cc') 
    if "`noisily'"!="" {
      char `gm'[varname] "Mean_`aff'"
      noi list `within' `gm' if `nw'==1, noobs subvarname
    }
    replace `cc' = `cc' - `gm'
    drop `gm'
    by `within': egen `gm' = sd(`cc') 
    count if `gm'==0 & `nw'==1
    if r(N)>0 {
      local novar = r(N)
      noi di "Warning: "  `novar' _continue
      noi di " strata are uninformative (no variation in `aff')"
      count if `gm'==0
      noi di "(" r(N) " observations lost in this analysis)"
      replace `use' = `use' & (`gm'>0)
      replace `cc' = . if `gm'==0
      gsort `within' - `use'
    }
    else {
      local novar 0
    }
    if "`aov'"!="" {
      noi {
        di
        di "Analysis of variance for variation of allele " _continue
        di "frequencies between strata"
        di
        di %10s "SNP" %15s "SSq Between" %15s "SSq Within" _continue
        di %15s "DF(num,den)" %10s "F-ratio" %10s "p-value"
      }
      tokenize `htsnps'
    }
    foreach xvar of local xvars {
      drop `gm'
      by `within': egen `gm' = mean(`xvar') 
      replace `xvar' = `xvar' - `gm'
      if "`aov'"!="" {
        summ `gm' if `xvar'!=.
        local ssb = r(Var) * (r(N) - 1)
        summ `xvar'
        local ssw = r(Var) * (r(N) - 1)
        local df1 = `nstrata' - 1
        local df2 = r(N) - `nstrata'
        local F = (`ssb'*`df2')/(`ssw'*`df1')
        local p = Ftail(`df1', `df2', `F')
        noi {
          di %10s "`1'" _continue
          di %15.3f `ssb' %15.3f `ssw' _continue
          di _skip(5) "(" %2.0f `df1' "," %5.0f `df2' ")" _continue
          di %10.3f `F' %10.2g `p'
        }
        mac shift
      }
    }
    drop `gm' `nw'
  }
  local nstrata = `nstrata' - `novar'
}
else {
  local nstrata 1
}
di
di "Single locus tests for tagging SNPs"
di "(using all available data for each SNP)"
di
di %20s "SNP" %15s "N" %15s "t-value" %15s "p-value"

local loc 0
local minp 1
tokenize `xvars'
foreach locus of local htsnps {
  local loc = `loc' + 1
  quietly {
    corr `cc' `1' 
    local r2 = r(rho)^2
    local df = r(N)-`nstrata'-1
    local tv = sqrt(`df'*`r2'/(1-`r2'))
    local pv = 2*ttail(`df', `tv')
    mac shift
  }  
  di %20s "`locus'" %15.0f r(N) %15.3f `tv' %15.4g `pv'
  if `pv'<`minp' {
    local minp `pv'
    local locmin "`locus'"
  }
  if "`saving'"!="" {
    post `post' (1) (1) (`tv') (r(N)) (`pv')
    local vlab "`vlab' `loc' `""`locus'""'"
  }
}
di
di "Smallest p-value is " `minp' " (`locmin')"
return scalar p_min_ty = `minp'
return local locus_p_min_ty = "`locmin'"
quietly count if `use' 
local Nr = r(N)
di
di "Number of records with complete data = " `Nr' _continue
di " (`Nr2' cases and `Nr1' controls)"
if "`impute'"!="" | `nht'<2 {
  di "(remaining analyses will use only these records)"
  local tyu  2
}
else {
  di "Missing values will be imputed using the Stata {bf:impute} command: "/*
    */ _continue
  if "`noisily'"!="" {
    di
  }
  quietly {
    local newxs
    tokenize `htsnps'
    foreach x of local xvars {
      noi di "[`1']" _continue
      mac shift     
      tempvar newx
      local others : subinstr local xvars "`x'" "", word
      `noisily' impute `x' `others' , gen(`newx')
      local newxs "`newxs' `newx'"
    }
    tokenize `newxs'
    foreach x of local xvars {
      replace `x' = `1'
      drop `1'
      mac shift
    }
    replace `use' = 1
    markout `use' `cc' `within' `xvars'
  }
  di
  local tyu 3
}
di 
di %20s "SNP" %15s "N" %15s "t-value" %15s "p-value"

if "`using'"!="" {
  quietly {
    foreach var of local topred {
      local row = `row'+1
      matrix b = beta[`row',1...]
      tempvar newv
      matrix score `newv' = b if `use'
      local xvars "`xvars' `newv'"
    }
    matrix drop beta b
  }
}

local minp 1
local loc 0
local allsnps "`htsnps' `topred'"
tokenize `xvars'
foreach locus of local allsnps {
  local loc = `loc'+1
  quietly {
    corr `cc' `1' if `use'
    local r2 = r(rho)^2 
    local df = r(N)-2 
    local tv = sqrt(`df'*`r2'/(1-`r2'))
    local pv = 2*ttail(`df', `tv')
    mac shift
  }  
  di %20s "`locus'" %15.0f r(N) %15.3f `tv' %15.4g `pv'
  if "`saving'"!="" {
    if "`using'"!="" {
      local ifht = (`loc'<=`nht')
      post `post' (`ifht') (`tyu') (`tv') (r(N)) (`pv')
    }
    else {
      post `post' (1) (`tyu') (`tv') (r(N)) (`pv')
    }
    if `loc'>`nht' { 
      local vlab "`vlab' `loc' `""`locus'""'"
    }
  }
  if `pv'<`minp' {
    local minp `pv'
    local locmin "`locus'"
  }
}
di
di "Smallest p-value (all loci) is " `minp' " (`locmin')"
di
quietly {
  reg `cc' `xvht' if `use'
  local rvar = e(rss)/e(df_r)
  local Fv = e(F)
  local df_n = e(df_m)
  local df_d = e(df_r) - `nstrata' + 1
  local FN = e(N)
  local Fp = Ftail(`df_n', `df_d', `Fv')
}
di "Global test for all `nht' htSNP loci, F(`df_n',`df_d') = " `Fv' _continue
di ", p = " `Fp'
if "`saving'"!="" {
  post `post' (.) (`tyu') (`Fv') (`FN') (`Fp')
  local code = `loc'+1
  local vlab "`vlab' `code' `""_global""'"
}
return scalar p_global = `Fp'
return scalar df_den = `df_d'
return scalar df_num = `df_n'
return scalar F_global = `Fv'
return scalar p_min_all = `minp'
return local locus_p_min_all = "`locmin'"
quietly {
  tempname U V
  tempvar ccd
  summ `cc' if `use'
  gen `ccd' = `cc' - r(mean)
  local yvar = r(Var)
  matrix accum `V' = `xvht' if `use', noconstant deviations
  matrix vecaccum `U' = `ccd' `xvht' if `use', noconstant
  matrix `V' = `yvar'*`V'
  matrix rownames `V' = `htsnps'
  matrix colnames `V' = `htsnps'
  matrix colnames `U' = `htsnps'
}
return matrix score_variance = `V'
return matrix score = `U'
if "`saving'"!="" {
  quietly {
    postclose `post'
    preserve
    use `save', clear
    gen int i = cond(_n<=`nht', _n, _n-`nht')
    label def loci `vlab'
    label values i loci
    decode i, generate(locus)
    label drop loci
    drop i
    if "`transpose'"!="" {
      local lsrt
      local nl = _N
      forvalues i=1/`nl' {
        local li = locus[`i']
        local li `li'
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
    else {
      label def ifhtlab 0 "Imputed SNP" 1 "Typed SNP"
      label values htSNP ifhtlab
      label define uselab 1 "Typed" 2 "Complete records" 3 "Typed + imputed"
      label values use uselab
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
  restore
} 
end 
