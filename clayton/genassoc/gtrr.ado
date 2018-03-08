*! version 1.8, Apr 28, 2003, DGC
program define gtrr, rclass 
  version 7.0
  syntax varlist(min=2 max=2) [, /*
               */ Ped(varname) Id(varname) Mother(varname) Father(varname) /*
               */ SEX(varname) AFfect(varname) ACode(numlist) /*
               */ PO SAVing(string) REPlace REF(string) EMin(real 5) /*
               */ noWARN noTABle ROBust CLuster(string) noANal]
  if "`cluster'"!="" {
    local robust "robust"
  }
  if "`robust'"!="" {
    if "`cluster'"=="" {
      local cl "n"
    }
    else {
      local cl = substr("`cluster'", 1, 1)
    }
    if "`cl'"!="n" & "`cl'"!="p" {
      di in red "Invalid cluster() option"
    }
  }
  if "`anal'"!="" & "`saving'"=="" {
    di in red "Use noanal option only when saving the case-control data"
    exit
  }
  if "`po'"=="" {
    local sep "/"
  }
  else {
    local sep "|"
  }
  preserve
  tokenize  `varlist'
  local v1 "`1'"
  local v2 "`2'"
  local allab : value label `v1'
  /* Put genotype into order --- smallest allele first */
  tempvar swap work error last sv1 sv2 
  quietly {
    gen byte `swap' = `v1'>`v2'
    gen `work' = `v1'
    replace `v1' = `v2' if `swap'
    replace `v2' = `work' if `swap'
    drop `swap' `work'
  }
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
  /* New variable names */
  local cc "case"
  local ma1 "mother_1"
  local ma2 "mother_2"
  local pa1 "father_1"
  local pa2 "father_2"
  local gtype "gt_child"
  local set "set"
  /* Other options */
  if "`ref'"!="" {
    local refopt "ref(`ref')"
  }
  local fmin = `emin'*4
  local fmin = "fmin(`fmin')" 
    
  keep `ped' `id' `mother' `father' `affect' `varlist' 
  tempfile idord
  sort `ped' `id'
  quietly save "`idord'"
  quietly recode `affect' `acode' = 2 * = 1 if `affect'!=. 
  quietly keep if `affect'==2 & `mother'!=. & `father'!=. /*
              */ & `v1'!=. & `v2'!=.
  local naff = _N
  di
  di "Total number of genotyped affected offspring = `naff'" 
  tempvar o1 o2 m1 m2 f1 f2 mid fid tid gtname res
  rename `v1' `o1'
  rename `v2' `o2'
  rename `id' `tid'
  rename `mother' `mid'
  rename `father' `fid'

  rename `mid' `id'
  sort `ped' `id'
  quietly merge `ped' `id' using "`idord'"
  quietly drop if _merge==2
  quietly drop _merge `mother' `father'
  rename `v1' `m1'
  rename `v2' `m2'
  rename `id'  `mid'

  rename `fid' `id'
  sort `ped' `id'
  quietly merge `ped' `id' using "`idord'", update
  quietly drop if _merge==2
  quietly drop _merge `mother' `father'
  rename `v1' `f1'
  rename `v2' `f2'
  rename `id'  `fid'
  rename `tid' `id'
  rename `o1' `v1'
  rename `o2' `v2'
  drop `affect'
  rename `mid' `mother'
  rename `fid' `father'
  quietly drop if `m1'==. | `m2'==. | `f1'==. | `f2'==. | `v1'==. | `v2'==.
  di "Number of trios with complete data = " _N
  quietly drop if `m1'==`m2' & `f1'==`f2'
  di "Number in which at least one parent is heterozygous = "  _N
  if "`po'"!="" {
    quietly drop if `m1'==`f1' & `m2'==`f2' & `v1'!=`v2'
    di "Number in which parental origin of child alleles is known = "  _N
  }
  di "(Any other trios have been dropped)"
  quietly {
    /* pseudo-controls: o1 and o2 contain genotype */ 
    expand 4
    sort `ped' `id'
    by `ped' `id': gen `o1' = cond(_n<3, `m1', `m2')
    by `ped' `id': gen `o2' = cond(_n==1 | _n==3, `f1', `f2')
    rename `m1' `ma1'
    rename `m2' `ma2'
    rename `f1' `pa1'
    rename `f2' `pa2'
    /* Only two controls for intercross parents if PO analysis */
    if "`po'"!="" {
      by `ped' `id': drop if (`ma1'==`pa1')&(`ma2'==`pa2')&(_n==2 | _n==3)
    }
    /* Sort pseudo-control alleles into f1 f2 */
    gen `f1' = cond(`o1'<`o2', `o1', `o2')
    gen `f2' = cond(`o1'<`o2', `o2', `o1')
    /* Case is first one which matches the offspring genotype */
    by `ped' `id': gen `cc' = sum(`v1'==`f1' & `v2'==`f2')
    by `ped' `id': gen `error' = `cc'[_N]==0
    by `ped' `id': gen `last' = (_n==_N)
    count if `error' & `last'
    if r(N)>0 {
      noisily if "`warn'"=="" {
        di
        di "Warning: " r(N) " genotyping errors or ex-paternities :"
        list `ped' `id' if `error' & `last' , noobs
        di
      }
      drop if `error'
    }
    drop `last'
    drop `error'
    by `ped' `id': replace `cc' = (sum(`cc')==1)
    drop `v1' `v2'
    if "`po'"=="" {
      drop `o1' `o2'
      local v1 child_1
      local v2 child_2
      rename `f1' `v1'
      rename `f2' `v2'
    }
    else {
      drop `f1' `f2'
      local v1 child_m
      local v2 child_p
      rename `o1' `v1'
      rename `o2' `v2'
    }
    /* Generate genotype codes from alleles */
    if "`allab'"!="" {
      label values `v1' `allab'
      label values `v2' `allab'
      decode `v1', gen(`sv1') maxlen(5)
      decode `v2', gen(`sv2') maxlen(5)
      gen str12 `gtname' = `sv1'+"`sep'"+`sv2'
      drop `sv1' `sv2'
    }
    else {
      gen str12 `gtname' = string(`v1')+"`sep'"+string(`v2')
    }
    encode `gtname', gen(`gtype')
    drop `gtname'
    /* Generate case-control set number */
    sort `ped' `id'
    by `ped' `id': gen `set' = (_n == 1)
    replace `set' = sum(`set')
  }
  if "`saving'"!="" {
    di
    save `saving', `replace'
  }
  if "`anal'"!="" {
    exit
  }
  di
  if "`po'"!="" {
    di "Genotypes are coded as maternal/paternal allele"
    di
  }
  label def capsco 0 "Control" 1 "Case"
  label values `cc' capsco
  if "`table'"=="" {
    di "Table of genotype distribution in cases and pseudo-controls"
    tabulate `gtype' `cc', col
  }
  di
  di
  di "Conditional logistic analysis"
  quietly tabulate `gtype', gen(_G)
  iv_rm _G*, `refopt' `fmin'
  local refgt "`s(ref)'"
  if "`robust'"!="" {
    if "`cl'"=="p" {
      local cluster "cluster(`ped')"
    }
    else {
      tempvar nuke
      egen `nuke' = group(`ped' `mother' `father')
      local cluster "cluster(`nuke')"
    }
  }
  quietly rclogit `cc' _G*, group(`set') `robust' `cluster'
  if "`robust'"!="" {
    di "Robust variance estimates have been used, " _continue
    di "taking clusters to be " _continue
    if "`cl'"=="p" {
      di "complete pedigrees"
    }
    else {
      di "nuclear families"
    }
    di "(Number of clusters = " `e(Nclus)' ")"
  }
  di
  di "`e(chi2type)' chi-squared test = " %9.2g = e(chi2) " on " _continue 
  di e(df_m) " df" "   (P =" %8.1e = chiprob(e(df_m),  e(chi2)) ")"
  di
  di "Genotype relative risk (GRR) estimates:"
  di 
  grrout
  di
  di "(Reference genotype: `refgt')"
  quietly drop _G*
  return local chi2type "`e(chi2type)'"
  return scalar p_val = chiprob(e(df_m),  e(chi2))
  return scalar df = e(df_m)
  return scalar chi2 = e(chi2)
  end

/* Slave function to remove indicator variables */

program def iv_rm, sclass
version 6.0
  syntax varlist [,REF(string) FMin(real 0)]
  local binvar
  local binfrq 0
  local least
  local minfnd _N
  local most
  local maxfnd 0
  tokenize `varlist'
  local found 0
  while "`1'"!="" {
    local label1 : var lab `1'
    local code1 = substr("`label1'", 2+index("`label1'", "=="), .)
    label var `1' "`code1'"
    quietly count if `1'
    local freq = r(N)
    if `freq'<`fmin' {
      if "`code1'"=="`ref'" {
        di  "Reference genotype insufficiently frequent: using most frequent"
        local ref 
      }
      if "`binvar'"=="" {
        local binvar `1'
        local binfrq `freq'
      }	
      else {
        quietly {
	  replace `binvar' = `binvar' + `1'
          local binfrq = `binfrq'+`freq'
          drop `1'
        }
      }
    }
    else {
      if "`code1'"=="`ref'" {
	quietly drop `1'
        local found 1
      }
      else {
        if `freq'<`minfnd' {
          local minfnd `freq'
          local least `1'
	}
        if `freq'>`maxfnd' {
          local maxfnd `freq'
          local most `1'
        }
      }
    }
    mac shift 
  }
  if "`binvar'"!="" {
    if `binfrq'<`fmin' {
      if "`least'"=="" {
        di in red "There are insufficient data to meet " _continue
        di in red "expected frequency threshold"
        exit
      }
      quietly {
        replace `binvar'=`binvar'+`least'
        local binfrq = `binfrq' + `minfnd'
        drop `least'
      }
    }
    label var `binvar' "Others"
    di "Insufficiently frequent genotypes have been pooled " /*
      */ " and labelled " in blue "Others"
    if `binfrq'>`maxfnd' {
      local most `binvar'
    }
  }
  if "`ref'"=="" {
    local ref : var label `most' 
    drop `most'
  }
  else if !`found' {
    di 
    di "Requested reference genotype not found: using most frequent"
    local ref : var label `most' 
    drop `most'
  }
  sreturn local ref "`ref'"
  end


/* Slave function to write GRRs  */

program define grrout
  version 6.0
  matrix bt=e(b)
  matrix Vt=e(V)
  local varlist: colnames bt
  di _skip(7) "Genotype" _skip(7) "GRR" _skip(7) "(z)" /*
    */ _skip(8) "P-value"  _skip(5) "[95% Conf. Interval]"
  tokenize `varlist'
  local i = 0
  while "`1'"!="" {
    local lab : var lab `1'
    di %15s = "`lab'" _continue
    local i = `i'+1
    local beta = el(bt, 1, `i')
    di %10.3f = exp(`beta') _continue
    local se = sqrt(el(Vt, `i', `i'))
    local z = `beta'/`se'
    di %10.2f = `z' _continue
    di %15.1e = chiprob(1, (`z')^2) _continue
    di %13.3f = exp(`beta'-1.96*`se') _continue
    di %12.3f = exp(`beta'+1.96*`se')
    mac shift
  }
  end

