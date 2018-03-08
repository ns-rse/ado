*! version 1.2 Nov 1, 2001 DGC
program define grprare, rclass
  version 7.0
  syntax varlist(min=1 max=2) [, GEN(string) MIN(real 0.05) /*
                 */ COde(string) LABel(string) FOrce]
  local nv: word count `varlist'
  if "`gen'"!="" {
    local ng : word count `gen'
    if `ng'!=`nv' {
      di in red "Number of generated variables does not match number input"
      exit
    }
  }
  if `min'<0 | `min'>100 {
    di in red "Minimum relative frequency must be a proportion or a percentage"
    exit
  }
  if `min'>1 {
    local min = `min'/100
  }
  quietly {
    tempname fr cfr cd cds
    tempvar freq vals last total 
    local add 0
    local fmt
    local typ
    foreach var of local varlist {
      local fmtv : value label `var'
      local typv : type `var'
      if "`fmtv'"!="" {
        if "`fmt'"=="" {
          local fmt "`fmtv'"
        }
        else if "`fmt'"!="`fmtv'" {
          noi di in red "Input variables have different value labels"
          exit
        }
      }
      if "`typ'"=="" {
        local typ "`typv'"
      }
      else if "`typ'"!="`typv'" {
        noi di in red "Input variables have different types"
        exit
      }
      if `add' {
        tabulate `var', matcell(`fr') matrow(`cd')
        matrix `cfr' = `cfr' \ `fr'
        matrix `cds' = `cds' \ `cd'
      }
      else {
        tabulate `var', matcell(`cfr') matrow(`cds')
      }
      local add 1
    }
    capture matrix drop `fr' `cd'
    svmat `cfr', names(`freq')
    svmat `cds', names(`vals')
    sort `vals'
    by `vals': gen `total' = sum(`freq')
    by `vals': gen `last' = (_n==_N)
    replace `vals' = . if !`last'
    replace `total' = . if `vals'==.
    sort `total'
    summarize `total'
    local tfrq = r(sum)
    local grpat = `tfrq'*`min'
    if "`code'"=="" {
      summarize `vals'
      local code = 1 + r(max)
    }
    local vgrp
    replace `freq' = sum(`total')
    local ngrp = 0
    local carryon = 1
    while `carryon' {
      local ngrp = `ngrp' + 1
      local vn = `vals'[`ngrp']
      local carryon = (`freq'[`ngrp']<`grpat') | (`total'[`ngrp']<`grpat')
      if `carryon' {
        local vgrp "`vgrp' `vn'"
      }
    }
    local grpf = `freq'[`ngrp']/`tfrq'
    drop `total' `freq' `vals'
  }
  if `grpf'==1 {
    di "ALL alleles would have to be grouped to achieve target frequency"
    di "Nothing done"
    local carryon 0
  }
  else if `ngrp'==1 {
    di "No grouping of codes necessary"
    if "`gen'"!="" & "`force'"!="" {
      di "New variables will be generated anyway"
      local carryon 1
    }
    else {
      di "Nothing done"
      local carryon 0 
    }
  }
  else {
    di "Codes `vgrp' will be recoded to `code'" _continue
    if "`fmt'"!="" {
      if "`label'"=="" {
        local label "Others"
      }
      di " with label `label'"
      label define `fmt' `code' "`label'", modify
    }
    else {
      di
    }
    di "Relative frequency of the new category = " %7.2f = 100*`grpf' "%"
    return local group `vgrp'
    return scalar new = `code'
    return scalar rf_new = `grpf'
    local carryon 1
  }
  if `carryon' {
    if "`gen'"!="" {
      tokenize `varlist'
      foreach var of local gen {
        quietly {
          gen `typ' `var' = `1'
          label variable `var' "`1' (recoded)"
        }
        mac shift
      }
      local varlist `gen'
    }
    foreach var of local  varlist {
      if "`vgrp'"!="" {
        quietly recode `var' `vgrp' = `code'
      }
      if "`fmt'"!="" {
        quietly label values `var' `fmt'
      }
    }
  } 
end

