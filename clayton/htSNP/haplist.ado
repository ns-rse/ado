*! version 1.1 , DGC, 21 November 2001

program define haplist
  version 7.0
  preserve
  syntax [varlist] [fw aw iw pw/], HTsnps(varlist)
  tempvar n sum wt cwt hwt use foll
  if "`weight'"=="" {
    gen `wt'= 1
  }
  else {
    quietly {
      gen `wt' = `exp'
      count if `wt'<0
      if r(N)>0 {
        di in red "Negative weights not allowed"
       exit
      }
    }
  }
  gen byte `use'=1
  markout `use' `varlist' `wt'
  quietly {
    keep if `use'
    sort `varlist'
    by `varlist': gen `cwt' = sum(`wt')
    by `varlist': keep if _n==_N
    summ `cwt' , mean
    local tot = r(sum)
    replace `cwt' = `cwt'/`tot'
    gsort `htsnps' -`cwt'
    by `htsnps': replace `wt' = sum(`cwt')
    by `htsnps': gen `hwt' = `wt'[_N]
    by `htsnps': gen byte `foll'=(_n>1)
    sort `htsnps' `foll'
    foreach locus of local varlist {
      local ht = index("`htsnps' ", "`locus' ")
      if `ht'==0 {
        sort `htsnps' `foll'
        by `htsnps': replace `locus'=. if `foll' & `locus'==`locus'[1]
      }
    gsort -`hwt' `htsnps' -`cwt'
    }
  }
  local last = _N
  di
  di _skip(12) _continue
  foreach locus of local varlist {
    local ht = index("`htsnps' ", "`locus' ")
    if `ht'==0 {
      di " " _continue
    }
    else {
      di "v" _continue
    }
  }
  di
  di "   %   Cum{sf}" _skip(2) _continue
  local nl : word count `varlist'
  foreach i of numlist 1/`nl' {
    di mod(`i',10) _continue
  }
  di
  foreach i of numlist 1/`last'{
    if !`foll'[`i'] {
      di
      di "{bf}" %5.1f = `hwt'[`i']*100 "{sf}" _skip(7) _continue
      foreach locus of local varlist {
        local ht = index("`htsnps' ", "`locus' ")
        if `ht'!=0 {
          di "{bf}" %1.0f = `locus'[`i'] "{sf}" _continue
        }
        else {
          di "*" _continue
        }
      }
      di
    }
    di %5.1f = `cwt'[`i']*100 %5.1f = `wt'[`i']*100 _skip(2) _continue
    foreach locus of local varlist {
      local ht = index("`htsnps' ", "`locus' ")
      if `ht'==0 {
        di %1.0f = `locus'[`i'] _continue
      }
      else if !`foll'[`i'] {
        di "{bf}" %1.0f = `locus'[`i'] "{sf}" _continue
      }
      else {
        di "." _continue
      }
    }
    di
  }
 end



