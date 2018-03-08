*! Version 1.0, DGC Sept 2002
program define snp2hap 
syntax varlist(min=4 max=4), [GENerate(string) RANdom ALL(string)]
if "`generate'"=="" {
  local generate _new
}
cap confirm new var `generate'_1 `generate'_2 
if _rc!=0 {
  di in red "Could not create new variables: `generate'_1 & `generate'_2"
  exit
}
if "`all'"!="" {
  if "`random'"!="" {
    di in red "random and all options cannot both be set"
    exit
  }
  cap confirm new var `all'
  if _rc!=0 {
    di in red "Could not create new variable `all'"
    exir
  }
}
tempvar g1 g2
tokenize `varlist'
quietly {
  gen `g1' = `1' + `2' - 2
  replace `g1' = . if `g1'<0 | `g1'>2
  gen `g2' = `3' + `4' - 2
  replace `g2' = . if `g2'<0 | `g2'>2
  gen  `generate'_1=.
  gen  `generate'_2=.
  replace `generate'_1 = 1 if `g1'<2 & `g2'<2
  replace `generate'_2 = 4 if `g1'>0 & `g2'>0
  replace `generate'_1 = 3 if `g1'==2 & `g2'<2         
  replace `generate'_1 = 2 if `g2'==2 & `g1'<2
  replace `generate'_2 = 3 if `g2'==0 & `g1'>0
  replace `generate'_2 = 2 if `g1'==0 & `g2'>0
  replace `generate'_2 = 1 if `g1'==0 & `g2'==0
  replace `generate'_1 = 4 if `g1'==2 & `g2'==2
  replace `generate'_1 = 0 if `g1'==1 & `g2'==1
  replace `generate'_2 = 0 if `g1'==1 & `g2'==1
  replace `generate'_1 = . if `g1'==. | `g2'==.
  replace `generate'_2 = . if `g1'==. | `g2'==.
  count if `generate'_1==1
  local f1 = r(N)
  count if `generate'_2==1
  local f1 = r(N) + `f1'
  count if `generate'_1==2
  local f2 = r(N)
  count if `generate'_2==2
  local f2 = r(N) + `f2'
  count if `generate'_1==3
  local f3 = r(N)
  count if `generate'_2==3
  local f3 = r(N) + `f3'
  count if `generate'_1==4
  local f4 = r(N)
  count if `generate'_2==4
  local f4 = r(N) + `f4'
  count if `generate'_1==0 & `generate'_2==0
  local f0 = r(N)
}
di "Number of subjects genotyped at both loci = " (`f1'+`f2'+`f3'+`f4')/2+`f0'
di "Number with uncertain haplotype phase     = " `f0'
local p 0
local q 0
local on 1
while (`on') {
  local ad = (`f1'+`p'*`f0')*(`f4'+`p'*`f0')
  local bc = (`f2'+`q'*`f0')*(`f3'+`q'*`f0')
  local pl = `p'
  local p = `ad'/(`ad'+`bc')
  local q = 1- `p'
  local on = (`p'-`pl')^2 > 1.0e-6
}
di "Probability 1.1/2.2 rather than 1.2/2.1   = " `p'
if "`random'"!="" {
  tempvar u
  gen `u' = uniform()
  di "Uncertain haplotypes have been assigned at random with this probability"
  quietly {
    replace `generate'_1 = cond(`u'<`p', 1, 2) if `generate'_1 == 0
    replace `generate'_2 = cond(`u'<`p', 4, 3) if `generate'_2 == 0
  }
}
else if "`all'"!="" & `p'!=0 & `p'!=1 {
  di "Uncertain assignments are expanded into two records"
  quietly {
    tempvar id
    gen `id'=_n
    expand 2 if `generate'_1==0 & `generate'_2==0
    sort `id'
    by `id': gen `all' = cond(`generate'_1!=0, 1, cond(_n==1, `p', `q'))
    by `id': replace `generate'_1 = cond(_n==1, 1, 2) if `generate'_1==0
    by `id': replace `generate'_2 = cond(_n==1, 4, 3) if `generate'_2==0
  }
}
else {
  if `p'>0.5 {
    di "Uncertain haplotypes have been assigned 1.1/2.2"
    quietly {
      replace `generate'_1 = 1 if `generate'_1 == 0
      replace `generate'_2 = 4 if `generate'_2 == 0
    }
  }
  else {
    di "Uncertain haplotypes have been assigned 1.2/2.1"
    quietly {
      replace `generate'_1 = 2 if `generate'_1 == 0
      replace `generate'_2 = 3 if `generate'_2 == 0
    }
  }
}
tempname labels
quietly {
  label define `labels' 1 "1.1" 2 "1.2" 3 "2.1" 4 "2.2"
  label values  `generate'_1 `labels'
  label values  `generate'_2 `labels'
}
end
  
