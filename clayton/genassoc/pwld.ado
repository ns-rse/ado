!* Version 1.5, David Clayton, Oct 2006
program define pwld
  version 7.0
  syntax varlist(min=2) [fw aw] [if] [in] /*
         */ [, HType GType Dp(integer 2) MATrix(string) noML QTol(real 0.0001) /*
         */  ITer(int 20) MEasure(string) TABle(string)  noLIst MIN(real 0) /*
         */  Pattern(numlist asc) Symbols(string) GRaph COLors(numlist) /*
         */  SAVing(string) REPlace]
  marksample include , novarlist
  if "`htype'"!="" & "`gtype'"!="" {
    di in red "Coding cannot be in both gtype and htype forms"
    exit
  }
  if "`ml'"=="" & "`htype'"!="" {
    di in red "ML estimation not relevant for htype data"
    exit
  }
  local defm "D'"
  local defcol 9 8 5 6 7 4 2 3
  if "`measure'"=="" {
    local measure "`defm'"
  }
  if "`measure'" == "D" | "`measure'" == "Cov"{
    local ldm 1
    local ldmeas "D (Covariance)"
  }
  else if "`measure'" == "D'" {
    local ldm 2
    local ldmeas "Lewontin's D'"
  }
  else if "`measure'" == "Delta" | "`measure'" == "r" {
    local ldm 3
    local ldmeas "Delta (correlation coefficient)"
  }
  else if "`measure'" == "rho" {
    local ldm 4
    local ldmeas "rho"
  }
  else if "`measure'" == "OR" {
    local ldm 5
    local ldmeas "Odds ratio"
  }
  else if "`measure'" == "Q" {
    local ldm 6
    local ldmeas "Yule's Q"
  }
  else if "`measure'" == "delta*" {
    local ldm 7
    local ldmeas "delta* (attributable fraction)"
  }
  else if "`measure'" == "R2" {
    local ldm 8
    local ldmeas "R-squared"
  }
  else {
    di in red "Unrecognized LD measure"
    exit
  }
  if "`pattern'"!="" {
    local ns : word count of `pattern'
    if "`graph'"=="" & "`symbols'"=="" {
      if `ns'>6 {
        di in red "`ns' display symbols should be defined"
        exit
      }
      local symbols ".:+*#@"
    }
  }
  if "`weight'"!="" {
    if "`ml'"=="" & "`weight'"!="fweight" {
      di in red "Only frequency weights are legal with ML estimation"
      exit
    }
    local weight "[`weight'`exp']"
  }
  local nv : word count `varlist'
  if "`gtype'"=="" & "`htype'"=="" & `nv'<4 {
    di in red "Too few variables in varlist"
    exit
  }
  if `ldm'==5 {
    local cld = -1.0
  }
  else {
    local cld = 1.0
  }
  if `min'!=0 {
    di "Loci with minor allele frequencies less than `min'% will be ignored"
    local min = `min'/100
  }
  local lname = 0
  tokenize `varlist'
  local newvlist
  local gnames
  while "`1'" != "" {
    if "`gtype'"!="" { /* Genotype coding */
      quietly count if `1'!=0 & `1'!=1 & `1'!=2 & `1'!=.
      if r(N)!=0 {
        di "Variable {bf:`1'} not coded 0, 1, or 2 --- omitted"
      }
      else {
        local l1 = length("`1'")
        if `l1'>`lname' {
          local lname = `l1'
        }
        local newvlist "`newvlist' `1'"
        local gnames "`gnames' `1'"
      }
      mac shift 
    }
    else if "`htype'"!="" { /* Haplotype coding */
      quietly count if (`1'!=1 & `1'!=2 & `1'!=.) 
      if r(N)!=0 {
        di  "Variable {bf:`1'}  not coded 1 or 2 --- omitted"
      }
      else {
        tempvar newv
        quietly gen `newv' = `1' - 1
        local newvlist "`newvlist' `newv'"
        local l1 = length("`1'")
        if `l1'>`lname' {
          local lname = `l1'
        }
        local gnames "`gnames' `1'"
      }
      mac shift 
    }
    else { /* pair of alleles coding */
      if "`2'" == "" {
        di in red "Uneven number of variables"
        exit
      }
      quietly count if (`1'!=1 & `1'!=2 & `1'!=.)|(`2'!=1 & `2'!=2 & `2'!=.)
      if r(N)!=0 {
        di "Variable {bf:`1'} and/or {bf:`2'} not coded 1 or 2 --- omitted"
      }
      else {
        tempvar newv
        quietly gen `newv' = `1' + `2' - 2
        local newvlist "`newvlist' `newv'"
        local l2 0
        while substr("`1'", `l2', 1)==substr("`2'", `l2', 1) {
          local l2 = `l2' + 1
        }
        local l2 = `l2' - 1
        while (`l2'>0) & (substr("`1'", `l2', 1)=="_") {
          local l2 = `l2' - 1
        }
        if `l2'==0 {
          local gnames "`gnames' `2'"
          local l2 = length("`2'")
        }
        else {
          local gn = substr("`1'", 1, `l2')
          local gnames "`gnames' `gn'"
        }
        if `l2'>`lname' {
          local lname = `l2'
        }
      }
      mac shift
      mac shift
    }
  }

  local nv : word count `newvlist'
  local msz : set matsize
  if "`table'"!="" {
    local nv2 = 2*(`nv'-1)
    if `nv2'>`msz' {
      di "Resetting matrix size to {bf:`nv2'}"
      set matsize `nv2'
    }
    matrix define `table' = J(`nv2', `nv2', 0)
    tokenize `gnames'
    local tcnam
    local trnam
    while "`2'"!="" {
      local tcnam "`tcnam' `1':2 `1':1"
      mac shift
      local trnam "`trnam' `1':2 `1':1"
    }
    matrix rownames `table' = `trnam'
    matrix colnames `table' = `tcnam'
  }
  else {
    if `nv'>`msz' {
      di "Resetting matrix size to {bf:`nv'}"
      set matsize `nv'
    }
  }
  tempname Q
  matrix define `Q' = J(`nv', `nv', `cld')
  matrix rownames `Q' = `gnames'
  matrix colnames `Q' = `gnames'

  /* Allele frequencies */

  tokenize `newvlist'
  local i 0
  while "`1'"!="" {
    local i = `i' + 1
    quietly  summ `1' if `include', meanonly
    if "`htype'"!="" {
      local freq = r(mean)
    }
    else {
      local freq = r(mean)/2
    }
    matrix `Q'[`i',`i'] = `freq'
    mac shift
  }

  /* Pairwise LD measures */
  
  tokenize `newvlist'
  mac shift
  local row 2
  while "`1'"!="" {
    local rvar "`1'"
    tokenize `newvlist'
    local col 1
    while "`1'"!="`rvar'" {
      tempvar use  xy  h2
      quietly {
        gen byte `use' = `include' & (`1'!=.) & (`rvar'!=.)
        gen byte `xy' = `1'*`rvar'
        count if `use'
      }
      if r(N)==0 {
        di "Warning: No valid data for {bf:`1'} versus {bf:`rvar'}"
        matrix def `Q'[`row',`col'] = .
        matrix def `Q'[`col',`row'] = .
      }
      else {
        local fm = min(`Q'[`row',`row'], 1-`Q'[`row',`row'],/*
           */ `Q'[`col',`col'], 1-`Q'[`col',`col'])
        if `fm' > `min' {
          quietly {
            if "`htype'"!="" {
              summ `1' `weight' if `use', meanonly
              local twon = r(sum_w)
              local csum = r(mean)
              summ `rvar' `weight' if `use', meanonly
              local rsum = r(mean)
              summ `xy' `weight' if `use', meanonly
              local c1 = r(mean)
            }
            else {
              summ `1' `weight' if `use', meanonly
              local twon = 2*r(sum_w)
              local csum = r(sum)/`twon'
              summ `rvar' `weight' if `use', meanonly
              local rsum = r(sum)/`twon'
              summ `xy' `weight' if `use', meanonly
              local sxy = r(sum)
              gen byte `h2' = (`1'==1) & (`rvar'==1)
              summ `h2' `weight' if `use', meanonly
              local sh2 = r(sum)
              drop `h2'
              local c1 = `sxy'/`twon' - `rsum'*`csum'
              if `c1'<0 {
                local c1 0
              }
              if `c1'>`rsum' {
                local c1 = `rsum'
              }
              if `c1'>`csum' {
                local c1 = `csum'
              }
              if `c1' < (`rsum'+`csum'-1) {
                local c1 = `rsum'+`csum'-1
              }
              if "`ml'"=="" {
                local q 2
                local cont 1
                local it 0
                while (`it'<`iter' & `cont') {
                  local it = `it' + 1
                  local c2 = `rsum' - `c1'
                  local c3 = `csum' - `c1'
                  local c4 = 1 - `rsum' - `csum' + `c1'
                  local ql `q'
                  local q = (`c1'*`c4'-`c2'*`c3')/(`c1'*`c4'+`c2'*`c3')
                  local c1 = (`sxy' + `q'*`sh2')/(2*`twon')
                  local cont = abs(`q'-`ql') > `qtol'
                  if `it'==`iter' & `cont' {
                    local nam1 : word  `row' of `gnames'
                    local nam2 : word  `col' of `gnames'
                    noi di "Warning: no convergence in `it' iterations "/*
                      */ "for {bf:`nam1'} vs {bf:`nam2'}"
                  }
                }
              }
            }
            local c2 = `rsum' - `c1'
            local c3 = `csum' - `c1'
            local c4 = 1 - `rsum' - `csum' + `c1'
            if "`table'"!="" {
              local row0 = 2*(`row'-1)
              local col0 = 2*`col'
              local row1 = `row0' - 1
              local col1 = `col0' - 1
              matrix `table'[`row1',`col1'] = `c1'
              matrix `table'[`row1',`col0'] = `c2'
              matrix `table'[`row0',`col1'] = `c3'
              matrix `table'[`row0',`col0'] = `c4'
            }
            local e1 = `rsum'*`csum'
            local e2 = `rsum' - `e1'
            local e3 = `csum' - `e1'
            local e4 = 1 - `rsum' - `csum' + `e1'
            local D = `c1' - `e1'
            if `ldm'==1  { /* Basic D measure */
              local q = `D'
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
            else if `ldm'==2 { /* D' measure */
              if `D'>0 {
                local q = `D'/min(`e2',`e3')
              }
              else {
                local q = -`D'/min(`e1',`e4')
              }
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
            else if `ldm'==3 { /* Delta */
              local q = `D'/sqrt(`e1'*`e4')
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
            else if `ldm'==4 { /* delta */
              if `D'>0 {
                local q = `D'/min(`e2',`e3')
              }
              else {
                local q = `D'/min(`e1',`e4')
              }
              matrix `Q'[`row', `col'] = abs(`q')
              matrix `Q'[`col', `row'] = abs(`q')
            }
            else if `ldm'==5 { /* Odds ratio */
              if `c2'>0 & `c3'>0 {
                local q = (`c1'*`c4')/(`c2'*`c3')
              }
              else {
                local q = `cld'
              }
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
            else if `ldm'==6 { /* Yule's Q */
              local q = (`c1'*`c4' - `c2'*`c3')/(`c1'*`c4' + `c2'*`c3')
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
            else if `ldm'==7 { /* Attributable fraction */
              if `D' > 0 {
                matrix `Q'[`row', `col'] = abs(`D')/(`rsum'*(1-`csum'))
                matrix `Q'[`col', `row'] = abs(`D')/((1-`rsum')*`csum')
              }
              else {
                matrix `Q'[`row', `col'] = abs(`D')/(`rsum'*`csum')
                matrix `Q'[`col', `row'] = abs(`D')/((1-`rsum')*(1-`csum'))
              } 
            }
            else if `ldm'==8 { /* R2 */
              local q = (`D')^2/(`e1'*`e4')
              matrix `Q'[`row', `col'] = `q'
              matrix `Q'[`col', `row'] = `q'
            }
          }
        }
        else {
          matrix  `Q'[`row', `col'] =  .
          matrix  `Q'[`col', `row'] = .
        }
      }
      mac shift
      local col = `col'+1
      drop `use' `xy'
    }
    mac shift
    local row = `row'+1
  }
  if "`list'"=="" {
    di
    di "Off-diagonal elements are estimates of `ldmeas'" _continue
    if "`htype'"=="" {
      di " (assuming H-W equilibrium)"
    }
    else {
      di
    }
    di "Diagonal elements are relative frequencies of allele 2"
    local wide = `dp'+3
    di
    matrix list `Q', format(%`wide'.`dp'f) noblank noheader
    if `ldm'==5 {
      di
      di "Infinite odds ratios have been denoted by " %`wide'.`dp'f = `cld'
    }
  }
  if "`symbols'"!="" {
    di
    di "Pattern of `ldmeas' matrix, using [cutpoints] symbols: " 
    local s = 0
    local w = `dp'+2
    tokenize `pattern'
    while "`1'"!="" {
      local s = `s' + 1
      local p = substr("`symbols'", `s', 1)
      di " [" %`w'.`dp'f = `1' "] " %1s = "`p'" _continue
      mac shift
    }
    di
    if `ldm'==5 {
      di "(all odds ratios expressed as >= 1.0)"
    }
    else {
      di "(all measures expressed as positive values)"
    }
    di 
    local pattern "`pattern' ."
    forvalues row = 1/`nv' {
      local vname : word `row' of `gnames'
      di %`lname's = "`vname'" " "  _continue
      forvalues col = 1/`nv'{
        if `row'==`col' {
          di %1s = "\" _continue
        }
        else {
          local q = `Q'[`row',`col']
          if `ldm'==5 {
            if `q'<1 {
              local q = 1/`q'
            }
          }
          else {
            if `q'<0 {
              local q = -`q'
            }
          }
          local s = 0
          tokenize `pattern'
          while `q'>`1' {
            local s = `s'+1
            mac shift
          }
          if `s'==0 {
            local p " "
          }
          else {
            local p = substr("`symbols'", `s', 1)
          }
          di %1s = "`p'" _continue
        }
      }
      di
    }   
  }

  /* Graphical display */
  
  if "`graph'"!="" {
    gph open, xsize(20.0) ysize(14.41)
    local step = 23063/`nv'
    if "`pattern'"=="" {
      if "`colors'"=="" {
        local colors  `defcol' 
      }
      local npen : word count `colors'
      local band = 1/`npen'
      local at 
      local pattern
      forvalues i = 2/`npen' {
        local at = `at' + `band'
        local pattern "`pattern' `at'"
      }
      local pattern "`pattern' ."
      local npat `npen'
    }
    else {
      local pattern "`pattern' ."
      local npat : word count `pattern'
      if "`colors'"=="" {
        forvalues i = 1/`npat' {
          local pen =  8 - `npat'  +`i'
          local pen : word `pen' of `defcol'
          local colors "`colors' `pen'"
        }
      }
      else {
        local npen : word count `colors'
        if `npen'<`npat' {
          di in red "Not enough colors specified"
          exit
        }
      }
      local npen `npat'
    }
      
    gph pen 1
    gph text  500 30250 0 0 `ldmeas' 
    local y 1400
    local yy 2000
    forvalues i = 1/`npen' {
      local pen : word `i' of `colors'
      gph pen `pen'
      gph box `y' 29500 `yy' 31000 4
      gph pen 1
      gph text `yy' 31250 0 -1 (`pen')
      if `i' < `npen' {
        local cut : word `i' of `pattern'
        local yt = `yy' + 1250
        gph text `yt' 30250 0 0 `cut'
      } 

      local y = `y' + 2700
      local yy = `yy' + 2700
    }
    local x  0
    forvalues i=1/`nv' {
      local y  0
      local xx = `x' + `step'
      forvalues j=1/`nv' {
        local yy = `y' + `step'
        if `i'==`j' {
          gph pen 1
          local xc = (`x' + `xx')/2
          local name : word `i' of `gnames'
          if `nv'<=40 {
            gph text `xc' `xc'  0 0 `i'
            gph text `xc' 23200 0 -1 `name'
          }
          else if (`nv'<=200 & mod(`i',5)==1) | mod(`i',10)==1 {
            gph text `xc' 23200 0 -1 `name'
          }
        }
        else {
          local qij = abs(`Q'[`i',`j'])
          if `qij'!=. {
            local k `npat'
            while `k'>0 {
              local high : word `k' of `pattern'
              if `qij'<`high' {
                local pen : word `k' of `colors'
              }
              local k = `k' - 1
            }
            gph pen `pen'
            gph box `y' `x' `yy' `xx' 4
          }
        }
        local y `yy'
      }
      local x = `xx'
    }
  }

  /* Save as file */

  if "`saving'"!="" {
    tempname pfile
    tempfile tfile
    postfile `pfile' row col LD_measure using "`tfile'"
    forvalues i=1/`nv' {
      local i1 = `i'-1
      forvalues j=1/`i1' {
        local q = `Q'[`i',`j']
        if `q'!=. {
          post `pfile' (`i') (`j') (`q')
        }
      }
    }
    postclose `pfile'
    preserve
    use "`tfile'"
    label define L_`pfile' 0 .
    forvalues i=1/`nv' {
      local name : word `i' of `gnames'
      label define L_`pfile' `i' "`name'", add
    }
    label values row L_`pfile'
    label values col L_`pfile'
    save "`saving'", `replace'
  }

  /* Save as matrix */
  
  if "`matrix'" == "" {
    matrix drop `Q'
  }
  else {
    matrix rename `Q' `matrix'
  }
    
  end
