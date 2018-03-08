*! version 3.0  dgc/mh Jan 29, 1996
program define tabrate
version 4.0
local varlist "req ex min(1) max(2)"
local if "opt"
local in "opt"
local weight "fweight aweight iweight"

#delimit ;
local options "Exposure(string) Graph Trend SMR Level(integer $S_level) 
Per(real 1000) *";
#delimit cr

parse "`*'"
tempvar touse
mark `touse' `if' `in' [`weight' `exp']
markout `touse' `varlist' `exposure'

parse "`varlist'", parse(" ")
local d `1'
recast float `d'
local level = `level'*0.005 +0.50

preserve

qui keep if `touse'

tempvar rate ci_low ci_high efac mu v e xp W
tempname mu dmean dtot exptot v chitr ptrend het one

if "`2'"!="" {
  local xp `2'
}
else {
  gen `one' = 1
  local xp `one'
}

* checks on pyears or expected numbers
 
if "`exposure'"==""{
        di in re "Needs person-years or expected numbers"
        exit
        }

qui count if `d' != float(0) & `d' != float(1) & `d' != .
    if _result(1)>0 {
        di in bl "WARNING: response `d' not coded 0/1"
    }

if "`weight'" != "" {
	gen `W' `exp'
	if "`weight'" == "aweight" {
		qui summ `W'
		qui replace `W' = `W'/_result(3)
	}
	qui replace `d' = `d'*`W'
	qui replace `exposure' = `exposure'*`W'
}
local stats = ("`weight'" == "" | "`weight'"=="fweight")


qui {
	sort  `xp'
	by  `xp': replace `d' = sum(`d')
	by  `xp': replace `exposure' = sum(`exposure')
	by  `xp': keep if _n==_N
        summ `xp' [aweight = `exposure']
        scalar `exptot' = _result(2)
        scalar `mu' = _result(3)
        scalar `v' = _result(4)
        summ `xp' [fweight = `d']
        scalar `dtot' = _result(2)
        scalar `dmean'= _result(3)

if "`trend'"!=""{
        scalar `v' = `v'/`dtot'
        scalar `chitr' = ( `dmean' - `mu')^2/`v'
        scalar `ptrend' = chiprob(1,`chitr')
}

else {
        gen `e' = `exposure' * `dtot' / `exptot'
        replace `e'=(`d'-`e')*(`d'-`e')/`e'
        replace `e'=sum(`e')
        scalar `het' = `e'[_N] 
}


if "`smr'"==""{

        gen `rate'=`per'* `d'/`exposure'
        gen `efac'  = exp(invnorm(`level')*sqrt(1/`d'))
        gen `ci_low'=`rate'/`efac'
        gen `ci_high'=`rate'*`efac'


        rename `d' _D
        rename `ci_low' ci_low
        rename `ci_high' ci_high
        rename `exposure' _Y
        rename `rate' _rate
        format _Y %8.1f
        format _rate ci_low ci_high %8.3f

                if "`graph'"!=""&"`2'"!=""{
                noi graph _rate `xp', c(l) `options'
                }
	#delimit ;
        noi di in gr _n "Table of cases (D), person-time (Y), and rates per " 
		`per' " units of person-time";
        #delimit cr
        if "`2'"!="" {
	  if `stats' {
            noi list  `xp' _D _Y _rate ci_low ci_high, noobs
	  }
	  else {
            noi list  `xp' _D _Y _rate, noobs
	  }
        }
        else {
	  if `stats' {
            noi list  _D _Y _rate ci_low ci_high, noobs
	  }
	  else {
            noi list  _D _Y _rate, noobs
	  }
        }
	if `stats' {
          if "`2'"!="" {
            if "`trend'"!=""{
                noi di in gr _n "Chi-squared for trend  " /*
                */ in ye %8.2f  = `chitr' " ( 1 df, p = " in ye %6.3f `ptrend' " )"
            }        

            else {
                noi di in gr _n "Chisq test for unequal rates = " /*
                */ in ye %8.2f = `het' " (" _N-1 " df, p = " %6.3f chiprob(_N-1, `het') " )"
            }
          }
	}
	else {
	  noi di in bl _n "No tests or confidence intervals with these weights"
        }
}

else    {
        gen `smr'=100 * `d'/`exposure'
        gen `efac'  = exp(invnorm(`level')*sqrt(1/`d'))
        gen `ci_low'=`smr'/`efac'
        gen `ci_high'=`smr'*`efac'

rename `d' _D
rename `ci_low' ci_low
rename `ci_high' ci_high
rename `exposure' _E
rename `smr' _SMR
format _E %8.1f
format _SMR ci_low ci_high %8.3f

        if "`graph'"!=""&"`2'"!=""{
                noi graph _SMR `xp', c(l) `options'
        }

noi di in gr _n "Table of failures (D), expected failures (E), and SMR's"
if "`2'"!="" {
  if `stats' {
    noi list  `xp' _D _E _SMR ci_low ci_high, noobs
  }
  else {
    noi list  `xp' _D _E _SMR, noobs
  }
}
else {
  if `stats' {
    noi list  _D _E _SMR ci_low ci_high, noobs
  }
  else {
    noi list  _D _E _SMR, noobs
  }
}
  if `stats' {
        if "`2'"!="" {
          if "`trend'"!=""{
                noi di in gr _n(2) "chi-squared for trend  " /*
                */ in ye %8.2f  = `chitr' " ( 1 df, p = " in ye %6.3f `ptrend' " )"
          }        

          else {
                noi di in gr _n "Chisq test for unequal SMRs = " /*
                */ in ye %8.2f = `het' " (" _N-1 " df, p = " %6.3f chiprob(_N-1, `het') " )"
          }

        }
  }
  else {
    noi di in bl _n "No tests or confidence intervals with these weights"
  }
}
}
end
