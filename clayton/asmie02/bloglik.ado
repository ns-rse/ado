*! version 6.0 mh 10/8/1996
*! revised 19/9/96 to make approx an option
cap program drop bloglik
program define bloglik
version 4.0

parse "`*'", parse(" ,")
confirm integer number `1'
confirm integer number `2'
local d = `1' 
mac shift
local h = `1'
mac shift
local options "Cut(real -1.353) Param(string) Null(string) Approx *"
parse "`*'"
preserve
clear

*if  (`d' == 0 | `h' == 0) {
*        di in red "zero events not allowed"
*        exit
*}


if "`param'"=="" {
di in re "must specify the parameter"
exit
}

if "`param'"=="odds" {
di in re "odds not allowed - please specify logodds"
exit
}

if "`null'"!="" {
   di in re "null not allowed"
   exit
}

if `cut' >0 {
   di in re "cut must be negative on the log likelihood scale"
   exit
}


local N=`d'+`h'
local R = 5.0
local length 1001
qui set obs `length'

local max = `d'*ln(`d')+`h'*ln(`h')-(`N')*ln(`N')

local yline "yline(`cut')"
  
local M =`d'/`N'
local S=sqrt(`M'*(1-`M')/`N')
scalar tol = 0.0001
* finds exact supported range
if `d' > 0 {
scalar r1=.00001
scalar r2 = `M'
local max=`d'*log(`M')+`h'*log(1-`M')
scalar f1=`d'*log(r1)+`h'*log(1-r1)-`max'-`cut'
scalar f2=`d'*log(r2)+`h'*log(1-r2)-`max'-`cut'
scalar f= 1

while abs(f)>tol {
  scalar r=(r1+r2)*0.5
  scalar f=`d'*log(r)+`h'*log(1-r)-`max' -`cut'
    if f*f1>0 {
      scalar r1=r
      scalar f1=f
    }
    else {
      scalar r2=r
      scalar f2=f
    }
  }
  local low=r
}
else {
local low = 0
local max = 0
}
  scalar r1=.99999
  scalar r2 = `M'
  scalar f1=`d'*log(r1)+`h'*log(1-r1)-`max' -`cut'
  scalar f2=`d'*log(r2)+`h'*log(1-r2)-`max' -`cut'
  scalar f= 1
  while abs(f)>tol {
    scalar r=(r1+r2)*0.5
    scalar f=`d'*log(r)+`h'*log(1-r)-`max'-`cut'
    if f*f1>0 {
      scalar r1=r
      scalar f1=f
    }
    else {
      scalar r2=r
      scalar f2=f
    }
  }
  local high=r

  global S_1 = `low'
  global S_2 = `high'

if "`param'"=="prob" {

* displays exact and approx supported range
di ""
di in gr "cut-point " in ye `cut'
di in gr "Most likely value for prob    " in ye %7.5f `M'
di in gr "Likelihood based limits for prob  " in ye %7.5f /*
*/`low' "  " %7.5f `high'
if "`approx'"!="" {
local plow = `M'- sqrt(-`cut'*2)*`S'
local phigh =  `M'+ sqrt(-`cut'*2)*`S'
di in gr "Approx quadratic limits for prob  " in ye %7.5f /*
*/`plow' "  " %7.5f `phigh'
}

* graphs exact and approx llr

        local start = max(1/`length',`M'-`R'*`S')
        local stop= min(1-1/`length',`M'+`R'*`S')

        qui gen prob =`start' + (`stop'-`start')*(_n-1)/`length' 
        qui gen true =`d'*ln(prob)+`h'*ln(1-prob)-`max'
        if "`approx'"!="" {
          qui gen approx =-0.5*((prob-`M')/`S')^2 
        }

        #delimit ;
        if "`approx'"!="" {;
        graph true approx prob if true>-5&approx>-4, 
	symbol(..) connect(s.)
        `yline' `xline' `options'
        l1title("log likelihood ratio") ;
        };
        else {;
        graph true prob if true>-5, symbol(.) connect(s)
        `yline' `xline' `options'
        l1title("log likelihood ratio") ;
        };
        #delimit cr

}
if "`param'"=="logodds" {
        
* displays exact and approx supported range on logodds scale

        local M =log(`d'/`h')
        local S=sqrt(1/`d'+1/`h')
        local low=log(`low'/(1-`low'))
        local high=log(`high'/(1-`high'))
di  ""
di in gr "cut-point " in ye `cut'
di in gr "Most likely value for logodds    " in ye %7.5f `M'
di in gr "Likelihood based limits for logodds  " in ye %7.5f /*
*/`low' "  " %7.5f `high'
if "`approx'"!="" {
local plow = `M'- sqrt(-`cut'*2)*`S'
local phigh =  `M'+ sqrt(-`cut'*2)*`S'
di in gr "Approx quadratic limits for logodds  " in ye %7.5f /*
*/`plow' "  " %7.5f `phigh'
}
*graphs exact and approx llr on logodds scale

        qui gen logodds =`M' - `R'*`S' + 2*`R'*`S'*(_n-1)/`length' 
        format logodds %4.2f        
        qui gen true =`d'*logodds-`N'*ln(1+exp(logodds))-`max'
        if "`approx'"!="" {
        qui gen approx =-0.5*((logodds-`M')/`S')^2
        }
        di ""
        di in bl "Back on original scale"
        di " "
        local low=exp(`low')/(1+exp(`low'))
        local high=exp(`high')/(1+exp(`high'))
        local M=exp(`M')/(1+exp(`M'))
        di in gr "Most likely value for param    " in ye %7.5f `M'
        di in gr "Likelihood based limits for param  " in ye %7.5f /*
        */`low' "  " %7.5f `high'
        if "`approx'"!="" {
        local plow=exp(`plow')/(1+exp(`plow'))
        local phigh=exp(`phigh')/(1+exp(`phigh'))
        di in gr "Approx quadratic limits for param  " in ye %7.5f /*
        */`plow' "  " %7.5f `phigh'
        }
        #delimit ;
        if "`approx'"!="" {;
        graph true approx logodds if true>-5&approx>-4, symbol(.i) connect(ss)
        `yline' `xline' `options'
        l1title("log likelihood ratio")
        title("llr curve is in yellow, approximation in red or white") ;
        };
        else {;
        graph true logodds if true>-5, symbol(.) connect(s)
        `yline' `xline' `options'
         l1title("log likelihood ratio") ;
        };
        #delimit cr

}

end

