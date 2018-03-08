*! version 6 mh 6/aug/1996
cap program drop blik
program define blik
version 4.0

parse "`*'", parse(" ,")
confirm integer number `1'
confirm integer number `2'
local d = `1' 
mac shift
local h = `1'
mac shift
local options " Cut(real 0.1465) Null(string) Param(string) * "
parse "`*'"
preserve
clear

if "`param'"=="" {
di in re "must specify the parameter"
exit
}

if "`null'"!="" {
  local nval=`null'
}

if "`null'"=="" {
  local yline "yline(`cut')"
  local xline ""
  }
else {
  local yline ""
  local xline "xline(`nval')"
  }



***************************************
*Checks that D and H not both zero
*Exchanges D and H if H zero
***************************************

if `h'==0 & `d'==0 {
 di in re "No data"
 exit
}
if `h'==0 {
 local h = `d'
 local d = 0
 if "`null'"!="" {
  local nval=1-`nval'
 }
di in bl "D and H have been interchanged"
}


********************************
*Headings etc.
********************************


local N=`d'+`h'
local R = 6.0
local length 1001
local M =`d'/`N'
qui set obs `length'




*********************************
*special case where D is zero
*********************************

if `d'==0 {
 if "`param'"=="logodds" {
  di in re "logodds not allowed with zero events"
  exit
  }
  gen prob=(_n-1)*0.5/`length'
  gen double lik = (1-prob)^`h'

  *The likelihood graph
  if "`param'"=="prob" {
   format prob %3.2f
   #delimit ;
   graph lik prob if lik>0.01, symbol(i) connect(s) 
   `yline' `xline' 
   `options'
   l1title("likelihood ratio");
   #delimit cr
  }
  if "`param'"=="odds" {
   gen odds=prob/(1-prob)
   format odds %4.2f
   #delimit ;
   graph lik odds if lik>0.01, symbol(i) connect(s) 
   `yline' `xline'  
   `options'
   l1title("likelihood ratio");
   #delimit cr
  }
  
*The supported range and area
  if "`null'"=="" {
    local M 0
    local low 0
    local high=1-exp(log(`cut')/`h')

   if "`param'"=="prob" {
    di " "
    di in gr "cut-point " in ye `cut'
    di in gr "Most likely value for prob    " in ye %7.5f `M'
    di in gr "Likelihood based limits for prob  " in ye %7.5f /*
    */`low' "  " %7.5f `high'
   }
   if "`param'"=="odds" {
    di " "
    di in gr "cut-point " in ye `cut'
    local M 0
    local low=`low'/(1-`low')
    local high=`high'/(1-`high')
    di in gr "Most likely value for odds  " in ye %7.5f `M'
    di in gr "Likelihood based limits for odds  " /*
    */in ye %7.5f `low' "  " %7.5f `high'
   }
  }

  *The likelihood for the null value

  else {
   if "`param'"=="prob" {
    local M 0
    local lrnull = (1-`nval')^`h'
    di " "
    di in ye `lrnull'
    di in gr "Most likely value for prob  " in ye %7.5f `M'
    di in gr "Null value for prob         " in ye %7.5f `nval'
    di in gr "Lik ratio for null value          " in ye %7.5f `lrnull'
   } 
   if "`param'"=="odds" {
    local M 0
    local lrnull = (1+`null')^(-`h')

    di in gr "Most likely value for odds  " in ye %7.5f `M'
    di in gr "Null value for odds         " in ye %7.5f `nval'
    di in gr "Lik ratio for null value          " in ye %7.5f `lrnull'
   } 
 }
}
*****************************
* The general case
*****************************

else {
  if `M'>0.01 {
   gen prob=(_n-0.5)/`length'
  }
  else {
   gen prob=(_n-0.5)*5*`M'/`length'
  }

  gen double lik =  `d'*log(prob/`M') + `h'*log((1-prob)/(1-`M'))
  qui replace lik=exp(lik)
  if "`param'"=="prob" {
   local M =`d'/`N'
  }
  if "`param'"=="odds" {
   local M =`d'/`h'
  }
  if "`param'"=="logodds" {
   local M = log(`d'/`h')
  }

  if "`param'"=="prob" {
   format prob %3.2f
   #delimit ;
   graph lik prob if lik>0.01, symbol(i) connect(s) 
   `yline' `xline' l1title("likelihood ratio")
   `options' ;
   #delimit cr
  }
  if "`param'"=="odds" {
   gen odds=prob/(1-prob)
   format odds %4.2f
   #delimit ;
   graph lik odds if lik>0.01, symbol(i) connect(s) 
   `yline' `xline' l1title("likelihood ratio")
   `options' ;
   #delimit cr
  }
  if "`param'"=="logodds" {
   gen logodds=log(prob/(1-prob))
   format logodds %4.2f
   #delimit ;
   graph lik logodds if lik>0.01, symbol(i) connect(s) 
   `yline' `xline' l1title("likelihood ratio")
   `options' ;
   #delimit cr
  }



*The supported range

  local M = `d'/`N'
  if "`null'"=="" {
   scalar r1=.00001
   scalar r2 = `M'
   local max=`d'*log(`M')+`h'*log(1-`M')
   scalar f1=`d'*log(r1)+`h'*log(1-r1)-`max'-log(`cut')
   scalar f2=`d'*log(r2)+`h'*log(1-r2)-`max'-log(`cut')
   scalar f= 1
   scalar tol = 0.0001
   while abs(f)>tol {
    scalar r=(r1+r2)*0.5
    scalar f=`d'*log(r)+`h'*log(1-r)-`max' -log(`cut')
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

   scalar r1=.99999
   scalar r2 = `M'
   scalar f1=`d'*log(r1)+`h'*log(1-r1)-`max'-log(`cut')
   scalar f2=`d'*log(r2)+`h'*log(1-r2)-`max'-log(`cut')
   scalar f= 1
   while abs(f)>tol {
    scalar r=(r1+r2)*0.5
    scalar f=`d'*log(r)+`h'*log(1-r)-`max'-log(`cut')
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


   if "`param'"=="prob" {
    di " "
    di in gr "cut-point " in ye `cut'
    di in gr "Most likely value for prob    " in ye %7.5f `M'
    di in gr "Likelihood based limits for prob  " in ye %7.5f /*
    */`low' "  " %7.5f `high'
   }
 
   if "`param'"=="odds" {
    di " "
    di in gr "cut-point " in ye `cut'
    local M=`M'/(1-`M')
    local low=`low'/(1-`low')
    local high=`high'/(1-`high')
    di in gr "Most likely value for odds  " in ye %7.5f `M'
    di in gr "Likelihood based limits for odds  " /*
    */in ye %7.5f `low' "  " %7.5f `high'
   }
   if "`param'"=="logodds" {
    di " "
    di in gr "cut-point " in ye `cut'
    local M=log(`M'/(1-`M'))
    local low=log(`low'/(1-`low'))
    local high=log(`high'/(1-`high'))
    di in gr "Most likely value for logodds  " in ye %7.5f `M'
    di in gr "Likelihood based limits for logodds  " /*
    */in ye %7.5f `low' "  " %7.5f `high'
   }
  }

  *The likelihood for the null value
 
  else {
   if "`param'"=="prob" {
    local M=`d'/`N'
    local lrnull = `d'*log(`nval'/`M') + `h'*log((1-`nval')/(1-`M'))
    local lrnull = exp(`lrnull')
    di in ye `lrnull'
    di in gr "Most likely value for prob  " in ye %7.5f `M'
    di in gr "Null value for prob         " in ye %7.5f `nval'
    di in gr "Lik ratio for null value          " in ye %7.5f `lrnull'
   } 
   if "`param'"=="odds" {
    local M = `d'/`h'
    local lrnull = `d'*log(`nval'/`M') - `N'*log((1+`nval')/(1+`M'))
    local lrnull = exp(`lrnull')
    di in gr "Most likely value for odds  " in ye %7.5f `M'
    di in gr "Null value for odds         " in ye %7.5f `nval'
    di in gr "Lik ratio for null value          " in ye %7.5f `lrnull'
   } 
   if "`param'"=="logodds" {
    local M = log(`d'/`h')
    local lrnull = `d'*(`nval'-`M') - `N'*log((1+exp(`nval'))/(1+exp(`M')))
    local lrnull = exp(`lrnull')
    di in gr "Most likely value for logodds  " in ye %7.5f `M'
    di in gr "Null value for logodds         " in ye %7.5f `nval'
    di in gr "Lik ratio for null value          " in ye %7.5f `lrnull'
   }
 
  *local area = ibeta(`d'+1,`h'+1,`nval')
  *di in gr "Area (%) to the left of null value   " in ye %5.1f 100*`area'
  *di in gr "Area (%) to the right of null value  " in ye %5.1f 100*(1-`area')


  }
}
end

