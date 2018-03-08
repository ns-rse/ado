program define binsim
set more off
clear
set obs `1'
rndbin `1' `2' `3'
local null=log(`2')-log(1-`2')
local i=1
local conf = 0
while `i'< `1' {
local d=xb[`i']
local h=`3' - `d'
bloglik `d' `h', param(logodds) xline(`null') cut(-1.353) xlab xscale(-6,2) ylab title("D= "`d' " H= " `h')
sleep 100
pause
local conf = `conf' + (($S_1 < `2') & (`2' < $S_2)) 
local i=`i' + 1
}
di ""
di "OUTPUT FROM BINSIM"
di ""
display "confidence level is " `conf' / `1'
set more on
end
