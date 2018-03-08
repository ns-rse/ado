*! Version 0.1, Sept 2003

program define labgen, rclass

version 8.0

local i 0
while 1 {
  local i = `i' + 1
  local labname L`i'
  capture label def `labname' 0 zero
  local myrc = _rc
  if (`myrc' == 0 ) {
    return local label `labname'
    label drop `labname'
    exit
   }
}

end
