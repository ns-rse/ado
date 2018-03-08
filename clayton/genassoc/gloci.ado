!* version 1.2, Dec 2003, DGC
program define gloci, rclass
  version 7.0
  syntax [varlist], [POSTfix(string) LBLdef(string)]
  if "`postfix'"=="" {
    local postfix "_1 _2"
  }
  local ac: word count `postfix'
  if `ac'!=2 {
    di in red "There must be 2 allele postfix strings"
    exit 1
  }
  tokenize `postfix'
  local pf1 "`1'"
  local ilen = length("`pf1'")
  local pf2 "`2'"
  local varies
  local remains `varlist'
  foreach var of local varlist {
    local end = substr("`var'", -`ilen', .)
    if "`end'"=="`pf1'" {
      local slen = length("`var'") - `ilen'
      local stub = substr("`var'", 1, `slen')
      local remains:subinstr local remains "`stub'`pf2'" "", /*
         */ word count(loc found)
      if `found'>0 {
        local remains:subinstr local remains "`stub'`pf1'" "", word
        local varies "`varies' `stub'"
      }
    }
  }
  di
  local nl : word count `varies'
  if `nl'>0 {
    di "The following loci were found: `varies'"
  }
  else {
    di "No variables found with specified postfix strings"
  }
  local remains `remains'
  if "`remains'"!="" {
    di "The following variables were not matched: `remains'"
  }
  return local loci `varies'
  if "`lbldef'"!="" {
    tokenize `varies'
    local i 0
    label define `lbldef' 0 "", modify
    while "`1'"!="" {
      local i = `i' + 1
      label define `lbldef' `i' "`1'", add
      mac shift
    }
  }
end
