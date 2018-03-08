!+ version 1.1, Oct 2001, DGC

program define greshape
  version 7.0
  syntax [varlist], [gen(string) ID(varname) POSTfix(string)]
  if "`postfix'"=="" {
    local postfix "_1 _2"
  }
  local ac: word count `postfix'
  if `ac'!=2 {
    di in red "There must be 2 allele postfix strings"
    exit
  }
  tokenize `postfix'
  if "`gen'"=="" {
    local gen "_ac`1'`2'"
  }
  confirm new variable `gen'
  cap gloci , postfix(`postfix')
  if _rc {
    di in red "Nothing done"
    exit
  }
  local varies = r(loci)
  /* If no id variable, try and find one */
  if "`id'" == "" {
    foreach var of local varlist {
      quietly inspect `var'
      if r(N_unique)==_N {
        local id "`var'"
        continue , break
      }
    }
  }
  di
  if "`id'"=="" {
    gen _id = _n
    local id "_id"
    di "A variable containing unique identifiers could not be found" _continue
    di ": `id' created"
  }
  else {
    di "Unique subject identifier: `id'"
  }
  di "The dataset will be reshaped so that one allele is contained on each record"
  local nl: word count `varies'
  if `nl'>0 {
    di "The following loci were found: `varies'"
    reshape long `varies', i(`id') j(`gen') string
  }
  else {
    di "No variables found with specified postfix strings; nothing done"
  }
end
    
