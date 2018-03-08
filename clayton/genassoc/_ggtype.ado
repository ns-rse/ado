*! Version 1.4, Sept 2003

program define _ggtype
    version 7.0
    gettoken type 0 : 0
    gettoken g 0 : 0
    gettoken eqs 0 : 0
    syntax varlist(min=2 max=2) [,MAXLen(int 17) po LABel(string)]
    if "`po'"=="" {
      local sep "/"
    }
    else {
      local sep "|"
    }
    local alen = int((`maxlen'-1)/2)
    tempvar one two a1 a2 gt pres
    tokenize `varlist'
    local labone :value label `1'
    local labtwo :value label `2'
    if "`labone'"!="" & "`labtwo'"!="" & "`labone'"!="`labtwo'" {
      di in red "Input variables have different value labels"
      exit
    }
    else if "`labone'"!="" {
      local alab `labone'
    }
    else if "`labtwo'"!="" {
      local alab `labtwo'
    }
    else {
      local alab
    }
    quietly {
      gen byte `pres' = (`1'==.)|(`2'==.)
      if "`po'"=="" {
        gen `one' = cond(`1'<=`2', `1', `2')
        gen `two' = cond(`1'<=`2', `2', `1')
      }
      else {
        gen `one' = `1'
        gen `two' = `2'
      }
      if "`alab'"=="" {
        gen str`alen' `a1' = string(`one')
      }
      else {
        label val `one' `alab'
        label val `two' `alab'
        decode `one', generate(`a1') maxl(`alen')
      }
      if "`alab'"=="" {
        gen str`alen' `a2' = string(`two')
      }
      else {
        decode `two', generate(`a2') maxl(`alen')
      }
      gen str`maxlen' `gt' = `a1'+"`sep'"+`a2' if !`pres'
      if "`label'"=="" {
        labgen         
        local label `r(label)'
      }
      encode `gt', gen(`g') label(`label')
    }
end

    
    
