*! Version 1.3, Sept 15 2003

program define _ghtype
    version 7.0
    gettoken type 0 : 0
    gettoken h 0 : 0
    gettoken eqs 0 : 0
    syntax varlist(min=2) [,COdeas(varname) MAXLen(int 80) NUM MISS]
    tempvar pres
    gen byte `pres' = 1
    if "`miss'" == "" {
      markout `pres' `varlist'
    }
    local alen 80
    if "`codeas'"!="" {
      local oldlab : value label `codeas'
      if "`oldlab'"=="" {
        di in red "Variable `codeas' does not have value labels"
        exit
      }
    }
    else {
      local oldlab 
    }
    tempvar hap
    tokenize `varlist'
    quietly {
      gen str`maxlen' `hap' = ""
      while "`1'" != "" {
        if "`codeas'"=="" {
          labgen
          local oldlab `r(label)' 
/*          if "`oldlab'"=="" {
            local oldlab "`1'"
          }
          else {
            local oldlab "`oldlab'_`1'"
          }         
JMMH edit */
        }
        tempvar loc
        local alab : value label `1'
        if ("`num'"!="") | ("`alab'"=="") {
          gen str`alen' `loc' = string(`1')
        }
        else {
          decode `1', generate(`loc') maxl(`alen')
        }
        if "`miss'" != "" {
          replace `loc' = "?" if `loc'=="" | `loc'=="."
        }
        replace `hap' = `hap'+"`sep'"+`loc' if `pres'
        local sep "."
        drop `loc'
        mac shift
      }

      encode `hap', gen(`h') label(`oldlab')
    }    
end
