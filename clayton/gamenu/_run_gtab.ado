*! Version 1.0, DGC Oct 2001
program define _run_gtab
  version 7.0
  syntax [varlist]
  global GA_gtab_vars "`varlist'"
  * Alleles
  global GA_gtab_selall "Select two alleles"
  window control static GA_gtab_selall 10 5 100 10
  window control msimple GA_gtab_vars 10 15 60 120 GA_gtab_locus
  * Generate
  global GA_gtab_iv "Generate indicator variables (optional)"
  global GA_gtab_px "Variable prefix string"
  global GA_gtab_ra "Omitted (ref) allele"
  window control static GA_gtab_iv 105 8 135 80 blackframe
  window control static GA_gtab_iv 118 5 109 10 center
  
  window control static GA_gtab_px 110 25 60 10
  window control edit 180 25 50 10 GA_gtab_gen
  * reference allele
  global GA_gtab_ra "Omitted indicator (defines reference allele)"
  global GA_gtab_ac "Allele code/label"
  window control static GA_gtab_ra 110 45 120 10
  * by code or by label
  global GA_gtab_rt 1
  window control radbegin "None" 110 55 30 10 GA_gtab_rt
  window control radio    "Auto" 140 55 30 10 GA_gtab_rt
  window control radio    "Code" 170 55 30 10 GA_gtab_rt
  window control radend   "Label" 200 55 30 10 GA_gtab_rt
  window control static GA_gtab_ac 110 70 68 10 right
  window control edit 180 70 50 10 GA_gtab_ref
  * if
  global GA_gtab_selif "Select if (optional)"
  window control static GA_gtab_selif 110 100 60 10
  window control edit 170 100 60 10 GA_gtab_if
  * OK/Cancel
  global GA_gtab_OK "exit 3001"
  global GA_gtab_cancel "exit 3000"
  window control button "OK" 150 120 40 20 GA_gtab_OK default
  window control button "Cancel" 200 120 40 20 GA_gtab_cancel
  * Do it all */
  cap noisily window dialog "gtab: Tabulate allele frequencies" /*
     */ 10 10 250 145
  if _rc==3001 {
    local na : word count $GA_gtab_locus
    if `na'!=2 { 
      _gamenu_error You must select two (and only two) allele variables
      _run_gtab
    }
    else {
      global GA_gtab_cmd "gtab $GA_gtab_locus"
      if "$GA_gtab_if"!="" {
        global GA_gtab_cmd "$GA_gtab_cmd if $GA_gtab_if"
      }
      if "$GA_gtab_gen"!="" {
        global GA_gtab_cmd "$GA_gtab_cmd, gen($GA_gtab_gen)"
        if $GA_gtab_rt!=1 {
          if $GA_gtab_rt==2 {
            global GA_gtab_cmd "$GA_gtab_cmd ref(_most)"
          }
          else if $GA_gtab_rt==3 {
            if "$GA_gtab_ref"=="" {
              di in red "You must specify a reference allele"
              exit
            }
            global GA_gtab_cmd "$GA_gtab_cmd ref($GA_gtab_ref)"
          }
          else {
            if "$GA_gtab_ref"=="" {
              di in red "You must specify a reference allele"
              exit
            }
            global GA_gtab_cmd "$GA_gtab_cmd ref(_$GA_gtab_ref)"
          }
        }
      }
      di "$GA_gtab_cmd"
      cap noisily $GA_gtab_cmd
      win push $GA_gtab_cmd
    }
  }
  macro drop GA_gtab_*
  end
