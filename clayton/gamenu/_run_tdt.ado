*! Version 1.0, DGC Oct 2001
* Run tdt from a dialog box
program def _run_tdt
  version 7.0
  syntax [varlist]
  global GA_tdt_vars "`varlist'"
  * Warnings
  global GA_tdt_warn = 1
  window control check "Show warnings" 120 5 60 10 GA_tdt_warn
  * Alleles
  global GA_tdt_selall "Select two alleles"
  window control static GA_tdt_selall 10 5 100 10
  window control msimple GA_tdt_vars 10 15 60 190 GA_tdt_locus
  * Weight
  global GA_tdt_whead "Weight expression (optional)"
  window control static GA_tdt_whead 10 220 100 10
  window control edit 10 230 80 10 GA_tdt_wt  
  * Parental origin
  global GA_tdt_or "Parental Origin"
  global GA_tdt_status "Status code(s)"
  global GA_tdt_opt "(optional)"
  window control static GA_tdt_or 115 28 125 45 blackframe
  window control static GA_tdt_or 150 25 60 12 center
  global GA_tdt_orval 3
  window control radbegin "Mother" 120 35 35 10 GA_tdt_orval
  window control radio    "Father" 160 35 35 10 GA_tdt_orval
  window control radend   "Both"   200 35 35 10 GA_tdt_orval
  window control static GA_tdt_status 120 55 60 10
  window control edit 162 55 40 10 GA_tdt_pa
  window control static GA_tdt_opt 207 55 30 10
  * Offspring
  global GA_tdt_aff "Affected offspring"
  global GA_tdt_ac = 2
  window control static GA_tdt_aff 115 88 125 45 blackframe
  window control static GA_tdt_aff 150 85 60 12 center
  global GA_tdt_sex 3
  window control radbegin "Male" 120 95 35 10 GA_tdt_sex
  window control radio    "Female" 160 95 35 10 GA_tdt_sex
  window control radend   "Both"   200 95 35 10 GA_tdt_sex
  window control static GA_tdt_status 120 115 60 10
  window control edit 162 115 40 10 GA_tdt_ac
  * Robust
  global GA_tdt_robust "Robust variance estimates"
  global GA_tdt_clus = 1
  window control static GA_tdt_robust 115 148 125 25 blackframe
  window control static GA_tdt_robust 140 145 80 12 center
  window control radbegin "No" 120 155 35 10 GA_tdt_clus
  window control radio    "Nuclear" 160 155 35 10 GA_tdt_clus
  window control radend   "Pedigree"   200 155 35 10 GA_tdt_clus
  * Minimum expected frequency
  global GA_tdt_emhead "Minimum Expected Frequency"
  global GA_tdt_emin = 5
  window control static  GA_tdt_emhead 115 195 100 10
  window control edit 230 195 10 10  GA_tdt_emin
  * OK/Cancel
  global GA_tdt_OK "exit 3001"
  global GA_tdt_cancel "exit 3000"
  window control button "OK" 150 220 40 20 GA_tdt_OK default
  window control button "Cancel" 200 220 40 20 GA_tdt_cancel
  *
  cap noisily window dialog "tdt: The Transmission/Disequilibrium Test" /*
     */ 10 10 250 250
  if _rc==3001 {
    local na : word count $GA_tdt_locus
    if `na' != 2 {
      _gamenu_error You must select two (and only two) allele variables
      _run_tdt
    }
    else {
      global GA_tdt_cmd "tdt $GA_tdt_locus"
      if "$GA_tdt_wt" != "" {
        global GA_tdt_cmd "$GA_tdt_cmd [iw=$GA_tdt_wt]$"
      }
      global  GA_tdt_cmd "$GA_tdt_cmd ,$GA_pedvars"
      if !$GA_tdt_warn {
        global  GA_tdt_cmd "$GA_tdt_cmd nowarn"
      }
      if $GA_tdt_orval == 1 {
        global  GA_tdt_cmd "$GA_tdt_cmd morigin"
      }
      if $GA_tdt_orval == 2 {
        global  GA_tdt_cmd "$GA_tdt_cmd porigin"
      }
      if "$GA_tdt_pa" != "" {
        global  GA_tdt_cmd "$GA_tdt_cmd porigin"
      }
      if $GA_tdt_sex == 1 {
        global  GA_tdt_cmd "$GA_tdt_cmd male"
      }
      if $GA_tdt_sex == 2 {
        global  GA_tdt_cmd "$GA_tdt_cmd female"
      }
      global  GA_tdt_cmd "$GA_tdt_cmd acode($GA_tdt_ac)"
      if $GA_tdt_clus==2 {
        global  GA_tdt_cmd "$GA_tdt_cmd robust cluster(nuclear)"
      }
      if $GA_tdt_clus==3 {
        global  GA_tdt_cmd "$GA_tdt_cmd robust cluster(pedigree)"
      }
      global  GA_tdt_cmd "$GA_tdt_cmd emin($GA_tdt_emin)"
      di "$GA_tdt_cmd"
      cap noisily $GA_tdt_cmd
      win push $GA_tdt_cmd
    }
  }
  macro drop GA_tdt_*
  end

