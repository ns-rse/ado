*! Version 1.0, DGC Oct 2001
program define _run_grprare
  version 7.0
  syntax [varlist]
  global GA_grprare_vars "`varlist'"
  * Alleles
  global GA_grprare_selall "Select allele(s)"
  global GA_grprare_locus ""
  window control static GA_grprare_selall 10 5 70 10
  window control msimple GA_grprare_vars 10 15 70 140 GA_grprare_locus
  * New variable name
  global GA_grprare_nv "New variable name(s)"
  global GA_grprare_var
  window control static GA_grprare_nv 90 5 70 10
  window control edit 90 15 70 10 GA_grprare_var
  * Force generation of new variables
  global GA_grprare_force 0
  window control check "Force generation" 90 25 70 10 GA_grprare_force right
  * Minimum relative frequency 
  global GA_grprare_fr "Min frequency (%)"
  global GA_grprare_mrf 5
  window control static GA_grprare_fr 90 40 55 10
  window control edit 150 40 10 10 GA_grprare_mrf
  if $GA_grprare_mrf<1 {
    global GA_grprare_mrf = $GA_grprare_mrf/100
  }
  * New code/label
  global GA_grprare_new "New category"
  global GA_grprare_opt "(optional)"
  window control static GA_grprare_new 90 65 70 40 blackframe
  window control static GA_grprare_new 105 55 40 7 center
  window control static GA_grprare_opt 105 62 40 10 center
  global GA_grprare_cd "Code"
  global GA_grprare_lb "Label"
  window control static GA_grprare_cd 95 75 20 10 
  window control static GA_grprare_lb 95 90 20 10 
  global GA_grprare_code
  global GA_grprare_lab
  window control edit 140 75 15 10  GA_grprare_code
  window control edit 115 90 40 10  GA_grprare_lab
  * OK/Cancel
  global GA_grprare_OK "exit 3001"
  global GA_grprare_cancel "exit 3000"
  window control button "OK" 120 110 40 20 GA_grprare_OK default
  window control button "Cancel" 120 135 40 20 GA_grprare_cancel
  * Do it all */
  cap noisily window dialog "grprare: Group rare alleles" /*
     */ 10 10 170 160
  if _rc==3001 {
    local na : word count $GA_grprare_locus
    local nn : word count $GA_grprare_var
    /* if (`na'!=2 | `na'!=1 ) {*/
    if (`na'!=2 & `na'!=1 ) {
      _gamenu_error You must select one or two allele variables
      _run_grprare
    }
    else if (`na'!=`nn') { 
      _gamenu_error Incorrect number of new variables
      _run_grprare
    }
    else {
      global GA_grprare_cmd "grprare $GA_grprare_locus, gen($GA_grprare_var)"
      if $GA_grprare_force {
        global GA_grprare_cmd "$GA_grprare_cmd force"
      }
      if "$GA_grprare_mrf"!=""{
        global GA_grprare_cmd "$GA_grprare_cmd min($GA_grprare_mrf)"
      }
      if "$GA_grprare_code"!=""{
        global GA_grprare_cmd "$GA_grprare_cmd code($GA_grprare_code)"
      }
      if "$GA_grprare_lab"!=""{
        global GA_grprare_cmd "$GA_grprare_cmd label($GA_grprare_lab)"
      }
      di "$GA_grprare_cmd"
      cap noisily $GA_grprare_cmd
      win push $GA_grprare_cmd
    }
  }
  macro drop GA_grprare_*
  end
