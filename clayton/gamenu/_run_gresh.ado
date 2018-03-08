*! Version 1.0, DGC Oct 2001
program define _run_gresh
  version 7.0

  global GA_gresh_post "Allele postfix strings"
  global GA_gresh_pf "_1 _2"
  window control static GA_gresh_post 10 5 60 10
  window control edit  75 5 35 10 GA_gresh_pf
  /* Variable options */
  global GA_gresh_ij "i/j variables (optional)"
  window control static  GA_gresh_ij 10 30 100 40 blackframe
  window control static  GA_gresh_ij 25 25 70 10 center
  global GA_gresh_it "Unique id (existing)"
  global GA_gresh_jt "Allele postfixes (new)"
  window control static GA_gresh_it 15 40 60 10
  window control static GA_gresh_jt 15 55 60 10
  window control edit 80 40 25 10 GA_gresh_i
  window control edit 80 55 25 10 GA_gresh_j
  /* OK/Cancel */
  global GA_gresh_OK "exit 3001"
  global GA_gresh_cancel "exit 3000"
  window control button "OK" 10 80 40 20 GA_gresh_OK default
  window control button "Cancel" 70 80 40 20 GA_gresh_cancel
  * Do it all */
  cap noisily window dialog "greshape: Expand g'type data" /*
     */ 10 10 120 110
  if _rc==3001 {
    local na : word count $GA_gresh_pf
    if `na'!=2 {
      _gamenu_error You must define two (and only two) postfix strings
      _run_gresh
    }
    else {
      global GA_gresh_cmd "greshape , post($GA_gresh_pf)"
      if "$GA_gresh_i"!="" {
        global GA_gresh_cmd "$GA_gresh_cmd id($GA_gresh_i)"
      }
      if "$GA_gresh_j"!="" {
        global GA_gresh_cmd "$GA_gresh_cmd gen($GA_gresh_j)"
      }
      di "$GA_gresh_cmd"
      cap noisily $GA_gresh_cmd
      win push $GA_gresh_cmd
    }
  }
  macro drop GA_gresh_*
  end
