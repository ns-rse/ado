*! Version 1.0, DGC Oct 2001
program define _run_htype
  version 7.0
  syntax [varlist]
  global GA_htype_fst "1st (maternal) h'type"
  global GA_htype_sec "2nd (paternal) h`type"
  window control static  GA_htype_fst 10 5 60 10
  window control static  GA_htype_sec 80 5 60 10
  * New variable name
  global GA_htype_new "New variable name"
  window control static GA_htype_new 10 20 60 10
  window control static GA_htype_new 80 20 60 10
  global GA_htype_nv1
  global GA_htype_nv2
  window control edit 10 30 60 10 GA_htype_nv1
  window control edit 80 30 60 10 GA_htype_nv2
  * Alleles
  global GA_htype_vars "`varlist'"
  global GA_htype_selall "Select alleles"
  window control static GA_htype_selall 10 50 60 10
  window control msimple GA_htype_vars 10 60 60 120 GA_htype_all1
  window control static GA_htype_selall 80 50 60 10
  window control msimple GA_htype_vars 80 60 60 120 GA_htype_all2
  * Allele coding
  global GA_htype_ac "Allele coding"
  window control static GA_htype_ac 10 195 130 20 blackframe
  window control static GA_htype_ac 55 190 40 10 center
  global GA_htype_num  0
  global GA_htype_miss 0
  window control check "Force numeric" 15 200 55 10  GA_htype_num
  window control check "Code missing" 80 200 55 10  GA_htype_miss
  * OK/Cancel
  global GA_htype_OK "exit 3001"
  global GA_htype_cancel "exit 3000"
  window control button "OK" 30 225 40 20 GA_htype_OK default
  window control button "Cancel" 100 225 40 20 GA_htype_cancel
  * Do it all */
  cap noisily window dialog "htype: Generate haplotype variables" /*
     */ 10 10 150 250
  if _rc==3001 {
    local n1 : word count $GA_htype_all1
    local n2 : word count $GA_htype_all2
    if `n1'!=`n2' { 
      _gamenu_error You must select equal numbers of loci for both haplotypes
      _run_htype
    }
    else if `n1'==1 {
      _gamenu_error Only one locus specified - nothing to do
      macro drop GA_htype_*
    }
    else {
      global GA_htype_cmd "egen $GA_htype_nv1 = htype($GA_htype_all1),"
      if $GA_htype_num {
        global GA_htype_cmd "$GA_htype_cmd num"
      }
      if $GA_htype_miss {
        global GA_htype_cmd "$GA_htype_cmd miss"
      }
      di "$GA_htype_cmd"
      $GA_htype_cmd
      win push $GA_htype_cmd
      global GA_htype_cmd "egen $GA_htype_nv2 = htype($GA_htype_all2),"
      global GA_htype_cmd "$GA_htype_cmd codeas($GA_htype_nv1)"
      if $GA_htype_num {
        global GA_htype_cmd "$GA_htype_cmd num"
      }
      if $GA_htype_miss {
        global GA_htype_cmd "$GA_htype_cmd miss"
      }
      di ". $GA_htype_cmd"
      cap noisily $GA_htype_cmd
      win push $GA_htype_cmd
    }
  }
  macro drop GA_htype_*
  end
