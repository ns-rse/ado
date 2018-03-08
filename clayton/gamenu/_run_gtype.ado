*! Version 1.0, DGC Oct 2001
program define _run_gtype
  version 7.0
  syntax [varlist]
  global GA_gtype_vars "`varlist'"
  * Alleles
  global GA_gtype_selall "Select two alleles"
  window control static GA_gtype_selall 10 5 60 10
  window control msimple GA_gtype_vars 10 15 60 120 GA_gtype_locus
  * New variable name
  global GA_gtype_nv "New variable name"
  global GA_gtype_var
  window control static GA_gtype_nv 80 5 60 10
  window control edit 80 15 60 10 GA_gtype_var
  * Parental origin
  global GA_gtype_pp "origin of alleles"
  global GA_gtype_po 0
  window control check  "Preserve parental" 80 35 60 10 GA_gtype_po right
  window control static GA_gtype_pp 88 43 60 10
  * Max label length
  global GA_gtype_ll "Maximum label length"
  global GA_gtype_mll 10
  window control static GA_gtype_ll 80 58 65 10
  window control edit 80 67 15 10 GA_gtype_mll
  * OK/Cancel
  global GA_gtab_OK "exit 3001"
  global GA_gtab_cancel "exit 3000"
  window control button "OK" 100 85 40 20 GA_gtab_OK default
  window control button "Cancel" 100 115 40 20 GA_gtab_cancel
  * Do it all */
  cap noisily window dialog "gtype: Generate genotype variable" /*
     */ 10 10 150 140
  if _rc==3001 {
    local na : word count $GA_gtype_locus
    if `na'!=2 { 
      _gamenu_error You must select two (and only two) allele variables
      _run_gtype
    }
    else {
      global GA_gtype_cmd "egen $GA_gtype_var = gtype($GA_gtype_locus)"
      global GA_gtype_cmd "$GA_gtype_cmd, maxlen($GA_gtype_mll)"
      if $GA_gtype_po {
        global GA_gtype_cmd "$GA_gtype_cmd po"
      }
      di "$GA_gtype_cmd"
      cap noisily $GA_gtype_cmd
      win push $GA_gtype_cmd
    }
  }
  macro drop GA_gtype_*
  end
