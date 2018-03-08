*! Version 1.0, DGC Oct 2001
program def _run_gtrr
  version 7.0
  syntax [varlist]
  global GA_gtrr_vars "`varlist'"
  * Warnings
  global GA_gtrr_warn = 1
  window control check "Show warnings" 90 5 60 10 GA_gtrr_warn
  * Table
  global GA_gtrr_table = 1
  window control check "Show table" 175 5 40 10 GA_gtrr_table
  * Alleles
  global GA_gtrr_selall "Select two alleles"
  window control static GA_gtrr_selall 10 5 80 10
  window control msimple GA_gtrr_vars 10 15 60 190 GA_gtrr_locus
  *Status
  global GA_gtrr_ast "Status code(s) for affected offspring"
  window control static GA_gtrr_ast 90 20 100 10
  global GA_gtrr_ac 2
  window control edit 195 20 20 10 GA_gtrr_ac
  * Parental Origin
  window control check "Preserve parental origin of alleles" /*
     */ 90 35 120 10 GA_gtrr_po
  * Reference allele
  global GA_gtrr_ra "Reference allele for relative risks"
  window control static GA_gtrr_ra 90 50 100 10
  window control edit 195 50 20 10 GA_gtrr_ref
  global GA_gtrr_opt "(Optional)"
  window control static GA_gtrr_opt 90 57 100 10 center
  * Minimum expected frequency
  global GA_gtrr_emhead "Minimum Expected Frequency"
  global GA_gtrr_emin = 5
  window control static  GA_gtrr_emhead 90 70 100 10
  window control edit 205 70 10 10  GA_gtrr_emin
  * Case/pseudo-control data
  global GA_gtrr_cc "Case/pseudo-control data"
  window control static GA_gtrr_cc 90 88 125 35 blackframe
  window control static GA_gtrr_cc 112 85 80 12 center
  global GA_gtrr_sv "Save C/C dataset as"
  window control static GA_gtrr_sv 95 97 60 10
  window control edit 160 97 45 10 GA_gtrr_file
  global GA_gtrr_rep 0
  window control check "Replace" 173 108 32 10 GA_gtrr_rep
  global GA_gtrr_noan 0
  window control check "Suppress analysis" 95 108 60 10 GA_gtrr_noan
  window control static GA_gtrr_opt 132 118 40 12 center 
  * Robust
  global GA_gtrr_robust "Robust variance estimates"
  global GA_gtrr_clus = 1
  window control static GA_gtrr_robust 90 138 125 25 blackframe
  window control static GA_gtrr_robust 115 135 80 12 center
  window control radbegin "No" 95 145 35 10 GA_gtrr_clus
  window control radio    "Nuclear" 135 145 35 10 GA_gtrr_clus
  window control radend   "Pedigree"   175 145 35 10 GA_gtrr_clus
  * OK/Cancel
  global GA_gtrr_OK "exit 3001"
  global GA_gtrr_cancel "exit 3000"
  window control button "OK" 125 180 40 20 GA_gtrr_OK default
  window control button "Cancel" 175 180 40 20 GA_gtrr_cancel
  *
  cap noisily window dialog "gtrr: Family-based genotype relative risks" /*
     */ 10 10 225 210
  if _rc==3001 {
    local na : word count $GA_gtrr_locus
    if `na' != 2 {
      _gamenu_error You must select two (and only two) allele variables
      _run_gtrr
    }
    else {
      global GA_gtrr_cmd "gtrr $GA_gtrr_locus,"
      if "$GA_pedvars"!="" {
        global  GA_gtrr_cmd "$GA_gtrr_cmd $GA_pedvars"
      }
      global  GA_gtrr_cmd "$GA_gtrr_cmd acode($GA_gtrr_ac)"
      if $GA_gtrr_clus==2 {
        global  GA_gtrr_cmd "$GA_gtrr_cmd robust cluster(nuclear)"
      }
      if $GA_gtrr_clus==3 {
        global  GA_gtrr_cmd "$GA_gtrr_cmd robust cluster(pedigree)"
      }
      global  GA_gtrr_cmd "$GA_gtrr_cmd emin($GA_gtrr_emin)"
      if "$GA_gtrr_file"!="" {
        global  GA_gtrr_cmd "$GA_gtrr_cmd saving($GA_gtrr_file)"
        if $GA_gtrr_rep {
          global  GA_gtrr_cmd "$GA_gtrr_cmd replace"
        }
      }
      if !$GA_gtrr_warn {
        global  GA_gtrr_cmd "$GA_gtrr_cmd nowarn"
      }
      if !$GA_gtrr_table {
        global  GA_gtrr_cmd "$GA_gtrr_cmd notab"
      }
      if $GA_gtrr_noan {
        global  GA_gtrr_cmd "$GA_gtrr_cmd noan"
      }
      if $GA_gtrr_po {
        global GA_gtrr_cmd "$GA_gtrr_cmd po"
      }
      if  "$GA_gtrr_ref" != "" {
        global GA_gtrr_cmd "$GA_gtrr_cmd ref($GA_gtrr_ref)"
      }
      di "$GA_gtrr_cmd"
      cap noisily $GA_gtrr_cmd
      win push $GA_gtrr_cmd
    }
  }
  macro drop GA_gtrr_*
  end

