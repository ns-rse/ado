*! Version 1.2, DGC May 2002
* Run pseudocc from a dialog box
program def _run_pscc
  version 7.0
  syntax [varlist]
  global GA_pscc_vars "`varlist'"
  * File
  global GA_pscc_sav "Save case-control data as"
  window control static GA_pscc_sav 115 5 80 10
  window control edit 200 5 40 10 GA_pscc_file
  global GA_pscc_rep 0
  window control check "Replace" 200 15 40 10 GA_pscc_rep 
  * Warnings
  global GA_pscc_warn = 1
  window control check "Show warnings" 10 220 60 10 GA_pscc_warn
  * Alleles
  global GA_pscc_selall "Select alleles (pairs)"
  window control static GA_pscc_selall 10 5 100 10
  window control msimple GA_pscc_vars 10 15 60 190 GA_pscc_locus
  * Select cases 
  global GA_pscc_cases "Select cases"
  window control static GA_pscc_cases 115 28 125 40 blackframe
  window control static GA_pscc_cases 150 25 60 12 center
  global GA_pscc_status "Status code(s)"
  global GA_pscc_ac 2
  window control static GA_pscc_status 120 35 60 10
  window control edit 182 35 40 10 GA_pscc_ac
  global GA_pscc_fst 0
  window control check "Only use first affected offspring" /*
    */ 120 50 100 10  GA_pscc_fst
  * Offspring
  global GA_pscc_aff "Control selection"
  window control static GA_pscc_aff 115 78 125 70 blackframe
  window control static GA_pscc_aff 150 75 60 12 center
  global GA_pscc_num 3
  global GA_pscc_ph 1
  global GA_pscc_p0 0
  global GA_pscc_ex 0
  window control radbegin "One" 120 85 35 10 GA_pscc_num
  window control radio    "Three" 160 85 35 10 GA_pscc_num
  window control radend   "Maximum"   200 85 35 10 GA_pscc_num
  window control check "Preserve haplotype phase" 120 100 100 10 GA_pscc_ph
  window control check "Preserve parental origin" 120 115 100 10 GA_pscc_po
  window control check "Exchangeable parental genotypes" 120 130 105 10 /*
    */ GA_pscc_ex
  * Parental variables
  global GA_pscc_pg "Include parental genotype"
  window control static GA_pscc_pg 115 158 125 25 blackframe
  window control static GA_pscc_pg 135 155 90 12 center
  global GA_pscc_mg 0
  global GA_pscc_gg 0
  window control check "Mother" 120 165 50 10 GA_pscc_mg
  window control check "Father" 185 165 50 10 GA_pscc_fg
  global GA_pscc_pv "Copy other parental variables"
  window control static GA_pscc_pv 115 193 125 70 blackframe
  window control static GA_pscc_pv 135 190 90 12 center
  global GA_psrr_ma "Mother"
  global GA_psrr_fa "Father"
  window control static GA_psrr_ma 120 200 50 10 center
  window control static GA_psrr_fa 185 200 50 10 center
  window control msimple GA_pscc_vars 120 210 50 45 GA_pscc_mv
  window control msimple GA_pscc_vars 185 210 50 45 GA_pscc_fv
  * OK/Cancel
  global GA_pscc_OK "exit 3001"
  global GA_pscc_cancel "exit 3000"
  window control button "OK" 10 243 40 20 GA_pscc_OK default
  window control button "Cancel" 55 243 40 20 GA_pscc_cancel
  *
  cap noisily window dialog "pseudocc: Generate pseudo-controls" /*
     */ 10 10 250 275
  if _rc==3001 {
    local na : word count $GA_pscc_locus
    local odd = mod(`na', 2)
    if `na' == 0 | `odd' {
      _gamenu_error You must select an even number of allele variables
      _run_pscc
    }
    else if "$GA_pscc_file"=="" {
      _gamenu_error You must specify a file in which to save the created dataset
      _run_pscc
    }
    else {
      global GA_pscc_cmd "pseudocc $GA_pscc_locus, saving($GA_pscc_file)"
      if $GA_pscc_rep {
        global  GA_pscc_cmd "$GA_pscc_cmd replace"
      }
      if "$GA_pedvars"!="" {
        global  GA_pscc_cmd "$GA_pscc_cmd $GA_pedvars"
      }
      if !$GA_pscc_warn {
        global  GA_pscc_cmd "$GA_pscc_cmd nowarn"
      }
      global  GA_pscc_cmd "$GA_pscc_cmd acode($GA_pscc_ac)"
      if $GA_pscc_fst {
        global  GA_pscc_cmd "$GA_pscc_cmd first"
      }
      if $GA_pscc_num==1 {
        global  GA_pscc_cmd "$GA_pscc_cmd one"
      }
      if $GA_pscc_num==2 {
        global  GA_pscc_cmd "$GA_pscc_cmd three"
      }
      if $GA_pscc_ph {
        global  GA_pscc_cmd "$GA_pscc_cmd phase"
      }
      if $GA_pscc_po {
        global  GA_pscc_cmd "$GA_pscc_cmd porigin"
      }
      if $GA_pscc_ex {
        global GA_pscc_cmd "$GA_pscc_cmd exch"
      }

      if $GA_pscc_mg {
        global GA_pscc_cmd "$GA_pscc_cmd mgen"
      }
      if $GA_pscc_fg {
        global GA_pscc_cmd "$GA_pscc_cmd fgen"
      }
      if "$GA_pscc_mv"!="" {
        global  GA_pscc_cmd "$GA_pscc_cmd mvar($GA_pscc_mv)"
      }
      if "$GA_pscc_fv"!="" {
        global  GA_pscc_cmd "$GA_pscc_cmd fvar($GA_pscc_fv)"
      }
      di "$GA_pscc_cmd"
      cap noisily $GA_pscc_cmd
      win push $GA_pscc_cmd
    }
  }
  macro drop GA_pscc_*
  end

