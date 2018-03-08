*! version 1.0, DGC Sept 2001
program define _def_pedvars
  syntax [varlist]
  global GA_defped_vars "`varlist'"
  tokenize `varlist'
  * Pedigree
  global GA_defped_ped "Pedigree"
  global GA_defped_1 "`1'"
  window control static GA_defped_ped 10 5 70 10
  window control ssimple GA_defped_vars 10 15 70 100 GA_defped_1
  * Member
  global GA_defped_mem "Member"
  global GA_defped_2 "`2'"
  window control static GA_defped_mem 10 125 70 10
  window control ssimple GA_defped_vars 10 135 70 100 GA_defped_2 
  * Father
  global GA_defped_pa "Father"
  global GA_defped_3 "`3'"
  window control static GA_defped_pa 90 5 70 10
  window control ssimple GA_defped_vars 90 15 70 100 GA_defped_3
  * Mother
  global GA_defped_ma "Mother"
  global GA_defped_4 "`4'"
  window control static GA_defped_ma 90 125 70 10
  window control ssimple GA_defped_vars 90 135 70 100 GA_defped_4 
  * Pedigree
  global GA_defped_sex "Sex"
  global GA_defped_5 "`5'"
  window control static GA_defped_sex 170 5 70 10
  window control ssimple GA_defped_vars 170 15 70 100 GA_defped_5
  * Member
  global GA_defped_aff "Disease status"
  global GA_defped_6 "`6'"
  window control static GA_defped_aff 170 125 70 10
  window control ssimple GA_defped_vars 170 135 70 100 GA_defped_6
  * OK/Cancel
  global GA_defped_restore "exit 3002"
  global GA_defped_OK "exit 3001"
  global GA_defped_cancel "exit 3000"
  window control button "OK" 80 240 40 20 GA_defped_OK default
  window control button "Restore defaults" 130 240 60 20 GA_defped_restore
  window control button "Cancel" 200 240 40 20 GA_defped_cancel
  * Do it
  cap noisily window dialog "Define pedigree/status variables" /*
     */ 10 10 250 265
  if _rc==3001 {
    global GA_pedvars ""
    if "$GA_defped_1"!="`1'" {
      global GA_pedvars "$GA_pedvars ped($GA_defped_1)"
    }
    if "$GA_defped_2"!="`2'" {
      global GA_pedvars "$GA_pedvars id($GA_defped_2)"
    }
    if "$GA_defped_3"!="`3'" {
      global GA_pedvars "$GA_pedvars fa($GA_defped_3)"
    }
    if "$GA_defped_4"!="`4'" {
      global GA_pedvars "$GA_pedvars mo($GA_defped_4)"
    }
    if "$GA_defped_5"!="`5'" {
      global GA_pedvars "$GA_pedvars sex($GA_defped_2)"
    }
    if "$GA_defped_6"!="`6'" {
      global GA_pedvars "$GA_pedvars aff($GA_defped_2)"
    }
  }
  else if _rc==3002 {
    global GA_pedvars
    _def_pedvars
  }
  mac drop GA_defped*
  end
