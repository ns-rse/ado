*! Version 1.0, Oct 2001, DGC
program define _run_ginsheet
  version 7.0
  cap window fopen GA_ginsh_file "Select file to open" /*
    */ ".dat files|*.dat|.txt files|*.txt|All files|*.*"
  global GA_ginsh_pp 0
  global GA_ginsh_ncv 0
  global GA_ginsh_ex "(Additional variables"
  global GA_ginsh_br ")"
  global GA_ginsh_zm 1
  global GA_ginsh_cl 0
  window control check "Clear" 10 5 25 10 GA_ginsh_cl
  window control check "Recode zero to missing"  40 5 75 10 GA_ginsh_zm
  window control check "Preped format" 120 5 50 10 GA_ginsh_pp
  window control static GA_ginsh_ex 120 17 60 10
  window control edit 180 17 13 10 GA_ginsh_ncv
  window control static GA_ginsh_br 197 17 3 10
  global GA_ginsh_vn "Other options passed to insheet"
  window control static GA_ginsh_vn 10 35 190 20 blackframe
  window control static GA_ginsh_vn 55 32 100 10 center
  global GA_ginsh_opt
  window control edit 20 40 170 10 GA_ginsh_opt
  * OK/Cancel
  global GA_ginsh_OK "exit 3001"
  global GA_ginsh_cancel "exit 3000"
  window control button "OK" 110 60 40 20 GA_ginsh_OK default
  window control button "Cancel" 160 60 40 20 GA_ginsh_cancel
  *
  cap noisily window dialog "$GA_ginsh_file" /*
     */ 10 10 210 85
  if _rc==3001 {
    global GA_ginsh_cmd `"ginsheet using "$GA_ginsh_file","'
    if $GA_ginsh_cl {
      global GA_ginsh_cmd "$GA_ginsh_cmd clear"
    }
    if $GA_ginsh_pp {
      global GA_ginsh_cmd "$GA_ginsh_cmd preped"
    }
    
    if $GA_ginsh_ncv>0 {
      global GA_ginsh_cmd "$GA_ginsh_cmd ncov($GA_ginsh_ncv)"
    }  
    if $GA_ginsh_zm {
      global GA_ginsh_cmd "$GA_ginsh_cmd zmiss"
    }
    if "$GA_ginsh_opt"!="" {
      global GA_ginsh_cmd "$GA_ginsh_cmd $GA_ginsh_opt"
    }
    di `"$GA_ginsh_cmd"'
    cap noisily $GA_ginsh_cmd
    win push $GA_ginsh_cmd
  }
  mac drop GA_ginsh_*
  end
 
