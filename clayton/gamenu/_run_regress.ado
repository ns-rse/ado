*! Version 1.0, DGC Oct 2001
program define _run_regress
  version 7.0
  syntax [varlist]
  global GA_regress_vars "`varlist'"
  * Command
  global GA_regress_rc "Regression command"
  global GA_regress_clist "logistic logit clogit regress"
  window control static GA_regress_rc 10 5 60 10 
  window control ssimple GA_regress_clist 75 5 40 50 GA_regress_cmd
  * Options
  global GA_regress_oo "Options"
  window control static  GA_regress_oo 130 5 25 10
  window control edit 160 5 80 10 GA_regress_opt
  * Conditions
  global GA_regress_co "Conditions"
  window control static  GA_regress_co 120 30 30 10
  window control edit 160 30 80 10 GA_regress_con
  * Depvar
  global GA_regress_rv "Response variable"
  window control static GA_regress_rv 10 75 65 10 center
  window control ssimple GA_regress_vars 10 85 65 100 GA_regress_dep
  * Explanatory variables
  global GA_regress_ev "Explanatory variables"
  window control static  GA_regress_ev 85 65  155 125 blackframe
  window control static  GA_regress_ev 117 62  80 125 center
  global GA_regress_cat "Categorical"
  global GA_regress_met "Metric"
  window control static GA_regress_cat 93 75 65 10 center
  window control static GA_regress_met 165 75 65 10 center
  window control msimple GA_regress_vars 93 85 65 100 GA_regress_cvar
  window control msimple GA_regress_vars 165 85 65 100 GA_regress_mvar
  * OK/Cancel
  global GA_regress_OK "exit 3001"
  global GA_regress_cancel "exit 3000"
  window control button "OK" 140 200 40 20 GA_regress_OK default
  window control button "Cancel" 200 200 40 20 GA_regress_cancel
  * Do it
  cap noisily window dialog "Run a regression command" /*
     */ 10 10 250 225
  if _rc==3001 {
    if "$GA_regress_cmd" == "" {
      _gamenu_error No regression command specified
      _run_regress
    }
    else if "$GA_regress_dep" == "" {
      _gamenu_error No response variable specified
      _run_regress
    }
    else {
      global GA_regress_cmd "$GA_regress_cmd $GA_regress_dep"
      if "$GA_regress_cvar" != "" {
        global GA_regress_cmd "xi:$GA_regress_cmd"
        foreach var of global GA_regress_cvar {
          global GA_regress_cmd "$GA_regress_cmd i.`var'"
        }
      }
      if "$GA_regress_mvar" != "" {
        global GA_regress_cmd "$GA_regress_cmd $GA_regress_mvar"
      }
      if "$GA_regress_opt"!="" {
        global GA_regress_cmd "$GA_regress_cmd, $GA_regress_opt"
      }
      if "$GA_regress_con"!="" {
        global GA_regress_cmd "$GA_regress_cmd $GA_regress_con"
      }
      di "$GA_regress_cmd"
      cap noisily $GA_regress_cmd
      window push $GA_regress_cmd
    }
  }
  macro drop GA_regress_*
  end


