*! Version 1.0, DGC Oct 2001
program define _run_rtest
  version 7.0
  * Coefficients
  global GA_rtest_co "Coefficients"
  matrix est = e(b)
  global GA_rtest_par : colnames est
  global GA_rtest_par : subinstr global GA_rtest_par "_cons" "", word
  window control static GA_rtest_co 10 5 50 10 center
  window control msimple GA_rtest_par 10 15 50 80 GA_rtest_test
  * Type
  global GA_rtest_tt "Type of test"
  window control static GA_rtest_tt 70 8 90 60 blackframe
  window control static GA_rtest_tt 90 5 50 60 center
  global GA_rtest_type 1
  window control radbegin "Wald" 80 22 70 10 GA_rtest_type
  window control radend "Likelihood ratio" 80 37 70 10 GA_rtest_type
  global GA_rtest_rm 0
  window control check "(Restore model)"  95 45 55 10 GA_rtest_rm
  * OK/Cancel
  global GA_rtest_OK "exit 3001"
  global GA_rtest_cancel "exit 3000"
  window control button "OK" 70 75 40 20 GA_rtest_OK default
  window control button "Cancel" 120 75  40 20 GA_rtest_cancel
  * Do it
  cap noisily window dialog "Test coefficients in a regression model" /*
     */ 10 10  170 100
  if _rc==3001 {
    if "$GA_rtest_par"=="" {
      _gamenu_error No e-class command has been run or no parameters to test
    }
    else if "$GA_rtest_test"=="" {
      _gamenu_error You must specify one or more coefficients to test
      _run_rtest
    }
    else {
      if $GA_rtest_type==1 {
        local cmd "testparm $GA_rtest_test"
        di "`cmd'"
        cap noisily `cmd'
        window push `cmd'
      }
      else {
        local cmd = e(cmd) + " " + e(depvar)
        local drop "$GA_rtest_par"
        foreach par of global GA_rtest_par {
          local drop : subinstr local drop "`par'" "", word
        }
        di "lrtest ,saving(0)"
        cap noisily lrtest ,saving(0)
        window push lrtest ,saving(0)
        if _rc==0 {
          di ". quietly `cmd' `drop' if e(sample)"
          cap `cmd' `drop' if e(sample)
          window push quietly `cmd' `drop' if e(sample)
          if _rc==0 {
            di ". lrtest"
            cap noisily lrtest
            window push lrtest
          }
          if $GA_rtest_rm {
            di
            di "Refitting previous model"
            di 
            di ". quietly `cmd' $GA_rtest_par"
            cap `cmd' $GA_rtest_par
            window push `cmd' $GA_rtest_par
          }
        
        
    }
  }
  macro drop GA_rtest_* 
  end


