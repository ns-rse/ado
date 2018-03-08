*! Version 1.0, DGC Oct 2001
program define _run_grename
  version 7.0

  global GA_grename_post "Allele postfix strings"
  global GA_grename_pf "_1 _2"
  window control static GA_grename_post 10 5 60 10
  window control edit  75 5 25 10 GA_grename_pf
  /* OK/Cancel */
  global GA_grename_OK "exit 3001"
  global GA_grename_cancel "exit 3000"
  window control button "Next >" 60 25 40 20 GA_grename_OK default
  window control button "Cancel" 10 25 40 20 GA_grename_cancel
  * Do it all */
  cap noisily window dialog "Rename loci" /*
     */ 10 10 110 50
  if _rc==3001 {
    local na : word count $GA_grename_pf
    if `na'!=2 { 
      _gamenu_error You must define two (and only two) postfix strings
      _run_grename
    }
    else {
      quietly gloci , post($GA_grename_pf)
      local loci = r(loci)
      foreach locus of local loci {
        _rename `locus' $GA_grename_pf
        if _rc==3002 {
          continue , break
        }
      }
    }
  }
  macro drop GA_grename_*
  end

program define _rename
  args locus pf1 pf2
  global GA_grename_new "New name for `locus'"
  window control static GA_grename_new 10 5 75 10
  global GA_grename_nl
  window control edit 90 5 30 10 GA_grename_nl
  /* OK/Skip */
  global GA_grename_ren "exit 3001"
  global GA_grename_skip "exit 3000"
  global GA_grename_back "exit 3002"
  window control button "Rename" 80 25 40 20 GA_grename_ren default
  window control button "Abort" 10 25 30 20 GA_grename_skip
  window control button "Skip" 45 25 30 20 GA_grename_skip
  * Do it all */
  cap noisily window dialog "Rename locus: `locus'" /*
     */ 20 20 130 50
  if _rc==3001 {
    cap confirm new var $GA_grename_nl`pf1' $GA_grename_nl`pf2'
    if "$GA_grename_nl"!="" & !_rc {
      di "rename `locus'`pf1' $GA_grename_nl`pf1'"
      cap rename `locus'`pf1' $GA_grename_nl`pf1'
      win push rename `locus'`pf1' $GA_grename_nl`pf1'
      if _rc!=0 {
        exit 3002
      }
      di ". rename `locus'`pf2' $GA_grename_nl`pf2'"
      cap rename `locus'`pf2' $GA_grename_nl`pf2'
      win push rename `locus'`pf2' $GA_grename_nl`pf2'
      if _rc!=0 {
        exit 3002
      }
    }
  }
  end 
    
    
