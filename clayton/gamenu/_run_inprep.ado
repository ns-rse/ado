*! Version 1.0, DGC Oct 2001
program define _run_inprep
  version 7.0
  global GA_inprep_fn "File:"
  global GA_inprep_clr 0
  window control static  GA_inprep_fn 10 5 15 10
  window control edit 25 5 70 10 GA_inprep_file
  window control check "Clear" 95 5 25 10 GA_inprep_clr
  * Options
  global GA_inprep_is "Insheet options"
  window control static GA_inprep_is 10 20 50 10
  window control edit 60 20 60 10 GA_inprep_opt
  * Loci
  global GA_inprep_OR "OR"
  global GA_inprep_lo "Loci"
  global GA_inprep_ln "Names"
  global GA_inprep_nl "Number of loci"
  global GA_inprep_pf "Allele postfix strings"
  global GA_inprep_post "_1 _2"
  window control static GA_inprep_lo 10 40 110 60 blackframe
  window control static GA_inprep_lo 55 37 20 10 center
  window control static GA_inprep_ln 15 45 25 10
  window control edit 40 45 75 10 GA_inprep_loci
  window control static GA_inprep_OR 15 57 100 10 center
  window control static GA_inprep_nl 15 65 60 10
  window control edit 100 65 15 10 GA_inprep_nloc
  window control static GA_inprep_pf 15 85 70 10
  window control edit 95 85 20 10 GA_inprep_post
  * Covariates
  global GA_inprep_co "Covariates"
  global GA_inprep_cn "Names"
  global GA_inprep_nc "Number of covariates"
  window control static GA_inprep_co 10 110 110 20 blackframe
  window control static GA_inprep_co 50 107 30 10 center
  window control static GA_inprep_cn 15 115 25 10
  window control edit 40 115 75 10 GA_inprep_cov
  global GA_inprep_mess1 "If no covariates, there is no need to"
  global GA_inprep_mess2 "give locus names or number of loci"
  window control static GA_inprep_mess1 10 140 100 10 center
  window control static GA_inprep_mess2 10 150 100 10 center
  * OK/Cancel
  global GA_inprep_OK "exit 3001"
  global GA_inprep_cancel "exit 3000"
  window control button "OK"   10 165 40 20 GA_gtab_OK default
  window control button "Cancel" 80 165 40 20 GA_gtab_cancel
  * Do it all */
  cap noisily window dialog "inprep: Read a preped format file" /*
     */ 10 10 130 190
  if _rc==3001 {
    local nl $GA_inprep_nloc
    if "`nl'"=="" {
      local nl 0
    }
    local nln : word count $GA_inprep_loci
    local ncn : word count $GA_inprep_cov
    if "$GA_inprep_file"=="" {
      _gamenu_error You must specify a file
      _run_inprep
    }
    else if `nln'!=`nl' & !(`nln'==0 | `nl'==0) {
      _gamenu_error Number of loci not equal to number of locus names
      _run_inprep
    }
    else {
      global GA_inprep_cmd "inprep using $GA_inprep_file, "
      if $GA_inprep_clr {
        global GA_inprep_cmd "$GA_inprep_cmd clear"
      }
      if "$GA_inprep_loci"!="" {
        global GA_inprep_cmd "$GA_inprep_cmd loci($GA_inprep_loci)"
      }
      else if "$GA_inprep_nloc"!="" {
        global GA_inprep_cmd "$GA_inprep_cmd nlocus($GA_inprep_nloc)"
      }
      global GA_inprep_cmd "$GA_inprep_cmd post($GA_inprep_post)"
      if "$GA_inprep_cov"!="" {
        global GA_inprep_cmd "$GA_inprep_cmd covar($GA_inprep_cov)"
      }
      if "$GA_inprep_opt"!="" {
        global GA_inprep_cmd "$GA_inprep_cmd $GA_inprep_opt"
      }
      di "$GA_inprep_cmd"
      cap noisily $GA_inprep_cmd
      win push  $GA_inprep_cmd
    }
  }
  mac drop GA_inprep_*
  end

