*! Version 2.1, 10 Oct 2003, DGC
program define ginsheet
  version 7.0
  syntax using [, ZMiss PREped NCov(int 0) NLoc(int 0) NId(int 0) /*
     */ LOCi(string) COVar(string) ID(string) * ]
  insheet `using', `options'
  unab vars :*
  local nv : word count `vars'
  local last
  local new
  forvalues i = 1/`nv' {
    local name : word `i' of `vars'
    if "`name'"!="v`i'" {
      local j 1
      local last `name'
      local new "`new' `name'"
    }
    else {
      if "`last'"!="" {
        if `j'==1 {
          local new "`new'_1"
        }
        local j = `j'+ 1
        local new "`new' `last'_`j'"
      }
      else {
        local new "`new' `name'"
      }
    }
  }
  if "`last'"!="" {
    di "Names have been read from the first line of the input file"
    if "`preped'`id'`loci'`covar'"!="" {
      di "preped, id, loci, and covar options are ignored"
    }
    local nto = `nv' - `ncov'
    forvalues i = 1/`nv' {
      local vari : word `i' of `vars'
      local newi : word `i' of `new'
      if "`vari'"!="`newi'" {
        rename `vari' `newi'
      }
      if "`zmiss'"!="" & `i' <= `nto' {
        quietly mvdecode `newi', mv(0)
      }
    }
  }
  else {      
    if "`preped'"!="" {
      if `nid'!=0 | "`id'"!="" {
        di in red "Options nid and id not legal for preped file"
        exit
      }
      local nid 6
      local id family id father mother sex affected
    }
    else {
      if "`id'"!="" {
        local nin : word count `id'
        if `nid'==0 {
          local nid `nin'
        }
        else if `nid'!=`nin' {
          di in red "id and nid options inconsistent"
          exit
        }
      }
      else {
        local id
        forvalues i = 1/`nid' {
          local id "`id'ID`i' "
        }
      }
    }
    if "`covar'"!="" {
      local ncn : word count `covar'
      if `ncov'==0 {
        local ncov `ncn'
      }
      else if `ncov'!=`ncn' {
        di in red "cov and ncov options inconsistent"
        exit
      }
    }
    else if `ncov'>0 {
      local covar
      forvalues i = 1/`ncov' {
        local covar "`covar'X`i' "
      }
    }
    if `nloc'==0 {
      local nv = `nv' - `ncov' - `nid'
      if `nv'<=0 {
        di in red "Insufficient variables in input file"
        exit
      }
      if mod(`nv', 2) {
        di in red "Inconsistent number of variables on input file"
      }
      local nloc = `nv'/2
    }
    if "`loci'"!="" {
      local nln : word count `loci'
      if `nln'!=`nloc' {
        di in red "Inconsitent number of locus names declared in loci option"
        exit
      }
    }
    else {
      local loci
      forvalues i = 1/`nloc' {
        local loci "`loci'L`i' "
      }
    }

    local iv=0
    forvalues i = 1/`nid' {
      local iv = `iv' + 1
      local name : word `i' of `id'
      rename v`iv' `name'
      if "`zmiss'"!="" {
        quietly mvdecode `name', mv(0)
      }
    }
    forvalues i = 1/`nloc' {
      local iv = `iv' + 1
      local name : word `i' of `loci'
      rename v`iv' `name'_1
      local iv = `iv' + 1
      rename v`iv' `name'_2
      if "`zmiss'"!="" {
        quietly mvdecode `name'_1 `name'_2, mv(0)
      }
    }
    if `ncov'>0 {
      forvalues i = 1/`ncov' {
        local iv = `iv' + 1
        local name : word `i' of `covar'
        rename v`iv' `name'
      }
    }
  }
  end
      
    
