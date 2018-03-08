*! Version 1.2, David Clayton 3 Apr, 2003
program define ibdmrg
version 7.0
syntax [varlist (default=none)] using/ /*
   */ [, Ped(varname) Mem(varname) Vars(string)]
/*
  Store varnames for pedigree and two members in ibd file in tp, m1 and m2

  Varnames for pedigree and member in phenotype file are ped and mem
*/
if "`varlist'" == "" {
local varlist "ped mem_1 mem_2"
  confirm variable `varlist'
}
local na : word count `varlist'
if `na' < 3 {
  error 102
}
if `na' > 3 {
  error 103
}
tokenize `varlist'
local tp `1'
local m1 `2'
local m2 `3'
if "`ped'" == "" {
  local ped "ped"
}
if "`mem'" == "" {
  local mem "mem"
}
di
di in whi "Merge IBD and phenotype files"
di in whi "Pedigree is variable " in ye "`tp'" in whi /*
*/ " (in the IBD file) and " in ye "`ped'" in whi " (in the phenotype file)"
di in whi "Members are variables " in ye "`m1'" in whi /*
*/ " and " in ye "`m2'" in whi " (IBD file) and " in ye "`mem'" in whi /*
*/ " (phenotype file)"

/*
  ibd is a scratch file to hold ibd data, phen a file to hold phenotype data

  uped, umem are temporary variable names for pedigree, and member
*/
tempfile ibd phen
tempvar uped umem
/*
  Rename pedigree and first member variables, sort, and store in scratch file
*/
ren `tp' `uped'
ren `m1' `umem'
sort `uped' `umem'
quietly save "`ibd'"
/*
  Read phenotype file. Find variables to be linked
*/
use `using'
confirm variable `ped' `mem'
if "`vars'" == "" {
  unab avars : _all
  tokenize `avars'
  while "`1'" != "" {
    if ("`1'" != "`ped'") & ("`1'" != "`mem'") {
      local vars "`vars' `1'"
    }
    mac shift
  }
}
else {
  confirm variable `vars'
}
di in whi "Variables extracted from phenotype file: " in ye "`vars'"
/*
  Rename pedigree and member variables, sort, and save to scratch file
*/
keep `ped' `mem' `vars'
ren `ped' `uped'
ren `mem' `umem'
sort `uped' `umem'
quietly save "`phen'"
/*
  Add _1 to names of variables to be linked, merge with ibd data
*/
di in whi "Adding variables: "  _continue
tokenize `vars'
while "`1'" != "" {
  ren `1' `1'_1
  di in ye "`1'_1 "  _continue
  mac shift
}
merge `uped' `umem' using "`ibd'"
quietly drop if  _merge == 1
drop _merge
/*
  Restore original name for member1 variable, rename member2 variable,
  sort, and store to scratch file
*/
ren `umem' `m1'
ren `m2' `umem'
sort `uped' `umem'
quietly save "`ibd'", replace
/*
  Reread phenotype file, and add _2 to names of linked variables
*/
use "`phen'"
di in whi "Adding variables: "  _continue
tokenize `vars'
while "`1'" != "" {
  ren `1' `1'_2
  di in ye "`1'_2 "  _continue
  mac shift
}
di
merge `uped' `umem' using "`ibd'"
quietly drop if  _merge == 1
drop _merge
/*
  Restore original variable names for member2 and for pedigree
*/
ren `umem' `m2'
ren `uped' `tp'
end


