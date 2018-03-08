*! version 1.1 , DGC,  3 Apr 2003

program define htboot, rclass
   version 7.0
   preserve
   syntax [varlist] using [fw aw iw pw/], HTsnps(varlist) [REPS(integer 100) /*
          */ RA(real 0.0) SAVing(string) REPlace CLuster(varname)  /*
          */ DOTS(integer 10) PERCentiles(numlist asc min=1 max=8 >0 <100)]
   if `ra'<0 | `ra'>1 {
     di in red "Illegal ra() option"
   }
   local pqmin = `ra'*(1-`ra')
   tempname beta
   tempvar wt bwt n sn use p first work
   if "`weight'"=="" {
     gen `wt'= 1
   }
   else {
     quietly {
       gen `wt' = `exp'
       count if `wt'<0
       if r(N)>0 {
         di in red "Negative weights not allowed"
        exit
       }
     }
   }
   foreach locus of local htsnps {
     local varlist : subinstr local varlist "`locus'" "", word
   }
   gen byte `use'=1
   markout `use' `varlist' `wt' `htsnps'
   foreach locus of local htsnps {
     quietly count if `use' & `locus'!=1 & `locus'!=2
     if r(N) > 0 {
       di in red "ht-SNP {bf:`locus'} is not coded 1 or 2"
       exit
     }
   }
   local nht : word count `htsnps'
   local vlist
   foreach locus of local varlist {
     quietly count if `use' & `locus'!=1 & `locus'!=2
     if r(N) > 0 {
       di "{bf:`locus'} is not coded 1 or 2 - not treated as a locus"
     }
     else {
       quietly {
         summ `locus' [aw=`wt'], meanonly
         local pq = (r(mean)-1)*(2-r(mean))
       }
       if `pq'<`pqmin' {
         di "Minor allele of {bf:`locus'} is too rare - ignored"
       }
       else {
         local vlist "`vlist' `locus'"
       }
     }
   }
   quietly {
     keep if `use'
     keep `cluster' `htsnps' `vlist' `wt'
   }
   sort `htsnps'
   merge `htsnps' `using'
   local newvlist
   foreach locus of local vlist {
     cap confirm var _hap_`locus'
     if _rc != 0 {
       di "Locus {bf:`locus'} could not be found in {bf:using} dataset"
     }
     else {
       local newvlist "`newvlist' `locus'"
     }
   }
   local vlist "`newvlist'"
   local nv: word count `vlist'
   if `nv'==0 {
     di in red "No loci other than ht-SNPs"
     exit
   }
   tempname beta
   local matsz: set matsize
   if `nv'>`matsz' {
     local matsz = `nv'+ 1
     di "Resetting matrix size to {bf:`matsz'}"
     quietly set matsize `matsz'
   }   
   matrix `beta' = J(`nv', 1+`nht', 0)
   matrix rownames `beta' = `vlist'
   matrix colnames `beta' = `htsnps' _cons
   forvalues i = 1/`nv' {
     local var : word `i' of `vlist'
     quietly reg _all_`var' `htsnps'
     matrix `beta'[`i',1] = e(b)
   }
   drop _all_*
   if "`cluster'"!="" {
     sort `cluster'
     by `cluster': gen `first' = (_n==1)
     quietly count if `first'
     local N = r(N)
   }
   else {
     local N = _N
   }
   local pvars
   foreach var of local vlist {
     local pvars "`pvars' _r2h_`var' _r2c_`var' _r2a_`var'"
   }
   if "`saving'"=="" {
     tempfile pfile
   }
   else {
     local pfile "`saving'"
   }
   local matsz: set matsize
   local npost: word count `pvars'
   if `npost'>`matsz' {
     local matsz = `npost'+ 1
     di "Resetting matrix size to {bf:`matsz'}"
     quietly set matsize `matsz'
   }
   tempname pname
   postfile `pname' `pvars' using "`pfile'", `replace'
   di "Generating {bf:`reps'} bootstrap estimates of R-squared values"
   if "`saving'"!="" {
     di "Results are saved in the the Stata dataset {bf:`saving'}"
   }
   forvalues rep = 1/`reps' {
     /* Generate bootstrap weights ... Bayesian for now */
     quietly {
       if "`cluster'"!="" {
         gen `work' = -log(uniform()) if `first'
         by `cluster': gen `bwt' = sum(`work')
         drop `work'
         summ `bwt' if `first', meanonly
         local sc = `N'/r(sum)
         replace `bwt' = `wt'*`bwt'*`sc'
       }
       else {
         gen `bwt' = -log(uniform())
         summ `bwt', meanonly
         local sc = `N'/r(sum)
         replace `bwt' = `wt'*`bwt'*`sc'
       }
     }
     /* For each variable in turn, calculate the 3 R^2 measures */
     local pline
     forvalues i = 1/`nv' {
       local var : word `i' of `vlist'
       quietly {
         reg `var' _hap_`var' [iw=`bwt']
         local r2 = e(r2)
         local pline "`pline' (`r2')"
         reg `var' _hcl_`var' [iw=`bwt']
         local r2 = e(r2)
         local pline "`pline' (`r2')"
         tempname bi
         matrix `bi' = `beta'[`i',1...]
         mat score `work' = `bi'
         reg `var' `work' [iw=`bwt']
         drop `work'
         local r2 = e(r2)
         local pline "`pline' (`r2')"
       }
     }
     post `pname' `pline' 
     drop `bwt'
     if `dots'>0 & mod(`rep',`dots')==0 {
       di "." _continue
     }
   }
   postclose `pname'
   di
   di "Analysis of bootstrap estimates:"
   di
   clear
   use "`pfile'"
   if "`percentiles'"=="" {
     local percentiles "5 10 25 50 75 90 95"
   }
   local npc : word count `percentiles'
   foreach index in r2h r2c r2a {
     di "R-squared for " _continue
     if "`index'"=="r2h" {
       di "ht-haplotype prediction:"
     }
     else if "`index'"=="r2c" {
       di "ht-haplotype clump prediction:"
     }
     else {
       di "ht-allelic prediction:"
     }
     di
     di "--------------------------------" _continue
     forvalues i = 1/`npc' {
       di %6s = "------" _continue
     }
     di
     di %16s = "Locus" %8s = "Mean" %8s = "SD" _continue
     foreach pct of local percentiles {
       di %6s = "`pct'%" _continue
     }
     di
     di "--------------------------------" _continue
     forvalues i = 1/`npc' {
       di %6s = "------" _continue
     }
     di
     foreach var of local vlist {
       quietly summ _`index'_`var'
       di %16s = "`var'" %8.2f = 100*r(mean) %8.2f = 100*r(sd) _continue
       if `npc'>0 {
         _pctile  _`index'_`var', percentiles(`percentiles')
         forvalues i = 1/`npc' {
           di %6.1f = 100*r(r`i') _continue
         }
       }
       di
     }
     di "--------------------------------"  _continue
     forvalues i = 1/`npc' {
       di %6s = "------" _continue
     }
     di
     di
   }       
   end
