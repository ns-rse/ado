*! version 1.3 , DGC,  24 Oct 2003

program define haptag, rclass
   version 7.0
   preserve
   syntax [varlist] [using] [fw aw iw pw/], HTsnps(varlist) /*
          */ [RA(real 0.0) SAVing(string) REPlace]
   if `ra'<0 | `ra'>1 {
     di in red "Illegal ra() option"
   }
   if "`using'"!="" & "`saving'"!="" {
     di in red "using ... and saving() cannot be used simultaneously"
     exit
   }
   local pqmin = `ra'*(1-`ra')
   tempname beta
   tempvar wt n sn use p var div pbar qbar res
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
   local vlist
   foreach locus of local varlist {
     quietly count if `use' & `locus'!=1 & `locus'!=2
     if r(N) > 0 {
       di "{bf:`locus'} is not coded 1 or 2 - ignored"
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
     keep `htsnps' `vlist' `wt'
   }
   sort `htsnps'
   quietly {
     by `htsnps': gen `n' = sum(`wt')
     foreach locus of local vlist {
       by `htsnps': replace `locus' = sum(`wt'*(`locus'-1))
     }
     by `htsnps': drop if _n!=_N
     local nhap =  _N
   }
   if "`using'"!="" {
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
   }
   local nv: word count `vlist'
   if `nv'==0 {
     di in red "No loci other than ht-SNPs"
     exit
   }
   quietly {
     summ `n', meanonly
     local N = r(sum)
     local MAF
     local PDE
     local R2
     local R2CL
     local R2A
     local mnpde = .
     local mnr2 = .
     local mnr2cl = .
     local mnr2a = .
     local mxpde = -0.1
     local mxr2 = -0.1
     local mxr2cl = -0.1
     local mxr2a = -0.1
     local swgp = 0
     local swcl = 0
     local sarl = 0
     local stot = 0
     local sdwg = 0
     local stdv = 0
     foreach locus of local vlist {
       summ `locus', meanonly
       local P = r(sum)/`N'
       local Q = 1-`P'
       if `P'< `Q' {
         local MAF "`MAF' `P'"
       }
       else {
         local MAF "`MAF' `Q'"
       }
       gen `p' = `locus'/`n'
       gen `var' = `n'*`p'*(1-`p')
       gen `div' = `n'*`var'
       /* PDE */
       summ `div', meanonly
       local dwg = r(sum)
       local tdv = `N'*`N'*`P'*`Q'
       local pde = 1 - `dwg'/`tdv'
       /* Within haplotype variance */
       summ `var' , meanonly
       local whss = r(sum)
       local tot = `N'*`P'*`Q'
       if "`using'" !="" {
         /* Haplotype r^2 */
         reg `p' _hap_`locus' [iw=`n']
         local  wgp = `whss' + e(rss)
         local r2 = 1 - `wgp'/`tot'
         /* Clump r^2 */
         reg `p' _hcl_`locus' [iw=`n']
         local  wcl = `whss' + e(rss)
         local r2cl = 1 - `wcl'/`tot'
         /* Allelic r^2 */
         reg _all_`locus' `htsnps'
         matrix `beta' = e(b)
         matrix score `pbar' = `beta' 
         reg `p' `pbar' [iw=`n']
         drop `pbar'
         local  arl = `whss' + e(rss)
         local r2a = 1 - `arl'/`tot'         
       }
       else {
         /* Haplotype r^2 */
         local wgp = `whss'
         local r2 = 1 - `wgp'/`tot'
         /* Clump r^2 */
         sort `p'
         summ `locus' , meanonly
         local suml = r(sum)
         gen `sn' = sum(`n')
         gen `pbar' = sum(`locus')/`sn'
         gen `qbar' = (`suml'-sum(`locus'))/(`N' - `sn')
         gen `res' = `sn'*`pbar'*(1-`pbar') + (`N' - `sn')*`qbar'*(1-`qbar')
         summ `res', meanonly
         local wcl = r(min)
         local r2cl = 1 - `wcl'/`tot'
         if "`saving'"!="" {
           tempvar clump
           gen `clump' = !((_n == 1) | (`res' < `res'[_n-1]))
         }
         /* Allelic r^2 */
         reg `p' `htsnps' [iw=`n']
         if "`saving'"!="" {
           matrix `beta' = e(b)
           tempvar allele
           matrix score `allele' = `beta'
         }
         local arl = `whss' + e(rss)
         local r2a = 1 - `arl'/`tot'
         drop `pbar' `qbar' `res' `sn'
       }
       local stdv = `stdv' + `tdv'
       local stot = `stot' + `tot'
       local sdwg = `sdwg' + `dwg'
       local swgp = `swgp' + `wgp'
       local swcl = `swcl' + `wcl'
       local sarl = `sarl' + `arl'
       local PDE "`PDE' `pde'"
       local R2 "`R2' `r2'"
       local R2CL "`R2CL' `r2cl'"
       local R2A "`R2A' `r2a'"
       if `pde'!=. & `pde'<`mnpde' {
         local mnpde `pde'
         local lnpde `locus'
       }
       if `pde'!=. & `pde'>`mxpde' {
         local mxpde `pde'
         local lxpde `locus'
       }
       if `r2'!=. & `r2'<`mnr2' {
         local mnr2 `r2'
         local lnr2 `locus'
       }
       if `r2'!=. & `r2'>`mxr2' {
         local mxr2 `r2'
         local lxr2 `locus'
       }
       if `r2a'!=. & `r2a'<`mnr2a' {
         local mnr2a `r2a'
         local lnr2a `locus'
       }       
       if `r2a'!=. & `r2a'>`mxr2a' {
         local mxr2a `r2a'
         local lxr2a `locus'
       }       
       if `r2cl'!=. & `r2cl'<`mnr2cl' {
         local mnr2cl `r2cl'
         local lnr2cl `locus'
       }
       if `r2cl'!=. & `r2cl'>`mxr2cl' {
         local mxr2cl `r2cl'
         local lxr2cl `locus'
       }
       if "`saving'"!="" {
         rename `allele' _all_`locus'
         rename `clump' _hcl_`locus'
         rename `p' _hap_`locus'
         drop `var' `div'
       }
       else {
         drop `p' `var' `div'
       }
     }
   }
   matrix define all = J(`nv',5,0)
   matrix rownames all = `vlist'
   matrix colnames all = MAF PDE R2 R2CL R2A
   di
   di "ht-SNPs: {bf:`htsnps'}"
   di "Additional SNPs:{bf:`vlist'}"
   di "Number of htSNP haplotypes: {bf:`nhap'}"
   if "`using'"!="" {
     di "R-squared values calculated using previously stored results"
   }
   else {
     di "R-squared values are maximum achievable for these data"
   }
   di
   di "--------------------------------------------------------------------------------"
   di _skip(15) %13s = "Min Allele" %13s = "PDE"  _continue
   di %26s = "R-squared (%)"
   di %15s = "Locus" %13s = "Freq (%)" %13s = "(%)" _continue
   di %13s = "Haplotype" %13s = "Clump" %13s = "Allelic"
   di "--------------------------------------------------------------------------------"
   forvalues i = 1/`nv' {
     local locus: word `i' of `vlist'
     local maf: word `i' of `MAF'
     local pde: word `i' of `PDE'
     local r2: word `i' of `R2'
     local r2cl: word `i' of `R2CL'
     local r2a: word `i' of `R2A'
     di %15s = "`locus'" %13.3f = 100*`maf' %13.3f = 100*`pde' _continue
     di %13.3f = 100*`r2'  %13.3f = 100*`r2cl' %13.3f = 100*`r2a'
     matrix all[`i',1] = `maf'
     matrix all[`i',2] = `pde'
     matrix all[`i',3] = `r2'
     matrix all[`i',4] = `r2cl'
     matrix all[`i',5] = `r2a'
   }
   di "--------------------------------------------------------------------------------"
   di _skip(13) %15s = "Mean" %13.3f = 100*(1 - `sdwg'/`stdv') _continue
   di %13.3f =  100*(1 - `swgp'/`stot') _continue
   di %13.3f =  100*(1 - `swcl'/`stot') _continue
   di %13.3f =  100*(1 - `sarl'/`stot')
   di _skip(13) %15s = "Minimum"  %13.3f = 100*`mnpde' _continue
   di %13.3f =  100*`mnr2' _continue
   di %13.3f =  100*`mnr2cl' _continue
   di %13.3f =  100*`mnr2a'
   di _skip(13) %15s = "Locus"  %13s = "`lnpde'" %13s = "`lnr2'" _continue
   di %13s = "`lnr2cl'" %13s ="`lnr2a'"
   di _skip(13) %15s = "Maximum"  %13.3f = 100*`mxpde' _continue
   di %13.3f =  100*`mxr2' _continue
   di %13.3f =  100*`mxr2cl' _continue
   di %13.3f =  100*`mxr2a'
   di _skip(13) %15s = "Locus"  %13s = "`lxpde'" %13s = "`lxr2'" _continue
   di %13s = "`lxr2cl'" %13s ="`lxr2a'"
   di "--------------------------------------------------------------------------------"
   if "`saving'"!="" {
     di
     di "Saving prediction scores"
     keep `htsnps' _hap_* _hcl_* _all_*
     sort `htsnps'
     save `saving', `replace'
   }
   return scalar pde_mean = 1 - `sdwg'/`stdv'
   return scalar pde_min = `mnpde'
   return local  pde_worst "`lnpde'"
   return scalar pde_max = `mxpde'
   return local  pde_best "`lxpde'"
   return scalar r2_mean = 1 - `swgp'/`stot'
   return scalar r2_min = `mnr2'
   return local r2_worst "`lnr2'"
   return scalar r2_max = `mxr2'
   return local r2_best "`lxr2'"
   return scalar r2c_mean = 1 - `swcl'/`stot'
   return scalar r2c_min = `mnr2cl'
   return local r2c_worst "`lnr2cl'"
   return scalar r2c_max = `mxr2cl'
   return local r2c_best "`lxr2cl'"
   return scalar r2a_mean = 1 - `sarl'/`stot'
   return scalar r2a_min = `mnr2a'
   return local r2a_worst "`lnr2a'"
   return scalar r2a_max = `mxr2a'
   return local r2a_best "`lxr2a'"
   return matrix cri all
   restore
   end

