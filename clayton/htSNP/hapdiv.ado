*! version 1.2 , DGC October 2001

program define hapdiv, rclass
   version 7.0
   preserve
   syntax [varlist] [fw aw iw pw/], HTsnps(varlist)
   tempvar n sum wt use 
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
   gen byte `use'=1
   markout `use' `varlist' `wt'
   quietly keep if `use'
   foreach locus of local varlist {
     quietly count if `use' & `locus'!=1 & `locus'!=2
     if r(N) > 0 {
       di in red "`locus' is not coded 1 or 2"
       exit 1
     }
   }
   sort `htsnps'
   quietly {
     summ `wt' , meanonly
     local N = r(sum)
     by `htsnps': gen `n' = sum(`wt')
     by `htsnps': replace `n'=. if _n!=_N
     count if `n'!=.
     local nhap = r(N)
     gen `sum' = sum(`n'^2) if `n'!=.
     local sn2 = `sum'[_N]
     local chance = 1 - `sn2'/(`N'^2)
     if "`weight'"=="" | "`weight'"=="fweight" {
       local chance = `chance'*`N'/(`N'-1)
     }
     drop `sum'
   }
   di
   di "htSNPs: `htsnps'"
   di "Number of htSNP haplotypes: `nhap'" 
   di
   di "--------------------------------------------------------------------"
   di _skip(26) "Mean diversity" _skip(7) "Percent"  _skip(8) "Chance"
   di %16s = "Locus" %13s = "Total"  %13s = "Residual" _continue
   di %13s = "captured" %13s = "corrected"
   di "--------------------------------------------------------------------"
   quietly {
     tokenize `varlist'
     local sres = 0
     local stot = 0
     local pmin = 100
     local rmax = 0
     while "`1'" != "" {
       summ `wt'  if `1'==2, meanonly
       local N2 = r(sum)
       tempvar div sum
       by `htsnps': gen `div' = sum(`wt'*(`1'==2))
       replace `div' = 2*`div'*(`n' - `div') if `n'!=.
       gen `sum' = sum(`div') if `n'!=.
       local resid = `sum'[_N]
       local total = 2*`N2'*(`N'-`N2')
       local pcap = (`total' - `resid')/`total'
       local pcor = (`pcap' - `chance')/(1-`chance')
       local rm = `resid'/(`sn2')
       if `pcap'<`pmin' {
         local pmin = `pcap'
         local smin "`1'"
       }
       if `rm'>`rmax' {
         local rmax = `rm'
         local smax "`1'"
       }
       noisily {
         di %16s = "`1'" %13.3f = `total'/(`N'*`N') _continue
         di %13.3f = `resid'/`sn2' %13.3f = 100*`pcap' %13.3f = 100*`pcor'
       }
       local sres = `sres' + `resid'
       local stot = `stot' + `total'
       drop `div' `sum'
       mac shift
     }
   }
   local pcap = (`stot'-`sres')/`stot'
   local pcor = (`pcap' - `chance')/(1-`chance')
   local rm = `sres'/(`sn2')
   di "--------------------------------------------------------------------"
   di %16s = "Total" %13.3f = `stot'/(`N'*`N') _continue
   di %13.3f = `rm' %13.3f = 100*`pcap' %13.3f = 100*`pcor'
   di "--------------------------------------------------------------------"
   di
   return scalar p_tot = 100*`pcap'
   return scalar p_min = 100*`pmin'
   return scalar p_chance = 100*`chance'
   return scalar rmd_tot = `rm'
   return scalar rmd_max = `rmax'
   return local p_worst "`smin'"
   return local rmd_worst "`smax'"
   end

