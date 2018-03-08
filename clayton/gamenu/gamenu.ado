*! Version 1.0, DGC Sept 2001
program def gamenu
  version 7.0
  * Initialization 
  global  GA_pedvars ""
  * Popout menu structure 
  window menu clear
  * window menu append popout "sysmenu" "GenAssoc"
  * Redefine standard menu and add GenAssoc
  window menu popout "MyTop"
  window menu append popout "MyTop" "sFile"
  window menu append popout "MyTop" "Edit"
  window menu append popout "MyTop" "Prefs"
  window menu append popout "MyTop" "Window"
  window menu append popout "MyTop" "sHelp"
  window menu append popout "MyTop" "GenAssoc"
  window menu append popout "GenAssoc" "Data management"
  window menu append string "Data management" "Read Stata dataset" "XEQ use"
  window menu append string "Data management" "Save dataset" "XEQ save"
  window menu append string "Data management" "Read spreadsheet" "_run_ginsheet"
  window menu append string "Data management" "Define Pedigree Variables" /*
    */ "_def_pedvars"
  window menu append string "Data management" "Rename loci" "_run_grename"
  window menu append string "Data management" "Reshape data file" "_run_gresh"
  window menu append popout "GenAssoc" "Recode"
  window menu append string "Recode" "Group rare alleles" "_run_grprare"
  window menu append string "Recode" "Create genotype variable" "_run_gtype"
  window menu append string "Recode" "Create haplotype variables" "_run_htype"
  window menu append popout "GenAssoc" "Tabulate" 
  window menu append string "Tabulate" "Allele frequencies" "_run_gtab"
  window menu append string "Tabulate" "Case-parent triads" "_run_trios"
  window menu append popout "GenAssoc" "TDT etc"
  window menu append string "TDT etc" "TDT" "_run_tdt"
  window menu append string "TDT etc" "Genotype RR" "_run_gtrr"
  window menu append string "TDT etc" "Pseudo-CC" "_run_pscc"
  window menu append string "TDT etc" "Parental origin" "_run_origin"
  window menu append popout  "GenAssoc" "Regression" 
  window menu append string "Regression" "Fit" "_run_regress"
  window menu append string "Regression" "Test (drop)" "_run_rtest"
  window menu append string "Regression" "Test (add)" /*
         */ "_gamenu_error Not implemented yet"
  * Start menu system 
  * window menu set "sysmenu"
  window menu set "MyTop"
  end
