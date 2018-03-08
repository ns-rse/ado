This directory contains David Clayton's Stata software. The software and 
associated files are grouped into a number of "packages". These are designed 
to be downloaded and installed from within Stata, using the "net" commands.

Packages currently available are:

genassoc	Programs for analysis of genetic association studies
gamenu		Rudimentary graphical user interface to the genassoc package 
exassoc		Exercises on genetic association studies
htSNP		Programs to aid choice of "haplotype tagging" SNPs
htSNP2		Mark2 version of htSNP package
ibdreg		Linkage analyis with covariates by IBD regression (regression
		of IBD state upon covariates)
permutation	Stata programs for Monte Carlo permutation tests
dropout_death	Stata code and for dropout & death in longitudinal studies
asmie02		Materials for course on "Advanced Statistical Models in 
		Epidemiology", Erasmus Summer Programme, Aug 26-30, 2002
rrrest		Chris Wallace's program for analysis of recurrence relative
		risks in family data.

Most commands require Stata 8.0 or later

For those unfamiliar with the Stata "net" command, the following sequence
of commands, typed into Stata itself, will download and install any package
(xxx represents the package name, e.g. genassoc):

net from http://www-gene.cimr.cam.ac.uk/clayton/software/stata
net install xxx
net get xxx

The second command (net install) installs the package, including any online 
help files. The third command (net get) dowloads any additional material 
(documentation files, example datasets etc) into the current directory.


David Clayton

Diabetes and Inflammation Laboratory		Tel: 44 (0)1223 762669 
Cambridge Institute for Medical Research	Fax: 44 (0)1223 762102
Wellcome Trust/MRC Building			david.clayton@cimr.cam.ac.uk
Addenbrooke's Hospital, Cambridge, CB2 2XY	www-gene.cimr.cam.ac.uk/clayton
