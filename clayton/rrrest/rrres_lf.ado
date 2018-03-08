*! Name    : rrres_lf.ado (method lf evaluator for rrrest)
*! Date    : 11 Oct 2002
*! Version : 0.1
*! Author  : Chris Wallace
*! Email   : chris.wallace@lshtm.ac.uk

	/* method lf evaluator for rrrest */
program define rrres_lf
version 6
	set more off
	args lnf e1 e2 eta
	tempvar lnfj theta theta1 s delta u1 u2
	qui {
		gen double `u1' = 1/(1+exp(-`e1'))
		gen double `u2' = 1/(1+exp(-`e2'))
    
		/* generate lnfj - double */
      gen double `theta' = exp(`eta')
		gen double `theta1' = `theta'-1
		gen double `s' = sqrt( (1+`theta1'*(`u1'+`u2'))^2 - /*
*/                            4*`u1'*`u2'*`theta1'*`theta' )
		gen double `delta' = ( 1+`theta1'*(`u1'+`u2') - `s' )/( 2*`theta1' )
		gen double `lnfj' = `delta' if $ML_y1==1 & $ML_y2==1
		replace `lnfj' = `u1'-`delta' if $ML_y1==1 & $ML_y2==0
		replace `lnfj' = `u2'-`delta' if $ML_y1==0 & $ML_y2==1
		replace `lnfj' = 1-`u1'-`u2'+`delta' if $ML_y1==0 & $ML_y2==0
		replace `lnf' = ln(`lnfj')
		}
end
