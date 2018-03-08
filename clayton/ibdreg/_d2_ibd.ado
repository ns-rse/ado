 /* Evaluator */

program def  _d2_ibd
  version 7.0
  args todo b lnf g negH g1
  quietly {
    tempvar r2 eta pi qi p02 p12 p22 f d1 d2
    mleval `eta' = `b'
    gen `pi' = 1/(1 + exp(-`eta'))
    gen `qi' = 1 -`pi'
    gen `p02' = `qi'^2
    gen `p12' = 2*`pi'*`qi'
    gen `p22' = `pi'^2
    /* likelihood */
    gen `f' = cond($rat2!=., $rat0*`p02' + $rat1*`p12' + $rat2*`p22', /*
          */  cond($rat1!=., $rat0*`qi' + $rat1*`pi', 0))
    /* Derivatives of likelihood wrt pi */
    gen `d1'= cond($rat2!=., 2*(`pi'*($rat2+$rat0-2*$rat1)+$rat1-$rat0),/*
          */  cond($rat1!=., $rat1-$rat0, 0))
    gen `d2'= cond($rat2!=., 2*($rat2 - $rat1*2 + $rat0), 0)

    /* log likelihood and derivatives */

    replace `d1' = `pi'*`qi'*`d1'/`f'
    replace `g1' = `d1' if $ML_samp
    replace `d2' = (`pi'*`qi')^2*`d2'/`f' - `d1'^2
    replace `f' = log(`f')
    /* Add over subjects */
    tempname d1s d2s
    mlsum `lnf' = `f' if $ML_samp
    if `todo'==0 | `lnf'==. {exit}
    mlvecsum `lnf' `d1s' = `d1'  if $ML_samp, eq(1)
    matrix `g' = (`d1s')
    mlmatsum `lnf' `d2s' = -(`d2') if $ML_samp, eq(1)
    matrix `negH' = (`d2s')
  }
  end
