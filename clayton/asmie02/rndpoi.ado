*!version 1.0.0   1993 Joseph Hilbe, Walter Linde-Zwirble       (sg44: STB-28)
* Poisson distribution random number generator
* Example: rndpoi 1000 4  [set obs 1000;  4 is the mean]

program define rndpoi
   version 3.1
   cap drop xp
   qui  {
      local cases `1'
      set obs `cases'
      mac shift
      local xm `1'
      mac shift
      tempvar em t ds sum1 ran1    
      local g = exp(-`xm')
      gen `em'= -1
      gen `t' = 1.0
      gen `ran1' = uniform()
      gen `ds' = 1
      egen `sum1' = sum(`ds')
      noi di in gr "( Generating " _c
      while  `sum1' > 0 {
          replace `em' = `em'+ 1 if (`ds'==1)
          replace `t' = `t' * `ran1' if (`ds'==1)
          replace `ds'=0 if (`g' > `t')
          replace `ran1' = uniform()
          drop `sum1'
          egen `sum1' = sum(`ds')
          noi di in gr "." _c
      }
      noi di in gr " )"
      gen xp = int(`em'+0.5)
      noi di in bl "Variable " in ye "xp " in bl "created." 
   }
end

