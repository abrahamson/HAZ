c ------------------------------------------------------------------  

      subroutine S31_Set_Rates ( sourceType, nParamVar, magRecur, rate, beta, minMag,           
     1           maxMag, iFlt, iWidth, faultArea, 
     1           RateParam, mpdf_param, magStep, RateType, charMeanMo, expMeanMo )                           

      implicit none
      include 'pfrisk.h'                       

      real magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH),                                
     1     beta(MAX_FLT,MAXPARAM,MAX_WIDTH),                                    
     1     minMag(MAX_FLT), maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH),                 
     1     rate(MAXPARAM,MAX_WIDTH),charMeanMo(MAXPARAM,MAX_WIDTH),    
     1     expMeanMo(MAXPARAM,MAX_WIDTH),                             
     2     RateParam(MAX_FLT,MAXPARAM,MAX_WIDTH),                                
     3     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6)                                     
      real magStep(MAX_FLT)             
      real mU, mL, rigidity, momentRate1, momentRate2                 
      real c, beta1, t1, t2, t3, c1, deltaM1, deltaM2, deltaM3             
      real mean, sigma, zmagL, zmagU, pmagL, pmagU, dd, mag                     
      real mU1, mL1, mch, mchr, mexp, aexp, achr, m1
      real*8 sum
      real term1, term2, term3, scale, pL1, pU1, contrmaxmag
      integer iParam, iFlt, i, nParamVar(MAX_FLT,MAX_WIDTH)                  
      integer nmstep, iWidth, i1, sourceType                    
      real RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), meanMoment2
      real maxmag1, magstep1, texp, contrexp, rate_M_gt_0
      real bAC, gamma, bGR, fGR, fAC, gAC, meanMoChar, pRatio
      real c2, c3, c4, bM2, bM1, scale1, faultArea, x, deltamac
      real meanMoRelease1, meanMoRelease2, meanMoRelease3, fix_step

c     No rate calculation needed for SourceType 7
      if (sourceType .eq. 7) then
        goto 100
      endif
      
      rigidity = 3.0e11       

      i = iWidth          
        do iParam=1,nParamVar(iFlt,i)          
     
          beta1 = beta(iFlt,iParam,i) 
          mU = maxMag(iFlt,iParam,i)           
          mL = minMag(iFlt)          

          deltaM1 = mpdf_param(iFlt,iParam,i,1)     
          deltaM2 = mpdf_param(iFlt,iParam,i,2) 
          deltaM3 = mpdf_param(iFlt,iParam,i,3) 
          
c         Is this a slip-rate or moment rate?
          if (RateType(iFlt,iParam,i) .eq. 1 .or. RateType(iFlt,iParam,i) .eq. 4) then

            if ( RateType(iFlt,iParam,i) .eq. 1 ) then
              momentRate2 = rateParam(iFlt,iParam,i)*faultArea * rigidity * 1.0e9
c              write (80,'( 3i5,2e12.4)') iflt, i, iParam,faultArea, rateParam(iFlt,iParam,i)
            else
              momentRate2 = rateParam(iFlt,iParam,i)
            endif          

            if ( magRecur(iFlt,iParam,i) .eq. 1 ) then                                  

c            EXPONENTIAL MODEL           
c            See Equation XX of Appendix #.
             c1 = 1.5*alog(10.)                 
             t1 = beta1 * exp(beta1*mL+c1*10.7) /                                       
     1           ( 1. - exp(-beta1*(mu-mL)))                                    
             t2 = -beta1+c1                     
             t3 = ( exp(t2*mU) - exp(t2*mL) )   
             momentRate1 = t1 / t2 * t3         

C            Compute the mean moment between Mag=0 and Mag=ML
             mL1 = 0.0
             c1 = 1.5*alog(10.)
             t1 = beta1 * exp(beta1*mL1+c1*10.7) /                                       
     1           ( 1. - exp(-beta1*(mL-mL1)))                                    
             t2 = -beta1+c1                     
             t3 = ( exp(t2*mL) - exp(t2*mL1) )   
             meanMoment2 = t1 / t2 * t3         

c            Find the fraction of the momentRate that is from Mag=0 and Mag=ML
             t1 = exp(-beta1*mL)-exp(-beta1*mU)
             t2 = 1 - t1
             scale = t1*momentRate1 / (t1*momentRate1 + t2*meanMoment2) 
            
             if (momentRate2 .eq. 0.) then
                rate(iParam,i) = 0.0
             else
                rate(iParam,i) = momentRate2 * scale / momentRate1
                
             endif

           elseif (magRecur(iFlt,iParam,i) .eq. 0.) then                        

c            CHARACTERISTIC MODEL         

c            calculate mean moment between Mag=5 and MMax   
             mch = mU - 0.5*deltam1
             c1 = 1.5*alog(10.)                 
             t1 = exp(-beta1*(mU-deltaM1-mL))   
             t2 = exp(-beta1*(mU-deltaM1-deltaM2-mL))                                   
             t3 = exp(-beta1*(mU-deltaM1-mL))                                   
             c = beta1*t2*deltam1/(1.0-t3)
             term1 = 1.0/(1.0+c)         
             term2 = c1 - beta1
             term3 = exp(beta1*mL+16.05*alog(10.0))
             texp = (term1*beta1*exp(beta1*mL+16.05*alog(10.)))        
             texp = texp/(1.0-t3)*(term2)        
             texp = texp*(exp(term2*(mch-0.5*deltam1))-exp(term2*mL))			 
             if (mL .gt. mU-deltaM1 ) then            
                mexp = 0.             
             else            
                aexp = (term1*beta1*exp(beta1*mL+16.05*alog(10.)))/
     1              (1.0-t3)     
                mexp = (aexp/(-beta1+c1))*(exp((mch-0.5*deltam1)*
     1              (-beta1+c1)) - exp(mL*(-beta1+c1)))     
             endif
             if (mL .gt. mU ) then             
                mchr = 0.
             else             
                achr = (term1*beta1*t2*exp(16.05*alog(10.))) /
     1              (1.0-t3)     
                mchr = (achr/c1)*(exp(c1*(mch+0.5*deltam1)) -
     1              exp(c1*(mch-0.5*deltam1)))     
             endif

             meanMoRelease1 = mexp + mchr

c            calculate mean moment between Mag=0 and Mag=5
             ML1=0.0
             c2 = (0.5*beta1*exp(-beta1*(mch-ML1-1.25)))/(1-exp(-beta1
     1            *(mch-ML1-0.25)))
             c3 = beta1/((1+c2)*(1-exp(-beta1*(mch-ML1-0.25))))
             c4 = c3/(-beta1+(alog(10.)*1.5))
             bM2 = exp((-beta1*(ML-ML1))+(alog(10.)*((1.5*ML)+16.05)))
             bM1 = exp((-beta1*(ML1-ML1))+(alog(10.)*((1.5*ML1)+16.05)))
             meanMoRelease2 = c4*(bM2-bM1)

c            calculate mean moment between Mag=0 and MMax   
             mch = mU - 0.5*deltam1
             c1 = 1.5*alog(10.)                 
             t1 = exp(-beta1*(mU-deltaM1-mL1))   
             t2 = exp(-beta1*(mU-deltaM1-deltaM2-mL1))                                   
             t3 = exp(-beta1*(mU-deltaM1-mL1))                                   
             c = beta1*t2*deltam1/(1.0-t3)
             term1 = 1.0/(1.0+c)         
             term2 = c1 - beta1
             term3 = exp(beta1*mL1+16.05*alog(10.0))
             texp = (term1*beta1*exp(beta1*mL1+16.05*alog(10.)))        
             texp = texp/(1.0-t3)*(term2)        
             texp = texp*(exp(term2*(mch-0.5*deltam1))-exp(term2*mL1))			 
             if (mL .gt. mU-deltaM1 ) then            
                mexp = 0.             
             else            
                aexp = (term1*beta1*exp(beta1*mL1+16.05*alog(10.)))/
     1              (1.0-t3)     
                mexp = (aexp/(-beta1+c1))*(exp((mch-0.5*deltam1)*
     1              (-beta1+c1)) - exp(mL1*(-beta1+c1)))     
             endif
             if (mL1 .gt. mU ) then             
                mchr = 0.
             else             
                achr = (term1*beta1*t2*exp(16.05*alog(10.))) /
     1              (1.0-t3)     
                mchr = (achr/c1)*(exp(c1*(mch+0.5*deltam1)) -
     1              exp(c1*(mch-0.5*deltam1)))     
             endif

             meanMoRelease3 = mexp + mchr

c        calculate the scale factor
	   if (momentRate2 .eq. 0.) then
             scale1 = 0.0
           else
             scale1 = (momentRate2/meanMoRelease3)/ (momentRate2/
     1                (meanMoRelease3-meanMoRelease2)) 
           endif 

c        calculate the rate
  
             if ( meanMoRelease1 .eq. 0. ) then              
                rate(iParam,i) = 0.                  
             else            
                rate(iParam,i) = momentRate2*scale1/meanMoRelease1  
                               
             endif		 
           elseif (magRecur(iFlt,iParam,i) .eq. 3.) then                              

c            SINGLE MAXIMUM MAGNITUDE MODEL (magRecur = 3)                       
             mean =  mU - mpdf_param(iFlt,iParam,i,1) 
             sigma = mpdf_param(iFlt,iParam,i,2)       

             if (sigma .eq. 0.0) then 
                sum = 10.**(1.5*mean+16.05)
             else
c            Use a fixed mag step of 0.01 for getting the moment balance
             fix_step = 0.01
             nmstep = int((mU - mL)/fix_step)  
             sum = 0.                           
             mag = mL + fix_step/2.        
             zmagU=(mU-mean)/sigma       
             zmagL=(mL-mean)/sigma       

             call S27_NDTR(zmagL,pmagL,dd)          
             call S27_NDTR(zmagU,pmagU,dd)   
       
             do i1 = 1,nmstep                   
               mL1 = mag - fix_step/2.    
               mU1 = mag + fix_step/2.    
               call S27_NDTR((mL1-mean)/sigma,pL1,dd)                                      
               call S27_NDTR((mU1-mean)/sigma,pU1,dd)                                       
               sum = sum + (pL1-pU1)/(1.-pmagU)*(10.**(1.5*mag+16.05))    
               mag = mag + fix_step       
             enddo                       
             endif

             rate(iParam,i) = momentRate2/sum  

c     WAACY Model             
          elseif (magRecur(iFlt,iParam,i) .eq. 10 ) then  
             call S31_calc_sum_waacy ( sum, mpdf_param, maxMag, beta, minMag, iFlt, iParam, iwidth, faultArea, pRatio )  
             rate_M_gt_0 = momentRate2/sum
 
c            Set the rate for M> Mmin (scale the rate for M>0 by pRatio,
c            where pRatio is the ratio of eqk with M>Mmin / M>0
             rate(iParam,i) = rate_M_gt_0 * pRatio   

C..........Working Group Model...................
           elseif (magRecur(iFlt,iParam,i) .eq. 4.) then     
             write (*,'( 2x,''WG99 model with slip-rates not working yet'')')
             stop 99

C..........Bi-exponential Distribution...................
           elseif (magRecur(iFlt,iParam,i) .eq. 5.) then     
             write (*,'( 2x,''Bi-Exponential Distribution with slip-rates not working yet'')')
             stop 99

C..........BC Hydro - Alternative Characteristic Model...................
           elseif (magRecur(iFlt,iParam,i) .eq. 6.) then     
             write (*,'( 2x,''BC Hydro Alt Characteristic model with slip-rates not working yet'')')
             stop 99

           else
              write (*,'( 2x,''Bad Mag Recurrence Model'')')                              
              write (*,'( 2x,''Check fault input file  '')') 
              write (*,'( 2x, ''for fault number = '',i5)') iFlt
              stop 99                          
           endif                               
         endif

C ***** ACTIVITY RATE OPTION FOR SOURCES *****           
c        Is this an activity rate?
         if (RateType(iFlt,iParam,i) .eq. 2) then

c........Working Group Case with activity rates..........................
            if (magRecur(iFlt,iParam,i) .eq. 4.) then
            
c            Find mean moment per event for the char part  
             mean =  mU - mpdf_param(iFlt,iParam,i,1) 
             sigma = mpdf_param(iFlt,iParam,i,2)       
             maxMag1 = mean + mpdf_param(iFlt,iParam,i,4)*sigma
             contrexp = mpdf_param(iFlt,iParam,i,5)
             contrmaxmag = 1.0 - contrexp

             magStep1 = 0.01

             if (sigma .eq. 0.0) then 
                meanMoChar = 10.**(1.5*mean+16.05)
             else
                nmstep = int((maxmag1 - (mean-2.0*sigma))/magStep1)  
                sum = 0.                           
                mag = mean - 2.0*sigma + magStep1/2.        
                zmagU=(maxMag1-mean)/sigma       
                zmagL=((mean-2.0*sigma) - mean)/sigma       
                call S27_NDTR(zmagL,pmagL,dd)          
                call S27_NDTR(zmagU,pmagU,dd)   
       
                do i1 = 1,nmstep                   
                  mL1 = mag - magStep1/2.    
                  mU1 = mag + magStep1/2.    
                  call S27_NDTR((mL1-mean)/sigma,pL1,dd)                                      
                  call S27_NDTR((mU1-mean)/sigma,pU1,dd)                                      
                  sum = sum + (pL1-pU1)/(pmagL-pmagU)*                                    
     1            (10.**(1.5*mag+16.05))    
                  mag = mag + magStep1      
                enddo                       
             endif

             charMeanMo(iParam,i) = sum   

c            EXPONENTIAL PART           
c            See Equation XX of Appendix #.
c Set upper magnitude for exponential part to 2Sigma less than mean mag
             mu = mean - 2.0*sigma

             c1 = 1.5*alog(10.)                 
             t1 = beta1 * exp(beta1*mL+c1*10.7) /                                       
     1           ( 1. - exp(-beta1*(mu-mL)))                                    
             t2 = -beta1+c1                     
             t3 = ( exp(t2*mU) - exp(t2*mL) )   
             momentRate1 = t1 / t2 * t3         

C            Compute the mean moment between Mag=0 and Mag=ML
             mL1 = 0.0
             c1 = 1.5*alog(10.)
             t1 = beta1 * exp(beta1*mL1+c1*10.7) /                                       
     1           ( 1. - exp(-beta1*(mL-mL1)))                                    
             t2 = -beta1+c1                     
             t3 = ( exp(t2*mL) - exp(t2*mL1) )   
             meanMoment2 = t1 / t2 * t3         

c            Find the fraction of the momentRate that is from Mag=0 and Mag=ML
             t1 = exp(-beta1*mL)-exp(-beta1*mU)
             t2 = 1 - t1
             scale = t1*momentRate1 / (t1*momentRate1 + t2*meanMoment2) 
            
             if (momentRate2 .eq. 0.) then
                expMeanMo(iParam,i) = 0.0
             else
                expMeanMo(iParam,i) = momentRate2 * scale / momentRate1
             endif
              ExpMeanMo(iParam,i) = momentRate1

C...........Set rate as the combination of the exponential and maxmagnitude cases......
c             rate(iParam,i) = rateexp(iParami,i) + ratechar(iParam,i)

            endif

c........BC Hydro Alternative Characteristic Model........................
            if (magRecur(iFlt,iParam,i) .eq. 6.) then
               deltamac = mpdf_param(iFlt,iParam,i,1)
               bAC = mpdf_param(iFlt,iParam,i,2)
               gamma = mpdf_param(iFlt,iParam,i,3)
               bGR = beta1/alog(10.0)
               fGR = 10**(bGR*(mL-mU+deltamac))
               fAC = 10**(-bAC*deltamac)
               gAC = 10**(-1.5*deltamac)
   
               term1 = bAC*(1.5-bGR)*(1.0-fGR)
               term2 = bGR*(1.5-bAC)*gAC*fGR*(1.0-fAC)

               rate(iParam,i)=RateParam(iFlt,iParam,i)*(1.0+((1.0-gamma)/(gamma)) *
     1              term1/term2)
            endif

c           Simply set the activity rate for all other cases.
            rate(iParam,i)=RateParam(iFlt,iParam,i)   

         endif

c        Is this a recurrence INterval?                                     
         if (RateType(iFlt,iParam,i) .eq. 3) then
c           write (*,'( i5)') magRecur(iFlt,iParam,i)
           
c          max mag model
           if (magRecur(iFlt,iParam,i) .eq. 3) then                           
              rate(iParam,i) =1./RateParam(iFlt,iParam,i) 
              
c          Y&C model                                    
           elseif (magRecur(iFlt,iParam,i) .eq. 0 ) then                        

             if ( mL .ge. mU-deltaM1 ) then     
               rate(iParam,i) = 1./RateParam(iFlt,iParam,i)                             
             else                               
               x = exp(beta1*(mU-deltaM1-deltaM2-mL))/deltaM1
               rate(iParam,i) = 1./RateParam(iFlt,iParam,i) *                           
     1             (1 + x/beta1*(1-exp(-beta1*(mU-deltaM1-mL))))       
             endif   
             
c          truncated exp model
           elseif (magRecur(iFlt,iParam,i) .eq. 1 ) then
             m1 = mu - deltaM2                        
             t1 = -exp(-beta1*(mu-mL)) + exp(-beta1*(m1-mL))
             t2 = 1. - exp(-beta1*(mU-mL))
             rate(iParam,i) = (1./RateParam(iFlt,iParam,i) ) / (t1/t2)

c........does not work with Working Group Case
           elseif (magRecur(iFlt,iParam,i) .eq. 4.) then
               write (*,*) 'Recurence interval does not currently work'
               write (*,*) 'with working group model!!!!'
               stop 99
c........does not work with Bi-Exponential Model
           elseif (magRecur(iFlt,iParam,i) .eq. 5.) then
               write (*,*) 'Recurence interval does not currently work'
               write (*,*) 'with Bi-Exponential model!!!!'
               stop 99
c........does not work with BC Hydro Alternative Characteristic Model
           elseif (magRecur(iFlt,iParam,i) .eq. 6.) then
               write (*,*) 'Recurence interval does not currently work'
               write (*,*) 'with BC Hydro Alternative model!!!!'
               stop 99
           endif
         endif           
        enddo                                  

 100  continue

      return                            
      end                               

c -----------------------------------------------------------

      subroutine S31_calc_sum_waacy ( sum, mpdf_param, maxMag, beta, minmag,
     1          iFlt, iParam, iwidth, faultArea, pRatio)  

      implicit none
      include 'pfrisk.h'                

      real  mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6),                                     
     1       beta(MAX_FLT,MAXPARAM,MAX_WIDTH), 
     1       maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real minMag(MAX_FLT), Mmin
      real MaxMagWA, Btail, SigM, Fract_Exp, mChar, b_value, stepM
      real*8 sum, moment, sum1, sum2
      real WA_PMag(10000)
      real cumProb(10000)
      real mag, M1
      integer iMag, nMag, iFlt, iParam, iWidth, iMag1
      real pratio
      real faultArea, areaRatio, area_rup


c     Set WAACY model parameters (haz45 version)
      if ( mpdf_param(iFlt,iParam,iwidth,5) .eq. 0. ) then
        MaxMagWA =  mpdf_param(iFlt,iParam,iwidth,2)
      else
        MaxMagWA = maxMag(iFlt,iParam,iWidth) + mpdf_param(iFlt,iParam,iwidth,5)
      endif
      if (MaxMagWA .gt. mpdf_param(iFlt,iParam,iwidth,6) ) then
        MaxMagWA = mpdf_param(iFlt,iParam,iwidth,6)
      endif
      Btail = mpdf_param(iFlt,iParam,iwidth,3)
      SigM = mpdf_param(iFlt,iParam,iwidth,1)
      Fract_Exp = mpdf_param(iFlt,iParam,iwidth,4)
      mChar = maxMag(iFlt,iParam,iWidth)
      b_value = beta(iFlt,iParam,iWidth)/alog(10.0)

c     start integration at Mag=0 for momment balance
      Mmin = 0
      stepM =0.001
      nMag = inT ( (MaxMagWA - Mmin) / stepM )
c      write (*,'( i5)') nMag
     
c     Compute the probabilities of magnitudes using WAACY 
c     without correction for rupture past the modelled fault      
      call S25_Calc_WA_Pmag2 ( mChar, sigM, b_value, bTail, Fract_Exp, 
     1                     MaxMagWA, Mmin, WA_Pmag, stepM, nMag )

c     Adjust the moment release to account for the part of the rupture
c     that is not modelled.

c     First, Initialize moment sum
      sum1 = 0.
      sum2 = 0.
      M1 = mChar - 0.25

c     Compute the moment * mag pdf for the part of the rupture this is modelled.
c     this is separated into the moment in the exp part (sum1) and the moment
c     in the char part (sum2)
      do iMag=1,nMag
        mag = Mmin + (iMag-0.5)*stepM
        moment = 10.**(1.5*mag+16.05)

c       Scale the moment from the eqk for the part that is released on the modelled fault
c       Just use log(A)= M-4 for now
        area_rup = 10.**(mag -4)
        areaRatio = area_rup / faultArea
        if (areaRatio .gt. 1. ) then
          moment = moment / areaRatio
        endif
        if ( mag .LT. M1 ) then
          sum1 = sum1 + moment*WA_Pmag(iMag)
        else
          sum2 = sum2 + moment*WA_Pmag(iMag)
        endif
      enddo

c     Find the fraction of the moment in the exp part after the correction for 
c     moment released in non-modelled ruptures
      Fract_Exp = (Fract_Exp / (sum1 / (sum1 + sum2)) ) * Fract_Exp  

c     Reset the value in the array so this this corrected fraction is 
c     later used for the mag prob
      mpdf_param(iFlt,iParam,iwidth,4)= Fract_Exp 

c     Recompute the probabilities of magnitudes using WAACY 
c     with the correction for rupture past the modelled fault      
      call S25_Calc_WA_Pmag2 ( mChar, sigM, b_value, bTail, Fract_Exp, 
     1                     MaxMagWA, Mmin, WA_Pmag, stepM, nMag )

c     Compute the moment * mag pdf for the part of the rupture this is modelled.
      sum = 0.
      do iMag=1,nMag
        mag = Mmin + (iMag-0.5)*stepM
        moment = 10.**(1.5*mag+16.05)

c       Scale the moment from the eqk for the part that is released on the modelled fault
c       Just use log(A)= M-4 for now
        area_rup = 10.**(mag -4)
        areaRatio = area_rup / faultArea
        if (areaRatio .gt. 1. ) then
          moment = moment / areaRatio
        endif
        sum = sum + moment*WA_Pmag(iMag)
      enddo

c     Find the cumulative rate which will be used to scale the rate from M>0 to M>Mmin
      cumProb(nMag) = WA_Pmag(nMag) 
      do iMag=nMag-1,1,-1
        cumProb(iMag) = cumProb(iMag+1) + WA_Pmag(iMag)
      enddo     

c     Set the ratio of Prob for M>Mmin to M>0
c     Note: cumProb(1) should be unity.
      iMag1 = Int( minMag(iFlt) / stepM ) 
      pRatio = cumProb(iMag1) / cumProb(1)

      return
      end
              
