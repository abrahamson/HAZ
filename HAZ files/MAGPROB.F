                           
      subroutine magProb ( sourceType, mag, maxMag, minMag, magStep, beta, 
     1                     iFlt, iFltWidth, pmag, nParam, nWidth, magRecur, 
     2                     mpdf_param, ExpMeanMo, CharMeanMo, Rup1_flag )

      implicit none
      include 'pfrisk.h'  
          
      integer rup1_flag, iFlt, nParam(MAX_FLT,MAX_WIDTH), nWidth(MAX_FLT),
     1        iWidth, iParam, nsummag, imagsum, sourceType, iFltWidth
      real mag, maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), beta(MAX_FLT,MAXPARAM,MAX_WIDTH),           
     1     minMag(MAX_FLT), magStep(MAX_FLT), magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH),
     2     pmag(MAXPARAM,MAX_WIDTH), mU, mL, magU, magL, mag1, mag2,           
     3     maxMagWA, mChar, zmL, zmU, zmagU, zmagL, meanMag, sigma, 
     4     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6), magcut, betaAC 
      real ExpMeanMo(MAXPARAM,MAX_WIDTH), CharMeanMo(MAXPARAM,MAX_WIDTH),
     1     magnorm, sumnorm, mUnorm, mLnorm, aknorm, pmagnorm(MAXPARAM,MAX_WIDTH),
     2     beta1, ak, deltaM1, deltaM2, p1, p2, pChar, d, pExp, c1, Btail,
     3     SigM, Fract_Exp, b_value, pMag1
      real*8 pmL, pmU, pmagL, pmagU                                         

c     set pmag for SourceType 7
      if (sourceType .eq. 7) then
        do iParam=1,nParam(iFlt,iFltWidth)
          pmag(iParam,iFltWidth) = 1.0
        enddo
        goto 100

c     all other SourceTypes       
      else
    
      rup1_flag = 0   

c     Compute the mag probability for each parameter variation                  
      do iWidth=1,nWidth(iFlt)          
        do iParam=1,nParam(iFlt,iWidth)
c         Set magnitude range                     
          mL = mag - magStep(iFlt)/2. 
          mU = mag + magStep(iFlt)/2.    
          beta1 = beta(iFlt,iParam,iWidth) 
          if ( magRecur(iFlt,iParam,iWidth) .ne. 10 ) then 
            magU = maxMag(iFlt,iParam,iWidth)                
          else
            if ( mpdf_param(iFlt,iParam,iwidth,5) .eq. 0. ) then
              magU =  mpdf_param(iFlt,iParam,iwidth,2)
            else
              magU = maxMag(iFlt,iParam,iWidth) + mpdf_param(iFlt,iParam,iwidth,5)
            endif
          endif
          
c         Lower magnitude equal to the minimum magnitude for each fault.		  
          magL = minMag(iFlt)        
          
c         CHECK IF MAGNITUDE STEP EXCEEDS THE MAX MAG                                      
          if ( mL .gt. magU ) then                              
            pmag(iParam,iWidth) = 0.           
          else                                 
      
c           WHICH MAG RECURRENCE RELATION?
c           TRUNCATED EXPONENTIAL MODEL   
            if ( magRecur(iFlt,iParam,iWidth) .eq. 1. ) then 
c             Reset mU to magU if necessary (when last mag is not a full step size)
              if (mU .gt. magU) then
                mU = magU
              endif                     
              if (beta(iFlt,iParam,iWidth) .ne. 0.0) then
                ak = 1. / (1. - exp(-beta(iFlt,iParam,iWidth)                            
     1          *(magU - magL)) )      
              else                               
                write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                stop 99                               
              endif                              

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE EXPONENTIAL MODEL.
              pmag(iParam,iWidth) = ak*(exp(-beta1*(mL - magL))                                
     1                   - exp(-beta1*(mU - magL)) )  
c             set flag if this is part of the large eqk                          
              if ( mag .gt. MagU-mpdf_param(iFlt,iParam,iwidth,3)) then
                Rup1_flag = 1
              else
                Rup1_flag = 0
              endif

c           MAXIMUM MAGNITUDE MODEL       
            elseif ( magRecur(iFlt,iParam,iWidth) .eq. 3. ) then 
              meanMag = magU - mpdf_param(iFlt,iParam,iwidth,1)                          
              sigma = mpdf_param(iFlt,iParam,iwidth,2)    
c           Check for Delta Function, i.e., mLSigma=0
              if (sigma .eq. 0.0) then
                 if (meanMag .ge. mL .and. meanMag .lt. mU) then
c            Set probability of Magnitude equal to 1.0 for this 
c            magnitude bin.
                    pmag(iParam,iWidth) = 1.0               
                 else
                    pmag(iParam,iWidth) = 0.0               
                 endif
              else
c               Reset mU to magU if necessary (when last mag is not a full step size)
                if (mU .gt. magU) then
                  mU = magU
                endif                                                  
                zmL = (mL-meanMag)/sigma    
                zmU = (mU-meanMag)/sigma    
                zmagU=(magU-meanMag)/sigma  
                zmagL=(magL-meanMag)/sigma 
                if ( zmU .lt. -4. .or. zmu .gt. 4. ) then
                  pmag(iParam,iWidth) = 0.
                else
                  call NDTR3(zmL,pmL)        
                  call NDTR3(zmU,pmU)        
                  call NDTR3(zmagL,pmagL)           
                  call NDTR3(zmagU,pmagU)        
                  pmag(iParam,iWidth) = (pmL - pmU)/(pmagL-pmagU)
                endif                        
              endif      

c             All ruptures are large mag ruptures for syn rup
              Rup1_flag = 1              
      
c           CHARACTERISTIC MODEL (Youngs and Coppersmith)                         
            elseif (magRecur(iFlt,iParam,iWidth) .eq. 0. ) then 
c             Reset mU to magU if necessary (when last mag is not a full step size)
              if (mU .gt. magU) then
                mU = magU
              endif                          
              deltaM1 = mpdf_param(iFlt,iParam,iwidth,1)                                 
              deltaM2 = mpdf_param(iFlt,iParam,iwidth,2)
              mag1 = magU - deltaM1 - deltaM2    
              mag2 = magU - deltaM1              
c             Mag range (mL to mU) is completely in the exponenital 
c             part of the pdf
              if ( mU .lt. mag2 ) then          
                pmag(iParam,iWidth) =            
     1             exp(-beta1*(mL - minMag(iFlt)))                           
     1            - exp(-beta1 *(mU - minMag(iFlt)))                        
              elseif (mL .gt. mag2 ) then                               
c               Mag range (mL to mU) is completely in the char mag 
c               part of pdf
                if ( mU .lt. magU ) then
                  pmag(iParam,iWidth) = beta1*(exp(-beta1                                  
     1               *(mag1-minMag(iFlt)))) * magStep(iFlt)   
                else
                  pmag(iParam,iWidth) = beta1*(exp(-beta1                                  
     1               *(mag1-minMag(iFlt)))) * (magU-mL)   
                endif
              else
c               Mag range (mL to mU) stradles the exponential/char 
c               part of the pdf
                p1 = exp(-beta1*(mL - minMag(iFlt)))                           
     1                - exp(-beta1 *(mag2 - minMag(iFlt)))                        
                p2 = beta1*(exp(-beta1                                  
     1               *(mag1-minMag(iFlt)))) * (mU-mag2)


c               COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c               FOR THE CHARACTERISTIC MAGNITUDE MODEL.
                pmag(iParam,iWidth) = p1 + p2
              endif
                        
c             NORMALIZE                   
              ak = (1-exp(-beta1*(magU-deltaM1-minMag(iFlt))))                           
     1         + beta1*exp(-beta1*(mag1-minMag(iFlt))) * deltaM1                
              pmag(iParam,iWidth) = pmag(iParam,iWidth) / ak   

c             Set large rup flag for syn rupture
              if ( mag .gt. mag2) then
                Rup1_flag = 1
              else
                Rup1_flag = 0
              endif
                          
       
c         WG99 model
          elseif ( magRecur(iFlt,iParam,iWidth) .eq. 4. ) then         

c           *** char eqk part - normal distribution
              meanMag = magU - mpdf_param(iFlt,iParam,iwidth,1)                          
              sigma = mpdf_param(iFlt,iParam,iwidth,2)    

c           Check for Delta Function, i.e., Sigma=0
              if (sigma .eq. 0.0) then
                 if (meanMag .ge. mL .and. meanMag .lt. mU) then
c            Set probability of Magnitude equal to 1.0 for this 
c            magnitude bin.
                    pChar = 1.0               
                 else
                    pChar = 0.0               
                 endif
              else
                 zmL = (mL-meanMag)/sigma    
                 zmU = (mU-meanMag)/sigma    
                 zmagU=(magU-meanMag)/sigma  
                 zmagL=(magL-meanMag)/sigma  
                 call NDTR(zmL,reaL(pmL),d)        
                 call NDTR(zmU,real(pmU),d)        
                 call NDTR(zmagL,real(pmagL),d)           
                 call NDTR(zmagU,real(pmagU),d)           

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE MAXIMUM MAGNITUDE MODEL.

                 if (mag .ge. meanMag-mpdf_param(iFlt,iParam,iwidth,4)*sigma 
     1         .and. mag .le. meanMag+mpdf_param(iFlt,iParam,iwidth,4)*Sigma) then                      
c   change to truncate only on the high end
                    pChar = (pmL - pmU)/(1.-pmagU)                        
                 else
                    pChar = 0.0
                 endif

              endif

c             *** Exp part
C Reset Maxmim magnitude for Exp part equal to mean mag minus 2.0 sigma
              magU = meanmag - mpdf_param(iFlt,iParam,iwidth,4)*sigma

              if (beta(iFlt,iParam,iWidth) .ne. 0.0) then
               ak = 1. / (1. - exp(-beta(iFlt,iParam,iWidth)                            
     1           *(magU - magL)) )      
              else                               
                write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                stop 99                               
              endif                              

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE EXPONENTIAL MODEL.
              if (mag .le. magU) then
                 pExp = ak*(exp(-beta1*(mL - magL))                                
     1                   - exp(-beta1*(mU - magL)) )                            
              else
                 pExp = 0.0
              endif

c             Combine the char and exp parts using the fraction of rate 
c             that goes into each model
              if (ExpMeanMo(iParam,iWidth)+CharMeanMo(iParam,iWidth) .ne. 0.0) then

                 c1 = mpdf_param(iFlt,iParam,iWidth,5)*charMeanMo(iParam,iWidth) /
     1                  ( (1.0-mpdf_param(iFlt,iParam,iWidth,5))*ExpMeanMo(iParam,iWidth) + 
     2                         mpdf_param(iFlt,iParam,iWidth,5)*CharMeanMo(iParam,iWidth) )

                 pmag(iParam,iWidth) = c1*pExp + (1.0-c1)*pChar

               else
                 pmag(iParam,iWidth) = 0.0             
               endif

c           BC Hydro Alternative Characteristic Model   
            elseif ( magRecur(iFlt,iParam,iWidth) .eq. 6. ) then                      
              betaAC = mpdf_param(iFlt,iParam,iWidth,2)*alog(10.0)
              magcut = magU - mpdf_param(iFlt,iParam,iWidth,1)
              if (mag .le. magcut) then 
                 if (beta(iFlt,iParam,iWidth) .ne. 0.0) then
                    ak = 1. / (1. - exp(-beta(iFlt,iParam,iWidth)                            
     1                   *(magU - magL)) )      
                 else                               
                   write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                   stop 99                               
                 endif                               
c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE EXPONENTIAL MODEL.
                 pmag(iParam,iWidth) = ak*(exp(-beta1*(mL - magL))                                
     1                   - exp(-beta1*(mU - magL)) )                            
              else
                 if (betaAC .ne. 0.0) then
                    ak = 1. / (1. - exp(-betaAC*(magU - magcut) ) )  
                 else                               
                   write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                   stop 99                               
                 endif                               
c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE SECOND EXPONENTIAL MODEL.
                 pmag(iParam,iWidth) = ak*(exp(-betaAC*(mL - magcut))                                
     1                   - exp(-betaAC*(mU - magcut)) )                            
              endif

C             Compute total pmag over all magnitudes to normalize model by such that Sum(mpdf)=1
              sumnorm = 0.0
              nsummag = (magU - magL)/magstep(iFlt) + 1
              do imagsum=1,nsummag

                 magNorm = real(magL + magstep(iFlt)/2. + (imagsum-1)*magstep(iFlt))
                 mUnorm = magnorm + magstep(iFlt)/2.
                 mLnorm = magnorm - magstep(iFlt)/2.
                 
                 betaAC = mpdf_param(iFlt,iParam,iWidth,2)*alog(10.0)
                 magcut = magU - mpdf_param(iFlt,iParam,iWidth,1)
                 if (magNorm .le. magcut) then 
                    if (beta(iFlt,iParam,iWidth) .ne. 0.0) then
                       aknorm = 1. / (1. - exp(-beta(iFlt,iParam,iWidth)                            
     1                      *(magU - magL)) )      
                    else                               
                      write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                      stop 99                               
                    endif                               
c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE EXPONENTIAL MODEL.
                    pmagnorm(iParam,iWidth) = aknorm*(exp(-beta1*(mLnorm - magL))                                
     1                      - exp(-beta1*(mUnorm - magL)) )                            
                 else
                    if (betaAC .ne. 0.0) then
                       aknorm = 1. / (1. - exp(-betaAC*(magU - magcut) ) )  
                    else                               
                      write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                      stop 99                               
                    endif                               
c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE SECOND EXPONENTIAL MODEL.
                    pmagnorm(iParam,iWidth) = aknorm*(exp(-betaAC*(mLnorm - magcut))                                
     1                      - exp(-betaAC*(mUnorm - magcut)) )                            
                 endif

                 sumnorm = sumnorm + pmagnorm(iParam,iWidth)

              enddo
              
C             Now normalize BC Hydro Model by total 
              pmag(iParam,iWidth) = pmag(iParam,iWidth)/sumnorm

c           WAACY Model
            elseif ( magRecur(iFlt,iParam,iWidth) .eq. 10. ) then  
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
              if ( MaxMagWA .lt. mChar ) then
                write (*,'( 2x,''bad parameter for WAACY model, maxmag is too small'')')
                stop 99
              endif  

c             Brute force calc of Pmag.  Fix this later to be faster.
              call Calc_WA_Pmag1 ( mChar, sigM, b_value, bTail, Fract_exp, MaxMagWA, 
     1                  minMag(iFlt), ML, MU, pMag1)
              pMag(iParam,iWidth) = pMag1

            endif                                
          endif

        enddo                               
      enddo  
      
      endif 

 100  continue

      return                            
      end                               
   
C------WAACY Subroutine for Magnitude Prob Values --------------

      subroutine Calc_WA_Pmag1 ( mChar, sigM, b_value, bTail, F, Mmax,
     1       Mmin, ML, MU, WA_Pmag)

      implicit none
      real mChar, sigM, b_value, btail, Mmax, F
      real M1, M2
      real Mmin, ML, MU, WA_Pmag
      real pdf1, PHI_15, PHI1, x, step1, mag1
      real beta, betaTail, c1, c2, c3, d1, d2, alpha, t1
      real Mo_bar_exp1, Mo_bar_exp2, Mo_bar_char
      integer nStep, i
      real temp1, temp2, oneMinusAlpha, x2, x3, x4

      betaTail = bTail * alog(10.)
      beta = b_value * alog(10.)

c     Calculate WAACY terms
      M1 = mChar - 0.25
      M2 = mChar + 1.5*sigM
      t1 = sqrt(2*3.1415926)
      c1 = (1./ (t1*sigM)) * exp( -((1.5*sigM)**2) / (2*sigM**2) )
      x = -0.25 / sigM
      phi1 = 0.159
c      call NDTR(x,PHI1)        
      PHI_15 = 0.933
      c2 = 1./ ( 1.-PHI1-(1.-PHI_15) )
      c3 = (1. - exp(-betaTail*(Mmax-M2)) ) / betaTail

      temp1 =  beta * 10.**(16.05) * ( exp((-beta+3.45)*M1)-1. ) 
      temp2 = ( 1.-exp(-beta*M1))*(-beta+3.45)
      Mo_bar_exp1 = temp1/temp2

      Mo_bar_exp2 = betaTail * 10**16.05 * exp( betaTail*M2 ) *
     1              ( exp(( -betaTail+3.45)*Mmax) - exp(( -betaTail+3.45)*M2) ) /
     2            ( ( 1. - exp(-betaTail*(Mmax-M2)) ) * ( -betaTail + 3.45 ) )
      Mo_bar_Char = 10.**(1.5*MChar+16.05) * (2.63*(sigM-0.2) + 1.19) 

      d1 = Mo_bar_exp1 / ( exp( -beta*Mmin) - exp(-beta*M1) )
      d2 = (1./ ( 1. + c1*c2*c3)) * Mo_bar_Char + c1*c2*c3*Mo_bar_exp2 / ( 1. + c1*c2*c3 ) 
      alpha = d2*F / ( d1*(1.-F)+d2*F )
      oneMinusAlpha = d1*(1.-F) / ( d1*(1.-F)+d2*F )

c     Compute the pdf with a small step size
      step1 = 0.001
      nStep = int ( (MU - ML ) / step1 )

c     Compute the probability of mag bewteen ML and MU
      WA_Pmag = 0.
      do i=1,nStep
        mag1 = ML + (i-0.5)*step1       
        if ( mag1 .lt. Mmin) then
          pdf1 = 0.
        elseif ( mag1 .lt. M1 ) then
          pdf1 = alpha * beta*exp(-beta*(mag1-Mmin)) / ( 1. - exp(-beta*(M1-Mmin)) )
        elseif ( mag1 .lt. M2 ) then
          x2 = c2 / (1. + c1*c2*c3)
          x3 = 1./  (t1*sigM)
          x4 = exp( -((mag1-mChar)**2) / (2*sigM**2) )
          pdf1 = oneMinusAlpha * x2 * x3 * x4

        elseif ( mag1 .lt. Mmax ) then
          pdf1 = oneMinusAlpha * c1*c2 / ( (1.+c1*c2*c3) ) * exp( -betaTail*(mag1-M2) ) 
        else
          pdf1 = 0.
        endif
        WA_Pmag = WA_Pmag + pdf1*step1   
      enddo

      return
      end

C------WAACY Subroutine for Magnitude Prob Values Used for Rates --------------

      subroutine Calc_WA_Pmag2 ( mChar, sigM, b_value, bTail, F, Mmax,
     1       Mmin, WA_Pmag, stepM, nMag)

      implicit none
      real mChar, sigM, b_value, btail, Mmax, F
      real M1, M2, d
      real Mmin, stepM, WA_Pmag(1)
      real pdf1(10000), PHI_15, PHI1, x, mag1
      real*8  beta, betaTail, c1, c2, c3, t1
      real*8 d1, d2, alpha, a1, a2, alpha1, sum1, sum2
      real*8  Mo_bar_exp1, Mo_bar_exp2, Mo_bar_char
      integer i, nMag
      real*8  temp1, temp2
      betaTail = bTail * alog(10.)
      beta = b_value * alog(10.)

c     Calculate WAACY terms
      M1 = mChar - 0.25
      M2 = mChar + 1.5*sigM
      t1 = sqrt(2*3.1415926)
      c1 = (1./ (t1*sigM)) * exp( -((1.5*sigM)**2) / (2*sigM**2) )
      x = -0.25 / sigM
      call NDTR(x,PHI1,d) 
      phi1 = 1-phi1
      PHI_15 = 0.933
      c2 = 1./ ( PHI_15-phi1) 
      c3 = (1. - exp(-betaTail*(Mmax-M2)) ) / betaTail

c      Mo_bar_exp1 = beta * 10.**(16.05) * ( exp((-beta+3.45)*M1)-1. ) /
c     1              ( 1.-exp(-beta*M1)) * (-beta+3.45)

      temp1 =  beta * 10.**(16.05) * ( exp((-beta+3.45)*M1)-1. ) 
      temp2 = ( 1.-exp(-beta*M1))*(-beta+3.45)
      Mo_bar_exp1 = temp1/temp2

      Mo_bar_exp2 = betaTail * 10**(16.05) * exp( betaTail*M2 ) *
     1              ( exp(( -betaTail+3.45)*Mmax) - exp(( -betaTail+3.45)*M2) ) /
     2            ( ( 1. - exp(-betaTail*(Mmax-M2)) ) * ( -betaTail + 3.45 ) )
      Mo_bar_Char = 10.**(1.5*MChar+16.05) * (2.63*(sigM-0.2) + 1.19) 

      d1 = Mo_bar_exp1 / ( exp( -beta*Mmin) - exp(-beta*M1) )
      temp1 = c1*c2*c3
      d2 = (1./ ( 1. + c1*c2*c3)) * Mo_bar_Char + c1*c2*c3*Mo_bar_exp2 / ( 1. + c1*c2*c3 ) 
      a1 = d2*F
      a2 = d1*(1.-F)+d2*F
      alpha = d2*F / ( d1*(1.-F)+d2*F )
      alpha1 = d1*(1-F) / (d1*(1.-F)+d2*F)

c     Compute the pdf with a small step size
      do i=1,nMag
        mag1 = Mmin + (i-1)*stepM 
        if ( mag1 .lt. Mmin) then
          pdf1(i) = 0.
        elseif ( mag1 .lt. M1 ) then
          pdf1(i) = alpha * beta*exp(-beta*(mag1-Mmin)) / ( 1. - exp(-beta*(M1-Mmin)) )

        elseif ( mag1 .lt. M2 ) then
          pdf1(i) = alpha1 * c2 / ( (1.+c1*c2*c3) * t1*sigM ) * 
     1             exp( -((mag1-mChar)**2) / (2*sigM**2) )

        elseif ( mag1 .lt. Mmax ) then
          pdf1(i) = alpha1 * c1*c2 / (1.+c1*c2*c3) * exp( -betaTail*(mag1-M2) ) 
          
        else
          pdf1(i) = 0.
        endif

      enddo

c     Compute the probability of mag 
c     Also compute the moment less than M1 (sum1) 
c          and momment greater than M1 (sum2)
       sum1 = 0.
       sum2 = 0.
      do i=1,nMag
        mag1 = Mmin + (i-1)*stepM
        WA_Pmag(i) = pdf1(i)*stepM  
        if (mag1 .le. M1 ) then
          sum1 = sum1 + 10.**(1.5*mag1+16.05)*WA_Pmag(i)
        else
          sum2 = sum2 + 10.**(1.5*mag1+16.05)*WA_Pmag(i)
        endif
      enddo
        
      return
      end                             

C------WAACY Subroutine for to correct for rupture outside of modelled flt -----

      subroutine Calc_F2_Waacy ( mChar, sigM, b_value, bTail, F, Mmax,
     1       Mmin, F2 )

      implicit none     
      
      real mChar, sigM, b_value, btail, Mmax, F, M1, M2, d, Mmin,
     1     F2, PHI_15, PHI1, x
      real*8  beta, betaTail, c1, c2, c3, t1, d1, d2, alpha, a1, a2, 
     1        alpha1, Mo_bar_exp1, Mo_bar_exp2, Mo_bar_char, temp1, 
     2        temp2
     
      betaTail = bTail * alog(10.)
      beta = b_value * alog(10.)

c     Calculate WAACY terms
      M1 = mChar - 0.25
      M2 = mChar + 1.5*sigM
      t1 = sqrt(2*3.1415926)
      c1 = (1./ (t1*sigM)) * exp( -((1.5*sigM)**2) / (2*sigM**2) )
      x = -0.25 / sigM
      call NDTR(x,PHI1,d) 
      phi1 = 1-phi1
      PHI_15 = 0.933
      c2 = 1./ ( PHI_15-phi1) 
      c3 = (1. - exp(-betaTail*(Mmax-M2)) ) / betaTail

      temp1 =  beta * 10.**(16.05) * ( exp((-beta+3.45)*M1)-1. ) 
      temp2 = ( 1.-exp(-beta*M1))*(-beta+3.45)
      Mo_bar_exp1 = temp1/temp2

      Mo_bar_exp2 = betaTail * 10**(16.05) * exp( betaTail*M2 ) *
     1              ( exp(( -betaTail+3.45)*Mmax) - exp(( -betaTail+3.45)*M2) ) /
     2            ( ( 1. - exp(-betaTail*(Mmax-M2)) ) * ( -betaTail + 3.45 ) )
      Mo_bar_Char = 10.**(1.5*MChar+16.05) * (2.63*(sigM-0.2) + 1.19) 

      d1 = Mo_bar_exp1 / ( exp( -beta*Mmin) - exp(-beta*M1) )
      temp1 = c1*c2*c3
      d2 = (1./ ( 1. + c1*c2*c3)) * Mo_bar_Char + c1*c2*c3*Mo_bar_exp2 / ( 1. + c1*c2*c3 ) 
      a1 = d2*F
      a2 = d1*(1.-F)+d2*F
      alpha = d2*F / ( d1*(1.-F)+d2*F )
      alpha1 = d1*(1-F) / (d1*(1.-F)+d2*F)
        
      return
      end                             
