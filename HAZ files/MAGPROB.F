c     magprob.f                             

      subroutine magProb ( mag, maxMag, minMag, magStep, beta, iFlt,            
     1                     pmag, nParam, nWidth,                                
     2                     magRecur, mpdf_param, ExpMeanMo, CharMeanMo )
      include 'pfrisk.h'            
      real mag, maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), beta(MAX_FLT,MAXPARAM,MAX_WIDTH),           
     1     minMag(MAX_FLT), magStep(MAX_FLT), magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH)                  
      real pmag(MAXPARAM,MAX_WIDTH), mU, mL, magU, magL, mag1, mag2         
      real maxMagWA, mMag, mChar                
      real zmL, zmU, zmagU, zmagL, d    
      real*8 pmL, pmU, pmagL, pmagU       
      real meanMag, sigma                                   
      real mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,5)
      integer iFlt, nParam(MAX_FLT,MAX_WIDTH), nWidth(MAX_FLT)
      real ExpMeanMo(MAXPARAM,MAX_WIDTH), CharMeanMo(MAXPARAM,MAX_WIDTH)                        
      real magcut, betaAC, magnorm, sumnorm, mUnorm, mLnorm
      real aknorm, pmagnorm(MAXPARAM,MAX_WIDTH) 
      
c     Set magnitude range               
      mU = mag + magStep(iFlt)/2.       
      mL = mag - magStep(iFlt)/2.   

c     Compute the mag probability for each parameter variation                  
      do iWidth=1,nWidth(iFlt)          
        do iParam=1,nParam(iFlt,iWidth)        
          beta1 = beta(iFlt,iParam,iWidth) 
          if ( magRecur(iFlt,iParam,iWidth) .ne. 10 ) then 
            magU = maxMag(iFlt,iParam,iWidth)                
          else
            magU = 8.5              
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
              if (beta(iFlt,iParam,iWidth) .ne. 0.0) then
              ak = 1. / (1. - exp(-beta(iFlt,iParam,iWidth)                            
     1           *(magU - magL)) )      
              else                               
                write (*,'( 2x,''Error - b-value = 0.0'')')                                   
                stop 99                               
              endif                              

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE EXPONENTIAL MODEL.
              pmag(iParam,iWidth) = ak*(exp(-beta1*(mL - magL))                                
     1                   - exp(-beta1*(mU - magL)) )                            

c           MAXIMUM MAGNITUDE MODEL       
            elseif ( magRecur(iFlt,iParam,iWidth) .eq. 3. ) then                         
cnjg              meanMag = magU - mpdf_param(iFlt,iParam,iwidth,3)                          
              meanMag = magU - mpdf_param(iFlt,iParam,iwidth,1)                          
              sigma = mpdf_param(iFlt,iParam,iwidth,2)    
c           Check for Delta Function, i.e., Sigma=0
              if (sigma .eq. 0.0) then
                 if (meanMag .ge. mL .and. meanMag .lt. mU) then
c            Set probability of Magnitude equal to 1.0 for this 
c            magnitude bin.
                    pmag(iParam,iWidth) = 1.0               
                 else
                    pmag(iParam,iWidth) = 0.0               
                 endif
              else
                 zmL = (mL-meanMag)/sigma    
                 zmU = (mU-meanMag)/sigma    
                 zmagU=(magU-meanMag)/sigma  
                 zmagL=(magL-meanMag)/sigma  
                 call NDTR(zmL,Real(pmL),d)        
                 call NDTR(zmU,Real(pmU),d)        
                 call NDTR(zmagL,Real(pmagL),d)           
                 call NDTR(zmagU,Real(pmagU),d)        

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE MAXIMUM MAGNITUDE MODEL.
                 pmag(iParam,iWidth) = (pmL - pmU)/(pmagL-pmagU)                        
              endif            
 
c           CHARACTERISTIC MODEL (Youngs and Coppersmith)                         
            elseif (magRecur(iFlt,iParam,iWidth) .eq. 0. ) then                         
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
                 call NDTR(zmL,pmL,d)        
                 call NDTR(zmU,pmU,d)        
                 call NDTR(zmagL,pmagL,d)           
                 call NDTR(zmagU,pmagU,d)           

c             COMPUTE PROBABILITY OF MAGNITUDE BETWEEN mL AND mU                
c             FOR THE MAXIMUM MAGNITUDE MODEL.

                 if (mag .ge. meanMag-mpdf_param(iFlt,iParam,iwidth,4)*sigma 
     1         .and. mag .le. meanMag+mpdf_param(iFlt,iParam,iwidth,4)*Sigma) then
c                    pChar = (pmL - pmU)/(pmagL-pmagU)                        
c   change to truncate only on th high end
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

c                 write (*,'(f8.3,2e12.4)') magnorm, pmagnorm(iParam,iWidth), sumnorm

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
              Btail = mpdf_param(iFlt,iParam,iwidth,3)
              SigM = mpdf_param(iFlt,iParam,iwidth,1)
              Fract_Exp = mpdf_param(iFlt,iParam,iwidth,4)
              mChar = maxMag(iFlt,iParam,iWidth)
              b_value = beta(iFlt,iParam,iWidth)/alog(10.0)

              call Calc_WA_Pmag1 ( mChar, sigM, b_value, bTail, Fract_exp, MaxMagWA, 
     1                  minMag(iFlt), ML, MU, pMag1)
              pMag(iParam,iWidth) = pMag1
          endif                                
         endif     

        enddo                                  

      enddo   

      return                            
      end                               
   
C------WAACY Subroutine for Magnitude Prob Values --------------

      subroutine Calc_WA_Pmag1 ( mChar, sigM, b_value, bTail, F, Mmax,
     1       Mmin, ML, MU, WA_Pmag)

      implicit none
      real mChar, sigM, b_value, btail, Fract_exp, Mmax, F
      real M1, M2
      real Mmin, stepM, ML, MU, WA_Pmag
      real pdf1(10000), PHI_15, PHI1, x, step1, mag1
      real beta, betaTail, c1, c2, c3, d1, d2, alpha, t1
      real Mo_bar_exp1, Mo_bar_exp2, Mo_bar_char
      integer nStep, i
      real temp1, temp2, x1, x2, x3, x4, x5

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

c     Compute the pdf with a small step size
      step1 = 0.01
      nStep = int ( (MU - ML ) / step1 )

c     Compute the probability of mag bewteen ML and MU
      WA_Pmag = 0.
      do i=1,nStep
        mag1 = ML + (i-0.5)*step1       
        if ( mag1 .lt. Mmin) then
          pdf1(i) = 0.
        elseif ( mag1 .lt. M1 ) then
          pdf1(i) = alpha * beta*exp(-beta*(mag1-Mmin)) / ( 1. - exp(-beta*(M1-Mmin)) )
        elseif ( mag1 .lt. M2 ) then
          x1 = (1.-alpha)
          x2 = c2 / (1. + c1*c2*c3)
          x3 = 1./  (t1*sigM)
          x4 = exp( -((mag1-mChar)**2) / (2*sigM**2) )
          pdf1(i) = x1 * x2 * x3 * x4

        elseif ( mag1 .lt. Mmax ) then
          pdf1(i) = (1.-alpha) * c1*c2 / ( (1.+c1*c2*c3) ) * exp( -betaTail*(mag1-M2) ) 
        else
          pdf1(i) = 0.
        endif
        WA_Pmag = WA_Pmag + pdf1(i)*step1   
      enddo

      return
      end

C------WAACY Subroutine for Magnitude Prob Values Used for Rates --------------

      subroutine Calc_WA_Pmag2 ( mChar, sigM, b_value, bTail, F, Mmax,
     1       Mmin, WA_Pmag, stepM, nMag)

      implicit none
      real mChar, sigM, b_value, btail, Fract_exp, Mmax, F
      real M1, M2
      real Mmin, stepM, ML, MU, WA_Pmag(1)
      real pdf1(10000), PHI_15, PHI1, x, step1, mag1
      real beta, betaTail, c1, c2, c3, d1, d2, alpha, t1
      real Mo_bar_exp1, Mo_bar_exp2, Mo_bar_char
      integer nStep, i, nMag, iBin, iMag
      real temp1, temp2
      betaTail = bTail * alog(10.)
      beta = b_value * alog(10.)

c     Calculate WAACY terms
      M1 = mChar - 0.25
      M2 = mChar + 1.5*sigM
      t1 = sqrt(2*3.1415926)
      c1 = (1./ (t1*sigM)) * exp( -((1.5*sigM)**2) / (2*sigM**2) )
      x = -0.25 / sigM
c      call NDTR(x,PHI1) 
      phi1 = 0.159       
      PHI_15 = 0.933
      c2 = 1./ ( PHI1-(1.-PHI_15) )
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
      d2 = (1./ ( 1. + c1*c2*c3)) * Mo_bar_Char + c1*c2*c3*Mo_bar_exp2 / ( 1. + c1*c2*c3 ) 
      alpha = d2*F / ( d1*(1.-F)+d2*F )
c      write (*,'( 2x,f10.4, 2x,''Mmax'')') Mmax

c     Compute the pdf with a small step size
      step1 = 0.01
      do i=1,nMag
        mag1 = Mmin + (i-1)*step1  
        if ( mag1 .lt. Mmin) then
          pdf1(i) = 0.
        elseif ( mag1 .lt. M1 ) then
          pdf1(i) = alpha * beta*exp(-beta*(mag1-Mmin)) / ( 1. - exp(-beta*(M1-Mmin)) )
c          write (*,'( ''bin 1'', 5e12.4)') pdf1(i), alpha,  beta*exp(-beta*(mag1-Mmin))
c     1           , ( 1. - exp(-beta*(M1-Mmin)) )


        elseif ( mag1 .lt. M2 ) then
          pdf1(i) = (1-alpha) * c2 / ( (1.+c1*c2*c3) * t1*sigM ) * 
     1             exp( -((mag1-mChar)**2) / (2*sigM**2) )
c          pdf1(i) =  c2 / ( (1.+c1*c2*c3) * t1*sigM ) * 
c     1             exp( -((mag1-mChar)**2) / (2*sigM**2) )

        elseif ( mag1 .lt. Mmax ) then
          pdf1(i) = (1.-alpha) * c1*c2 / (1.+c1*c2*c3) * exp( -betaTail*(mag1-M2) ) 
c          pdf1(i) = c1*c2 / (1.+c1*c2*c3) * exp( -betaTail*(mag1-M2) ) 
          
        else
          pdf1(i) = 0.
        endif

      enddo

c     Compute the probability of mag bewteen ML and MU
      do i=1,nMag
        WA_Pmag(i) = pdf1(i)*step1   
      enddo
        
      return
      end                             
