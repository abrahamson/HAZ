c     magprob.f                             

      subroutine magProb ( mag, maxMag, minMag, magStep, beta, iFlt,            
     1                     pmag, nParam, nWidth,                                
     2                     magRecur, mpdf_param, ExpMeanMo, CharMeanMo )
      include 'pfrisk.h'            
      real mag, maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), beta(MAX_FLT,MAXPARAM,MAX_WIDTH),           
     1     minMag(MAX_FLT), magStep(MAX_FLT), magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH)                  
      real pmag(MAXPARAM,MAX_WIDTH), mU, mL, magU, magL, mag1, mag2                     
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
          magU = maxMag(iFlt,iParam,iWidth)                

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
                 call NDTR3(zmL,pmL)        
                 call NDTR3(zmU,pmU)        
                 call NDTR3(zmagL,pmagL)           
                 call NDTR3(zmagU,pmagU)       

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
       
       
c               write (*,'(i4,2x,f8.3,3e10.3)') nsummag, mag, pmag(iParam,iWidth), sumnorm, pmagnorm(iParam,iWidth)*sumnorm
c               write (*,*)
c               pause 
              
          endif                                
         endif     

        enddo                                  

      enddo   

      return                            
      end                               
