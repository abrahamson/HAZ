c ntdr.f                                                                        
                                                                                
      SUBROUTINE NDTR( X, P, D)                                                 
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                                
C          D IS DENSITY VALUE (OUTPUT).
      IF (X) 1,2,2                                                              
    1 AX = -X                                                                   
      GOTO 3                                                                    
    2 AX = X                                                                    
    3 IF ( AX-6.0 ) 5,4,4                                                       
    4 P = 1.                                                                    
      D = 0.                                                                    
      GOTO 6                                                                    
    5 T = 1. / (1.0 + 0.2316419 * AX)                                           
      D = 0.3989423 * EXP(-X*X / 2.0)                                           
      P = 1.0 - D*T*( (((1.330274*T - 1.821256)*T + 1.781478) * T-              
     1    0.3565638) * T + 0.3193815)                                           
    6 IF (X) 8,7,7                                                              
    7 P = 1.0 - P                                                               
    8 continue                                                                  
      RETURN                                                                    
      END                                                                       

c ------------------------------------------------------------------
      SUBROUTINE NDTR2( X, P)                                                 
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                                
C          D IS DENSITY VALUE (OUTPUT).
      real*8 p, x1
	  real*8 d1, d2, d3, d4, d5, d6
	  real x
	  data d1, d2, d3, d4, d5, d6 / 0.0498673470, 0.0211410061,
     1     0.0032776263, 0.0000380036, 0.0000488906, 0.0000053830 /
      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
	  endif
	  
	  p = 1. - 0.5*(  (((((d6*x1+d5)*x1+d4)*x1+d3)*x1+d2)*x1+d1)*
     1         x1+1.)**(-16.)
      if ( x .gt. 0. ) then
	    p = 1. - p
	  endif

      RETURN                                                                    
      END                                                                       

c ------------------------------------------------------------------
      SUBROUTINE NDTR3( X, P)    
C          Reference: Abramowitz and Stegan equation 7.1.26                                             
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                               
      real*8 p, x1, x2
	  real*8 p1, a1, a2, a3, a4, a5, t
	  real x
	  data p1, a1, a2, a3, a4, a5 / 0.3275911, 0.254829592, 
     1     -0.284496736, 1.421413741, -1.453152027, 1.061405429 /
      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
      endif
	  
        x2 = x1/(sqrt(2.))
        t = 1/(1+(p1*x2))
        p = 1-0.5*((a1*t)+(a2*(t**2))+(a3*(t**3))+(a4*(t**4))+(a5*(t**
     1      5)))*(exp(-(x2**2)))
     
      if ( x .gt. 0. ) then
	    p = 1. - p
      endif

      RETURN                                                                    
      END       
c ------------------------------------------------------------------

      real function pxceed ( eti, ti, siga, jj, i, sigTrunc )                   
      include 'pfrisk.h'                                                        
      real eti, ti(MAX_PROB,1), siga, sigTrunc 
                                                                                      
      W = (ti(JJ,I)-ETI)/SIGA                                                   
      CALL NDTR(W,G,D)                                                          
                                                                                
c     TRUNCATE DISTRIBUTION AT SigTrunc                                         
      if (w .gt. sigTrunc) then                                                 
        pxceed = 0.0                                                                   
      else                                                                      
        CALL NDTR(sigTrunc,g1,D)                                                       
        pxceed = (g-g1)*(1.+g1)                                                        
      endif 
      return                                                                    
      end                                                                       

c ------------------------------------------------------------------

      real*8 function pxceed2 ( eti, ti, siga, jj, i, sigTrunc )                   
      include 'pfrisk.h'                                                        
      real eti, ti(MAX_PROB,1), siga, sigTrunc, w
	  real*8 g, g1
                                                                                
      W = (ti(JJ,I)-ETI)/SIGA                                                   
      CALL NDTR2(W,G)     
                                                                                
c     TRUNCATE DISTRIBUTION AT SigTrunc                                         
      if (w .gt. sigTrunc) then                                                 
        pxceed2 = 0.0                                                                   
      else                                                                      
        CALL NDTR2(sigTrunc,g1)                                                       
        pxceed2 = (g-g1)*(1.+g1)                                                        
      endif                                                                     
      return                                                                    
      end  

c ------------------------------------------------------------------

      real*8 function pxceed3 ( eti, ti, siga, jj, i, sigTrunc )                   
      include 'pfrisk.h'                                                        
      real eti, ti(MAX_PROB,1), siga, sigTrunc, w
	  real*8 g, g1
                                                                                
      W = (ti(JJ,I)-ETI)/SIGA                                                   
      CALL NDTR3(W,G)     
                                                                                
c     TRUNCATE DISTRIBUTION AT SigTrunc                                         
      if (w .gt. sigTrunc) then                                                 
        pxceed3 = 0.0                                                                   
      else                                                                      
        CALL NDTR3(sigTrunc,g1)                                                       
        pxceed3 = (g-g1)*(1.+g1)                                                             
      endif                                                                     
      return                                                                    
      end                                                                       
c ------------------------------------------------------------------

      real*8 function probOccur ( eps1, eps2, sTrunc )                   
      include 'pfrisk.h'                                                        
      real eps1, eps2, sTrunc
      real*8 g1, g2, gt
                                                                                
      CALL NDTR2(-eps1,g1)                                                       
      CALL NDTR2(-eps2,g2)                                                       
      CALL NDTR2(-sTrunc,gt)                                                       

c     TRUNCATE DISTRIBUTION AT SigTrunc                                         
      if (eps1 .ge. sTrunc) then                                                 
        probOccur = 0.0                                                                   
      else                                                                      
        probOccur = (g2-g1)/(gt)                                                        
      endif                                                                     

      return                                                                    
      end                                                                       


c ------------------------------------------------------------------
                                                                                
      subroutine rupDimProb ( mag, coef, sigma, step, sigmaMax,                 
     1           rupDim, prob, iFlt, idim )                                     
                                                                                
c      include 'pfrisk.h'                                                       
                                                                                
      real mag, coef(2,1), sigma(1), rupDim, prob                                  
      real nSigma, nSigma_plus, nSigma_minus                                    
                                                                                
      nSigma = -sigmaMax + (idim-0.5)*step                                      
      nSigma_plus = (nSigma + step/2.)                                          
      nSigma_minus = (nSigma - step/2.)                                         
      rupDim = 10.**(coef(1,iflt)+coef(2,iflt)*mag+nSigma*sigma(iflt))              
                                                                                
c     Compute probability that (log) rupture dimension is between               
c     dim_log_minus and dim_log_plus                                            
      CALL NDTR ( sigmaMax, F0, D )                                             
      CALL NDTR ( nSigma_minus, F1, D )                                         
      CALL NDTR ( nSigma_plus, F2, D )                                          
      prob = (F1-F2)/(1-2*f0)                                                   
                                                                                
      return                                                                    
      end                                                                       
