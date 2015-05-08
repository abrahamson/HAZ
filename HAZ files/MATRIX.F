c matrix.f

      subroutine ztilde (zrow,m,n,z)
c     This subroutine converts a matrix from partitioned format
c     to tilde format.

      integer zrow
      real*8 z(zrow,1)

      do 100 i=1,m/2
        do 90 j=1,n
          z(i,j+n) = -z(i+m/2,j)
          z(i+m/2,j+n) = z(i,j)
  90    continue
 100  continue
      n = 2*n
      return
      end

c-------------------------------------------------------------
      subroutine mult (a,arow,m1,n1,b,brow,m2,n2,c,crow)
c     This subroutine multiplies two matrices.
c       C = A B
c     where A is m1 x n1  and  B is m2 x n2

      INTEGER arow, brow, crow
      double precision a(arow,1), b(brow,1), c(crow,1)

      if (n1 .ne. m2) stop 99 
      do 100 i=1,m1
        do 90 j=1,n2
          c(i,j) = 0.0
          do 80 ii=1,n1
            c(i,j) = c(i,j) + a(i,ii)*b(ii,j)
  80      continue
  90    continue
  100 continue
      return
      end

c-------------------------------------------------------------
      SUBROUTINE simul( N, A, X, EPS, INDIC, NRC , DETER )
C
C     WHEN INDIC IS NEGATIVE, SIMUL COMPUTES THE INVERSE OF THE N BY
C     N MATRIX A IN PLACE.  WHEN INDIC IS ZERO, SIMUL COMPUTES THE
C     N SOLUTIONS X(1)...X(N) CORRESPONDING TO THE SET OF LINEAR
C     EQUATIONS WITH AUGMENTED MATRIX OF COEFFICIENTS IN THE N BY
C     N+1 ARRAY A AND IN ADDITION COMPUTES THE INVERSE OF THE
C     COEFFICIENT MATRIX IN PLACE AS ABOVE.  IF INDIC IS POSITIVE,
C     THE SET OF LINEAR EQUATIONS IS SOLVED BUT THE INVERSE IS NOT
C     COMPUTED IN PLACE.  THE GAUSS-JORDAN COMPLETE ELIMINATION METHOD
C     IS EMPLOYED WITH THE MAXIMUM PIVOT STRATEGY.  ROW AND COLUMN
C     SUBSCRIPTS OF SUCCESSIVE PIVOT ELEMENTS ARE SAVED IN ORDER IN
C     THE IROW AND JCOL ARRAYS RESPECTIVELY.  K IS THE PIVOT COUNTER,
C     PIVOT THE ALGEBRAIC VALUE OF THE PIVOT ELEMENT, MAX
C     THE NUMBER OF COLUMNS IN A AND DETER THE DETERMINANT OF THE
C     COEFFICIENTS MATRIX.  THE SOLUTIONS ARE COMPUTED IN THE (N+1)TH
C     COLUMN OF A AND THEN UNSCRAMBLED AND PUT IN PROPER ORDER IN
C     X(1)...X(N) USING THE PIVOT SUBSCRIPT INFORMATION AVAILABLE
C     IN THE IROW AND JCOL ARRAYS.  THE SIGN OF THE DETERMINANT IS
C     ADJUSTED, IF NECESSARY, BY DETERMINING IF AN ODD OR EVEN NUMBER
C     OF PAIRWISE INTERCHANGES IS REQUIRED TO PUT THE ELEMENTS OF THE 
C     JORD ARRAY IN ASCENDING SEQUENCE WHERE JORD(IROW(I)) = JCOL(I).
C     IF THE INVERSE IS REQUIRED, IT IS UNSCRAMBLED IN PLACE USING
C     Y(1)...Y(N) AS TEMPORARY STORAGE.  THE VALUE OF THE DETERMINANT
C     IS RETURNED AS THE VALUE OF THE FUNCTION.  SHOULD THE POTENTIAL
C     PIVOT OF LARGEST MAGNITUDE BE SMALLER IN MAGNITUDE THAN EPS,
C     THE MATRIX IS CONSIDERED TO BE SINGULAR AND A TRUE ZERO IS
C     RETURNED AS THE VALUE OF THE FUNCTION.
C
      IMPLICIT REAL*8(A-H, O-Z)
      real*8 A(nrc,nrc), X(N), Y(200), eps, deter
      integer IROW(200), JCOL(200), JORD(200)
C
      MAX = N
      IF ( INDIC.GE.0 )  MAX = N + 1
C
C     ..... IS N LARGER THAN 78 .....
      IF ( N.LE.78 )  GO TO 5
      WRITE (6,200)
      DETER = 0.
      RETURN
C
C     ..... BEGIN ELIMINATION PROCEDURE .....
 5    DETER = 1.
      DO 18 K = 1, N
      KM1 = K - 1
C     ..... SEARCH FOR THE PIVOT ELEMENT .....
      PIVOT = 0.
      DO 11 I = 1, N
      DO 11 J = 1, N
C     ..... SCAN IROW AND JCOL ARRAYS FOR INVALID PIVOT SUBSCRIPTS ....
      IF ( K.EQ.1 ) GO TO 9
      DO 8 ISCAN = 1, KM1
      DO 8 JSCAN = 1, KM1
      IF ( I.EQ.IROW(ISCAN) ) GO TO 11
      IF ( J.EQ.JCOL(JSCAN) ) GO TO 11
 8    CONTINUE
 9    CONTINUE
      PIVOT = A(I,J)
      IROW(K) = I
      JCOL(K) = J
 11   CONTINUE
C
C     ..... INSURE THAT SELECTED PIVOT IS LARGER THAN EPS .....
      IF ( DABS(PIVOT).GT.EPS ) GO TO 13
      DETER = 0.
      RETURN
C
C     ..... UPDATE THE DETERMINANT VALUE .....
 13   IROWK = IROW(K)
      JCOLK = JCOL(K)
      DETER = DETER*PIVOT
C
C     ..... NORMALIZE PIVOT ROW ELEMENTS .....
      DO 14 J = 1, MAX
 14   A(IROWK,J) = A(IROWK,J)/PIVOT
C
C     ..... CARRY OUT ELIMINATION AND  DEVELOP INVERSE .....
      A(IROWK,JCOLK) = 1./PIVOT
      DO 18 I = 1, N
      AIJCK = A(I,JCOLK)
      IF ( I.EQ.IROWK ) GO TO 18
      A(I,JCOLK) = - AIJCK/PIVOT
      DO 17 J = 1, MAX
 17   IF ( J.NE.JCOLK ) A(I,J) = A(I,J) - AIJCK*A(IROWK,J)
 18   CONTINUE
C
C     ..... ORDER SOLUTION VALUES (IF ANY) AND CREATE JORD ARRAY .....
      DO 20 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
      JORD(IROWI) = JCOLI
 20   IF ( INDIC.GE.0 ) X(JCOLI) = A(IROWI,MAX)
C
C     ..... ADJUST SIGN OF DETERMINANT .....
      INTCH = 0
      NM1 = N - 1
      DO 22 I = 1, NM1
      IP1 = I + 1
      DO 22 J = IP1,N
      IF ( JORD(J).GE.JORD(I) ) GO TO 22
      JTEMP = JORD(J)
      JORD(J) = JORD(I)
      JORD(I) = JTEMP
      INTCH = INTCH + 1
 22   CONTINUE
      IF( INTCH/2*2.NE.INTCH ) DETER = - DETER
C
C     ..... IF INDIC IS POSITIVE RETURN WITH RESULTS .....
      IF ( INDIC.LE.0 ) GO TO 26
c     DETER = DETER
      RETURN
C
C     ..... IF INDIC IS NEGATIVE OR ZERO, UNSCRAMBLE THE INVERSE
C           FIRST BY ROWS .....
 26   DO 28 J = 1, N
      DO 27 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
 27   Y(JCOLI) = A(IROWI,J)
      DO 28 I = 1, N
 28   A(I,J) = Y(I)
C     ..... THEN BY COLUMNS .....
      DO 30 I = 1, N
      DO 29 J = 1, N
      IROWJ = IROW(J)
      JCOLJ = JCOL(J)
 29   Y(IROWJ) = A(I,JCOLJ)
      DO 30 J = 1, N
 30   A(I,J) = Y(J)
C
C     ..... RETURN FOR INDIC NEGATIVE OR ZERO .....
c     DETER = DETER
      RETURN
C
C     ..... FORMAT FOR OUTPUT STATEMENT .....
 200  FORMAT( 10H0N TOO BIG )
C
      END
	  
C-------------------------------------------------------------------
      subroutine SetHypocenter (iH, nH, h, probH )

c     Use a uniform distribution for the location of the hypocenter
      h = (iH-0.5)/nH
      probH = 1./nH

      return
      end
