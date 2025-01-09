!source code dgemm obtained from https://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f_source.html

!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
! \par Purpose:
!  =============
!
! \verbatim
!
! DGEMM  performs one of the matrix-matrix operations
!
!    C := alpha*op( A )*op( B ) + beta*C,
!
! where  op( X ) is one of
!
!    op( X ) = X   or   op( X ) = X**T,
!
! alpha and beta are scalars, and A, B and C are matrices, with op( A )
! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
! \endverbatim
!
!  Arguments:
!  ==========
!
! \param[in] TRANSA
! \verbatim
!          TRANSA is CHARACTER*1
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A**T.
!
!              TRANSA = 'C' or 'c',  op( A ) = A**T.
! \endverbatim
!
! \param[in] TRANSB
! \verbatim
!          TRANSB is CHARACTER*1
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B**T.
!
!              TRANSB = 'C' or 'c',  op( B ) = B**T.
! \endverbatim
!
! \param[in] M
! \verbatim
!          M is INTEGER
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
! \endverbatim
!
! \param[in] N
! \verbatim
!          N is INTEGER
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
! \endverbatim
!
! \param[in] K
! \verbatim
!          K is INTEGER
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
! \endverbatim
!
! \param[in] ALPHA
! \verbatim
!          ALPHA is DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
! \endverbatim
!
! \param[in] A
! \verbatim
!          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
! \endverbatim
!
! \param[in] LDA
! \verbatim
!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
! \endverbatim
!
! \param[in] B
! \verbatim
!          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
! \endverbatim
!
! \param[in] LDB
! \verbatim
!          LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
! \endverbatim
!
! \param[in] BETA
! \verbatim
!          BETA is DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
! \endverbatim
!
! \param[in,out] C
! \verbatim
!          C is DOUBLE PRECISION array, dimension ( LDC, N )
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
! \endverbatim
!
! \param[in] LDC
! \verbatim
!          LDC is INTEGER
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
! \endverbatim
!
! Authors:
! ========
!
! \author Univ. of Tennessee
! \author Univ. of California Berkeley
! \author Univ. of Colorado Denver
! \author NAG Ltd.
!
! \ingroup gemm
!
! \par Further Details:
! =====================
!
! \verbatim
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
! \endverbatim
!
! =====================================================================
Module dgemm_module
    implicit none
    contains

    subroutine dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC,Nmult)
        !
        ! -- Reference BLAS level3 routine --
        ! -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        ! -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !    .. Scalar Arguments ..
        REAL(kind=8), intent(in) :: ALPHA,BETA
        INTEGER, intent(in) ::  K,LDA,LDB,LDC,M,N
        INTEGER, intent(out) :: Nmult
        CHARACTER,intent(in) :: TRANSA,TRANSB
!    ..
!    .. Array Arguments ..
        REAL(kind=8),intent(in) :: A(LDA,*), B(LDB,*)
        REAL(kind=8),intent(out) ::  C(LDC,*)

!    ..
!
! =====================================================================
!
!    .. External Functions ..
      !log10ICAL LSAME
      !EXTERNAL lsame
!    ..
!    .. External Subroutines ..
!    ..
!    .. Intrinsic Functions ..
      INTRINSIC max
!    ..
!    .. Local Scalars ..
      REAL(kind=8) TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      log10ICAL NOTA,NOTB
!    ..
!    .. Parameters ..
      REAL(kind=8) ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!    ..
!
!    Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!    transposed and set  NROWA and NROWB  as the number of rows of  A
!    and  B  respectively.
!
      if (transa.eq.'N') then
              nota=.true.
      else 
              nota=.false.
      end if 
      if (transb.eq.'N') then 
              notb=.true.
      else 
              notb=.false.
      end if

      if (nota) then
          nrowa = m
      else
          nrowa = k
      end if 
      if (notb) then
          nrowb = k
      else
          nrowb = n
      end if
!
!    Test the input parameters.
!
!      info = 0
      !IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND. &
      !   (.NOT.lsame(transa,'T'))) THEN
      !    info = 1
      !ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND. &
      !        (.NOT.lsame(transb,'T'))) THEN
      !    info = 2
      !ELSE IF (m.LT.0) THEN
      !    info = 3
      !ELSE IF (n.LT.0) THEN
      !    info = 4
      !ELSE IF (k.LT.0) THEN
      !    info = 5
      !ELSE IF (lda.LT.max(1,nrowa)) THEN
      !    info = 8
      !ELSE IF (ldb.LT.max(1,nrowb)) THEN
      !    info = 10
      !ELSE IF (ldc.LT.max(1,m)) THEN
      !    info = 13
      !END IF
      !IF (info.NE.0) THEN
      !    CALL xerbla('DGEMM ',info)
      !    RETURN
      !END IF
!
!    Quick return if possible.
!
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
         (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
!
!    And if  alpha.eq.zero.
!
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
                      Nmult = Nmult + 1
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
!
!    Start the operations.
!
      IF (notb) THEN
          IF (nota) THEN
!
!          Form  C := alpha*A*B + beta*C.
!
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
                          Nmult = Nmult + 1 
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      Nmult = Nmult + 1 
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
                          Nmult = Nmult + 1 
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
!
!          Form  C := alpha*A**T*B + beta*C
!
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
                          Nmult = Nmult + 1 
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                          Nmult = Nmult + 1
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                          Nmult = Nmult + 1 
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (nota) THEN
!
!          Form  C := alpha*A*B**T + beta*C
!
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
  130                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
                          Nmult = Nmult + 1 
  140                 CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      Nmult = Nmult + 1
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
                          Nmult = Nmult + 1 
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
!
!          Form  C := alpha*A**T*B**T + beta*C
!
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
                          Nmult = Nmult + 1 
  180                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                          Nmult = Nmult + 1 
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                          Nmult = Nmult + 1 
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
!
      RETURN
!
!    End of DGEMM
    end SUBROUTINE dgemm
end Module dgemm_module
