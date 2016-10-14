        program readmc
c
c       This subroutine reads the monte carlo data file specified
c       by file.  
c       

        implicit none
        integer i1,i2,i,j,ier,ikernel,ixvgr,nl,run,nf,nuse,
     1          nneg,l
        double precision r1,r2,temp,beta,ta,pi,bwidth
        integer nfmax,ndatamax
        parameter(ndatamax=400,nfmax=500)
        double precision gtau(ndatamax),tau(ndatamax),data(ndatamax),
     1       R(ndatamax,nfmax),Rp(ndatamax,nfmax),w(nfmax),dw(nfmax),     
     2       work(2*ndatamax),umat(ndatamax,ndatamax),
     3       cov(ndatamax,ndatamax),eigs(ndatamax),
     4       umatp(ndatamax,ndatamax),eigsp(ndatamax),
     5       Bf(nfmax,nfmax)
       character*72 linemc
        integer symques
        parameter(pi=3.1415927)
        external ssvdc
        double precision range

c       unit    file
c       5       standard input
c       6       standard output
c       92      QMC data and covariance
c       10      image default
c       9       eigenvalue output

      write(6,11) 
 11   format('   enter the input data filename:  ',$)
      read(5,10) linemc 
      open(unit=92,file=linemc,status='old')

      write(6,12) 
 12   format('  enter the output data filename:  ',$)
      read(5,10) linemc 
      open(unit=96,file=linemc,status='unknown')

      write(6,13) 
 13   format('      enter the default filename:  ',$)
      read(5,10) linemc 
      open(unit=10,file=linemc,status='old')

      write(6,14) 
 14   format('enter the eig. spectrum filename:  ',$)
      read(5,10) linemc 
      open(unit=9,file=linemc,status='unknown')
      if(index(linemc(1:72),'xvgr').ne.0) ixvgr=1

      write(6,15) 
 15   format('  enter the range of eigenvalues:  ',$)
      read(5,*) range
      if(range.le.0.0) 
     &write(6,*) 'Only the diagonal elements of cov will be used'

 16   write(6,17) 
 17   format('        enter the type of kernel:  ',$)
      read(5,*) ikernel

      if(ikernel.le.0) then
        write(6,*) '1   symmetric fermion kernel'
        write(6,*) '2   asymmetric fermion kernel'
        write(6,*) '3   symmetric boson kernel (for chipp(w)/w)'
        write(6,*) '4   symmetric boson structure factor kernel'
        write(6,*) '5   symmetric matsubara frequency kernel'
        goto 16
      end if

      symques=0
 401  read(92,10) linemc
 10   format(a72)
c
      if(index(linemc(1:72),'sym').ne.0) symques=1
      if(index(linemc(1:72),'Results').ne.0 .and.
     1   index(linemc(1:72),' nl').ne.0   ) then
         read(92,*) nl,run,beta
         temp=1.0/beta
      else
         goto 401
      endif
c
        if(symques.eq.1)then
c       the dataset is symmetric and is stored in reduced form
          nuse=nl/2+1
        else
          nuse=nl
        end if

        do i=1,10000
          read(10,*,end=122) w(i),r1,dw(i)
        end do
 122    nf=i-1

c       check for inconsistencies
        if(w(1).gt.w(nf)) then
          write(6,*) 'ERROR, w(1) > w(nf) '
          stop
        end if
        if(w(1).lt.0.0.and.ikernel.gt.2) then
          write(6,*) 'ERROR, w(1) < 0 for boson kernel '
          stop
        end if

c
c     Input the raw data
 404  read(92,10) linemc
      if(index(linemc(1:72),'Gtau(n)').ne.0) then
c         read(92,10) linemc
         do 20 i=1,nuse
           read(92,*) tau(i),gtau(i)
           gtau(i)=abs(gtau(i))
  20     continue
      else
        goto 404
      end if


c       Now form the kernel of the transform.  Recall that the
c       chi(tau) data is the -+ data, which is 1/2 of the zz data.

       if(ikernel.eq.1) then
c        The kernel should reflect the symmetry of the data
         if(w(1).lt.0.0d0)then 
           do i=1,nuse
             ta=tau(i)
             do j=1,nf
               if(w(i).ge.0.0d0) R(i,j)=
     &             0.5d0*(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     &            /(1.0d0+exp(-beta*w(j)))
               if(w(i).lt.0.0d0) R(i,j)=
     &             0.5d0*(exp(ta*w(j))+exp((beta-ta)*w(j)))
     &            /(1.0d0+exp(beta*w(j)))
             end do
           end do
         else  ! run w(i) over positive frequencies only
           do i=1,nuse
             ta=tau(i)
             do j=1,nf
               R(i,j)=(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     1            /(1.0d0+exp(-beta*w(j)))
             end do
           end do
         end if
       else if(ikernel.eq.2) then  ! asymmetric fermion kernel
         do i=1,nuse
           ta=tau(i)
           do j=1,nf
             if(w(i).ge.0.0d0) R(i,j)=exp(-ta*w(j))/
     &                             (1.0+exp(-beta*w(j)))
             if(w(i).lt.0.0d0) R(i,j)=exp((beta-ta)*w(j))/
     &                             (1.0+exp(beta*w(j)))
           end do
         end do
       else if (ikernel.eq.3) then  ! symmetric boson kernel (for chipp(w)/w)
         do i=1,nuse
           ta=tau(i)
           do j=1,nf
             if(w(j).gt.1.0d-5) then
               R(i,j)=w(j)*(exp(-ta*w(j)) + exp(-(beta-ta)*w(j)))/
     &                (1.0-exp(-beta*w(j)))
             else
               R(i,j)=temp*(exp(-ta*w(j)) + exp(-(beta-ta)*w(j)))
             end if
           end do
         end do
       else if (ikernel.eq.4) then  ! symmetric boson structure factor kernel
         do i=1,nuse
           ta=tau(i)
           do j=1,nf
             R(i,j) = exp(-ta*w(j)) + exp(-(beta-ta)*w(j))
           end do
         end do
       else if (ikernel.eq.5) then  ! symmetric matsubara freq. kernel (kinda)
         do i=1,nuse
           ta=tau(i)
           do j=1,nf
             R(i,j) = 1.0/(ta**2+w(j)**2)
           end do
         end do
       end if

c      Input the covariance
 405   read(92,10) linemc
       if(index(linemc(1:72),'Gij-GiGj').ne.0) then
c          read(92,10) linemc
          do 30 i=1,nuse
          do 30 j=i,nuse
            read(92,*) i1,i2,cov(i,j)
            cov(j,i)=cov(i,j)
  30     continue
       else
         goto 405
       end if
        if(cov(1,1).lt.0.000001) cov(1,1)=cov(2,2)

      if(range.gt.0.0) then
c       Now diagonalize cov with SSVDC
c       SUBROUTINE SSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
        call ssvdc(cov,ndatamax,nuse,nuse,eigs,eigsp,umat,ndatamax,
     1             umatp,ndatamax,work,10,ier)
        if(ier.ne.0) write(6,*) 'ssvdc, info= ',ier
c       call DGESVD( JOBU, JOBVT, M, N, cov, ndatamax, eigs, U, LDU, VT, LDVT,
c     1                   WORK, LWORK, INFO )
c       call DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
c     1                   WORK, LWORK, INFO )
c       The singular values (eigenvalues) are in eigs arranged in 
c       decending order.  The corresponding eigenvectors are the columns
c       of umat.

c       find out how many eigenvalues to discard.
c       range must be set to the numerical precision of the computer
c       or of the covariance data; whichever is desired.
        nneg=0
          do 33 i=1,nuse
             if(eigs(i).lt.range*eigs(1)) nneg=nneg+1
 33       continue

c              T                
c       Rp=Umat * R   
c       R(nuse-nneg,nf)  Rp(nuse-nneg,nf)  Umat(nuse-nneg,nuse-nneg)
c
        do i=1,nuse
        do j=1,nuse
          Umatp(i,j)=Umat(j,i)
        end do
        end do
        call mxma(Umatp,1,ndatamax,R,1,ndatamax,
     &            Rp,1,ndatamax,nuse-nneg,nuse-nneg,nf)

c                                          T
c       form the transformed data data=umat *Gtau 
c       data(nuse-nneg)  Umat(nuse-nneg,nuse-nneg)  Gtau(nuse-nneg)
c
        call mxva(Umatp,1,ndatamax,Gtau,1,data,1,nuse-nneg,nuse-nneg)

        do i=1,nuse
          eigs(i)=sqrt(abs(eigs(i))/float(run-1))
        end do
      else
c       Use only the diagonal elements of the covariance
        nneg=0
        do i=1,nuse
          data(i)=Gtau(i)
          eigs(i)=sqrt(abs(cov(i,i))/float(run-1))
          do j=1,nf
            Rp(i,j)=R(i,j)
          end do
        end do
      end if
c       
        write(96,*) nuse-nneg,nf,run
        if(ixvgr.eq.1) then
          write(9,*) '@g0 type logy'
          write(9,*) '@    default font 0'
          write(9,*) '@    yaxis  label "eigenvalue"'
          write(9,*) '@    xaxis  label "index"'
          write(9,*) '@    s0 symbol 2'
          write(9,*) '@    s0 linestyle 0'
          write(9,*) '@    s0 symbol size 0.5'
          write(9,*) '@TYPE xy'
        end if
        do 35 i=1,nuse-nneg
           write(96,*) data(i),eigs(i)
           write(9,*) i,eigs(i)
 35     continue

        do 40 i=1,nuse-nneg
        do 40 j=1,nf
          write(96,*) Rp(i,j)
 40     continue

      stop
      end


      subroutine mxva(a,ira,ica,b,irb,c,irc,na,nb)
c     This should emulate the SCILIB routine MXVA.  It should form the
c     matrix-vector product C=AB, where NA is the number of rows in A and the
c     number of elements in C, NB is the number of elements in B and columns
c     in A.  To get between two elements in the same row but different columns,
c     add IC. To get between two elements in the same column but different
c     rows, add IR. 
      
      implicit double precision (a-h,o-z)
      double precision a(*),b(*),c(*)
      
      iaddc=1
      do 100 irow=1,na
         iadda=1+irow*ira-ira
         iaddb=1
         xxx=0.
         do 50 icol=1,nb
            xxx=xxx+a(iadda)*b(iaddb)
            iadda=iadda+ica
            iaddb=iaddb+irb
 50      continue
         c(iaddc)=xxx
         iaddc=iaddc+irc
 100  continue
      
      return
      end
      
      subroutine mxma(a,ira,ica,b,irb,icb,c,irc,icc,na,nb,nc)
c     This should emulate the SCILIB routine MXMA.  It should form the
c     matrix product AB=C, where NA is the number of rows in A and C,
c     NB is the number of rows in B and columns in A, and NC is the
c     number of columns of C and B.  To get from one element to the
c     next column, add IC.  To get from one element to the next row,
c     add IR.
      
      implicit double precision (a-h,o-z)
      double precision a(*),b(*),c(*)
      
      do 110 icol=1,nc
         do 100 irow=1,na
            iadda=1+irow*ira-ira-ica
            iaddb=1+icol*icb-irb-icb
            xxx=0.
            do 50 i=1,nb
               iadda=iadda+ica
               iaddb=iaddb+irb
               xxx=xxx+a(iadda)*b(iaddb)
 50         continue
            iaddc=1+irow*irc+icol*icc-irc-icc
            c(iaddc)=xxx
 100     continue
 110  continue
      return
      end
      
      


c     This is the corrected version of the LINPAK code
      SUBROUTINE SSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
C***BEGIN PROLOGUE  SSVDC
C***DATE WRITTEN   790319   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D6
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=SINGLE PRECISION(SSVDC-S DSVDC-D CSVDC-C),
C             LINEAR ALGEBRA,MATRIX,SINGULAR VALUE DECOMPOSITION
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Perform the singular value decomposition of a real NXP
C            matrix
C***DESCRIPTION
C
C     SSVDC is a subroutine to reduce a real NxP matrix X by
C     orthogonal transformations U and V to diagonal form.  The
C     diagonal elements S(I) are the singular values of X.  The
C     columns of U are the corresponding left singular vectors,
C     and the columns of V the right singular vectors.
C
C     On Entry
C
C         X         REAL(LDX,P), where LDX .GE. N.
C                   X contains the matrix whose singular value
C                   decomposition is to be computed.  X is
C                   destroyed by SSVDC.
C
C         LDX       INTEGER
C                   LDX is the leading dimension of the array X.
C
C         N         INTEGER
C                   N is the number of rows of the matrix X.
C
C         P         INTEGER
C                   P is the number of columns of the matrix X.
C
C         LDU       INTEGER
C                   LDU is the leading dimension of the array U.
C                   (See below).
C
C         LDV       INTEGER
C                   LDV is the leading dimension of the array V.
C                   (See below).
C
C         WORK      REAL(N)
C                   work is a scratch array.
C
C         JOB       INTEGER
C                   JOB controls the computation of the singular
C                   vectors.  It has the decimal expansion AB
C                   with the following meaning
C
C                        A .EQ. 0  Do not compute the left singular
C                                  vectors.
C                        A .EQ. 1  Return the N left singular vectors
C                                  in U.
C                        A .GE. 2  Return the first MIN(N,P) singular
C                                  vectors in U.
C                        B .EQ. 0  Do not compute the right singular
C                                  vectors.
C                        B .EQ. 1  Return the right singular vectors
C                                  in V.
C
C     On Return
C
C         S         REAL(MM), where MM=MIN(N+1,P).
C                   The first MIN(N,P) entries of S contain the
C                   singular values of X arranged in descending
C                   order of magnitude.
C
C         E         REAL(P).
C                   E ordinarily contains zeros.  However, see the
C                   discussion of INFO for exceptions.
C
C         U         REAL(LDU,K), where LDU .GE. N.  If JOBA .EQ. 1, then
C                                   K .EQ. N.  If JOBA .GE. 2 , then
C                                   K .EQ. MIN(N,P).
C                   U contains the matrix of right singular vectors.
C                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
C                   or if JOBA .EQ. 2, then U may be identified with X
C                   in the subroutine call.
C
C         V         REAL(LDV,P), where LDV .GE. P.
C                   V contains the matrix of right singular vectors.
C                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
C                   then V may be identified with X in the
C                   subroutine call.
C
C         INFO      INTEGER.
C                   the singular values (and their corresponding
C                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
C                   are correct (here M=MIN(N,P)).  Thus if
C                   INFO .EQ. 0, all the singular values and their
C                   vectors are correct.  In any event, the matrix
C                   B = TRANS(U)*X*V is the bidiagonal matrix
C                   with the elements of S on its diagonal and the
C                   elements of E on its super-diagonal (TRANS(U)
C                   is the transpose of U).  Thus the singular
C                   values of X and B are the same.
C
C     LINPACK.  This version dated 03/19/79 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     ***** Uses the following functions and subprograms.
C
C     External drot
C     BLAS daxpy,ddot,dscal,dswap,dnrm2,drotg
C     Fortran ABS,dMAX1,MAX0,MIN0,MOD,SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  daxpy,ddot,dnrm2,drot,drotg,dscal,dswap
C***END PROLOGUE  SSVDC
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
c      REAL X(LDX,P),S(P),E(P),U(LDU,P),V(LDV,P),WORK(N)
      double precision X(LDX,P),S(P),E(P),U(LDU,N),V(LDV,P),WORK(N)
C
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     1        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      double precision ddot,T,R
      double precision B,C,CS,EL,EMM1,F,G,dnrm2,SCALE,SHIFT,SL,SM,SN,
     1                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
C***FIRST EXECUTABLE STATEMENT  SSVDC
      MAXIT = 300
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = dnrm2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0d0) GO TO 10
               IF (X(L,L) .NE. 0.0d0) S(L) = SIGN(S(L),X(L,L))
               CALL dscal(N-L+1,1.0d0/S(L),X(L,L),1)
               X(L,L) = 1.0d0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0d0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -ddot(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL daxpy(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = dnrm2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0d0) GO TO 80
               IF (E(LP1) .NE. 0.0d0) E(L) = SIGN(E(L),E(LP1))
               CALL dscal(P-L,1.0d0/E(L),E(LP1),1)
               E(LP1) = 1.0d0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0d0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0d0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL daxpy(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL daxpy(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0d0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0d0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0d0
  180       CONTINUE
            U(J,J) = 1.0d0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0d0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -ddot(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL daxpy(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL dscal(N-L+1,-1.0d0,U(L,L),1)
               U(L,L) = 1.0d0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0d0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0d0
  260          CONTINUE
               U(L,L) = 1.0d0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0d0) GO TO 320
               DO 310 J = LP1, P
                  T = -ddot(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL daxpy(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0d0
  330       CONTINUE
            V(L,L) = 1.0d0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
C     ......EXIT
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 400
            TEST = ABS(S(L)) + ABS(S(L+1))
            ZTEST = TEST + ABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0d0
C        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0d0
               IF (LS .NE. M) TEST = TEST + ABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + ABS(E(LS-1))
               ZTEST = TEST + ABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0d0
C           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (490,520,540,570), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0d0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL drotg(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL drot(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0d0
            DO 530 K = L, M
               T1 = S(K)
               CALL drotg(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL drot(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = dmax1(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),ABS(S(L)),
     1                    ABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0d0
            C = (SM*EMM1)**2
            SHIFT = 0.0d0
            IF (B .EQ. 0.0d0 .AND. C .EQ. 0.0d0) GO TO 550
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0d0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL drotg(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL drot(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL drotg(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     1            CALL drot(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0d0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL dscal(P,-1.0d0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
C           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     1            CALL dswap(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     1            CALL dswap(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END




