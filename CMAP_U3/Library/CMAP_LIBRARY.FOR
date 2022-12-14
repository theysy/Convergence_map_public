      SUBROUTINE PRE_STRAIN(ANG,BC,STRAN,PROPS,STATEV,DDSDDE)
C     # THIS SUBROUTINE SIMULATE PRE-STRAIN IMPOSED UNDER AN UNIAXIAL TENSION.
C     # THE SIMULATION IS COMPOSED OF 1. PRE-STRAIN AND 2. UNLOADING
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION STRESS(NTENS), STRAN(NTENS),
     1          PROPS(NPROPS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     2          DDEDDS(NTENS,NTENS), DDSDDE0(NTENS,NTENS)
      DIMENSION STRESS0(NDATA,NTENS), STRAN0(NDATA,NTENS),
     1          DSTRAN(NDATA,NTENS)
      DIMENSION DFDS(NTENS), SN(NTENS), TERM(NTENS)
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      COMMON /KCAMP2/ INDX
C     [1] BOUNDARY CONDITION
      DT=1.D0/NDATA
      STRAN0=0.D0
      ENU=PROPS(8) ! POISSON'S RATIO
      CALL LINSPACE(0.D0,BC,NDATA,STRAN0(:,1))
      STRAN0=DLOG(STRAN0+1.D0)
      STRAN0(:,2)=-ENU*STRAN0(:,1)
C     [2] ROTATE STRAIN TENSOR
      DO I=2,NDATA
          CALL ROTMAT(ANG,STRAN0(I,:),2.D0)
      END DO
C     [3] STRAIN INCREMENT
      DSTRAN=0.D0
      DO I=2,NDATA
          DSTRAN(I,:)=STRAN0(I,:)-STRAN0(I-1,:)
      END DO
C     [4] REFERENCE STRESS STATE
      SN=0.D0
      SN(1)=1.D0
      CALL ROTMAT(ANG,SN,1.D0)
      DO I=1,NTENS
          IF(SN(I).EQ.1.D0) INDX=I
      END DO
      CALL UMAT(PROPS,STATEV,SN,DSTRAN(1,:),DDSDDE0)
      DFDS=STATEV(51:50+NTENS)
C     [5] RUN SIMULATION
      STRESS=0.D0
      STATEV=0.D0
      STRAN=0.D0
      DO N=1,NDATA
          IF(STATEV(4).NE.0.D0) THEN
C              CALL INVERSE(DDSDDE,NTENS,DDEDDS)
C              DFDS=MATMUL(DDEDDS,SN)
C              DFDS=DFDS/DFDS(INDX)
              DSTRAN(N,:)=MAXVAL(DSTRAN(N,:))*DFDS
          END IF
          CALL UMAT(PROPS,STATEV,STRESS,DSTRAN(N,:),DDSDDE)
          DO I=1,NTENS
              STRAN(I)=STRAN(I)+DSTRAN(N,I)
          END DO
      END DO
C     [6] UNLOADING
      DDSDDE=DDSDDE0 ! EXCHANGE CEP->CE
      CALL INVERSE(DDSDDE,NTENS,DDEDDS)
      STRAN=STRAN-MATMUL(DDEDDS,STRESS)
      RETURN
      END SUBROUTINE PRE_STRAIN
C-----------------------------------------------------------------------
      SUBROUTINE INIT_STRESS(BC,HARD_VAR,SIG,STRAN0,EQPLAS,DDSDDE)
C     # THIS SUBROUTINE OBTAINS THE INITIAL STRESS TENSOR
C       RIGHT BEFORE PLASTIC DEFORMATION STARTS
C     # STRAIN INCREMENT IS CALULATED FROM A DEFORMED MATERIAL
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3), STRAN0(NDIM3),
     1 DDSDDE(NDIM3,NDIM3), DDEDDS(NDIM3,NDIM3)
      DIMENSION EIG_VAL(3), EIG_VEC(3,3), DSTRESS(NDIM3),
     1 DSTRAN(NDIM3), STRAN1(NDIM3)
      DIMENSION DISP(NDIM3)

C      DIMENSION STRESS0(NDIM3),ELAS(NDIM3),PLAS(NDIM3),STRAN(NDIM3)
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KCAMP2/ INDX
      PARAMETER(TOL=1.D-16)
      CALL INVERSE(DDSDDE,NTENS,DDEDDS)
      DSTRAN=MATMUL(DDEDDS,SIG)
C     [0] DETERMINE DISPLACEMENT
      DT=1.D0/NDATA
C      SIGN=1.D0
      INDX3=NDIM1+1
      VAL=SIG(1)*SIG(2)
      IF(VAL.LE.TOL.AND.SIG(INDX3).EQ.0.D0) THEN
          FLAG=1.D0
      ELSEIF(VAL.LE.TOL.AND.SIG(INDX3).NE.0.D0) THEN
          FLAG=2.D0
      ELSE
          FLAG=0.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN ! UNIAXIAL CASE
        DO I=1,NDIM3
          VAL=MAXVAL(DABS(SIG))
          IF(DABS(SIG(I)).EQ.VAL) INDX2=I
        END DO
        SIGN=DSTRAN(INDX2)/DABS(DSTRAN(INDX2))
        IF(INDX2.EQ.INDX.AND.SIGN.EQ.1.D0) BC0=2*BC
        BC1=DLOG(SIGN*BC0+1.D0)
        DISP2=BC1-STRAN0(INDX2)
        DISP2=DEXP(DISP2)-1.D0
        SCALE=DABS(DISP2*DT)/DABS(DSTRAN(INDX2))
        DSTRAN=SCALE*DSTRAN
      ELSEIF(FLAG.EQ.2.D0) THEN ! UNIAXIAL + SHEAR
C     UNIAXIAL
        DO I=1,NDIM1
          VAL=MAXVAL(DABS(SIG(1:NDIM1)))
          IF(DABS(SIG(I)).EQ.VAL) INDX2=I
        END DO
        SIGN=DSTRAN(INDX2)/DABS(DSTRAN(INDX2))
C          IF(INDX2.EQ.INDX.AND.SIGN.EQ.1.D0) BC=2*BC
        BC1=DLOG(SIGN*BC+1.D0)
        DISP2=BC1-STRAN0(INDX2)
        DISP2=DEXP(DISP2)-1.D0
        SCALE=DABS(DISP2*DT)/DABS(DSTRAN(INDX2))
        DSTRAN(1:NDIM1)=SCALE*DSTRAN(1:NDIM1)
C     SHEAR
        INDX2=NDIM1+1
        SIGN=DSTRAN(INDX2)/DABS(DSTRAN(INDX2))
        BC1=DLOG(SIGN*BC+1.D0)
        DISP2=BC1-STRAN0(INDX2)
        DISP2=DEXP(DISP2)-1.D0
        SCALE=DABS(DISP2*DT)/DABS(DSTRAN(INDX2))
        DSTRAN(INDX2)=SCALE*DSTRAN(INDX2)
      ELSE ! MULTIAXIAL CASE
        DO I=1,NDIM3
            IF(DABS(DSTRAN(I)).GE.TOL) THEN
                SIGN=DSTRAN(I)/DABS(DSTRAN(I))
            ELSE
                SIGN=0.D0
            END IF
            BC1=DLOG(SIGN*BC+1.D0)
            DISP(I)=BC1-STRAN0(I)
            DISP(I)=DEXP(DISP(I))-1.D0
            SCALE=DABS(DISP(I)*DT)/DABS(DSTRAN(I))
            IF(ISNAN(SCALE)) SCALE=0.D0
            IF(DISP(I).NE.0.D0) SIGN2=DISP(I)/DABS(DISP(I))
            DSTRAN(I)=SCALE*SIGN2*DABS(DSTRAN(I))
        END DO
      END IF
      SIG_BAR=0.D0
C     [2] FLOW STRESS
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      DSTRESS=MATMUL(DDSDDE,DSTRAN)
C..................................................................
      DO WHILE(SIG_BAR .LT. FLOW_SIG)
C     [3] UPDATE STRESS
        SIG=SIG+DSTRESS
C     [4] EQUIVALENT STRESS
        CALL EQV_SIG(HARD_VAR, SIG, SIG_BAR)
      END DO
      RETURN
      END SUBROUTINE INIT_STRESS
C-----------------------------------------------------------------------
      SUBROUTINE INIT_STRESS2(BC,HARD_VAR,SIG,STRAN0,EQPLAS,DDSDDE)
C     # THIS SUBROUTINE OBTAINS THE INITIAL STRESS TENSOR
C       RIGHT BEFORE PLASTIC DEFORMATION STARTS
C     # STRAIN INCREMENT IS CALULATED FROM AN UNDEFORMED MATERIAL
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3), STRAN0(NDIM3),
     1 DDSDDE(NDIM3,NDIM3), DDEDDS(NDIM3,NDIM3)
      DIMENSION EIG_VAL(3), EIG_VEC(3,3), DSTRESS(NDIM3),
     1 DSTRAN(NDIM3), STRAN1(NDIM3)
      DIMENSION DISP(NDIM3)

C      DIMENSION STRESS0(NDIM3),ELAS(NDIM3),PLAS(NDIM3),STRAN(NDIM3)
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KCAMP2/ INDX
      PARAMETER(TOL=1.D-16)
      CALL INVERSE(DDSDDE,NTENS,DDEDDS)
      DSTRAN=MATMUL(DDEDDS,SIG)
C     [0] MAJOR STRESS COMPONENT
      VAL=MAXVAL(DABS(SIG))
      DO I=1,NDIM3
          IF(SIG(I).EQ.VAL) THEN
              INDX2=I
              EXIT
          END IF
      END DO
C     [1] DETERMINE DISPLACEMENT
      DT=1.D0/NDATA
      SIGN=DSTRAN(INDX2)/DABS(DSTRAN(INDX2))
      IF(INDX2.EQ.INDX.AND.SIGN.EQ.1.D0) THEN
          BC0=2*BC
      ELSE
          BC0=BC
      END IF
      BC1=DLOG(SIGN*BC0+1.D0)
      DISP2=BC1-STRAN0(INDX2)
      DISP2=DEXP(DISP2)-1.D0
      SCALE=DABS(DISP2*DT)/DABS(DSTRAN(INDX2))
      DSTRAN=SCALE*DSTRAN
      SIG_BAR=0.D0
C     [2] FLOW STRESS
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      DSTRESS=MATMUL(DDSDDE,DSTRAN)
C..................................................................
      DO WHILE(SIG_BAR .LT. FLOW_SIG)
C     [3] UPDATE STRESS
        SIG=SIG+DSTRESS
C     [4] EQUIVALENT STRESS
        CALL EQV_SIG(HARD_VAR, SIG, SIG_BAR)
      END DO
      RETURN
      END SUBROUTINE INIT_STRESS2
C-----------------------------------------------------------------------
      SUBROUTINE SHEAR_STRESS(HARD_VAR, STRESS)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), STRESS(NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      TAU=STRESS(NDIM1+1)
      BETA1= 0.D0
      BETA2= 2.D0*DTAN(1.D0)
      BETA3= BETA1
      ITER=0
      ITER_MAX= 1000
      DO I=1, ITER_MAX
          STRESS(NDIM1+1)= DTAN(BETA3)
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          STRESS=STRESS/SIG_BAR
          IF(STRESS(NDIM1+1) .GT. TAU) THEN
              BETA2=BETA3
          ELSE
              BETA1=BETA3
          END IF
          BETA3=(BETA1+BETA2)/2.D0
          ITER=ITER+1
          IF(DABS(BETA1-BETA2) .LT. 1.D-10) THEN
              EXIT
          END IF
      END DO

      RETURN
      END SUBROUTINE SHEAR_STRESS
C-----------------------------------------------------------------------
      SUBROUTINE PI_PLANE(STRESS, PI_SIG)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION STRESS(NDIM3), PI_SIG(2), DEV_SIG(NDIM6)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C      PI_SIG(1)= 2.D0*STRESS(1)/DSQRT(6.D0)-STRESS(2)/DSQRT(6.D0)
C      PI_SIG(2)= STRESS(2)/DSQRT(2.D0)
      CALL DEVIATORIC(STRESS, DEV_SIG)

      PI_SIG(1)= 2.D0*DEV_SIG(1)/DSQRT(6.D0)-DEV_SIG(2)/DSQRT(6.D0)
     1          -DEV_SIG(3)/DSQRT(6.D0)
      PI_SIG(2)= DEV_SIG(2)/DSQRT(2.D0)-DEV_SIG(3)/DSQRT(2.D0)


      RETURN
      END SUBROUTINE PI_PLANE
C-----------------------------------------------------------------------
      SUBROUTINE PI_PLANE2(DEV_H, PI_SIG)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_H(NDIM6), PI_SIG(2), DEV_SIG(NDIM6)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      DEV_SIG=DEV_H

      PI_SIG(1)= 2.D0*DEV_SIG(1)/DSQRT(6.D0)-DEV_SIG(2)/DSQRT(6.D0)
     1          -DEV_SIG(3)/DSQRT(6.D0)
      PI_SIG(2)= DEV_SIG(2)/DSQRT(2.D0)-DEV_SIG(3)/DSQRT(2.D0)


      RETURN
      END SUBROUTINE PI_PLANE2
C-----------------------------------------------------------------------
      SUBROUTINE PI_PLANE3(DFDS, PI_SIG)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DFDS(NDIM3), PI_SIG(2), DEV_SIG(NDIM6)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      IF(NDIM3.EQ.3) THEN
          DEV_SIG(1:NDIM1)=DFDS(1:NDIM1)
          DEV_SIG(NDIM3)=-(DFDS(1)+DFDS(2))
          DEV_SIG(NDIM6)=DFDS(NDIM3)
      ELSE
          DEV_SIG=DFDS
      END IF

      PI_SIG(1)= 2.D0*DEV_SIG(1)/DSQRT(6.D0)-DEV_SIG(2)/DSQRT(6.D0)
     1          -DEV_SIG(3)/DSQRT(6.D0)
      PI_SIG(2)= DEV_SIG(2)/DSQRT(2.D0)-DEV_SIG(3)/DSQRT(2.D0)


      RETURN
      END SUBROUTINE PI_PLANE3
C-----------------------------------------------------------------------
      SUBROUTINE ROTMAT(ANG,TENS1,FACT)
C     # THIS SUBROUTINE ROTATES A TENSOR
C     # FACT=1: STRESS | FACT=2: STRAIN
C     RZ=[COS(Z) -SIN(Z) 0]  |  RY=[ COS(Y) 0 SIN(Y)]
C        [SIN(Z)  COS(Z) 0]  |     [     0  1     0 ]
C        [    0       0  1]  |     [-SIN(Y) 0 COS(Y)]
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION TENS1(NTENS)
      DIMENSION R(3,3), RT(3,3), TENS2(3,3)
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      PARAMETER(RADIAN=DATAN(1.D0)/45.D0)
      TH=ANG*RADIAN
C     [1] ROTATION MATRIX
      R(1,1)= DCOS(TH)
      R(1,2)=-DSIN(TH)
      R(1,3)= 0.D0
      R(2,1)= DSIN(TH)
      R(2,2)= DCOS(TH)
      R(2,3)= 0.D0
      R(3,1)= 0.D0
      R(3,2)= 0.D0
      R(3,3)= 1.D0
      RT=TRANSPOSE(R)
C     [2] TRANSFORM VECTOR INTO TENSOR
      CALL TRANS_VOIGT(TENS1,TENS2,FACT,1)
C     [3] ROTATION THE SECOND-ORDER TENSOR
      TENS2=MATMUL(RT,TENS2)
      TENS2=MATMUL(TENS2,R)
C     [4] TRANSFORM TENSOR INTO VECTOR
      CALL TRANS_VOIGT(TENS1,TENS2,FACT,2)
      RETURN
      END SUBROUTINE ROTMAT
C-----------------------------------------------------------------------
      SUBROUTINE TRANS_VOIGT(TENS1,TENS2,FACT,IOPT)
      IMPLICIT REAL*8(A-H, O-Z)
C     # THIS SUBROTUINE TRANSFORM TENSORS INTO VOIGT NOTATION.
C     S=[S11 S22 S12] OR S=[S11 S22 S33 S12 S13 S23]
C     S=[S11 S22 S12] OR S=[S11 S22 S33 S12 S13 S23]
C     # FACT=1: STRESS TENSOR | FACT=2: STRAIN TENSOR
C     # IOPT=1: VECTOR -> TENSOR | IOPT=2: TENSOR -> VECTOR
      DIMENSION TENS1(NTENS), TENS2(3,3)
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA

      IF(IOPT.EQ.1) THEN ! VECTOR -> TENSOR
        TENS2=0.D0
        IF(NTENS.EQ.3) THEN
          TENS2(1,1)=TENS1(1)
          TENS2(2,2)=TENS1(2)
          TENS2(1,2)=TENS1(3)/FACT
          TENS2(2,1)=TENS1(3)/FACT
        ELSE
          TENS2(1,1)=TENS1(1)
          TENS2(2,2)=TENS1(2)
          TENS2(3,3)=TENS1(3)
          TENS2(1,2)=TENS1(4)/FACT
          TENS2(2,1)=TENS1(4)/FACT
          TENS2(1,3)=TENS1(5)/FACT
          TENS2(3,1)=TENS1(5)/FACT
          TENS2(2,3)=TENS1(6)/FACT
          TENS2(3,2)=TENS1(6)/FACT
        END IF
      ELSEIF(IOPT.EQ.2) THEN ! TENSOR -> VECTOR
        TENS1=0.D0
        IF(NTENS.EQ.3) THEN
          TENS1(1)=TENS2(1,1)
          TENS1(2)=TENS2(2,2)
          TENS1(3)=TENS2(1,2)*FACT
        ELSE
          TENS1(1)=TENS2(1,1)
          TENS1(2)=TENS2(2,2)
          TENS1(3)=TENS2(3,3)
          TENS1(4)=TENS2(1,2)*FACT
          TENS1(5)=TENS2(1,3)*FACT
          TENS1(6)=TENS2(2,3)*FACT
        END IF
      END IF

      RETURN
      END SUBROUTINE TRANS_VOIGT
C-----------------------------------------------------------------------
      SUBROUTINE EIG(TENS1, EIG_VAL, EIG_VEC, FACT, IOPT)
C     #   THIS SUBROUTINE CALCULATES EIGEN VALUES AND
C         VECTORS OF THE 2ND-ORDER TENSOR.
C     IOPT=1: NORMAL SPACE -> PRINCIPAL SPACE
C     IOPT=2: PRINCIPAL SPACE -> NORMAL SPACE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TENS11(NDIM3), TENS2(3,3), EIG_VAL(3), EIG_VEC(3,3)
      DIMENSION R(3,3), RT(3,3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      IF(IOPT.EQ.1) THEN
          CALL TRANS_VOIGT(TENS1,TENS2,FACT,1)
          CALL JACOBI(TENS2, 3, 3, EIG_VAL, EIG_VEC)
C          CALL EIGSRT(EIG_VAL, EIG_VEC, 3, 3)
      ELSEIF(IOPT.EQ.2)THEN
          R=EIG_VEC
          RT=TRANSPOSE(R)
          CALL TRANS_VOIGT(EIG_VAL,TENS2,FACT,1)
          TENS2=MATMUL(R,TENS2)
          TENS2=MATMUL(TENS2,RT)
          CALL TRANS_VOIGT(TENS1,TENS2,FACT,2)
      END IF

      RETURN
      END SUBROUTINE EIG
C-----------------------------------------------------------------------
      SUBROUTINE JACOBI(A,N,NP,D,V)
C     #   THIS SUBROUTINE RETURNS EIGENVALUE AND EIGENVECTOR OF
C         A SYMMTERIC MATRIX.
C     #   THIS SUBROUTINE WAS COPIED FROM THE NUMERICAL RECIPIES.
      INTEGER N,NP,NROT,NMAX
      REAL*8 A(NP,NP),D(NP),V(NP,NP)
      PARAMETER (NMAX=500)
      INTEGER I,IP,IQ,J
      REAL*8 C,G,H,S,SM,T,TAU,THETA,TRESH,B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+
     * G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE

      RETURN
      END SUBROUTINE JACOBI
C-----------------------------------------------------------------------
      SUBROUTINE EIGSRT(D,V,N,NP)
C     #   THIS SUBROUTINE SORTS THE EIGENVALUES AND VECTORS CALCULATED FROM
C         SUB. JACOBI.
      INTEGER N,NP
      REAL*8 A(NP),D(NP),V(NP,NP)
      INTEGER I,J,K
      REAL*8 P
C      D=DABS(A)
      DO I=1,N-1
        K=I
        P=D(I)
        DO J=I+1,N
          IF(D(J).GE.P)THEN
            K=J
            P=D(J)
          ENDIF
        END DO
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
          END DO
        ENDIF
      END DO
      RETURN
      END SUBROUTINE EIGSRT
C-----------------------------------------------------------------------
      SUBROUTINE LINSPACE(FROM, TO, N, ARRAY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(N)

      RANGE = TO - FROM

      IF (N == 0) RETURN

      IF (N == 1) THEN
          ARRAY(1) = FROM
          RETURN
      END IF
      DO I=1, N
          ARRAY(I) = FROM + RANGE*(I-1)/(N-1)
      END DO
      END SUBROUTINE LINSPACE
C-----------------------------------------------------------------------
