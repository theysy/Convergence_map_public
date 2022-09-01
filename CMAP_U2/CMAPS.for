      INCLUDE 'LIBRARY\MML_U2SA.FOR'
      INCLUDE 'LIBRARY\CMAP_LIBRARY.FOR'
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H, O-Z)
      REAL*8, DIMENSION(:), ALLOCATABLE :: PROPS, STATEV, STRESS,
     1 STRAN, DFDS
      REAL*8, DIMENSION(:), ALLOCATABLE :: HARD_VAR, DEV_H0,
     1 DEV_SIG, DUMMY11, DUMMY12, HARD_VAR0
      REAL*8, DIMENSION(:), ALLOCATABLE :: DEV_H, DEV_HA, DFDH,
     1 DEV_EH0, DEV_EH
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DDFDDH, T
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: DDSDDE, DDFDDS, DDSDDE0
      DIMENSION PI_SIG(2), YL_SIG(2), DFDS0(2)
      CHARACTER :: FILEPATH*128, FILEPROPS*128, FILESTATV*128
C
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      PARAMETER(PI=DACOS(-1.D0), RADIAN=DATAN(1.D0)/45.D0)
2     FORMAT(1X,F8.4,',',F8.4,',',F8.4,',',F8.4)
21     FORMAT(1X,F8.4,',',F8.4,',',F8.4,',',F8.4)
22     FORMAT(1X,F8.4,',',F8.4)
3     FORMAT(1X,F8.2,',',F8.2,',',F8.2)
4     FORMAT(1X,F8.5,',',F8.5,',',F8.2)
5     FORMAT(1X,I8,',',I8,',',F8.2',',F8.2,',',F8.2,',',F8.2,',',F8.2)
C-----------------------------------------------------------------------
      STATEV_PAR=1.D0 ! 0: STATEV FROM SIMULATION | 1: READ STATEV FROM CSV FILE
      PLOT_PAR=0.D0 ! 0: PI-PLANE + CMAP | 1: PI-PLANE
C      FILEPROPS='PROPS_AA2090_HERSHEY.csv'
C      FILEPROPS='PROPS_AA2090_YLD2K.csv'
      FILEPROPS='PROPS_AA6022_YLD2K.csv'

C      FILESTATV='STATEV_AA2090_HERSHEY.CSV'
      FILESTATV='STATEV_AA6022_YLD2K_HAH20.csv'
C-----------------------------------------------------------------------
C     1.0   READ USER MATERIAL PROPERTIES
      FILEPATH='UMAT_PROPS\'
      OPEN(1, FILE=TRIM(FILEPATH)//TRIM(FILEPROPS), STATUS='UNKNOWN')
C
      I=0
      IOS=0
      DO WHILE(IOS .NE. -1)
            READ(1, *, IOSTAT=IOS)
            IF(IOS .EQ. -1) THEN
                  NPROPS=I
                  EXIT
            ELSE
                  I=I+1
            END IF
      END DO
      REWIND(1)
      ALLOCATE(PROPS(NPROPS))
      DO I=1, NPROPS
          READ(1, *) PROPS(I)
      END DO
      CLOSE(1)
C
      IF(STATEV_PAR.EQ.0.D0) THEN
            NSTATV=60
      ELSE
            NSTATV=40
      END IF
      ALLOCATE(STATEV(NSTATV))
C
      IF(PROPS(2) .NE. 4.) THEN
          NTENS=3
      ELSE
          NTENS=6
      END IF
      ALLOCATE(STRESS(NTENS),STRAN(NTENS),DDSDDE(NTENS,NTENS))
      ALLOCATE(DUMMY11(NTENS),DUMMY12(NTENS), DFDS(NTENS))
      ALLOCATE(DDFDDS(NTENS,NTENS))
C
C-----------------------------------------------------------------------
C     1.1  DEFINE DIMENSION OF MATRICE
      DUMMY11=0.D0
      DUMMY12=0.D0
      CALL UMAT(PROPS, STATEV, DUMMY11, DUMMY12, DDSDDE)
C     # REUTRN: 1. SIMULATION OPTION PARAMETERS
C               2. DIMENSIONS OF MATRICE
      ALLOCATE(HARD_VAR(NDIM7), DEV_H0(NDIM6), DEV_SIG(NDIM6))
      ALLOCATE(HARD_VAR0(NDIM7), DDSDDE0(NDIM3,NDIM3))
C-----------------------------------------------------------------------
C     1.2 STATE VARIABLES AFTER PRE-STRAIN SIMULATION
      IF(STATEV_PAR.EQ.0.D0) THEN
C     [1] PRE-STRAIN AND UNLOADING
      ANG=0.D0
      BC=6.D-2 !BOUNDARY CONDITION
      DT=1.D-3 !TIME INCREMENT
      NDATA=1.D0/DT
      CALL PRE_STRAIN(ANG,BC,STRAN,PROPS,STATEV,DDSDDE)
C     # RETURNS:  1. UNLOADED STRAIN AFTER PRE-STRAIN
C                 2. EVOLVED STATE VARIABLES
      ELSE
C     [2] READ STATE VARIABLE FROM A FILE
      FILEPATH='UMAT_STATEV\'
      OPEN(2, FILE=TRIM(FILEPATH)//TRIM(FILESTATV), STATUS='UNKNOWN')
      DO I=1, NSTATV
          READ(2, *) STATEV(I)
      END DO
      IF(STATEV(2).EQ.0.D0) THEN
          PRINT *, '#ERROR: STATEV FILE IS NOT DETECTED.'
      END IF
      CLOSE(2)
      END IF
C-----------------------------------------------------------------------
C     1.3 READ STATE VARIABLES OF HARDENING MODEL AFTER PRE-STRAIN
      HARD_VAR=0.D0
      EQPLAS=STATEV(1)
      IF(EQPLAS .EQ. 0.D0) THEN
          IF(HARD_PAR .EQ. 1.) THEN ! CHABOCHE
              HARD_VAR= 0.D0
          ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN ! YOSHIDA-UEMORI
              HARD_VAR= 0.D0
          ELSEIF(HARD_PAR .EQ. 3.) THEN ! HAH11
C         #   G-VALUES
              HARD_VAR(1:NDIM5)=1.D0
C         #   DEV_H
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM6
                  HARD_VAR(NDIM5+I)=DEV_H0(I)
              END DO
          ELSEIF(HARD_PAR .EQ. 4.) THEN ! HAH14
C         #   G-VALUES
              HARD_VAR(1:NDIM5)=1.D0
C         #   DEV_H
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM6
                  HARD_VAR(NDIM5+I)=DEV_H0(I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
C         #   G-VALUES
              HARD_VAR(1:NDIM5)= 1.D0
              CALL DEVIATORIC(STRESS, DEV_SIG)
              CALL NORM_HAH(DEV_SIG, DEV_H0)
              DO I=1, NDIM6
C         #   DEV_H
                  HARD_VAR(NDIM5+I)=DEV_H0(I)
C         #   DEV_HP
                  HARD_VAR(NDIM5+NDIM6+I)=DEV_H0(I)
C         #   DEV_HS
                  HARD_VAR(NDIM5+2*NDIM6+I)=DEV_H0(I)
              END DO
          END IF
      ELSE
          IF(HARD_PAR .EQ. 1.) THEN     ! CHABOCHE
              DO I=1, NDIM3
                  HARD_VAR(I)=STATEV(4+I)
                  HARD_VAR(NDIM3+I)=STATEV(10+I)
              END DO
          ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN ! YOSHIDA-UEMORI
              DO I=1, NDIM3
                  HARD_VAR(I)=STATEV(4+I)
                  HARD_VAR(NDIM3+I)=STATEV(10+I)
                  HARD_VAR(2*NDIM3+I)=STATEV(16+I)
                  HARD_VAR(3*NDIM3+I)=STATEV(22+I)
              END DO
              HARD_VAR(4*NDIM4+1)=STATEV(29)
              HARD_VAR(4*NDIM4+2)=STATEV(30)
          ELSEIF(HARD_PAR .EQ. 3.) THEN ! HAH11
C         #   G-VALUES
              HARD_VAR(1:NDIM5)=STATEV(6:9)
C         #   DEV_H
              DO I=1, NDIM6
                  HARD_VAR(NDIM5+I)=STATEV(9+I)
              END DO
          ELSEIF(HARD_PAR .EQ. 4.) THEN ! HAH14
C         #   G-VALUES
              HARD_VAR(1:NDIM5)=STATEV(6:11)
C         #   DEV_H
              DO I=1, NDIM6
                  HARD_VAR(NDIM5+I)=STATEV(11+I)
              END DO
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
C         #   G-VALUES
              HARD_VAR(1:NDIM5)= STATEV(6:15)
              DO I=1, NDIM6
C         #   DEV_H
                  HARD_VAR(NDIM5+I)=STATEV(15+I)
C         #   DEV_HP
                  HARD_VAR(NDIM5+NDIM6+I)=STATEV(21+I)
C         #   DEV_HS
                  HARD_VAR(NDIM5+2*NDIM6+I)=STATEV(27+I)
              END DO
          END IF
      END IF
C     [3] READ STATE VARIABLES OF DISLOCATION HARDENING MODEL
      IF(EQPLAS .EQ. 0.D0) THEN
          IF(FLOW_PAR .EQ. 7.) THEN ! RGBV
              RHO0=P10
              HARD_VAR(NDIM7-NDIM8+1)=RHO0 !RHO_F
              HARD_VAR(NDIM7-NDIM8+2)=0.D0 !RHO_R1
              HARD_VAR(NDIM7-NDIM8+3)=0.D0 !RHO_R2
              HARD_VAR(NDIM7-NDIM8+4)=0.D0 !RHO_L
              HARD_VAR(NDIM7-NDIM8+5)=RHO0 !RHO_F01
              HARD_VAR(NDIM7-NDIM8+6)=RHO0 !RHO_F02
              HARD_VAR(NDIM7-NDIM8+7)=RHO0 !RHO_TOT
              HARD_VAR(NDIM7-NDIM8+8)=0.D0 !DRDE
              HARD_VAR(NDIM7-NDIM8+9)=0.D0 !RHO_RTOT (TEMPORARY)
          END IF
      ELSE
          IF(FLOW_PAR .EQ. 7.) THEN ! RGBV
              DO I=1, NDIM8
                  HARD_VAR(NDIM7-NDIM8+I)=STATEV(39+I)
              END DO
          END IF
      END IF
C-----------------------------------------------------------------------
C     1.4     YIELD SURFACE PLOTTING
      OPEN(100, FILE='OUT\PI_PLANE.CSV', STATUS='UNKNOWN')
      OPEN(101, FILE='OUT\YLD_LOCUS.CSV', STATUS='UNKNOWN')
      OPEN(102, FILE='OUT\DEV_H.CSV', STATUS='UNKNOWN')
      OPEN(103, FILE='OUT\DEV_S.CSV', STATUS='UNKNOWN')
c      HARD_VAR(8)=1.1D0
      YL_SIG=0.D0
      PI_SIG=0.D0
      DO I=0, 360
          ANG=I*RADIAN
C       [1] DISTORTED YIELD SURFACE
          STRESS=0.D0
          STRESS(1)= DCOS(ANG)
          STRESS(2)= DSIN(ANG)
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          STRESS=STRESS/SIG_BAR
          CALL PI_PLANE(STRESS, PI_SIG)
C       [2] ORIGINAL YIELD SURFACE
          STRESS=0.D0
          STRESS(1)= DCOS(ANG)
          STRESS(2)= DSIN(ANG)
          HARD_PAR=0.D0
          CALL YLD(STRESS, SIG_BAR0)
          STRESS=STRESS/SIG_BAR0
          DUMMY12=STRESS
          CALL PI_PLANE(STRESS, YL_SIG)
          HARD_PAR=PROPS(1)
          WRITE(100, 2) PI_SIG(1), PI_SIG(2), YL_SIG(1), YL_SIG(2)
c          WRITE(101, 2) HARD_VAR(8), HARD_VAR0(8),
c     1    COS_X, DUMMY12(2)
          WRITE(103, 2) STRESS
      END DO
      ANG=0.D0*RADIAN
      STRESS=0.D0
      STRESS(1)= DCOS(ANG)
      STRESS(2)= DSIN(ANG)
      CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
      CALL GRAD(1.D0, HARD_VAR, STRESS, SIG_BAR, DFDS, DDFDDS)
      STRESS=STRESS/SIG_BAR
      CALL DEVIATORIC(STRESS, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)
      CALL PI_PLANE2(DEV_H0, PI_SIG)
      CALL PI_PLANE(STRESS, YL_SIG)

      WRITE(102, 2) PI_SIG(1), PI_SIG(2), YL_SIG(1), YL_SIG(2)
      CALL PI_PLANE2(DEV_SIG, YL_SIG)
      WRITE(103, 22) YL_SIG(1), YL_SIG(2)
      CLOSE(100)
      CLOSE(101)
      CLOSE(102)
      CLOSE(103)
C-----------------------------------------------------------------------
      OPEN(105, FILE='OUT\PI_PLANE2.CSV', STATUS='UNKNOWN')
      OPEN(106, FILE='OUT\YLD_LOCUS2.CSV', STATUS='UNKNOWN')
      OPEN(107, FILE='OUT\DEV_H2.CSV', STATUS='UNKNOWN')
      OPEN(108, FILE='OUT\DFDS.CSV', STATUS='UNKNOWN')
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      RATIO=FLOW_SIG/FLOW_SIG0
      YL_SIG=0.D0
      PI_SIG=0.D0
      HARD_PAR=PROPS(1)
      DO I=0, 720
          ANG=I/2.D0
          STRESS=0.D0
          STRESS(1)= DCOS(ANG*DATAN(1.D0)/45.D0)
          STRESS(2)= DSIN(ANG*DATAN(1.D0)/45.D0)
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          CALL GRAD(0.D0, HARD_VAR, STRESS, SIG_BAR, DFDS, DDFDDS)
          STRESS=STRESS/SIG_BAR
          DUMMY11=STRESS
          CALL PI_PLANE(STRESS, PI_SIG)
          STRESS=0.D0
          STRESS(1)= DCOS(ANG*DATAN(1.D0)/45.D0)
          STRESS(2)= DSIN(ANG*DATAN(1.D0)/45.D0)
          HARD_PAR=0.D0
          CALL YLD(STRESS, SIG_BAR0)
          STRESS=STRESS/SIG_BAR0
          DUMMY12=STRESS
          CALL PI_PLANE(STRESS, YL_SIG)
C          HARD_PAR=5.5D0
          WRITE(105, 2) PI_SIG(1), PI_SIG(2), YL_SIG(1), YL_SIG(2)
          WRITE(106, 2) DUMMY11(1), DUMMY11(2), DUMMY12(1), DUMMY12(2)
          WRITE(108, 2) DFDS(1), DFDS(2), DFDS(3)
      END DO
      ALLOCATE(DEV_H(NDIM6), DEV_HA(NDIM3), DFDH(NDIM3))
      ALLOCATE(DEV_EH(NDIM6), DEV_EH0(NDIM6))
      ALLOCATE(T(NDIM3,NDIM3), DDFDDH(NDIM3,NDIM3))
      ANG=0.D0
      STRESS=0.D0
      STRESS(1)= DCOS(ANG)
      STRESS(2)= DSIN(ANG)
      CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
      STRESS=STRESS/SIG_BAR
      CALL DEVIATORIC(STRESS, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)
      CALL DEVIATORIC_TRANS(DEV_H0, DEV_HA)
      CALL GRAD(0.D0, HARD_VAR, DEV_HA, SIG_BAR_H, DFDH, DDFDDH)
      CALL DEVIATORIC_TENS(T)
      DFDH=MATMUL(T, DFDH)
      IF(NDIM3 .EQ. 3) THEN
          DEV_EH0(1)=DFDH(1)
          DEV_EH0(2)=DFDH(2)
          DEV_EH0(3)=-(DFDH(1)+DFDH(2))
          DEV_EH0(4)=DFDH(3)
      ELSE
          DEV_EH0=DFDH
      END IF
      CALL NORM_HAH(DEV_EH0, DEV_EH)
      CALL PI_PLANE2(DEV_EH, PI_SIG)
      CALL PI_PLANE(STRESS, YL_SIG)
      WRITE(107, 2) PI_SIG(1), PI_SIG(2), YL_SIG(1), YL_SIG(2)
      CLOSE(105)
      CLOSE(106)
      CLOSE(107)
      CLOSE(108)

      HARD_PAR=PROPS(1)
      IF(PLOT_PAR.EQ.1.D0) THEN
          PRINT *, '#MESSAGE: PI-PLANE COMPLETED!'
          STOP
      END IF
C-----------------------------------------------------------------------
C     1.5     CONVERGENCE MAP
      YL_SIG=0.D0
      PI_SIG=0.D0
      OPEN(110, FILE='OUT\CMAP.CSV', STATUS='UNKNOWN')
      OPEN(120, FILE='DEBUG\DEBUG.CSV', STATUS='UNKNOWN')
      OPEN(130, FILE='OUT\COSX.CSV', STATUS='UNKNOWN')
C      OPEN(121, FILE='DEBUG\DEBUG2.CSV', STATUS='UNKNOWN')
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
C     #   STRAIN INCREMENT AND STRESS AMPLITUDE
      DO I=0, 720
      DO J=0, 360
          ANG1=FLOAT(I)/2.D0*RADIAN
          ANG2=FLOAT(J)/4.D0*RADIAN
          STRESS=0.D0
          STRESS(1)= DCOS(ANG1)*DCOS(ANG2)
          STRESS(2)= DSIN(ANG1)*DCOS(ANG2)
          STRESS(NDIM1+1)= DSIN(ANG2)
C     # OBTAIN CONVERGENCE BEHAVIOR
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          STRESS=STRESS/SIG_BAR
          CALL PI_PLANE(STRESS, PI_SIG)
          STRESS=STRESS*FLOW_SIG*(1.D0+1.D0)
          IF(HARD_PAR.EQ.1. .OR. HARD_PAR.EQ.2.) THEN
            DO K=1, NDIM3
              STRESS(K)=STRESS(K)-HARD_VAR(K)
            END DO
          END IF
          IF(I.EQ.0 .AND. J.EQ.0) THEN
              CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
              PRINT *, '-------------------------------'
              PRINT *, 'RESIDUAL:', SIG_BAR-FLOW_SIG
              PRINT *, '-------------------------------'
          END IF
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          WRITE(120, 3) FLOAT(I), FLOAT(J), SIG_BAR-FLOW_SIG
C          WRITE(121, 5) I, J, STRESS(1), STRESS(2), STRESS(NDIM1+1)
C-----------------------------------------------------------------------
          HARD_VAR0=HARD_VAR
          DDSDDE0=DDSDDE
          EQPLAS0=EQPLAS
          DFDS=0.D0
          CALL STRESS_UPDATE
     1      (HARD_VAR0,STRESS,EQPLAS0,SIG_BAR0,DFDS,DDSDDE0,ITER)
          IF(HARD_PAR.EQ.1. .OR. HARD_PAR.EQ.2.) THEN
            DO K=1, NDIM3
              STRESS(K)=STRESS(K)+HARD_VAR(K)
            END DO
          END IF
C     # STRAIN PATH CHANGE PARAMETER
          CALL COSX(HARD_VAR,STRESS,COS_X)
C     # DEVIATORIC PLANE COMPONENT
          STRESS=0.D0
          STRESS(1)= DCOS(ANG1)*DCOS(ANG2)
          STRESS(2)= DSIN(ANG1)*DCOS(ANG2)
          STRESS(NDIM1+1)= DSIN(ANG2)
          CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
          STRESS=STRESS/SIG_BAR
          CALL PI_PLANE(STRESS, PI_SIG)
          WRITE(110, 4) PI_SIG(1), PI_SIG(2), FLOAT(ITER)
          WRITE(130, 4) PI_SIG(1), PI_SIG(2), COS_X
      END DO
          PRINT *, I/720.*100, '%'
      END DO
      PRINT *, '#MESSAGE CMAP COMPLETED!!'
      CLOSE(110)
      CLOSE(120)
      CLOSE(130)
      END
C-----------------------------------------------------------------------