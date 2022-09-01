CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     THE ORIGINAL UMAT CODE WAS WRITTEN BY J.-W. LEE (2010. 11)              C
C                                                                             C
C     NAME    : SEONG-YONG, YOON                                              C
C     ADVISOR : FREDERIC BARLAT                                               C
C     E-MAIL  : theysy@postech.ac.kr                                          C
C     GIT-HUB : https://github.com/theysy                                     C
C     AFFILIATION:                                                            C
C         GRADUATE INSTITUTE OF FERROUS TECHNOLOGY (GIFT)                     C
C         POHANG UNIVERSITY OF SCIENCE AND TECHNOLOGY (POSTECH)               C
C                                                                             C
C     VERSION INFORMATION                                                     C
C         1. MML_U0: ONLY FOR PLANE STRESS (SHELL ELEMENT)                    C
C         2. MML_U1: SHELL OR SOLID ELEMENT                                   C
C         3. MML_VF: MML_V1 + FRACTURE MODEL                                  C
C         4. MML_U2: HAH20 AND CHABOCHE ARE UPDATED.                          C
C         5. MML_V2F: MML_V2 + FRACTURE MODEL                                 C
C         6. MML_U3: PRESSURE EFFECT AND HAH20 ARE UPDATED.                   C
C             DEVELOPMENT LOG:                                                C
C             -LINE SEARCH EULER BACKWARD IS UPDATED.(190702)                 C
C             -HAH20E IS UPDATED. (200318)                                    C
C             -MODIFIED YOSHIDA-UEMORI IS UPDATED. (200318)                   C
C             -PERMANENT SOFTENING OF HAH20 IS ENHANCED. (20200331)           C
C             -HAH14 IS UPDATED. (20200610)                                   C
C             -HAH20 LATENT HARDENING EQUATION IS UPDATE (20200810)           C
C             -RGBV MODEL IS UPDATED (20201102)                               C
C             -BROYDEN METHOD IS UPDATED (20210113)                           C
C                                                                             C
C     REFERENCES                                                              C
C         [1] F.BARLAT  ET AL. IJP     (2003): YLD2000-2D                     C
C         [2] F.BARLAT  ET AL. IJP     (2012): HAH11                          C
C         [3] F.BARLAT  ET AL. IJP     (2014): HAH14                          C
C         [4] J.W. LEE  ET AL. IJP     (2012): CUTTING PLANE METHOD           C
C         [5] J.W. LEE  ET AL. CMAM    (2012): CLOSEST POINT PROJECTION       C
C         [6] J.W. YOON ET AL. IJP     (2004): MULTI-STAGE EULER BACKWARD     C
C                                            : THE ANALYTICAL DERIVATIVE OF   C
C                                              YLD2000_2D MODEL               C
C         [7] F. BARLAT ET AL.  IJP     (2005): YLD2004_18P                   C
C         [8] H. ARETZ          ESAFORM (2007): SCALED FDM                    C
C         [9] T.J. PARK ET AT.  IJSS    (2012): CHABOCHE1-SIMPLER MODEL       C
C         [10]J.Y. LEE  ET AL.  IJSS    (2012): CHABOCHE2                     C
C         [11]F. YOSHIDA ET AL. IJP     (2002): YOSHIDA-UEMORI                C
C         [12]M. ORTIZ  ET AL.  IJNME   (1985): TANGENT MODULUS METHOD        C
C         [13]W.M. SCHERZINGER  CMAME   (2017): LINE SEARCH EULER BACKWARD    C
C         [14]F. BARLAT ET AL.  IJSS    (2020): HAH20H & HAH20E               C
C         [15]E. RAUCH ET AL.   MSMSE   (2011): RGBV                          C
C         [16]K. KITAYAMA ET AL.IJP     (2013): CRYSTALLOGRAPHIC RGBV         C
C         [17]S.Y. YOON ET AL.  CMAME   (2020): LINE-SEARCH FOR HARDENING     C
C         [18]C.G. BROYDEN      MOC     (1965): BROYDEN METHOD                C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C-----------------------------------------------------------------------
C     STATE VARIABLES FOR UMAT
C     STATEV(1)   : EQUIVALENT PLASTIC STRAIN
C     STATEV(2)   : EQUIVALENT STRESS
C     STATEV(3)   : YOUNG'S MODULUS
C     STATEV(4)   : YIELD CRITERION(0: ELASTIC / 1: PLASTIC)
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF CHABOCHE
C     STATEV(5:10)    : BACK-STRESS: ALPHA
C     STATEV(11:16)   : BACK-STRESS INCREMENT: DEL_ALPHA
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF YOSHIDA-UEMORI
C     STATEV(5:10)    : BACK-STRESS: ALPHA
C     STATEV(11:16)   : BACK-STRESS INCREMENT: DEL_ALPHA
C     STATEV(17:22)   : CENTER OF BOUNDING SURFACE: BETA
C     STATEV(23:28)   : CENTER OF NON-IH SURFACE  : Q
C     STATEV(29)      : RADIUS OF BOUNDING SURFACE: R
C     STATEV(30)      : RADIUS OF NON-IH SURFACE  : r: STAG_R
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH11 MODELS
C     [1] G-VALUES
C     STATEV(5)   : STRAIN-PATH CHANGE PARAMETER
C     STATEV(6) : G(1): G1
C     STATEV(7) : G(2): G2
C     STATEV(8) : G(3): G3
C     STATEV(9) : G(4): G4
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(10) : H11
C     STATEV(11) : H22
C     STATEV(12) : H33
C     STATEV(13) : H12
C     STATEV(14) : H23
C     STATEV(15) : H31
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH14 MODELS
C     [1] G-VALUES
C     STATEV(5)   : STRAIN-PATH CHANGE PARAMETER
C     STATEV(6) : G(1): G1
C     STATEV(7) : G(2): G2
C     STATEV(8) : G(3): G3
C     STATEV(9) : G(4): G4
C     STATEV(10): G(5): GS
C     STATEV(11): G(6): GL
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(12) : H11
C     STATEV(13) : H22
C     STATEV(14) : H33
C     STATEV(15) : H12
C     STATEV(16) : H23
C     STATEV(17) : H31
C-----------------------------------------------------------------------
C     # STATE VARIABLES IN HAH20H & HAH20E MODELS
C     [1] G-VALUES
C     STATEV(5)   : STRAIN-PATH CHANGE PARAMETER
C     STATEV(6) : G(1): G-    *
C     STATEV(7) : G(2): G+    *
C     STATEV(8) : G(3): G3    *
C     STATEV(9) : G(4): G4    *
C     STATEV(10): G(5): GP    *
C     STATEV(11): G(6): GS    *
C     STATEV(12): G(7): GC    *
C     STATEV(13): G(8): GL    *
C     STATEV(14): G(9): G3*   *
C     STATEV(15): G(10): GP*  *
C     [2] MICROSTRUCTURE DEVIATOR
C     STATEV(16) : H11        *
C     STATEV(17) : H22        *
C     STATEV(18) : H33        *
C     STATEV(19) : H12        *
C     STATEV(20) : H23
C     STATEV(21) : H31
C     STATEV(22) : HP11       *
C     STATEV(23) : HP22       *
C     STATEV(24) : HP33       *
C     STATEV(25) : HP12       *
C     STATEV(26) : HP23
C     STATEV(27) : HP31
C     STATEV(28) : HS11       *
C     STATEV(29) : HS22       *
C     STATEV(30) : HS33       *
C     STATEV(31) : HS12       *
C     STATEV(32) : HS23
C     STATEV(33) : HS31
C-----------------------------------------------------------------------
C     # STATE VARIABLES OF RGBV
C     STATEV(40): RHO_F
C     STATEV(41): RHO_R1(+)
C     STATEV(42): RHO_R2(-)
C     STATEV(43): RHO_L
C     STATEV(44): RHO_F01
C     STATEV(45): RHO_F02
C     STATEV(46): RHO_TOT
C     STATEV(47): DRDE
C     STATEV(48): BACK-STRESS(OPTIONAL)
C-----------------------------------------------------------------------
C     #   UMAT PARAMETERS
C     PROPS(1): OPTIONS FOR ANISOTROPIC HARDENING
C                 0- ISOTROPIC HARDENING
C                 1- CHABOCHE [9][10]
C                 2- YOSHIDA-UEMORI [11] | 2.5- MODIFIED YOSHIDA-UEMORI
C                 3- HAH11 [2]
C                 4- HAH14 [3]
C                 5- HAH20H [14] | 5.5-HAH20E [14]
C     PROPS(2): OPTIONS FOR YTELD FUNCTION MODEL
C                 1- VON_MISES
C                 2- HILL 1948
C                 3- YLD2000_2D   [1]
C                 4- YLD2004_18P  [7]
C     PROPS(3): OPTIONS FOR HARDENING INPUT: 
C                 1- SWIFT: SIGMA=P1*(P2+EQPLAS)**P3
C                 2- VOCE: SIGMA=P1-P2*EXP(-P3*EQPLAS)
C                 3- MODIFIED VOCE: SIGMA= P1 + P2*EQPLAS + P3*(1-EXP(-P4*EQPLAS))
C                 4- SWIFT+MODIFIED VOCE: P1*(P2+EQPLAS)**P3 + P4 + P5*EQPLAS + P6*(1-EXP(-P7*EQPLAS))
C                 5- HOCKETT-SHERBY: P1-(P1-P2)*EXP(-P3*EQPLAS**P4)
C                 6- COMBINED SIWFT-VOCE: P1*(P2*(P3+EQPLAS)**P4)+(1.D0-P1)*(P5-P6*DEXP(-P7*EQPLAS))
C                 7- DISLOCATION-BASED HARDENING MODEL: P1*[P2+P3*P4*P5*SQRT(RHO(EQPLAS))
C     PROPS(4): OPTIONS FOR STRESS-UPDATE ALGORITHM
C                 1- CUTTING PLANE MEHTOD     (EXPLICIT)
C                 2- CLOSEST POINT PROJECTION (FULLY-IMPLICIT)
C                 3- CLOSEST POINT PROJECTION (SEMI-IMPLICIT)
C                 4- TANGENT MODULUS METHOD   (IMPLICIT)
C                 5- BROYDEN METHOD           (FULLY-IMPLICIT)
C                 SUA_PAR=UDT_PAR+SSC_PAR
C                 EX) 1.55 = CTM + UPDATE METHOD + LINE-SEARCH
C     PROPS(5): OPTIOSN FOR DIFFERENTIATION METHOD
C                 1- ANAYLTICAL DIFFERENTIATION
C                 2- NUMERICAL DIFFERENTIATION
C-----------------------------------------------------------------------
C     #   MATERIAL PARAMETERS
C     PROPS(6) : PRESSURE-EFFECT COEFFICIENT
C     PROPS(7) : YOUNG'S MODULUS
C     PROPS(8) : POISSON'S RATIO
C-----------------------------------------------------------------------
C     #   HARDENING LAW
C     PROPS(9:18)=P(1:10): HARDENING LAW PARAMETERS
C-----------------------------------------------------------------------
C     #   ELASTIC MODULUS DEGRADATION
C     PROPS(19): EMOD_A (FOR CHORD-E) EX. EMOD_A=175700.
C     PROPS(20): EMOD_B (FOR CHORD-E) EX. EMOD_B=49.6
C-----------------------------------------------------------------------
C     ##  PROPS(1)=1: CHABOCHE
C     #   PROPS(21:29): COEFFICIENTS FOR CHABOCHE
C     #   CHABOCHE1 [9]: TABLE6
C     PROPS(21):      H1: H1
C     PROPS(22):      H2: H2
C     PROPS(23):      H3: RATIO OF ISO- TO KINEMATIC HARDENING
C     #   CHABOCHE2 [10]
C     PROPS(21):      H1: A
C     PROPS(22):      H2: B
C     PROPS(23:29):   H3:H9
C
C     ##  PROPS(1)=2: YOSHIDA-UEMORI
C     #   PROPS(3)=0: ORIGINAL | >0: MODIFIED
C     #   PROPS(21:31): COEFFICIENTS FOR YOSHIDA-UEMORI MODEL [11]
C     PROPS(21):      H1: B
C     PROPS(22):      H2: Y
C     PROPS(23):      H3: C
C     PROPS(24):      H4: m
C     PROPS(25):      H5: b
C     PROPS(26):      H6: h
C     PROPS(27):      H7: YH1: R_sat
C     PROPS(28):      H8: YH2
C     PROPS(29):      H9: YH3
C     PROPS(30):      H10:YH4
C     PROPS(31):      H11:YH5
C
C     ##  PROPS(1)=3: HAH11 MODEL
C     #   PROPS(21:27): COEFFICIENTS USED FOR HAH11 MODEL [2]
C     PROPS(21):      QM: EXPONENT OF HAH MODEL
C     PROPS(22):      K
C     PROPS(23:27):   K1~K5
C
C     ##  PROPS(1)=4: HAH14 MODEL
C     #   PROPS(21:31): K-VALUES USED FOR HAH14 MODEL [3]
C     PROPS(21):      QM: EXPONENT OF HAH MODEL
C     PROPS(22):      K
C     PROPS(23:27):   K1~K5
C     PROPS(28:31):   CROSS-LOADING / LATENT-HARDENING COEFFICIENTS
C
C     ##  PROPS(1)=4: HAH20 MODEL
C     #   PROPS(21:31): K-VALUES USED FOR HAH20 MODEL [2]
C     PROPS(21:22):   QM, PM: EXPONENT OF HAH MODEL
C     PROPS(23:27):   K1~K5
C     PROPS(28:31):   CROSS-LOADING / LATENT-HARDENING COEFFICIENTS
C-----------------------------------------------------------------------
C     ##  PROPS(2)=2: HILL1948
C     #   PROPS(32:37): COEFFICIENTS FOR HILL1948
C     PROPS(32):      F
C     PROPS(33):      G
C     PROPS(34):      H
C     PROPS(35):      L
C     PROPS(36):      M
C     PROPS(37):      N
C
C     ##  PROPS(2)=3: YLD2000_2D MODEL [1]
C     #   PROPS(32:40): COEFFICIENTS FOR YLD2000-2D
C     PROPS(32):      AM -> EXPONENT FOR YLD2000-2D MODEL
C     PROPS(33:40):   A1~A8
C
C     ##  PROPS(2)=4: YLD2004_18P MODEL
C     #   PROPS(32:50): COEFFICIENTS FOR YLD2004-18P
C     PROPS(32):      CM -> EXPONENT FOR YLD2004-18P MODEL
C     PROPS(33:41):   C1~C9
C     PROPS(42:50):   C10~C18
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.1   UMAT                                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UMAT(PROPS, STATEV, STRESS, DSTRAN, DDSDDE)
C
      IMPLICIT REAL*8(A-H, O-Z)

      CHARACTER*8 CMNAME
C
      DIMENSION PROPS(NPROPS), STATEV(NSTATV), STRESS(NTENS),
     1          DSTRAN(NTENS), DDSDDE(NTENS,NTENS)
      DIMENSION DFDS(NTENS)
C
      REAL*8, DIMENSION(:), ALLOCATABLE :: HARD_VAR, DEV_H0, DEV_SIG
C
      COMMON /KUMAT/ NPROPS, NSTATV, NTENS, NDATA
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KISO_HARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10
      COMMON /KHILL1948/ HF, HG, HH, HL, HM, HN
      COMMON /KYLD2K/ AM, A1, A2, A3, A4, A5, A6, A7, A8
      COMMON /KYLD2K4_1/CM, C1, C2, C3, C4, C5, C6, C7, C8, C9
      COMMON /KYLD2K4_2/C10,C11,C12,C13,C14,C15,C16,C17,C18
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KHAH2/ XC, XKC, XL, XKL, XI_L
      COMMON /KPRESS/ PEC0
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
C
C     1.0.    UMAT PARAMETERS
C
C     [1]     HAH MODEL PARAMETERS
      HARD_PAR=PROPS(1)
C     [2]     CONSITUTIVE MODEL PARAMETER
      YLD_PAR=PROPS(2)
C     [3]     HARDENING LAW PARAMETER
      FLOW_PAR=PROPS(3)
C     [4]     STRESS UPDATE ALGORITHM PARAMETER
      SUA_PAR=PROPS(4)
C     [5]     DIFFERENTIATION METHOD
      GRAD_PAR=PROPS(5)
      IF(YLD_PAR .NE. 3.) GRAD_PAR=2.
      IF(HARD_PAR .GE. 4.) GRAD_PAR=2.
C     [7]     DIMENSION OF ARRAY
      IF(NTENS .EQ. 3) THEN
          NDIM1=2
          NDIM2=1
      ELSE
          NDIM1=3
          NDIM2=3
      END IF
      NDIM3=NTENS
C     #   NDIM4: DIMENSION OF JACOBIAN MATRIX
      NDIM4=NTENS+1
C     #   NDIM5: DIMENSION OF G-VALUES
      IF(HARD_PAR .EQ. 3.) THEN   ! HAH11
          NDIM5= 4
      ELSEIF(HARD_PAR .EQ. 4.) THEN   ! HAH14
          NDIM5= 6
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
          NDIM5= 10
      ELSE
          NDIM5= 1
      END IF
C     #   NDIM6: SIZE OF DEVIATORIC STRESS TENSOR AND
C         MICROSTRUCTURAL DEVIATOR.
      IF(NDIM3 .EQ. 3) THEN
C     #   FOR PLANE STRESS AND STRAIN
          NDIM6=NDIM4
      ELSEIF(NDIM3 .EQ. 6) THEN
C     #   FOR 3-D STRESS AND STRAIN
          NDIM6=NDIM3
      ELSE
          WRITE(*,*) '# ERROR: CHECK THE INPUT FILE'
          STOP
      END IF
C     #   NDIM7: DIMENSION FOR HARDENING MODEL VARIABLES
      IF(HARD_PAR.EQ.1.) THEN ! CHABOCHE
C     #   HARD_VAR(1:6) = ALPHA
C     #   HARD_VAR(7:12)= DEL_ALPHA
          NDIM7= 2*NDIM3
      ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN ! YOSHIDA-UEMORI
C     #   HARD_VAR(1:6)  = ALPHA
C     #   HARD_VAR(7:12) = DEL_ALPHA
C     #   HARD_VAR(13:18)= BETA
C     #   HARD_VAR(19:24)= Q
C     #   HARD_VAR(25)   = R
C     #   HARD_VAR(26)   = r
          NDIM7= 4*NDIM3+2
      ELSEIF(HARD_PAR .EQ. 3.) THEN ! HAH11
C     #   HARD_VAR(1:4)= G
C     #   HARD_VAR(5:10)= DEV_H
          NDIM7= NDIM5+NDIM6
      ELSEIF(HARD_PAR .EQ. 4.) THEN ! HAH14
C     #   HARD_VAR(1:6)= G
C     #   HARD_VAR(7:12)= DEV_H
          NDIM7= NDIM5+NDIM6
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
C     #   HARD_VAR(1:10)= G
C     #   HARD_VAR(11:16)= DEV_H
C     #   HARD_VAR(17:22)= DEV_HP
C     #   HARD_VAR(23:28)= DEV_HS
C          NDIM7= NDIM5+NDIM6+NDIM6
          NDIM7= NDIM5+3*NDIM6
      ELSE
          NDIM7= 1
      END IF
C     THE NUMBER OF STATE VARIABLE FOR DISLOCATION HARDENING MODEL
      IF(FLOW_PAR.EQ.7) THEN ! RGBV (K. KITAYAMA ET AL. [16])
C     #   HARD_VAR(NDIM7+1)=RHO_F
C     #   HARD_VAR(NDIM7+2)=RHO_R1
C     #   HARD_VAR(NDIM7+3)=RHO_R2
C     #   HARD_VAR(NDIM7+4)=RHO_L
C     #   HARD_VAR(NDIM7+5)=RHO_F01
C     #   HARD_VAR(NDIM7+6)=RHO_F02
C     #   HARD_VAR(NDIM7+7)=RHO_RTOT
C     #   HARD_VAR(NDIM7+8)=RHO_TOT
C     #   HARD_VAR(NDIM7+9)=DRDE
          NDIM8=9
      ELSE
          NDIM8=0
      END IF
      NDIM7=NDIM7+NDIM8

      ALLOCATE(DEV_H0(NDIM6), DEV_SIG(NDIM6), HARD_VAR(NDIM7))
C     [8]     PARAMETERS FOR DEGUBING
      KNOEL=NOEL
      KNPT=NPT
      KKSTEP=KSTEP
      KKINC=KINC

C
C     1.1.    READ MATERIAL PROPERTIES FORM INPUT FILE
C
C     [1]     PRESSURE EFFECT COEFFICIENT
      PEC0= PROPS(6)
C     [2]     YOUNG'S MODULUS
      EMOD=STATEV(3)
      EMOD0=PROPS(7)
      IF(EMOD.EQ.0.D0) THEN
          EMOD=EMOD0
      END IF
C     [3]     POISSON'S RATIO
      ENU=PROPS(8)
C
C     1.2.    READ PARAMETERS
C
C     [1]     ISOTROPIC HARDENING LAW PARAMETERS
      P1=PROPS(9)
      P2=PROPS(10)
      P3=PROPS(11)
      P4=PROPS(12)
      P5=PROPS(13)
      P6=PROPS(14)
      P7=PROPS(15)
      P8=PROPS(16)
      P9=PROPS(17)
      P10=PROPS(18)

C     [2]     ELASTIC MODULUS DEGRADATION
      IF(PROPS(19) .EQ. 0.) THEN
          EMOD_A=EMOD0
          EMOD_B=0.D0
      ELSE
          EMOD_A=PROPS(19)
          EMOD_B=PROPS(20)
      END IF

C     [3]     COEFFICIENTS FOR HARDENING MODELS
      IF(HARD_PAR.EQ.1.) THEN       ! CHABOCHE
          H1=PROPS(21)
          H2=PROPS(22)
          H3=PROPS(23)
          H4=PROPS(24)
          H5=PROPS(25)
          H6=PROPS(26)
          H7=PROPS(27)
          H8=PROPS(28)
          H9=PROPS(29)
      ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN   ! YOSHIDA-UEMORI
          H1=PROPS(21)
          H2=PROPS(22)
          H3=PROPS(23)
          H4=PROPS(24)
          H5=PROPS(25)
          H6=PROPS(26)
          H7=PROPS(27)
          H8=PROPS(28)
          H9=PROPS(29)
          H10=PROPS(30)
          H11=PROPS(31)
      ELSEIF(HARD_PAR .EQ. 3.) THEN   ! HAH11
          QM =PROPS(21)
          PM =0.D0
          XK =PROPS(22)
          XK1=PROPS(23)
          XK2=PROPS(24)
          XK3=PROPS(25)
          XK4=PROPS(26)
          XK5=PROPS(27)
      ELSEIF(HARD_PAR .EQ. 4.) THEN   ! HAH14
          QM =PROPS(21)
          PM =0.D0
          XK =PROPS(22)
          XK1=PROPS(23)
          XK2=PROPS(24)
          XK3=PROPS(25)
          XK4=PROPS(26)
          XK5=PROPS(27)
          XC =PROPS(28) !S
          XKC=PROPS(29) !kS
          XL =PROPS(30) !L
          XKL=PROPS(31) !KL
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN   ! HAH20
          QM = PROPS(21)
          PM = PROPS(22)
          XK = PROPS(23)
          XK1= XK       ! XK=XK1
          XK2= PROPS(24)
          XK3= PROPS(25)
          XK4= PROPS(26)
          XK5= PROPS(27)
          XC=  PROPS(28)
          XKC= PROPS(29)
          XL=  PROPS(30)
          XKL= PROPS(31)
          XI_L=5.D-1
      END IF

C     [4]     YIELD FUNCTION PARAMETERS
      IF(YLD_PAR .EQ. 2.) THEN
C     #   HILL 1948
          HF=PROPS(32)
          HG=PROPS(33)
          HH=PROPS(34)
          HL=PROPS(35)
          HM=PROPS(36)
          HN=PROPS(37)
      ELSEIF(YLD_PAR .EQ. 3.) THEN
C     #   YLD2000_2D
          AM= PROPS(32)
          A1= PROPS(33)
          A2= PROPS(34)
          A3= PROPS(35)
          A4= PROPS(36)
          A5= PROPS(37)
          A6= PROPS(38)
          A7= PROPS(39)
          A8= PROPS(40)
      ELSEIF(YLD_PAR .EQ. 4) THEN
C     #   YLD2004_18P
          CM= PROPS(32)
          C1= PROPS(33)
          C2= PROPS(34)
          C3= PROPS(35)
          C4= PROPS(36)
          C5= PROPS(37)
          C6= PROPS(38)
          C7= PROPS(39)
          C8= PROPS(40)
          C9= PROPS(41)
          C10=PROPS(42)
          C11=PROPS(43)
          C12=PROPS(44)
          C13=PROPS(45)
          C14=PROPS(46)
          C15=PROPS(47)
          C16=PROPS(48)
          C17=PROPS(49)
          C18=PROPS(50)
      END IF
C-----------------------------------------------------------------------
C
C     1.3.    INITIALIZE THE TANGENT MODULUS
C
C     #       INITIALIZATION OF ELASTIC STIFFNESS MATRIX
C             ISOTROPIC ELASTIC STIFFNESS MATRIX FOR PLANE STRESS
      DDSDDE=0.D0
      XMU= EMOD/(2.D0*(1.D0+ENU))                 ! #SHEAR MODULUS:   μ=E/[2*(1+ν)]
      XLM= ENU*EMOD/((1.D0+ENU)*(1.D0-2.D0*ENU))  ! #LAME'S CONSTANT: λ=Eν/[(1+ν)(1-2ν)]
      XVAL=EMOD/(1.D0-ENU**2)
      IF(NDIM3 .EQ. 3) THEN
C     #   DDSDDE(1:2, 1:2)
          DO I=1, NDIM1
          DO J=1, NDIM1
              DDSDDE(I,J)= XVAL*ENU
          END DO
              DDSDDE(I,I)= XVAL*1.D0
          END DO
C     #   DDSDDE(3,3)
          DDSDDE(NDIM3,NDIM3)= 1.D0*XMU
      ELSE
C     #   DDSDDE(1:3, 1:3)
          DO I=1, NDIM1
          DO J=1, NDIM1
              DDSDDE(I,J)=XLM
          END DO
          DDSDDE(I,I)=XLM+2.D0*XMU
          END DO
C     #   DDSDDE(4:6, 4:6)
          DO I=NDIM1+1, NDIM3
              DDSDDE(I,I)=XMU
          END DO
      END IF

C
C     1.4.    TIRAL STRESS
C
C     #   THE UPDATED STRESS IS INITIALLY ASSUMED TO BE PURELY ELASTIC
C         FOR A GIVEN STRAIN INCREMENT, WHILE PLASTIC STATE VARIABLES
C         ARE KEPT CONSTANT [4].

C     [1]     TOTAL STRAIN
      DO I=1, NDIM3
C          STRAN(I)=STRAN(I)+DSTRAN(I)
      END DO

C     [2]     ELASTIC TRIAL STRESS
C     C = DDSDDE
C     SIGT(N+1) = SIG(N) + C:DSTRAN(N+1)  [EQ.(23)][4]
      STRESS= STRESS + MATMUL(DDSDDE,DSTRAN)
C
C     1.5.    READ STATE VARIABLES
C
C     [1] EQUIVALENT PLASTIC STRAIN AT A PREVIOUS TIME-STEP
      EQPLAS=STATEV(1)

C     [2] STATE VARIABLES OF HARDENING MODELS
C     #   WHEN EQPLAS=0, THE STATE VARIABLES WOULD BE INITIALIZED
      DFDS=0.D0
      HARD_VAR=0.D0
      IF(EQPLAS .EQ. 0.D0) THEN
          IF(HARD_PAR .EQ. 1.) THEN ! CHABOCHE
              HARD_VAR= 0.D0
          ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN ! YOSHIDA-UEMORI
              HARD_VAR= 0.D0
C              HARD_VAR(NDIM7)= 1.D-9
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
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5) THEN ! HAH20
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
              HARD_VAR(4*NDIM3+1)=STATEV(29)
              HARD_VAR(4*NDIM3+2)=STATEV(30)
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
          ELSEIF(FLOOR(HARD_PAR) .EQ. 5) THEN ! HAH20
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
C     [3] READ STATE VARIABLES OF DISLOCATION-BASED HARDENING MODEL
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
C
C     1.6.    INITIAL EQUIVALENT STRESS
C
C     EQUIVALENT STRESS
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              STRESS(I)=STRESS(I)-HARD_VAR(I)
          END DO
      END IF
c      CALL UPDATE_HARD_VAR(HARD_VAR, STRESS, EQPLAS, 0.D0)
      CALL EQV_SIG(HARD_VAR, STRESS, SIG_BAR)
C     RETURN:
C         SIG_BAR

C
C     1.7.    FLOW STRESS FROM HARDENING LAW
C
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      FLOW_SIG= FLOW_SIG*GP
C
C     1.8.    STRESS UPDATE ALGORITHM
C
C     [1]     YIELD CONDITION
      CALL PRESSURE_EFFECT(PEC, STRESS, FLOW_SIG, SIG_KK, BVAL, CVAL)
      SIG_BAR=SIG_BAR+PEC*CVAL*SIG_KK
      F=SIG_BAR - CVAL
      IF(F .GE. 1.D-6) THEN
          STATE=1.D0
C     [2]     PLASTIC CORREECTOR
          CALL STRESS_UPDATE(HARD_VAR,STRESS,EQPLAS,SIG_BAR,DDSDDE,ITER)
      ELSE
          STATE=0.D0
      END IF
C     [3]     GRADIENT OF YIELD SURFACE
      CALL FLOW_RULE(HARD_VAR, STRESS, EQPLAS, DFDS)
C     [4]     STRESS=STRESS+BACKSTRESS
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              STRESS(I)=STRESS(I)+HARD_VAR(I)
          END DO
      END IF

C     RETURN:
C             (1) HARD_VAR: UPDATED STATE VARIALBES OF HARDENING MODELS
C             (2) STRESS: UPDATED STRESS TENSOR
C             (3) EQPLAS: EQUIVALENT PLASTIC STRAIN
C             (4) DDSDDE: CONSISTENT TANGENT MODULUS

C
C     1.9.    UPDATE STATE VARIABLES
C
C     [1]     UPDATE EQUIVALENT PLASTIC STRAIN
      D_EQPLAS=EQPLAS-STATEV(1)
      STATEV(1)=EQPLAS
C     [2]     UPDATE EQUIVALENT STRESS
      STATEV(2)=SIG_BAR
C     [3]     YOUNG'S MODULUS DEGRADATION
      STATEV(3)=EMOD0-(EMOD0-EMOD_A)*(1.D0-DEXP(-EMOD_B*EQPLAS))
C     [4]     DEFORMATION STATE
      STATEV(4)=STATE
C     [5]     STATE VARIALBES OF HARDENING MODELS
      IF(HARD_PAR .EQ. 1.) THEN ! CHABOCHE
C         #STATEV(5:10)    : ALPHA
C         #STATEV(11:16)   : DEL_ALPHA
          DO I=1, NDIM3
              STATEV(4+I)= HARD_VAR(I)
          END DO
      ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN ! YOSHIDA-UEMORI
C         #STATEV(5:10)    : ALPHA
C         #STATEV(11:16)   : DEL_ALPHA
C         #STATEV(17:22)   : BETA
C         #STATEV(23:28)   : Q
C         #STATEV(29)      : R
C         #STATEV(30)      : STAG_R
          DO I=1, NDIM3
              STATEV(4+I)=HARD_VAR(I)
              STATEV(10+I)=HARD_VAR(NDIM3+I)
              STATEV(16+I)=HARD_VAR(2*NDIM3+I)
              STATEV(22+I)=HARD_VAR(3*NDIM3+I)
          END DO
          STATEV(29)=HARD_VAR(4*NDIM3+1)
          STATEV(30)=HARD_VAR(4*NDIM3+2)
      ELSEIF(HARD_PAR .EQ. 3.) THEN   ! HAH11
C         #STATEV(5)      : COS_X
          CALL COSX(HARD_VAR,STRESS,DEV_H0,COS_X)
          STATEV(5)= COS_X
C         #STATEV(6:9)    : G-VALUES
          DO I=1, NDIM5
              STATEV(5+I)=HARD_VAR(I)
          END DO
C         #STATEV(10:15)  : DEV_H
C         #STATEV(16:21)  : DEV_H0
          DO I=1, NDIM6
              STATEV(9+I)=HARD_VAR(NDIM5+I)
              STATEV(15+I)=DEV_H0(I)
          END DO
      ELSEIF(HARD_PAR .EQ. 4.) THEN   ! HAH14
C         #STATEV(5)      : COS_X
          CALL COSX(HARD_VAR,STRESS,DEV_H0,COS_X)
          STATEV(5)= COS_X
C          STATEV(5)=1.D0*COS_X
C         #STATEV(6:11)    : G-VALUES
          DO I=1, NDIM5
              STATEV(5+I)=HARD_VAR(I)
          END DO
C         #STATEV(12:17)  : DEV_H
C         #STATEV(18:23)  : DEV_H0
          DO I=1, NDIM6
              STATEV(11+I)=HARD_VAR(NDIM5+I)
              STATEV(17+I)=DEV_H0(I)
          END DO
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5) THEN
C         #STATEV(5)      : COS_X
          CALL COSX(HARD_VAR,STRESS,DEV_H0,COS_X)
          STATEV(5)= COS_X
C         #STATEV(6:15)   : G-VALUES
          DO I=1, NDIM5
              STATEV(5+I)=HARD_VAR(I)
          END DO
C         #STATEV(16:21)  : DEV_H
C         #STATEV(22:27)  : DEV_HP
C         #STATEV(28:33)  : DEV_HS
C         #STATEV(34:39)  : DEV_H0
          DO I=1, NDIM6
              STATEV(15+I)=HARD_VAR(NDIM5+I)
              STATEV(21+I)=HARD_VAR(NDIM5+NDIM6+I)
              STATEV(27+I)=HARD_VAR(NDIM5+2*NDIM6+I)
              STATEV(33+I)=DEV_H0(I)
          END DO
      ELSEIF(HARD_PAR .EQ. 0.) THEN
          STATEV(5:39)=0.D0
      END IF
      IF(FLOW_PAR .EQ. 7) THEN
C         # STATE VARIABLE OF RGBV #
C         STATEV(40): RHO_F
C         STATEV(41): RHO_R1
C         STATEV(42): RHO_R2
C         STATEV(43): RHO_L
C         STATEV(44): RHO_F01
C         STATEV(45): RHO_F02
C         STATEV(46): RHO_TOT
C         STATEV(47): DRDE
C         STATEV(48): RHO_RTOT
          DO I=1,NDIM8
              STATEV(39+I)=HARD_VAR(NDIM7-NDIM8+I)
          END DO
      END IF
      STATEV(49)=ITER
      STATEV(50)=D_EQPLAS
      DO I=1, NTENS
          STATEV(50+I)=DFDS(I)
      END DO
      DEALLOCATE(DEV_H0, DEV_SIG, HARD_VAR)

      RETURN
      END SUBROUTINE UMAT
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.2   STRESS UPDATE ALGORITHMS                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STRESS_UPDATE(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION HARD_VAR(NDIM7),SIG(NDIM3),CE(NDIM3,NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KITER/ MAX_ITER
C     [1] MAXIMUM ITERATION
      MAX_ITER=200
C     [2] SELECTION OF STRESS UPDATE ALGORITHM
      IF(FLOOR(SUA_PAR) .EQ. 1.) THEN
          CALL CTM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 2.) THEN
          CALL CPPM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 3.) THEN
          CALL CPPM2(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 4.) THEN
          CALL TMM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      ELSEIF(FLOOR(SUA_PAR) .EQ. 5.) THEN
          CALL BROYDEN(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      END IF

      RETURN
      END SUBROUTINE STRESS_UPDATE
C-----------------------------------------------------------------------
C
C     EULER FORWARD METHOD= CUTTING PLANE METHOD (EXPLICIT SCHEME)#
C     REF: [4] J.W. LEE ET AL. IJP  (2012): CUTTING PLANE METHOD  #
      SUBROUTINE CTM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), SIG_K(NDIM3), CE(NDIM3,NDIM3), DFDS(NDIM3),
     1          DDFDDS(NDIM3,NDIM3),CE_DFDS_CE_DGDS(NDIM3,NDIM3),
     2          CE_DFDS(NDIM3),DEL_SIG(NDIM3), D_SIG(NDIM3)
      DIMENSION HARD_VAR(NDIM7), DADE(NDIM3), DELTA(NDIM3)
      DIMENSION DGDS(NDIM3), DDGDDS(NDIM3, NDIM3), CE_DGDS(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      COMMON /KITER/ MAX_ITER
      PARAMETER(TOL=1.D-6)

C     1.0     INITIALIZATION
C     Ε_BAR(0,N+1) = Ε_BAR(N)
      EQPLAS_K = EQPLAS
C     SIG(0, N+1) = SIG_T(N+1)
      SIG_K = SIG
      DEL_EQPLAS=0.D0
      DEL_SIG=0.D0
C     KRONECKER DELTA
      DELTA= 0.D0
      NDI= NDIM3-NDIM2
      DELTA(1:NDI)= 1.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
210   CONTINUE
C     1.1     CALCULATE FLOW STRESS FROM HARDENING LAW
C     #   CALCULATE (1) FLOW_SIG (2) DHDE AT THE CURRENT SUB-STEP: ITER    #
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS_K, FLOW_SIG_K, DHDE)
      FLOW_SIG_K= FLOW_SIG_K*GP
      DHDE= DHDE*GP
C     1.2     PRESSURE-EFFECT
C         SIG_BAR(SY)= SIG_BAR_K
C         SIG_BAR(SP)= SIG_BAR_P
C         (LARGE) SIG_BAR(SY)= SIG_BAR_L
      CALL PRESSURE_EFFECT(PEC, SIG_K, FLOW_SIG_K, SIG_KK, BVAL, CVAL)
C     1.3     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE EQUIVALENT STRESS AND PLASTIC POTENTIAL
      CALL EQV_SIG(HARD_VAR, SIG_K, SIG_BAR_K)
      CALL EQV_SIG(HARD_VAR, SIG_K/BVAL, SIG_BAR_P)
C     #   CALCULATE THE FIRST AND SECOND DERIVATIVES OF THE YIELD SURFACE
      CALL GRAD(1.D0, HARD_VAR, SIG_K, SIG_BAR_K, DFDS, DDFDDS)
      CALL GRAD(1.D0, HARD_VAR, SIG_K/BVAL, SIG_BAR_P, DGDS, DDGDDS)
      DFDS= DFDS+PEC*CVAL*DELTA
      DDGDDS= DDGDDS/(BVAL*BVAL)
      SIG_BAR_L= SIG_BAR_K+PEC*CVAL*SIG_KK
      DCDE= DHDE/((1.D0-PEC*FLOW_SIG_K)**2)
      DEDL= SIG_BAR_K/SIG_BAR_L
      DFDL= (1.D0-PEC*SIG_KK)*DCDE*DEDL
C     1.4     RESIDUAL FUNCTION: YIELD CONDITION
      F_K=SIG_BAR_L-CVAL
      IF(ITER.EQ.0) R0=F_K
C     1.5     DENOMIATOR OF TANGENT STIFFNESS TENSOR
C     [1]     DFDS_C_DGDS
      DFDS_CE_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_CE_DGDS= DFDS_CE_DGDS+DFDS(I)*CE(I,J)*DGDS(J)
      END DO
      END DO
C     [2]     DENOMINATOR
      XDENOM=DFDS_CE_DGDS+DFDL
C     [3]     KINEMATIC HARDENING
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              DADE(I)= HARD_VAR(NDIM3+I)
              XDENOM=XDENOM+DFDS(I)*DADE(I)
          END DO
      END IF
C-----------------------------------------------------------------------
C     1.6     CONVERGENCE CHECK
      IF(DABS(F_K) .LT. TOL) THEN
          GOTO 215
      END IF
      IF(ITER.GT.MAX_ITER .OR. ISNAN(SIG_BAR_L)) THEN
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          ITER=MAX_ITER
          GOTO 215
      END IF
C-----------------------------------------------------------------------
C     1.7     STEP-SIZE CONTROL
      VAL=F_K/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,0.D0,F_K,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     1.8     STEP SIZE OF SOLUTION VARIABLE INCREMENTS [EQ. (28)][2]
C     [1]     PLASTIC MULTIPLIER (CONSISTENCY PARAMETER) AT A SUB-STEP
      DEL_LAMBDA=ALPHA*VAL
C     [2]     PLASTIC CORRECTOR AT A SUB-STEP
      DEL_SIG=0.D0
      DO I=1,NDIM3
      DO J=1,NDIM3
          DEL_SIG(I)= DEL_SIG(I)-CE(I,J)*DGDS(J)*DEL_LAMBDA
      END DO
      END DO
C     [3]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS=DEL_LAMBDA*DEDL
C     1.9     INTEGRATE SOLUTION VARIABLES
C     [1]     EQUIVALENT PLASTIC STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS
C     [2]     STRESS TENSOR
      SIG_K=SIG_K+DEL_SIG

C     1.10    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG_K, EQPLAS_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 210
C-------------------------------------------------------------------------
215   CONTINUE
C     1.11    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     [3]     EQUIVALENT STRESS: SIG_BAR(N+1)=SIG_BAR(K+1)
      SIG_BAR=SIG_BAR_L
C     1.12    UPDATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, D_EQPLAS)
      END IF
C     1.13    CONTINUUM TANGENT MODULUS
C     [1]     C_DFDS & C_DGDS
      CE_DFDS=MATMUL(CE,DFDS)
      CE_DGDS=MATMUL(CE,DGDS)
C     [2]     C_DFDS_C_DFDS
      CE_DFDS_CE_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE_DFDS_CE_DGDS(I,J)=CE_DFDS(I)*CE_DGDS(J)
      END DO
      END DO
C     [3]     ELASTOPLASTIC TANGENT STIFFNESS
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE(I,J)=CE(I,J)-CE_DFDS_CE_DGDS(I,J)/XDENOM
      END DO
      END DO
      RETURN
      END SUBROUTINE CTM
C------------------------------------------------------------------------------
C
C     EULER BACKWARD METHOD= CLOSETST POINT PROJECTION METHOD (FULLY IMPLICIT)
      SUBROUTINE CPPM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), SIG_K(NDIM3), CE(NDIM3,NDIM3), DFDS(NDIM3),
     1          DDFDDS(NDIM3,NDIM3), R1(NDIM3), E1(NDIM3,NDIM3),
     2          INDX(NDIM3), INDX2(NDIM4), DEL_SIG(NDIM3), 
     3          D_SIG(NDIM3), D_PLAS(NDIM3), E1_DFDS(NDIM3), 
     4          E1_DGDS_E1_DFDS(NDIM3,NDIM3)
      DIMENSION HARD_VAR(NDIM7), DADE(NDIM3)
      REAL*8    INV_CE(NDIM3,NDIM3), LU_CE(NDIM3,NDIM3),
     1          INV_E1(NDIM3,NDIM3), LU_INV_E1(NDIM3,NDIM3),
     2          JACOB(NDIM4,NDIM4), LHS(NDIM4), RHS(NDIM4),
     3          JACOBI(NDIM4,NDIM4), LU_J(NDIM4,NDIM4)

      DIMENSION DGDS(NDIM3), DDGDDS(NDIM3, NDIM3), E1_DGDS(NDIM3),
     1          DELTA(NDIM3)
      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      COMMON /KITER/ MAX_ITER
C     2.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_PLAS      = 0.D0
      D_LAMBDA    = 0.D0
      DEL_EQPLAS  = 0.D0
C     KRONECKER DELTA
      DELTA= 0.D0
      NDI= NDIM3-NDIM2
      DELTA(1:NDI)= 1.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C------------------------------------------------------------------------------
220   CONTINUE
C     2.1     CALCULATE FLOW STRESS FROM HARDENING LAW
C     #   CALCULATE (1) FLOW_SIG (2) DHDE AT THE CURRENT SUB-STEP: ITER    #
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS_K, FLOW_SIG_K, DHDE)
      FLOW_SIG_K= FLOW_SIG_K*GP
C     2.2     PRESSURE-EFFECT
C         SIG_BAR(SY)= SIG_BAR_K
C         SIG_BAR(SP)= SIG_BAR_P
C         LSIG_BAR(SY)= SIG_BAR_L
      CALL PRESSURE_EFFECT(PEC, SIG_K, FLOW_SIG_K, SIG_KK, BVAL, CVAL)
C     2.3     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE EQUIVALENT STRESS AND PLASTIC POTENTIAL
      CALL EQV_SIG(HARD_VAR, SIG_K, SIG_BAR_K)
      CALL EQV_SIG(HARD_VAR, SIG_K/BVAL, SIG_BAR_P)
C     #   CALCULATE THE FIRST AND SECOND DERIVATIVES OF THE YIELD SURFACE
      CALL GRAD(1.D0, HARD_VAR, SIG_K, SIG_BAR_K, DFDS, DDFDDS)
      CALL GRAD(1.D0, HARD_VAR, SIG_K/BVAL, SIG_BAR_P, DGDS, DDGDDS)

      DHDE= DHDE*GP
      DFDS= DFDS+PEC*CVAL*DELTA
      DDGDDS= DDGDDS/(BVAL*BVAL)
      SIG_BAR_L= SIG_BAR_K+PEC*CVAL*SIG_KK
      DCDE= DHDE/((1.D0-PEC*FLOW_SIG_K)**2)
      DEDL= SIG_BAR_K/SIG_BAR_L
      DFDL= (1.D0-PEC*SIG_KK)*DCDE*DEDL
C     2.4     RESIDUAL FUNCTIONS
C     [1]     NON-ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_LAMBDA*DGDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)*R1(I)
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     PRESSURE-DEPENDENT YIELD CONDITION
      R2= SIG_BAR_L - CVAL
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     2.5    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 225
      END IF
      IF(ITER.GT.MAX_ITER .OR. ISNAN(SIG_BAR_L)) THEN
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          ITER=MAX_ITER
          GOTO 225
      END IF
C-----------------------------------------------------------------------
C     2.6     GRADIENT OF RESIDUAL FUNCTION
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, INV_CE)
C     [2]     DR1DSP
C     #       INV_E1= INV_C+ΔΕ_BAR*[∂/∂Σ(∂F/∂Σ)]
      DO I=1, NDIM3
      DO J=1, NDIM3
          INV_E1(I,J)= INV_CE(I,J)+D_LAMBDA*DDGDDS(I,J)
      END DO
      END DO
C-----------------------------------------------------------------------
C     2.7     JACOBIAN MATRIX
      JACOB(1:NDIM3,1:NDIM3)=-INV_E1
      JACOB(NDIM4,NDIM4)=DFDL
      DO I=1, NDIM3
          JACOB(I,NDIM4)=-DGDS(I)
          JACOB(NDIM4,I)=-DFDS(I)
      END DO
C     #       INVERSE JACOBIAN MATRIX
      CALL INVERSE(JACOB, NDIM4, JACOBI)
C     2.8     SOLUTION VARIABLES: LHS & RHS
C     [1] LEFT-HAND SIDE: LHS
      DO I=1, NDIM4
          IF(I .LE. NDIM3) THEN
              LHS(I)=R1(I)
          ELSE
              LHS(I)=R2
          END IF
      END DO
C     [2] RIGHT-HAND SIDE: RHS
      RHS=MATMUL(JACOBI,LHS)
C-----------------------------------------------------------------------
C     2.9     STEP-SIZE CONTROL
      VAL=RHS(NDIM4)
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     2.10    STEP SIZE OF SOLUTION VARIABLE INCREMENTS [EQ. (28)][2]
C     [1]     PLASTIC MULTIPLIER (CONSISTENCY PARAMETER) AT A SUB-STEP
      DEL_LAMBDA=ALPHA*VAL
C     [2]     PLASTIC CORRECTOR AT A SUB-STEP
      DEL_SIG(1:NDIM3)=ALPHA*RHS(1:NDIM3)
C     [3]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS= DEL_LAMBDA*DEDL
C     2.11    INCREMENT AT A TIME-STEP
C     [1]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [2]     PLASTIC MULTIPLIER INCREMENT AT A TIME-STEP
      D_LAMBDA=D_LAMBDA+DEL_LAMBDA
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-INV_CE,D_SIG)

C     2.12    UPDATE SOLUTION VARIABLES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS

C     2.13    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG_K, EQPLAS_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 220
C------------------------------------------------------------------------------
225   CONTINUE
C     2.14    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     [3]     EQUIVALENT STRESS: SIG_BAR(N+1)=SIG_BAR(K+1)
      SIG_BAR=SIG_BAR_L
C     2.15    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, D_EQPLAS)
      END IF
C     2.16    CONSISTENT TANGENT MODULUS
C     [0]     INVERSE DR1DSP
      E1= 0.D0
      CALL INVERSE(INV_E1, NDIM3, E1)
C     [1]     DENOMINATOR
      DFDS_E1_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_DGDS=DFDS_E1_DGDS+DFDS(I)*E1(I,J)*DGDS(J)
      END DO
      END DO
      XDENOM=DFDS_E1_DGDS+DFDL
C     [2]     KINEMATIC HARDENING
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              DADE(I)= HARD_VAR(NDIM3+I)
              XDENOM=XDENOM+DFDS(I)*DADE(I)
          END DO
      END IF
C     [3]     E1_DGDS & E1_DFDS
      E1_DGDS=MATMUL(E1,DGDS)
      E1_DFDS=MATMUL(E1,DFDS)
C     [4]     NUMERATOR:E_DGDS_E_DFDS
      DO I=1, NDIM3
      DO J=1, NDIM3
          E1_DGDS_E1_DFDS(I,J)=E1_DGDS(I)*E1_DFDS(J)
      END DO
      END DO
C     [5]     ELASTO TANGENT STIFFNESS
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE(I,J)=E1(I,J)-E1_DGDS_E1_DFDS(I,J)/XDENOM
      END DO
      END DO
      
      RETURN
      END SUBROUTINE CPPM
C------------------------------------------------------------------------------
C
C     EULER BACKWARD METHOD= CLOSETST POINT PROJECTION METHOD (SEMI-IMPLICIT)
      SUBROUTINE CPPM2(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), SIG_K(NDIM3), CE(NDIM3,NDIM3), DFDS(NDIM3),
     1          DDFDDS(NDIM3,NDIM3), R1(NDIM3), E1(NDIM3,NDIM3),
     2          INDX(NDIM3), INDX2(NDIM4), DEL_SIG(NDIM3), 
     3          D_SIG(NDIM3), D_PLAS(NDIM3), E1_DFDS(NDIM3), 
     4          E1_DGDS_E1_DFDS(NDIM3,NDIM3)
      DIMENSION HARD_VAR(NDIM7), DADE(NDIM3)
      DIMENSION DGDS(NDIM3), DDGDDS(NDIM3, NDIM3), E1_DGDS(NDIM3),
     1          DELTA(NDIM3)
      REAL*8    INV_CE(NDIM3,NDIM3), LU_CE(NDIM3,NDIM3),
     1          INV_E1(NDIM3,NDIM3), LU_INV_E1(NDIM3,NDIM3),
     2          JACOB(NDIM4,NDIM4), LHS(NDIM4), RHS(NDIM4),
     3          JACOBI(NDIM4,NDIM4), LU_J(NDIM4,NDIM4)

      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      COMMON /KITER/ MAX_ITER
C     3.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_LAMBDA    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_SIG     = 0.D0
      ITER        = 0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
230   CONTINUE
C     3.1     CALCULATE FLOW STRESS FROM HARDENING LAW
C     #   CALCULATE (1) FLOW_SIG (2) DHDE AT THE CURRENT SUB-STEP: ITER    #
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS_K, FLOW_SIG_K, DHDE)
      FLOW_SIG_K= FLOW_SIG_K*GP
C     3.2     PRESSURE-EFFECT
C         SIG_BAR(SY)= SIG_BAR_K
C         SIG_BAR(SP)= SIG_BAR_P
C         LSIG_BAR(SY)= SIG_BAR_L
      CALL PRESSURE_EFFECT(PEC, SIG_K, FLOW_SIG_K, SIG_KK, BVAL, CVAL)
C     3.3     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE EQUIVALENT STRESS AND PLASTIC POTENTIAL
      CALL EQV_SIG(HARD_VAR, SIG_K, SIG_BAR_K)
      CALL EQV_SIG(HARD_VAR, SIG_K/BVAL, SIG_BAR_P)
C     #   CALCULATE THE FIRST AND SECOND DERIVATIVES OF THE YIELD SURFACE
      CALL GRAD(1.D0, HARD_VAR, SIG_K, SIG_BAR_K, DFDS, DDFDDS)
      CALL GRAD(1.D0, HARD_VAR, SIG_K/BVAL, SIG_BAR_P, DGDS, DDGDDS)

      DHDE= DHDE*GP
      DFDS= DFDS+PEC*CVAL*DELTA
      DDGDDS= DDGDDS/(BVAL*BVAL)
      SIG_BAR_L= SIG_BAR_K+PEC*CVAL*SIG_KK
      DCDE= DHDE/((1.D0-PEC*FLOW_SIG_K)**2)
      DEDL= SIG_BAR_K/SIG_BAR_L
      DFDL= (1.D0-PEC*SIG_KK)*DCDE*DEDL

C     3.4     RESIDUAL FUNCTIONS
C     [1]     NON-ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_LAMBDA*DGDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)*R1(I)
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     PRESSURE-DEPENDENT YIELD CONDITION
      R2= SIG_BAR_L - CVAL
      IF(ITER.EQ.0) R0=R2
C     3.5     GRADIENT OF RESIDUAL FUNCTION:DR1DS
C     [0]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, INV_CE)
C     [1]     DENOMINATOR
C     #       INV_E1= INV_C+ΔΕ_BAR*[∂/∂Σ(∂F/∂Σ)]
      DO I=1, NDIM3
      DO J=1, NDIM3
          INV_E1(I,J)= INV_CE(I,J)+D_LAMBDA*DDGDDS(I,J)
      END DO
      END DO
C     #       E1
      CALL INVERSE(INV_E1, NDIM3, E1)
C     #   DFDS_E1_DGDS
      DFDS_E1_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_DGDS=DFDS_E1_DGDS+DFDS(I)*E1(I,J)*DGDS(J)
      END DO
      END DO
      XDENOM=DFDS_E1_DGDS+DFDL
C     #   KINEMATIC HARDENING
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              DADE(I)= HARD_VAR(NDIM3+I)
              XDENOM=XDENOM+DFDS(I)*DADE(I)
          END DO
      END IF
C     [2]     NUMERATOR
C     #   DFDS_E1_R1
      DFDS_E1_R1=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_R1= DFDS_E1_R1 + DFDS(I)*E1(I,J)*R1(J)
      END DO
      END DO
      XNUMER= R2 - DFDS_E1_R1
C-----------------------------------------------------------------------
C     3.6    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 235
      END IF
      IF(ITER.GT.MAX_ITER .OR. ISNAN(SIG_BAR_L)) THEN
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          ITER=MAX_ITER
          GOTO 235
      END IF
C-----------------------------------------------------------------------
C     3.7     STEP-SIZE CONTROL
      VAL=XNUMER/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     3.8     STEP SIZE OF SOLUTION VARIABLE INCREMENTS [EQ. (28)][2]
C     [1]     PLASTIC MULTIPLIER (CONSISTENCY PARAMETER) AT A SUB-STEP
      DEL_LAMBDA=ALPHA*VAL
C     [2]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS= DEL_LAMBDA*DEDL
C     [3]     PLASTIC CORRECTOR INCREMENT
      DEL_SIG= MATMUL(-E1,R1+DGDS*DEL_LAMBDA)

C     3.9     INCREMENTS AT A TIME-STEP
C     [1]     EQUIVALENT PALSTIC STRAIN INCREMENT AT A TIME-STEP
      D_LAMBDA=D_LAMBDA+DEL_LAMBDA
      D_EQPLAS=D_EQPLAS+DEL_EQPLAS
C     [2]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-INV_CE,D_SIG)

C     3.10    UPDATE STATE VARIALBES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS

C     3.11    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG_K, EQPLAS_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 230
C----------------------------------------------------------------------
235   CONTINUE
C     3.12    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     [3]     EQUIVALENT STRESS: SIG_BAR(N+1)=SIG_BAR(K+1)
      SIG_BAR=SIG_BAR_L
C     3.13    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, D_EQPLAS)
      END IF
C     3.14    CONSISTENT TANGENT MODULUS
C     [1]     E1_DGDS & E1_DFDS
      E1_DGDS=MATMUL(E1,DGDS)
      E1_DFDS=MATMUL(E1,DFDS)
C     [2]     NUMERATOR:E_DGDS_E_DFDS
      DO I=1, NDIM3
      DO J=1, NDIM3
          E1_DGDS_E1_DFDS(I,J)=E1_DGDS(I)*E1_DFDS(J)
      END DO
      END DO
C     [3]     ELASTOPLASTIC TANGENT STIFFNESS METHOD
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE(I,J)=E1(I,J)-E1_DGDS_E1_DFDS(I,J)/XDENOM
      END DO
      END DO
      
      RETURN
      END SUBROUTINE CPPM2
C------------------------------------------------------------------------------
C
C     TANGENT MODULUS METHOD (SEMI-IMPLICIT)
C     REF: [12] M. ORTIZ  ET AL.  IJNME   (1985)
      SUBROUTINE TMM(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), SIG_K(NDIM3), CE(NDIM3,NDIM3), DFDS(NDIM3),
     1          DDFDDS(NDIM3,NDIM3), R1(NDIM3), E1(NDIM3,NDIM3),
     2          INDX(NDIM3), INDX2(NDIM4), DEL_SIG(NDIM3),
     3          D_SIG(NDIM3), D_PLAS(NDIM3), E1_DFDS(NDIM3),
     4          E1_DGDS_E1_DFDS(NDIM3,NDIM3), DEL_SIG0(NDIM3),
     5          DEL_SIG1(NDIM3), CE_DGDS(NDIM3)
      DIMENSION HARD_VAR(NDIM7), DADE(NDIM3)
      DIMENSION DEV_SIG(NDIM6), DEV_H0(NDIM6)
      DIMENSION DGDS(NDIM3), DDGDDS(NDIM3, NDIM3), E1_DGDS(NDIM3),
     1          DELTA(NDIM3)
      REAL*8    INV_CE(NDIM3,NDIM3), LU_CE(NDIM3,NDIM3),
     1          INV_E1(NDIM3,NDIM3), LU_INV_E1(NDIM3,NDIM3),
     2          JACOB(NDIM4,NDIM4), LHS(NDIM4), RHS(NDIM4),
     3          JACOBI(NDIM4,NDIM4), LU_J(NDIM4,NDIM4)

      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      COMMON /KITER/ MAX_ITER
C     4.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_LAMBDA    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_LAMBDA  = 0.D0
      DEL_SIG     = 0.D0
      DEL_SIG0    = 0.D0
      DEL_SIG1    = 0.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C     KRONECKER DELTA
      DELTA= 0.D0
      NDI= NDIM3-NDIM2
      DELTA(1:NDI)= 1.D0
C     NUMERICAL SCHEME PARAMETER
C     θ=1    : IMPLICIT INTEGRATION
C     θ=0    : EXPLICIT INTEGRATION
C     θ=0.5  : TRAPEZOIDAL RULE
      TH=0.5D0
C----------------------------------------------------------------------
240   CONTINUE
C     4.1     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE (1) FLOW_SIG (2) DHDE AT THE CURRENT SUB-STEP: ITER    #
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS_K, FLOW_SIG_K, DHDE)
      FLOW_SIG_K= FLOW_SIG_K*GP
C     4.2     PRESSURE-EFFECT
C         SIG_BAR(SY)= SIG_BAR_K
C         SIG_BAR(SP)= SIG_BAR_P
C         LSIG_BAR(SY)= SIG_BAR_L
      CALL PRESSURE_EFFECT(PEC, SIG_K, FLOW_SIG_K, SIG_KK, BVAL, CVAL)
C     4.3     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE EQUIVALENT STRESS AND PLASTIC POTENTIAL
      CALL EQV_SIG(HARD_VAR, SIG_K, SIG_BAR_K)
      CALL EQV_SIG(HARD_VAR, SIG_K/BVAL, SIG_BAR_P)
C     #   CALCULATE THE FIRST AND SECOND DERIVATIVES OF THE YIELD SURFACE
      CALL GRAD(1.D0, HARD_VAR, SIG_K, SIG_BAR_K, DFDS, DDFDDS)
      CALL GRAD(1.D0, HARD_VAR, SIG_K/BVAL, SIG_BAR_P, DGDS, DDGDDS)

      DHDE= DHDE*GP
      DFDS= DFDS+PEC*CVAL*DELTA
      DDGDDS= DDGDDS/(BVAL*BVAL)
      SIG_BAR_L= SIG_BAR_K+PEC*CVAL*SIG_KK
      DCDE= DHDE/((1.D0-PEC*FLOW_SIG_K)**2)
      DEDL= SIG_BAR_K/SIG_BAR_L
      DFDL= (1.D0-PEC*SIG_KK)*DCDE*DEDL
C     4.4     RESIDUAL FUNCTIONS
C     [1]     NON-ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_LAMBDA*DGDS(I) !(1,3)
          R1_NORM=R1_NORM+(TH*R1(I))**2
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     PRESSURE-DEPEDENT YIELD CONDITION
      R2= SIG_BAR_L - CVAL
      IF(ITER.EQ.0) R0=R2
C     4.5     GRADIENT OF RESIDUAL FUNCTION:DR1DS
C     [0]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, INV_CE)
C     [1]     DENOMINATOR
      DO I=1, NDIM3
      DO J=1, NDIM3
          INV_E1(I,J)= INV_CE(I,J)+D_LAMBDA*DDGDDS(I,J)
      END DO
      END DO
C     #       E1
      CALL INVERSE(INV_E1, NDIM3, E1)
C     #   DFDS_E1_DGDS
      DFDS_E1_DGDS=0.D0
      DFDS_CE_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_DGDS=DFDS_E1_DGDS+DFDS(I)*E1(I,J)*DGDS(J)
          DFDS_CE_DGDS=DFDS_CE_DGDS+DFDS(I)*CE(I,J)*DGDS(J)
      END DO
      END DO
      XDENOM=TH*DFDS_E1_DGDS+(1.D0-TH)*DFDS_CE_DGDS+DFDL
C     #   KINEMATIC HARDENING
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              DADE(I)= HARD_VAR(NDIM3+I)
              XDENOM=XDENOM+DFDS(I)*DADE(I)
          END DO
      END IF
C     [2]     NUMERATOR
C     #   DFDS_E1_R1
      DFDS_E1_R1=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_R1= DFDS_E1_R1 + DFDS(I)*E1(I,J)*R1(J)
      END DO
      END DO
      XNUMER= R2 - TH*DFDS_E1_R1
C-----------------------------------------------------------------------
C     4.6    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 245
      END IF
      IF(ITER.GT.MAX_ITER .OR. ISNAN(SIG_BAR_L)) THEN
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          ITER=MAX_ITER
          GOTO 245
      END IF
C-----------------------------------------------------------------------
C     4.7     STEP-SIZE CONTROL
      VAL=XNUMER/XDENOM
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C-----------------------------------------------------------------------
C     4.8     STEP-SIZE OF SOLUTION VARIABLES AT A SUB-STEP
C     [1]     EQUIVALENT PLASTIC STRAIN INCREMENT
C     DEL_EQPLAS(K+1/2)
      DEL_LAMBDA= ALPHA*VAL
      DEL_EQPLAS= DEL_LAMBDA*DEDL
C     [2]     PLASTIC CORRECTOR INCREMENT
      DEL_SIG1= MATMUL(-E1,R1+DGDS*DEL_LAMBDA)
      CE_DGDS=MATMUL(CE,DGDS)
      DEL_SIG0=-CE_DGDS*DEL_LAMBDA
      DEL_SIG=TH*DEL_SIG1+(1.D0-TH)*DEL_SIG0
C     4.9     INCREMENTS AT A TIME-STEP
C     [1]     EQUIVALENT PALSTIC STRAIN INCREMENT AT A TIME-STEP
      D_LAMBDA= D_LAMBDA+DEL_LAMBDA
      D_EQPLAS= D_EQPLAS+DEL_EQPLAS
C     [2]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-INV_CE,D_SIG)

C     4.10    UPDATE SOLUTION VARIALBES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS

C     4.11    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG_K, EQPLAS_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 240
C----------------------------------------------------------------------
245   CONTINUE
C     4.12    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     [3]     EQUIVALENT STRESS: SIG_BAR(N+1)=SIG_BAR(K+1)
      SIG_BAR=SIG_BAR_L
C     4.13    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, D_EQPLAS)
      END IF
C     4.14    CONSISTENT TANGENT MODULUS
C     [1]     E1_DGDS & E1_DFDS
      E1_DGDS=MATMUL(E1,DGDS)
      E1_DFDS=MATMUL(E1,DFDS)
C     [2]     NUMERATOR:E_DGDS_E_DFDS
      DO I=1, NDIM3
      DO J=1, NDIM3
          E1_DGDS_E1_DFDS(I,J)=E1_DGDS(I)*E1_DFDS(J)
      END DO
      END DO
C     [3]     ELASTOPLASTIC TANGENT STIFFNESS
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE(I,J)=E1(I,J)-E1_DGDS_E1_DFDS(I,J)/XDENOM
      END DO
      END DO

      RETURN
      END SUBROUTINE TMM
C-----------------------------------------------------------------------
C     EULER BACKWARD METHOD SOLVED WITH QUASI-NEWTON METHOD [18]
      SUBROUTINE BROYDEN(HARD_VAR,SIG,EQPLAS,SIG_BAR,CE,ITER)
C     # THIS SUBROUTINE UPDATE JACOBIAN MATRIX BASED ON THE CLASSICAL 
C       BROYDEN METHOD.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), SIG_K(NDIM3), CE(NDIM3,NDIM3), DFDS(NDIM3),
     1          DDFDDS(NDIM3,NDIM3), R1(NDIM3), E1(NDIM3,NDIM3),
     2          INDX(NDIM3), INDX2(NDIM4), DEL_SIG(NDIM3),
     3          D_SIG(NDIM3), D_PLAS(NDIM3), E1_DFDS(NDIM3),
     4          E1_DGDS_E1_DFDS(NDIM3,NDIM3), DEL_SIG0(NDIM3),
     5          DEL_SIG1(NDIM3), CE_DGDS(NDIM3)
      DIMENSION HARD_VAR(NDIM7), DADE(NDIM3)
      DIMENSION DEV_SIG(NDIM6), DEV_H0(NDIM6)
      DIMENSION DGDS(NDIM3), DDGDDS(NDIM3, NDIM3), E1_DGDS(NDIM3),
     1          DELTA(NDIM3)
      REAL*8    INV_CE(NDIM3,NDIM3), LU_CE(NDIM3,NDIM3),
     1          INV_E1(NDIM3,NDIM3), LU_INV_E1(NDIM3,NDIM3),
     2          JACOB(NDIM4,NDIM4), LHS(NDIM4), RHS(NDIM4),
     3          JACOBI(NDIM4,NDIM4), LU_J(NDIM4,NDIM4)
C     #   BROYDEN ARRAY
      REAL*8    JACOBI_LHS(NDIM4), WNUMER(NDIM4,NDIM4),
     1          DEL_JACOBI(NDIM4,NDIM4), VNUMER(NDIM4,NDIM4)

      PARAMETER(TOL=1.D-6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      COMMON /KITER/ MAX_ITER
C     5.0     INITIALIZATION
      SIG_K       = SIG
      EQPLAS_K    = EQPLAS
      D_SIG       = 0.D0
      D_EQPLAS    = 0.D0
      D_LAMBDA    = 0.D0
      D_PLAS      = 0.D0
      DEL_EQPLAS  = 0.D0
      DEL_SIG     = 0.D0
      ITER        = 0
C     KRONECKER DELTA
      DELTA= 0.D0
      NDI= NDIM3-NDIM2
      DELTA(1:NDI)= 1.D0
C     LINE-SEARCH
      ALPHA=1.D0
C     SUB_STEP
      ITER=0
C     STRESS UDPATE PARAMETER
C     SUA_PAR=SUA_PAR+UDT_PAR+SSC_PAR
      VAL1=SUA_PAR*10.D0
      VAL2=SUA_PAR*100.D0
      UDT_PAR=FLOOR(MOD(VAL1,10.D0)) ! UPDATE(INTEGRATION) PARAMETER
      SSC_PAR=MOD(VAL2,10.D0) ! STEP-SIZE CONTROL PARAMETER
C-----------------------------------------------------------------------
250   CONTINUE
C     5.1     CALCULATE FLOW STRESS FROM HARDENING LAW
C     #   CALCULATE (1) FLOW_SIG (2) DHDE AT THE CURRENT SUB-STEP: ITER    #
      IF(FLOOR(HARD_PAR) .EQ. 5.) THEN
          GP= HARD_VAR(5)
      ELSE
          GP= 1.D0
      END IF
      CALL FLOW_STRESS(HARD_VAR, EQPLAS_K, FLOW_SIG_K, DHDE)
      FLOW_SIG_K= FLOW_SIG_K*GP
C     5.2     PRESSURE-EFFECT
C         SIG_BAR(SY)= SIG_BAR_K
C         SIG_BAR(SP)= SIG_BAR_P
C         LSIG_BAR(SY)= SIG_BAR_L
      CALL PRESSURE_EFFECT(PEC, SIG_K, FLOW_SIG_K, SIG_KK, BVAL, CVAL)
C     5.3     CALCAULTE EQUIVALENT STRESS AND DERIVATIVE
C     #   CALCULATE EQUIVALENT STRESS AND PLASTIC POTENTIAL
      CALL EQV_SIG(HARD_VAR, SIG_K, SIG_BAR_K)
      CALL EQV_SIG(HARD_VAR, SIG_K/BVAL, SIG_BAR_P)
C     #   CALCULATE THE FIRST AND SECOND DERIVATIVES OF THE YIELD SURFACE
      CALL GRAD(1.D0, HARD_VAR, SIG_K, SIG_BAR_K, DFDS, DDFDDS)
      CALL GRAD(1.D0, HARD_VAR, SIG_K/BVAL, SIG_BAR_P, DGDS, DDGDDS)

      DHDE= DHDE*GP
      DFDS= DFDS+PEC*CVAL*DELTA
      DDGDDS= DDGDDS/(BVAL*BVAL)
      SIG_BAR_L= SIG_BAR_K+PEC*CVAL*SIG_KK
      DCDE= DHDE/((1.D0-PEC*FLOW_SIG_K)**2)
      DEDL= SIG_BAR_K/SIG_BAR_L
      DFDL= (1.D0-PEC*SIG_KK)*DCDE*DEDL
C     5.4     RESIDUAL FUNCTIONS
C     [1]     NON-ASSOCIATED FLOW RULE
      R1_NORM=0.D0
      DO I=1, NDIM3
          R1(I)= -D_PLAS(I) +D_LAMBDA*DGDS(I) !(1,3)
          R1_NORM=R1_NORM+R1(I)*R1(I)
      END DO
      R1_NORM=DSQRT(R1_NORM)
C     [2]     PRESSURE-DEPENDENT YIELD CONDITION
      R2= SIG_BAR_L - CVAL
      IF(ITER.EQ.0) R0=R2
C-----------------------------------------------------------------------
C     5.5     GRADIENT OF RESIDUAL FUNCTION
C     [1]     ELASTIC COMPLIANCE TENSOR
      CALL INVERSE(CE, NDIM3, INV_CE)
C     [2]     DR1DSP
      DO I=1, NDIM3
      DO J=1, NDIM3
          INV_E1(I,J)= INV_CE(I,J)+D_LAMBDA*DDGDDS(I,J)
      END DO
      END DO
C-----------------------------------------------------------------------
C     5.6    CONVERGENCE CHECK
      IF(DABS(R2).LT.TOL .AND. R1_NORM .LT. TOL) THEN
          GOTO 255
      END IF
      IF(ITER.GT.MAX_ITER .OR. ISNAN(SIG_BAR_L)) THEN
          FLAG=1.D0
      END IF
      IF(FLAG.EQ.1.D0) THEN
          ITER=MAX_ITER
          GOTO 255
      END IF
C-----------------------------------------------------------------------
C     5.7     JACOBIAN TENSOR (QUASI-NEWTON METHOD)
      IF(ITER .EQ. 0) THEN ! INITIAL JACOBIAN MATRIX
C     [1]     JACOBIAN MATRIX
      JACOB(1:NDIM3,1:NDIM3)=-INV_E1
      JACOB(NDIM4,NDIM4)=DFDL
      DO I=1, NDIM3
          JACOB(I,NDIM4)=-DGDS(I)
          JACOB(NDIM4,I)=-DFDS(I)
      END DO
C     [2]     INVERSE JACOBIAN MATRIX
      CALL INVERSE(JACOB, NDIM4, JACOBI)
      ELSE ! UPDATE JACOBIAN TENSOR BASED ON BROYDEN METHOD
      JACOBI_LHS=MATMUL(JACOBI,LHS)
      WDENOM=0.D0
      DO I=1,NDIM4
          WDENOM=WDENOM+JACOBI_LHS(I)*RHS(I)
      END DO
      DO I=1,NDIM4
      DO J=1,NDIM4
          WNUMER(I,J)=(RHS(I)-JACOBI_LHS(I))*RHS(J)
      END DO
      END DO
      WNUMER=MATMUL(WNUMER,JACOBI)
      DEL_JACOBI=WNUMER/WDENOM
C     #   UPDATE JACOBIAN MATRIX
      JACOBI=JACOBI+DEL_JACOBI
      END IF
C----------------------------------------------------------------------
C     5.8     SOLUTION VARIABLES: LHS & RHS
C     [1] LEFT-HAND SIDE: LHS
      DO I=1, NDIM4
          IF(I .LE. NDIM3) THEN
              LHS(I)=R1(I)
          ELSE
              LHS(I)=R2
          END IF
      END DO
C     [2] RIGHT-HAND SIDE: RHS
      RHS=MATMUL(JACOBI,LHS)
C----------------------------------------------------------------------
C     5.9     STEP-SIZE CONTROL
      VAL=RHS(NDIM4)
      IF(DABS(SSC_PAR-5.D0) .LT. TOL) THEN ! LINE-SEARCH
          CALL LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,VAL)
      END IF
C----------------------------------------------------------------------
C     5.10    STEP SIZE OF SOLUTION VARIABLE INCREMENTS [EQ. (28)][2]
C     [1]     PLASTIC MULTIPLIER (CONSISTENCY PARAMETER) AT A SUB-STEP
      DEL_LAMBDA=ALPHA*RHS(NDIM4)
C     [2]     PLASTIC CORRECTOR AT A SUB-STEP
      DEL_SIG(1:NDIM3)=ALPHA*RHS(1:NDIM3)
C     [3]     EQUIVALENT PLASTIC STRAIN
      DEL_EQPLAS= DEL_LAMBDA*DEDL
C     5.11    INCREMENT AT A TIME-STEP
C     [1]     PLASTIC CORRECTOR AT A TIME-STEP
      D_SIG=D_SIG+DEL_SIG
C     [2]     PLASTIC MULTIPLIER INCREMENT AT A TIME-STEP
      D_LAMBDA=D_LAMBDA+DEL_LAMBDA
C     [3]    PLASTIC STRAIN INCREMENT AT A TIME-STEP
      D_PLAS=MATMUL(-INV_CE,D_SIG)

C     5.12    UPDATE SOLUTION VARIABLES
C     [1]     STRESS
      SIG_K=SIG_K+DEL_SIG
C     [2]     EQUIVALENT STRAIN
      EQPLAS_K=EQPLAS_K+DEL_EQPLAS

C     5.13    UPDATE STATE VARIALBES OF HARDENING MODELS
      IF(DABS(UDT_PAR-0.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG_K, EQPLAS_K, DEL_EQPLAS)
      END IF

      ITER=ITER+1
      GOTO 250
C------------------------------------------------------------------------------
255   CONTINUE
C     5.14    UPDATE SOLUTION VARIABLES
C     [0]     EQUIVALENT STRAIN INCREMENT
      D_EQPLAS=EQPLAS_K-EQPLAS
C     [1]     EQUIVALENT STRAIN: E_BAR(N+1)= E_BAR(K+1)
      EQPLAS=EQPLAS_K
C     [2]     STRESS: SIG(N+1)= SIG(K+1)
      SIG=SIG_K
C     [3]     EQUIVALENT STRESS: SIG_BAR(N+1)=SIG_BAR(K+1)
      SIG_BAR=SIG_BAR_L
C     5.15    INTEGRATE STATE VARIABLES OF HARDENING MODELS II
      IF(DABS(UDT_PAR-5.D0) .LT. TOL) THEN
          CALL UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, D_EQPLAS)
      END IF
C     5.16    CONSISTENT TANGENT MODULUS
C     [0]     INVERSE DR1DSP
      E1= 0.D0
      CALL INVERSE(INV_E1, NDIM3, E1)
C     [1]     DENOMINATOR
      DFDS_E1_DGDS=0.D0
      DO I=1, NDIM3
      DO J=1, NDIM3
          DFDS_E1_DGDS=DFDS_E1_DGDS+DFDS(I)*E1(I,J)*DGDS(J)
      END DO
      END DO
      XDENOM=DFDS_E1_DGDS+DFDL
C     [2]     KINEMATIC HARDENING
      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          DO I=1, NDIM3
              DADE(I)= HARD_VAR(NDIM3+I)
              XDENOM=XDENOM+DFDS(I)*DADE(I)
          END DO
      END IF
C     [3]     E1_DGDS & E1_DFDS
      E1_DGDS=MATMUL(E1,DGDS)
      E1_DFDS=MATMUL(E1,DFDS)
C     [4]     NUMERATOR:E_DGDS_E_DFDS
      DO I=1, NDIM3
      DO J=1, NDIM3
          E1_DGDS_E1_DFDS(I,J)=E1_DGDS(I)*E1_DFDS(J)
      END DO
      END DO
C     [5]     ELASTO TANGENT STIFFNESS
      DO I=1, NDIM3
      DO J=1, NDIM3
          CE(I,J)=E1(I,J)-E1_DGDS_E1_DFDS(I,J)/XDENOM
      END DO
      END DO

      RETURN
      END SUBROUTINE BROYDEN
C-----------------------------------------------------------------------
      SUBROUTINE LINESEARCH(ALPHA,PSI0,PSIK,R1,R2,CE,ITER,DEL_EQPLAS)
C     # THIS SUBROUTINE CONTROLS THE STEP SIZE OF THE SOLUTION VARIABLE
C       BASED ON THE LINE-SEARCH TECHNIQUE.
C     REF: [13] W. M. SCHERZINGER  CMAME   (2017)
C          [17] S. Y. YOON ET AL. CMAME (2020)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION R1(NDIM3), CE(NDIM3,NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
      PARAMETER(TOL=1.D-6)
C     [0] COEFFICIENTS OF LINE-SEARCH
      BETA= 1.D-4
      ETA= 1.D-1
      EPS_MAX=1.D-4
      EM=2.D0
      GAM1= 0.4D0
      GAM2= 9.D0
      GAM3= 0.5D0
C     [1] MERIT FUNCTION: PSI  [13][EQ. (46)]
      CALL DOUBLE_DOT(R1,R1,NDIM3,DOUBLE_R1)
      XMU= CE(NDIM3,NDIM3)
      PSI= (5.D-1)*(DOUBLE_R1**(EM/2.D0)+(DABS(R2)/(2.D0*XMU))**EM)
C     [2] GOLDSTEIN'S CONDITION   [13][EQ. (51)]
      IF(ITER .EQ. 0) THEN
          PSI0= PSI
          ALPHA1= ALPHA
      ELSE
          PSI2= (1.D0-2.D0*BETA*ALPHA)*PSI0
          IF(PSI-PSI2.GT.0.D0) THEN
              ALPHA1= PSI0/(PSI0+PSI)
          ELSE
              ALPHA1= ALPHA
          END IF
      END IF
      VAL2=ALPHA*DEL_EQPLAS
C     [3] ADDITIONAL STEP SIZE CONTROL CONDITION
      IF(VAL2.LT.0.D0 .OR. VAL2.GT.EPS_MAX) THEN
C      IF(DABS(VAL2).GT.EPS_MAX) THEN
          ALPHA1= GAM1/((1.D0+GAM2*(PSI/PSI0)**GAM3))
      ELSE
          ALPHA1=1.D0
      END IF
C-----------------------------------------------------------------------
      PSIK=PSI
      ALPHA= MAX(ETA, ALPHA1)
      RETURN
      END SUBROUTINE LINESEARCH
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.3  HARDENING LAWS                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRESSURE_EFFECT(PEC, SIG, FLOW_SIG, SIG_KK, BVAL, CVAL)
C     #   THIS SUBROUTINE RETURNS SOME VALUES NEEDED IN PRESSURE-EFFECT.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KPRESS/ PEC0

C     [1] FIRST INVARIANT OF STRESS
      SIG_KK= 0.D0
      NDI= NDIM3-NDIM2
      DO I=1, NDI
          SIG_KK= SIG_KK + SIG(I)
      END DO
C      PEC=PEC0/(1.D0+DEXP(SIG_KK/FLOW_SIG))
      PEC=PEC0
C     [2] C-VALUE
      CVAL= FLOW_SIG/(1.D0-PEC*FLOW_SIG)
C     [3] B-VALUE
      BVAL= (CVAL/FLOW_SIG)*(1.D0-PEC*SIG_KK)
      RETURN
      END SUBROUTINE PRESSURE_EFFECT
C-----------------------------------------------------------------------
      SUBROUTINE FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION HARD_VAR(NDIM7)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

      IF(HARD_PAR.EQ.1. .OR. FLOOR(HARD_PAR).EQ.2.) THEN
          IF(FLOW_PAR .EQ. 7.) THEN
              WRITE(*,*) '#ERROR: RGBV IS NOT AVAILABLE IN KINEMATIC HAR
     1                    DENING'
              STOP
          END IF
      END IF

      IF(HARD_PAR.EQ.1.) THEN
          IF(H4.EQ.0.D0) THEN
C     #   REF: T.J. PARK, IJSS, 2012
              CALL ISO_HARD(HARD_VAR, 0.D0, FLOW_SIG0, DHDE)
              CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
              FLOW_SIG=H3*(FLOW_SIG+FLOW_SIG0)
              DHDE=H3*DHDE
          ELSE
              CALL CB_HARD(EQPLAS, FLOW_SIG, DHDE)
          END IF
          
      ELSEIF(FLOOR(HARD_PAR).EQ.2.) THEN
          FLOW_SIG= H2
          DHDE= 0.D0
      ELSE
          CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      END IF

      RETURN
      END SUBROUTINE FLOW_STRESS
C-----------------------------------------------------------------------
      SUBROUTINE ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS OF ISOTROPIC HARDENING LAW.
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION HARD_VAR(NDIM7)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KISO_HARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10

      IF(FLOW_PAR.EQ.1.) THEN
C     #   SWIFT HARDENING LAW #
          FLOW_SIG=P1*(P2+EQPLAS)**P3             ! H(Ε_BAR) = P1*(P2 + Ε_BAR)^P3
          DHDE    =P1*P3*(P2+EQPLAS)**(P3-1.D0)    ! ∂H/∂(Ε_BAR) = P1*P3*(P2 + Ε_BAR)^(P3-1)
      ELSEIF(FLOW_PAR.EQ.2.) THEN
C     #   VOCE HARDENING LAW #
          FLOW_SIG=P1-P2*DEXP(-P3*EQPLAS)    ! H(Ε_BAR) = P1-P2*EXP(-P3*Ε_BAR)]
          DHDE    =P2*P3*DEXP(-P3*EQPLAS)    ! ∂H/∂(Ε_BAR) = P2*P3*EXP(-P3*Ε_BAR)
      ELSEIF(FLOW_PAR.EQ.3.) THEN
C     #   MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=P1+P2*EQPLAS +P3*(1.D0-DEXP(-P4*EQPLAS))
          DHDE    =P2+ P3*P4*DEXP(-P4*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.4.) THEN
C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
          FLOW_SIG=P1*(P2+EQPLAS)**P3+P4+P5*EQPLAS
     1             +P6*(1.D0-DEXP(-P7*EQPLAS))
          DHDE    =P5+P1*P3*(P2 + EQPLAS)**(P3-1.D0)
     1             +P6*P7*DEXP(-P7*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.5.) THEN
C     #   HOCKETT-SHERBY HARDENING LAW #
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 1.D-9
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=P1-(P1-P2)*DEXP(-P3*EQPLAS0**P4)
          DHDE=(P1-P2)*P3*P4*EQPLAS0**(P4-1.D0)
     1        *DEXP(-P3*EQPLAS0**P4)
      ELSEIF(FLOW_PAR.EQ.6.) THEN
C     #   COMBINED SWIFT-VOCE #
          FLOW_SIG= P1*(P2*(P3+EQPLAS)**P4)
     1             +(1.D0-P1)*(P5-P6*DEXP(-P7*EQPLAS))
          DHDE= P1*P2*P4*(P3+EQPLAS)**(P4-1.D0)
     1         +(1.D0-P1)*P6*P7*DEXP(-P7*EQPLAS)
      ELSEIF(FLOW_PAR.EQ.7.) THEN
C     #   DISLOCATION-BASED HARDENING: RGBV [15]#
          IF(EQPLAS .EQ. 0.D0) THEN
              RHO_TOT=P10
              DRDE=0.D0
          ELSE
              RHO_TOT=HARD_VAR(NDIM7-NDIM8+7)
              DRDE=DABS(HARD_VAR(NDIM7-NDIM8+8))
          END IF
C     #   RESOLVED SHEAR STRESS
          XAL=0.5D0 ! DISLOCATION INTERACTION
          TAU=P2+XAL*P3*P4*DSQRT(RHO_TOT)
C     #   WEI WEN BACK-STRESS
C          FLOW_SIG=P1*(TAU*(1.D0-P10*XX))
C          DHDE=0.5D0*P1*P3*P4*P5/DSQRT(RHO_TOT)*DRDE
C     #   ORIGINAL EQUATION
          FLOW_SIG=P1*(TAU)
          DHDE=0.5D0*P1*XAL*P3*P4/DSQRT(RHO_TOT)*DRDE
      END IF

      RETURN
      END SUBROUTINE ISO_HARD
C-----------------------------------------------------------------------
      SUBROUTINE CB_HARD(EQPLAS, FLOW_SIG, DHDE)
C     #   THIS SUBROUTINE RETURNS STATE VARIABLES OF KINEMATIC HARDENING MODELS.
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   REF: J.Y. LEE, IJSS, 2012
C     #   SWIFT HARDENING LAW #
      IF(FLOW_PAR.EQ.1.) THEN
          FLOW_SIG=H3*(H4+EQPLAS)**H5             ! H(Ε_BAR) = H3*(H4 + Ε_BAR)^H5
          DHDE    =H3*H5*(H4+EQPLAS)**(H5-1.D0)    ! ∂H/∂(Ε_BAR) = H3*H5*(H4 + Ε_BAR)^(H5-1)
C     #   VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.2.) THEN
          FLOW_SIG=H3-H4*DEXP(-H5*EQPLAS)    ! H(Ε_BAR) = H3-H4*EXP(-H5*Ε_BAR)]
          DHDE    =H4*H5*DEXP(-H5*EQPLAS)    ! ∂H/∂(Ε_BAR) = H4*H5*EXP(-H5*Ε_BAR)
C     #   MODIFIED-VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.3.) THEN
          FLOW_SIG=H3+H4*EQPLAS +H5*(1.D0-DEXP(-H6*EQPLAS))
          DHDE    =H4+ H5*H6*DEXP(-H6*EQPLAS)
C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.4.) THEN
          FLOW_SIG=H3*(H4+EQPLAS)**H5+H6+H7*EQPLAS
     1             +H8*(1.D0-DEXP(-H9*EQPLAS))
          DHDE    =H7+H3*H5*(H4 + EQPLAS)**(H5-1.D0)
     1             +H8*H9*DEXP(-H9*EQPLAS)
C     #   HOCKETT-SHERBY HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.5.) THEN
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 1.D-9
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=H3-(H3-H4)*DEXP(-H5*EQPLAS0**H6)
          DHDE=(H3-H4)*H5*H6*EQPLAS0**(H6-1.D0)
     1            *DEXP(-H5*EQPLAS0**H6)
C     #   COMBINED SWIFT-VOCE #
      ELSEIF(FLOW_PAR.EQ.6.) THEN
          FLOW_SIG= H1*(H2*(H3+EQPLAS)**H4)
     1             +(1.D0-H1)*(H5-H6*DEXP(-H7*EQPLAS))
          DHDE= H1*H2*H4*(H3+EQPLAS)**(H4-1.D0)
     1             +(1.D0-H1)*H6*H7*DEXP(-H7*EQPLAS)
      END IF
      RETURN
      END SUBROUTINE CB_HARD
C-----------------------------------------------------------------------
      SUBROUTINE YU_HARD(EQPLAS, FLOW_SIG, DHDE)
C     #   THIS SUBROUTINE RETURNS FLOW STRESS OF ISOTROPIC HARDENING TERM OF YU.
      IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11
C     #   ORITINAL HARDENING LAW#
      IF(FLOW_PAR.EQ.0.) THEN
          FLOW_SIG= H7*(1.D0-DEXP(-H4*EQPLAS))
          DHDE    = H7*H4*DEXP(-H4*EQPLAS)
C     #   SWIFT HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.1.) THEN
          FLOW_SIG=H7*((H8+EQPLAS)**H9-H8**H9)             ! H(?_BAR) = H7*(H8 + ?_BAR)^H9
          DHDE    =H7*H9*(H8+EQPLAS)**(H9-1.D0)    ! ?H/?(?_BAR) = H7*H9*(H8 + ?_BAR)^(H9-1)

C     #   VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.2.) THEN
          FLOW_SIG=H7-H8*DEXP(-H9*EQPLAS)    ! H(?_BAR) = H7-H8*EXP(-H9*?_BAR)]
          DHDE    =H8*H9*DEXP(-H9*EQPLAS)    ! ?H/?(?_BAR) = H8*H9*EXP(-H9*?_BAR)

C     #   MODIFIED-VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.3.) THEN
          FLOW_SIG=H7+H8*EQPLAS +H9*(1.D0-DEXP(-H10*EQPLAS))
          DHDE    =H8+ H9*H10*DEXP(-H10*EQPLAS)

C     #   SWIFT+MODIFIED-VOCE HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.4.) THEN
          FLOW_SIG= H7*(EQPLAS)**H8- H7**H8
     1             +H9*EQPLAS+H10*(1.D0-DEXP(-H11*EQPLAS))
          DHDE= H7*H8*(EQPLAS)**(H8-1.D0)
     1         +H9+H10*H11*DEXP(-H11*EQPLAS)

C     #   HOCKETT-SHERBY HARDENING LAW #
      ELSEIF(FLOW_PAR.EQ.5.) THEN
          IF(EQPLAS .EQ. 0.D0) THEN
              EQPLAS0= 1.D-9
          ELSE
              EQPLAS0= EQPLAS
          END IF
          FLOW_SIG=H7-(H7-H8)*DEXP(-H9*EQPLAS0**H10)
          DHDE=(H7-H8)*H9*H10*EQPLAS0**(H10-1.D0)
     1        *DEXP(-H9*EQPLAS0**H10)
C     #   COMBINED HOLLOMON-VOCE #
      ELSEIF(FLOW_PAR.EQ.6.) THEN
          FLOW_SIG= H7*(H8*(EQPLAS)**H9- H8**H9)
     1             +(1.D0-H7)*(H10*(1.D0-DEXP(-H11*EQPLAS)))
          DHDE= H7*H8*H9*(EQPLAS)**(H9-1.D0)
     1         +(1.D0-H7)*H10*H11*DEXP(-H11*EQPLAS)
      END IF

      RETURN
      END SUBROUTINE YU_HARD
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.4  EQUIVLANET STRESS                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EQV_SIG(HARD_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS EQUIVALENT STRESS FROM CONSTITUTIVE MODEL.
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION SIG(NDIM3), HARD_VAR(NDIM7)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      XSIG=0.D0
      DO I=1, NDIM3
          XSIG= XSIG + SIG(I)*SIG(I)
      END DO
      IF(XSIG .EQ. 0.D0) THEN
          SIG_BAR=0.D0
          RETURN                      ! -> RETURN(END THIS SUBROUTINE) HERE!!
      END IF

      IF(HARD_PAR .LT. 3.) THEN ! HAH: OFF
          CALL YLD(SIG, SIG_BAR)
      ELSEIF(HARD_PAR .EQ. 3.) THEN ! HAH11
          CALL HAH11(HARD_VAR, SIG, SIG_BAR)
      ELSEIF(HARD_PAR .EQ. 4.) THEN ! HAH14
          CALL HAH14(HARD_VAR, SIG, SIG_BAR)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
          CALL HAH20(HARD_VAR, SIG, SIG_BAR)
      END IF
      RETURN
      END SUBROUTINE EQV_SIG
C-----------------------------------------------------------------------
      SUBROUTINE YLD(SIG, SIG_BAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      XSIG=0.D0
      DO I=1, NDIM3
          XSIG= XSIG + SIG(I)*SIG(I)
      END DO
      IF(XSIG .EQ. 0.D0) THEN
          SIG_BAR=0.D0
          RETURN                      ! -> RETURN(END THIS SUBROUTINE) HERE!!
      END IF

      IF(YLD_PAR .EQ. 1.) THEN        ! VON_MISES
          CALL VON_MISES(SIG, SIG_BAR)
      ELSEIF(YLD_PAR .EQ. 2.) THEN    ! HILL48
          CALL HILL_1948(SIG, SIG_BAR)
      ELSEIF(YLD_PAR .EQ. 3.) THEN    ! YLD2000_2D
          IF(NDIM3 .NE. 3) THEN
              PRINT *, '#ERROR: CHECK THE YIELD FUNCTION OPTION!'
              STOP
          ELSE
              CALL YLD2000_2D(SIG, SIG_BAR)
          END IF
      ELSEIF(YLD_PAR .EQ. 4.) THEN    ! YLD2004/_18P
          IF(NDIM3 .NE. 6) THEN
              PRINT *, '#ERROR: CHECK THE YIELD FUNCTION OPTION!'
              STOP
          ELSE
              CALL YLD2004_18P(SIG, SIG_BAR)
          END IF
      END IF
      RETURN
      END SUBROUTINE YLD
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.5   ANISOTROPIC HARDENING MODELS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HAH11(HARD_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EQUIVALENT STRESS OF HAH11 MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM6), DEV_SIG(NDIM6)
      DIMENSION HARD_VAR(NDIM7)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5

C     #   STATEV VARIABLES FOR HAH11 MODEL    #
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO

C     #   EFFECTIVE STRESS OF YIELD_FUNCTION  #
      CALL YLD(SIG, PHI)

C     #   CALCULATE DEVIATORIC STRESS AND DEVIATOR    #
      CALL DEVIATORIC(SIG,DEV_SIG)

C     #   DOUBLE = ^H:S  IN (EQ. 9)   #
      CALL DOUBLE_DOT(DEV_H, DEV_SIG, NDIM6, DOUBLE)
C
C     APPENDIX A: [EQ. A.28][2]
C
C     [3] STATE VARIABLES IN HAH11
C     (1) F-FUNCTION
      IF(DOUBLE .LT. 0.D0) THEN
          F=((G(1)**(-QM))-1.D0)**(1.D0/QM)
      ELSE
          F=((G(2)**(-QM))-1.D0)**(1.D0/QM)
      END IF

C     (2) THE FLUCTUAING COMPONENT: Φ_H
      PHI_H= F*DABS(2.D0*DOUBLE)
C
C     APPENDIX A: [EQ. A.27][2]
C
C     Σ_BAR(Σ) = [Φ^Q + (Φ_H)^Q]^(1/Q) = SIG_BAR
C     PHI = Φ = YLD2000-2D WHICH IS FROM SUB. YLD2000_2D
C     PHI_H: Φ_H = {(F1**Q)*|H:S-|H:S||**Q + (F2**Q)*|H:S+|H:S||**Q}**(1/Q)
      SIG_BAR = (PHI**QM + PHI_H**QM)**(1.D0/QM)

      RETURN
      END SUBROUTINE HAH11
C-----------------------------------------------------------------------
      SUBROUTINE GRAD_HAH11(HARD_VAR, SIG, SIG_BAR, DFDS, DDFDDS)
C     # THIS SUBROUTINE RETURNS ANALYTICAL DERIVATIVE OF HAH11 MODEL.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3)
      DIMENSION G(NDIM5), DEV_H(NDIM6), DEV_SIG(NDIM6), 
     1          PHI_SIG(NDIM3), PHI_SIG_SIG(NDIM3,NDIM3)
      DIMENSION DEV_HK(NDIM3), DEV_HK0(NDIM6), PHI_H_SIG(NDIM3), 
     2          PSI(NDIM3), PSI_SIG(NDIM3,NDIM3)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5

C     [0] STATE VARIABLES FOR HAH11 MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO

C     [1] CALCULATE DOUBLE VALUE
      CALL DEVIATORIC(SIG,DEV_SIG)
      CALL DOUBLE_DOT(DEV_H, DEV_SIG, NDIM6, DOUBLE)
C     [2] VALUES FROM YIELD FUNCTIONS
C     (1) THE EQUIVALENT STRESS OF YIELD FUNCTION
      CALL YLD(SIG, PHI)
C     (2) THE DERIVATIVES OF YIELD FUNCTIONS
      CALL GRAD_YLD(SIG, PHI, PHI_SIG, PHI_SIG_SIG)

C     [3] STATE VARIABLES IN HAH11
C     (1) F-FUNCTION
      IF(DOUBLE .LT. 0.D0) THEN
          F=((G(1)**(-QM))-1.D0)**(1.D0/QM)
          SIGN=-1.D0
      ELSE
          F=((G(2)**(-QM))-1.D0)**(1.D0/QM)
          SIGN=1.D0
      END IF

C     (2) THE FLUCTUATING COMPONENT: Φ_H
      PHI_H= F*DABS(2.D0*DOUBLE)

C     (3) DEV_HK: H:∂Σ'/∂Σ
      NDI=NDIM3-NDIM2
      DO I=1, NDI
          DEV_HK(I)=DEV_H(I)
      END DO

      IF(NDIM3 .EQ. 3) THEN !PLANE STRESS
          DEV_HK(3)=DEV_H(4)
      ELSE IF(NDIM3 .EQ. 6) THEN !3D-STRESS
          DO I=NDI+1,NDIM6
              DEV_HK(I)=DEV_H(I)
          END DO
      END IF
C     (4) ∂Φ_H/∂Σ=2*F1*H'
      PHI_H_SIG=2.D0*F*DEV_HK*SIGN

C     (5) Ψ=Φ^(1-Q)*∂Φ/∂Σ + SIGN(DOUBLE)*Φ_H^(1-Q)*∂Φ_H/∂Σ
      PSI=PHI**(QM-1.D0)*PHI_SIG + PHI_H**(QM-1.D0)*PHI_H_SIG

C     [2] THE FIRST-ORDER DERIVATIVE: DFDS=Σ_BAR^(1-Q)*Ψ
      DFDS= SIG_BAR**(1.D0-QM)*PSI
C-----------------------------------------------------------------------
C     [3] THE SECOND-ORDER DERIVATIVE OF HAH11
      IF(SUA_PAR .GT. 1.) THEN
      DO 900 I=1, NDIM3
      DO 900 J=1, NDIM3
          PSI_SIG(I,J)= (QM-1.D0)*PHI**(QM-2.D0)*PHI_SIG(I)*PHI_SIG(J)
     1                  +PHI**(QM-1.D0)*PHI_SIG_SIG(I,J)
     2                  +(QM-1.D0)*PHI_H**(QM-2.D0)
     3                  *PHI_H_SIG(I)*PHI_H_SIG(J)
900   CONTINUE

      DO 910 I=1, NDIM3
      DO 910 J=1, NDIM3
          DDFDDS(I,J)= (1.D0-QM)*SIG_BAR**(-QM)*DFDS(J)*PSI(I)
     1                +SIG_BAR**(1.D0-QM)*PSI_SIG(I,J)
910   CONTINUE
      ELSE
          DDFDDS=0.D0
      END IF

      RETURN
      END SUBROUTINE GRAD_HAH11
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH11(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
C     #   THIS SUBROUTINE CALCULATE THE INCREMENTS OF HAH11 STATE VARIABLES
C         [EQ. 20 - 24][5]
C         K(1)=K,K(2)=K1,..,K(6)=K5
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM6), G(NDIM5), DEL_G(NDIM5),
     1          DEV_H(NDIM6), DEL_DEV_H(NDIM6), DEV_S(NDIM6)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC

C     #   STATE VARIABLES FOR HAH11 MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL ISO_HARD(HARD_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
C     #   CALCULATE THE INCREMENTS OF INTERNAL PARAMETERS OF HAH MODEL
C     S= DEV_SIG(I+1, N+1)
      CALL DEVIATORIC(SIG,DEV_SIG)
C     S^= DEV_S(I+1, N+1)
      CALL NORM_HAH(DEV_SIG,DEV_S)
C     H^(N):S(I+1, N+1)= DOUBLE
      CALL DOUBLE_DOT(DEV_H, DEV_SIG, NDIM6, DOUBLE)
C     [H^:S^] = COS_X
      CALL DOUBLE_DOT(DEV_S, DEV_H, ,NDIM6, COS_X)
C     COS_X=(8/3)*(S^:H^)
      COS_X=(8.D0/3.D0)*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF
C     #   MICROSTRUCTURE DEVIATOR
      DO I=1, NDIM6
          DEL_DEV_H(I)=XLM*XK*(DEV_S(I)-DEV_H(I)*COS_X)
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
C     #   REVERSE LOADING VARIABLES
C     [1] BAUSCHINGER EFFECT
      DEL_G(1)= ((1.D0-XLM)/2.D0)*XK1*(G(4)/G(1)-1.D0)
     1         +((1.D0+XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(1))
      DEL_G(2)= ((1.D0-XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(2))
     1         +((1.D0+XLM)/2.D0)*XK1*(G(3)/G(2)-1.D0)
C     [2] PERMANENT SOFTENING
      DEL_G(3)= ((1.D0-XLM)/2.D0)*XK5*(XK4-G(3))
      DEL_G(4)= ((1.D0+XLM)/2.D0)*XK5*(XK4-G(4))

      G=G+DEL_G*DEL_EQPLAS
      IF(G(1) .GT. 1.D0) G(1)=1.D0
      IF(G(2) .GT. 1.D0) G(2)=1.D0

C     #   UPDATE STATE VARIABLES FOR HAH11 MODEL
      DO I=1, NDIM5
          HARD_VAR(I)=G(I)
      END DO
      DO I=1, NDIM6
          HARD_VAR(NDIM5+I)=DEV_H(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH11
C-----------------------------------------------------------------------
      SUBROUTINE HAH14(HARD_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EFFECTIVE STRESS GENERATED BY HAH (HAH20H) MODEL.
C     #   BOTH HAH20H & HAH20E ARE IMPLEMENTED IN THIS SUBROUTINE.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM6), DEV_H0(NDIM3),
     1          DEV_SIG(NDIM6)
      DIMENSION SIG_C(NDIM6), SIG_O(NDIM6), SIG_L(NDIM6), SIG_X(NDIM6)
      DIMENSION SIG_L0(NDIM3), SIG_X0(NDIM3)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KHAH2/ XC, XKC, XL, XKL, XI_L
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO

C     [1] PREREQUIRED VARIALBES
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)
      CALL DOUBLE_DOT(DEV_SIG, DEV_H, NDIM6, DOUBLE)
      CALL DOUBLE_DOT(DEV_H0, DEV_H, NDIM6, COS_X)
      HVAL=(8.D0/3.D0)
C      HVAL=1.D0
      COS_X=HVAL*COS_X

C     [2] STRESS DEVIATOR DECOMPOSITION
C     2.1 COLLINEAR STRESS COMPONENTS
      SIG_C=HVAL*DOUBLE*DEV_H
C     2.2 ORTHOGONAL STRESS COMPONENTS
      SIG_O=DEV_SIG-SIG_C
C     2.3 STRESS COMPONENTS FOR LATENT HARDENING
      IF(G(6) .EQ. 0.D0) THEN
          WRITE(*,*) '#ERROR: GL= 0.D0'
          STOP
      END IF
      SIG_L=SIG_C+(1.D0/G(6))*SIG_O
C     2.4 STRESS COMPONENTS FOR CROSS-LOADING CONTRACTION
      SIG_X=(4.D0*(1.D0-G(5)))*SIG_O

C     [3] STABLE COMPONENT: XI(Ξ)
      CALL DEVIATORIC_TRANS(SIG_L, SIG_L0)
      CALL DEVIATORIC_TRANS(SIG_X, SIG_X0)
      CALL YLD(SIG_L0, PHI_L)
      CALL YLD(SIG_X0, PHI_X)
      
      XI=(PHI_L**2.D0 + PHI_X**2.D0)**(1.D0/2.D0)
C-----------------------------------------------------------------------
C     [4] STATE VARIABLES IN HAH11
C     (1) F-FUNCTION
C      VAL=(DSQRT(6.D0*HVAL)/4.D0)
      IF(DOUBLE .LT. 0.D0) THEN
          F=((G(1)**(-QM))-1.D0)**(1.D0/QM)
      ELSE
          F=((G(2)**(-QM))-1.D0)**(1.D0/QM)
      END IF

C     (2) THE FLUCTUAING COMPONENT: Φ_H
      PHI_H= F*DABS(2.D0*DOUBLE)
C-----------------------------------------------------------------------
C     [6] EQUIVALENT STRESS
      SIG_BAR=(XI**QM+PHI_H**QM)**(1.D0/QM)
C      SIG_BAR=SIG_BAR/G(5)

      RETURN
      END SUBROUTINE HAH14
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH14(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION G(NDIM5), DEV_H(NDIM6), DEV_HP(NDIM6), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM6), DEV_H0(NDIM6), DEV_HS(NDIM6),
     1          DEL_DEV_H(NDIM6), DEL_DEV_HP(NDIM6), DEL_G(NDIM5)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KHAH2/ XC, XKC, XL, XKL, XI_L
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL ISO_HARD(HARD_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
C     [1] PREREQUISITE
C     1.1 MICROSCTURE DEVIATOR COEFFICIENTS
      XZ  =5.D0
      XKR =0.1D0
      XKRP=0.02D0
      HVAL=8.D0/3.D0
C      HVAL=1.D0
C     1.2 STRAIN-PATH CHANGE PARAMETER
C     1.2.1
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)

C     1.2.2
      CALL DOUBLE_DOT(DEV_H, DEV_H0, ,NDIM6, COS_X)
      COS_X=HVAL*COS_X
      IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
      IF(DABS(COS_X) .LT. 1.D-7) COS_X=0.D0
      COS_X2=COS_X*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF

C     [2] MICROSTRUCTURE DEVIATOR EVOLUTION
C      VAL=XK*XLM*((DABS(COS_X))**(1.D0/XZ)+G(7))
      VAL=XK*XLM*((DABS(COS_X))**(XKR)+XKRP)
      DO I=1, NDIM6
          DEL_DEV_H(I)=VAL*(DEV_H0(I)-DEV_H(I)*COS_X)
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
C-----------------------------------------------------------------------
C     [3] REVERSE LOADING VARIABLE EVOLUTION
      DEL_G=0.D0
      DEL_G(1)= ((1.D0-XLM)/2.D0)*XK1*(G(4)/G(1)-1.D0)
     1         +((1.D0+XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(1))
      DEL_G(2)= ((1.D0-XLM)/2.D0)*XK2*(XK3*FLOW_SIG0/FLOW_SIG-G(2))
     1         +((1.D0+XLM)/2.D0)*XK1*(G(3)/G(2)-1.D0)
      DEL_G(3)= ((1.D0-XLM)/2.D0)*XK5*(XK4-G(3))
      DEL_G(4)= ((1.D0+XLM)/2.D0)*XK5*(XK4-G(4))
C     [5] CROSS-LOADING CONTRACTION EVOLUTION
      DEL_G(5)=XKC*(1.D0+(XC-1.D0)*COS_X2-G(5))

C     [6] LATENT HARDENING EVOLUTION
      VAL=(FLOW_SIG-FLOW_SIG0)/FLOW_SIG
      VAL2=XL*(1.D0-COS_X2)+COS_X2
      DEL_G(6)= XKL*(VAL*(DSQRT(VAL2)-1.D0)+1.D0-G(6))

C     [8] UPDTAE G-VALUES
      G=G+DEL_G*DEL_EQPLAS
      VAL=XK3*(FLOW_SIG0/FLOW_SIG)
      IF(G(1) .GT. 1.D0) G(1)=1.D0
      IF(G(2) .GT. 1.D0) G(2)=1.D0

C     [9] UPDATE STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          HARD_VAR(I)=G(I)
      END DO
      DO I=1, NDIM6
          HARD_VAR(NDIM5+I)=DEV_H(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH14
C-----------------------------------------------------------------------
      SUBROUTINE HAH20(HARD_VAR, SIG, SIG_BAR)
C     #   THIS SUBROUTINE RETURNS THE EFFECTIVE STRESS GENERATED BY HAH (HAH20H) MODEL.
C     #   BOTH HAH20H & HAH20E ARE IMPLEMENTED IN THIS SUBROUTINE.
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), G(NDIM5), DEV_H(NDIM6), DEV_HP(NDIM6), 
     1          DEV_SIG(NDIM6)
      DIMENSION SIG_C(NDIM6), SIG_O(NDIM6), SIG_L(NDIM6), SIG_X(NDIM6)
      DIMENSION SIG_L0(NDIM3), SIG_X0(NDIM3)
      DIMENSION DEV_H0(NDIM3), DEV_EH(NDIM6), DEV_EH0(NDIM6),
     1          DFDH(NDIM3), DDFDDH(NDIM3, NDIM3),
     2          T(NDIM3, NDIM3), DEV_SIG0(NDIM6)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KHAH2/ XC, XKC, XL, XKL, XI_L
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC
	  
C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
          DEV_HP(I)=HARD_VAR(NDIM5+NDIM6+I)
      END DO

C     [1] PREREQUIRED VARIALBES
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL DOUBLE_DOT(DEV_SIG, DEV_H, NDIM6, DOUBLE)

C     [2] STRESS DEVIATOR DECOMPOSITION
C     2.1 COLLINEAR STRESS COMPONENTS
      SIG_C=DOUBLE*DEV_H
C     2.2 ORTHOGONAL STRESS COMPONENTS
      SIG_O=DEV_SIG-SIG_C
C     2.3 STRESS COMPONENTS FOR LATENT HARDENING
      IF(G(8) .EQ. 0.D0) THEN
          WRITE(*,*) '#ERROR: GL= 0.D0'
          STOP
      END IF
      VAL=1.D0/(XI_L*(G(8)-1.D0)+1.D0)
      SIG_L=VAL*SIG_C+(1.D0/G(8))*SIG_O
C     2.4 STRESS COMPONENTS FOR CROSS-LOADING CONTRACTION
      SIG_X=(4.D0*(1.D0-G(7))/G(8))*SIG_O

C     [3] STABLE COMPONENT: XI(Ξ)
      CALL DEVIATORIC_TRANS(SIG_L, SIG_L0)
      CALL DEVIATORIC_TRANS(SIG_X, SIG_X0)
      CALL YLD(SIG_L0, PHI_L)
      CALL YLD(SIG_X0, PHI_X)
      
      XI=(PHI_L**PM + PHI_X**PM)**(1.D0/PM)
C-----------------------------------------------------------------------
C     [4] HAH20H
      IF(HARD_PAR .EQ. 5.D0) THEN ! HAH20H
C     [4.1] F-FUNCTION
          IF(DOUBLE .LT. 0.D0) THEN
              F=DSQRT(3.D0/8.D0)*(((G(1)**(-QM))-1.D0)**(1.D0/QM))
          ELSE
              F=DSQRT(3.D0/8.D0)*(((G(2)**(-QM))-1.D0)**(1.D0/QM))
          END IF
C     [4.1] FLUCTUATING COMPONENT: Φ_H
          PHI_H= F*DABS(2.D0*DOUBLE)
C-----------------------------------------------------------------------
C     [5] HAH20E
      ELSEIF(HARD_PAR .EQ. 5.5D0) THEN ! HAH20E
C     [5.1] GRADIENT OF YIELD SURFACE AT DEV_H: ε_h
          CALL DEVIATORIC_TRANS(DEV_H, DEV_H0)
          CALL GRAD(0.D0, HARD_VAR, DEV_H0, SIG_BAR_H, DFDH, DDFDDH)
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

C     [5.2] F-FUNCTION
          CALL DOUBLE_DOT(DEV_EH, DEV_SIG, NDIM6, DOUBLE2)
          IF(DOUBLE2 .LT. 0.D0) THEN
              F=DSQRT(3.D0/8.D0)*(((G(1)**(-QM))-1.D0)**(1.D0/QM))
          ELSE
              F=DSQRT(3.D0/8.D0)*(((G(2)**(-QM))-1.D0)**(1.D0/QM))
          END IF

C     [5.3] FLUCTUATING COMPONENT: Φ_H
          CALL DOUBLE_DOT(DEV_EH, DEV_H, NDIM6, DOUBLE3)
          PHI_H= F*DABS(2.D0*DOUBLE2/DOUBLE3)
      END IF
C-----------------------------------------------------------------------
C     [6] EQUIVALENT STRESS
      SIG_BAR=(XI**QM+PHI_H**QM)**(1.D0/QM)
C      SIG_BAR=SIG_BAR/G(5)

      RETURN
      END SUBROUTINE HAH20
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HAH20(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION G(NDIM5), DEV_H(NDIM6), DEV_HP(NDIM6), SIG(NDIM3)
      DIMENSION DEV_SIG(NDIM6), DEV_H0(NDIM6), DEV_HS(NDIM6),
     1          DEL_DEV_H(NDIM6), DEL_DEV_HP(NDIM6), DEL_G(NDIM5)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KHAH1/ QM, PM, XK, XK1, XK2, XK3, XK4, XK5
      COMMON /KHAH2/ XC, XKC, XL, XKL, XI_L
      COMMON /KDEBUG/ KNOEL, KNPT, KKSTEP, KKINC

C     [0] STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          G(I)=HARD_VAR(I)
      END DO
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
          DEV_HP(I)=HARD_VAR(NDIM5+NDIM6+I)
          DEV_HS(I)=HARD_VAR(NDIM5+2*NDIM6+I)
      END DO
C     # DEFINE YIELD STRENGHTH AND FLOW STRESS
      CALL ISO_HARD(HARD_VAR, 0.D0, FLOW_SIG0, DHDE)
      CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)

C     [1] PREREQUIRED VARIABLES
C     1.1 THE OTHER COEFFICIENTS
C     1.1.1   MICROSTRUCTURE DEVIATOR
C      XKP= XK1/2.D0
      XKP= 15.D0
      XI_R= 8.D0
C     1.1.2   BAUSCHINGER RELATED EFFECTS
      XI_B =4.D0
      XI_BP=1.5D0
C     1.1.3   PERMANENT SOFTENING
      XKS =XK1
      XI_S =3.D0
C     1.1.4   CROSS-LOADIN CONTRACTION
      XKCP =XKC
      XI_CP=6.D0
C     1.1.5   LATENT HARDENING
      XI_LP =0.5D0

C     1.2 STRAIN-PATH CHANGE PARAMETER
C     1.2.1
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL NORM_HAH(DEV_SIG, DEV_H0)

C     1.2.2
      CALL DOUBLE_DOT(DEV_H, DEV_H0, NDIM6, COS_X)
      IF(DABS(COS_X) .GT. 1.D0) COS_X=COS_X/DABS(COS_X)
      COS_X2=COS_X*COS_X
      IF(COS_X .LT. 0.D0) THEN
          XLM= -1.D0
      ELSE
          XLM= 1.D0
      ENDIF

C     1.2.3
      CALL DOUBLE_DOT(DEV_HP, DEV_H0, NDIM6, COS_XP)
      IF(DABS(COS_XP) .GT. 1.D0) COS_XP=COS_XP/DABS(COS_XP)
      COS_XP2=COS_XP*COS_XP
      IF(COS_XP .LT. 0.D0) THEN
          XLMP= -1.D0
      ELSE
          XLMP= 1.D0
      ENDIF

C     [2] MICROSTRUCTURE DEVIATOR EVOLUTION
      DO I=1, NDIM6
C     2.1 DEV_HP
          DEL_DEV_HP(I)= XLMP*XKP*(DEV_H0(I)-COS_XP*DEV_HP(I))

C     2.2 DEV_H
C     #   DEV_H EVOLVES ONLY WHEN ITSELF OR DEV_HP ARE CLOSE ENOUGH FROM
C         THE ACTIVE STRESS STATE DEV_SIG AS CONTROLLED BY THE COEFFICIENT
C         XI_R.
          DEL_DEV_H(I)= XLM*XK*((COS_XP2**XI_R)
     1                 +(COS_X2**XI_R))*(DEV_H0(I)-COS_X*DEV_H(I))
      END DO
      DEV_H=DEV_H+DEL_DEV_H*DEL_EQPLAS
      DEV_HP=DEV_HP+DEL_DEV_HP*DEL_EQPLAS
      CALL NORM_HAH(DEV_H, DEV_H)
      CALL NORM_HAH(DEV_HP, DEV_HP)
C-----------------------------------------------------------------------
C     [3] REVERSE LOADING VARIABLE EVOLUTION
      DEL_G=0.D0
C     3.1 G(1)= G-
      DEL_G(1)= ((1.D0+XLM)/2.D0)*(XK1*(1.D0-G(1))/(G(1)**XI_BP)
     1          *(1.D0-(COS_X2**XI_B))
     2          +XK2*(XK3*(FLOW_SIG0/FLOW_SIG)-G(1))*(COS_X2**XI_B))
     3          +((1.D0-XLM)/2.D0)*(XK1*(1.D0-G(1))/(G(1)**XI_BP))

C     3.2 G(2)= G+
      DEL_G(2)= ((1.D0-XLM)/2.D0)*(XK1*(1.D0-G(2))/(G(2)**XI_BP)
     1          *(1.D0-(COS_X2**XI_B))
     2          +XK2*(XK3*(FLOW_SIG0/FLOW_SIG)-G(2))*(COS_X2**XI_B))
     3          +((1.D0+XLM)/2.D0)*(XK1*(1.D0-G(2))/(G(2)**XI_BP))

C     [4] PERMANENT SOFTENING EVOLUTION
      IF(G(4)*XLM .LT. 0.D0) THEN
          G(4) =XLM
          G(9) =G(3)  ! G3*
          G(10)=G(5)  ! GP*
          DEV_HS=DEV_H! H*
          CALL DOUBLE_DOT(DEV_H0, DEV_HS, NDIM6, DOUBLE2) !DOUBLE2= S:H*
          IF(DABS(DOUBLE2).GT.1.D0) DOUBLE2=DOUBLE2/DABS(DOUBLE2)
          DOUBLE22=DOUBLE2*DOUBLE2
          GPP=G(10)-(G(10)-G(9))*DABS(DOUBLE22)**XI_S   !
          IF((GPP-G(5)) .LT. 0.D0) G(5)=GPP
          G(3)= G(5)
      ELSE
          CALL DOUBLE_DOT(DEV_H0, DEV_HS, NDIM6, DOUBLE2) !DOUBLE2= S:H*
          IF(DABS(DOUBLE2).GT.1.D0) DOUBLE2=DOUBLE2/DABS(DOUBLE2)
          DOUBLE22=DOUBLE2*DOUBLE2
          GPP=G(10)-(G(10)-G(9))*DABS(DOUBLE22)**XI_S   !
          IF((GPP-G(5)) .LT. 0.D0) G(5)=GPP
C     4.1 DEL_G3P= ΔG3'
          DEL_G3P= XK5*G(6)*(XK4-G(3))
          G3P= G(3)+DEL_G3P*DEL_EQPLAS
C     4.2 DEL_G(6)= ΔGS
          DEL_G(6)=-XKS*(1.D0-DABS(DOUBLE2))*G(6)
          IF((G3P-G(5)) .LE. 0.D0) THEN
              G(3)= G3P
          ELSE
              G(3)= G(5)
          END IF
          
      END IF

C     [5] CROSS-LOADING CONTRACTION EVOLUTION
      IF(DABS(COS_X) .GE. 0.996) THEN
          DEL_G(7)=XKC*(XC-G(7))
      ELSE
          DEL_G(7)=XKCP*(1.D0-G(7))/(G(7)**XI_CP)
      END IF

C     [6] LATENT HARDENING EVOLUTION
      VAL1= (XL-G(8))*(1.D0-DABS(COS_XP)**XI_LP)
      VAL2= (1.D0-G(8))*(1.D0+DABS(COS_XP)**XI_LP)
      DEL_G(8)=XKL*(VAL1+VAL2)

C     [7] UPDTAE G-VALUES
      G=G+DEL_G*DEL_EQPLAS
      VAL=XK3*(FLOW_SIG0/FLOW_SIG)
      IF(G(1) .GT. 1.D0) THEN
          G(1)=1.D0
      ELSEIF(G(1) .LT. VAL) THEN
          G(1)= VAL
      END IF
      IF(G(2) .GT. 1.D0) THEN
          G(2)=1.D0
      ELSEIF(G(2) .LT. VAL) THEN
          G(2)= VAL
      END IF
      IF(G(3) .GT. 1.D0) THEN
          G(3)=1.D0
      ELSEIF(G(3) .LT. XK4) THEN
          G(3)=XK4
      END IF

C     [9] UPDATE STATE VARIABLES FOR HAH MODEL
      DO I=1, NDIM5
          HARD_VAR(I)=G(I)
      END DO
      DO I=1, NDIM6
          HARD_VAR(NDIM5+I)=DEV_H(I)
          HARD_VAR(NDIM5+NDIM6+I)=DEV_HP(I)
          HARD_VAR(NDIM5+2*NDIM6+I)=DEV_HS(I)
      END DO

      RETURN
      END SUBROUTINE UPDATE_HAH20
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_CHABOCHE(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES STATE VARIABLES OF KINEMATIC HARDENING.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3)
      DIMENSION ALPHA(NDIM3), DEL_ALPHA(NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

      DO I=1, NDIM3
          ALPHA(I)= HARD_VAR(I)
      END DO
      IF(H4 .EQ. 0.D0) THEN   ! CHABOCHE1[9]
          HP1= H1
          HP2= H2
          CALL CB_HARD(EQPLAS, FLOW_SIG, DHDE)
      ELSE                    ! CHABOCHE2[10]
          CALL CB_HARD(EQPLAS, FLOW_SIG, DHDE)
          CALL ISO_HARD(HARD_VAR, EQPLAS, FLOW_ISO, DSDE)
          AP=FLOW_ISO-FLOW_SIG
          DAPDE=DSDE-DHDE
          HP1= H1*EQPLAS+ H2
          HP2= (HP1-DAPDE)/AP
      END IF
      DEL_ALPHA= HP1*SIG/FLOW_SIG-HP2*ALPHA
      DO I=1, NDIM3
          HARD_VAR(I)=ALPHA(I)+DEL_ALPHA(I)*DEL_EQPLAS
          HARD_VAR(NDIM3+I)=DEL_ALPHA(I)
      END DO
      SIG= SIG-DEL_ALPHA*DEL_EQPLAS
      RETURN
      END SUBROUTINE UPDATE_CHABOCHE
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_YU(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
C     #   THIS SUBROUTINE UPDATES BACK-STRESS VALUES OF YOSHIDA-UEMORI.
      IMPLICIT REAL*8(A-H, O-Z)
      
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3)
      DIMENSION ALPHA(NDIM3), DEL_ALPHA(NDIM3), ALPHAS(NDIM3), 
     1          DEL_ALPHAS(NDIM3), BETA(NDIM3), DEL_BETA(NDIM3), 
     2          Q(NDIM3), DEL_Q(NDIM3)
      DIMENSION DFDS(NDIM3), DDFDDS(NDIM3,NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KKIN_HARD/ H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11

C     #   ISOTORPIC-KINEMATIC HARDENING
C     [1] STATE VARIALBES
      DO I=1, NDIM3
          ALPHA(I)= HARD_VAR(I)
          BETA(I) = HARD_VAR(2*NDIM3+I)
          Q(I)=HARD_VAR(3*NDIM3+I)
      END DO
      R= HARD_VAR(4*NDIM3+1) !BOUNDING SURFACE

C     [2] STATE VARIALBE INCREMENTS
C     #   PREREQUISITES
      ALPHAS=ALPHA-BETA
      A=DABS(H1+R-H2) ! [EQ. 7b][10]: A=B+R-Y
      CALL YLD(ALPHAS, ALPHAS_BAR) ! [EQ. 7b][10]: ALPHAS_BAR
      CALL YLD(SIG+ALPHA, SIG_BAR)
      CALL GRAD(0.D0, HARD_VAR, SIG+ALPHA, SIG_BAR, DFDS, DDFDDS)
      IF(ALPHAS_BAR .EQ. 0.D0) THEN
          DEL_ALPHAS=H3*(A/H2*SIG)
      ELSE
          DEL_ALPHAS= H3*((A/H2)*SIG-DSQRT(A/ALPHAS_BAR)
     1                      *ALPHAS)
      END IF
      DEL_BETA=H4*((2.D0/3.D0)*H5*DFDS-BETA)
      DEL_ALPHA= DEL_ALPHAS+DEL_BETA

C     #   ISOTROPIC HARDENING FOR BOUNDING SURFACE
      CALL YU_HARD(EQPLAS, R, DEL_R)
      DEL_R=0.D0

C     #   WORK-HARDENING STAGNATION
C     [1] STAGNATION SURFACE: VON-MISES YIELD FUNCTION
      CALL YLD(BETA-Q, STAG_BAR)
C     [2] STATE VARIALBE INCREMENTS
      CALL DOUBLE_DOT(BETA-Q,DEL_BETA*DEL_EQPLAS,NDIM3,STAG_DOUBLE)! [EQ. 23a][10]
      IF(STAG_DOUBLE.GT.0.D0) THEN
          ! #ASSUMPTION: STAG_R= STAG_BAR
          STAG_GAMMA=(3.D0/2.D0)*STAG_DOUBLE/STAG_BAR
          STAG_MU=(1.D0-H6)*STAG_GAMMA/STAG_BAR         ! [EQ. 27][10]
      ELSE
C         DEL_R=0.D0
          STAG_GAMMA=0.D0
          STAG_MU=0.D0
      END IF
      DEL_Q=STAG_MU*(BETA-Q)
c     [3] UPDATE STATE VARIALBES
      ALPHA=ALPHA+DEL_ALPHA*DEL_EQPLAS
      BETA=BETA+DEL_BETA*DEL_EQPLAS
      DO I=1, NDIM3
          HARD_VAR(I)=    ALPHA(I)
          HARD_VAR(NDIM3+I)=  DEL_ALPHA(I)
          HARD_VAR(2*NDIM3+I)= BETA(I)
          HARD_VAR(3*NDIM3+I)= Q(I)+DEL_Q(I)*DEL_EQPLAS
      END DO
      HARD_VAR(4*NDIM3+1)= R
      HARD_VAR(4*NDIM3+2)= ALPHAS_BAR

      SIG= SIG-DEL_ALPHA*DEL_EQPLAS
      RETURN
      END SUBROUTINE UPDATE_YU
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_RGBV(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
C     # THIS SUBROUTINE FOLLOWS THE DISLOCATION EVOLUTION LAWS BY E. RAUCH'S WORK [15].
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3), DEV_SIG(NDIM6), 
     1          DEV_H0(NDIM6), H(3), PSIG(3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KISO_HARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10

C     # DETERMINE STRAIN-PATH CHANGE PARAMETER
      IF(HARD_PAR.GE.3.) THEN ! HAH MODEL
          CALL COSX(HARD_VAR,SIG,DEV_H0,COS_X)
      ELSE
          CALL PRIN_YLD2004_18P(SIG, H, PSIG)
          VAL=SUM(PSIG)
          IF(VAL .GE. 0.D0) THEN
              COS_X=1.D0
          ELSE
              COS_X=-1.D0
          END IF
      END IF
C     # READ STATE VARIABLE
      RHO_F1 =HARD_VAR(NDIM7-NDIM8+1)
      RHO_R1 =HARD_VAR(NDIM7-NDIM8+2)
      RHO_F01=HARD_VAR(NDIM7-NDIM8+3)
      RHO_F2 =HARD_VAR(NDIM7-NDIM8+4)
      RHO_R2 =HARD_VAR(NDIM7-NDIM8+5)
      RHO_F02=HARD_VAR(NDIM7-NDIM8+6)
      RHO_MAX=1.D16
      RHO_MIN=1.D6
      XM=1.D0
      XM2=0.5D0
C     # UPDATE DISLOCATION DENSITY
      IF(COS_X .GE. 0.D0) THEN
          XLAM=DSQRT(RHO_F1)/(P5*P6)+1.D0/(P5*P7)
          XLAM=P1*XLAM
          DEL_RHOF=XLAM-P8*RHO_F1
          DEL_RHOR=-XLAM*(RHO_R1/(RHO_F01))**XM
          RHO_F1=RHO_F1+DEL_RHOF*DABS(DEL_EQPLAS)
          RHO_R1=RHO_R1+DEL_RHOR*DABS(DEL_EQPLAS)
          RHO_F2=(1.D0-P9)*RHO_F1!+RHO_R1
          RHO_R2=P9*RHO_F1
          RHO_F02=RHO_F1
          RHO_TOT=RHO_F1+RHO_R1
          DRDE=DEL_RHOF+DEL_RHOR
          XX=(RHO_R1/RHO_TOT)**1.D0
      ELSE
          XLAM=DSQRT(RHO_F2)/(P5*P6)+1.D0/(P5*P7)
          XLAM=P1*XLAM
          DEL_RHOF=XLAM-P8*RHO_F2
          DEL_RHOR=-XLAM*(RHO_R2/(RHO_F02))**XM
          RHO_F2=RHO_F2+DEL_RHOF*DABS(DEL_EQPLAS)
          RHO_R2=RHO_R2+DEL_RHOR*DABS(DEL_EQPLAS)
          RHO_F1=(1.D0-P9)*RHO_F2!+RHO_R2
          RHO_R1=P9*RHO_F2
          RHO_F01=RHO_F2
          RHO_TOT=RHO_F2+RHO_R2
          DRDE=DEL_RHOF+DEL_RHOR
          XX=(RHO_R2/RHO_TOT)**1.D0
      END IF
C     # UPDATE BACK-STRESS
C      XX=HARD_VAR(NDIM7-NDIM8+9)
C      TAU=P2+P3*P4*P5*DSQRT(RHO_TOT)
C      XXS=P10*(TAU-P2)
C      DEL_XX=P11*(XXS-XX)
C      XX=XX+DEL_XX*DABS(DEL_EQPLAS)
C     # UPDATE STATE VARIABLE
      HARD_VAR(NDIM7-NDIM8+1)=RHO_F1
      HARD_VAR(NDIM7-NDIM8+2)=RHO_R1
      HARD_VAR(NDIM7-NDIM8+3)=RHO_F01
      HARD_VAR(NDIM7-NDIM8+4)=RHO_F2
      HARD_VAR(NDIM7-NDIM8+5)=RHO_R2
      HARD_VAR(NDIM7-NDIM8+6)=RHO_F02
      HARD_VAR(NDIM7-NDIM8+7)=DABS(RHO_TOT)
      HARD_VAR(NDIM7-NDIM8+8)=DRDE
      HARD_VAR(NDIM7-NDIM8+9)=XX

      RETURN
      END SUBROUTINE UPDATE_RGBV
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_RGBV2(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
C     # THIS SUBROUTINE FOLLOWS THE DISLOCATION EVOLUTION LAWS BY K. KITAYAMA'S WORK [16].
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3), DEV_SIG(NDIM6), 
     1          DEV_H0(NDIM6), H(3), PSIG(3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KISO_HARD/ P1, P2, P3, P4, P5, P6, P7, P8, P9, P10

C     # DETERMINE STRAIN-PATH CHANGE PARAMETER
      IF(HARD_PAR.GE.3.) THEN ! HAH MODEL
          CALL COSX(HARD_VAR,SIG,DEV_H0,COS_X)
      ELSE
          CALL PRIN_YLD2004_18P(SIG, H, PSIG)
          VAL=SUM(PSIG)
          IF(VAL .GE. 0.D0) THEN
              COS_X=1.D0
          ELSE
              COS_X=-1.D0
          END IF
      END IF
C     # READ STATE VARIABLE
      RHO_F=HARD_VAR(NDIM7-NDIM8+1)
      RHO_R1=HARD_VAR(NDIM7-NDIM8+2)
      RHO_R2=HARD_VAR(NDIM7-NDIM8+3)
      RHO_L=HARD_VAR(NDIM7-NDIM8+4)
      RHO_F01=HARD_VAR(NDIM7-NDIM8+5)
      RHO_F02=HARD_VAR(NDIM7-NDIM8+6)
      RHO_TOT=HARD_VAR(NDIM7-NDIM8+7)
      RHO_RTOT=HARD_VAR(NDIM7-NDIM8+9)
C     # VARIOUS REVERSIBLE PARAMETER
      VAL1=(RHO_MAX+RHO_MIN)/(RHO_MAX-RHO_MIN)
      VAL2=(RHO_MAX+RHO_MIN)/2.D0
      VAL3=(DSQRT(RHO_TOT)-VAL2)/VAL2
      RP=0.5D0-0.5D0*DTANH(3.D0*VAL1*VAL3)
C     # UPDATE DISLOCATION DENSITY
      XLAM=DSQRT(RHO_TOT)/(P4*P5)+1.D0/(P4*P6) ! MEAN FREE-PATH
      XLAM=P1*XLAM
      DEL_RHOF=(1.D0-P8)*XLAM-P7*RHO_F1
      XM=1.D0
      IF(COS_X .GE. 0.D0) THEN
          DEL_RHOR1=P8*XLAM-P7*RHO_R1
          DEL_RHOR2=-XLAM*(DABS(RHO_R2/RHO_F02))**XM
          RHO_F01=RHO_TOT
          XX=(RHO_R1/RHO_TOT)**0.5D0
      ELSE
          DEL_RHOR1=-XLAM*(DABS(RHO_R1/RHO_F01))**XM
          DEL_RHOR2=P8*XLAM-P7*RHO_R2
          RHO_F02=RHO_TOT
          XX=(RHO_R2/RHO_TOT)**0.5D0
      END IF
      RHO_F=RHO_F+DEL_RHOF*DABS(DEL_EQPLAS)
      RHO_R1=RHO_R1+DEL_RHOR1*DABS(DEL_EQPLAS)
      RHO_R2=RHO_R2+DEL_RHOR2*DABS(DEL_EQPLAS)
      RHO_TOT=RHO_F+RHO_R1+RHO_R2
      RHO_RTOT=RHO_R1+RHO_R2
      DRDE=DEL_RHOF+DEL_RHOR1+DEL_RHOR2
C     # UPDATE BACK-STRESS
C      XX=HARD_VAR(NDIM7-NDIM8+9)
C      TAU=P2+P3*P4*P5*DSQRT(RHO_TOT)
C      XXS=P10*(TAU-P2)
C      DEL_XX=P11*(XXS-XX)
C      XX=XX+DEL_XX*DABS(DEL_EQPLAS)
C     # UPDATE STATE VARIABLE
      HARD_VAR(NDIM7-NDIM8+1)=RHO_F
      HARD_VAR(NDIM7-NDIM8+2)=RHO_R1
      HARD_VAR(NDIM7-NDIM8+3)=RHO_R2
      HARD_VAR(NDIM7-NDIM8+4)=RHO_L
      HARD_VAR(NDIM7-NDIM8+5)=RHO_F01
      HARD_VAR(NDIM7-NDIM8+6)=RHO_F02
      HARD_VAR(NDIM7-NDIM8+7)=RHO_TOT
      HARD_VAR(NDIM7-NDIM8+8)=DRDE
      HARD_VAR(NDIM7-NDIM8+9)=RHO_RTOT
      RETURN
      END SUBROUTINE UPDATE_RGBV2
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_HARD_VAR(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     # UPDATE ANISOTROPIC HARDENING MODEL
      IF(HARD_PAR .EQ. 1.) THEN ! CHABOCHE
          CALL UPDATE_CHABOCHE(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 2.) THEN ! YOSHIDA-UEMORI
          CALL UPDATE_YU(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      ELSEIF(HARD_PAR .EQ. 3.) THEN ! HAH11
          CALL UPDATE_HAH11(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      ELSEIF(HARD_PAR .EQ. 4.) THEN ! HAH14
          CALL UPDATE_HAH14(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN ! HAH20
          CALL UPDATE_HAH20(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      END IF
      
C     # UPDATE DISLOCATION HARDENING MODEL
      IF(FLOW_PAR.EQ.7.) THEN
          CALL UPDATE_RGBV2(HARD_VAR, SIG, EQPLAS, DEL_EQPLAS)
      END IF

      RETURN
      END SUBROUTINE UPDATE_HARD_VAR
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.6 ANISOTROPIC YIELD FUNCTIONS                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VON_MISES(SIG, PHI)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM6)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      IF(FLOOR(HARD_PAR) .GE. 4.) THEN
          IF(NDIM3 .EQ. 3.) THEN
              DEV_SIG(1)= SIG(1)
              DEV_SIG(2)= SIG(2)
              DEV_SIG(3)= -(DEV_SIG(1)+DEV_SIG(2))
              DEV_SIG(4)= SIG(3)
          ELSE
              DEV_SIG= SIG
          END IF
      ELSE
          CALL DEVIATORIC(SIG, DEV_SIG)
      END IF

      CALL DOUBLE_DOT(DEV_SIG, DEV_SIG, NDIM6, PSI)
      PHI= DSQRT((3.D0/2.D0)*PSI)

      RETURN
      END SUBROUTINE VON_MISES
C-----------------------------------------------------------------------
      SUBROUTINE GRAD_VM(SIG, PHI, PHI_SIG, PHI_SIG_SIG)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM6), DEV_SIG0(NDIM3),
     1          PHI_SIG(NDIM3), PHI_SIG_SIG(NDIM3, NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      
      PHI_SIG=0.D0
      PHI_SIG_SIG=0.D0
      
      CALL DEVIATORIC(SIG, DEV_SIG)
      CALL DEVIATORIC_TRANS(DEV_SIG, DEV_SIG0)
      PHI_SIG=(3.D0/2.D0)*DEV_SIG0/PHI

      RETURN
      END SUBROUTINE GRAD_VM
C-----------------------------------------------------------------------
      SUBROUTINE HILL_1948(SIG, PHI)
C     #   THIS SUBROUTINE CALCULATES THE EQUIVALENT STRESS THROUGH HILL 1948
C         YIELD FUNCTION UNDER PLANE STRESS ASSUMPTION.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), STRESS(6)

      COMMON /KHILL1948/ HF, HG, HH, HL, HM, HN
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
C     S(6)=[S11 S22 S33 S12 S13 S23]
      NDI=NDIM3-NDIM2
      STRESS=0.D0
      DO I=1, NDI
          STRESS(I)= SIG(I)
      END DO

      IF(NDIM3 .EQ. 3) THEN
          STRESS(4)= SIG(3)
      ELSE
          STRESS(4:6)=SIG(4:6)
      END IF

      PSI= HF*(STRESS(2)-STRESS(3))**2+HG*(STRESS(3)-STRESS(1))**2
     1   +HH*(STRESS(1)-STRESS(2))**2
     2   +2.D0*HL*STRESS(6)**2+2.D0*HM*STRESS(5)**2+2.D0*HN*STRESS(4)**2

      PHI=DSQRT(PSI)

      RETURN
      END SUBROUTINE HILL_1948
C-----------------------------------------------------------------------
      SUBROUTINE YLD2000_2D(SIG, PHI)
C     #   THIS SUBROUTINE RETURNS THE EQUIVALENT STRESS OF YLD2000_2D MODEL
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(3), PHI_SIG(3), X1P(3), X2P(3), PX1P(2), PX2P(2),
     &          X1P_SIG(3,3), X2P_SIG(3,3), PX1P_X1P(2,3),PX2P_X2P(2,3),
     &          PHI1P_PX1P(2), PHI2P_PX2P(2), PSI_SIG(3), ST(3),
     &          XL1P(3,3), XL2P(3,3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KYLD2K/ AM, A1, A2, A3, A4, A5, A6, A7, A8

C     [1] LINEAR TRANSFORMATION
      CALL LTM_YLD2000_2D(SIG, XL1P, XL2P, X1P, X2P)
C     [2] PRINCIPAL STRESS
      CALL PRIN_YLD2000_2D(X1P, XR1P, PX1P)
      CALL PRIN_YLD2000_2D(X2P, XR2P, PX2P)

C     [3] ANISOTROPIC YILED FUNCTION [EQ. 10 & EQ. 15][1]
      ST(1)=DABS(PX1P(1)-PX1P(2))            !-> |X'_1 - X'_2|
      ST(2)=DABS(2.D0*PX2P(2)+PX2P(1))        !-> |2*X''_2 + X''_1|
      ST(3)=DABS(2.D0*PX2P(1)+PX2P(2))        !-> |2*X''_1 + X''_2|
C
      PHI1P=ST(1)**AM                ! -> Φ'
      PHI2P=ST(2)**AM + ST(3)**AM    ! -> Φ''

C     Ψ = Φ' + Φ''
C     Φ = {0.5*[Φ' + Φ'']}**(1/A) IN [1]
C       = {0.5*Ψ}**(1/AM)
      PSI= PHI1P + PHI2P
      PHI=((1.D0/2.D0)*(PSI))**(1.D0/AM) ! -> Σ_BAR

      RETURN
      END SUBROUTINE YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE LTM_YLD2000_2D(SIG, XL1P, XL2P, X1P, X2P)
C     #   THIS SUBROUTINE RETURNS THE LINEAR TRANSFORMATION MATIRX REQUIRED
C         IN YLD2000_2D MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(3), XL1P(3,3), XL2P(3,3), X1P(3), X2P(3),
     1          C1P(3,3), C2P(3,3), T(3,3)

      COMMON /KYLD2K/ AM, A1, A2, A3, A4, A5, A6, A7, A8
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      C1P= 0.D0
      C1P(1,1)= A1
      C1P(2,2)= A2
      C1P(3,3)= A7

      C2P= 0.D0
      C2P(1,1)= (4.D0/3.D0)*A5 - (1.D0/3.D0)*A3
      C2P(1,2)= (2.D0/3.D0)*(A6-A4)
      C2P(2,1)= (2.D0/3.D0)*(A3-A5)
      C2P(2,2)= (4.D0/3.D0)*A4 - (1.D0/3.D0)*A6
      C2P(3,3)= A8

      CALL DEVIATORIC_TENS(T)

C     LINEAR TRANSFORMATION [EQ. 14][1]
      IF(FLOOR(HARD_PAR) .GE. 4.) THEN
          X1P=MATMUL(C1P, SIG)
          X2P=MATMUL(C2P, SIG)
      ELSE
          XL1P=MATMUL(C1P,T)
          XL2P=MATMUL(C2P,T)
          X1P=MATMUL(XL1P, SIG)
          X2P=MATMUL(XL2P, SIG)
      END IF

      RETURN
      END SUBROUTINE LTM_YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE PRIN_YLD2000_2D(XP, XRP, PXP)
C     # THIS SUBROUTINE RETURNS PRINCIPAL VALUES OF STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION XP(3), PXP(2)

C     THE PRINCIPAL VALUES OF X' & X'' [EQ. 17][1]
      XRP=DSQRT(((XP(1)-XP(2))/2.D0)**2.D0 + XP(3)**2.D0)    ! -> SQRT(Δ'/4)
C
      PXP(1)= (XP(1)+XP(2))/2.D0+XRP     !-> X'_1
      PXP(2)= (XP(1)+XP(2))/2.D0-XRP     !-> X'_2


      RETURN
      END SUBROUTINE PRIN_YLD2000_2D
C-----------------------------------------------------------------------
      SUBROUTINE GRAD_YLD2000(SIG, PHI, PHI_SIG, PHI_SIG_SIG)
C     #   THIS SUBROUTINE RETURNS ANALYTICAL DERIVATIVES OF YLD2000_2D MODEL
C         REF: [6] J. W. YOON ET AL. IJP (2004)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(3), XL1P(3,3), XL2P(3,3), X1P(3), X2P(3), PX1P(2),
     1          PX2P(2), PHI1P_PX1P(2), PHI2P_PX2P(2), PX1P_X1P(2,3),
     2          PX2P_X2P(2,3), X1P_SIG(3,3), X2P_SIG(3,3), PSI_SIG(3),
     3          PHI_SIG(3)
      
      DIMENSION PHI_SIG_SIG(3,3), PHI_SIG_PSI(3), PSI_SIG_SIG(3,3),
     1          PHI1P_PX1P_PX1P(2,2), PHI2P_PX2P_PX2P(2,2),
     2          PX1P_X1P_X1P(2,3,3), PX2P_X2P_X2P(2,3,3),
     3          PX1P_SIG(2,3), PX2P_SIG(2,3),
     4          PHI1P_PX1P_SIG(2,3), PHI2P_PX2P_SIG(2,3),
     5          PHI1P_SIG_SIG(3,3), PHI2P_SIG_SIG(3,3),
     6          PHI1P_X1P_X1P(3,3), PHI2P_X2P_X2P(3,3),
     7          PHI1P_X1P_SIG(3,3), PHI2P_X2P_SIG(3,3)
      
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KYLD2K/ AM, A1, A2, A3, A4, A5, A6, A7, A8

      CALL LTM_YLD2000_2D(SIG, XL1P, XL2P, X1P, X2P)
      CALL PRIN_YLD2000_2D(X1P, XR1P, PX1P)
      CALL PRIN_YLD2000_2D(X2P, XR2P, PX2P)

      PHI_SIG=0.D0
      PHI_SIG_SIG=0.D0

C     #   [1] THE FIRST DERIVATIVE OF YLD2000_2D MODEL

C     (1) ∂Φ/∂Ψ
      PHI_PSI=(1.D0/(2.D0*AM))*PHI**(1.D0-AM)

C     (2) ∂Ψ/∂*X'=∂Φ'/∂*X'
C     APPENDIX A: [EQ. A1.4][1]
      PHI1P_PX1P(1)=AM*(PX1P(1)-PX1P(2))**(AM-1.D0)
      PHI1P_PX1P(2)=-PHI1P_PX1P(1)

C     (3) ∂Ψ/∂*X''=∂Φ''/∂*X''
C     APPENDIX A: [EQ. A1.8][1]
      PHI2P_PX2P(1)=AM*(2.D0*PX2P(2)+PX2P(1))**(AM-1.D0)
     1            +2.D0*AM*(2.D0*PX2P(1)+PX2P(2))**(AM-1.D0)
      
      PHI2P_PX2P(2)=2.D0*AM*(2.D0*PX2P(2)+PX2P(1))**(AM-1.D0)
     1            +AM*(2.D0*PX2P(1)+PX2P(2))**(AM-1.D0)

C     (4) ∂*X'/∂X'
C     APPENDIX A: [EQ. A1.5][1]
C     SINGULAR CASE: Δ' = 0 / X'_XX = X'_YY / X'_XY = 0
      IF(XR1P .EQ. 0.D0) THEN
          PX1P_X1P(1,1)= 1.D0/2.D0        ! -> D(X'_1)/D(X'_XX)
          PX1P_X1P(1,2)= 1.D0/2.D0        ! -> D(X'_2)/D(X'_XX)
          PX1P_X1P(1,3)= 0.D0             ! -> D(X'_1)/D(X'_XY)
          PX1P_X1P(2,1)= 1.D0/2.D0        ! -> D(X'_1)/D(X'_YY)
          PX1P_X1P(2,2)= 1.D0/2.D0        ! -> D(X'_2)/D(X'_XX)
          PX1P_X1P(2,3)= 0.D0             ! -> D(X'_2)/D(X'_XY)

C     GENERAL CASE = Δ' ≠ 0 / X'_XX ≠ X'_YY / X'_XY ≠ 0
      ELSE
          PX1P_X1P(1,1)= 1.D0/2.D0+(X1P(1)-X1P(2))/(4.D0*XR1P)
          PX1P_X1P(1,2)= 1.D0/2.D0-(X1P(1)-X1P(2))/(4.D0*XR1P)
          PX1P_X1P(1,3)= X1P(3)/XR1P
          PX1P_X1P(2,1)= PX1P_X1P(1,2)
          PX1P_X1P(2,2)= PX1P_X1P(1,1)
          PX1P_X1P(2,3)=-PX1P_X1P(1,3)
      END IF

C     (5) ∂*X''/∂X''
C     APPENDIX A: [EQ. A1.9][1]
C     SINGULAR CASE: Δ'' = 0 / X''_XX = X''_YY / X''_XY = 0
      IF(XR2P .EQ. 0.0) THEN
          PX2P_X2P(1,1)= 1.D0/2.D0        ! -> D(X''_1)/D(X''_XX)
          PX2P_X2P(1,2)= 1.D0/2.D0        ! -> D(X''_2)/D(X''_XX)
          PX2P_X2P(1,3)= 0.D0             ! -> D(X''_1)/D(X''_XY)
          PX2P_X2P(2,1)= 1.D0/2.D0        ! -> D(X''_1)/D(X''_YY)
          PX2P_X2P(2,2)= 1.D0/2.D0        ! -> D(X''_2)/D(X''_XX)
          PX2P_X2P(2,3)= 0.D0             ! -> D(X''_2)/D(X''_XY)

C     GENERAL CASE = Δ'' ≠ 0 / X''_XX ≠ X''_YY / X''_XY ≠ 0
      ELSE
          PX2P_X2P(1,1)= 1.D0/2.D0+(X2P(1)-X2P(2))/(4.D0*XR2P)
          PX2P_X2P(1,2)= 1.D0/2.D0-(X2P(1)-X2P(2))/(4.D0*XR2P)
          PX2P_X2P(1,3)= X2P(3)/XR2P
          PX2P_X2P(2,1)= PX2P_X2P(1,2)
          PX2P_X2P(2,2)= PX2P_X2P(1,1)
          PX2P_X2P(2,3)=-PX2P_X2P(1,3)
      ENDIF


C     (6) ∂X'/∂Σ AND ∂X''/∂Σ
      DO I=1,3
          DO J=1,3
              X1P_SIG(I,J)=XL1P(I,J)
              X2P_SIG(I,J)=XL2P(I,J)
          END DO
      END DO

C     (7) ∂Ψ/∂Σ = [∂Φ'/∂*X']*[∂*X'/∂X']*[∂X'/∂Σ] + [∂Φ''/∂*X'']*[∂*X''/∂X'']*[∂X''/∂Σ]
C     -> *X: PRINCIPAL VALUES OF X
      PSI_SIG=0.D0
      DO I=1,3
          DO J=1,2
              DO K=1,3
                  PSI_SIG(I)= PSI_SIG(I)
     &               +(PHI1P_PX1P(J)*PX1P_X1P(J,K)*X1P_SIG(K,I))
     &               +(PHI2P_PX2P(J)*PX2P_X2P(J,K)*X2P_SIG(K,I))
              END DO
          END DO
C     (8) ∂Φ/∂Σ = ∂Φ/∂Ψ * ∂Ψ/∂Σ = DFDS
          PHI_SIG(I)=PHI_PSI*PSI_SIG(I)
      END DO
C-----------------------------------------------------------------------
C     #   THE SECOND DERIVATIVE OF YLD2000_2D MODEL
C     #   ∂∂Φ/∂Σ∂Σ= [∂∂Φ/∂Σ∂Ψ]*[∂Ψ/∂Σ]+[∂Φ/∂Ψ]*[∂∂Ψ/∂Σ∂Σ]
      IF(SUA_PAR .NE. 1.) THEN
C     #   [1] PART-1
C         (1) ∂∂Φ/∂Σ∂Ψ
          PHI_SIG_PSI=((1.D0-AM)/PHI)*PHI_PSI*PHI_SIG
C     #   [2] PART-2
C         [A]   ∂∂Φ'/∂Σ∂*X'=[∂∂Φ'/∂*X'∂*X']*[∂*X'/∂X']*[∂X'/∂Σ]
C         (A.1) ∂∂Ψ/∂*X'∂*X' = ∂∂Φ'/∂*X'∂*X'
          PHI1P_PX1P_PX1P(1,1)= AM*(AM-1.D0)*(PX1P(1)-PX1P(2))
     1                          **(AM-2.D0)
          PHI1P_PX1P_PX1P(2,1)=-PHI1P_PX1P_PX1P(1,1)
          PHI1P_PX1P_PX1P(1,2)= PHI1P_PX1P_PX1P(2,1)
          PHI1P_PX1P_PX1P(2,2)= PHI1P_PX1P_PX1P(1,1)

C         (A.2) ∂∂Ψ/∂*X''∂*X'' = ∂∂Φ''/∂*X''∂*X''
          X2_X1=AM*(AM-1.D0)*(2.D0*PX2P(2)+PX2P(1))**(AM-2.D0)
          X1_X2=AM*(AM-1.D0)*(2.D0*PX2P(1)+PX2P(2))**(AM-2.D0)

          PHI2P_PX2P_PX2P(1,1)=1.D0*X2_X1+4.D0*X1_X2
          PHI2P_PX2P_PX2P(2,1)=2.D0*X2_X1+2.D0*X1_X2
          PHI2P_PX2P_PX2P(1,2)=PHI2P_PX2P_PX2P(2,1)
          PHI2P_PX2P_PX2P(2,2)=4.D0*X2_X1+1.D0*X1_X2

C         [B]   ∂∂*X'/∂Σ∂X'=[∂∂*X'/∂X'∂X']*[∂X'/∂Σ]
C         (B.1) ∂∂*X1'/∂X'∂X'
          PX1P_X1P_X1P(1,1,1)= (1.D0/(4.D0*XR1P))
     1                   -(X1P(1)-X1P(2))**(2.D0)/(16.D0*XR1P**(3.D0))
          PX1P_X1P_X1P(1,2,1)=-PX1P_X1P_X1P(1,1,1)
          PX1P_X1P_X1P(1,3,1)=-(X1P(1)-X1P(2))*X1P(3)
     1                         /(4.D0*XR1P**(3.D0))
          PX1P_X1P_X1P(1,1,2)= PX1P_X1P_X1P(1,2,1)
          PX1P_X1P_X1P(1,2,2)= PX1P_X1P_X1P(1,1,1)
          PX1P_X1P_X1P(1,3,2)=-PX1P_X1P_X1P(1,3,1)

          PX1P_X1P_X1P(1,1,3)= PX1P_X1P_X1P(1,3,1)
          PX1P_X1P_X1P(1,2,3)= PX1P_X1P_X1P(1,3,2)
          PX1P_X1P_X1P(1,3,3)= (1.D0/XR1P)-(X1P(3)**2.D0)
     1                        /(XR1P**(3.D0))

C         (B.2) ∂∂*X2'/∂X'∂X'
          PX1P_X1P_X1P(2,:,:)=-PX1P_X1P_X1P(1,:,:)

C         (B.3) ∂∂*X1''/∂X''∂X''
          PX2P_X2P_X2P(1,1,1)= (1.D0/(4.D0*XR2P))
     1                   -(X2P(1)-X2P(2))**(2.D0)/(16.D0*XR2P**(3.D0))
          PX2P_X2P_X2P(1,2,1)=-PX2P_X2P_X2P(1,1,1)
          PX2P_X2P_X2P(1,3,1)=-(X2P(1)-X2P(2))*X2P(3)
     1                         /(4.D0*XR2P**(3.D0))
          PX2P_X2P_X2P(1,1,2)= PX2P_X2P_X2P(1,2,1)
          PX2P_X2P_X2P(1,2,2)= PX2P_X2P_X2P(1,1,1)
          PX2P_X2P_X2P(1,3,2)=-PX2P_X2P_X2P(1,3,1)

          PX2P_X2P_X2P(1,1,3)= PX2P_X2P_X2P(1,3,1)
          PX2P_X2P_X2P(1,2,3)= PX2P_X2P_X2P(1,3,2)
          PX2P_X2P_X2P(1,3,3)= (1.D0/XR2P)-(X2P(3)**2.D0)
     1                        /(XR2P**(3.D0))

C         (B.4) ∂∂*X2''/∂X''∂X''
          PX2P_X2P_X2P(2,:,:)=-PX2P_X2P_X2P(1,:,:)

C         ##  ∂∂Ψ/(∂Σ_I)(∂Σ_J)
          PX1P_SIG=MATMUL(PX1P_X1P,X1P_SIG)
          PX2P_SIG=MATMUL(PX2P_X2P,X2P_SIG)
C         #   PART-A
          PHI1P_PX1P_SIG=MATMUL(PHI1P_PX1P_PX1P,PX1P_SIG)
          PHI2P_PX2P_SIG=MATMUL(PHI2P_PX2P_PX2P,PX2P_SIG)

          PHI1P_SIG_SIG=MATMUL(TRANSPOSE(PX1P_SIG),PHI1P_PX1P_SIG)
          PHI2P_SIG_SIG=MATMUL(TRANSPOSE(PX2P_SIG),PHI2P_PX2P_SIG)

          PSI_SIG_SIG=PHI1P_SIG_SIG+PHI2P_SIG_SIG

C         #   PART-B
          PHI1P_X1P_X1P=0.D0
          PHI2P_X2P_X2P=0.D0
          DO 1000 IA=1, 2
          DO 1000 IR=1, 3
          DO 1000 IB=1, 3
              PHI1P_X1P_X1P(IR,IB)=PHI1P_X1P_X1P(IR,IB)
     1       +PHI1P_PX1P(IA)*PX1P_X1P_X1P(IA,IR,IB)

              PHI2P_X2P_X2P(IR,IB)=PHI2P_X2P_X2P(IR,IB)
     1       +PHI2P_PX2P(IA)*PX2P_X2P_X2P(IA,IR,IB)
1000       CONTINUE

          PHI1P_X1P_SIG=MATMUL(PHI1P_X1P_X1P,X1P_SIG)
          PHI1P_SIG_SIG=MATMUL(TRANSPOSE(X1P_SIG), PHI1P_X1P_SIG)

          PHI2P_X2P_SIG=MATMUL(PHI2P_X2P_X2P,X2P_SIG)
          PHI2P_SIG_SIG=MATMUL(TRANSPOSE(X2P_SIG), PHI2P_X2P_SIG)

          PSI_SIG_SIG= PSI_SIG_SIG
     1                +PHI1P_SIG_SIG+PHI2P_SIG_SIG

C         (11) ∂∂Φ/∂Σ∂Σ = DDFDDS
          PHI_SIG_SIG=0.D0
          DO 1001 I=1, 3
          DO 1001 J=1, 3
              PHI_SIG_SIG(I,J)= PHI_SIG_SIG(I,J)
     1       +PHI_SIG_PSI(J)*PSI_SIG(I)+PHI_PSI*PSI_SIG_SIG(I,J)
1001       CONTINUE

      END IF
      RETURN
      END SUBROUTINE GRAD_YLD2000
C-----------------------------------------------------------------------
      SUBROUTINE YLD2004_18P(SIG, PHI)
C     # THIS SUBROUTINES RETURNS THE EQUIVALENT STRESS OF YLD2004-18P MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), S1P(NDIM3), S2P(NDIM3),  
     1          PS1P(NDIM1), PS2P(NDIM1), H1P(NDIM1), H2P(NDIM1),
     2          XL1P(NDIM3, NDIM3), XL2P(NDIM3, NDIM3)

      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KYLD2K4_1/CM, C1, C2, C3, C4, C5, C6, C7, C8, C9
      COMMON /KYLD2K4_2/C10,C11,C12,C13,C14,C15,C16,C17,C18

C     [1] LINEAR TRANSFORMATION
      CALL LTM_YLD2004_18P(SIG, XL1P, XL2P, S1P, S2P)

C     [2] PRINCIPAL STRESS
      CALL PRIN_YLD2004_18P(S1P, H1P, PS1P)
      CALL PRIN_YLD2004_18P(S2P, H2P, PS2P)

C     [3] ANISOTROPIC YIELD FUNCTION
C         Φ   =(1/4)*Ψ**(1/AM)
      PSI=0.D0
      DO 1100 IP=1, NDIM1
      DO 1100 IQ=1, NDIM1
          PSI=PSI+DABS(PS1P(IP)-PS2P(IQ))**CM
1100  CONTINUE

      PHI= (PSI/4.D0)**(1.D0/CM)
      RETURN
      END SUBROUTINE YLD2004_18P
C-----------------------------------------------------------------------
      SUBROUTINE LTM_YLD2004_18P(SIG, XL1P, XL2P, S1P, S2P)
C     # THIS SUBROUTINE RETURNS LIEAR-TRANSFORMED STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), S1P(NDIM3), S2P(NDIM3)
      DIMENSION C1P(NDIM3, NDIM3), C2P(NDIM3, NDIM3), T(NDIM3, NDIM3),
     1          XL1P(NDIM3, NDIM3), XL2P(NDIM3, NDIM3)
      DIMENSION DEV_SIG(NDIM6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KYLD2K4_1/CM, C1, C2, C3, C4, C5, C6, C7, C8, C9
      COMMON /KYLD2K4_2/C10,C11,C12,C13,C14,C15,C16,C17,C18

      C1P=0.D0
      C1P(1,2)=-C1
      C1P(1,3)=-C2
      C1P(2,1)=-C3
      C1P(2,3)=-C4
      C1P(3,1)=-C5
      C1P(3,2)=-C6
      C1P(4,4)= C9
      C1P(5,5)= C8
      C1P(6,6)= C7

      C2P=0.D0
      C2P(1,2)=-C10
      C2P(1,3)=-C11
      C2P(2,1)=-C12
      C2P(2,3)=-C13
      C2P(3,1)=-C14
      C2P(3,2)=-C15
      C2P(4,4)= C18
      C2P(5,5)= C17
      C2P(6,6)= C16

      CALL DEVIATORIC_TENS(T)

      XL1P=MATMUL(C1P,T)
      XL2P=MATMUL(C2P,T)
      S1P=MATMUL(XL1P,SIG)
      S2P=MATMUL(XL2P,SIG)

      RETURN
      END SUBROUTINE LTM_YLD2004_18P
C-----------------------------------------------------------------------
      SUBROUTINE PRIN_YLD2004_18P(S0, H, PS)
C     #   THIS SUBROUTINE RETURNS THE PRINCIPAL STRESSES USED IN YLD2004 MODEL.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION S0(NDIM3), S(6), PS(3), H(3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      PARAMETER(PI=DACOS(-1.D0))

      S=0.D0
      IF(NDIM3 .EQ. 3) THEN
          S(1:2)= S0(1:2)
          S(4)= S0(3)
      ELSE
          S=S0
      END IF
          
C     #   S(6)={Σ_XX, Σ_YY, Σ_ZZ, Σ_XY, Σ_XZ, Σ_YZ)
C     [1] CALCULATE INVARIANTS
      H(1)=(S(1)+S(2)+S(3))/3.D0
      H(2)=(S(4)**2+S(5)**2+S(6)**2-S(1)*S(2)-S(1)*S(3)-S(2)*S(3))/3.D0
      H(3)=(2.D0*S(4)*S(5)*S(6)+S(1)*S(2)*S(3)-S(1)*S(6)**2-S(2)*S(5)**2
     1     -S(3)*S(4)**2)/2.D0

C     [2] PARAMETERS FROM INVARIANTS
      P=H(1)**2+H(2)
      Q=(2.D0*H(1)**3+3.D0*H(1)*H(2)+2.D0*H(3))/2.D0

      IF(DABS(Q/(P**(3.D0/2.D0))) .GT. 1.D0) THEN
          TH=0.D0
      ELSE
          TH=DACOS(Q/(P**(3.D0/2.D0)))
      END IF

C     [3] PRINCIPAL STRESSES
      PS(1)=2.D0*DSQRT(H(1)**2+H(2))*DCOS(TH/3.D0)+H(1)
      PS(2)=2.D0*DSQRT(H(1)**2+H(2))*DCOS((TH+4.D0*PI)/3.D0)+H(1)
      PS(3)=2.D0*DSQRT(H(1)**2+H(2))*DCOS((TH+2.D0*PI)/3.D0)+H(1)

      RETURN
      END SUBROUTINE PRIN_YLD2004_18P
C-----------------------------------------------------------------------
      SUBROUTINE GRAD_YLD2004(SIG, PHI, PHI_SIG, PHI_SIG_SIG)
C     # THIS SUBROUTINE RETURNS ANALYTICAL DERIVATIVES OF YLD2004-18P MODEL.
C       THIS PART IS UNDER DEVELOPMENT (18.06.18).
      IMPLICIT REAL*8(A-H, O-Z)

      DIMENSION SIG(NDIM3), PHI_SIG(NDIM3), PHI_SIG_SIG(NDIM3,NDIM3),
     1          PSI_SIG_SIG(NDIM3,NDIM3), PHI_SIG_PSI(NDIM3)
      DIMENSION S1P(NDIM3), S2P(NDIM3), PS1P(NDIM1), PS2P(NDIM1),
     1          XL1P(NDIM3,NDIM3), XL2P(NDIM3,NDIM3), H1P(NDIM1),
     2          H2P(NDIM1), PSI_PS1P(NDIM1), PSI_PS2P(NDIM1),
     3          PS1P_H1P(NDIM1,NDIM1), PS2P_H2P(NDIM1,NDIM1),
     4          H1P_S1P(NDIM1,NDIM3), H2P_S2P(NDIM1,NDIM3),
     5          S1P_SIG(NDIM3,NDIM3), S2P_SIG(NDIM3,NDIM3),
     6          PSI_SIG(NDIM3)
      DIMENSION PSI_PS1P_PS1P(NDIM1,NDIM1), PSI_PS2P_PS2P(NDIM1,NDIM1),
     1         PSI_PS2P_PS1P(NDIM1,NDIM1), PSI_PS1P_PS2P(NDIM1,NDIM1),
     2         PS1P_SIG(NDIM1,NDIM3), PS2P_SIG(NDIM1,NDIM3),
     3         PSI_SIG_PS1P(NDIM3,NDIM1), PSI_SIG_PS2P(NDIM3,NDIM1),
     4         PS1P_H1P_H1P(NDIM1,NDIM1,NDIM1),
     5         PS2P_H2P_H2P(NDIM1,NDIM1,NDIM1),H1P_SIG(NDIM1,NDIM3),
     6         H2P_SIG(NDIM1,NDIM3),PS1P_SIG_H1P(NDIM1,NDIM3,NDIM1),
     7         PS2P_SIG_H2P(NDIM1,NDIM3,NDIM1)
      DIMENSION H1P_S1P_S1P(NDIM1,NDIM3,NDIM3),
     1 H2P_S2P_S2P(NDIM1,NDIM3,NDIM3), H1P_SIG_S1P(NDIM1,NDIM3,NDIM3),
     2 H2P_SIG_S2P(NDIM1,NDIM3,NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      COMMON /KYLD2K4_1/CM, C1, C2, C3, C4, C5, C6, C7, C8, C9
      COMMON /KYLD2K4_2/C10,C11,C12,C13,C14,C15,C16,C17,C18

C     [0] VALUES FROM YLD2004_18P
      CALL LTM_YLD2004_18P(SIG, XL1P, XL2P, S1P, S2P)
      CALL PRIN_YLD2004_18P(S1P, H1P, PS1P)
      CALL PRIN_YLD2004_18P(S2P, H2P, PS2P)

C     [1] INITIALIZATION
      PHI_SIG=0.D0
      PHI_SIG_SIG=0.D0

C     [2] THE FIRST DERIVATIVE OF YLD2004_18P MODEL
C     (1) ∂Φ/∂Ψ
      PHI_PSI=1.D0/(4.D0*CM)*PHI**(1.D0-CM)
C     (2) ∂Ψ/∂Σ
C     (2.1) ∂Ψ/∂PS1P & ∂Ψ/∂PS2P
      PSI_PS1P=0.D0
      PSI_PS2P=0.D0
      DO 1101 IP=1, NDIM1
      DO 1101 IK=1, NDIM1
          PSI_PS1P(IP)= PSI_PS1P(IP)+CM*(PS1P(IP)-PS2P(IK))
     1                 *DABS(PS1P(IP)-PS2P(IK))**(CM-2.D0)
          PSI_PS2P(IP)= PSI_PS2P(IP)-CM*(PS1P(IK)-PS2P(IP))
     1                 *DABS(PS1P(IK)-PS2P(IP))**(CM-2.D0)
1101   CONTINUE
C     (2.2) ∂PS1P/∂H1P & ∂PS2P/∂H2P
      CALL YLD2K4_PS_H(PS1P, H1P, PS1P_H1P)
      CALL YLD2K4_PS_H(PS2P, H2P, PS2P_H2P)

C     (2.3) ∂H1P/∂S1P & ∂H2P/∂S2P
      CALL YLD2K4_H_S(S1P, H1P_S1P)
      CALL YLD2K4_H_S(S2P, H2P_S2P)

C     (2.4) ∂S1P/∂Σ, ∂S2P/∂Σ
      S1P_SIG=XL1P
      S2P_SIG=XL2P

C     (3) ∂Ψ/∂Σ
      PSI_SIG=0.D0
      DO 1102 II=1,NDIM3
      DO 1102 IR=1,NDIM3
      DO 1102 IQ=1,NDIM1
      DO 1102 IP=1,NDIM1
          PSI_SIG(II)= PSI_SIG(II)
     1                +PSI_PS1P(IP)*PS1P_H1P(IP,IQ)*H1P_S1P(IQ,IR)
     2                *S1P_SIG(IR,II)
     1                +PSI_PS2P(IP)*PS2P_H2P(IP,IQ)*H2P_S2P(IQ,IR)
     2                *S2P_SIG(IR,II)
1102   CONTINUE
C     (6.4) ∂Φ/∂Σ
      PHI_SIG=PHI_PSI*PSI_SIG
C-----------------------------------------------------------------------
C     [3] THE SECOND DERIVATIVE OF YLD2004_18P MODEL
      IF(SUA_PAR .NE. 1.) THEN
C     (1) ∂∂Φ/∂Σ∂Ψ
      PHI_SIG_PSI=((1.D0-CM)/PHI)*PHI_PSI*PHI_SIG

C     (2) ∂∂Ψ/∂Σ∂Σ
C     (2.A)
C     (2.A.1) ∂∂Ψ/∂PS1P∂PS1P & ∂∂Ψ/∂PS2P∂PS2P
      PSI_PS1P_PS1P=0.D0
      PSI_PS2P_PS2P=0.D0
      DO 1103 IP=1, NDIM1
      DO 1103 IM=1, NDIM1
      DO IK=1, NDIM1
          IF(IP .EQ. IM) THEN
              PSI_PS1P_PS1P(IM,IP)= PSI_PS1P_PS1P(IM,IP)+CM*(CM-1.D0)
     1                         *DABS(PS1P(IM)-PS2P(IK))**(CM-2.D0)
              PSI_PS2P_PS2P(IM,IP)= PSI_PS2P_PS2P(IM,IP)+CM*(CM-1.D0)
     1                         *DABS(PS1P(IK)-PS2P(IM))**(CM-2.D0)
          ELSE
              PSI_PS1P_PS1P(IM,IP)= 0.D0
              PSI_PS2P_PS2P(IM,IP)= 0.D0
          END IF
      END DO
C     (2.A.2) ∂∂Ψ/∂PS2P∂PS1P & ∂∂Ψ/∂PS1P∂PS2P
      PSI_PS2P_PS1P(IM,IP)= CM*(1.D0-CM)
     1                     *DABS(PS1P(IM)-PS2P(IP))**(CM-2.D0)
      PSI_PS1P_PS2P(IM,IP)= PSI_PS2P_PS1P(IM,IP)
1103   CONTINUE
C     (2.A.3) ∂PS1P/∂Σ & ∂PS2P/∂Σ
      PS1P_SIG=0.D0
      PS2P_SIG=0.D0
      DO 1104 IP=1,NDIM1
      DO 1104 IQ=1,NDIM1
      DO 1104 IR=1,NDIM3
      DO 1104 II=1,NDIM3
          PS1P_SIG(IP,II)= PS1P_SIG(IP,II)
     1                    +PS1P_H1P(IP,IQ)*H1P_S1P(IQ,IR)*S1P_SIG(IR,II)
          PS2P_SIG(IP,II)= PS2P_SIG(IP,II)
     1                    +PS2P_H2P(IP,IQ)*H2P_S2P(IQ,IR)*S2P_SIG(IR,II)
1104   CONTINUE
C     (2.A.4) ∂∂Ψ/∂Σ∂PS1P & ∂∂Ψ/∂Σ∂PS2P
      PSI_SIG_PS1P= MATMUL(PSI_PS1P_PS1P, PS1P_SIG)
     1             +MATMUL(PSI_PS2P_PS1P, PS2P_SIG)
      PSI_SIG_PS2P= MATMUL(PSI_PS2P_PS2P, PS2P_SIG)
     1             +MATMUL(PSI_PS1P_PS2P, PS1P_SIG)

      PSI_SIG_SIG= MATMUL(PSI_SIG_PS1P, PS1P_SIG)
     1              +MATMUL(PSI_SIG_PS2P, PS2P_SIG)

C     (2.B)
C     (2.B.1) ∂∂PS1P/∂H1P∂H1P, ∂∂PS2P/∂H2P∂H2P
      CALL YLD2K4_PS_H_H(PS1P, H1P, PS1P_H1P, PS1P_H1P_H1P)
      CALL YLD2K4_PS_H_H(PS2P, H2P, PS2P_H2P, PS2P_H2P_H2P)
C     (2.B.2) ∂H1P/∂Σ, ∂H2P/∂Σ
      H1P_SIG=MATMUL(H1P_S1P,S1P_SIG)
      H2P_SIG=MATMUL(H2P_S2P,S2P_SIG)
C     (2.B.3) ∂∂PS1P/∂Σ∂H1P, ∂∂PS2P/∂Σ∂H2P
      PS1P_SIG_H1P=0.D0
      PS2P_SIG_H2P=0.D0
      DO 1105 IP=1,NDIM1
      DO 1105 IJ=1,NDIM3
      DO 1105 IQ=1,NDIM1
      DO 1105 IM=1,NDIM1
          PS1P_SIG_H1P(IP,IJ,IQ)= PS1P_SIG_H1P(IP,IJ,IQ)
     1                           +PS1P_H1P_H1P(IP,IM,IQ)*H1P_SIG(IM,IJ)
          PS2P_SIG_H2P(IP,IJ,IQ)= PS2P_SIG_H2P(IP,IJ,IQ)
     1                           +PS2P_H2P_H2P(IP,IM,IQ)*H2P_SIG(IM,IJ)
1105   CONTINUE
C     (2.B.4) ∂∂Ψ/∂Σ∂Σ
      PSI_SIG_SIG=0.D0
      DO 1106 IJ=1,NDIM3
      DO 1106 II=1,NDIM3
      DO 1106 IR=1,NDIM3
      DO 1106 IQ=1,NDIM1
      DO 1106 IP=1,NDIM1
          PSI_SIG_SIG(IJ,II)= PSI_SIG_SIG(IJ,II)
     1                       +PSI_PS1P(IP)*PS1P_SIG_H1P(IP,IJ,IQ)
     2                       *H1P_S1P(IQ,IR)*S1P_SIG(IR,II)
          PSI_SIG_SIG(IJ,II)= PSI_SIG_SIG(IJ,II)
     1                       +PSI_PS2P(IP)*PS2P_SIG_H2P(IP,IJ,IQ)
     2                       *H2P_S2P(IQ,IR)*S2P_SIG(IR,II)
1106   CONTINUE

C     (2.C) ∂∂H1P/∂Σ∂S1P & ∂∂H2P/∂Σ∂S2P
C     (2.C.1) ∂∂H1P/∂S1P∂S1P & ∂∂H2P/∂S2P∂S2P
      CALL YLD2K4_H_S_S(S1P, H1P_S1P_S1P)
      CALL YLD2K4_H_S_S(S2P, H2P_S2P_S2P)
C     (2.C.2) ∂∂H1P/∂Σ∂S1P & ∂∂H2P/∂Σ∂S2P
      H1P_SIG_S1P=0.D0
      H2P_SIG_S2P=0.D0
      DO 1107 IQ=1,NDIM1
      DO 1107 IM=1,NDIM3
      DO 1107 IR=1,NDIM3
      DO 1107 IJ=1,NDIM3
          H1P_SIG_S1P(IQ,IJ,IR)= H1P_SIG_S1P(IQ,IJ,IR)
     1                          +H1P_S1P_S1P(IQ,IM,IR)*S1P_SIG(IM,IJ)
          H2P_SIG_S2P(IQ,IJ,IR)= H2P_SIG_S2P(IQ,IJ,IR)
     1                          +H2P_S2P_S2P(IQ,IM,IR)*S2P_SIG(IM,IJ)
1107   CONTINUE
C     (2.C.3) ∂∂Ψ/∂Σ∂Σ
      PSI_SIG_SIG=0.D0
      DO 1108 IJ=1,NDIM3
      DO 1108 II=1,NDIM3
      DO 1108 IR=1,NDIM3
      DO 1108 IQ=1,NDIM1
      DO 1108 IP=1,NDIM1
          PSI_SIG_SIG(IJ,II)= PSI_SIG_SIG(IJ,II)
     1                       +PSI_PS1P(IP)*PS1P_H1P(IP,IQ)
     2                       *H1P_SIG_S1P(IQ,IJ,IR)*S1P_SIG(IR,II)
          PSI_SIG_SIG(IJ,II)= PSI_SIG_SIG(IJ,II)
     1                       +PSI_PS2P(IP)*PS2P_H2P(IP,IQ)
     2                       *H2P_SIG_S2P(IQ,IJ,IR)*S2P_SIG(IR,II)
1108   CONTINUE

C     (2.D) ∂∂Φ/∂Σ∂Σ
C     (2.D.1) ∂∂Φ/∂Σ∂Σ
      DO 1109 II=1, NDIM3
      DO 1109 IJ=1, NDIM3
          PHI_SIG_SIG(II,IJ)= PHI_SIG_PSI(IJ)*PHI_SIG(II)
     1                       +PHI_PSI*PSI_SIG_SIG(II,IJ)
1109   CONTINUE
      END IF
      RETURN
      END SUBROUTINE GRAD_YLD2004
C-----------------------------------------------------------------------
      SUBROUTINE YLD2K4_PS_H(PS, H, PS_H)
C     #   THIS SUBROUTINE RETURNS PS1P_H1P OR PS2P_H2P.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION PS(NDIM1), H(NDIM1), PS_H(NDIM1,NDIM1), XI(NDIM1)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      P=DSQRT(H(1)**2+H(2))
      IF(DABS(PS(2)-PS(3)) .EQ. 0.D0) THEN         ! SINGULAR CASE I
          PS_H(1,1)=1.D0+2.D0*H(1)/P
          PS_H(1,2)=1.D0/P
          PS_H(1,3)=0.D0
          PS_H(2,1)=1.D0-H(1)/P
          PS_H(2,2)=-1.D0/(2.D0*P)
          PS_H(2,3)=0.D0
          PS_H(3,:)=PS_H(2,:)
      ELSEIF(DABS(PS(2)-PS(1)) .EQ. 0.D0) THEN     ! SINGULAR CASE II
          PS_H(3,1)=1.D0-2.D0*H(1)/P
          PS_H(3,2)=-1.D0/P
          PS_H(3,3)=0.D0
          PS_H(1,1)=1.D0+H(1)/P
          PS_H(1,2)=1.D0/(2.D0*P)
          PS_H(1,3)=0.D0
          PS_H(2,:)=PS_H(1,:)
      ELSE
      DO IP=1, NDIM1
          XI(IP)=PS(IP)**2-2.D0*H(1)*PS(IP)-H(2)  ! [EQ. A.23][7]
          IF(XI(IP) .EQ. 0.D0) PRINT *, '## XI .EQ. 0.D0'
          PS_H(IP,1)=PS(IP)**2/XI(IP)
          PS_H(IP,2)=PS(IP)/XI(IP)
          PS_H(IP,3)=(2.D0/3.D0)*(1.D0/XI(IP))
      END DO
      END IF
      RETURN
      END SUBROUTINE YLD2K4_PS_H
C-----------------------------------------------------------------------
      SUBROUTINE YLD2K4_H_S(S, H_S)
C     #   THIS SUBROUTINE RETURNS H1P_S1P & H2P_S2P
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION S(NDIM3), H_S(NDIM1,NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
C     #   SIG(6)={Σ_XX, Σ_YY, Σ_ZZ, Σ_YZ, Σ_XZ, Σ_XY)
      H_S=0.D0
C     #   FOR IQ=1
      H_S(1,1:3)=1.D0/3.D0
      H_S(1,4:6)=0.D0
C     #   FOR IQ=2
      H_S(2,1)=-(S(2)+S(3))/3.D0
      H_S(2,2)=-(S(1)+S(3))/3.D0
      H_S(2,3)=-(S(1)+S(2))/3.D0
      H_S(2,4:6)=(2.D0/3.D0)*S(4:6)
C     #   FOR IQ=3
      H_S(3,1)=(S(2)*S(3)-S(4)**2)/2.D0
      H_S(3,2)=(S(1)*S(3)-S(5)**2)/2.D0
      H_S(3,3)=(S(1)*S(2)-S(6)**2)/2.D0
      H_S(3,4)=S(5)*S(6)-S(1)*S(4)
      H_S(3,5)=S(4)*S(6)-S(2)*S(5)
      H_S(3,6)=S(4)*S(5)-S(3)*S(6)
      RETURN
      END SUBROUTINE YLD2K4_H_S
C-----------------------------------------------------------------------
      SUBROUTINE YLD2K4_PS_H_H(PS, H, PS_H, PS_H_H)
C     #   THIS SUBROUTINE RETURNS PS1P_H1P_H1P, PS2P_H2P_H2P
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION PS(NDIM1), XI(NDIM1), XI_H(NDIM1,NDIM1),
     1          PS_H(NDIM1,NDIM1), PS_H_H(NDIM1,NDIM1,NDIM1), H(NDIM1)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      PS_H_H=0.D0
      DO 1111 IP=1, NDIM1
          XI(IP)=PS(IP)**2-2.D0*H(1)*PS(IP)-H(2)
          XI_H(IP,1)=-2.D0*PS(IP)
          XI_H(IP,2)=-1.D0
          XI_H(IP,3)=0.D0
      DO 1111 IQ=1, NDIM1
      DO 1111 IM=1, NDIM1
      IF(PS(2).EQ.PS(3)) THEN         ! SINGULAR CASE I
          PS_H_H(2,IM,IQ)=-(5.D-1)*PS_H_H(1,IM,IQ)
          PS_H_H(3,IM,IQ)= PS_H_H(2,IM,IQ)
      ELSEIF(PS(2).EQ.PS(1)) THEN     ! SINGULAR CASE II
          PS_H_H(2,IM,IQ)=-(5.D-1)*PS_H_H(3,IM,IQ)
          PS_H_H(1,IM,IQ)= PS_H_H(2,IM,IQ)
      ELSE                            ! GENERAL CASE
          PS_H_H(IP,IM,IQ)=-XI(IP)**(-2)*XI_H(IP,IM)*PS_H(IP,IQ)
      END IF
          PS_H_H(IP,IQ,IM)=PS_H_H(IP,IM,IQ)
1111   CONTINUE
      RETURN
      END SUBROUTINE YLD2K4_PS_H_H
C-----------------------------------------------------------------------
      SUBROUTINE YLD2K4_H_S_S(S, H_S_S)
C     #   THIS SUBROUTINE RETURNS H1P_S1P_S1P & H2P_S2P_S2P
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION S(NDIM3), H_S_S(NDIM1,NDIM3,NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      H_S_S=0.D0

      DO 1112 IM=1, NDIM1
      DO 1112 IR=1, NDIM1
          IF(IM .EQ. IR) THEN
              H_S_S(2,IM,IR)=0.D0
          ELSE
              H_S_S(2,IM,IR)=-1.D0/3.D0
          END IF
1112   CONTINUE

      H_S_S(3,2,1)=(5.D-1)*S(3)
      H_S_S(3,3,1)=(5.D-1)*S(1)
      H_S_S(3,4,1)=0.D0
      H_S_S(3,5,1)=0.D0
      H_S_S(3,6,1)=-S(6)

      H_S_S(3,1,2)=H_S_S(3,2,1)
      H_S_S(3,3,2)=(5.D-1)*S(1)
      H_S_S(3,4,2)=0.D0
      H_S_S(3,5,2)=-S(5)
      H_S_S(3,6,2)=0.D0

      H_S_S(3,1,3)=(5.D-1)*S(2)
      H_S_S(3,2,3)=H_S_S(3,3,2)
      H_S_S(3,4,3)=-S(4)
      H_S_S(3,5,3)=0.D0
      H_S_S(3,6,3)=0.D0

      H_S_S(3,1,4)=0.D0
      H_S_S(3,2,4)=0.D0
      H_S_S(3,3,4)=H_S_S(3,4,3)
      H_S_S(3,4,4)=-S(3)
      H_S_S(3,5,4)=S(6)
      H_S_S(3,6,4)=S(5)

      H_S_S(3,1,5)=0.D0
      H_S_S(3,2,5)=-S(5)
      H_S_S(3,3,5)=H_S_S(3,5,3)
      H_S_S(3,4,5)=S(6)
      H_S_S(3,5,5)=-S(2)
      H_S_S(3,6,5)=S(4)

      H_S_S(3,1,6)=H_S_S(3,6,1)
      H_S_S(3,2,6)=H_S_S(3,6,2)
      H_S_S(3,3,6)=H_S_S(3,6,3)
      H_S_S(3,4,6)=H_S_S(3,6,4)
      H_S_S(3,5,6)=H_S_S(3,6,5)
      H_S_S(3,6,6)=-S(1)

      RETURN
      END SUBROUTINE YLD2K4_H_S_S
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.7  DIFFERENTIATION                                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRAD(FLAG, HARD_VAR, SIG, SIG_BAR, DFDS, DDFDDS)
C     #   THIS SUBROUTINE CALCULATES THE DERIVATIVES OF CONSTITUTIVE MODELS.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     #   CHECK STRESS TENSOR
      XSIG=0.D0
      DO I=1,NDIM3
          XSIG=XSIG+SIG(I)*SIG(I)
      END DO
      IF(XSIG.EQ.0.D0) THEN
          DFDS=0.D0
          DDFDDS=0.D0
          RETURN
      END IF

      IF(GRAD_PAR .EQ. 1.) THEN
C     #   ANALYTICAL DIFFERENTIATION
          IF(HARD_PAR .LT. 3.) THEN
              CALL GRAD_YLD(SIG, SIG_BAR, DFDS, DDFDDS)
          ELSEIF(HARD_PAR .EQ. 3.) THEN
              CALL GRAD_HAH11(HARD_VAR, SIG, SIG_BAR, DFDS, DDFDDS)
          ELSE
              WRITE(*,*) 'WRONG OPTION!'
          END IF
      ELSE
C     #   NUMERICAL DIFFERENTIATION
          CALL FDM(FLAG, HARD_VAR, SIG, DFDS, DDFDDS)
      END IF

      RETURN
      END SUBROUTINE GRAD
C-----------------------------------------------------------------------
      SUBROUTINE GRAD_YLD(SIG, SIG_BAR, DFDS, DDFDDS)
C     #   THIS SUBROUTINE CALCULATES THE DERIVATIVE OF YIELD FUNCTION
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION SIG(NDIM3), DFDS(NDIM3), DDFDDS(NDIM3,NDIM3)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      IF(YLD_PAR .EQ. 3.) THEN        ! YLD2000_2D
          CALL GRAD_YLD2000(SIG, SIG_BAR, DFDS, DDFDDS)
      ELSEIF(YLD_PAR .EQ. 4.) THEN    ! YLD2004_18P
          CALL GRAD_YLD2004(SIG, SIG_BAR, DFDS, DDFDDS)
      ELSE
          WRITE(*,*) 'WRONG OPTION!'
      END IF

      RETURN
      END SUBROUTINE GRAD_YLD
C-----------------------------------------------------------------------
      SUBROUTINE FDM(FLAG, HARD_VAR, SIG, DFDS, DDFDDS)
C     # THIS SUBROUTINE CALCULATES PARTIAL DERIVATIVE THROUGH NUMERICAL
C       DIFFERENTIATION WHICH IS FINITE DIFFERENCE METHOD.
C       THREE-POINT MIDPOINT RULE WAS USED.
C     # FLAG=0: YIELD FUNCTION DERIVATIVE
C           =1: EQUIVALENT STRESS DERIVATIVE (HAH OR YLD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3),SIGP(NDIM3),DFDS(NDIM3),DDFDDS(NDIM3,NDIM3),
     1          ASIG(NDIM3,3,NDIM3),ASIG_BAR(NDIM3,3),HM(NDIM3,3,NDIM3),
     2          H(NDIM3),DGDS(3)
      DIMENSION HARD_VAR(NDIM7)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
C     SCALING FACTOR [8]
      XI=0.D0
      DO I=1, NDIM3
          XI=XI+SIG(I)*SIG(I)
      END DO
      XI=1.D0/DSQRT(XI)
C     SCALED STRESS COMPONENTS
      SIGP=SIG*XI

C     # DEFINE THE STEP SIZE
      H=1.D-6

      HM= 0.D0
      DFDS=0.D0
      DDFDDS=0.D0
      DO I=1, NDIM3
          DO J=1, 3
C         HM: MULTIPLIEROF STEP-SIZE, H
              HM(I, J, I)= 2.D0-J*1.D0

C         INPUT MATRIX FOR EQUIVALENT STRESS
C         ASIG=   {[SIG(1)+H1     SIG(2)  SIG(3)] [SIG(1) SIG(2)+H2   SIG(3)] [SIG(1) SIG(2)  SIG(3)+H3   ]}
C                 {[SIG(1)        SIG(2)  SIG(3)] [SIG(1) SIG(2)      SIG(3)] [SIG(1) SIG(2)  SIG(3)      ]}
C                 {[SIG(1)-H1     SIG(2)  SIG(3)] [SIG(1) SIG(2)-H2   SIG(3)] [SIG(1) SIG(2)  SIG(3)-H3   ]}
              ASIG(I, J, :) = SIGP
              DO K=1, NDIM3
                  ASIG(I, J, K)= SIGP(K)+HM(I,J,K)*H(K)
              END DO
C         CALCULATE EQUIVALENT STRESS
              IF(FLAG .EQ. 0.D0) THEN
                  CALL YLD(ASIG(I,J,:), ASIG_BAR(I,J))
              ELSEIF(FLAG .EQ. 1.D0) THEN
                  CALL EQV_SIG(HARD_VAR, ASIG(I,J,:), ASIG_BAR(I,J))
              END IF
          END DO
C         CALCULATE THE FIRST-ORDER PARTIAL DERIVATIVE
          CALL FDM_CAL(H(I), ASIG_BAR(I,:), DFDS(I))
      END DO
C-----------------------------------------------------------------------
C         CALCULATE THE SECOND-ORDER PARTIAL DERIVATIVE
      IF(SUA_PAR.NE.1. .AND. FLAG.EQ.1.D0) THEN
      INDX=0.D0
      DO I=1, NDIM3
          DO J=1, NDIM3
              IF(I .EQ. J) THEN
C         #   I=J: (1,1), (2,2), (3,3)
                  DDFDDS(I, J)=
     1             (ASIG_BAR(I,1)-2.D0*ASIG_BAR(I,2)
     2            +ASIG_BAR(I,3))/(H(I)*H(I))
              ELSEIF(I .LT. J) THEN
C         #   I<J: (1,2), (1,3), (2,3)
                  INDX=INDX+1.D0
                  CALL FDM_NUMER(HARD_VAR,INDX,I,J,H,HM,SIGP,DGDS)
                  CALL FDM_CAL(H(I), DGDS, DDFDDS(I, J))
                  DDFDDS(J, I)= DDFDDS(I, J)
              END IF
          END DO
      END DO
      END IF
C-----------------------------------------------------------------------
      DDFDDS=XI*DDFDDS

      RETURN
      END SUBROUTINE FDM
C-----------------------------------------------------------------------
      SUBROUTINE FDM_NUMER(HARD_VAR, INDX, I, J, H, HM, SIG, DGDS)
C     #   THIS SUBROUTINE RETURNS G-FUNCTION REQUIRED TO CALCULATE
C         SECOND-ORDER DERIVATIVE.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIG(NDIM3),ASIG(9,NDIM3),ASIG_BAR(9),H(NDIM3),
     1          HM(NDIM3,3,NDIM3),HM2ND(9,NDIM3),DGDS(3)
      DIMENSION HARD_VAR(NDIM7)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     #   INITIALIZATION
      HM2ND=0.D0
      ASIG=0.D0
      ASIG_BAR=0.D0
      DGDS=0.D0

C     #   H-VALUE MATRIX
      DO J0=1, 3
          N0=3*J0-2
          N1=3*J0
          HM2ND(N0:N1,I)=HM(I,J0,I)
          HM2ND(N0:N1,J)=HM(J,:,J)
      END DO

      DO J1=1, 9
C     INPUT MATRIX FOR EQUIVALENT STRESS FOR THREE-POINT RULE
C     ASIG=   {[SIG(1)+H1  SIG(2)+H2  SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)+H2  SIG(3)+H3]} 1
C             {[SIG(1)+H1  SIG(2)     SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)+H2  SIG(3)   ]} 2
C             {[SIG(1)+H1  SIG(2)-H2  SIG(3)] [SIG(1)+H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)+H2  SIG(3)-H3]} 3
C             {[SIG(1)     SIG(2)+H2  SIG(3)] [SIG(1)     SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)     SIG(3)+H3]} 4
C             {[SIG(1)     SIG(2)     SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)     SIG(3)   ]} 5
C             {[SIG(1)     SIG(2)-H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)     SIG(3)-H3]} 6
C             {[SIG(1)-H1  SIG(2)+H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)+H3] [SIG(1)  SIG(2)-H2  SIG(3)+H3]} 7
C             {[SIG(1)-H1  SIG(2)     SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)   ] [SIG(1)  SIG(2)-H2  SIG(3)   ]} 8
C             {[SIG(1)-H1  SIG(2)-H2  SIG(3)] [SIG(1)-H1  SIG(2)  SIG(3)-H3] [SIG(1)  SIG(2)-H2  SIG(3)-H3]} 9
          DO K1=1, NDIM3
              ASIG(J1,K1)= SIG(K1) + HM2ND(J1,K1)*H(K1)
          END DO

C     CALCULATE EQUIVALENT STRESS
          CALL EQV_SIG(HARD_VAR, ASIG(J1,:), ASIG_BAR(J1))
      END DO

C     CALCULATE G-FUNCTION FOR THE SECOND-DERIVATIVE
      DO N2=1, 3
          N3=3*N2-2
          N4=3*N2
          CALL FDM_CAL(H(J), ASIG_BAR(N3:N4), DGDS(N2))
      END DO

      RETURN
      END SUBROUTINE FDM_NUMER
C-----------------------------------------------------------------------
      SUBROUTINE FDM_CAL(H, F, DF)
C     # THIS SUBROUTINE CALCULATES THREE MID-POINT RULE OF FDM.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(3)

      DF=(F(1)-F(3))/(2.D0*H)

      RETURN
      END SUBROUTINE FDM_CAL
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUB.8   OPERATORS                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLOW_RULE(HARD_VAR, SIG, EQPLAS, DGDS)
C     # THIS SUBROUTINE CALCULATES THE GRADIENT YIELD SURFACE IN TERMS OF
C       STRESS TENSOR.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HARD_VAR(NDIM7), SIG(NDIM3), DGDS(NDIM3)
      DIMENSION DEV_SIG(NDIM4), DEV_SIG0(NDIM3)
      DIMENSION T(NDIM3,NDIM3), DDGDDS(NDIM3,NDIM3)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     [1] GRADIENT OF YIELD SURFACE
      CALL FLOW_STRESS(HARD_VAR, EQPLAS, FLOW_SIG, DHDE)
      CALL PRESSURE_EFFECT(PEC, SIG, FLOW_SIG, SIG_KK, BVAL, CVAL)
      IF(FLOOR(HARD_PAR) .GE. 4.) THEN
          CALL DEVIATORIC(SIG/BVAL, DEV_SIG)
          CALL DEVIATORIC_TRANS(DEV_SIG, DEV_SIG0)
          CALL GRAD(0.D0, HARD_VAR, DEV_SIG0, SIG_BAR, DGDS, DDGDDS)
          CALL DEVIATORIC_TENS(T)
          DGDS=MATMUL(T, DGDS)
      ELSE
          CALL EQV_SIG(HARD_VAR, SIG/BVAL, SIG_BAR_P)
          CALL GRAD(0.D0, HARD_VAR, SIG, SIG_BAR_P, DGDS, DDGDDS)
      END IF

      RETURN
      END SUBROUTINE FLOW_RULE
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC(SIG, DEV_SIG)
C     #   THIS SUBROUTINE CALCULATE DEVIATORIC STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION SIG(NDIM3), DEV_SIG(NDIM6), DELTA(NDIM6)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      DEV_SIG=0.D0
      NDI=NDIM6-NDIM2
C     INITIAL DEVIATORIC STRESS TENSOR
      IF(NDIM3 .EQ. 3) THEN ! PLANE STRESS
          DEV_SIG(1)= SIG(1)
          DEV_SIG(2)= SIG(2)
          DEV_SIG(3)= 0.D0
          DEV_SIG(4)= SIG(3)
      ELSE
          DEV_SIG= SIG
      END IF

      SIG_KK=0.D0
      DELTA=0.D0
      DO I=1, NDI
          SIG_KK=SIG_KK+DEV_SIG(I)/3.D0
          DELTA(I)=1.D0
      END DO
      DEV_SIG=DEV_SIG-SIG_KK*DELTA

      RETURN
      END SUBROUTINE DEVIATORIC
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC_TRANS(DEV_SIG, DEV_SIG0)
C     #   THIS SUBROUTINE REMOVES DEV_SIG33.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION DEV_SIG(NDIM6), DEV_SIG0(NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8
      
      NDI=NDIM3-NDIM2
C     SPHERICAL COMPONENTS
      DO I=1, NDI
          DEV_SIG0(I)=DEV_SIG(I)
      END DO
C     SHEAR COMPONENTS
      IF(NDIM3 .EQ. 3) THEN !PLANE STRESS
          DEV_SIG0(3)=DEV_SIG(4)
      ELSE IF(NDIM3 .EQ. 6) THEN
          DO I=NDI+1,NDIM3
              DEV_SIG0(I)=DEV_SIG(I)
          END DO
      ELSE
          WRITE(*,*) '#ERROR: DIMENSION OF STRESS TENSOR'
      END IF

      RETURN
      END SUBROUTINE DEVIATORIC_TRANS
C-----------------------------------------------------------------------
      SUBROUTINE DEVIATORIC_TENS(T)
C     #   THIS SUBROUTINE REMOVES DEV_SIG33.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION T(NDIM3, NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      T=0.D0
      IF(NDIM3 .EQ. 3) THEN
          T(1,1)= 2.D0
          T(1,2)= -1.D0
          T(2,1)= -1.D0
          T(2,2)= 2.D0
          T(3,3)= 3.D0
      ELSE
          T(1,1)= 2.D0
          T(1,2)=-1.D0
          T(1,3)=-1.D0
          T(2,1)=-1.D0
          T(2,2)= 2.D0
          T(2,3)=-1.D0
          T(3,1)=-1.D0
          T(3,2)=-1.D0
          T(3,3)= 2.D0
          T(4,4)= 3.D0
          T(5,5)= 3.D0
          T(6,6)= 3.D0
      END IF
      T=T/3.D0

      RETURN
      END SUBROUTINE DEVIATORIC_TENS
C-----------------------------------------------------------------------
      SUBROUTINE NORM_HAH(H,DEV_H)
C     #   THIS SUBROUTINE RETURNS THE NORMALIZED MICRO STRUCTURE DEVIATOR.
C         [EQ.(2)][3] WAS USED.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION H(NDIM6), DEV_H(NDIM6)
      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     #   NUMERATOR
C     XDENOM = H_KL*H_KL
      CALL DOUBLE_DOT(H,H,NDIM6,XDENOM)

      IF (XDENOM .EQ. 0.D0) THEN
          DEV_H=0.D0
          RETURN
      END IF
      HVAL=8.D0/3.D0
      IF(HARD_PAR .EQ. 3.) THEN !HAH11
C         XDENOM = SQRT[(8/3)*H_KL*H_KL]
          XDENOM= DSQRT(HVAL*XDENOM)
      ELSEIF(HARD_PAR .EQ. 4.) THEN !HAH14
C          XDENOM= DSQRT((8.D0/3.D0)*XDENOM)
          XDENOM= DSQRT(HVAL*XDENOM)
      ELSEIF(FLOOR(HARD_PAR) .EQ. 5.) THEN !HAH20
          XDENOM= DSQRT(XDENOM)
      END IF

C     #   NORMALIZED DEVIATOR
C     DEV_H = H^_IJ = H_IJ / SQRT[(8/3)*H_KL*H_KL]
      DEV_H=H/XDENOM

      RETURN
      END SUBROUTINE NORM_HAH
C-----------------------------------------------------------------------
      SUBROUTINE COSX(HARD_VAR, STRESS, DEV_S, COS_X)
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION HARD_VAR(NDIM7), STRESS(NDIM3), DEV_SIG(NDIM6)
      DIMENSION DEV_S(NDIM6), DEV_H(NDIM6)

      COMMON /KOPTION/ HARD_PAR, YLD_PAR, FLOW_PAR, SUA_PAR, GRAD_PAR
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      IF(FLOOR(HARD_PAR).EQ.3.) THEN ! HAH11
          VAL=8.D0/3.D0
      ELSEIF(FLOOR(HARD_PAR).EQ.4.) THEN ! HAH14
          VAL=8.D0/3.D0
      ELSEIF(FLOOR(HARD_PAR).EQ.5.) THEN ! HAH20
          VAL=1.D0
      ELSE
          RETURN
      END IF

C     NORMALIZED DEVIATORIC STRESS
      CALL DEVIATORIC(STRESS, DEV_SIG)
      CALL NORM_HAH(DEV_SIG,DEV_S)
C     MICROSTRUCTURE DEVIATOR
      DO I=1, NDIM6
          DEV_H(I)=HARD_VAR(NDIM5+I)
      END DO

C     STRAIN-PATH CHANGE PARAMETER
      CALL DOUBLE_DOT(DEV_S, DEV_H, NDIM6, COS_X)
      COS_X=VAL*COS_X
      IF(DABS(COS_X).GE.1.D0) COS_X=COS_X/DABS(COS_X)

      RETURN
      END SUBROUTINE COSX
C-----------------------------------------------------------------------
      SUBROUTINE DOUBLE_DOT(A,B,NDIM,DOUBLE)
C     #   A:B = A_IJ*B_IJ = TR(A^T*B)
C         A:B = A11*B11 + A12*B12 + A13*B13 + A21*B21 + A22*B22 + A23*B23 + A31*B31 + A32*B32 + A33*B33
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION A(NDIM), B(NDIM)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

C     DOUBLE = Σ DEV_H(I)*DEV_SIG(I)
      DOUBLE = 0.D0
      NDI= NDIM-NDIM2
      DO I=1,NDI
          DOUBLE = DOUBLE + A(I)*B(I)
      END DO
      DO I=NDI+1, NDIM
          DOUBLE= DOUBLE + 2.D0*A(I)*B(I)
      END DO

      RETURN
      END SUBROUTINE DOUBLE_DOT
C-----------------------------------------------------------------------
      SUBROUTINE DOUBLE_DOT21(A,B,DOUBLE)
C     #   THIS SUBROUTINE CALCULATE DOUBLE DOT PRODUCT FOR NON-DEVIATORIC VARAIBLES
C         IN PLANE STRESS.
      IMPLICIT REAL*8(A-H, O-Z)
      DIMENSION A(NDIM3,NDIM3), B(NDIM3), DOUBLE(NDIM3)
      COMMON /KSIZE/ NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8

      DOUBLE = 0.D0
      NDI= NDIM3-NDIM2
      DO I=1, NDIM3
      DO J=1, NDIM3
          IF(J .LE. NDI) THEN
              DOUBLE(I) = DOUBLE(I) + A(I,J)*B(J)
          ELSE
              DOUBLE(I) = DOUBLE(I) + 2.D0*A(I,J)*B(J)
          END IF
      END DO
      END DO

      RETURN
      END SUBROUTINE DOUBLE_DOT21
C-----------------------------------------------------------------------
      SUBROUTINE INVERSE(A, N, INV_A)
C     #   THIS SUBROUTINE RETURNS THE INVERSE MATRIX A(NXN) USING LU-DECOMPOSITION.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N), INDX(N)
      REAL*8 INV_A(N,N), LU_A(N,N)

      LU_A= A
      INV_A= 0.D0
      DO I=1,N
          INV_A(I,I)= 1.D0
      END DO
      CALL LUDCMP(LU_A, N, N, INDX, DDCMP)
      DO J=1,N
          CALL LUBKSB(LU_A, N, N, INDX, INV_A(1,J))
      END DO

      RETURN
      END SUBROUTINE INVERSE
C-----------------------------------------------------------------------
      SUBROUTINE LUDCMP (A, N, NP, INDX, D)
C-----  LU decomposition
C-----  Use single precision on cray
C-----  THIS CODE RETURN THE MIXED LU DECOMPOSED MATRIX
C-----  A = L*U = LU'
C-----  LU' = [U11 U12 U13]
C             [L21 U22 U23]
C             [L31 L32 L33]
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.D0
      DO I=1,N
         AAMAX=0.D0

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.D0) THEN
              RETURN
         END IF
         VV(I)=1.D0/AAMAX
      END DO
C     
      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO
      RETURN
      END SUBROUTINE LUDCMP
C-----------------------------------------------------------------------
      SUBROUTINE LUBKSB (A, N, NP, INDX, B)
C-----  Linear equation solver based on LU decomposition
C-----  Use single precision on cray
C----   IN THIS CODE, WE CAN GET THE INVERSE MATRIX OF A BY CALCULATING A MATRIX X
C----   A*X = I -> X = A^-1
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)
C     ## A = L*U
C     ## A*X = (L*U)*X = L*(U*X) = B
C     ## FORWARD SUBSTITUTION ##
C     ## L*Y = B: GET Y(I)!!
C     L(I,I) = 1
      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF
C     ## HERE: B(I) = SUM -> Y(I)
         B(I)=SUM 
      END DO
C     ## BACKWARD SUBSTITUTION ##
C     ## U*X = Y: GET X(I)!!
      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF
C     ## HERE: B(I) = SUM/A(I,I) -> X(I)
         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END SUBROUTINE LUBKSB
C-----------------------------------------------------------------------