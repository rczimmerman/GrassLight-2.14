	SUBROUTINE E05NM(ED5NM,THETAZ)
C     *******************************************************************
C                 Subroutine E05NM Ver 2.14
C
C         Computes from astronomic data provided by the user:
C             1) downwelling irradiance spectrum just below the sea surface
C                 at 5 nm resolution
C             2) photoperiod
C
C         Usage: CALL E05NM(ED5NM,THETAZ)
C
C         Definitions of passed variables
C             ED5NM -  Downwelling spectral irradiance below the water surface (5 nm resolution)
C             THETAZ - Sun angle in air, as reqired for implementation of the Lee et al. (2005) model for Kd
C
C     Subroutines required
C     --------------------------
C         LINFIT - Linear regression calculation, appended to this file
C
C     Files required by E05NM:
C     -----------------------
C         None
C
C     Files created by E05NM:
C     ------------------------
C         SFCIRR.TXT - Ssurface irradiance spectrum in water, at 5 nm intervals, for use by OWQINIT - CHUCK: IS THIS IN WATER OR ABOVE WATER?
C
C
C
C		Ver 1.0 Created: 1 May 2012                  by: CLG
C	Date		     Modification                       Version   By
C	-------------- ----------------------------------  --------  --------
C       23 Jul 2012   DAYLEN now calculated from          2.04         RCZ
C                     ND and returned to OWQINIT
C       21 May 2015   OZONE calculation corrected         2.12         RCZ
C                     to operate on radians for all
C                     terms
C        4 Jun 2015   Modified OZONE calculation to       2.13         CLG
C                     calculate climatology for all
C                     hemispheres.  Removed disclaimer
C                     comments limiting calculation to
C                     northern hemisphere
C       22 Dec 2015   Corrected calculation of MTHETA     2.14        RCZ
C                     replacing "-" with "+"
C
C     *********************************************************************
	INCLUDE 'CommonE0.f'
	REAL ND, MTHETA, MOZTHETA, MPTHETA
      CHARACTER*64 COPYRITE,VER
	DIMENSION X(3),Y(3)
	DIMENSION F0(301),TR(301),TOZ(301),TOX(301),TW(301),TAUA(301),
     1     TA(301),TAA(301),TAS(301),EDD(301),EDS(301),ED(301)
      REAL IR(301),IA(301)
      DIMENSION ED5NM(61),ED10NM(31)
      DATA COPYRITE/'Copyright (c) 2016 by CL Gallegos & RC Zimmerman'/
      DATA VER/'Version 2.14'/
C   S/R to calculate incident irradiance spectrum (W m^-2 nm^-1)
C   just below the water surface.  Implements Gregg and Carder (1990)
C   
C   First section calculates constants used to determine sun position
	ECC=-0.0167
	PI=3.14156265
	PLANK=6.626068E-34
	CLIGHT=299800000.
	GAMMA=2*PI*(JD-1)/365.
	
	DEL=(0.006918-0.399912*COS(GAMMA)+0.070257*SIN(GAMMA)-0.006758*
     &COS(2*GAMMA)+0.000907*SIN(2*GAMMA)-0.002697*COS(3*GAMMA)+0.00148*
     &SIN(3*GAMMA))
        ET=(0.000075+0.001868*COS(GAMMA)-0.032077*SIN(GAMMA)-0.014615*
     &COS(2*GAMMA)-0.04089*SIN(2*GAMMA))*229.18
     	RPHI=PI*XLAT/180.
     	RLONG=PI*XLON/180.
     	OMEGASR=(ACOS(-TAN(RPHI)*TAN(DEL)))
     	OMEGASS=-OMEGASR
     	ND=(180./PI)*2.*OMEGASR/15.   !Hours of sun
     	HSR=12.-ND/2.                !Time of sunrise
     	HSS=12.+ND/2.                !Time of sunset
      DAYLEN=ND               !Photoperiod - pass back to GL
C      WRITE(*,999)DAYLEN           !debug photoperiod
C999   FORMAT('In s/r E05NM, DAYLEN = ',F10.3) !debug photoeriod
     	OMEGA=12.*OMEGASR/(12.-HSR)+(OMEGASR/(HSR-12.))*HR
     	THETAZ=(ACOS(SIN(DEL)*SIN(RPHI)+COS(DEL)*COS(RPHI)*COS(OMEGA)))
     	MTHETA=1/(COS(THETAZ)+0.50572*(96.07995-180.*THETAZ/PI)**
     1     (-1.6364))
      MOZTHETA=1.0035/SQRT((COS(THETAZ))**2.+0.007)
      MPTHETA=MTHETA*PRESS/1013.25
C   Ozone climatology section.  Default coefficients to N. America
C   Coefficients from Table 1 in Van Heuklon (1979, Solar Energy, p. 66)
      OZ_A=150.0
      OZ_BETA=1.28
      OZ_C=40.0
      OZ_F=-30.0
      OZ_G=20.0
      OZ_H=3.0
      IF(XLON.GT.0.) THEN 
          OZ_I=20.0 
      ELSE OZ_I=0.0
      ENDIF
      IF(XLAT.GT.0.0) GO TO 40
C   Coefficients for Southern hemisphere
      OZ_A=100.0
      OZ_BETA=1.5
      OZ_C=30.0
      OZ_F=152.625
      OZ_H=2.0
      OZ_I=-75.0           
 40     OZONE=(235+(OZ_A+OZ_C*SIN(0.9856*(JD+OZ_F)*PI/180.)+OZ_G*
     1    SIN(OZ_H*PI*(XLON+OZ_I)/180.))*(SIN(OZ_BETA*RPHI))**2)/1000.
C   Section to compute aerosol size distribution
	F=((2.-RH/100.)/(6.*(1-RH/100.)))**(1./3.)
	A1COEF=2000*AM**2.
	Y(1)=LOG(A1COEF)-LOG((0.1/(F*0.03))**2.)-LOG(F)
	A2COEF=5.866*(WM-2.2)
	Y(2)=LOG(A2COEF)-LOG((1/(F*0.24))**2.)-LOG(F)
	A3COEF=0.01527*(W-2.2)*0.05
	Y(3)=LOG(A3COEF)-LOG((10/(F*2.))**2.)-LOG(F)
	X(1)=LOG(0.1)
	X(2)=0.
	X(3)=LOG(10.)
C   Get aerosol Junge coef. and amplitude from slope and intercept 
	CALL LINFIT(X,Y,3,CAMPL,GAMAER)
	CAMPL=EXP(CAMPL)
	CA=3.91/V
	ALPHAE0=-(GAMAER+3.)
	TAUA550=CA*HA
	BETA=TAUA550/(0.55**(-ALPHAE0))
	IF(ALPHAE0.LT.0.) THEN 
	  AVCOSTHETA=0.82
	  ELSEIF(ALPHAE0.GT.1.2) THEN 
	    AVCOSTHETA=0.65
	  ELSE 
	    AVCOSTHETA=-0.1417*ALPHAE0+0.82
	  ENDIF
	B3COEF=LOG(1-AVCOSTHETA) 
	B2COEF=B3COEF*(0.0783+B3COEF*(-0.3824-0.5874*B3COEF))
	B1COEF=B3COEF*(1.459+B3COEF*(0.1595+0.4129*B3COEF))
	FA=1-0.5*EXP((B1COEF+B2COEF*COS(THETAZ))*COS(THETAZ))
	OMEGA_A=(-0.0032*AM+0.972)*EXP(0.000306*RH)
C   Reflectance section
	IF(W.GT.7.) THEN 
	  CD=(0.49+0.065*W)*0.001  !Drag coefficient
	  ELSE 
	    CD=(0.62+1.56*W)*0.001
	  ENDIF
	RHOA=1200.
	IF(W.GT.7.) THEN
	  RHOF=(0.000045*RHOA*CD-0.00004)*W**2.
	  ELSEIF(W.LT.4.) THEN
	  RHOF=0.
	  ELSE
	  RHOF=0.000022*RHOA*CD*W**2.-0.0004
	  ENDIF
	RHODSP=0.0253*EXP((-0.000714*W+0.0618)*((180.*THETAZ/PI)-40.))
	IF(W.GT.4.) THEN 
	  RHOSSP=0.057
	  ELSE 
	  RHOSSP=0.066
	  ENDIF
	RHOD=RHODSP+RHOF
	RHOS=RHOSSP+RHOF
C   Begin calculation of flux and transmittance arrays at 1 nm intervals
	DO 100 I=1,301
	WL=399+FLOAT(I)
	F0(I)=H0(I)*(1+ECC*COS(2*PI*(JD-3)/365.))**2
	TR(I)=EXP(-(MPTHETA)/(0.0000000001156406*WL**4
     &-0.000001335*WL**2))
        TOZ(I)=EXP(-AOZ(I)*OZONE*MOZTHETA)
        TOX(I)=EXP(-((1.41*AOX(I)*MPTHETA)/(1+118.3*AOX(I)*
     &MPTHETA)**0.45))
        TW(I)=EXP(-((0.238*AW_INC(I)*WV*MTHETA)/(1+20.07*AW_INC(I)*WV*
     &MTHETA)**0.45))
        TAUA(I)=BETA*(WL/1000.)**(-ALPHAE0)
        TA(I)=EXP(-TAUA(I)*MTHETA)
        TAA(I)=EXP(-(1.-OMEGA_A)*TAUA(I)*MTHETA)
        TAS(I)=EXP(-OMEGA_A*TAUA(I)*MTHETA)
        IR(I)=F0(I)*COS(THETAZ)*TOZ(I)*TOX(I)*TW(I)*TAA(I)*
     &(1-TR(I)**0.95)*0.5
        IA(I)=F0(I)*COS(THETAZ)*TOZ(I)*TOX(I)*TW(I)*TAA(I)*TR(I)**1.5*
     &(1-TAS(I))*FA
        EDD(I)=F0(I)*COS(THETAZ)*TR(I)*TA(I)*TOZ(I)*TOX(I)*TW(I)*
     &(1-RHOD)
        EDS(I)=(IR(I)+IA(I))*(1-RHOS)
        ED(I)=EDD(I)+EDS(I)
c        write(*,8000)i,edd(i),eds(i)
c8000    format(1x,i5,2e12.6)
100	CONTINUE
C  Average ED spectrum at 5 nm intervals from 403 to 697 nm
	ED5NM(1)=ED(1)
	ED5NM(61)=ED(301)
      XWAV=400.
      OPEN(UNIT=99,FILE='SFCIRR.TXT',STATUS='REPLACE')
      WRITE(99,5000)XLAT,XLON,JD
5000  FORMAT(2G15.3,I10)
      WRITE(99,9000)XWAV,ED5NM(1)
	DO 200 I=2,60
          XWAV=XWAV+5.
	    SUMED=0.
	    JSTART=5*(I-1)-2
          JEND=JSTART+4
	    DO 150 J=JSTART,JEND,1
	        SUMED=SUMED+ED(J)
 150	    CONTINUE
 	    ED5NM(I)=SUMED/5.
          WRITE(99,9000)XWAV,ED5NM(I)
9000      FORMAT(1X,G12.4,1X,G12.4)          
 200  CONTINUE 
      XWAV=XWAV+5
      WRITE(99,9000)XWAV,ED5NM(61)
      CLOSE(99)
 1000	FORMAT(10X,2G12.4)
 1001	FORMAT(1X,3G12.4)	
	RETURN
      END
C  End of subroutine 
	SUBROUTINE LINFIT(X,Y,NPTS,A,B)
	DIMENSION X(3),Y(3)
 11	SUM=0.
 	SUMX=0.
 	SUMY=0.
 	SUMX2=0.
 	SUMXY=0.
 	SUMY2=0.
 21	DO 50 I=1,3
 	X1=X(I)
 	Y1=Y(I)
 	SUM=SUM+1.
 	SUMX=SUMX+X1
 	SUMY=SUMY+Y1
 	SUMX2=SUMX2+X1*X1
 	SUMXY=SUMXY+X1*Y1
 	SUMY2=SUMY2+Y1*Y1
 50	CONTINUE
 	DELTA=SUM*SUMX2-SUMX*SUMX
 	A=(SUMX2*SUMY-SUMX*SUMXY)/DELTA  !Intercept
 	B=(SUMXY*SUM-SUMX*SUMY)/DELTA    !Slope	
 	RETURN
 	END 	 		
