       SUBROUTINE OWQINIT(IRRFILE)
C     *********************************************************************
C                     Subroutine OWQINIT Ver 2.14
C
C     Subroutine OWINIT - Optical Water Quality Initialization
C     derived from CallingShell original program for
C     incorporation into GrassLight (Ver 2).  Originally created by Charles L.
C     Gallegos, modified by Richard C. Zimmerman for manual input of some
C     variables and input file names.
C
C     Computes spectral Kd and irradiance file for use in GrassLight
C     from user-supplied inputs. Current implementation provides 5 nm and
C     10nm resolution files.
C     
C     Data passed to OWQINIT from GL via call stmt:
C     ----------------------------------
C         IRRFILE - Name of new ".IRR" file to be created
C     
C     Subroutines and Functions Required:
C     -----------------------------------
C         CommonE0 - Defines common variables between CGWQ and other subroutines
C         CommonKd - Defines common variables between CGWQ and other subroutines
C         E05NM    - Calculates incident irradiance spectrum (W m^-2 nm^-1)
C                     just below the water surface using Gregg and Carder (1990)
C         SPECKD   - Calculates absorption components, scattering, and attenuation
C         ATTEN    - Attenuates Ed to top of seagrass canopy using Beer's Law
C
C		Originally Created: 1 May 2012                  by: CLG
C         Modified for GrassLight: 25 May 2012            by: RCZ
C
C
C
C         Files required by OWQINIT_210:
C         ----------------------------------
C             SiteData.txt - provides information on geographic position (latitude, longitude), 
C                 date, time, water depth, CDOM absorption at ___nm (m-1), Chl a concentration (µg m-3) 
C                 and turbidity (NTU)
C             IOPCnsts.txt - provides constants required to determine Kd() from water quality parameters
C             Atmos.txt - provides atmospheric constants for determining downwelling plane irradiance at the surface of the earth
C             E0InputArrays.txt - radiant flux at top of atmosphere, and absorption spectra of atmospheric components
C             CANHTZ.TXT - canopy height and depth data; created in GL main routine
C
C         Files created by OWQINIT_210:  
C         --------------------------------
C             CGIRR05.IRR - Top of canopy irradiance at 5 nm intervals.  Not otherwise used by GL 2.10
C             OWQ10NM.IRR - Top of canopy irradiance at 10 nm intervals; used by EDCLC.  Name can be 
C                 modified by the user in this routine
C
C         Files modified by OWQINIT_210:
C         --------------------------------
C             SiteData.txt - based on user input
C             IOPCnsts.txt - based on user input
C             Atmos.txt - based on user input
C
C             
C
C	Date		     Modification		      Version   By
C	-------------- ----------------------------------  --------  --------
C       15 Jun 2012   Reinstated write/read of CANHT,DEPTH   2.0    CLG
C       15 Jun 2012   real declaration for NUCDOM, etc.      2.0    CLG
C       10 Jul 2012   Fixed some small bugs in 2.0           2.02   RCZ
C       14 Jul 2012   Format of ".TXT" input files changed   2.03   RCZ
C                     so labels follow data; see FORMAT stmt
C                     label 1001
C       18 Jul 2012   Added statements to write file of      2.03   CLG
C                     as WC_IOPs.txt in SPECKD.f
C       24 Jul 2012   Changed structure of CANHTZ.TXT        2.04   RCZ
C                     and added photoperiod, computed by
C                     E05NM, and sent to GL via CANHTZ.TXT
C       26 Jul 2012   Changed FORMAT stmt for input of WQ    2.04   RCZ
C                     parameters from G12.4 to F12.0 to 
C                     eliminate need to enter data with 
C                     decimal points
C       28 Jun 2012   Streamlined and similified menu        2.05    RCZ
C                     options
C          1 Aug 12   New hierarchical menu created          2.06   RCZ
C                     and ability to change IOP constants added
C                     to OWQINIT
C         8 Aug 12    All driving data values are now        2.07  RCZ
C                     read from and written to data files.
C                     Structure for SiteData.txt changed to
C                     anticipate time-series calculations and
C                     is no longer compatible with older versions
C        19 Aug 12    Name of ".IRR" file now user-definable  2.07 RCZ
C        21 Aug 12    Writing ZTOTAL to .IRR file header      2.08 RCZ
C         7 Sep 12    Correctred OWQINIT ver info in ".irr"   2.08 RCZ
C                     header files
C         1 Nov 12    Fixed spacing of header information     2.09 RCZ
C                     when writing irradiance (.IRR) files
C        29 Jun 13    Modified .IRR file output for Ver 2.10  2.10 RCZ
C        12 Aug 13    Commented out section to change         2.11 CLG
C                     atmospheric parameters   
C        13 Nov 13    Installed trap for incorrect            2.11 RCZ
C                      entry of option number for the
C                      main menu
C        21 Dec 15    Adjusted menu to indicate Calendar      2.14 RCZ
C                       date, not Julian Date
C        12 Nov 16    Cleaned up .IRR filenaming steps        2.14 RCZ
C                     to prevent overwriting files
C
C     ********************************************************************************
      INCLUDE 'CommonE0.f'
	INCLUDE 'CommonKd.f'
	DIMENSION LAMBDA1(301)
	DIMENSION EDCANTOP(61),SPECK5NM(61),LAM5(61)
      CHARACTER ANS
      CHARACTER*16 IRRFILE,SITEFILE,IOPCNST,ATMCNST,WQCONC,NEWFIL
      CHARACTER*64 COPYRITE,VER,HEADER
      REAL NUCDOM,NUCHLA,NUTURB		!clg added 15 June 2012
      REAL NULAT,NULON,NUDAY,NUHR,NUZ  !RCZ added 28 july 2012
      DATA COPYRITE/'Copyright (c) 2012 by CL Gallegos & RC Zimmerman'/
      DATA VER/'Version 2.14'/SITEFILE/'SiteData.txt'/
      DATA IOPCNST/'IOPcnsts.txt'/ATMCNST/'Atmos.txt'/
C      DATA WQCONC/'WQCONC.TXT'/  !not necessary
      DATA CDOM/0.15/CHLA/6.2/TURB/5.0/   !DEFAULT M-1, MG/M3, NTU
      DATA DUMMY/0./ISITE/0/IIOPC/0/IATMC/0/
      IRRFILE='OWQ10NM.IRR'
C
C   ******  Begin section to load required arrays and constants  **********
C
C  Load IOP arrays, AW, BW, and ASTCHL
	OPEN(UNIT=1,FILE='IOPArrays.txt',STATUS='OLD')
c      write(*,3)
c3     format('reading iopArrays')
	READ(1,*) (LAMBDA,AW(I),BW(I),ASTCHL(I),I=1,61)
	CLOSE(1)
C   Load IOP constants
10	OPEN(UNIT=1,FILE=IOPCNST,STATUS='OLD')
c      write(*,11)
c11    format(' reading IOP constants')
	READ(1,1001) SG              !Spectral slope of CDOM absorption
	READ(1,1001) SIGATRB         !NAP absorb. cross section at 440 nm
	READ(1,1001) STRB            !Spectral slope of NAP absorption
	READ(1,1001) BLTRB           !Baseline NAP absorption
	READ(1,1001) SIGBTRB     !Scattering cross section of turbidity
	READ(1,1001) ETA             !Scattering spectral exponent
	READ(1,1001) BB2B            !Particulate backscattering ratio
	READ(1,1001) APHST675	 !Chl-specific absorption at 675 nm
      IF(ISITE.EQ.1)GO TO 30
	CLOSE(1)
C  Load canopy height data
      OPEN(UNIT=1,FILE='CANHTZ.TXT',STATUS='OLD')
c      write(*,12)
c12    format('reading canopy height data')
      READ(1,15)CANOPYHT	  !Had to re-instate read cg 6/15/12
15    FORMAT(F10.3)
      CLOSE(UNIT=1)
C  Load site location, date, time and depth & WATER QUALITY data      
20    OPEN(UNIT=1,FILE=SITEFILE,STATUS='OLD') !,ERR=120) - REMOVE ERROR TRAP FOR NOW
c      write(*,21)
c21    format('reading site file')
	READ(1,1001)XLAT             !LATITUDE (DEG.) POSITIVE N      
      READ(1,1001)XLON             !Longitude (deg.) negative W
	READ(1,22)JD,HR,TOTALZ,CDOM,CHLA,TURB   !Day,Time,Depth,CDOM,Chl a,Turbidity
22    FORMAT(I4,5F10.4)   
	CLOSE(1)
      IF(ISITE.EQ.1)GO TO 30
C  Load atmospheric parameters
25	OPEN(UNIT=1,FILE=ATMCNST,STATUS='OLD')
c      write(*,26)
c26    format('reading atmospheric constants')      
	READ(1,1001) AM               ! Air mass type, 1-10
	READ(1,1001) WM	   ! Windspeed averaged over prev 24 h, 1-10 m/s
      READ(1,1001) W                !Instantaneous windspeed, 1-20 m/s
      READ(1,1001) RH               !Relative humidity
      READ(1,1001) PRESS            !Atmospheric pressure, mb
      READ(1,1001) WV               !Precipitable water vapor, cm
      READ(1,1001) HA               !Aerosol scale height
      READ(1,1001) V                !Visibility, km, 5-25
      CLOSE(1)
 1000	FORMAT(I4)
 1001	FORMAT(F10.4)
C   *********  End of section to load required arrays and constants   ******
C
C         **** Main menu to change irradiance file parameters ***
C
30    ISITE=0
      NEWSITE=0
      NEWIOP=0
      NEWATM=0
      NEWCWQ=0
32    WRITE(*,35)SITEFILE,IOPCNST        !,ATMCNST removed 8/12/2013 clg
35    FORMAT(//'               Irradiance File Creation Parameters'//
     1    '   1) Location, date, depth & WQ concentrations'
     2    ' from ',A16/
     3    '   2) IOP parameters from 'A16///
C     4    '   3) Atmospheric parameters from 'A16///
     5    ' Enter MENU number to change parameters, <RET> to continue:')
      READ(*,40)MENU
40    FORMAT(I1)
      IF(MENU.EQ.0)GO TO 600
      IF(MENU.NE.1    .AND.    MENU.NE.2)GO TO 32 !ERROR TRAP FOR INCORRECT ENTRY
C
C     *** Now, act on the chosen menu option
C
      IF(MENU.EQ.1)THEN       !Get new location, date time and depth
          WRITE(*,41)
41        FORMAT('New location & WQ concentrations from a different',
     1        ' file (Y/N)?')
          READ(*,42)ANS
42        FORMAT(A1)          
          IF(ANS.EQ.'Y'   .OR.    ANS.EQ.'y')THEN
              WRITE(*,44)
44            FORMAT('Name of file:')
              READ(*,46)NEWFIL
46            FORMAT(A16)
              IF(NEWFIL.NE.' ')THEN
                  SITEFILE=NEWFIL
                  ISITE=1
                  GO TO 20
              ELSE
                  GO TO 30
              END IF
          ELSE
              WRITE(*,50)XLAT   !,ADVANCE = 'NO')
50            FORMAT(/' Latitude now ',F10.4,' deg, new value (decimal'
     1            ' deg, postitive N):')
              READ(*,60)NULAT
60            FORMAT(F10.0)
              IF(NULAT.NE.0.)THEN
                  XLAT=NULAT
                  NEWSITE=1
              END IF
              WRITE(*,70)XLON   !,ADVANCE='NO')
70            FORMAT(/' Longitude now ',F10.4,' deg, new value (decimal'
     1            ' deg, positive E):')
              READ(*,60)NULON
              IF(NULON.NE.0.)THEN
                  XLON=NULON
                  NEWSITE=1
              END IF
              WRITE(*,80)JD   !,ADVANCE='NO')
80            FORMAT(/' Calendar date now ',I4,
     1            ', new value (1 TO 365):')
              READ(*,90)NUJD
90            FORMAT(I4)
              IF(NUJD.NE.0)THEN
                  JD=NUJD
                  NEWSITE=1
              END IF
              WRITE(*,95)HR   !,ADVANCE='NO')
95            FORMAT(/' Time of day now',F10.2,
     1            ' h, new value (decimal hr):')
              READ(*,60)NUHR
              IF(NUHR.GT.0.)THEN
                  HR=NUHR
                  NEWSITE=1
              END IF
              WRITE(*,100)TOTALZ
100           FORMAT(/' Water depth now',F10.2,' m, new value (m):')
              READ(*,60)NUZ
              IF(NUZ.GT.0)THEN
                  TOTALZ=NUZ
                  NEWSITE=1
              END IF
              WRITE(*,110)CDOM   !,ADVANCE='NO')
110           FORMAT(/' CDOM abs now ',F10.4 ' per m, new value (m-1):')
              READ(*,120)NUCDOM
120           FORMAT(F10.0)
              IF(NUCDOM.GT.0.)THEN
                  CDOM=NUCDOM
                  NEWSITE=1
              END IF
              WRITE(*,130)CHLA   !,ADVANCE='NO')
130           FORMAT(/' [Chl a] now',F10.4,' ug/L, new value:')
              READ(*,120)NUCHLA
              IF(NUCHLA.GT.0.)THEN
                  CHLA=NUCHLA
                  NEWSITE=1
              END IF
140           WRITE(*,150)TURB   !,ADVANCE='NO')
150           FORMAT(/' Turbidity now',F10.4,' NTU, new value:')
              READ(*,120)NUTURB
              IF(NUTURB.LT.0)GO TO 140
              IF(NUTURB.GT.0)THEN
                  TURB=NUTURB
                  NEWSITE=1
              END IF
              IF(NEWSITE.EQ.1) THEN           !save new Site DATA
                  WRITE(*,160)
160               FORMAT(/'Save Site Data to file (Y/N)?')
                  READ(*,42)ANS
                  IF(ANS.EQ.    'y'   .OR.    ANS.EQ.'Y')THEN
                      WRITE(*,170)SITEFILE
170                   FORMAT('New file name for site data (CURRENTLY ',
     1                    A16,'):')
                      READ(*,46)NEWFIL
                      IF(NEWFIL.NE.' ')THEN
                          SITEFILE=NEWFIL
                      ELSE
                          WRITE(*,180)
180                       FORMAT('Overwrite existing file (Y/N)?')
                          READ(*,110)ANS
                          IF(ANS.EQ.'y'   .OR.    ANS.EQ.'Y')THEN
                              SITEFILE=NEWFIL
                          ELSE
                              GO TO 30
                          END IF
                      END IF
                      OPEN(UNIT=2,FILE=SITEFILE,STATUS='REPLACE')
                      WRITE(2,190)XLAT,XLON,JD,HR,TOTALZ,CDOM,CHLA,TURB
190                   FORMAT(F10.4,'  LAT'/F10.4,'  LON'/I4,5F10.4/
     1                    'JDAY      HOUR     DEPTH      CDOM',
     2                    '   Chl a      Turb')
                      WRITE(*,195)
195                   FORMAT(' Enter a descriptor line to save '
     1                    'with the data (64 char max):')
                      READ(*,200)HEADER
200                   FORMAT(A64)                      
                      WRITE(2,200)HEADER
                      CLOSE(UNIT=2)
                  END IF
              END IF
          END IF
          GO TO 30
      ELSEIF(MENU.EQ.2)THEN       !Get new IOP parameters
          WRITE(*,210)
210       FORMAT('New IOP parameters from existing file (Y/N)?')
          READ(*,42)ANS
          IF(ANS.EQ.'Y'   .OR.    ANS.EQ.'y')THEN
              WRITE(*,220)
220           FORMAT('Name of IOP params file:')
              READ(*,230)NEWFIL
230           FORMAT(A16)
              IF(NEWFIL.NE.' ')THEN
                  IOPCNST=NEWFIL
                  ISITE=1
                  GO TO 10
              ELSE
                  GO TO 30
              END IF
          ELSE
              WRITE(*,240)SG
240           FORMAT(/' Spectral slope of CDOM abs now',F10.4,
     1            ' per nm, new value:')
              READ(*,250)XVAL
250           FORMAT(F10.0)
              IF(XSVAL.GT.0)THEN
                  SG=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,260)SIGATRB
260           FORMAT(/' Abs x-sec for NAP now',F10.4,
     1            ' m^2/g, new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  SIGATRB=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,270)STRB
270           FORMAT(/' Spectral slope for NAP-a now',F10.4,' per nm,',
     1            'new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  STRB=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,280)BLTRB
280           FORMAT(/' Baseline for NAP-a now',F10.4,
     1            ' m^2 g^-1, new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  BLTRB=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,290)SIGBTRB
290           FORMAT(/' Turbidity scattering x-sec now',F10.4,
     1            ' m^2 g^-1, new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  SIGBTRB=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,300)ETA
300           FORMAT(/' Spectral exponent for particulate scattering ',
     1            'now',F10.4,' new value:')   
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  ETA=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,310)BB2B
310           FORMAT(/' bb/b ratio now',F10.4,' new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  BB2B=XVAL
                  NEWIOP=1
              END IF
              WRITE(*,320)APHST675
320           FORMAT(/' Chl a*(675 nm) now',F10.4,
     1            ' m^2 mg^-1, new value:')
              READ(*,250)XVAL
              IF(XVAL.GT.O)THEN
                  APHST675=XVAL
                  NEWIOP=1
              END IF
              IF(NEWIOP.EQ.1) THEN        !Save new IOP constants
                  WRITE(*,330)
330               FORMAT(/'Save IOP constants to file (Y/N)?')
                  READ(*,42)ANS
                  IF(ANS.EQ.    'y'   .OR.    ANS.EQ.'Y')THEN
                      WRITE(*,332)IOPCNST
332                   FORMAT('New file name for IOP constants '
     1                    '(CURRENTLY ',A16,'):')
                      READ(*,46)NEWFIL
                      IF(NEWFIL.NE.' ')THEN
                          IOPCNST=NEWFIL
                      ELSE
                          WRITE(*,180)
                          READ(*,42)ANS
                          IF(ANS.EQ.'y'   .OR.    ANS.EQ.'Y')THEN
                              IOPCNST=NEWFIL
                          ELSE
                              GO TO 30
                          END IF
                      END IF
                      OPEN(UNIT=2,FILE=IOPCNST,STATUS='REPLACE')
                      WRITE(2,340)SG,SIGATRB,STRB,BLTRB,SIGBTRB,ETA,
     1                    BB2B,APHST675,DUMMY
340                   FORMAT(F10.4,' SG'/F10.4,' SIGATRB'/F10.3,' STRB'/
     1                    F10.4,' BLTRB'/F10.3,' SIGBTRB'/F10.1,' ETA'/
     2                    F10.4,' BB2B'/F10.4,' APHST675'/,F10.4,
     3                    ' DUMMY')
                      WRITE(*,195)
                      READ(*,360)HEADER
360                   FORMAT(A64)
                      WRITE(2,360)HEADER
                      CLOSE(UNIT=2)
                  END IF
              END IF
C              GO TO 30
              END IF
C   Atmospheric parameters input section commented out 8/12/2013, CLG
C      ELSEIF(MENU.EQ.3)THEN       !Get new atmospheric parameters
C          WRITE(*,370)
C370       FORMAT('New atmos parameters from existing file (Y/N)?')
C          READ(*,42)ANS
C          IF(ANS.EQ.'Y'   .OR.    ANS.EQ.'y')THEN
C              WRITE(*,380)
C380           FORMAT('Name of file:')
C              READ(*,390)NEWFIL
C390           FORMAT(A16)
C              IF(NEWFIL.NE.' ')THEN
C                  IOPCNST=NEWFIL
C                  ISITE=1
C                  GO TO 10
C              ELSE
C                  GO TO 30
C              END IF
C          ELSE
C              WRITE(*,500)AM
C500           FORMAT(/' Air mass type now,'F10.1,
C     1            ', new value (1-10):')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  AM=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,510)WM
C510           FORMAT(/'24-h mean windspeed now',F10.1,
C     1            ' m/s, new value (1-10)')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  WM=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,520)W
C520           FORMAT(/' Instantaneous windspeed now',F10.1,
C     1            ' m/s, new value (1-20):')          
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  W=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,530)RH
C530           FORMAT(/' Relative humidity now',F10.1,
C     1            '%, new value (%):')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  RH=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,540)PRESS
C540           FORMAT(/' Atm press now',F10.1,' mb, new value (mb):')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  PRESS=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,550)WV
C550           FORMAT(/' Precipitable wate vapor now',F10.1,
C     1            ' cm, new value (cm):')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  WV=XVAL
C                  NEWATM=1
C              END IF
C              WRITE(*,560)V
C560           FORMAT(/' Visibility now',F10.1,
C     1            ' km, new value (km):')
C              READ(*,250)XVAL
C              IF(XVAL.GT.O)THEN
C                  V=XVAL
C                  NEWATM=1
C              END IF
C              IF(NEWATM.EQ.1) THEN        !Save new atmos params
C                  WRITE(*,570)
C570               FORMAT('Save Atmospheric params to file (Y/N?') 
C                  READ(*,42)ANS
C                  IF(ANS.EQ.    'y'   .OR.    ANS.EQ.'Y')THEN
C                      WRITE(*,575)ATMCNST
C575                   FORMAT('New file name for atmos params ',
C     1                    '(CURRENTLY ',A16,'):')
C                      READ(*,46)NEWFIL
C                      IF(NEWFIL.NE.' ')THEN
C                          ATMCNST=NEWFIL
C                      ELSE
C                          WRITE(*,180)
C                          READ(*,110)ANS
C                          IF(ANS.EQ.'y'   .OR.    ANS.EQ.'Y')THEN
C                              ATMCNST=NEWFIL
C                          ELSE
C                              GO TO 30
C                          END IF
C                      END IF
C                      OPEN(UNIT=2,FILE=ATMCNST,STATUS='REPLACE')
C                      WRITE(2,580)AM,WM,W,RH,PRESS,WV,HA,V,DUMMY
C580                   FORMAT(F10.1,' AM'/F10.1,' WM'/F10.1,' W'/
C     1                    F10.1,' RH'/F10.1,' PRESS'/F10.1,' WV'/
C     2                    F10.1,' HA'/F10.1,' V'/,F10.1,' DUMMY')
C                      WRITE(*,195)
C                      READ(*,195)HEADER
C                      WRITE(2,195)HEADER
C                      CLOSE(UNIT=2)
C                  END IF
C              END IF
C          END IF
      END IF
C      
C     ******  Section to read input and call incident irradiance routine  **** 	
C  Read input arrays for radiant flux at top of atmosphere, and
C  absorption spectra of atmospheric components.  This will be 
C  transparent to user
C      pause 'Opening E0InputArrays.txt'
600	OPEN(UNIT=1,FILE='E0InputArrays.txt',STATUS='OLD')
C      pause 'file E0InputArrays.txt opened OK'
	DO 620 I=1,301
C          pause 'getting ready to read first line'
          READ(1,610)LAMBDA1(I),H0(I),AOZ(I),AW_INC(I),AOX(I)
610       FORMAT(F8.0,4(G10.3))
C           pause 'line read'
C8010      FORMAT(1x,i4,5(1x,f8.3))
620   CONTINUE
C        WRITE(*,8010)I-1,LAMBDA1(I-1),H0(I-1),AOZ(I-1),AW_INC(I-1),
C     &AOX(I-1)
	CLOSE(1)
 	CALL E05NM(ED5NM,THETAZ)   !Returns subsurface Ed at 5 nm intervals
      OPEN(UNIT=2,FILE='CANHTZ.TXT',STATUS='REPLACE') !now write depth & DAYLENGTH to canhtz
      WRITE(2,630)CANOPYHT
630   FORMAT(G10.3,5X,'Canopy Height, m')
      WRITE(2,640)TOTALZ
640   FORMAT(G10.3,5X,'Water Depth, m')
      WRITE(2,650)DAYLEN
650   FORMAT(G10.3,5X,'Photoperiod, h'//)
      CLOSE(UNIT=2)
C   ********* End of call to incident irradiance spectrum  ***********
C
C  ***** Begin section to call routine to calculate spectral Kd  ******
C          Routine returns Kd spectrum at 5nm intervals, SPECK5NM
C
      CALL SPECKD(CDOM,CHLA,TURB,SPECK5NM)      
      CALL ATTEN(ED5NM,SPECK5NM,EDCANTOP)	
C   EDCANTOP is Ed at top of seagrass canopy, W m^-2 nm^-1, 5 nm 
C   intervals
      K=1
      DO 700 I=400,700,5
          LAM5(K)=I
C          WRITE(*,1550)K, LAM5(K)
C1550      FORMAT(' LAM5(',I3,') ='I5)         
          K=K+1
700   CONTINUE          
C     write 5 nm resolution file	
      OPEN(UNIT=2,FILE='CGIRR05.IRR',STATUS='REPLACE')
C      write(*,1700)
C1700  format('Unit 2 file cgirr05.irr opened')
      WRITE(2,800)VER,XLAT,XLON,JD,HR,ATMCNST,CDOM,CHLA,TURB,TOTALZ,
     1 ZTOC,DAYLEN
800   FORMAT(' Spectral Ed & Kd for GL input from OWQINIT ',A64/
     1 ' Latitude',F10.3,' Longitude',F10.3/
     2 ' Day',6X,I4,' Time',F10.3/
     3 ' Atmospheric parameters from ',A16/
     4 ' ag(440) (per m) =',F10.3,' Chl a (ug/L) = ',F10.3,
     5 ' TSM (mg/L or NTU) ='F10.3/
     6 ' Water depth (m) =',F10.3,
     7 ' Depth to top of canopy (m) =',F10.3/
     8 ' Photoperiod (h) =',F10.3/
     9 ' Blank line 9 for future use'/
     1 ' WAVEL Ed Kd: line 10 last variable header line')
      DO 810 I=1,61
          WRITE(2,805)LAM5(I),EDCANTOP(I),SPECK5NM(I)
805       FORMAT(I5,3F10.3)
810   CONTINUE    
      CLOSE(UNIT=2)
C       write(*,2025)
c2025  format('Unit 2 file cgirr05.irr closed')
C
C     Now write 10 nm resolution file
815   WRITE(*,820)IRRFILE
820   FORMAT('Name of new .IRR file (curently ',A16,'):')
      READ(*,46)NEWFIL
      IF(NEWFIL.NE.' ')THEN
          IRRFILE=NEWFIL
      ELSE
          WRITE(*,180)
          READ(*,42)ANS
          IF(ANS.NE.'y'   .OR.    ANS.NE.'Y')THEN
              GO TO 815
          ELSE
              IRRFILE=NEWFIL           
          END IF
      END IF
      OPEN(UNIT=2,FILE=IRRFILE,STATUS='REPLACE')
      WRITE(2,900)VER,XLAT,XLON,JD,HR,ATMCNST,CDOM,CHLA,TURB,TOTALZ,
     1 ZTOC,DAYLEN
900   FORMAT(' Spectral Ed & Kd for GL input from OWQINIT 'A64/
     1 ' Latitude',F10.3,' Longitude'F10.3/
     2 ' Day',6X,I4,' Time',F10.3/
     3 ' Atmospheric parameters from ',A16/
     4 ' ag(440) (per m) =',F10.3,' Chl a (ug/L) = ',F10.3,
     5 ' TSM (mg/L or NTU) ='F10.3/
     6 ' Water Depth (m) =',F10.3/
     7 ' Depth to top of canopy (m) =',F10.3/
     8 ' Photoperiod (h) =',F10.3/
     9 ' Blank line 9 for future use'/
     1 ' WAVEL Ed Kd: line 10 last variable header line')
      DO 910 I=1,61,2
          WRITE(2,805)LAM5(I),EDCANTOP(I),SPECK5NM(I)
910   CONTINUE
      CLOSE (UNIT=2)
C      WRITE(*,2045)      !debug
C2045  FORMAT('UNIT 2 FILE CGIRR10.IRR CLOSED')   !debug
 	RETURN
      END
