	SUBROUTINE EDCALC(NLAYER,OUTINC,LATHIK,DENSTY,BETA,
     1    MUBARD,IRRFILE,SEDREF,LEAFIOP,AEPI,HR)
C	*************************************************************************
C                 Subroutine EDCALC Ver. 2.14
C
C		Calculates the propagation of downwelling irradiance
C		through a seagrass canopy of biomass distribution determined
C		by BIODIS
C		
C		Usage: CALL EDCALC(NLAYER,OUTINC,LATHIK,DENSTY,BETA,IRRFILE,SEDREF,LEAFIOP,AEPI)
C		
C		Definitions of Passed Variables:
C		---------------------------------
C		NLAYER - Number of canopy layers determined by subroutine BIODIS
C		LATHIK - Thickness of each canopy layer (mm) deterined by 
C				subroutine BIODIS.
C		LTHIKM - Thickmess of each canopy layer (m) used here to apply water
C				column attenuation coefficients (K's).
C		DENSTY - Shoot denisty (shoots/m2).
C		BETA - Canopy Bending Angle (radians) determined by subroutine BIODIS.
c		MUBARD - Average cosine of downwelling irradiance at top of canopy
C		PIZ - Spectral array of photosynthetic coefficient for absorbed photons.
C 		IRRFIL - Input file containing downwelling spectral irradiance
C				(Ed, W/m2/nm) at canopy surface and water column attenuation
C				coefficients (Kd and Ku).  
C		SEDREF - Input file containing sediment spectral reflectance data
C				(proportional values).  Starting wavelength value and
C				wavelength intervals MUST match those of IRRFIL.
C		LEAFIOP - Input file containing leaf spectral absorption (ln) and 
C				reflectance (relative) data.  Starting wavelength and 
C				wavelength intervals MUST match those of IRRFIL and SEDREF
C				but total file length may be longer to allow correction for
C				non-photosynthetic absorption using long wavelengths (>700 nm).
c		OUTINC - Output increment (in cm) for printing results of Ed and Eu calculations.
C         AEPI - Epiphyte Absorptance, set in GL
C         HR - Time of day from .IRR file.  If HR =  12 (noon), then calculate
C                 daily production estimates using sin model
C
C		COMMOM variables are array variables shared with other subroutines and the 
C				main program
C		Definitions of common variables:
C		-------------------------------
C		BIOMAS(I) - Relative biomass present in canopy layer (I).  All
C				biomass elemenct should sum to 1.  Calculated in
C				subroutine BIODIS
C		LAI(I) - Leaf area index (m2 leaf/m2 area) in each canopy layer (I).
C				Determined in subroutine BIODIS
C		LAP(I) - Leaf area projected toward the zenith.  Calculated as a
C				function of bending angle in subroutine BIODIS 
C		WAVLED(J) - Wavelength of spectral element j, determined in 
C					subroutine EDCALC from irradiance inupt file
C		PIZ(I) - Spectral production coefficient Pi(z) calculated in
C				subroutine EDCALC
C		EO(I) - Scalar quantum irradiance (Ed + Eu) in canopy layer I calculated in
C				subroutine EDCALC
C         EK(I) - Spectrally integrated value of Ek for each canopy layer; used in
C                 subroutine PRODUC to integrate daily P as f(Ek/Enoon)
C
C		INTERNAL variables and arrays are used only in Edcalc
C		Definitions of important internal variables
C		-------------------------------------------
C		LABSL(I) - Leaf absorption spectrum, ln units, per leaf
C		LREFL(I) - Leaf reclectance spectrum, relative units, per leaf
C		aSTAR - Leaf absorption spectrum corrected for non-photosynthetic
C				absorption, ln units
C		LABSB(I) - Leaf photosynthetic absorptance spectrum, relative units
C
C
C         Files required by EDCLC_210:
C         ----------------------------------
C             IRRFILE - data file containing spectral irradiance and diffuse attenuatio coefficient.
C                 Name of the specific file is passed to this subroutine from GL.  
C             SEDREF - data file containing the reflectance spectrum of the sediment underlying the
C                 plant canopy.  Name of the specific file is passed to this subroutine from GL.
C             LEAFIOP - data file containing the absorption and reflectance spectra for the leaves.
C                 Name of the specific file is passed to this subroutine from GL.
C
C
C
C         Files created by EDCLC_210:
C         ----------------------------------
C             OUTFIL (Unit = 3) - opened in GL, this subroutine provides Ed and Eu spectra at
C                 specified depth intervals through the plant canopy
C
C
C			Created: 9 Mar 2000                  by:  RCZ
C			Version 1.0 Finalized 30 Dec 2000    by:  RCZ
C
C			Date				Modification		                Ver     By
C		---------------  ----------------------------------------  ------ ------	
C          5 Sep 01    Eu calculation modified to include             1.1     RCZ
C                           effects of lambertian (cos 45 deg)
C						submarine light field
C          6 Sep 01	Irradiances now output at user-specified        1.15    RCZ
C						intervals (whole cms)
C          8 Sep 01	MUBARD added to Ed attenuation and 		        1.2     RCZ
C						reflectance terms
C         11 Sep 01	EDREF replaced by EU(I,J) in both		        1.3     RCZ
C						ED and EU attenuation calcs
C         11 Sep 01	Productivity calcs revised to allow		        1.3     RCZ
C						variable wavelength resolution
C         12 Sep 01	MUBARD now passed from GL as user-		        1.4     RCZ
C						modifiable variable 
C         29 Oct 01	MUBAR asymptotes to 0.5 not cos(45)		        1.5     RCZ
C         17 Dec 01	KD(z) now depends on MUBAR 			            1.6     RCZ
C         17 Dec 01	KU(z) now calculated from Kd(z) & MUBAR	        1.6     RCZ
C         31 Jan 02	Added double quotes (" ") to text written       1.61    RCZ
C					to output files to improve parsing by
C						Excel, etc.
C         28 Jun 07	Corrected calculation of MUBARD to 		        1.62    RCZ
C					increment with each layer, not wavelength	
C         14 Jul 07	Corrected calculation of EU(I,J) by		        1.63    RCZ
C					inserting missing parentheses in exponential
C					decay function
C         15 May 12    Revised to accomodate Intel                    2.00    RCZ
C                     FORTRAN 95 Compiler Requirements for
C                     Win7 x64 OS
C         14 Jul 12    Added /CONVRG/ to common block and             2.03    RCZ
C                        DONE to integer variable list
C         22 Aug 12  File name dimensions now CHARACTER*16            2.08    RCZ   
C                       Fixed formatting bug reading IRRFILE
C         22 Oct 12  Eu absorbed now affected by epiphyte load        2.09    RCZ
C         13 Jun 13  ALPHA dimensions reduced to (100), and           2.10    RCZ
C                       made available to other subroutines
C                       in common block for future use
C                       better determination of integrated daily
C                       production
C         3 Nov 16    Reads HR variable and passes it to subroutine   2.14    RCZ 
C                     PRODUC to calculate daily integral P when HR = 12
C
C	**************************************************************************
C
	IMPLICIT REAL*8 (A-Z)
	INTEGER I,J,K,L,M,N,NLAYER,NUMINT,NWAVED,NWAVRB,NWAVLF,
     1	NCOL,MROW,OUTINC,INTHT,ITER,DONE,FINL
	CHARACTER ANS,RUN
	CHARACTER*16 IRRFILE,SEDREF,LEAFIOP
	CHARACTER*64 COPYRITE,VER
	COMMON /MASS/ASYMP,INFLEC,SHAPE,BETA0,BIOMAS(100)
     1	/AREA/LAI(100)
     2    /PROJEC/LAP(100)
     3	/WAVE/WAVLED(100)
     4    /PROD/PIZ(100),ALPHA(100)
     5    /IRR/EO(100)
     6    /CONVRG/ITER,DONE,FINL,RUN
	DIMENSION WAVLRB(100),KD(100),KU(100),RB(100),WAVLLF(100),
     1	LABSL(100),LREFL(100),LABSM(100,100),LREFM(100,100),
     2	ED(101,100),EU(101,100),ASTAR(100),LABSB(100)
	DATA COPYRITE/'Copyright (c) 2013 by Richard C. Zimmerman'/
	DATA VER/'2.14'/
C
C				*** INITIAL PARAMETERIZATION ***
C
c	WRITE(*,10)VER
c10	FORMAT(' Subroutine EdCalc Ver. ',A64)
	LTHIKM=LATHIK/1000.
	I=0
C
C				*** OPEN IRRADIANCE DATA FILE ***
C
	OPEN(2,FILE=IRRFILE)
C
C			*** READ 10 HEADER RECORD LINES ***
C
C*****THIS SECTION JUST BURNS THE HEADER LINES
C	DO 40 I=1,10
C		READ(2,30)HDR
C30		FORMAT(A64)			
C40    CONTINUE
C******END SECTION TO BURN HEADER LINES
C
C             *** READ 3 HEADER LINES TO GET TIME OF DAY ***
C
      READ(2,30)HDR
30    FORMAT(A64)
      READ(2,30)HDR
      READ(2,40)HR
40    FORMAT(22X,F10.3)
C
C             *** NOW BURN THE REMAINING HEADER LINES ***
      DO 45 I=1,7
		READ(2,30)HDR
45    CONTINUE
C     
C
C			***LOOP THROUGH IRRADIANCE DATA FILE***
C
C	READ(2,*)NCOL,MROW  !No longer needed, removed these data from all .IRR files in Ver 2.10 and higher
	NWAVED=1
50    READ(2,60,END=90)WAVLED(NWAVED),ED(NLAYER+1,NWAVED),
     1	KD(NWAVED)
60    FORMAT(F10.0,2F10.3)
c         write(*,60)WAVLED(NWAVED),ED(NLAYER+1,NWAVED),KD(NWAVED)
c          pause
		NWAVED=NWAVED+1
	GO TO 50
90	CLOSE(2)
	NWAVED=NWAVED-1
	WAVDIF=ABS(WAVLED(1)-WAVLED(2))
C
C			*** GET BOTTOM REFLECTANCE FILE ***
C
	OPEN(2,FILE=SEDREF)
C
C			*** READ 10 HEADER RECORDS ***
C
	DO 140 I=1,10
		READ(2,30,ERR=143)HDR
140	CONTINUE
	GO TO 150
C			*** MISSING SEDREF ERROR TRAP ***
143	WRITE(*,145)JFILE
145	FORMAT(' ......can not find file ',A15)
	GO TO 1000
C
C			***LOOP THROUGH RB DATA FILE***
C
150	NWAVRB=1
155		READ(2,*,END=190)WAVLRB(NWAVRB),RB(NWAVRB)
C		IF(WAVLRB(NWAVRB).EQ.-1)GO TO 190   !no longer needed
		NWAVRB=NWAVRB+1
		GO TO 155
190	CLOSE(2)
	NWAVRB=NWAVRB-1
C
C			*** GET LEAF IOPs FILE***
C
	OPEN(2,FILE=LEAFIOP)
C
C			*** READ 10 HEADER RECORDS ***
C
	DO 240 I=1,10
		READ(2,30,ERR=243)HDR
240	CONTINUE
	GO TO 250
C			*** MISSING LEAFIOP ERROR TRAP ***
243	WRITE(*,245)LEAFIOP
245	FORMAT(' ......can not find file ',A15)
	GO TO 1000
C
C			***LOOP THROUGH LEAF IOP DATA FILE***
C
250	NWAVLF=1
255	READ(2,*,END=290)WAVLLF(NWAVLF),LABSL(NWAVLF),LREFL(NWAVLF)
		IF(LREFL(NWAVLF).GT.1.0)LREFL(NWAVLF)=LREFL(NWAVLF)/100.
		NWAVLF=NWAVLF+1
	GO TO 255
270	WRITE(*,280)KFILE
280	FORMAT(' .....can not find ',A12)
	GO TO 1000
290	CLOSE(2)
	NWAVLF=NWAVLF-1
C
C			*** COMPUTE LEAF aSTAR AND PHOTOSYNTHETIC ABSORPTANCE ***
C			*** BY SUBTRACTING LABSL(NWAVLF) FROM EACH LABSL(NWAVED) ***
	DO 295 J=1,NWAVED
	    ASTAR(J)=LABSL(J)-LABSL(NWAVLF)
	    LABSB(J)=1.-EXP(-ASTAR(J))
c	write(*,298)wavllf(j),labsl(j),aSTAR(j),labsb(j)
c298	format(' Wavel=',f5.0,' Leaf a=',f7.3,' a*=',F7.3,' A=',f7.3)
295	CONTINUE
C
C		*** SCALE LEAF TOTAL ABSORPTION AND REFLECTANCE COEFFICIENTS ***
C			*** BY PROJECTED LEAF AREA WITHIN EACH LAYER***
C
	DO 300 I=1,NLAYER
		DO 300 J=1,NWAVED
			LABSM(I,J)=LABSL(J)*LAP(I)
			LREFM(I,J)=LREFL(J)*LAP(I)
300	CONTINUE
C
C			*** CHECK WAVELENGTH SPACING ON IRRADIANCE, ***
C				    REFLECTANCE AND LEAF IOPS
C
	IF(NWAVLF.GE.NWAVRB	.AND.	NWAVLF.GE.NWAVED)GO TO 320
310		STOP 'Error: Input file wavelengths do not match '
320	IF(WAVLLF(1).NE.WAVLRB(1).AND.WAVLLF(1).NE.WAVLED(1))GO TO 310
C
C			*** future version should interpolate wavelength data ***
C			    to provide a common wavelength array for all calcs
C
C			*** NOW ATTENUATE Ed DOWN THROUGH THE CANOPY ***
C				AND DETERMINE THE PHOTOSYNTHETIC PARAMETER PI(z)
C				THE CONSTANT 0.008351 CONVERTS ENERGY (W/m2/nm)
C				TO QUANTA (uE/m2/s) PER NM
C				THE CONSTANT 0.5 is the avg cosine for Lambertian scattering
C				REFLECTED ED NOW ACCUMULATED IN EU ARRAY   11 SEP 01
400	DO 411 I=NLAYER,1,-1
		L=I+1
		PIZ(I)=0.
		DO 410 J=1,NWAVED
              IF(LAP(I).LT.LAI(I))THEN
                  ALPHA(J)=0.1*LABSB(J)*SIN(BETA)   !Leaves are not horizontal
              ELSE
                  ALPHA(J)=0.1*LABSB(J) !Leaves are horizontal
              END IF
			EU(I,J)=ED(L,J)*LREFM(I,J)/MUBARD
			EXPONE=EXP(-(KD(J)*LTHIKM+LABSM(I,J)/MUBARD))
			ED(I,J)=(ED(L,J)-EU(I,J))*EXPONE
c              write(*,405)i,j,ed(i,j),kd(j)
c405           format('layer ='i5,' wavel='i5,'Ed =',f10.3,' Kd='f10.3)
			PIZ(I)=PIZ(I)+ALPHA(J)*((ED(L,J)-EU(I,J))-
     1            (AEPI*(ED(L,J)-EU(I,J))))*          !SUBTRACT EPIPHYTE ABSORBTION
     2            WAVLED(J)*0.008351*WAVDIF           !CONVERT ENERGY TO QUANTA          
              EO(I)=EO(I)+WAVLED(J)*ED(I,J)*0.008351*WAVDIF
c              write(*,406)EO(i),wavled(j),ed(i,j),wavdif
c406           format('e0(i)=',4f10.3)              
			OLDMU=MUBARD
		    MUBARD=MUBARD-(MUBARD-0.5)*LAP(I)
		    KD(J)=KD(J)*OLDMU/MUBARD
c          write(*,4105)i,mubard,oldmu,lap(i),kd(j)
c4105      format('layer='i4,'Mu-bar='f10.3'old Mu='f10.3,'lap(i)='
c     1    f10.3,' Kd(j)'f10.3)          
410       CONTINUE
c          write(*,4050)piz(i),ed(l,j),eu(i,j),aepi
c4050      format(' piz =',f10.3,' ed=',f10.3,' eu=',f10.3,' aepi='f10.3) 
c          pause
c	write(*,4110)i,eo(i)
c4110	format(' Layer ',i3,' PAR =',f10.0,' umol/m2/s')
c	if(mod(i,10).eq.0)pause
411	CONTINUE
c	pause
C
C			*** WRITE THE ED ARRAY TO A TEXT FILE ***
C
	WRITE(3,418)DENSTY
418	FORMAT(' "Shoot density (shoots/m2) = "',F6.0/
     1	' "Ed (W/m2/nm) as a function of canopy height:"')
	WRITE(3,420)(WAVLED(J),J=1,NWAVED)
420	FORMAT('"Height (m)" "Wavel (nm)"',31F15.0)
	HT=LTHIKM*FLOAT(NLAYER+1)
	WRITE(3,430)HT,(ED(NLAYER+1,J),J=1,NWAVED)
	DO 440 I=NLAYER,1,-1
		HT=LTHIKM*FLOAT(I)
		INTHT=INT(HT*100.)
		M=MOD(INTHT,OUTINC)
c	write(*,421)i,ht,intht,outinc,m
c421	format(' checking layer ',i5,' height =',f10.3,' intht =',
c     1	i10,' outinc=',i5,' m =',i5)
c		pause
		IF(M.EQ.0)WRITE(3,430)HT,(ED(I,J),J=1,NWAVED)
430		FORMAT(F10.2,31F15.4)
440	CONTINUE
	HT=0.
	WRITE(3,430)HT,(ED(1,J),J=1,NWAVED)
C
C			*** NOW REFLECT Ed OFF THE BOTTOM TO MAKE Eu ***
C
	DO 500 J=1,NWAVED
		EU(1,J)=ED(1,J)*RB(J)
		KU(J)=KD(J)*MUBARD/0.5
500	CONTINUE
C
C			*** NOW CALCULATE Eu UPWARD THROUGH THE CANOPY ***
C				AND DETERMINE THE AMOUNT OF Eu ABSORBED
C					IN EACH LAYER ON THE WAY UP
C				INITIAL VALUE OF Eu DETERMINED BY Ed REFLECTED BY EACH
C					   LAYER ON THE WAY DOWN 
C				the constant 0.5 = Mubar for lambertian scattering
	DO 512 I=2,NLAYER
		DO 510 J=1,NWAVED
              IF(LAP(I).LT.LAI(I))THEN
                  ALPHA(J)=0.1*LABSB(J)*SIN(BETA)   !Leaves are not horizontal
              ELSE
                  ALPHA(J)=0.1*LABSB(J) !Leaves are horizontal
              END IF
			EUREF=EU(I-1,J)*LREFM(I,J)/0.5
			EU(I,J)=(EU(I,J)+(EU(I-1,J)-EUREF))
     1			*EXP(-(KU(J)*LTHIKM+LABSM(I,J)/0.5))
		    PIZ(I)=PIZ(I)+ALPHA(J)*(EU(I,J)-AEPI*EU(I,J))
     2             *0.008351*WAVLED(J)*WAVDIF
			EO(I)=EO(I)+EU(I,J)*WAVLED(J)*0.008351*WAVDIF
510		CONTINUE
512	CONTINUE
C
C			*** WRITE Eu DATA TO A TEXT FILE ***
C
	WRITE(3,515)
515	FORMAT(///' "Eu (W/m2/nm) as a function of depth:"')
	WRITE(3,420)(WAVLED(I),I=1,NWAVED)
	HT=LTHIKM*FLOAT(NLAYER)
	WRITE(3,430)HT,(EU(NLAYER,J),J=1,NWAVED)
		DO 520 I=NLAYER,1,-1
		HT=LTHIKM*FLOAT(I)
		INTHT=INT(HT*100.)
		M=MOD(INTHT,OUTINC)
		IF(M.EQ.0)WRITE(3,430)HT,(EU(I,J),J=1,NWAVED)
520	CONTINUE
	HT=0.
	WRITE(3,430)HT,(EU(1,J),J=1,NWAVED)
C      canr555=((eu(nlayer,16)+eu(nlayer,17))/2.)/  !canopy reflectance @ 555
C     1 ((ed(nlayer,16)+ed(nlayer,17))/2.)              !for Victoria's paper
C      write(*,600)canr555
C600   format('Canopy Reflectance (555)='f10.3//)
1000	RETURN
	END
