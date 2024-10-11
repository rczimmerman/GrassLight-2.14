	SUBROUTINE SPECKD(CDOM,CHLA,TURB,SPECK5NM)
	DIMENSION SPECK5NM(61)
	INCLUDE 'CommonKd.f'
C......../........./........./........./........./........./........./.S	
C   Calculate absorption components, scattering, and attenuation
	THETADEG=THETAZ*180./3.14159265  !In-air solar angle, (deg.) as required by Lee et al. (2005)	
c      write(*,10)sg,bltrb,bb2b
c10    format('sg= 'e16.4,' bltrb= ',e16.4,' bb2b= ',e16.4/)
      OPEN(UNIT=2,FILE="WC_IOPs.txt",STATUS="UNKNOWN")
      WRITE(2,1002)
      DO 100 I=1,61
	    WL=395.+5.*FLOAT(I)
	    AG=CDOM*EXP(-SG*(WL-440.))  ! CDOM  absorption
	    ACHL=CHLA*APHST675*ASTCHL(I)  ! Chlorophyll absorption
	    ANAP=TURB*(BLTRB+SIGATRB*EXP(-STRB*(WL-440.)))  ! NAP absorption
	    AT=AW(I)+AG+ACHL+ANAP   !Total absorption
	    BP=TURB*SIGBTRB*(555./WL)**ETA  !Particulate scattering
	    BB=0.5*BW(I)+BB2B*BP    !Total backscattering
	    WRITE(2,1001) WL,AG,ACHL,ANAP,BP,BB
C   Next statement implements Eq. 11 of Lee et al. (2005)	
	    SPECK5NM(I)=(1.+0.005*THETADEG)*AT
     &    +4.18*(1.-0.52*EXP(-10.8*AT))*BB
c          write(*,50)ag,achl,anap
c50          format('ag= 'e16.4,'achl= ',e16.4,'anap= ',e16.4)      
 100  CONTINUE        
       CLOSE(2)
 1001 FORMAT(F6.0,5F12.4)
 1002 FORMAT("WL	AG	ACHL	ANAP	BP	BB")
 	RETURN
 	END	
	
