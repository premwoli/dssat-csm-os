C=======================================================================
C  FOR_GRAZE, Subroutine, P. Woli
C-----------------------------------------------------------------------
C  Governing Model Equations
C  Integrates all grazing related variables
C-----------------------------------------------------------------------
C  REVISION       HISTORY
C  11/07/2016 PW  Written
C-----------------------------------------------------------------------
C  Called by:  forage
C  Calls:      CATTLEGROWTH, GRAZINGDATA, GETLUN, FIND, IGNORE, ERROR
C=======================================================================
!
      SUBROUTINE FOR_GRAZE(YDOY,RHUM,WINDSP,Rain,Tmin,Tmean,FILECC,
     &                RHOL,RHOS,PCNL,PCNST,SLA,           !Input
     &                WTLF,STMWT,TOPWT,TOTWT,WCRLF,WCRST, !Input/Output
     &                WTNLF,WTNST,WNRLF,WNRST,WTNCAN,     !Input/Output
     &                AREALF,XLAI,XHLAI,VSTAGE,CANHT,     !Input/Output/CHTCM
     &                Gtot,Gtotn,Gpctlf,Gpctn)
      
      USE MODULEDEFS

      IMPLICIT NONE

      INTEGER :: ISECT,ERR
      INTEGER :: LUNCRP,LNUM,FOUND
      INTEGER :: YDOY,J,IYR,i,DOP,DAY,age
      !INTEGER :: GYB=-99,GYE=-99,GDB=-99, GDE=-99
      INTEGER :: FLAGRW=0,FLAGVF=0,FLAGVP=0
      INTEGER :: STB,STE,LCB,LCE,OLB,OLE,PLB,PLE,RSB,RSE,RDB,RDE,GDB,GDE
      REAL    :: RHUM,WINDSP,Rain,Tmin,Tmean,Tc,Tp,RHc,RHp,WSc,WSp
      REAL    :: RHOL,RHOS,PCNL,PCNST,SLA,LSR
      REAL    :: AREALF,XLAI,XHLAI,VSTAGE,vstagp,CANHT,CHTCM
      REAL    :: WTLF,STMWT,TOPWT,TOTWT,WCRLF,WCRST,RTWT,GP,WTLFa,STMWTa
      REAL    :: WTNLF,WTNST,WNRLF,WNRST,WTNCAN,WTLFb,STMWTb
      REAL    :: TOPb,TOPa,TOPWTb=0,STBLE=0,RSPLF=-1.0E25
      REAL    :: GLEAF,GSTEM,GVSTG,AREAH,PROLFF,PROSTF,pliglf,pligst
      REAL    :: Gtot,Gtotn,Gpctlf,Gpctn,Glfn,Gstn
      REAL    :: Gcrlf,Gcrst,Gplig,Gpcho
      REAL    :: SW =-1.0E25,SWF =-1.0E25,ADG =0,SRW,EQSBW,DFA,CFA
      REAL    :: AFA,AFB,AFM,HAF,HMF,AFF,AFT,FAF,ABF,PTAF,PETI,EBF
      REAL    :: INTAKE=0,DMIP,PDMI,DMI=0.001,DMIm=0.001,GRAZED=0,HayIn
      REAL    :: CP,NDF,ASH,LIG,TDN,NDIP,IADIP,NDFn,dNFC,dCP,dNDF,tTDN
      REAL    :: aTDN,dTDN,wTDN,DMIR,Phi
      REAL    :: NEmI,NEmr,NEmrB,NEmrPA,NEmrCS,NEmrHS,NEmrCHS,NEmrAGE
      REAL    :: DEd,MEd,NEd,NEtr,NEr,NEgr,NEmd,NEgd,NEgI,NEAG
      REAL    :: TdEE,WdEE,ADIP,TdNDF,WdNDF,kdcp,discount
      CHARACTER*6  SECTION,ERRKEY,SEASON,GROUP,HAY
      CHARACTER*30 FILECC
      CHARACTER*80 MOW80
      CHARACTER (LEN=300) :: FMT=' '
      !
      REAL,    PARAMETER :: FBGL   = 0.80    !Fraction of biomass grazed as leaves USE ROTH SEASONAL VALUES           **
      REAL,    PARAMETER :: MVS    = 11.0    !leaf no./stem stubble VERY SENSITIVE
      REAL,    PARAMETER :: DEdBGH = 2.48    !DE,bermudagrass hay
      REAL,    PARAMETER :: SR     = 7.0     !stockers/ha (11=high, mid=8.5, 7=low) for 2008  ** SR=f(FAL)as an objective for the subsequent paper  
      REAL,    PARAMETER :: GTH    = 0.0     !Min.biomass to be retained (kg/ha) (in LEAVES?) **
      REAL,    PARAMETER :: SWI    = 250.0   !initial stocker wt. (kg/stocker)                **
      REAL,    PARAMETER :: EE     = 2.4     !ether extract (%) assumed b/o literature data
      REAL,    PARAMETER :: PAF    = 1.0     !processing adjustment factor (-); grasses = 1
      INTEGER, PARAMETER :: FSTB   = 135     !fall stocker begin (doy):       MAY 15
      INTEGER, PARAMETER :: FSTE   = 300     !fall stocker end (doy):         Oct 27
      INTEGER, PARAMETER :: FLCB   = 258     !fall lac.cow begin:             Sep 15
      INTEGER, PARAMETER :: FLCE   = 152     !fall lac.cow end:               Jun  1
      INTEGER, PARAMETER :: FOLB   = 258     !fall open lac.cow begin:        Sep 15
      INTEGER, PARAMETER :: FOLE   = 015     !fall open lac.cow end:          Jan 15
      INTEGER, PARAMETER :: FPLB   = 016     !fall preg.lac.cow begin:        Jan 16
      INTEGER, PARAMETER :: FPLE   = 152     !fall preg.lac.cow end:          Jun  1
      INTEGER, PARAMETER :: FRSB   = 153     !fall rep.heif.w/steers begin:   Jun 10
      INTEGER, PARAMETER :: FRSE   = 288     !fall rep.heif.w/steers end:     Oct 15
      INTEGER, PARAMETER :: FRDB   = 289     !fall rep.heif.w/o steers begin: Oct 16
      INTEGER, PARAMETER :: FRDE   = 001     !fall rep.heif.w/o steers end:   Jan  1
      INTEGER, PARAMETER :: WSTB   = 283     !wint stocker begin (doy):       Oct 10
      INTEGER, PARAMETER :: WSTE   = 135     !wint stocker end (doy):         May 15
      INTEGER, PARAMETER :: WLCB   = 015     !wint lac.cow begin:             Jan 15
      INTEGER, PARAMETER :: WLCE   = 274     !wint lac.cow end:               Oct  1
      INTEGER, PARAMETER :: WOLB   = 015     !wint open lac.cow begin:        Jan 15
      INTEGER, PARAMETER :: WOLE   = 105     !wint open lac.cow end:          Apr 15
      INTEGER, PARAMETER :: WPLB   = 106     !wint preg.lac.cow begin:        Apr 16
      INTEGER, PARAMETER :: WPLE   = 274     !wint preg.lac.cow end:          Oct  1
      INTEGER, PARAMETER :: WRSB   = 283     !wint rep.heif.w/steers begin:   Oct 10
      INTEGER, PARAMETER :: WRSE   = 091     !wint rep.heif.w/steers end:     Apr  1
      INTEGER, PARAMETER :: WRDB   = 092     !wint rep.heif.w/o steers begin: Apr  2
      INTEGER, PARAMETER :: WRDE   = 001     !wint rep.heif.w/o steers end:   Jan  1
      !
      SAVE FLAGRW, FLAGVF
      !SAVE GYB,GDB,GYE,GDE
      SAVE SW,SWF,DOP,DAY
      !SAVE GTH,SWI,SW,SWF,SR,DOP
      !
      J      = MOD(YDOY,1000)
      IYR    = (YDOY-J)/1000
      !
      HAY    = 'NO'                           !Feeding hay: 'YES/NO'?
      SEASON = 'FALL'                         !Calving season: FALL or WINT?
      GROUP  = 'ST'                           !Animal group: ST=stocker,LC=lac.cow,
                                              !OL=open lac.cow,PL=preg.lac.cow,
                                              !RS=rep.heif.w/steers,RD=rep.heif.w/o steers
      IF (SEASON .EQ. 'FALL') THEN
          STB = FSTB
          STE = FSTE
          LCB = FLCB
          LCE = FLCE
          OLB = FOLB
          OLE = FOLE
          PLB = FPLB
          PLE = FPLE
          RSB = FRSB
          RSE = FRSE
          RDB = FRDB
          RDE = FRDE
      ELSE IF (SEASON .EQ. 'WINT') THEN
          STB = WSTB
          STE = WSTE
          LCB = WLCB
          LCE = WLCE
          OLB = WOLB
          OLE = WOLE
          PLB = WPLB
          PLE = WPLE
          RSB = WRSB
          RSE = WRSE
          RDB = WRDB
          RDE = WRDE
      ENDIF
      IF      (GROUP .EQ. 'ST') THEN
          GDB = STB
          GDE = STE
      ELSE IF (GROUP .EQ. 'LC') THEN
          GDB = LCB
          GDE = LCE
      ELSE IF (GROUP .EQ. 'OL') THEN
          GDB = OLB
          GDE = OLE
      ELSE IF (GROUP .EQ. 'PL') THEN
          GDB = PLB
          GDE = PLE
      ELSE IF (GROUP .EQ. 'RS') THEN
          GDB = RSB
          GDE = RSE
      ELSE IF (GROUP .EQ. 'RD') THEN
          GDB = RDB
          GDE = RDE
      ENDIF
      IF (GDE .LT. GDB) GDE = GDE + 365
      !    
      !Open GRAZE.RPT to report output and BG_GRAZCHK.DAT for checking
      IF (FLAGRW .EQ. 0) THEN                                 !Read input, write output
          OPEN (56,FILE='C:\DSSAT48\BG_Grazing\BG_GRAZE.DAT') 
          OPEN (58,FILE='C:\DSSAT48\BG_Grazing\BG_GRAZCHK.DAT')
          FLAGRW = 1
      ENDIF
      !
      !Calculate grazed biomass and stocker weight
      TOPWT = 10.0*TOPWT                                      !g/m2 to kg/ha
      WTLF  = 10.0*WTLF                                       !g/m2 to kg/ha
      DAY = DAY + 1
      IF ((DAY .GE. GDB) .AND. (DAY .LE. GDE))  THEN	        !GRAZING PERIOD
          DOP = DOP + 1
          IF (DOP .EQ. 1) THEN                                !First day of season 
              SW = SWI                      	  
          ENDIF
          NDF = (-0.00095*J**2 + 0.4713*J + 13.9)             !%
          CP  = (0.0004*J**2 - 0.2095*J + 39.5)               !%
          ASH = (-0.007*J + 5.9)                              !%
          LIG = (-0.00011*J**2 + 0.0564*J - 2.2)              !%
          Phi = (-0.0051*J - 614.0/J + 5.79)
          IF (EE .GT. 1) THEN
              TdEE = 1-0.0314*(EE-1)*(EE-1)                   !% TdEE = Tedeski dEE
          ELSE
              TdEE = EE                                       !%
          ENDIF
          WdEE = EE-1                                         !% WdEE = Weiss dEE
          NDIP = -8.77+0.33*CP+0.143*NDF                      !%
          ADIP = 85/1000*CP                                   !%
          IADIP= 0.7*ADIP                                     !%
          NDFn = NDF-NDIP+IADIP                               !%
          dNFC = 0.98*PAF*(100-NDFn-CP+IADIP-ASH-EE)          !%
          kdcp = EXP(-0.012*(100*ADIP/CP))                    !g/g
          dCP = CP*kdcp                                       !%
          TdNDF = 0.75*(NDFn-LIG)*(1-Phi*(LIG/NDFn)**0.667)   !% TdNDF = Tedeschi dNDF
          WdNDF = 0.75*(NDFn-LIG)*(1-(LIG/NDFn)**0.667)       !% WdNDF = Weiss dNDF
          tTDN = dNFC+dCP+2.25*TdEE+TdNDF                     !% Total TDN
          aTDN = tTDN-7                                       !% Apparent TDN (non-discounted)
          DMIR = DMI/DMIm                                     
          discount = 0.15-17.371*(DMIR-1)+0.21*aTDN*(DMIR-1)
     &             + 0.341*NDF*(DMIR-1)
     &             - 0.0031*aTDN*NDF*(DMIR-1)
          dTDN = aTDN*(100-discount)/100                      ! Discounted TDN
          wTDN = dNFC+dCP+2.25*WdEE+WdNDF-7                   ! Weiss(1992) TDN, %
          
          ! Pick one: aTDN, dTDN, or wTDN
          TDN = aTDN;
          !TDN = 59.0
          DEd  = 4.409*TDN/100                                !Mcal/kg
          CALL CATTLEGROW (IYR,J,RHUM,WINDSP,Rain,Tmin,Tmean,SW,SR,DEd,
     &         TOPWT,WTLF,GTH,DMI,DMIm,ADG,NEmr,NEAG,NEmI,NEmd,
     &         PDMI,AFA,AFB,AFM,HAF,HMF,AFF,AFT,FAF,ABF,DOP,age,
     &         NEmrB,NEmrPA,NEmrAGE,NEmrCHS,EBF,DFA,CFA)
          
          IF (TOPWT .GT. GTH) THEN                            !GRAZZING ON B/G OCCURS **
              SW = SW + ADG/0.96                              !non-shrunk wt.
	        !INTAKE = DMI*SR                                 !kg/ha
              INTAKE = 0
		    GRAZED = INTAKE * 1.1                           !10% waste, kg/ha
              TOPWTb = TOPWT                                  !kg/ha
              GP = GRAZED/TOPWTb
	        TOPWT = 0.1*TOPWT                               !g/m2
	        WTLF  = 0.1*WTLF                                !g/m2
              GRAZED = 0.1*GRAZED                             !g/m2
              !
              !Get PROLFF and PROSTF from the .SPE file 
              ERRKEY = 'FRGRZ'
              CALL GETLUN('FILEC', LUNCRP)
              OPEN (LUNCRP,FILE = FILECC, STATUS = 'OLD',IOSTAT=ERR)
              IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,0)
              LNUM = 1
              SECTION = '!*PLAN'
              CALL FIND(LUNCRP, SECTION, LNUM, FOUND)
              IF (FOUND .EQ. 0) THEN
                  CALL ERROR(ERRKEY, 1, FILECC, LNUM)
              ELSE
                  CALL IGNORE(LUNCRP,LNUM,ISECT,MOW80)
                  READ(MOW80,'(12X,F6.3,12X,F6.3)',IOSTAT=ERR) 
     &                PROLFF,PROSTF
                  do i=1,5; CALL IGNORE(LUNCRP,LNUM,ISECT,MOW80); end do
                  READ(MOW80,'(2F6.0)',IOSTAT=ERR)
     &                pliglf, pligst
                  CLOSE(LUNCRP)
                  IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILECC,LNUM)
                  ENDIF
              !Update variable values for the calling subroutine
              !GLEAF = FBGL*GRAZED                         !g/m2
              !GSTEM = (1-FBGL)*GRAZED                     !g/m2
              GLEAF = WTLF*GP                            !g/m2
              GSTEM = STMWT*GP                           !g/m2
              GLEAF = MAX(GLEAF,0.0)
              GSTEM = MAX(GSTEM,0.0)
              GVSTG = MAX(MVS,0.0)
		  
              Gtot = GLEAF+GSTEM
              Glfn = GLEAF*pcnl/100
              Gstn = GSTEM*pcnst/100
              Gtotn = Glfn+Gstn
              
              Gcrlf = GLEAF*RHOL
              Gcrst = GSTEM*RHOS

              Gpctn = Gtotn/Gtot*100
              Gplig = (GLEAF*pliglf+GSTEM*pligst)/Gtot*100
              Gpcho = (Gcrlf+Gcrst)/Gtot*100
              Gpctlf = GLEAF/Gtot*100

              WTLFb  = WTLF                              !g/m2
              STMWTb = STMWT                             !g/m2
              WTLF   = MAX((WTLF-GLEAF),0.0)             !g/m2
              STMWT  = MAX((STMWT-GSTEM),0.0)            !g/m2
              WTLFa  = WTLF                              !g/m2
              STMWTa = STMWT                             !g/m2
              TOPWT  = MAX((TOPWT - GLEAF - GSTEM),0.0)  !g/m2
              TOTWT  = MAX((TOTWT - GLEAF - GSTEM),0.0)  !g/m2
              WCRLF  = WTLF*RHOL
              WCRST  = STMWT*RHOS
              WTNLF  = WTLF*PCNL/100.
              WTNST  = STMWT*PCNST/100.
              WTNCAN = WTNCAN - GLEAF*PCNL/100. - GSTEM*PCNST/100.
              IF ((WTLF - WCRLF) .GT. 0.0) THEN
                  WNRLF = MAX (WTNLF-PROLFF*0.16*(WTLF-WCRLF),0.0)
              ELSE
                  WNRLF = 0.0
              ENDIF
              IF ((STMWT - WCRST) .GT. 0.0) THEN
                  WNRST = MAX (WTNST-PROSTF*0.16*(STMWT-WCRST),0.0)
              ELSE
                  WNRST = 0.0
              ENDIF
              AREALF = WTLF*SLA
              XLAI   = AREALF/10000.
              XHLAI  = XLAI
              VSTAGE = GVSTG 

              TOPWT  = 10*TOPWT                           !kg/ha
              WTLF   = 10*WTLF                            !kg/ha
              GRAZED = 10*GRAZED                          !kg/ha
              WTLFb  = 10*WTLFb                           !kg/ha
              WTLFa  = 10*WTLFa                           !kg/ha
              STMWTb = 10*STMWTb                          !kg/ha
              STMWTa = 10*STMWTa                          !kg/ha
              TOPb = WTLFb + STMWTb
              TOPa = WTLFa + STMWTa
              HayIn  = 0.
          ELSE IF (HAY .EQ. 'YES') THEN                   !No grazing but feeding hay
              DEd = DEdBGH                                
              CALL CATTLEGROW (IYR,J,RHUM,WINDSP,Rain,Tmin,Tmean,SW,SR,
     &             DEd,TOPWT,WTLF,GTH,DMI,DMIm,ADG,NEmr,NEAG,NEmI,NEmd,
     &             PDMI,AFA,AFB,AFM,HAF,HMF,AFF,AFT,FAF,ABF,DOP,age,
     &             NEmrB,NEmrPA,NEmrAGE,NEmrCHS,EBF,DFA,CFA)
              SW = SW + ADG/0.96                          !wt.gain
              GRAZED = 0.0
              HayIn  = DMI*SR
              TOPWTb = TOPWT
              GLEAF = FBGL*GRAZED                        !g/m2
              GSTEM = (1-FBGL)*GRAZED                    !g/m2
          ELSE                                            !No grazing, no feeding
              CALL CATTLEGROW (IYR,J,RHUM,WINDSP,Rain,Tmin,Tmean,SW,SR,
     &             DEd,TOPWT,WTLF,GTH,DMI,DMIm,ADG,NEmr,NEAG,NEmI,NEmd,
     &             PDMI,AFA,AFB,AFM,HAF,HMF,AFF,AFT,FAF,ABF,DOP,age,
     &             NEmrB,NEmrPA,NEmrAGE,NEmrCHS,EBF,DFA,CFA)
              ADG = - NEmr/5                              !wt.loss
              SW = SW + ADG/0.96
              INTAKE = 0.0
              GRAZED = 1.1*INTAKE
              DMI = INTAKE
              TOPWTb = TOPWT
              GLEAF = FBGL*GRAZED                         !g/m2
              GSTEM = (1-FBGL)*GRAZED                     !g/m2
              HayIn  = 0.0
          ENDIF                                           !End of grazing
      ELSE                                                !Outside of the grazing period
          ADG = 0.0
          SW = SW + ADG/0.96									
          GRAZED = 0.0                                    !g/m2
          DMI = 1.1*GRAZED
          TOPWTb = TOPWT
          GLEAF = FBGL*GRAZED                             !g/m2
          GSTEM = (1-FBGL)*GRAZED                         !g/m2
      ENDIF
 
      STMWT = 10*STMWT                                    !kg/ha
      GLEAF = 10*GLEAF
      GSTEM = 10*GSTEM
      DMIP = INTAKE/SR/SW*100
      LSR = WTLF/STMWT
      !bGRAZED = TOPWT + INTAKE*1.1
      !aGRAZED = TOPWT
      !
      !Verifying input data by writing to BG_GRAZE.RPT file
      FLAGVP = 1
      IF (FLAGVP .EQ. 1) THEN                             !verif. print out
          IF (FLAGVF .EQ. 0) THEN                         !write verif. output
!              FMT='(/,1X,''SW:     Stocker weight [kg/stocker]'',/,
!     &        1X,''GRAZED: Grazed biomass [kg/ha]'',/,
!     &        1X,''STBLE:  Stubble remained after grazing [kg/ha]'',/,
!     &        1x,''CWADb:  Canopy weight, before grazing [kg/ha]'',/,
!     &        1x,''CWADa:  Canopy weight, after grazing [kg/ha]'',/)'
!             WRITE (56,FMT)
              !
              !Writing daily (grazing dates) outputs to BG_GRAZE.RPT file  
      FMT='(1X,''  YR DOY AGE     SW    TOPWT     WTLF    STMWT   GRAZED
     &    ADG   DMIP   DMIm    DMI   PDMI    AFA    AFB    ABF    HAF
     &    HMF    AFM    AFT   NEmI   NEmr  NEmrB NEmrPA
     & NEmrAG NEmrCH   NEAG    TDN    DFA    CFA'')'
              
              WRITE (56,FMT)
      FMT='(1X,''   -   - [m]   [kg]  [kg/ha]  [kg/ha]  [kg/ha]  [kg/ha]
     &     kg      %   kg/s   kg/s   kg/s      -      -      -      -   
     &  -      -      -    Mcal   Mcal   Mcal   Mcal   Mcal   Mcal 
     &  Mcal      -      -      -'')'
      
             WRITE (56,FMT)
      FMT='(1X,'' (1) (2) (3)    (4)      (5)      (6)      (7)
     &      (8)    (9)   (10)   (11)   (12)   (13)   (14)   (15)   (16)
     &   (17)   (18)   (19)   (20)   (21)   (22)   (23)   (24)   (25)
     &   (26)   (27)   (28)   (29)   (30)'')'
      
              WRITE (56,FMT)
              FLAGVF = 1
          ENDIF
      FMT='(1X,I4,1X,I3,1X,I3,1X,F6.1,1X,F8.1,1X,F8.1,1X,F8.1,1X,F8.1,
     & 1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,
     & 1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,
     & 1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.1,1X,F6.1,1X,F6.1)'
      WRITE (56,FMT) IYR,J,age,SW,TOPWT,WTLF,STMWT,INTAKE,ADG,DMIP,
     & DMIm,DMI,PDMI,AFA,AFB,ABF,HAF,HMF,AFM,AFT,
     & NEmI,NEmr,NEmrB,NEmrPA,NEmrAGE,NEmrCHS,NEAG,TDN,DFA,CFA
      !FMT='(1X,I4,1X,I3,1X,F6.1,1X,F7.1,1X,F6.1,1X,F7.1,1X,F6.2)'  
      !WRITE (56,FMT) IYR,J,SW,TOPWT,WTLF,STMWT,ADG
      !WRITE (56,FMT) IYR,J,SW,TOPWTb,INTAKE,DMIP,ADG,NEAG,NEmI,NEmr,NEgI
      !WRITE (56,FMT) IYR,J,SW,TOPWTb,GRAZED,TOPWT,INTAKE,ADG,DMIP,AFT
      !WRITE (56,FMT) IYR,J,SW,TOPWTb,GRAZED,TOPWT,DEd,NEmd,ADG,dTDN,TDN
      !WRITE (56,FMT) IYR,J,TOPWTb,GRAZED,TOPWT,NEmd,SW,INTAKE,ADG,AFF
          !WRITE (56,FMT) IYR,J,SW,TOPWTb,GRAZED,TOPWT,HayIn,ADG
      !WRITE (56,FMT) IYR,J,SW,DMI,NEmr,NEmd,NEAG,ADG
      !FMT='(1X,I4,1X,I3,1X,F6.1,1X,F6.1,1X,F6.1,1X,F6.1,1X,F6.1)'  
          !WRITE (56,FMT) IYR,J,SW,TOPWTb,GRAZED,Tp,RHp
      ENDIF
      !
      !Checking values for certain variables
      FMT='(1X,''Min. leaf protein compo after N mining:'',1X,F5.3)'
      WRITE (58,FMT) PROLFF
      FMT='(1X,''LCT:'',1X,I3)'
      WRITE (58,FMT) PROSTF
      !
      TOPWT = 0.1*TOPWT							    !g/m2
      !GRAZED = 0.1*GRAZED							!g/m2
      WTLF = 0.1*WTLF                                 !g/m2
      STMWT = 0.1*STMWT                               !g/m2
      GLEAF = 0.1*GLEAF                               !g/m2
      GSTEM = 0.1*GSTEM                               !g/m2
      !
      END SUBROUTINE FOR_GRAZE
      !
!***********************************************************************
!     Variable listing
!***********************************************************************
! Required I/O files:
! GRAZEIN      : Input file
! GRAZE.RPT    : Output file; opened/closed within this subroutine
! GRAZCHK.DAT  : Report file to check; opened/closed within this s/r
!
! Definition of variables passed from/to program forage.for:
! AREALF  real : Total leaf area (one side) (cm2 leaf/m2 ground)
! canht   real : Canopy height (m)
! DEd     real : Digestible energy (Mcal/kg diet)      
! FILECC  char : Path plus filename for species file (*.spe)
! PCNL    real : Percentage of N in leaf tissue (100 g N/g leaf)
! PCNST   real : Percent N in stem tissue (100 g N/g stem)
! RHOL    real : Fraction of leaf which is carbohydrate (g CH20/g leaf)
! RHOS    real : Fraction of stem which is carbohydrate (g CH2O/g stem)
! SLA     real : Specific leaf area (cm2/g leaf)        
! STMWT   real : Dry mass of stem, including C and N (g stem/m2 ground)
! Tmean   real : Average daily temperature (°C) 
! TOPWT   real : Total above-ground weight, including pods (g tissue/m2)
! TOTWT   real : Total weight of crop (g tissue/m2)
! VSTAGE  real : Number of nodes on main stem of plant
! WCRLF   real : Mass of CH2O reserves in leaves (g leaf CH2O/m2 ground)
! WCRST   real : Mass of CH2O reserves in stems (g stem CH2O/m2 ground)
! WNRLF   real : N available for mobil. from leaves above lower mining limit(g N/m2)
! WNRST   real : N available for mobil. from stems above lower mining limit(g N/m2)
! WTLF    real : Dry mass of leaf tissue including C and N (g leaf/m2 ground)
! WTNCAN  real : Mass of N in canopy (g N/m2 ground)
! WTNLF   real : Mass of N in leaves (g leaf N/m2 ground)
! WTNST   real : Mass of N in stems (g stem N/m2 ground)
! XHLAI   real : Leaf area index (m2 leaf/m2 ground)
! XLAI    real : Leaf area (one side) (m2 leaf/m2 ground)
! YDOY    real : Current day of simulation (YYYYDDD) 
!
! Definition of user variables read from file the GRAZIN.DAT:
! DMT1    real : Threshold defining transition bet. avg. & high biomass (kg/ha) 
! DMT2    real : Threshold defining transition bet. avg. & low biomass (kg/ha)
! GTH     real : Threshold for initialization/termination of grazing (kg/ha)
! SR      real : Stocking rate (stockers/ha)
! SWI     real : Initial stocker weight at the beginning of grazing (kg)
!
! Local variable definitions:
! ADG     real : Average daily weight gain per steer (kg carcass/steer)
! AREAH   real :
! DOP      int : Days on pasture
! ERR      int :
! ERRKEY  char :
! GLEAF  real  : The amount of leaf grazed
! GSTEM  real  : The amount of stem grazed
! GVSTG  real  : The vegetative stage after grazing
! FLAGOW   int : Flag to open and write the output data file once
! FLAGRW   int : Flag to read input file and write output report file once
! FLAFVF   int : Flag to write verification output titles one time
! FLAGVP   int : Flag for verification print out
! FM      char : Format statement
! FOUND    int :
! GDB      int : The day of year grazing begins
! GDE      int : The day of year grazing ends
! GRAZED  real : Grazed plus wasted biomass per day (kg/ha) 
! GP      real : Grazing fraction (GRAZED/TOTWTb)
! GYB      int : The year grazing begins
! GYE      int : The year grazing ends
! HayIn   real : Hay intake (kg/ha)
! INTAKE  real : Biomass intake of cattle (kg/ha)
! ISECT    int :
! IYR      int : The current year of simulation
! J        int : The current day of year of simulation
! LNUM     int :
! LUNCRP   int :
! MOW80   char :
! MVS     real : Parameter, number of leaves per stem in the stubble
! PROLFF  real : Minimum leaf protein composition after N mining
! PROSTF  real : Minimum stem protein composition after N mining
! RSPLF   real : Percent of leaf in the stubble (%)
! RTWT    real : Dry mass of root, including C and N, (g root/m2 ground)
! SECTION char :
! STBLE   real : Stubble biomass after grazing (kg/ha)
! SW      real : Stocker weight (kg/steer)
! SWF     real : Final stocker weight at the end of grazing (kg/steer)
! TDN     real : Total digestible nutrients (% DM)
! TOPWTb  real : Top weight before grazing (kg/ha)
!*****************************************************************************