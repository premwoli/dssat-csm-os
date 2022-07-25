C=======================================================================
C  CATTLEGROW, Subroutine, P. Woli, F.M. Rouquette Jr.
C-----------------------------------------------------------------------
C  Governing Model Equations
C  Integrates all cattle growth related variables
C-----------------------------------------------------------------------
C  REVISION       HISTORY
C  11/07/2016 PW  Written
C-----------------------------------------------------------------------
C  Called by:  for_graze
C  Calls:      none
C=======================================================================
!
      SUBROUTINE CATTLEGROW (YR,DOY,RH,WS,Rain,Tmin,Tmean,FBW,SR,DEd, !input
     &           TOPWT,WTLF,GTH,DMI,DMIm,WG,NEmr,NEAG,NEmI,NEmd,
     &           PDMI,AFA,AFB,AFM,HAF,HMF,AFF,AFT,FAF,ABF,DOP,age,
     &           NEmrB,NEmrPA,NEmrAGE,NEmrCHS,EBF,DFA,CFA) !output 
      !          
      IMPLICIT NONE
      !
      INTEGER :: YR,DOY,days,storm,DOP,age
      REAL :: RH,WS,Rain,Tmin,Tmean,FBW,TOPWT,WTLF,SR,Tc,Tp,RHc,RHp
      REAL :: WSc,WSp,WG,DMI,AFA,AFB,HAF,HMF,AFF,AFM,FAF,ABF
      REAL :: NEmrB,NEmrPA,NEmrPA1,NEmrPA2,NEmrPA3,NEmrPA4
      REAL :: MEmrCS,NEmrCS,NEmrHS,NEmrCHS,NEmrAGE,PNA
      Real :: PTAF,PTAF1,PTAF2,PTAF3,DAYL,sinDEC,cosDEC,shift,amplitude
      REAL :: CFA,DFA,UGI,EBF,EQSBW,SRW,DMINC,CETI,PETI,EXI,TSI,MEI,FFM
      Real :: TDN,DEd,MEd,NEmd,NEmdA,NEmr,NEgd,NEAG,LCT,TOI,Km,HE
      Real :: SA,SA1,SA2,GTH,NEd,NEgr,NEr,DMIm,SBW,EBW,a,b,c,d,e
      Real :: PDMI1,PDMI2,PDMI,AFT1,AFT2,AFT,NEgI,NEmI,k_TL,k_LB,k_UB
      !
      !User inputs
      CHARACTER*4, PARAMETER :: AI = 'Y'      !anabolic implant:yes='Y', no='N'
      CHARACTER*4, PARAMETER :: FA = 'NO'     !feed additive:Monensin='MO',none/others='NO'
      CHARACTER*4, PARAMETER :: AT = 'GR'     !animal type:growing='GR',lactating/dry='LD'
      CHARACTER*4, PARAMETER :: BT = 'TA'     !breed type:taurus=TA,indicus=IN,dual=DU,dairy/Holstein=HO,all combined=AL
      INTEGER, PARAMETER   :: iAge =   8      !months old (6-9)
      INTEGER, PARAMETER   :: BCS  =   5      !body cond.score(1-9);normal=5
      INTEGER, PARAMETER   :: HSG  =   2      !heat stress group: Indicus=1,Taurus/cross-bred=2
      REAL, PARAMETER      :: BAF  =   1.0    !breed adj. fac., Brahman: 0.9,Taurus=1,braford=0.95 (Table 19.1, NRoBC 8th ed.)
      REAL, PARAMETER      :: BMC  =   0.074  !basic met. coef.; Brahman: 0.07 w/phiscal activity(Eq.19.2, NRoBC),TAURUS= 0.077, intermediate=0.074
      REAL, PARAMETER      :: CW   =   0.0    !conceptus weight(kg/animal)
      REAL, PARAMETER      :: GAU  =   1.0    !grazing area unit (ha)
      REAL, PARAMETER      :: hair =   0.5    !hair depth(cm);prev=0.5
      REAL, PARAMETER      :: hide =   1.1    !hide thickness(-):thin=.8,avg=1,thick=1.2
      REAL, PARAMETER      :: LAF  =   1.0    !lactation adj. factor(-);lactating=1.2
      REAL, PARAMETER      :: MD   =   0.0    !mud depth (cm)
      REAL, PARAMETER      :: mud  =   1.0    !mud coat(-):no=1,some lower=.8,lower/side=.5,heavy=.2
      REAL, PARAMETER      :: MSBW = 600.0    !mature shrunk body weight(kg/animal) at the targeted empty body fat content; Angus=550(Table 20-1,NASEM) 
      REAL, PARAMETER      :: PC   =   6.0    !position changes (#/d) (Table 12.1; T&F)
      REAL, PARAMETER      :: SEX  =   1.0    !gender effect factor(-); bulls=1.15,others=1
      REAL, PARAMETER      :: TEBF =  28.0    !Targeted (final) empty body fat (%)
      REAL, PARAMETER      :: TRF  =   1.2    !terrain(leveled=1,undulating=1.5,hilly=2)
      REAL, PARAMETER      :: TS   =  18.0    !time standing (h/d) (Table 12.1; T&F)
      REAL, PARAMETER      :: WFL  =   2.0    !walking flat(km/d) (Table 12.1; T&F)
      REAL, PARAMETER      :: WSL  = 0.001    !walking slope(km/d) (Table 12.1; T&F)
      REAL, PARAMETER      :: GTH1 = 2500     !Graz. th.(WTLF?) above which no decrease in consumption (kg DM/ha)  **
      REAL, PARAMETER      :: GTH2 =    0     !Graz. th.(WTLF?) below which animal goes to survival mode (kg DM/ha)**
      REAL, PARAMETER      :: pi   = 3.1415926
      REAL, PARAMETER      :: Lat  = 32.3     !Latitude of the location (deg)
      !
      !Estimating energy values in the diet:
      TDN = 100*DEd/4.409
      MEd  = 0.9611*DEd - 0.2999                             !Galyean et al., 2016
      NEmd = 1.1104*MEd-0.0946*MEd**2+0.0065*MEd**3-0.7783   !Galyean et al., 2016
      NEgd = 1.1376*MEd-0.1198*MEd**2+0.0076*MEd**3-1.2979   !Galyean et al., 2016
      
      !Estimating the dynamic value of the correction factor (KAPPA); 100 ≤ KAPPA ≤ 285
      k_LB = -0.000000206*DOY**3 +0.000078*DOY**2   -0.005801*DOY -0.018
      k_TL = -0.000000206*DOY**3 +0.00008643*DOY**2 -0.009171*DOY +0.349
      k_UB = -0.000000206*DOY**3 +0.0000933*DOY**2  -0.012055*DOY +0.693
      
      !IF (DOY .LE. 226) THEN
          !k_TL =0.0000645*DOY**2 -0.01335*DOY + 1.69      !if w/o physical activity
          !k_TL = 0.00000044*DOY**3 -0.0001053*DOY**2 +0.009117*DOY +0.74     !OLD
          !k_LB = 0.00000022*DOY**3 -0.000016*DOY**2  -0.002496*DOY +1.15
          !k_UB = 0.00000063*DOY**3 -0.00018*DOY**2   +0.0187*DOY   +0.42
          !k_TL = 0.00000154*DOY**3 -0.0005558*DOY**2 +0.070656*DOY -2.00
          !k_TL = 0.0000008*DOY**3  -0.000218*DOY**2 + 0.0209*DOY + 0.36
          !k_LB = 0.0000012*DOY**3  -0.000413*DOY**2  +0.05133*DOY  -1.24
          !k_UB = 0.00000187*DOY**3 -0.000688*DOY**2  +0.0883*DOY   -2.70
      !ELSE
          !k_TL =0.000010025*DOY**3 -0.0080203*DOY**2 +2.1062*DOY-180.16 !if w/o PA
          !k_TL = 0.000018861*DOY**3 -0.014965*DOY**2 +3.91*DOY   -334.60 !OLD
          !k_LB = 0.000018922*DOY**3 -0.0148845*DOY**2+3.85466*DOY-327.09
          !k_UB = 0.00001944*DOY**3  -0.0155288*DOY**2+4.08455*DOY-351.71
          !k_TL = 0.00004277*DOY**3 -0.033705*DOY**2 +8.77075*DOY -751.10
          !k_TL = 0.0000266*DOY**3 -0.021373*DOY**2 +5.648328*DOY -488.79
          !k_LB = 0.00004095*DOY**3 -0.032075*DOY**2 +8.29355*DOY -705.72
          !k_UB = 0.00004559*DOY**3 -0.03609*DOY**2  +9.4345*DOY  -811.53
      !ENDIF
      !
      !Updating age
      age = iAge + INT(DOP/30)
      !
      !Predicting dry matter intake:
      !=========================================================================
      !Potential/maximum/basal intake-
      IF (NEmd .LT. 1.0) THEN
          NEmdA = 0.95
      ELSE
          NEmdA = NEmd
      ENDIF
      IF (iAge .LE. 12) THEN                               !Eq.19.89, NRC,2018
          PDMI1= FBW**0.75*(0.2435*NEmd-0.0466*NEmd**2-0.1128)/NEmdA
      ELSE
          PDMI1= FBW**0.75*(0.2435*NEmd-0.0466*NEmd**2-0.0869)/NEmdA
      ENDIF												
      PDMI2 =  FBW*(0.012425+0.019218*NEmd-0.007259*NEmd**2)!Eq.19.92, NRC,2016
      PDMI = PDMI2                                          !pick PDMI1 or PDMI2
      !=========================================================================
      !
      !Adjustment factor for feed additive (AFA)
      IF (AI .EQ. 'N') THEN
          IF (FA .EQ. 'MO') THEN
              AFA = 0.94*0.97
          ELSE
              AFA = 0.94
          ENDIF
      ELSE IF (AI .EQ. 'Y') THEN
          IF (FA .EQ. 'MO') THEN
              AFA = 0.97
          ELSE
              AFA = 1.0
          ENDIF
      ELSE
          AFA = 1.0
      ENDIF
      !
      !Adjustment factor for breed type (AFB)
      BREED1: SELECT CASE (BT)
      CASE ('TA') 
          AFB = 1.0  !under grazing conditions
      CASE ('IN')
          AFB = 1.0 !1.05
      CASE ('DU') 
          AFB = 1.04 !1.10
      CASE ('HO') 
          AFB = 1.08 !0.9
      END SELECT BREED1
      !
      !Adjustment factor for mud depth (AFM)
      AFM = 1-0.01*MD
      !
      !Adjustment factor for body fat (AFF)
      SBW = 0.96*FBW                                       !shrunk body wt.
      EBW = 0.891*SBW                                      !Eq.19-43; NRC(2016)
      BREED2: SELECT CASE (BT)
      CASE ('TA') 
          EBF = -1.8372+0.061*EBW+0.00047*EBW**2           !Tab.15-1; taurus;T&F(2018)
      CASE ('IN')
          EBF = 13.0228-0.096*EBW+0.00072*EBW**2           !Tab.15-1; indicus
      CASE ('HO') 
          EBF = 6.282-0.0365*EBW+0.00055*EBW**2            !Tab.15-1; Holstein
      CASE ('AL') 
          EBF = 0.91524+0.01286*EBW+0.00054*EBW**2         !Fig.15.1; all combined.
      END SELECT BREED2
      SRW = 399.89 - 10.159*TEBF + 0.4621*TEBF**2          !Eq. 19.44, NRC(2016); Eq. 15.16; T&F (2018)
      EQSBW = (SBW-CW)*SRW/MSBW                            !Eq. 19.46, NRC(2016); Eq. 15.15, T&F (2018)
      IF (EQSBW .GE. 350) THEN                             !Eq. 19.100
          AFF = 0.7714+0.00196*EQSBW-0.00000371*EQSBW**2
      ELSE
          AFF = 1.0
      ENDIF
      !
      !--------------------------------------------------------------------------
      !Adjustment factor for temperature (AFT)
      !Compute daylength (DAYL)
      IF (MOD(YR,4) .EQ. 0) THEN 
          days = 366
      ELSE
          days = 365
      ENDIF
      sinDEC = -SIN(23.45*pi/180)*COS(2*pi*(DOY+10)/days)
      cosDEC = SQRT(1-sinDEC*sinDEC)
      shift =  sinDEC * SIN(Lat*pi/180)
      amplitude = cosDEC * COS(Lat*pi/180)
      DAYL = 12*(1+(2/PI)*ASIN(shift/amplitude))
      
      WS = MIN(32.0,WS/24.0)  		            !km/d to km/h
      IF ((Rain .GE. 100) .OR. (WS .GE. 32)) THEN !
          storm = 1
      ELSE
          storm = 0
      ENDIF
      CALL cattle_PrevWTH(DOY,Tc,Tp,RHc,RHp,WSc,WSp)
      CETI = (27.88-0.456*Tc+0.010754*Tc**2
     &       -0.4905*RHc+0.00088*RHc**2
     &       +1.1507*WSc-0.126447*WSc**2
     &       +0.019876*Tc*RHc-0.046313*Tc*WSc+0.4167*DAYL) 
      !The NRC (2016) method; Eq.19.101:103)
      DMINC = (119.62-0.9708*CETI)/100            !DMI for night cooling ad.fac.
      IF (Tc .LE. -20) THEN                       !is Tmin = least nt. temp.?
          AFT1 = 1.16
      ELSE IF ((Tc .GT. -20) .AND. (Tc .LE. 20)) THEN
          AFT1 = 1.043-0.0044*Tc+0.0001*Tc**2
      ELSE IF ((Tc .GT. 20) .AND. (Tmin .GT. 20)) THEN !Tmin=lowest n.t.?
          AFT1 = DMINC
      ELSE IF ((Tc .GT. 20) .AND. (Tmin .LE. 20)) THEN
          AFT1 = (1-DMINC)*0.75 + DMINC
      ELSE IF ((storm .EQ. 1) .AND. (Tc .LE. 13)) THEN
          AFT1 = 1.0
      ENDIF
      !Another method (Eq.10.3; T&F, 2018)
       IF (CETI .LT. -10) THEN
          a = 1.07
          b = 0.09
          c = -14.4
          d = 0.0046
          e = -1.04
      ELSE IF ((CETI .GE. -10) .AND. (CETI .LT. 0)) THEN
          a = 1.05
          b = 0.02
          c = -5.5
          d = -9.7*10**(-5)
          e = -0.85
      ELSE IF ((CETI .GE. 0) .AND. (CETI .LT. 10)) THEN
          a = 1.03
          b = 0.02
          c = 4.5
          d = -0.0019
          e = -0.9
      ELSE IF ((CETI .GE. 10) .AND. (CETI .LT. 20)) THEN
          a = 1.0
          b = 0.03
          c = 14.5
          d = -0.0019
          e = -0.9
      ELSE IF ((CETI .GE. 20) .AND. (CETI .LT. 29)) THEN
          a = 0.9
          b = 0.1
          c = 24.7
          d = -0.00145
          e = -0.96
      ELSE
          IF (Tmin .GE. 20) THEN
              a = 0.9
              b = 0.1
              c = 24.7
              d = -0.00145
              e = -0.96
          ELSE
              a = 0.65
              b = 0.25
              c = 35.5
              d = 0.33
              e = -1.21
          ENDIF
      ENDIF
      AFT2 = a+b*(2*e*(LOG(EXP((CETI+d/2)/e)+EXP(c/e)) -
     &                 LOG(EXP((c+d/2)/e)+EXP(CETI/e)))+d)/(2*d)     
      AFT = AFT2                       !pick AFT1 or AFT2 
      
      !--------------------------------------------------------------------------
      !Adjustment factor for forage availability
      !The Loewer (1987) approach (herbage allowance factor; HAF):
      CFA = 10.454 - 0.1353*TDN        ! Critical herbage allowance(kgDM/kgBW)
      DFA = TOPWT/SR/SBW               ! Daily herbage allowance(kgDM/kgBW)
      IF (DFA .LT. CFA) THEN
          HAF = 2*DFA/CFA - (DFA**2)/(CFA**2)
      ELSE
          HAF = 1.0
      ENDIF
      !
      !The Rayburn (1992) approach (herbage mass factor; HMF):         
      DFA = TOPWT/(SR*SBW)                                 !kgDM/kgSBW
      UGI = PDMI/SBW                                       !unadjusted grazing intake (kgDM/kgSBW)
      IF ((DFA .GT. 4*UGI) .OR. (TOPWT .GT. 1150)) THEN
          HMF = 1.0                                        !Eq.19.97;NRC(2016)
      ELSE
          HMF = (0.17*TOPWT - 0.0000764*TOPWT**2 + 2.4)/100!if TOPWT > 0.
      ENDIF
      
      !The Zhang et al. (2008) approach:
      !IF (TOPWT .LT. GTH2) THEN
      !    FAF = 0.0001*SBW/PDMI
      ! ELSE IF ((TOPWT .GE. GTH2) .AND. (TOPWT .LE. GTH1)) THEN
      !    FAF = 1-(1-0.0001*SBW/PDMI)*(GTH1-TOPWT)/(GTH1-GTH2)
      !ELSE
      !    FAF = 1.0
      !ENDIF
      
      !
      !Now computing actual intake-
      DMI = PDMI*AFA*AFB*AFM*AFT*HAF       
      !DMI= PDMI*AFA*AFB*AFM*AFT*HAF  !use HAF, HMF, or FAF
      
      !**************************************************************************************************
      !Computing animal maintenance requirement
      !==================================================================================================
      !Basal maintenance requirement (NEmrB)-
      PNA = 0.8+0.05*(BCS-1)
      IF (Tp .GT. 20) THEN
          PETI = (27.88-0.456*Tp+0.010754*Tp**2
     &           -0.4905*RHp+0.00088*RHp**2
     &           +1.1507*WSc-0.126447*WSc**2
     &           +0.019876*Tp*RHp-0.046313*Tp*WSc+0.4167*DAYL)
      ELSE
          PETI = Tp
      ENDIF
      PTAF1 = 0.0007*(20-Tp)								            !Eq.12.2, T%F
      PTAF2 = ((88.426-0.785*Tp+0.0116*Tp**2)-77)/1000		        !Eq.12.5, T%F
      PTAF3 = ((88.426-0.785*PETI+0.0116*PETI**2)-77)/1000	        !Eq.12.7, T%F
      PTAF  = PTAF3                                                   !pick PTAF1,PTAF2,or PTAF3
      NEmrB = (BMC*BAF*LAF*PNA*SEX+PTAF)*SBW**0.75                    !Eq.19.1:4; NRC,2016
      !==================================================================================================
      !
      !Adjustment factor for physical activity (NEmrPA):Eq.12.12:14,T&F
      IF (TOPWT .GT. GTH) THEN                                        !grazing occurs
          NEmrPA1 = NEmrB*0.06                                        !6% assumed for stockers
          NEmrPA2 = NEmrB*((1.05+0.074*GAU)-1)                        ! smallest (+12.4%)
          NEmrPA3 = ((0.006*DMI*(0.9-TDN/100)+
     &               (0.05*TRF)/(TOPWT/1000+3))*SBW)/4.184            !greatest (+22-27%)
          NEmrPA4 = ((0.1*TS+0.062*PC+0.621*WFL+6.69*WSL)*FBW)/1000 	!intermediate (+18-21%)
          NEmrPA  = NEmrPA1                                          !pick NEmrPA1,2,3,or 4
      ELSE
          NEmrPA  = 0.0
      ENDIF
      
      !Adjustment factor for cold stress (NEmrCS); hold on until for winter pasture	  
      EXI = (6.1816-0.5575*WSc+0.0152*WSc**2+5.298*hair-0.4297*hair**2 -
     &      0.1029*WS*hair)*mud*hide
      IF (iAge .LE. 1) THEN
          TSI = 2.5
      ELSE IF ((iAge .GT. 1) .AND. (iAge .LE. 6)) THEN
          TSI = 6.5
      ELSE IF ((iAge .GT. 6) .AND. (iAge .LE. 12)) THEN
          TSI = 5.1875+0.3125*BCS
      ELSE IF (iAge .GT. 12) THEN
          TSI = 5.25+0.75*BCS
      ENDIF
      TOI = EXI + TSI
      MEI = DMI*MEd					
      FFM = NEmrB/NEmd
      IF (AT .EQ. 'GR') THEN
          HE = MEI-(DMI-FFM)*NEgd
      ELSE IF (AT .EQ. 'LD') THEN
          HE = MEI-(DMI-FFM)*NEmd
      ENDIF
      SA1 = 0.09*FBW**0.67				!Eq.12.16;T&F
      SA2 = 0.14*FBW**0.57			    !Eq.12.24;T&F 
      SA  = SA1                       !pick SA1 or SA2
      LCT = 39-0.85*TOI*HE/SA
      IF (LCT .GT. Tc) THEN           
          MEmrCS = SA*(LCT-Tc)/TOI
      ELSE
          MEmrCS = 0.0
      ENDIF
      Km = NEmd/MEd					
      NEmrCS = MEmrCS * Km
      !-------------------------------------------------------------------------------------
      !Adjustment factor for heat stress (NEmrHS) play with various values!
      IF (HSG .EQ. 1) THEN                          !indicus T&F Eq.12.25:26
          IF ((Tc .GT. 35) .AND. (Tc .LT. 38)) THEN !CETI/Tc/T?
              NEmrHS = NEmrB*0.07                   !rapid/shallow panting; Temp threshold?
          ELSE IF (Tc .GT. 38) THEN                 !if Temp >= 30C (NRC, 1996)
              NEmrHS = NEmrB*0.18                   !open mouth panting; Temp threshold? 
          ELSE
              NEmrHS = 0.0
          ENDIF
      ELSE IF (HSG .EQ. 2) THEN                     !taurus
          IF (CETI .GT. 20) THEN                    
             NEmrHS = NEmrB*(0.09857-0.01343*CETI+0.000457*CETI**2)
          ELSE
             NEmrHS = 0.0
          ENDIF
      ENDIF
      !
      IF (NEmrCS .GT. 0.0) THEN
          NEmrCHS = NEmrCS
      ELSE IF (NEmrHS .GT. 0.0) THEN
          NEmrCHS = NEmrHS
      ELSE
          NEmrCHS = 0.0
      ENDIF
      !------------------------------------------------------------------------------------
      !AAdjustment factor for age (NEmrAGE)
      IF (iAge .LT. 72) THEN		   		!Eq.12.49;T&M,2016
          NEmrAGE = NEmrB*(1-EXP(-0.08*(iAge/12.)))
      ELSE
          NEmrAGE = NEmrB*(1-EXP(-0.08*6))
      ENDIF
      !
      !Adjusted maintenance requirement (NEmr)-
      NEmr = NEmrB + NEmrPA + NEmrCHS - NEmrAGE 
      !NEmr = NEmrB + NEmrPA + NEmrCHS - NEmrAGE
      !
      !====================================================================================	  
      !Computing animal daily weight gain (WG)
      DMIm = NEmr/NEmd                    !DMI at maintenance(kg/hd/d) (Eq. 19-5, NASEM)
      NEgI = DMI*NEgd                     !NE intake for growth
      NEmI = DMI*NEmd                     !NE intake for maintenance
      NEAG = (DMI - DMIm)*NEgd            !NEld,NEcd not available!NE required for growth?
                                          !Zhang et al., 2008; Rodriguez et al, 1990
                                          !Eq. 15.12,T&F(2018);Eq.12.7.NRC(2016)
      !
      !WG adjustment factor for iAge-BCS effects (ABF)
      IF (iAge .LE. 5) THEN
          ABF = -0.035*BCS + 1.073
      ELSE IF ((iAge .GT. 5) .AND. (iAge .LE. 10)) THEN
          ABF = -0.053*BCS + 1.269
      ELSE                                                  !but should not exceed 15
          ABF = -0.073*BCS + 1.469
      ENDIF
      !
      IF (NEmI .LT. NEmr) THEN            !If NEtaken<NEreq, WL
          WG = -(NEmr - NEmI)/4.85        !assumed 5 Mcal/kg beef (Tedeschi suggestion: 4.85)
      ELSE
          WG = k_TL + 15.54*(NEAG/SBW**0.75)**0.9116  !Eq.12.7, NRC(2016);for shrunk wt. gain
          !WG = k_TL*ABF*15.54*(NEAG/SBW**0.75)**0.9116  !Eq.12.7, NRC(2016);for shrunk wt. gain
          !WG = 15.54*(NEAG/SBW**0.75)**0.9116  !for large-frame steers;
                                               !13.91 for med-frame steers
      ENDIF
      WS = WS*24.                         !back to km/d
      !
      END SUBROUTINE CATTLEGROW
!***********************************************************************
!     Variable listing
!***********************************************************************
! Definition of variables passed from/to program for_graze.for:
! FBW     real : Animal full body weight (kg carcass/animal)
! DMI     real : Actual dry matter intake (kg forage/animal/d)
! DOG     int  : Days on pasture grazing
! GTH     real : Grazing threshold (kg DM/ha)
! Rain    real : Daily total rainfall (mm)
! RH      real : Daily average relative humidity (%)
! SR      real : Animal stocking rate (#/ha)
! Tmean   real : Daily average temperature (°C)
! Tmin    real : Daily minimum temperature (°C) 
! TOPWT   real : Total above-ground weight, including pods (g/m2)
! WTLF    real : Leaf weight (g/m2)
! WG      real : Average daily weight gain (kg/animal/d)
! WS      real : Wind speed (km/d)
! YR      int  : the year of simulation
!
! Local variable definitions:
! AFA     real : Adjustment factor for feed additive (-)
! AFB     real : Adjustment factor for breed type (-)
! ABF     real : Age-BCS adjustment factor for weight gain (WG)
! AFF     real : Adjustment factor for body fat (-)
! HAF     real : Herbage allowance adjustment factor for grazing animals (-)
! HMF     real : Herbage mass adjustment factor for grazing animals (-)
! AFM     real : Adjustment factor for body mud (-)
! AFT     real : Adjustment factor for temperature (-)
! age      int : Age of an animal (months)
! AI       int : Anabolic implant (yes/no)
! AT       int : Animal type (1-2)
! BAF     real : Breed adjustment factor(-)
! BCS      int : Body condition score (scale:1-9)
! BMC      int : Basic metabolism coefficient (Mcal/FBW**75/d)
! BT       int : Breed type (1=Bos taurus/indicus,2=dual purpose animals,3=dairy animals)
! CETI    real : Current effective temperature index (°C)
! CFA     real : Critical forage allowance below which potential intake decreases (kgDM/kgFBW)
! CW      real : Conceptus weight (kg/animal)
! days     int : The total number of days in a year
! DAYL    real : Day length hours (#)
! DFA     real : Daily forage allowance (gDM/kgFBW)
! DEd     real : Digestible energy in the diet (Mcal/kg DM)
! DMINC   real : DMI for night cooling adjustment factor (-)
! TEBF    real : Targeted (final) empty body fat (%)
! EBW     real : Empty body weight (kg/animal)
! EQSBW   real : Equivalent shrunk body weight (kg/animal)
! EXI     real : External insulation (°C.m2.d/Mcal)
! FA       int : Feed additive (1=monensin, 2=none/lasalocid/laidlomycin)
! FFM     real : Feed for maintanance (kg/d)
! GAU      real : Grazing unit area (ha)
! F     real : Green forage available (kg/ha) = TOPWT
! hair    real : Hair depth (cm)
! HE      real : Heat production (Mcal/d)
! hide    real : Hide thickness factor (-)
! Km      real : Partial efficiency of utilization of MEd for NEmd (-)
! LAF     real : Lactation adjustment factor (-)
! LCT     real : Lower critical temperature (°C)
! MD      real : Mud depth on the animal body (cm)
! MDMI    real : Auxiliary variable for PDMI (kg forage/animal/d)
! MEd     real : Metabolizable energy in the diet (Mcal/kg DM)
! MEI     real : Metabolizable energy intake (Mcal/d)
! MEmrCS  real : MEmr due to cold stress (Mcal/animal/d)
! mud     real : Hair mud coat factor (-)
! MSBW    real : Mature shrunk body wt.(kg/animal), typically 96% of mature full body wt.
! NEAG    real : Net energy available for weight gain (Mcal/animal/d)
! NEgd    real : Net energy in the diet for weight gain (Mcal/kg DM)
! NEmd    real : Net energy in the diet for maintenance (Mcal/kg DM)
! NEmdA   real : Auxiliary variable for NEmd (Mcal/kg DM)
! NEmr    real : Net energy required by the cattle for maintenance (Mcal/animal/d)
! NEmrACT real : NEmr due to physical activity (Mcal/animal/d); 3 options
! NEmrAGE real : NEmr due to age (Mcal/animal/d)
! NEmrB   real : Basic net energy maintenance requirement (Mcal/animal/d)
! NEmrCS  real : NEmr due to cold stress (Mcal/animal/d)
! NEmrHS  real : NEmr due to heat stress (Mcal/animal/d)
! NEmrCHS real : NEmr due to cold or heat stress (Mcal/animal/d)
! PANT     int : Panting due to heat stress (breed-specific)
! PC      real : Position changes (no./d)
! PDMI    real : Potential voluntary dry matter intake (kg forage/animal/d)
! PETI    real : Previous effective temperature index (°C)
! PNA     real : Previous nutrition adjustment factor (-)
! PTAF    real : Previous temp. acclimatization fac. (Mcal/FBW**.75/d);3 options
! RHp     real : Daily average relative humidity in the previous month(%)
! SA      real : Surface area (m2): 2 options
! SBW     real : Shrunk body weight (kg/animal)
! SEX     real : Gender effect factor (-)
! SRW     real : Standard reference weight (kg/animal)
! storm    int : Storm exposure (0=no,1=yes)
! TOI     real : Total insulation (°C.m2.d/Mcal)
! Tp      real : Daily average temperature in the previous month (°C)
! TRF     real : Terrain factor (-)
! TS      real : Time standing (h/d)
! TSI     real : Tissue insulation (°C.m2.d/Mcal)
! UGI     real : Unadjusted grazing DMI (g DM/kg BW)
! WFL     real : Walking flat (km/d)
! WSL     real : Walking sloped (km/d)
!*****************************************************************************