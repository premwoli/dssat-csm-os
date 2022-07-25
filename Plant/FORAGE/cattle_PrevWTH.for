      SUBROUTINE cattle_PrevWTH (DOY,                             !input
     &           Tc,Tp,RHc,RHp,WSc,WSp)                           !output
      !
      implicit none
      INTEGER :: DOY,AS,i,j,k,WTHCOUNT,WTHLUN,ISECT               !indices
      INTEGER,PARAMETER  :: MWm=30, MWw=7                         !moving window size
      INTEGER,ALLOCATABLE,DIMENSION(:) :: DATE,WIND               !arrayS
      REAL,ALLOCATABLE,DIMENSION(:) :: TMAX,TMIN,TEMP,RHUM,WSPD   !arrays
      REAL,DIMENSION(1:DOY) :: maTEMPm,maRHUMm,maWSPDm            !arrays
      REAL,DIMENSION(1:DOY) :: maTEMPw,maRHUMw,maWSPDw            !arrays
      REAL               :: SumTEMPm,SumRHUMm,SumWSPDm,Tp,RHp,WSp     
      REAL               :: SumTEMPw,SumRHUMw,sumWSPDw,Tc,RHc,WSc      
      CHARACTER*80  WTH80
      CHARACTER*200 FMT
      !
      IF (.NOT. ALLOCATED(TMAX)) THEN
      WTHLUN=999
      OPEN (UNIT=WTHLUN,FILE='C:\DSSAT48\BG_Grazing\TXOV2001.WTH',
     & STATUS='old')                                              !assign WTH year file
      REWIND(WTHLUN)
      !
      ISECT = 0
      WTHCOUNT = 0
      DO WHILE (ISECT.EQ.0)
        READ (WTHLUN,'(A80)',IOSTAT=ISECT) WTH80
        IF (WTH80(1:1).NE."@" .AND. 
     &      WTH80(1:1).NE."*" .AND.
     &      WTH80(1:1).NE."!" .AND.
     &      WTH80(1:2).NE."  ".AND.
     &      ISECT.EQ.0)THEN
            WTHCOUNT = WTHCOUNT + 1
        ENDIF
      END DO
      !
      REWIND(WTHLUN)
      IF (WTHCOUNT.GT.0) THEN
        ALLOCATE(DATE(WTHCOUNT),TMAX(WTHCOUNT))
        ALLOCATE(TMIN(WTHCOUNT),WIND(WTHCOUNT),RHUM(WTHCOUNT))
        ALLOCATE(TEMP(WTHCOUNT),WSPD(WTHCOUNT))
      ELSE
        ALLOCATE(TMAX(1))
        TMAX (1) = -99
      ENDIF
      !
      k = 0
      ISECT = 0
      DO WHILE (ISECT.EQ.0)
        READ (WTHLUN,'(A80)',IOSTAT=ISECT) WTH80
        IF (WTH80(1:1).NE."@" .AND. 
     &      WTH80(1:1).NE."*" .AND.
     &      WTH80(1:1).NE."!" .AND.
     &      WTH80(1:2).NE."  ".AND.
     &      ISECT.EQ.0)THEN
            k = k + 1
            FMT = '(1I5,6X,2F6.1,12X,1I6,12X,1F6.1)'
            READ (WTH80,FMT,IOSTAT=ISECT)
     &            DATE(k),TMAX(k),TMIN(k),WIND(k),RHUM(k)
            TEMP(k) = (TMAX(k)+TMIN(k))/2
            WSPD(k) = WIND(k)/24.0        !km/d to km/h
        ENDIF
      END DO
      ENDIF
      !
      DO i = 1, DOY
        !compute 30-day moving averages
        IF (i .LT. MWm) THEN
            SumTEMPm = 0.0
            SumRHUMm = 0.0
            SumWSPDm = 0.0
            DO j = 1,i                   
                SumTEMPm = SumTEMPm + TEMP(j)
                SumRHUMm = SumRHUMm + RHUM(j)
                SumWSPDm = SumWSPDm + WSPD(j)
            END DO
            maTEMPm(i) = SumTEMPm/i
            maRHUMm(i) = SumRHUMm/i
            maWSPDm(i) = SumWSPDm/i
        ELSE
            SumTEMPm = 0.0
            SumRHUMm = 0.0
            SumWSPDm = 0.0
            DO j = i+1-MWm, i                   
                SumTEMPm = SumTEMPm + TEMP(j)
                SumRHUMm = SumRHUMm + RHUM(j)
                SumWSPDm = SumWSPDm + WSPD(j)
            END DO
            maTEMPm(i) = SumTEMPm / MWm
            maRHUMm(i) = SumRHUMm / MWm
            maWSPDm(i) = SumWSPDm / MWm
        ENDIF
        !compute 7-day moving averages
        IF (i .LT. MWw) THEN
            SumTEMPw = 0.0
            SumRHUMw = 0.0
            SumWSPDw = 0.0
            DO j = 1,i                   
                SumTEMPw = SumTEMPw + TEMP(j)
                SumRHUMw = SumRHUMw + RHUM(j)
                SumWSPDw = SumWSPDw + WSPD(j)
            END DO
            maTEMPw(i) = SumTEMPw/i
            maRHUMw(i) = SumRHUMw/i
            maWSPDw(i) = SumWSPDw/i
        ELSE
            SumTEMPw = 0.0
            SumRHUMw = 0.0
            SumWSPDw = 0.0
            DO j = i+1-MWw, i                   
                SumTEMPw = SumTEMPw + TEMP(j)
                SumRHUMw = SumRHUMw + RHUM(j)
                SumWSPDw = SumWSPDw + WSPD(j)
            END DO
            maTEMPw(i) = SumTEMPw / MWw
            maRHUMw(i) = SumRHUMw / MWw
            maWSPDw(i) = SumWSPDw / MWw
        ENDIF
      END DO
      Tc  = maTEMPw(DOY)
      Tp  = maTEMPm(DOY)
      RHc = maRHUMw(DOY)
      RHp = maRHUMm(DOY)
      WSc = maWSPDw(DOY)
      WSp = maWSPDm(DOY)
      !RETURN                                
      END SUBROUTINE cattle_PrevWTH