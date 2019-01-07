SUBROUTINE MAKE_LEVEL2B(year,month,iday,itime_start,itime_end,&
     EOT,options,previous,M,time_of_day)  
  ! Offline function that samples the model data so that the
  ! local time of the model data exactly match the satellite overpass time
  !
  ! Reads in a whole day at a time, plus the first UTC of the next
  ! day. This extra time step is kept around for the next day to avoid
  ! rereading the same data. The time steps are used to interpolate
  ! the data.
  !
  !       For each day and the correspond UTC
  ! itime 1  2  3  4  5  6  7  8  9 10 11 12 13
  ! UTC:  6 12 18 24  6 12 18 24  6 12 18 24  6
  ! day1: x  x  x  x  x 
  ! day2:             x  x  x  x  x
  ! day3:                         x  x  x  x  x
  ! 
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,        ONLY: WP
  USE HANDY,             ONLY: &
       BUILD_FILENAME,         &
       CHECK_FILE
  USE MODEL_INPUT,       ONLY: &
       ALLOCATE_MODEL_MATRIX,  &
       DEALLOCATE_MODEL_MATRIX,&
       INITIALISE_MODEL_MATRIX,&
       MODEL_AUX,              &
       MODEL_TYPE,             &
       READ_MODEL
  USE MY_MATHS,            ONLY: DAY_SINCE_DAY, DAYSINMONTH
  USE NAMELIST_INPUT,      ONLY: NAME_LIST

  IMPLICIT NONE

  TYPE(name_list), INTENT(in)        :: options
  INTEGER, INTENT(in)                :: year, month, iday, itime_start,itime_end
  REAL(wp),INTENT(in)                :: EOT
  TYPE(model_type), INTENT(inout)    :: previous       ! from previous day
  TYPE(model_type), INTENT(inout)    :: M              ! sampled model
  REAL(wp),INTENT(out)               :: time_of_day(M%aux%ngrids)

  !INTERNAL
  ! one for every timestep used in for this day (1+ntimes)
  TYPE(model_type),ALLOCATABLE       :: fullday(:)

  TYPE(model_aux)      :: aux
  REAL(wp)             :: UTC1, t1, t2 
  REAL(wp),ALLOCATABLE :: T0(:,:),T0_b(:,:),g_c2u1(:,:),g_c2u2(:,:)
  REAL(wp),ALLOCATABLE :: c2u1(:),c2u2(:),c2u1_2(:,:),c2u2_2(:,:)
  INTEGER              :: nlat,nlon,nlev,ngrids
  INTEGER              :: i
  REAL(wp)             :: deltaTime
  INTEGER              :: tmpD,tmpM,tmpY,tmpItime,nsteps
  LOGICAL,ALLOCATABLE  :: isL2b(:),g_isL2b(:,:),isL2b_2(:,:)

  ! equatorial overpass time converted to UTC
  LOGICAL :: interpolate
  INTEGER,PARAMETER :: step = 10
  CHARACTER(1000)  :: rx,file,fileDay

  aux   = M%aux
  ngrids= aux%ngrids
  nlev  = aux%nlev
  nlon  = aux%nlon
  nlat  = aux%nlat

  ALLOCATE(T0 (nlon  ,nlat),&
       T0_b   (nlon  ,nlat),&
       g_c2u1 (nlon  ,nlat),&
       g_c2u2 (nlon  ,nlat),&
       g_isL2b(nlon  ,nlat),&
       c2u1   (ngrids     ),&
       c2u2   (ngrids     ),&
       isL2b  (ngrids     ),&
       c2u1_2 (ngrids,nlev),&
       c2u2_2 (ngrids,nlev),&
       isL2b_2(ngrids,nlev))

  interpolate = options%L2b%interpolate
  rx          = options%paths%model_input_regexp  
  deltaTime   = aux%ref%deltaTime
  nsteps      = itime_end-itime_start+1

  ALLOCATE(fullday(nsteps+1)) ! +1 for the additional step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE (0,'(a)') '-------------- enter  make_level2B ------------------'
  WRITE (0,'(a,5i8)') ' year, month, iday, itime_start:itime_end: ',&
       year, month, iday, itime_start,itime_end 
  WRITE (0,'(a,f8.2)') ' EOT : ',EOT
  WRITE (0,'(a)') '-----------------------------------------------------'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1 FORMAT('Reading file from the ',a,' day:',a,', time step =',i3)

  ! ==============================
  ! 
  !     Put together the data
  !
  ! ==============================
  fileDay = BUILD_FILENAME(options%paths%model_input_regexp,&
       year,month,iday, model=options%model)

  ! --------------------
  ! Allocate
  !
  DO i = 1,nsteps+1
     CALL ALLOCATE_MODEL_MATRIX(fullday(i),ngrids,nlev)
     CALL INITIALISE_MODEL_MATRIX(fullday(i),ngrids,nlev)
     fullday(i)%aux=M%aux
  END DO

  ! --------------------
  ! Time step from PREVIOUS run
  !
  tmpItime=itime_start
  file=fileday
  i=1
  IF (.NOT.ALLOCATED(previous%T)) THEN
     ! then I need to read the first time step
     PRINT 1,'given',TRIM(file),tmpItime
     fullday(i)%aux%netcdf_file=TRIM(file)
     CALL READ_MODEL(fullday(i),tmpItime,options)
  ELSE
     PRINT '("Loading data from time step =",i3)',tmpItime

     fullday(i)=previous
  END IF

  ! --------------------
  ! Read all the model data for the GIVEN day
  !
  file=fileday
  DO i= 2,nsteps
     tmpItime=tmpItime+1
     PRINT 1,'given',TRIM(file),tmpItime
     fullday(i)%aux%netcdf_file=TRIM(file)
     CALL READ_MODEL(fullday(i),tmpItime,options,0)
  END DO


  IF (interpolate) THEN

     tmpItime=tmpItime+1
     file=fileday

     ! --------------------
     ! Read data from the NEXT day
     !
     IF (itime_end+1.GT.aux%ntlen) THEN
        tmpD=iday+1
        tmpM=month
        tmpY=year
        IF (tmpD.GT.DAYSINMONTH(year,month)) THEN
           tmpM=tmpM+1
           tmpD=1
        END IF
        IF (tmpM.GT.12) THEN 
           tmpY=tmpY+1
           tmpM=1
        END IF
        ! find the next file
        file = BUILD_FILENAME(rx,tmpY,tmpM,tmpD,model=options%model)

        IF (CHECK_FILE(TRIM(file))) THEN
           tmpItime=1 
        ELSE
           ! There is no next file. "You are probably at the end of your
           ! dataset", so choose the first time step of your current day
           file=fileDay
           tmpItime=itime_start
        END IF
     END IF
     fullday(i)%aux%netcdf_file=TRIM(file)
     ! --------------------

     PRINT 1,'next',TRIM(file),tmpItime

     CALL READ_MODEL(fullday(nsteps+1),tmpItime,options,0)
     ! ------------------------------
     ! SAVE output to be used in the next day
     previous=fullday(nsteps+1)
  END IF

  ! end gathered all the model data
  ! ==============================

  ! ==============================
  !
  ! Get the UTC corresponding to the EOT
  !
  !
  ! If the local time is 'EOT', what is that time in UTC at that longitude?
  T0 = MODULO(EOT - &
       24._wp/360._wp*aux%lon,&
       24._wp )

  IF (options%L2b%node.EQ."all") THEN
     PRINT *, "Simulating both asc and dec satellite overpasses"

     T0_b = MODULO( (EOT-12) - 24._wp/360._wp * aux%lon, 24._wp )
  END IF

  time_of_day = RESHAPE(T0,(/ ngrids /))

  !
  ! ============================== end EOT

  ! ====================
  ! Loop over the time steps
  ! conversion to avoid if UTC1+delta > 24
  ! One node at a time if you are doing both asc and dec
  !
  
  ! init the output data 
  CALL INITIALISE_MODEL_MATRIX(M,ngrids,nlev,mv=0._wp)

  UTC1 = MOD(aux%time(itime_start),24._wp) + aux%ref%hour

  DO i=1,nsteps ! one less than n*UTC interfaces used
     isL2b   = .FALSE.
     g_isL2b = .FALSE.
     c2u1    = 0._wp
     c2u2    = 0._wp
     c2u1_2  = 0._wp
     c2u2_2  = 0._wp

     t1 = MOD(UTC1,24._wp)
     t2 = MOD(UTC1,24._wp)+deltatime !utc2

     WRITE(0,'(2(a,1x,i2,1x))') &
          "Finding matching EOT data between UTC1=",&
          INT(t1),"and UTC2 =",INT(t2)

     ! -------------------- 
     ! find the longitudes where the EOT is contained between t1 and
     ! t2

     g_isL2b = (T0 .GE. t1) .AND. (T0 .LT. t2)
     isL2b= RESHAPE(g_isL2b,(/ ngrids /))
     isL2b_2 = SPREAD(isL2b,2,nlev)
     !--------------------

     ! ==============================
     ! Read the data
     ! ==============================

     IF (interpolate) THEN

        ! --------------------
        ! Find the weights between to 2 steps to match the EOT
        !
        g_c2u1=0
        g_c2u2=0
        WHERE(g_isL2b)
           g_c2u1  = 1-ABS(t1-T0)/deltaTime
           g_c2u2  = 1-ABS(t2-T0)/deltaTime
        END WHERE
        c2u1 = RESHAPE(g_c2u1, (/ ngrids /))
        c2u2 = RESHAPE(g_c2u2, (/ ngrids /))
     ELSE
        ! Emulate only using the time step at hand
        c2u1=1
        c2u2=0
     END IF
     c2u1_2=SPREAD(c2u1,2,nlev)
     c2u2_2=SPREAD(c2u2,2,nlev)
     ! --------------------

     !WRITE(rx,'(a,I0.2,".txt")') "testUTC=",INT(t1)
     !PRINT *,TRIM(rx)
     !OPEN(UNIT=8, FILE=TRIM(rx)) 
     !WRITE(0,'(a3,1x,5(a7,1x),a)') 'ind',' lon ',' T0 ','c2u1','c2u2','sum','isL2b'
     !DO ln=1,nlon
     !   WRITE(0,'(i3,1x,f7.3,1x,4(f7.4,1x),L)') &
     !        ln,aux%lon(ln,50),T0(ln,50),g_c2u1(ln,50),g_c2u2(ln,50),g_c2u1(ln,50)+g_c2u2(ln,50),g_isL2b(ln,50)
     !END DO
     !WRITE(0,'(a3,1x,5(a7,1x),a)') 'ind',' lon ',' T0 ','c2u1','c2u2','sum','isL2b'
     !CLOSE(8)

     ! ======================================
     ! Do the linear interpolation
     !
     !  Loop over longitudes and save data for the longitudes where UTC1
     !  and UTC2 encompass the sought local time. The data saved
     !  will be a linearly weighted averaged between the data
     !  at UTC1 and UTC2


     ! ------
     ! 2D fields
     WHERE (isL2b)
        ! sea ice
        M%CI    = M%CI    + c2u1*fullday(i)%CI + c2u2*fullday(i+1)%CI
        ! surface pressure
        M%PSURF = M%PSURF + EXP(c2u1*LOG(fullday(i)%PSURF) + c2u2*LOG(fullday(i+1)%PSURF))
        ! skin temperature
        M%SKT   = M%SKT   + c2u1*fullday(i)%SKT + c2u2*fullday(i+1)%SKT
        ! 2m temperature
        M%T2M   = M%T2M   + c2u1*fullday(i)%T2M + c2u2*fullday(i+1)%T2M
        ! Total cloud cover
        M%TCC   = M%TCC   + c2u1*fullday(i)%TCC + c2u2*fullday(i+1)%TCC
        ! Total column water vapor
        M%TCWV  = M%TCWV  + c2u1*fullday(i)%TCWV + c2u2*fullday(i+1)%TCWV

     END WHERE

     ! ------
     ! 3D fields
     WHERE (isL2b_2)
        ! layer cloud cover
        M%CC   = M%CC   + c2u1_2*fullday(i)%CC   + c2u2_2*fullday(i+1)%CC
        ! layer cloud liquid water content
        M%CLWC = M%CLWC + c2u1_2*fullday(i)%CLWC + c2u2_2*fullday(i+1)%CLWC
        ! layer cloud ice water content
        M%CIWC = M%CIWC + c2u1_2*fullday(i)%CIWC + c2u2_2*fullday(i+1)%CIWC
        ! layer humidity
        M%Q    = M%Q    + c2u1_2*fullday(i)%Q    + c2u2_2*fullday(i+1)%Q
        ! layer temperature
        M%T    = M%T    + c2u1_2*fullday(i)%T    + c2u2_2*fullday(i+1)%T

     END WHERE
     ! set the new time 
     UTC1=t2
  END DO

  IF (options%L2b%node.EQ."all") THEN
     !--------------------
     ! Do it again for the other node
     !
     UTC1 = MOD(aux%time(itime_start),24._wp) + aux%ref%hour
     DO i=1,nsteps
        g_isL2b=.FALSE.
        c2u1  =0._wp
        c2u2  =0._wp
        c2u1_2=0._wp
        c2u2_2=0._wp

        t1 = MOD(UTC1,24._wp)
        t2 = MOD(UTC1,24._wp)+deltatime !utc2
        
        WRITE(0,'(2(a,1x,i2,1x))') &
             "Finding matching EOT +12 data between UTC1 =",&
             INT(t1)," and UTC2 = ",INT(t2)

        ! --------------------
        ! find the longitudes where the EOT is contained between utc1 and utc2
        !
        g_isL2b = (T0_b .GE. t1) .AND. (T0_b .LT. t2)
        isL2b= RESHAPE(g_isL2b,(/ ngrids /))
        isL2b_2 = SPREAD(isL2b,2,nlev)
        !--------------------

        ! ==============================
        ! Read the data
        ! ==============================
        
        IF (interpolate) THEN

           ! --------------------
           ! Find the weights between to 2 steps to match the EOT
           !
           g_c2u1=0
           g_c2u2=0
           WHERE(g_isL2b)
              g_c2u1  = 1-ABS(t1-T0_b)/deltaTime
              g_c2u2  = 1-ABS(t2-T0_b)/deltaTime
           END WHERE
           c2u1 = RESHAPE(g_c2u1, (/ ngrids /))
           c2u2 = RESHAPE(g_c2u2, (/ ngrids /))
        ELSE
           ! Emulate only using the time step at hand
           c2u1=1
           c2u2=0
        END IF
        c2u1_2=SPREAD(c2u1,2,nlev)
        c2u2_2=SPREAD(c2u2,2,nlev)

        ! ======================================
        ! Do the linear interpolation

        ! ------
        ! 2D fields
        WHERE (isL2b)
           M%CI    = M%CI    + c2u1*fullday(i)%CI + c2u2*fullday(i+1)%CI
           M%PSURF = M%PSURF + EXP(c2u1*LOG(fullday(i)%PSURF) + c2u2*LOG(fullday(i+1)%PSURF))
           M%SKT   = M%SKT   + c2u1*fullday(i)%SKT + c2u2*fullday(i+1)%SKT
           M%T2M   = M%T2M   + c2u1*fullday(i)%T2M + c2u2*fullday(i+1)%T2M
           M%TCC   = M%TCC   + c2u1*fullday(i)%TCC + c2u2*fullday(i+1)%TCC
           M%TCWV  = M%TCWV  + c2u1*fullday(i)%TCWV + c2u2*fullday(i+1)%TCWV
        END WHERE

        ! ------
        ! 3D fields
        WHERE (isL2b_2)
           M%CC   = M%CC   + c2u1_2*fullday(i)%CC   + c2u2_2*fullday(i+1)%CC
           M%CLWC = M%CLWC + c2u1_2*fullday(i)%CLWC + c2u2_2*fullday(i+1)%CLWC
           M%CIWC = M%CIWC + c2u1_2*fullday(i)%CIWC + c2u2_2*fullday(i+1)%CIWC
           M%Q    = M%Q    + c2u1_2*fullday(i)%Q    + c2u2_2*fullday(i+1)%Q
           M%T    = M%T    + c2u1_2*fullday(i)%T    + c2u2_2*fullday(i+1)%T
        END WHERE


        ! set the new time 
        UTC1=t2
     END DO

     ! divide by 2 for two nodes since I just added data from this
     ! node to data from the previous node
     M%CI    = M%CI   /2._wp
     M%PSURF = M%PSURF/2._wp
     M%SKT   = M%SKT  /2._wp
     M%T2M   = M%T2M  /2._wp
     M%TCC   = M%TCC  /2._wp
     M%TCWV  = M%TCWV /2._wp
     M%CC    = M%CC   /2._wp
     M%CLWC  = M%CLWC /2._wp
     M%CIWC  = M%CIWC /2._wp
     M%Q     = M%Q    /2._wp
     M%T     = M%T    /2._wp
  END IF

  ! --------------------
  ! DEALLOCATING
  !
  DEALLOCATE(T0 ,&
       T0_b   ,&
       g_c2u1 ,&
       g_c2u2 ,&
       g_isL2b,&
       c2u1   ,&
       c2u2   ,&
       isL2b  ,&
       c2u1_2 ,&
       c2u2_2 ,&
       isL2b_2)
  
  DO i= 1,nsteps+1
      CALL DEALLOCATE_MODEL_MATRIX(fullday(i))
  END DO
  DEALLOCATE(fullday)
  !
  !--------------------

  PRINT *,"---------- sampled model data for one day"

END SUBROUTINE MAKE_LEVEL2B
