MODULE handy

  ! Several small handy functions
  !
  ! Salomon.Eliasson@smhi.se

  PUBLIC ::            &
       build_filename, &
       check_file,     &
       get_lun,        &
       replace_text,   &
       str2int,        &
       tally,          &
       time_keeper,    &
       to_upper

CONTAINS

  FUNCTION build_filename(formstr,&
       y,m,d,utc,&
       dir,model,sim,dataset,&
       sat,node,version,string,check) &
       RESULT(file)

    ! INPUT: 'formstr' should contaimn the fullpath regexp if dir is not provided

    IMPLICIT NONE

    CHARACTER(*), INTENT(in)           :: formstr 
    CHARACTER(*), INTENT(in), OPTIONAL :: dir,model,sim,dataset
    CHARACTER(*), INTENT(in), OPTIONAL :: sat,version
    CHARACTER(*),INTENT(in),OPTIONAL   :: node, string
    INTEGER, INTENT(in),OPTIONAL       :: y,m,d,utc ! year
    CHARACTER(LEN(formstr)+1000)       :: file
    LOGICAL, OPTIONAL :: check
    CHARACTER(4) :: s4
    CHARACTER(2) :: s2
    CHARACTER(1) :: s1
    LOGICAL   :: exists,kolla
    exists=.FALSE.

    IF (.NOT.PRESENT(check)) THEN
       kolla=.FALSE.
    ELSE
       kolla=check
    END IF
    ! --------
    ! FILE NAME
    !
    ! wildcards (all optional):
    ! #SIM='sim'                  
    ! #MODEL= options%model (A)   
    ! #Y4 = year  (i4)            
    ! #M2 = month (i2)            
    ! #D2 = day   (i2)            
    ! #Y2 = year  (i2)            
    ! #M1 = month (i1)            
    ! #D1 = day   (i1)            
    ! #UTC= utc   (i2)            
    ! #NODE = string              
    ! #VERSION = string           
    ! #STRING = string            
    !
    file=formstr
    IF (TALLY(file,'#SIM').GT.0) THEN
       file=Replace_Text (file,'#SIM',trim(sim))
    END IF
    IF (TALLY(file,'#DS').GT.0) THEN
       file=Replace_Text (file,'#DS',trim(dataset))
    END IF
    IF (TALLY(file,'#MODEL').GT.0) THEN
       file=Replace_Text (file,'#MODEL',trim(model))
    END IF
    IF (TALLY(file,'#Y4').GT.0) THEN
       WRITE(s4,'(i0.4)'),y
       file=Replace_Text (file,'#Y4',s4)
    ENDIF
    IF (TALLY(file,'#M2').GT.0) THEN
       WRITE(s2,'(i0.2)'),m
       file=Replace_Text (file,'#M2',s2)
    ENDIF
    IF (TALLY(file,'#D2').GT.0) THEN
       WRITE(s2,'(i0.2)'),d
       file=Replace_Text (file,'#D2',s2)
    ENDIF
    IF (TALLY(file,'#Y2').GT.0) THEN
       WRITE(s2,'(i0.2)'),y
       file=Replace_Text (file,'#Y2',s2)
    ENDIF
    IF (TALLY(file,'#M1').GT.0) THEN
       WRITE(s1,'(i1)'),m
       file=Replace_Text (file,'#M1',s1)
    ENDIF
    IF (TALLY(file,'#D1').GT.0) THEN
       WRITE(s1,'(i1)'),d
       file=Replace_Text (file,'#D1',s1)
    ENDIF
    IF (TALLY(file,'#UTC').GT.0) THEN
       WRITE(s2,'(i0.2)'),utc
       file=Replace_Text (file,'#D1',s2)
    ENDIF
    IF (TALLY(file,'#SAT').GT.0) THEN
       file=Replace_Text (file,'#SAT',TRIM(sat))
    ENDIF
    IF (TALLY(file,'#NODE').GT.0) THEN
       IF (node.EQ.'all') THEN
          file=Replace_Text (file,'_#NODE','')
       ELSE
          file=Replace_Text (file,'#NODE',TRIM(node))
       END IF
    ENDIF
    IF (TALLY(file,'#VERSION').GT.0) THEN
       file=Replace_Text (file,'#VERSION',TRIM(version))
    ENDIF
    IF (TALLY(file,'#STRING').GT.0) THEN
       file=Replace_Text (file,'#STRING',TRIM(string))
    ENDIF
    IF (PRESENT(dir)) THEN
       file=TRIM(dir)//TRIM(file)
    END IF
    IF (kolla) THEN
       IF (.NOT.CHECK_FILE(TRIM(file))) THEN
          STOP "stopped in build_filename() in model_input.F90"
       END IF
       PRINT '("File:",A)',TRIM(file)
    END IF
  END FUNCTION build_filename

  ELEMENTAL SUBROUTINE str2int(str,int,stat)
    ! convert string to integer
    !
    ! Code by Alexander Vogt found on
    ! http://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90

    IMPLICIT NONE
    ! Arguments
    CHARACTER(len=*),INTENT(in) :: str
    INTEGER,INTENT(out)         :: int
    INTEGER,INTENT(out)         :: stat

    READ(str,*,iostat=stat)  int
  END SUBROUTINE str2int

  ! ------

  FUNCTION check_file(file,directory) RESULT(exists)
    ! Check to see if the file or directory exists

    IMPLICIT NONE

    CHARACTER(*), OPTIONAL, INTENT(in) :: file
    CHARACTER(*), OPTIONAL, INTENT(in) :: directory
    LOGICAL                            :: exists

    exists = .TRUE.
    IF (PRESENT(directory)) INQUIRE(FILE=TRIM(directory), EXIST=exists)
    IF (PRESENT(file)) INQUIRE(FILE=TRIM(file), EXIST=exists)
    IF (.NOT.exists) THEN
       IF (PRESENT(file)) PRINT '(A,A,A)',"file '",TRIM(file),"' doesn't exist"
       IF (PRESENT(directory)) PRINT '(A,A,A)',"directory '",TRIM(directory),"' doesn't exist"
    END IF

  END FUNCTION check_file

  FUNCTION get_lun (lu_max)
    !
    !   get_file_unit returns a unit number that is not in use
    INTEGER, INTENT(in) :: lu_max

    ! out
    INTEGER :: get_lun

    ! INTERNAL
    INTEGER :: lu, m, iostat
    LOGICAL :: opened
    !
    m = lu_max  ;  IF (m < 1) m = 97
    DO lu = m,1,-1
       INQUIRE (unit=lu, opened=opened, iostat=iostat)
       IF (iostat.NE.0) CYCLE
       IF (.NOT.opened) EXIT
    END DO
    !
    get_lun = lu

  END FUNCTION get_lun

  FUNCTION to_upper(strIn) RESULT(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: strIn
    CHARACTER(len=LEN(strIn)) :: strOut
    INTEGER :: i,j

    DO i = 1, LEN(strIn)
       j = IACHAR(strIn(i:i))
       IF (j>= IACHAR("a") .AND. j<=IACHAR("z") ) THEN
          strOut(i:i) = ACHAR(IACHAR(strIn(i:i))-32)
       ELSE
          strOut(i:i) = strIn(i:i)
       END IF
    END DO

  END FUNCTION to_upper

  SUBROUTINE time_keeper(elapsed)
    ! print time keeping for poor-mans profiling

    IMPLICIT NONE

    REAL(8), INTENT(in) :: elapsed
    REAL(8) :: minute, hour, second

    minute                   = 0.
    hour                     = 0.
    second                   = 0.

    hour=FLOOR(REAL(elapsed)/3600);
    minute=FLOOR(60*(REAL(elapsed)/3600-hour));
    second=REAL(elapsed)-60*minute-3600*hour;
    IF (second>=60) THEN
       minute = minute+1;
       second = second-60;
    END IF
    IF (minute>=60) THEN
       hour = hour+1;
       minute = minute-60;
    END IF
    PRINT '(A,F3.0,1x,"h",1x,F3.0,1x,"min",1x,F6.3,1x,"s")',&
         "Elapsed time = ",hour,minute,second

  END SUBROUTINE time_keeper

  !
  ! String functions
  !
  ! Replace_Text   in all occurances in string with replacement string
  ! Tally          occurances in string of text arg

  FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
    ! by David Frank  dave_frank@hotmail.com
    ! http://home.earthlink.net/~dave_gemini/strings.f90
    ! (http://fortranwiki.org/fortran/show/String_Functions)
    ! copyright:
    ! All text on the Fortran Wiki is available under the terms of the
    ! GNU Free Documentation License (GFDL) with no invariant sections,
    ! front-cover texts, or back-cover texts.  $ $If you contribute to the
    ! Fortran Wiki, you thereby agree to license the contributed material
    ! to the public under the GFDL (version 1.2 or any later version
    ! published by the Free Software Foundation, with no invariant
    ! sections, front-cover texts, or back-cover texts).

    CHARACTER(*),INTENT(in)  :: s,text,rep
    CHARACTER(LEN(s)+100)    :: outs     ! provide outs with extra 100 char len
    INTEGER                  :: i, nt, nr

    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
    DO
       i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
       outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    END DO
  END FUNCTION Replace_Text

  INTEGER FUNCTION Tally (s,text)
    CHARACTER(*) :: s, text
    INTEGER :: i, nt

    Tally = 0 ; nt = LEN_TRIM(text)
    DO i = 1,LEN(s)-nt+1
       IF (s(i:i+nt-1) == text(:nt)) Tally = Tally+1
    END DO
  END FUNCTION Tally

END MODULE handy
