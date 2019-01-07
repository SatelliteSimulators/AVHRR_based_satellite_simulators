MODULE my_maths
  ! This module contains subroutine that:
  !     sort an array
  !     finds the medium in an array
  !
  !
  ! Salomon.Eliasson@smhi.se

  IMPLICIT NONE

  PUBLIC :: get_median,    &
       get_sorted,         &
       daysinmonth,        &
       isleapyear,         &
       dayofyear,          &
       day_since_year,     &
       day_since_day

CONTAINS

  SUBROUTINE get_sorted(DATA,len,sorted)
    ! Function that sorts a data vector and finds the median value. If the
    ! length of the data vector is even then median is 
    ! ( sorted(len/2) +sorted(len/2+1) )/2.0

    implicit none

    ! IN
    integer, intent(in) :: len
    real, intent(in),dimension(len) :: data

    ! OUT
    real, intent(out), dimension(len) :: sorted

    ! internal
    integer, dimension(len) :: indi
    logical, dimension(len) :: not_sorted
    integer :: ii, k
    real :: mini

    ! first sort the data into ascending order
    !
    ! 1) Loop over every element and flag the elements that are already
    !    copied to the sorted array (especially if there are non unique values)
    ! 2) Once the lowest value is found, Restart. keep going until
    !    everything is sorted


    indi(:) = 0
    mini = MINVAL(data)
    not_sorted(:) = .true. 
    k = 1
    do while (k .le. len) 
       do ii = 1, len
          if ( (data(ii) .eq. mini) .and. not_sorted(ii) ) then

             indi(k) = ii
             not_sorted(ii) = .false.
             mini = MINVAL(data,not_sorted)
             k = k+1
             exit
          end if
       end do
    end do

    do ii = 1,len
       sorted(ii) = data(indi(ii))
    end do

  END SUBROUTINE get_sorted

  !  ----------------------------------------

  SUBROUTINE get_median(DATA,len,median)
    ! Function that sorts a data vector and finds the median value. If the
    ! length of the data vector is even then median is 
    ! ( sorted(len/2) +sorted(len/2+1) )/2.0

    implicit none

    ! IN
    integer, intent(in) :: len
    real, intent(in),dimension(len) :: data

    ! OUT
    real, intent(out) :: median

    ! INTERNAL
    real, dimension(len) :: sorted

    call get_sorted(data,len,sorted)

    ! Find the median
    if (mod(len,2) .eq. 1) then
       median = sorted((len-1)/2 +1)
    else
       median = (sorted(len/2)  + sorted(len/2+1))/ 2.0
    end if

  END SUBROUTINE get_median

  ELEMENTAL FUNCTION daysInMonth(year, month)
    !DAYSINMONTH Number of days in a month.
    !
    !   DAYSINMONTH(YEAR, MONTH) returns the number of days in the given month.
    !
    !   Gregorian calendar is assumed.
    !
    !   Author:      Peter J. Acklam
    !   Time-stamp:  2002-03-03 12:52:00 +0100
    !   E-mail:      pjacklam@online.no
    !   URL:         http://home.online.no/~pjacklam
    !
    ! Made into fortran by Salomon Eliasson


    IMPLICIT NONE

    ! In 
    INTEGER, INTENT(in) :: year, month

    ! Out
    INTEGER             :: daysInMonth

    ! internal
    INTEGER :: nDays(12)

    ! Now get the number of days in the month.
    nDays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    daysInMonth  = nDays(month)

    ! Add leap day as necessary.
    IF (month .EQ. 2) daysInMonth=daysInMonth+ISLEAPYEAR(year)

  END FUNCTION daysinmonth

  ELEMENTAL FUNCTION isleapyear(year)
    !ISLEAPYEAR True for leap years.
    !
    !ISLEAPYEAR(YEAR) returns .true. if YEAR is a leap
    !year and .false. if it is not. Gregorian calendar is assumed.
    !
    !A year is a leap year if the following returns true
    !
    !    ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400)
    !
    !A year is not a leap year if the following returns true
    !
    !   rem(year, 4) | ( ~rem(year, 100) & rem(year, 400) )
    !
    !Author:      Peter J. Acklam
    !Time-stamp:  2002-03-03 12:51:45 +0100
    !E-mail:      pjacklam@online.no
    !URL:         http://home.online.no/~pjacklam
    !
    ! Turned into fortran by Salomon Eliasson
    implicit none

    ! In
    INTEGER, INTENT(in) :: year

    ! OUT
    INTEGER :: isleapyear

    ! out
    LOGICAL :: ily ! is leap year

    ily = ( (year/4 .EQ. REAL(year)/REAL(4)) .AND. ( (year/100) .NE. REAL(year)/REAL(100) ) )&
         .OR. ( (year/4 .EQ. REAL(year)/REAL(4)) )    

    IF (ily) THEN
       isleapyear = 1
    ELSE
       isleapyear = 0
    END IF

  END FUNCTION isleapyear

  ELEMENTAL FUNCTION dayofyear(year,month,day,hour,minute,second)
    !DAYOFYEAR Ordinal number of day in a year.
    !
    !   DAYOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
    !   day number in the given year plus a fractional part depending on the
    !   time of day.
    !
    !   Any missing MONTH or DAY will be replaced by 1.  HOUR, MINUTE or SECOND
    !   will be replaced by zeros.
    !
    !   If no date is specified, the current date and time is used.  Gregorian
    !   calendar is assumed.
    !
    !   Author:      Peter J. Acklam
    !   Time-stamp:  2002-03-03 12:52:04 +0100
    !   E-mail:      pjacklam@online.no
    !   URL:         http://home.online.no/~pjacklam
    !*
    ! Turned into fortran by Salomon Eliasson


    IMPLICIT NONE

    INTEGER, INTENT(in)           :: year, month, day
    INTEGER, OPTIONAL, INTENT(in) :: hour, minute, second

    REAL                          :: dayofyear

    INTEGER                       :: days_in_prev_months(12)
    INTEGER                       :: h,m,s

    days_in_prev_months = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /);

    dayofyear = days_in_prev_months(month)  ! days in prev. months

    ! day in month
    dayofyear = dayofyear +day 

    ! leap day
    if (month .gt. 2) dayofyear=dayofyear+ISLEAPYEAR(year)

    ! part of day
    if (.not.present(hour)) then 
       h   = 0
    else
       h = hour
    end if
    if (.not.present(minute)) then 
       m = 0
    else
       m = minute
    end if
    if (.not.present(second)) then
       s = 0
    else
       s = second
    end if
    dayofyear = dayofyear + ( s + 60.0*m + 3600.0*h )/86400.0;  

  END FUNCTION dayofyear

  ELEMENTAL FUNCTION day_since_year(year,iyr,imon,iday) RESULT(day)

    INTEGER, INTENT(in) :: year
    INTEGER, INTENT(in) :: iyr,imon,iday

    INTEGER :: y
    INTEGER :: day

    day = 0
    y = year

    ! cycle years
    DO WHILE (y .LT. iyr)
       day = day + 365 + ISLEAPYEAR(y)
       y  = y+1
    END DO

    !add day of year (1st day of the month)
    !e  add offset: 19700101 is 0 days since 19700101
    day = day + int(DAYOFYEAR(iyr,imon,iday)) - 1

  END FUNCTION day_since_year

  ELEMENTAL FUNCTION day_since_day(ref_year,ref_month,ref_day,&
       iyr,imon,iday) RESULT(day)
    ! counts the number of days since a reference date. This works as long
    ! as both dates are ahead of 500BC
    !
    !

    INTEGER, INTENT(in) :: ref_year,ref_month,ref_day
    INTEGER, INTENT(in) :: iyr,imon,iday

    INTEGER :: day

    ! years (from 500BC)
    day=DAY_SINCE_YEAR(-500,iyr,imon,iday)&                !day compared to 500BC
         - DAY_SINCE_YEAR(-500,ref_year,ref_month,ref_day) !ref compared to 500BC


  END FUNCTION day_since_day
END MODULE my_maths
