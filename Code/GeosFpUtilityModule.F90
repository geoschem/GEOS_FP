!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpUtilityModule
!
! !DESCRIPTION: Module GeosFpUtilModule contains several utility routines
!  for the GEOS-FP regridding package.
!\\
!\\
! !INTERFACE: 
!
MODULE GeosFpUtilityModule
! 
! !USES:
!
  USE CharpakModule

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GetNFields
  PUBLIC :: NotDir
  PUBLIC :: SystemDateTime
  PUBLIC :: SystemTimeStamp
  PUBLIC :: TimeStampString
  PUBLIC :: UnitsForTime
!
! !REVISION HISTORY:
!  04 Jan 2012 - R. Yantosca - Initial version

!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetNFields
!
! !DESCRIPTION: Returns the list of fields and number of fields to regrid
!  for each GEOS-FP raw data file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetNFields( dataList, nFields, fields )
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN)  :: dataList   ! Comma-sep'd field name list
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: nFields    ! Number of fields
    CHARACTER(LEN=*), INTENT(OUT) :: fields(:)  ! Array of field names
! 
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: F

    ! Split the data field list by commas into an array
    CALL makeCharArrayFromCharList( dataList, ',', fields )
    
    ! Compute the number of data fields we will process
    nFields = 0

    DO F = 1, SIZE( fields )
       IF ( TRIM( fields(F) ) /= ''      .and. &
            TRIM( fields(F) ) /= 'none' ) THEN
          nFields = nFields + 1
       ENDIF
    ENDDO

  END SUBROUTINE GetNFields
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NotDir
!
! !DESCRIPTION: 
!  Returns the file name part of a full directory path.  Akin to the GNU
!  Make "notdir" function.For each GEOS-FP raw data file 
!\\
!\\
! !INTERFACE:
!
  FUNCTION NotDir( pathName ) RESULT( fileName )
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*),   INTENT(IN) :: pathName   ! Directory path
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)             :: fileName   ! Just the filename part
! 
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  20 Sep 2013 - R. Yantosca - Add NotDir function
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: F
    INTEGER            :: nFields
    CHARACTER(LEN=255) :: fields(50)  ! Allow for 50 subdirs

    ! Split the data field list by commas into an array
    CALL makeCharArrayFromCharList( pathName, '/', fields )
    
    ! Find the last occupied field
    DO F = 1, SIZE( fields )
       IF ( TRIM( fields(F) ) /= '' ) nFields = nFields + 1
    ENDDO

    ! File name
    fileName = TRIM ( fields(nFields ) )

  END FUNCTION NotDir
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SystemDateTime
!
! !DESCRIPTION: Subroutine SystemDateTime returns the onboard system date and 
!  time.  If the optional timeZone argument is passed, then SystemDateTime
!  will return the GMT date and time.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SystemDateTime( sysDate, sysTime, timeZone )
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)           :: sysDate   ! Sys date (YYYYMMDD
    INTEGER,          INTENT(OUT)           :: sysTime   ! Sys time (hhmmss)
    CHARACTER(LEN=5), INTENT(OUT), OPTIONAL :: timeZone  ! Time zone offset
!
! !REMARKS:
!  Uses the F90 intrinsic function DATE_AND_TIME.  The VALUES argument of 
!  DATE_AND_TIME returns the following quantities as a 1-D vector:  
!  (/YYYY, MM, DD, GMT_MIN, HH, MM, SS, MSEC/)
!
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Local variables
    INTEGER              :: V(8)
    CHARACTER(LEN=5)     :: Z
    CHARACTER(LEN=8)     :: D
    CHARACTER(LEN=10)    :: T

    !=================================================================
    ! SystemDateTime begins here!
    !=================================================================

    ! Initialize
    D = 'ccyymmdd'
    T = 'hhmmss.sss'

    ! Call F90 intrinsic function DATE_AND_TIME to get the onboard time
    IF ( PRESENT( timeZone ) ) THEN
       CALL Date_And_Time( DATE=D, TIME=T, ZONE=Z, VALUES=V )  ! GMT time
    ELSE
       CALL Date_And_Time( DATE=D, TIME=T,         VALUES=V )  ! Local std time
    ENDIF
     
    ! Save to YYYYMMDD and HHMMSS format
    sysDate = ( V(1) * 10000 ) + ( V(2) * 100 ) + V(3) 
    sysTime = ( V(5) * 10000 ) + ( V(6) * 100 ) + V(7)
    
    ! Return optional time zone
    IF ( PRESENT( timeZone ) ) timeZone = Z

  END SUBROUTINE SystemDateTime
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SystemTimeStamp
!
! !DESCRIPTION: Function SystemTimeStamp returns a string with 
!  the onboard system GMT time in YYYY/MM/DD HH:MM format.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SystemTimeStamp() RESULT( timeStamp )
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255) :: timeStamp
! 
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER          :: sysDate, sysTime
    CHARACTER(LEN=5) :: timeZone

    !=================================================================
    ! SystemTimeStamp begins here!
    !=================================================================

    ! Get system date and time
    CALL SystemDateTime( sysDate, sysTime, timeZone )

    ! Create a string w/ system date & time
    timeStamp = TimeStampString( sysDate, sysTime, timeZone )

  END FUNCTION SystemTimeStamp
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TimeStampString
!
! !DESCRIPTION: Function TimeStampString returns a formatted string 
!  "YYYY/MM/DD hh:mm" for the a date and time specified by YYYYMMDD and 
!  hhmmss.
!\\
!\\
! !INTERFACE:
!
  FUNCTION TimeStampString( yyyymmdd, hhmmss, timeZone ) RESULT( timeStamp )
!
! !INPUT PARAMETERS: 
!
    INTEGER,          INTENT(IN)           :: yyyymmdd   ! YYYY/MM/DD date
    INTEGER,          INTENT(IN)           :: hhmmss     ! hh:mm:ss time
    CHARACTER(LEN=5), INTENT(IN), OPTIONAL :: timeZone   ! Time zone offset
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)                     :: timeStamp
!
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on GEOS-Chem
!  05 Jan 2012 - R. Yantosca - Print time zone offset (e.g. GMT-0500) if
!                              the timeZone argument is passed
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: year, month, day, hour, minute, second

    !=================================================================
    ! TimeStampString begins here!
    !=================================================================

    ! Split up date and time into individual variables
    CALL YMD_Extract( yyyymmdd, year, month,  day    )
    CALL YMD_Extract( hhmmss,   hour, minute, second )

    ! Use FORTRAN internal write to create the string
    IF ( PRESENT( timeZone ) ) THEN
       WRITE( timeStamp, 10 ) year, month, day, hour, minute, second, timeZone
    ELSE
       WRITE( timeStamp, 20 ) year, month, day, hour, minute, second
    ENDIF

    ! Format statements
 10 FORMAT( i4.4, '/', i2.2, '/', i2.2, ' ',        &
            i2.2, ':', i2.2, ':', i2.2, ' GMT', a5 )
 20 FORMAT( i4.4, '/', i2.2, '/', i2.2, ' ',        &
            i2.2, ':', i2.2, ':', i2.2             )

  END FUNCTION TimeStampString
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: UnitsForTime
!
! !DESCRIPTION: Function UnitsForTime returns the units string 
!  for the "time" index array.  The string takes the form of
!  "minutes since YYYY-MM-DD hh:mm:ss.0".
!\\
!\\
! !INTERFACE:
!
  FUNCTION UnitsForTime( yyyymmdd ) RESULT( units )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: yyyymmdd   ! YYYY/MM/DD date
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)  :: units
!
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: year, month, day

    !=================================================================
    ! UnitsForTime begins here!
    !=================================================================

    ! Split up date and time into individual variables
    CALL YMD_EXTRACT( yyyymmdd, year, month, day )
      
    ! Use FORTRAN internal write to create the string
    WRITE( units, 10 ) year, month, day

    ! Format statement
 10 FORMAT( 'minutes since ', i4.4, '-', i2.2, '-', i2.2, ' 00:00:00.0' )

  END FUNCTION UnitsForTime
!EOC
END MODULE GeosFpUtilityModule
