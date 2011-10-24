!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: MerraBinaryModule
!
! !DESCRIPTION: MerraBinaryModule contains routines to write to GEOS-Chem
!  binary output format.
!\\
!\\
! !INTERFACE: 

MODULE MerraBinaryModule
! 
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: OpenBinary
  PUBLIC  :: WriteBinary

  INTERFACE  WriteBinary
     MODULE PROCEDURE WriteBinary2D
     MODULE PROCEDURE WriteBinary3D
  END INTERFACE  
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: WriteBinary2D
  PRIVATE :: WriteBinary3D
!
! !REVISION HISTORY:
!  23 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
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
! !IROUTINE: OpenBinary
!
! !DESCRIPTION: Subroutine OpenBinary opens the binary file into which the 
!  regridded met field data will be written to.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OpenBinary( unit, ident, fileName )
!
! !INPUT PARAMETERS: 
!
    INTEGER,          INTENT(IN) :: unit       ! Logical unit #
    CHARACTER(LEN=8), INTENT(IN) :: ident      ! top-of-file ID str
    CHARACTER(LEN=*), INTENT(IN) :: fileName   ! File name
!
! !REVISION HISTORY: 
!  23 Jul 2010 - R. Yantosca - Initial Version, based on GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC  
    ! Open the file for writing
    OPEN( unit,               FILE=TRIM( fileName ), &
          FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )

    ! Place the file ident string at the top of the f`ile
    ! This helps GAMAP identify the modeltype and resolution
    WRITE( unit ) ident

  END SUBROUTINE openBinary
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WriteBinary2D
!
! !DESCRIPTION: Subroutine WriteBinary2D writes a 2-D data block to the 
!  GEOS-Chem binary format met field file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WriteBinary2D( unit, name, yyyymmdd, hhmmss, data, IX, JX )
!
! !INPUT ARGUMETNTS:
!
    INTEGER,          INTENT(IN) :: unit            ! Logical unit #
    INTEGER,          INTENT(IN) :: yyyymmdd        ! Date in YYYY/MM/DD form
    INTEGER,          INTENT(IN) :: hhmmss          ! Time in hh:mm:ss form
    INTEGER,          INTENT(IN) :: IX              ! Lon dimension of DATA
    INTEGER,          INTENT(IN) :: JX              ! Lat dimension of DATA
    REAL*4,           INTENT(IN) :: data(IX,JX)     ! Array w/ data to write 
    CHARACTER(LEN=8), INTENT(IN) :: name            ! Name of the data field
!
! !REVISION HISTORY: 
!  23 Jul 2010 - R. Yantosca - Initial Version, based on GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC  
    ! Write field name 
    WRITE( unit ) name

    ! Write YYYYMMDD, HHMMSS, data
    WRITE( unit ) yyyymmdd, hhmmss, data

  END SUBROUTINE WriteBinary2D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WriteBinary3D
!
! !DESCRIPTION: Subroutine WriteBinary3D writes a 3-D data block to the 
!  GEOS-Chem binary format met field file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WriteBinary3D( unit, name, yyyymmdd, hhmmss, data, IX, JX, LX )
!
! !INPUT ARGUMETNTS:
!
    INTEGER,          INTENT(IN) :: unit            ! Logical unit #
    INTEGER,          INTENT(IN) :: yyyymmdd        ! Date in YYYY/MM/DD form
    INTEGER,          INTENT(IN) :: hhmmss          ! Time in hh:mm:ss form
    INTEGER,          INTENT(IN) :: IX              ! Lon dimension of DATA
    INTEGER,          INTENT(IN) :: JX              ! Lat dimension of DATA
    INTEGER,          INTENT(IN) :: LX              ! Alt dimension of DATA
    REAL*4,           INTENT(IN) :: data(IX,JX,LX)  ! Array w/ data to write 
    CHARACTER(LEN=8), INTENT(IN) :: name            ! Name of the data field
!
! !REVISION HISTORY: 
!  23 Jul 2010 - R. Yantosca - Initial Version, based on GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Write field name 
    WRITE( unit ) name

    ! Write YYYYMMDD, HHMMSS, data
    WRITE( unit ) yyyymmdd, hhmmss, data

  END SUBROUTINE WriteBinary3D
!EOC
END MODULE MerraBinaryModule

