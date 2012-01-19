!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver5
!
! !DESCRIPTION: Program Geos57Driver5 is a top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver
!
! !USES:
!
  USE Geos57I3Module
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Geos57Driver5 creates the I3 (3hr instantaneous) data files for input 
!  into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  19 Jan 2012 - R. Yantosca - Initial version, created from Geos57Driver
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the 3-hour instantaneous data file
  CALL Geos57MakeI3

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver
!EOP
