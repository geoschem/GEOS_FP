!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver4
!
! !DESCRIPTION: Program Geos57Driver4 is a top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver4
!
! !USES:
!
  USE Geos57A3MstEModule
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Geos57Driver4 creates the A3mstE (3hr time-average moist parameters, on
!  level edges) data files for input into GEOS-Chem.
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

  ! Create the 3-hour average data files
  CALL Geos57MakeA3MstE

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver4
!EOP
