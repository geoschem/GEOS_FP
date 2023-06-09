!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosFpDriver0
!
! !DESCRIPTION: Program GeosFpDriver0 is a top-level driver for the 
!  GEOS-FP regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosFpDriver0
!
! !USES:
!
  USE GeosFpA1Module
  USE GeosFpCnModule
  USE GeosFpI3Module
  USE GeosFpInputsModule
  USE GeosFpRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosFpDriver1 creates the CN (constant), A1 (1hr time average), and
!  I3 (3hr instantaneous) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Now renamed Geos57 to GeosFP in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosFpInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosFpRegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL GeosFpMakeCn

  ! Create the 1-hour average data file
  CALL GeosFpMakeA1

  ! Create the 3-hour instantaneous files
  CALL GeosFpMakeI3

  ! Cleanup and quit 
  CALL GeosFpCleanup

END PROGRAM GeosFpDriver0
!EOP
