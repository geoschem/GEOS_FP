!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosFpDriver1
!
! !DESCRIPTION: Program GeosFpDriver1 is a top-level driver for the 
!  GEOS-5.7.x regridding programs. 
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosFpDriver1
!
! !USES:
!
  USE GeosFpA3CldModule
  USE GeosFpA3DynModule
  USE GeosFpInputsModule
  USE GeosFpRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosFpDriver1 creates the A3cld (3hr time-averaged cloud parameters) and
!  A3dyn  (3hr time-averaged dynamics parameters) data files for input into 
!  GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Now renamed Geos57 to GeosFp in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosFpInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosFpRegridInit

  ! Create the 3-hour average data files
  CALL GeosFpMakeA3Cld
  CALL GeosFpMakeA3Dyn

  ! Cleanup and quit 
  CALL GeosFpCleanup

END PROGRAM GeosFpDriver1
!EOP
