!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosFpDriver2
!
! !DESCRIPTION: Program GeosFpDriver2 is the top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosFpDriver2
!
! !USES:
!
  USE GeosFpA3MstCModule
  USE GeosFpA3MstEModule
  USE GeosFpInputsModule
  USE GeosFpRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  GeosFpDriver1 creates the A3mstC (3hr time-averaged moist parameters, on
!  level centers) and A3MstE (3hr time-averaged moist parameters on level
!  edges) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  20 Sep 2013 - R. Yantosca - Rename Geos57 to GeosFp in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosFpInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosFpRegridInit

  ! Create the 3-hour average data file
  CALL GeosFpMakeA3MstC
  CALL GeosFpMakeA3MstE

  ! Cleanup and quit 
  CALL GeosFpCleanup

END PROGRAM GeosFpDriver2
!EOP
