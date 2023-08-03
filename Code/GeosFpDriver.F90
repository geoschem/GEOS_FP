!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GeosFpDriver
!
! !DESCRIPTION: Program GeosFpDriver is the top-level driver for the 
!  GEOS-FP regridding programs.  GeosFpDriver will call routines to 
!  extract, regrid, and save the GEOS-FP met data to files for 
!  input to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
PROGRAM GeosFpDriver
!
! !USES:
!
  USE GeosFpA1Module
  USE GeosFpA3CldModule
  USE GeosFpA3DynModule
  USE GeosFpA3MstCModule
  USE GeosFpA3MstEModule
  USE GeosFpCnModule
  USE GeosFpI3Module
  USE GeosFpInputsModule
  USE GeosFpRegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Renamed Geos57 to GeosFp in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL GeosFpInitialize

  ! Initialize GEOS-5 regridding code
  CALL GeosFpRegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL GeosFpMakeCn

  ! Create the 1-hour average data file
  CALL GeosFpMakeA1

  ! Create the 3-hour average data files
  CALL GeosFpMakeA3Cld
  CALL GeosFpMakeA3Dyn
  CALL GeosFpMakeA3MstC
  CALL GeosFpMakeA3MstE

  ! Create the 6-hour instantaneous data file
  CALL GeosFpMakeI3

  ! Cleanup and quit 
  CALL GeosFpCleanup

END PROGRAM GeosFpDriver
!EOP
