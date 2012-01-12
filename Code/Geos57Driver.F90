!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver
!
! !DESCRIPTION: Program Geos57Driver is the top-level driver for the 
!  GEOS-5.7.x regridding programs.  Geos57Driver will call routines to 
!  extract, regrid, and save the GEOS-5.7.x met data to files for 
!  input to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver
!
! !USES:
!
  USE Geos57A1Module
  USE Geos57A3CldModule
  USE Geos57A3DynModule
  USE Geos57A3MstCModule
  USE Geos57A3MstEModule
  USE Geos57CnModule
  USE Geos57I3Module
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  26 Oct 2011 - R. Yantosca - Initial Version, based on MerraDriver
!  03 Jan 2012 - R. Yantosca - Activate calls to Geos57MakeCn, Geos57MakeI3
!  12 Jan 2012 - R. Yantosca - Activate call to Geos57MakeA3* routines
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL Geos57MakeCn

  ! Create the 6-hour instantaneous data file
  CALL Geos57MakeI3

  ! Create the 1-hour average data file
  CALL Geos57MakeA1

  ! Create the 3-hour average data files
  CALL Geos57MakeA3Cld
  CALL Geos57MakeA3Dyn
  CALL Geos57MakeA3MstC
  CALL Geos57MakeA3MstE

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver
!EOP
