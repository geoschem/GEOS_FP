!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver1
!
! !DESCRIPTION: Program Geos57Driver1 is a top-level driver for the 
!  GEOS-5.7.x regridding programs. 
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver1
!
! !USES:
!
  USE Geos57A3CldModule
  USE Geos57A3DynModule
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Geos57Driver1 creates the A3cld (3hr time-averaged cloud parameters) and
!  A3dyn  (3hr time-averaged dynamics parameters) data files for input into 
!  GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  19 Jan 2012 - R. Yantosca - Initial version, based on Geos57Driver
!  23 Jan 2012 - R. Yantosca - Now create the A3dyn file

!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the 3-hour average data files
  CALL Geos57MakeA3Cld
  CALL Geos57MakeA3Dyn

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver1
!EOP
