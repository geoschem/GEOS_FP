!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver2
!
! !DESCRIPTION: Program Geos57Driver2 is the top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver2
!
! !USES:
!
  USE Geos57A3MstCModule
  USE Geos57A3MstEModule
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Geos57Driver1 creates the A3mstC (3hr time-averaged moist parameters, on
!  level centers) and A3MstE (3hr time-averaged moist parameters on level
!  edges) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  19 Jan 2012 - R. Yantosca - Initial version, created from Geos57Driver
!  23 Jan 2012 - R. Yantosca - Now create the A3MstC and A3MstE files
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the 3-hour average data file
  CALL Geos57MakeA3MstC
  CALL Geos57MakeA3MstE

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver2
!EOP
