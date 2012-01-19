!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Geos57Driver0
!
! !DESCRIPTION: Program Geos57Driver0 is a top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Geos57Driver0
!
! !USES:
!
  USE Geos57A1Module
  USE Geos57CnModule
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Geos57Driver1 creates the CN (constant) and A1 (1hr time average)
!  data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  19 Jan 2012 - R. Yantosca - Initial version, created from Geos57Driver
!EOP
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL Geos57MakeCn

  ! Create the 1-hour average data file
  CALL Geos57MakeA1

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver0
!EOP
