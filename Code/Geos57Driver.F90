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
  USE Geos57A3Module
  USE Geos57CnModule
  USE Geos57I6Module
  USE Geos57InputsModule
  USE Geos57RegridModule

  IMPLICIT NONE
!
! !REVISION HISTORY: 
!  30 Aug 2011 - R. Yantosca - Initial Version, based on MerraDriver
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Geos57Initialize

  ! Initialize GEOS-5 regridding code
  CALL Geos57RegridInit

  ! Create the constant data file
  IF ( doMakeCn ) CALL Geos57MakeCn

  ! Create the 1-hour average data file
  CALL Geos57MakeA1

  ! Create the 3-hour average data file
  CALL Geos57MakeA3

  ! Create the 6-hour instantaneous data file
  CALL Geos57MakeI6

  ! Cleanup and quit 
  CALL Geos57Cleanup

END PROGRAM Geos57Driver
!EOP
