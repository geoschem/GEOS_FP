!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: MerraCnModule
!
! !DESCRIPTION: Module MerraCnModule contains routines to create the 
!  GEOS-Chem average constant data files from the MERRA raw data.
!\\
!\\
! !INTERFACE: 

MODULE MerraCnModule
! 
! !USES:
!
  ! MERRA data modules
  USE CharpakModule
  USE MerraInputsModule
  USE MerraBinaryModule
  USE MerraRegridModule

  ! HDF-EOS modules
  USE He4IncludeModule
  USE He4GridModule
  USE He4ErrorModule

  IMPLICIT NONE
  PRIVATE

# include "He4Define.h"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: MerraMakeCn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GetNFields
  PRIVATE :: ProcessCn2dAsmNx
  PRIVATE :: OpenCnBinary
!
! !REVISION HISTORY:
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  11 Aug 2010 - R. Yantosca - Now get constant data from const_2d_chm_Nx file
!  12 Aug 2010 - R. Yantosca - Renamed ProcessCn2dChmNx to ProcessCn2dAsmNx
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MerraMakeCnFiles 
!
! !DESCRIPTION: Routine MerraMakeConstFiles is the the driver routine for 
! \begin{enumerate}
! \item Extracting constant data fields (surface values) from 
!       the MERRA raw data files (HDF4-EOS format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data in a format that GEOS-Chem can read.
! \end{enumerate}
! This routine is called directly from the main program MerraDriver.F90
!\\
!\\
! !INTERFACE:
  SUBROUTINE MerraMakeCn
!
! !REVISION HISTORY: 
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  11 Aug 2010 - R. Yantosca - Now get constant data from const_2d_chm_Nx file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: F
    INTEGER                 :: nFields
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: msg

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE MerraMakeConst %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Return the list of fields and number of fields to process
    ! from each of the MERRA raw met data files
    CALL GetNFields( const_2d_asm_Nx_data, nFields, fields )

    ! Total number of fields that we will process
    nAllFields = nFields

    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( const_2d_asm_Nx_file ), nFields
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Process data
    !=======================================================================

    ! Open Const file for output
    CALL OpenCnBinary( nAllFields )

    ! Regrid fields from the various raw data files
    CALL ProcessCn2dAsmNx( nFields, fields )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing CN output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    IF ( do2x25  ) CLOSE( CN_2x25  )
    IF ( do4x5   ) CLOSE( CN_4x5   )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE MerraMakeConst %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE MerraMakeCn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetNFields
!
! !DESCRIPTION: Returns the list of fields and number of fields to regrid
!  for each MERRA raw data file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetNFields( dataList, nFields, fields )
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN)  :: dataList   ! Comma-sep'd field name list
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: nFields    ! Number of fields
    CHARACTER(LEN=*), INTENT(OUT) :: fields(:)  ! Array of field names
! 
! !REVISION HISTORY: 
!  26 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: F

    ! Split the data field list by commas into an array
    CALL makeCharArrayFromCharList( dataList, ',', fields )
    
    ! Compute the number of data fields we will process
    nFields = 0

    DO F = 1, SIZE( fields )
       IF ( TRIM( fields(F) ) /= ''      .and. &
            TRIM( fields(F) ) /= 'none' ) THEN
          nFields = nFields + 1
       ENDIF
    ENDDO

  END SUBROUTINE GetNFields
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessCn2dChmNx
!
! !DESCRIPTION: Subroutine ProcessCn2dChmNx regrids the MERRA met fields from 
!  the "const\_2d\_chm\_Nx" file and saves to the GEOS-Chem file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ProcessCn2dAsmNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  26 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  03 Aug 2010 - R. Yantosca - Now use name8 to save name to binary file
!  11 Aug 2010 - R. Yantosca - Renamed to ProcessCn2dChmNx
!  11 Aug 2010 - R. Yantosca - Now get constant data from const_2d_asm_Nx file
!  12 Aug 2010 - R. Yantosca - Remove pointers
!  12 Aug 2010 - R. Yantosca - Renamed to ProcessCn2dAsmNx
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: F, DD, HH

    ! Arrays
    REAL*4                  :: Q    ( I05x0666, J05x0666, 1 )
    REAL*4                  :: Q2x25( I2x25,    J2x25       )
    REAL*4                  :: Q4x5 ( I4x5,     J4x5        )

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fileHDF
    CHARACTER(LEN=MAX_CHAR) :: fileName
    CHARACTER(LEN=MAX_CHAR) :: file2x25
    CHARACTER(LEN=MAX_CHAR) :: file4x5
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE ProcessCn2dAsmNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Create filename from the template
    fileHDF = TRIM( dataDirHDF ) // TRIM( const_2d_asm_Nx_file )
    CALL expandDate( fileHDF, yyyymmdd, 000000 )

    ! Echo info
    msg = '%%% Reading ' // TRIM( fileHDF )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Open the HDF4-EOS file for input
    CALL He4SetVerbose( VERBOSE )
    CALL He4GridOpen( fileHDF )
    CALL He4GridGetDimInfo
    CALL He4GridReadX
    CALL He4GridReadY
    CALL He4GridReadZ
    CALL He4GetNymdNhms

    !=======================================================================
    ! Process data
    !=======================================================================

    ! Loop over data fields
    DO F = 1, nFields

       ! Save field name into an 9-char variable. 
       ! This will truncate field names longer than 8 chars.
       name  = TRIM( fields(F) )

       ! Skip if fieldname is empty
       IF ( name == '' ) CYCLE

       ! Save to 8-character variable for output
       name8 = name

       !-----------------------------
       ! Read data 
       !-----------------------------

       ! Zero data arrays
       Q     = 0e0
       Q2x25 = 0e0
       Q4x5  = 0e0

       ! Read data
       msg = '%%% Reading     ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL He4GridReadData( name, Q )

       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0

       !-----------------------------
       ! Regrid data to 2x25, 4x5
       !-----------------------------
       msg = '%%% Regridding  ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Regrid to 2 x 2.5
       IF ( do2x25 ) THEN
          CALL RegridMerraNTo2x25( 0, Q(:,:,1), Q2x25 )
       ENDIF

       ! Regrid to 4x5 
       IF ( do4x5 ) THEN
          CALL RegridMerraNTo4x5( 0, Q(:,:,1), Q4x5 )
       ENDIF

       !----------------------------
       ! Write binary output
       !----------------------------
       msg = '%%% Writing     ' // name // ' to disk'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Save date in local variables (pick 2000/01/01)
       DD  = 20000101
       HH  = 000000
                 
       ! 2 x 2.5 output
       IF ( do2x25 ) THEN
          CALL WriteBinary( CN_2x25, name8, DD, HH, Q2x25, I2x25, J2x25 )
       ENDIF
          
       ! 4 x 5 output
       IF ( do4x5 ) THEN
          CALL WriteBinary( CN_4x5,  name8, DD, HH, Q4x5,  I4x5,  J4x5  )
       ENDIF
    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Detach from grid and close HDF file
    msg = '%%% Closing ' // TRIM( fileHDF )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL He4GridClose( fileHDF )
    CALL He4CleanUpIndexFields

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE ProcessCn2dAsmFx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE ProcessCn2dAsmNx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: OpenCnBinary
!
! !DESCRIPTION: Subroutine OpenCmBinary opens the constant ("CN") files for 
!  output in GEOS-Chem binary format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OpenCnBinary( nFields )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: nFields   ! # of Const fields in the file
!
! !REVISION HISTORY: 
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  06 Aug 2010 - R. Yantosca - Change MERRA identifier from "G5" to "ME"
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=2)        :: numFields
    CHARACTER(LEN=8)        :: ident
    CHARACTER(LEN=MAX_CHAR) :: fileName

    ! Save to character
    WRITE( numFields, '(i2.2)' ) nFields

    !-------------------------------
    ! Open 2 x 2.5 file for output
    !-------------------------------
    IF ( do2x25 ) THEN

       ! Replace date in the output filename
       fileName = dataTmpl2x25
       CALL expandDate( fileName, 20000101, 000000 )
       CALL StrRepl( fileName, '%%', 'cn' )

       ! binary
       ident = 'ME 22 ' // numFields
       CALL OpenBinary( CN_2x25, ident, filename )

    ENDIF

    !-------------------------------
    ! Open 4 x 5 file for output
    !-------------------------------
    IF ( do4x5 ) THEN

       ! Replace tokens in the output filename
       fileName = dataTmpl4x5
       CALL expandDate( fileName, 20000101, 000000 )
       CALL StrRepl( fileName, '%%', 'cn' )

       ! Binary
       ident = 'ME 45 ' // numFields
       CALL OpenBinary( CN_4x5, ident, filename )

    ENDIF

  END SUBROUTINE OpenCnBinary
!EOC
END MODULE MerraCnModule

