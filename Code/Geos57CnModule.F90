!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57CnModule
!
! !DESCRIPTION: Module Geos57CnModule contains routines to create the 
!  GEOS-Chem average constant data files from the GEOS-5.7.x raw data.
!\\
!\\
! !INTERFACE: 

MODULE Geos57CnModule
! 
! !USES:
!
  ! GEOS-5.7.x data modules
  USE CharpakModule
  USE Geos57InputsModule
  USE Geos57RegridModule

  ! Modules for writing netCDF
  USE m_netcdf_io_create
  USE m_netcdf_io_define
  USE m_netcdf_io_write
  USE m_netcdf_io_close
  
  ! Modules for reading netCDF
  USE m_netcdf_io_open
  USE m_netcdf_io_close     
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read

  IMPLICIT NONE
  PRIVATE
  
  ! Include files
  INCLUDE "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Geos57MakeCn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: GetNFields
  PRIVATE :: ProcessCn2dAsmNx
!
! !REVISION HISTORY:
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
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
! !IROUTINE: NcOutFileDef
!
! !DESCRIPTION: Subroutine NcOutFileDef pre-defines variable names and 
!  attributes that will be added to the netCDF output files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcOutFileDef( XDim, YDim, outFileName, fOut )
!
! !INPUT PARAMETERS:
! 
    INTEGER,          INTENT(IN)    :: XDim          ! Longitude dimension
    INTEGER,          INTENT(IN)    :: YDim          ! Latitude dimension
    CHARACTER(LEN=*), INTENT(IN)    :: outFileName   ! Output file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fOut          ! Output netCDF file ID
!
! !REVISION HISTORY: 
!  25 Oct 2011 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: lName,   units
    INTEGER            :: idLon,   idLat,   idTime
    INTEGER            :: vId,     omode

    ! Arrays
    INTEGER            :: var1(1), var2(2), var3(3)

    !=========================================================================
    ! %%% BEGINNING OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Defining netCDF file variables & attributes ==='

    ! Open netCDF file for writing
    CALL NcCr_Wr( fOut, TRIM( outFileName ) )

    ! Turn filling off
    CALL NcSetFill( fOut, NF_NOFILL, omode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------
 
    ! Title string
    lName = 'GEOS-5.7.2 Constant Fields for GEOS-Chem'
    CALL NcDef_Glob_Attributes( fOut, 'title',       TRIM( lName ) )

    ! Version history
    lName = 'Version: 25 Oct 2011'
    CALL NcDef_Glob_Attributes( fOut, 'history',     TRIM( lName ) )

    ! Conventions
    lName = 'COARDS'
    CALL NcDef_Glob_Attributes( fOut, 'Conventions', TRIM( lName ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and variables
    !-------------------------------------------------------------------------

    ! Define geospatial dimensions for netCDF file
    CALL NcDef_Dimension( fOut, 'XDim', XDim, idLon  )
    CALL NcDef_Dimension( fOut, 'YDim', YDim, idLat  )
    CALL NcDef_Dimension( fOut, 'TDim', 1,    idTime )


    ! Longitude
    vId   = 0
    var2  = (/ idLon, idTime /)
    lName = 'Longitude'
    units = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 2, var2, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )         )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )         )
    CALL NcDef_Var_Attributes( fOut, vId, '_FillValue', '1.0e+15f'           )
  
    ! Latitude
    var2  = (/ idLat, idTime /)
    vId   = vId + 1
    lName = 'Latitude'
    units = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 2, var2, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name', TRIM( lName )         )
    CALL NcDef_Var_attributes( fOut, vId, 'units',     TRIM( units )         ) 

    ! FRLAKE
    IF ( StrPos( 'FRLAKE', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of lake type in grid box' 
       units = 'fraction'
       CALL NcDef_Variable      ( fOut, 'FRLAKE', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )      )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )      )
    ENDIF

    ! FRLAND
    IF ( StrPos( 'FRLAND', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land in grid box' 
       units = 'fraction'
       CALL NcDef_Variable      ( fOut, 'FRLAND', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )      )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )      )
    ENDIF

    ! FRLANDICE
    IF ( StrPos( 'FRLAND', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land ice in grid box' 
       units = 'fraction'
       CALL NcDef_Variable      ( fOut, 'FRLANDICE', NF_FLOAT, 3, var3, vId  )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )      )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )      )
    ENDIF

    ! FROCEAN
    IF ( StrPos( 'FROCEAN', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of ocean in grid box' 
       units = 'fraction'
       CALL NcDef_Variable      ( fOut, 'FROCEAN', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )      )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )      )
    ENDIF

    ! PHIS
    IF ( StrPos( 'PHIS', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface geopotential' 
       units = 'm2 s-2'
       CALL NcDef_Variable      ( fOut, 'PHIS', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name', TRIM( lName )      )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',     TRIM( units )      )
    ENDIF

    !=========================================================================
    ! %%% END OF NETCDF DEFINITION SECTION %%%
    !=========================================================================
    CALL NcEnd_def( fOut )

  END SUBROUTINE NcOutFileDef
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Geos57MakeCnFiles 
!
! !DESCRIPTION: Routine Geos57MakeConstFiles is the the driver routine for 
! \begin{enumerate}
! \item Extracting constant data fields (surface values) from 
!       the MERRA raw data files (HDF4-EOS format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data in a format that GEOS-Chem can read.
! \end{enumerate}
! This routine is called directly from the main program Geos57Driver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57MakeCn
!
! !REVISION HISTORY: 
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
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
    CHARACTER(LEN=MAX_CHAR) :: fName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57MakeConst %%%%%%%%%%'
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
    ! Open files for  output; define variables and attributes
    !=======================================================================

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = dataTmpl2x25
       CALL ExpandDate  ( fName, yyyymmdd, 000000          )
       CALL StrRepl     ( fName, '%%',     'cn'            )
       CALL NcOutFileDef( I2x25, J2x25,    fName, fOut2x25 )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = dataTmpl4x5
       CALL ExpandDate  ( fName, yyyymmdd, 000000          )
       CALL StrRepl     ( fName, '%%',     'cn'            )
       CALL NcOutFileDef( I4x5,  J4x5,     fName, fOut4x5  )
    ENDIF

    ! Regrid fields from the various raw data files
    CALL ProcessCn2dAsmNx( nFields, fields )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing CN output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close 2 x 2.5 output file
    IF ( do2x25 ) THEN
       CALL NcCl( fOut2x25 )
    ENDIF

    ! Close 4 x 5 output file
    IF ( do2x25 ) THEN
       CALL NcCl( fOut4x5  )
    ENDIF

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57MakeConst %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Geos57MakeCn
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
    ! Loop variables
    INTEGER                 :: F,       DD,     HH

    ! Variables for netCDF I/O
    INTEGER                 :: xDim,    YDim,   TDim
    INTEGER                 :: st2d(2), st3d(3)
    INTEGER                 :: ct2d(2), ct3d(3)

    ! Index arrays
    REAL*4                  :: lon  ( I025x03125                )
    REAL*4                  :: lat  ( J025x03125                )

    ! Data arrays
    REAL*4                  :: Q    ( I025x03125, J025x03125, 1 )
    REAL*4                  :: Q2x25( I2x25,      J2x25         )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5          )

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: fNameNested
    CHARACTER(LEN=MAX_CHAR) :: fName2x25
    CHARACTER(LEN=MAX_CHAR) :: fName4x5
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Initialization & open input file
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE ProcessCn2dAsmNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Create input filename from the template
    fNameInput = TRIM( inputDataDir ) // TRIM( const_2d_asm_Nx_file )
    CALL expandDate( fNameInput, yyyymmdd, 000000 )

    ! Echo info
    msg = '%%% Reading ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Open the netCDF4 file for input
    CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
    
    ! Get the dimensions from the netCDF file
    CALL NcGet_DimLen( fIn, 'lon',  XDim )
    CALL NcGet_DimLen( fIn, 'lat',  YDim )
    CALL NcGet_DimLen( fIn, 'time', TDim )

    ! Create index arrays for netCDF
    st3d = (/ 1,    1,    1    /)
    ct3d = (/ XDim, YDim, TDim /)

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
       CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )

       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0
       
       print*, '### min, max: ', minval( Q ), maxval( Q )

!       !-----------------------------
!       ! Regrid data to 2x25, 4x5
!       !-----------------------------
!       msg = '%%% Regridding  ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Regrid to 2 x 2.5
!       IF ( do2x25 ) THEN
!          CALL RegridGeos57to2x25( 0, Q(:,:,1), Q2x25 )
!       ENDIF
!
!       ! Regrid to 4x5 
!       IF ( do4x5 ) THEN
!          CALL RegridGeos57To4x5( 0, Q(:,:,1), Q4x5 )
!       ENDIF
!
!       !----------------------------
!       ! Write binary output
!       !----------------------------
!       msg = '%%% Writing     ' // name // ' to disk'
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Save date in local variables (pick 2000/01/01)
!       DD  = 20000101
!       HH  = 000000
!              
!       ! 2 x 2.5 output
!       IF ( do_nested ) THEN
!          CALL WriteBinary(  name8, DD, HH, Q2x25, I2x25, J2x25 )
!       ENDIF
!          
!       ! 2 x 2.5 output
!       IF ( do2x25 ) THEN
!          CALL WriteBinary( CN_2x25, name8, DD, HH, Q2x25, I2x25, J2x25 )
!       ENDIF
!          
!       ! 4 x 5 output
!       IF ( do4x5 ) THEN
!          CALL WriteBinary( CN_4x5,  name8, DD, HH, Q4x5,  I4x5,  J4x5  )
!       ENDIF
    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================


!    IF ( doNested ) THEN
!       CALL NcOutFileDef()
!    ENDIF

    ! Close 2 x 2.5 output file (if necessary)
    IF ( do2x25 ) THEN
       msg = '%%% Closing ' // TRIM( fName2x25 )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fOut2x25 )
    ENDIF

    ! Close 4 x 5 output file (if necessary)
    IF ( do4x5 ) THEN
       msg = '%%% Closing ' // TRIM( fName4x5 )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fOut4x5 )
    ENDIF

    ! Close input file
    msg = '%%% Closing ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL NcCl( fIn )

    ! Echo info    
    msg = '%%%%%% EXITING ROUTINE ProcessCn2dAsmNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE ProcessCn2dAsmNx
!EOC
END MODULE Geos57CnModule

