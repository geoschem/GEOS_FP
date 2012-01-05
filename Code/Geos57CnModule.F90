!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57CnModule
!
! !DESCRIPTION: Module Geos57CnModule contains routines to create the 
!  GEOS-Chem constant data files from the GEOS-5.7.x raw data.
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
  USE Geos57UtilityModule

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
! !REMARKS:
!  NOTE: Hardwire the constant data file to 00:00 GMT on 2011/01/01.
!
! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version, based on MERRA
!  20 Dec 2011 - R. Yantosca - Updates to achieve COARDS netCDF compliance
!  04 Jan 2012 - R. Yantosca - Updated comments, cosmetic changes
!  04 Jan 2012 - R. Yantosca - Add extra global attributes
!  04 Jan 2012 - R. Yantosca - Now reference Geos57UtilityModule
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
  SUBROUTINE NcOutFileDef( X,        Y,           T,    &
                           xMid,     YMid,        time, &
                           gridName, outFileName, fOut )
!
! !INPUT PARAMETERS:
! 
    INTEGER,          INTENT(IN)    :: X             ! Longitude dimension
    INTEGER,          INTENT(IN)    :: Y             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: T             ! Time dimension
    REAL*4,           INTENT(IN)    :: xMid(X)       ! Array of lon centers
    REAL*4,           INTENT(IN)    :: yMid(Y)       ! Array of lat centers
    INTEGER,          INTENT(IN)    :: time(T)       ! Array of times
    CHARACTER(LEN=*), INTENT(IN)    :: gridName      ! Name of the grid
    CHARACTER(LEN=*), INTENT(IN)    :: outFileName   ! Output file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fOut          ! Output netCDF file ID
!
! !REMARKS:
!  NOTE: Hardwire the constant data file to date 2011/01/01; 00:00 GMT.

! !REVISION HISTORY: 
!  25 Oct 2011 - R. Yantosca - Initial version
!  21 Dec 2011 - R. Yantosca - Modified for COARDS compliance
!  21 Dec 2011 - R. Yantosca - Also now write index arrays
!  22 Dec 2011 - R. Yantosca - Added gridName argument
!  04 Jan 2012 - R. Yantosca - Now use separate attributes "begin_date" and
!                              "begin_time" for the "time" index array
!  04 Jan 2012 - R. Yantosca - Now use all lowercase for index array names
!  04 Jan 2012 - R. Yantosca - Add extra global attributes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,   DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr, msg
    INTEGER            :: idLon,   idLat,   idTime,  vId,  oMode

    ! Arrays
    INTEGER            :: var1(1), var3(3)

    !=========================================================================
    ! %%% BEGINNING OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Echo info
    WRITE( 6, 100 ) TRIM( gridName )
100 FORMAT ( '%%% Defining netCDF file vars & attrs for ', a' grid' )

    ! Open netCDF file for writing
    CALL NcCr_Wr( fOut, TRIM( outFileName ) )

    ! Turn filling off
    CALL NcSetFill( fOut, NF_NOFILL, oMode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------
  
    ! Title string
    lName = 'GEOS-5.7.2 constant (CN) fields for GEOS-Chem'
    CALL NcDef_Glob_Attributes( fOut, 'Title',       TRIM( lName )   )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',     TRIM( lName )   )

    ! References
    lName = "www.geos-chem.org; wiki.geos-chem.org"
    CALL NcDef_Glob_Attributes( fOut, 'References',  TRIM( lName )   )

    ! Filename
    lName = outFileName
    CALL NcDef_Glob_Attributes( fOut, 'Filename',    TRIM( lName )   )
    
    ! History
    sysTime = SystemTimeStamp()
    lName = 'File generated on: ' // TRIM( sysTime )
    CALL NcDef_Glob_Attributes( fOut, 'History' ,    TRIM( lName )   )

    ! Format
    lName = "NetCDF-3" ;
    CALL NcDef_Glob_Attributes( fOut, 'Format' ,     TRIM( lName )   )

    ! Conventions
    lName = 'COARDS'
    CALL NcDef_Glob_Attributes( fOut, 'Conventions', TRIM( lName )   )

    ! Version
    lName = 'GEOS-5,7.2'
    CALL NcDef_Glob_Attributes( fOut, 'Version',     TRIM( lName )   )

    ! Model
    lName = 'GEOS5'
    CALL NcDef_Glob_Attributes( fOut, 'Model',       TRIM( lName )   )

    ! NLayers
    lName = '72'
    CALL NcDef_Glob_Attributes( fOut, 'Nlayers',     TRIM( lName )   )

    ! Start Date (hardwire to 2011/01/01)
    lName = '20110101' 
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',  TRIM( lName )   )

    ! Start Time
    lName = '00:00:00.0'
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',  TRIM( lName )   )

    ! End Date (hardwire to 2011/01/01)
    lName = '20110101' 
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',    TRIM( lName )   )

    ! End Time
    lName = '00:00:00.0'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',    TRIM( lName )   )

    ! Delta-time
    lName = '000000'
    CALL NcDef_Glob_Attributes( fOut, 'Delta_time',  TRIM( lName )   )

    ! Pick DI and DJ attributes based on the grid
    SELECT CASE ( TRIM( gridName ) )
       CASE( 'native', 'SEA4CRS', 'nested China' )
          DI = '0.3125'
          DJ = '0.25'
       CASE ( 'nested 0.5 x 0.625' )
          DI = '0.625'
          DJ = '0.5'
       CASE( '2 x 2.5 global' )
          DI = '2.5'
          DJ = '2'
       CASE( '4 x 5 global' )
          DI = '5'
          DJ = '4'
    END SELECT

    ! Delta-lon
    CALL NcDef_Glob_Attributes( fOut, 'Delta_lon',   TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_lat',   TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'lon',  X, idLon  )
    CALL NcDef_Dimension( fOut, 'lat',  Y, idLat  )
    CALL NcDef_Dimension( fOut, 'time', 1, idTime )

    ! Longitude index array
    vId     = 0
    var1    = (/ idLon /)
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )
  
    ! Latitude index array
    var1    = (/ idLat /)
    vId     = vId + 1
    lName   = 'latitude'
    units   = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    ! Time index array (hardwire date to 2011/01/01)
    var1    = (/ idTime /)
    vId     = vId + 1
    lName   = 'time'
    units   = UnitsForTime( 20110101 )
    delta_t = '0000-00-00 00:00:00'
    begin_d = '20110101'
    begin_t = '000000'
    incr    = '000000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! FRLAKE
    IF ( StrPos( 'FRLAKE', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of lake type in grid box' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAKE', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FRLAND
    IF ( StrPos( 'FRLAND', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land in grid box' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAND', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FRLANDICE
    IF ( StrPos( 'FRLANDIC', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land ice in grid box' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLANDIC', NF_FLOAT, 3, var3, vId  )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FROCEAN
    IF ( StrPos( 'FROCEAN', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of ocean in grid box' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FROCEAN', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! PHIS
    IF ( StrPos( 'PHIS', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface geopotential' 
       units = 'm2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PHIS', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    !=========================================================================
    ! %%% END OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! End the definition section
    CALL NcEnd_def( fOut )

    ! Write index arrays
    CALL NcWr( xMid, fOut, 'lon',  (/ 1 /), (/ X /) )
    CALL NcWr( yMid, fOut, 'lat',  (/ 1 /), (/ Y /) )
    CALL NcWr( time, fOut, 'time', (/ 1 /), (/ T /) )

    ! Echo info    
    msg = '%%%%%% LEAVING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE NcOutFileDef
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Geos57MakeCn 
!
! !DESCRIPTION: Routine Geos57MakeCn is the the driver routine for 
! \begin{enumerate}
! \item Extracting constant data fields (surface values) from 
!       the MERRA raw data files (HDF4-EOS format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
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
!  04 Jan 2012 - R. Yantosca - Updated comments
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
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    INTEGER                 :: time(1)
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57MakeCn %%%%%%%%%%'
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
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Hours in the day
    time = (/ 0 /)
    
    ! Open nested China output file
    IF ( doNestCh ) THEN
       fName = dataTmplNestCh
       gName = 'SEA4CRS'
       CALL ExpandDate  ( fName,     20110101,  000000      )      
       CALL StrRepl     ( fName,     '%%',     'cn'         )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  1,           &
                          xMid_025x03125(I0_ch:I1_ch),       &
                          yMid_025x03125(J0_ch:J1_ch),       &
                          time,      gName,    fName,        &
                          fOutNestCh                        )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = dataTmpl2x25
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     20110101,  000000      )      
       CALL StrRepl     ( fName,     '%%',      'cn'        )
       CALL NcOutFileDef( I2x25,     J2x25,     1,           &
                          xMid_2x25, yMid_2x25, time,        &
                          gName,     fName,     fOut2x25    )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = dataTmpl4x5
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     20110101,  000000      )      
       CALL StrRepl     ( fName,     '%%',     'cn'         )
       CALL NcOutFileDef( I4x5,      J4x5,      1,           &
                          xMid_4x5,  yMid_4x5,  time,        &
                          gName,     fName,     fOut4x5     )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================

    ! Regrid fields from the various raw data files
    CALL ProcessCn2dAsmNx( nFields,    fields,               &
                           fOutNestCh, fOut2x25, fOut4x5    )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing CN output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57MakeCn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Geos57MakeCn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessCn2dAsmNx
!
! !DESCRIPTION: Subroutine ProcessCn2dAsmNx regrids the GEOS-5.7.x met fields 
!  from the "const\_2d\_asm\_Nx" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ProcessCn2dAsmNx( nFields,    fields,  &
                               fOutNestCh, fOut2x25, fOut4x5 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
    INTEGER,          INTENT(IN) :: fOutNestCh  ! NestCh  netCDF file ID
    INTEGER,          INTENT(IN) :: fOut2x25    ! 2 x 2.5 netCDF file ID
    INTEGER,          INTENT(IN) :: fOut4x5     ! 4 x 5   netCDF file ID
!
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop variables
    INTEGER                 :: F

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: st2d(2),  st3d(3)
    INTEGER                 :: ct2d(2),  ct3d(3)

    ! Data arrays
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, 1 )
    REAL*4                  :: Q2x25( I2x25,      J2x25         )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5          )

    ! Pointers
    REAL*4,  POINTER        :: QNest(:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! 2 x 2.5 global grid       
    IF ( do2x25 ) THEN
       CALL NcGet_DimLen( fOut2x25,   'lon',  X2x25   )
       CALL NcGet_DimLen( fOut2x25,   'lat',  Y2x25   ) 
       CALL NcGet_DimLen( fOut2x25,   'time', T2x25   )
    ENDIF

    ! 4x5 global grid
    IF ( do4x5 ) THEN
       CALL NcGet_DimLen( fOut4x5,    'lon',  X4x5    )
       CALL NcGet_DimLen( fOut4x5,    'lat',  Y4x5    )   
       CALL NcGet_DimLen( fOut4x5,    'time', T4x5    )
    ENDIF
    
    !=======================================================================
    ! Open input file
    ! NOTE: For constant file, hardwire date to 2011/01/01
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
    CALL NcGet_DimLen( fIn, 'lon',  X )
    CALL NcGet_DimLen( fIn, 'lat',  Y ) 
    CALL NcGet_DimLen( fIn, 'time', T )

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

       ! Index arrays for netCDF
       st3d = (/ 1, 1, 1 /)
       ct3d = (/ X, Y, T /)

       ! Read data
       msg = '%%% Reading     ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )

       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0
       
       !-----------------------------
       ! Regrid data
       !-----------------------------
       msg = '%%% Regridding  ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Regrid to 2 x 2.5
       IF ( do2x25 ) THEN
          CALL RegridGeos57to2x25( 0, Q(:,:,1), Q2x25 )
       ENDIF

       ! Regrid to 4x5 
       IF ( do4x5 ) THEN
          CALL RegridGeos57To4x5( 0, Q(:,:,1), Q4x5 )
       ENDIF

       !-----------------------------
       ! Write netCDF output
       !-----------------------------

       msg = '%%% Archiving   ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Special handing
       IF ( TRIM( name ) == 'FRLANDICE' ) name ='FRLANDIC'

       ! Nested China
       IF ( doNestCh ) THEN
          QNest => Q( I0_ch:I1_ch, J0_ch:J1_ch, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestCh, YNestCh, TNestCh /)
          CALL NcWr( QNest, fOutNestCh, TRIM( name ), st3d, ct3d )          
       ENDIF

       ! Write 2 x 2.5 data
       IF ( do2x25 ) THEN
          st3d = (/ 1,     1,     1     /)
          ct3d = (/ X2x25, Y2x25, T2x25 /)
          CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )          
       ENDIF
       
       ! Write 4x5 data
       IF ( do4x5 ) THEN
          st3d = (/ 1,    1,    1    /)
          ct3d = (/ X4x5, Y4x5, T4x5 /)
          CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )          
       ENDIF

    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

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

