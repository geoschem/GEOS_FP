!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpCnModule
!
! !DESCRIPTION: Module GeosFpCnModule contains routines to create the 
!  GEOS-Chem constant data files from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE: 

MODULE GeosFpCnModule
! 
! !USES:
!
  ! GEOS-FP data modules
  USE CharpakModule
  USE GeosFpInputsModule
  USE GeosFpRegridModule
  USE GeosFpUtilityModule

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
# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GeosFpMakeCn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: GetNFields
  PRIVATE :: ProcessCn2dAsmNx
!
! !REMARKS:
!  NOTE: Hardwire the constant data file to 00:00 GMT on 2011/01/01.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version, based on MERRA
!  20 Dec 2011 - R. Yantosca - Updates to achieve COARDS netCDF compliance
!  04 Jan 2012 - R. Yantosca - Updated comments, cosmetic changes
!  04 Jan 2012 - R. Yantosca - Add extra global attributes
!  04 Jan 2012 - R. Yantosca - Now reference GeosFpUtilityModule
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Renamed to GeosFpCnModule; adjusted for COARDS
!  08 Oct 2013 - R. Yantosca - Now save CH, EU, NA, SE nested grids in one pass
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
!  01 Feb 2012 - R. Yantosca - Make all global attribute names lowercase
!  19 Sep 2013 - R. Yantosca - Change and/or add attributes for COARDS standard
!  23 Sep 2013 - R. Yantosca - Add calendar attribute to time
!  24 Sep 2013 - R. Yantosca - Now write dims in order: time, lat, lon
!  08 Oct 2013 - R. Yantosca - Updated CASE statement for gridName
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,   DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr, msg,  cal
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
    lName = 'GEOS-FP constant parameters (CN), processed for GEOS-Chem input'
    CALL NcDef_Glob_Attributes( fOut, 'Title',                TRIM( lName ) )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',              TRIM( lName ) )

    ! References
    lName = "www.geos-chem.org; wiki.geos-chem.org"
    CALL NcDef_Glob_Attributes( fOut, 'References',           TRIM( lName ) )

    ! Filename
    lName = NotDir( outFileName )
    CALL NcDef_Glob_Attributes( fOut, 'Filename',             TRIM( lName ) )
    
    ! History
    sysTime = SystemTimeStamp()
    lName = 'File generated on: ' // TRIM( sysTime )
    CALL NcDef_Glob_Attributes( fOut, 'History' ,             TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ProductionDateTime',   TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ModificationDateTime', TRIM( lName ) )

    ! Format
    lName = "NetCDF-3" ;
    CALL NcDef_Glob_Attributes( fOut, 'Format' ,              TRIM( lName ) )
                                                              
    ! Format                                                  
    lName = "global" ;                                        
    CALL NcDef_Glob_Attributes( fOut, 'SpatialCoverage',      TRIM( lName ) )
                                                              
    ! Conventions                                             
    lName = 'COARDS'                                          
    CALL NcDef_Glob_Attributes( fOut, 'Conventions',          TRIM( lName ) )
                                                              
    ! Version                                                 
    lName = 'GEOS-FP'                                        
    CALL NcDef_Glob_Attributes( fOut, 'Version',              TRIM( lName ) )
                                                              
    ! Model                                                   
    lName = 'GEOS-5'                                         
    CALL NcDef_Glob_Attributes( fOut, 'Model',                TRIM( lName ) )
                                                              
    ! NLayers                                                 
    lName = '72'                                              
    CALL NcDef_Glob_Attributes( fOut, 'Nlayers',              TRIM( lName ) )
                                                              
    ! Start Date (hardwire to 2011/01/01)                     
    lName = '20110101'                                        
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',           TRIM( lName ) )
                                                              
    ! Start Time                                              
    lName = '00:00:00.0'                                      
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',           TRIM( lName ) )
                                                              
    ! End Date (hardwire to 2011/01/01)                       
    lName = '20110101'                                        
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',             TRIM( lName ) )
                                                              
    ! End Time                                                
    lName = '00:00:00.0'                                      
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',             TRIM( lName ) )
                                                              
    ! Delta-time                                              
    lName = '000000'                                          
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Time',           TRIM( lName ) )

    ! Pick DI and DJ attributes based on the grid
    SELECT CASE ( TRIM( gridName ) )
       CASE( 'native', 'nested CH', 'nested EU', 'nested NA', 'nested SE' )
          DI = '0.3125'
          DJ = '0.25'
!       CASE ( 'nested 0.5 x 0.625' )
! (lzh,06/21/2014)
       CASE( 'nested CH 05', 'nested EU 05', 'nested NA 05', 'nested SE 05' )
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
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lon',            TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lat',            TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'time', 1, idTime )
    CALL NcDef_Dimension( fOut, 'lat',  Y, idLat  )
    CALL NcDef_Dimension( fOut, 'lon',  X, idLon  )

    ! Time index array (hardwire date to 2011/01/01)
    var1    = (/ idTime /)
    vId     = 0
    cal     = 'gregorian'
    lName   = 'time'
    units   = UnitsForTime( 20110101 )
    delta_t = '0000-00-00 00:00:00'
    begin_d = '20110101'
    begin_t = '000000'
    incr    = '000000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',       TRIM( cal     )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Latitude index array
    var1    = (/ idLat /)
    vId     = vId + 1
    lName   = 'latitude'
    units   = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    ! Longitude index array
    var1    = (/ idLon /)
    vId     = vId + 1
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )
  
    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! FRLAKE
    IF ( StrPos( 'FRLAKE', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of lake type in grid box' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAKE', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRLAND
    IF ( StrPos( 'FRLAND', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land in grid box' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLAND', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRLANDICE
    IF ( StrPos( 'FRLANDIC', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of land ice in grid box' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRLANDIC', NF_FLOAT, 3, var3, vId  )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FROCEAN
    IF ( StrPos( 'FROCEAN', const_2d_asm_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of ocean in grid box' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FROCEAN', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
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
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
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
! !IROUTINE: GeosFpMakeCn 
!
! !DESCRIPTION: Routine GeosFpMakeCn is the the driver routine for 
! \begin{enumerate}
! \item Extracting constant data fields (surface values) from 
!       the GEOS-FP raw data files (netCDF-4 format)
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeCn
!
! !REVISION HISTORY: 
!  27 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  04 Jan 2012 - R. Yantosca - Updated comments
!  11 Jan 2012 - R. Yantosca - Now call StrCompress to remove white space
!                              in the input file name.
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Now save output to nested Europe grid
!  23 Sep 2013 - R. Yantosca - Now define netCDF latitude such that the poles
!                              are at -90/+90.  This facilitates the GIGC
!                              using ESMF/MAPL.
!  08 Oct 2013 - R. Yantosca - Now save output to the nested SE Asia grid (SE)
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
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeCn %%%%%%%%%%'
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
    
    ! Open nested CH output file
    IF ( doNestCh ) THEN
       fName = TRIM( tempDirTmplNestCh ) // TRIM( dataTmplNestCh )
       gName = 'nested CH'
       CALL ExpandDate  ( fName,     20110101,     000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,     1,           &
                          xMid_025x03125(I0_ch:I1_ch),          &
                          yMid_025x03125(J0_ch:J1_ch),          &
                          time,      gName,        fName,       &
                          fOutNestCh                           )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu ) THEN
       fName = TRIM( tempDirTmplNestEu ) // TRIM( dataTmplNestEu )
       gName = 'nested EU'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',     'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestEu,  J_NestEu,  TIMES_A1,       &
                          xMid_025x03125(I0_eu:I1_eu),          &
                          yMid_025x03125(J0_eu:J1_eu),          &
                          a1Mins,    gName,        fName,       &
                          fOutNestEu                           )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa ) THEN
       fName = TRIM( tempDirTmplNestNa ) // TRIM( dataTmplNestNa )
       gName = 'nested NA'
       CALL ExpandDate  ( fName,     20110101,     000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestNa,  J_NestNa,  1,              &
                          xMid_025x03125(I0_na:I1_na),          &
                          yMid_025x03125(J0_na:J1_na),          &
                          time,      gName,        fName,       &
                          fOutNestNa                            )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe ) THEN
       fName = TRIM( tempDirTmplNestSe ) // TRIM( dataTmplNestSe )
       gName = 'nested SE'
       CALL ExpandDate  ( fName,     20110101,     000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestSe,  J_NestSe,  1,              &
                          xMid_025x03125(I0_se:I1_se),          &
                          yMid_025x03125(J0_se:J1_se),          &
                          time,      gName,        fName,       &
                          fOutNestSe                            )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     20110101,     000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.           )
       CALL NcOutFileDef( I2x25,     J2x25,        1,           &
                          xMid_2x25, nc_yMid_2x25, time,        &
                          gName,     fName,        fOut2x25    )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     20110101,  000000         )      
       CALL StrRepl     ( fName,     '%%%%%%',  'CN    '       )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I4x5,      J4x5,      1,              &
                          xMid_4x5,  nc_yMid_4x5,  time,        &
                          gName,     fName,     fOut4x5        )
    ENDIF
    !----- (lzh,06/21/2014)------------
    ! Open nested CH output file
    IF ( doNestCh05 ) THEN
       fName = TRIM( tempDirTmplNestCh05 ) // TRIM( dataTmplNestCh05 )
       gName = 'nested CH 05'
       CALL ExpandDate  ( fName,     20110101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestCh05,  J_NestCh05,     1,           &
                          xMid_05x0625(I0_ch05:I1_ch05),          &
                          yMid_05x0625(J0_ch05:J1_ch05),          &
                          time,      gName,        fName,       &
                          fOut05NestCh                           )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     20110101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05,     1,       &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          time,    gName,        fName,       &
                          fOut05NestEu                           )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     20110101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05,  1,              &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          time,      gName,        fName,       &
                          fOut05NestNa                            )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe05 ) THEN
       fName = TRIM( tempDirTmplNestSe05 ) // TRIM( dataTmplNestSe05 )
       gName = 'nested SE 05'
       CALL ExpandDate  ( fName,     20110101,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'CN    '    )
       CALL StrCompress ( fName, RemoveAll=.TRUE.              )
       CALL NcOutFileDef( I_NestSe05,  J_NestSe05,  1,              &
                          xMid_05x0625(I0_se05:I1_se05),          &
                          yMid_05x0625(J0_se05:J1_se05),          &
                          time,      gName,        fName,       &
                          fOut05NestSe                            )
    ENDIF
    !------(finish edit)---------------
    
    !=======================================================================
    ! Process data
    !=======================================================================

    ! Regrid fields from the various raw data files
    CALL ProcessCn2dAsmNx( nFields, fields )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing CN output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( doNestEu ) CALL NcCl( fOutNestEu )
    IF ( doNestNa ) CALL NcCl( fOutNestNa )
    IF ( doNestSe ) CALL NcCl( fOutNestSe )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )
    ! (lzh, 06/20/2014) add nested 0.5x0.625
    IF ( doNestCh05 ) CALL NcCl( fOut05NestCh )
    IF ( doNestEu05 ) CALL NcCl( fOut05NestEu )
    IF ( doNestNa05 ) CALL NcCl( fOut05NestNa )
    IF ( doNestSe05 ) CALL NcCl( fOut05NestSe )
   
    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeCn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeCn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessCn2dAsmNx
!
! !DESCRIPTION: Subroutine ProcessCn2dAsmNx regrids the GEOS-FP met fields 
!  from the "const\_2d\_asm\_Nx" file and saves output to netCDF format.
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
!  04 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  21 Jun 2012 - R. Yantosca - Bug fix: remove 2nd instance of doNestCh
!  08 Oct 2013 - R. Yantosca - Now save out to SE Asia nested grid
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
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: XNestSe,  YNestSe,  TNestSe
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: st2d(2),  st3d(3)
    INTEGER                 :: ct2d(2),  ct3d(3)

    ! Data arrays
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, 1 )
    REAL*4                  :: Q2x25( I2x25,      J2x25         )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5          )
    ! (lzh,06/20/2014) 0.5x0.625
    INTEGER                 :: XNestCh05,  YNestCh05,  TNestCh05
    INTEGER                 :: XNestEu05,  YNestEu05,  TNestEu05
    INTEGER                 :: XNestNa05,  YNestNa05,  TNestNa05
    INTEGER                 :: XNestSe05,  YNestSe05,  TNestSe05
    REAL*4, TARGET          :: Q05  ( I05x0625, J05x0625, 1     )

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

    ! Nested CH grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested EU grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested NA grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
    ENDIF

    ! Nested SE grid
    IF ( doNestSe ) THEN
       CALL NcGet_DimLen( fOutNestSe, 'lon',  XNestSe )
       CALL NcGet_DimLen( fOutNestSe, 'lat',  YNestSe ) 
       CALL NcGet_DimLen( fOutNestSe, 'time', TNestSe )
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
    
    ! (lzh, 06/20/2014) 0.5x0.625
    ! Nested CH grid 0625
    IF ( doNestCh05 ) THEN
       CALL NcGet_DimLen( fOut05NestCh, 'lon',  XNestCh05 )
       CALL NcGet_DimLen( fOut05NestCh, 'lat',  YNestCh05 )
       CALL NcGet_DimLen( fOut05NestCh, 'time', TNestCh05 )
    ENDIF

    ! Nested EU grid 0625
    IF ( doNestEu05 ) THEN
       CALL NcGet_DimLen( fOut05NestEu, 'lon',  XNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lat',  YNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'time', TNestEu05 )
    ENDIF

    ! Nested NA grid 0625
    IF ( doNestNa05 ) THEN
       CALL NcGet_DimLen( fOut05NestNa, 'lon',  XNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lat',  YNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'time', TNestNa05 )
    ENDIF

    ! Nested SE grid 0625
    IF ( doNestSe05 ) THEN
       CALL NcGet_DimLen( fOut05NestSe, 'lon',  XNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lat',  YNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'time', TNestSe05 )
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
    msg = '%%% Opening ' // TRIM( fNameInput )
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
       Q05   = 0e0        ! (lzh,06/21/2014)

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
          CALL RegridGeosFpto2x25( 0, Q(:,:,1), Q2x25 )
       ENDIF

       ! Regrid to 4x5 
       IF ( do4x5 ) THEN
          CALL RegridGeosFpTo4x5( 0, Q(:,:,1), Q4x5 )
       ENDIF
       ! Regrid to 0.5 x 0.625   ! (lzh,06/21/2014)
       IF ( do05x0625 ) THEN
          CALL RegridGeosFpto05x0625( 0, Q(:,:,1), Q05 )
       ENDIF
       
       !-----------------------------
       ! Write netCDF output
       !-----------------------------

       msg = '%%% Archiving   ' // name
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Special handing
       IF ( TRIM( name ) == 'FRLANDICE' ) name ='FRLANDIC'

       ! Nested CH
       IF ( doNestCh ) THEN
          QNest => Q( I0_ch:I1_ch, J0_ch:J1_ch, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestCh, YNestCh, TNestCh /)
          CALL NcWr( QNest, fOutNestCh, TRIM( name ), st3d, ct3d )   
          NULLIFY( QNest )
       ENDIF

       ! Nested EU
       IF ( doNestEu ) THEN
          QNest => Q( I0_eu:I1_eu, J0_eu:J1_eu, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestEu, YNestEu, TNestEu /)
          CALL NcWr( QNest, fOutNestEu, TRIM( name ), st3d, ct3d )   
          NULLIFY( QNest )
       ENDIF

       ! Nested NA
       IF ( doNestNa ) THEN
          QNest => Q( I0_na:I1_na, J0_na:J1_na, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestNa, YNestNa, TNestNa /)
          CALL NcWr( QNest, fOutNestNa, TRIM( name ), st3d, ct3d )   
          NULLIFY( QNest )
       ENDIF

       ! Nested SE
       IF ( doNestSe ) THEN
          QNest => Q( I0_se:I1_se, J0_se:J1_se, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestSe, YNestSe, TNestSe /)
          CALL NcWr( QNest, fOutNestSe, TRIM( name ), st3d, ct3d )   
          NULLIFY( QNest )
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
    !-------(lzh, 06/21/2014) add 0.5x0.625------------
       ! Nested CH
       IF ( doNestCh05 ) THEN
          QNest => Q05( I0_ch05:I1_ch05, J0_ch05:J1_ch05, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestCh05, YNestCh05, TNestCh05 /)
          CALL NcWr( QNest, fOut05NestCh, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

       ! Nested EU
       IF ( doNestEu05 ) THEN
          QNest => Q05( I0_eu05:I1_eu05, J0_eu05:J1_eu05, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestEu05, YNestEu05, TNestEu05 /)
          CALL NcWr( QNest, fOut05NestEu, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

       ! Nested NA
       IF ( doNestNa05 ) THEN
          QNest => Q05( I0_na05:I1_na05, J0_na05:J1_na05, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestNa05, YNestNa05, TNestNa05 /)
          CALL NcWr( QNest, fOut05NestNa, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF

       ! Nested SE
       IF ( doNestSe05 ) THEN
          QNest => Q05( I0_se05:I1_se05, J0_se05:J1_se05, 1 )  ! Point to proper slice

          st3d  = (/ 1,       1,       1       /)
          ct3d  = (/ XNestSe05, YNestSe05, TNestSe05 /)
          CALL NcWr( QNest, fOut05NestSe, TRIM( name ), st3d, ct3d )
          NULLIFY( QNest )
       ENDIF
    !-------(finish edit)------------------------------
       
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
END MODULE GeosFpCnModule

