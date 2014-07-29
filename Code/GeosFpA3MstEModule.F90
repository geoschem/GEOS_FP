!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpA3MstEModule
!
! !DESCRIPTION: Module GeosFpA3MstEModule contains routines to create the 
!  GEOS-Chem average 3-hr data files from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE: 

MODULE GeosFpA3MstEModule
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

# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GeosFpMakeA3MstE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dMstNe
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  12 Jan 2012 - R. Yantosca - Now just save out the fields on level edges
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Renamed to GeosFpA3MstEModule
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
  SUBROUTINE NcOutFileDef( X,        Y,           Z,     T,    &
                           xMid,     yMid,        zEdge, time, &
                           gridName, outFileName, fOut )
!
! !INPUT PARAMETERS:
! 
    INTEGER,          INTENT(IN)    :: X             ! Longitude dimension
    INTEGER,          INTENT(IN)    :: Y             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: Z             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: T             ! Time dimension
    REAL*4,           INTENT(IN)    :: xMid(X)       ! Array of lon centers
    REAL*4,           INTENT(IN)    :: yMid(Y)       ! Array of lat centers
    REAL*4,           INTENT(IN)    :: zEdge(Z)      ! Array of alt edges
    INTEGER,          INTENT(IN)    :: time(T)       ! Array of times
    CHARACTER(LEN=*), INTENT(IN)    :: gridName      ! Name of the grid
    CHARACTER(LEN=*), INTENT(IN)    :: outFileName   ! Output file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fOut          ! Output netCDF file ID
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version
!  01 Feb 2012 - R. Yantosca - Make all global attribute names lowercase
!  01 Feb 2012 - R. Yantosca - Bug fix: do not use Z+1 for vertical dimension,
!                              as the proper value is passed in the arg list
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  23 Sep 2013 - R. Yantosca - Add calendar attribute to time
!  24 Sep 2013 - R. Yantosca - Bug fix: now use correct start & end dates
!  24 Sep 2013 - R. Yantosca - Now save dims in order: time, lev, lat, lon
!  08 Oct 2013 - R. Yantosca - Update CASE statement for gridName
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,    DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr,  msg, cal
    INTEGER            :: idLon,   idLat,   idLev,   idAp
    INTEGER            :: idBp,    idTime,  vId,     omode

    ! Arrays
    INTEGER            :: var1(1), var4(4)

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
    CALL NcSetFill( fOut, NF_NOFILL, omode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------
  
    ! Title string
    lName = 'GEOS-FP time-averaged 3-hour moist parameters on level edges (A3mstE), processed for GEOS-Chem input'
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
                                                              
    ! Start Date
    lName = yyyymmdd_string
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',           TRIM( lName ) )
                                                              
    ! Start Time                                              
    lName = '00:00:00.0'                                      
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',           TRIM( lName ) )
                                                              
    ! End Date
    lName = yyyymmdd_string
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',             TRIM( lName ) )
                                                              
    ! End Time                                                
    lName = '23:59:59.99999'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',             TRIM( lName ) )
                                                              
    ! Delta-time                                              
    lName = '030000'                                          
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
    CALL NcDef_Dimension( fOut, 'time', T,   idTime )
    CALL NcDef_Dimension( fOut, 'lev',  Z,   idLev  )
    CALL NcDef_Dimension( fOut, 'lat',  Y,   idLat  )
    CALL NcDef_Dimension( fOut, 'lon',  X,   idLon  )

    ! Time index array
    var1    = (/ idTime /)
    vId     = 0
    cal     = 'gregorian'
    lName   = 'time'
    units   = UnitsForTime( yyyymmdd )
    delta_t = '0000-00-00 03:00:00'
    begin_d = yyyymmdd_string
    begin_t = '000000'
    incr    = '030000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',      TRIM( cal     )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Level index array
    var1    = (/ idLev /)
    vId     = vId + 1
    lName   = 'levels'
    units   = '1'
    CALL NcDef_Variable      ( fOut, 'lev', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

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

    IF ( StrPos( 'CMFMC', tavg3_3d_mst_Ne_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Upward moist convective mass flux'
       units = 'kg m-2 s-2'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'CMFMC', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PFICU
    IF ( StrPos( 'PFICU', tavg3_3d_mst_Ne_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Downward flux of ice precipitation (convective)'
       units = 'kg m-2 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PFICU', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PFILSAN
    IF ( StrPos( 'PFILSAN', tavg3_3d_mst_Ne_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Downward flux of ice precipitation (large scale + anvil)'
       units = 'kg m-2 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PFILSAN', NF_FLOAT, 4, var4, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PFLCU
    IF ( StrPos( 'PFLCU', tavg3_3d_mst_Ne_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Downward flux of liquid precipitation (convective)'
       units = 'kg m-2 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PFLCU', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PFLLSAN
    IF ( StrPos( 'PFLLSAN', tavg3_3d_mst_Ne_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Downward flux of liquid precipitation (large scale + anvil)'
       units = 'kg m-2 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PFLLSAN', NF_FLOAT, 4, var4, vId    )
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
    CALL NcWr( xMid,  fOut, 'lon',  (/ 1 /), (/ X /) )
    CALL NcWr( yMid,  fOut, 'lat',  (/ 1 /), (/ Y /) )
    CALL NcWr( zEdge, fOut, 'lev',  (/ 1 /), (/ Z /) )
    CALL NcWr( time,  fOut, 'time', (/ 1 /), (/ T /) )

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
! !IROUTINE: GeosFpMakeA3MstE
!
! !DESCRIPTION: Routine GeosFpMakeA3MstE is the the driver routine for 
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (moist parameters) from 
!       the GEOS-FP raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeA3MstE()
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version, based on MERRA
!  12 Jan 2012 - R. Yantosca - Now process only fields on level edges
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Change and/or add attributes for COARDS
!  20 Sep 2013 - R. Yantosca - Now save CMFMC to the A3mstE file
!  20 Sep 2013 - R. Yantosca - Now save out nested Europe grid
!  23 Sep 2013 - R. Yantosca - Now define netCDF latitude such that the poles
!                              are at -90/+90.  This facilitates the GIGC
!                              using ESMF/MAPL.
!  08 Oct 2013 - R. Yantosca - Now save output to nested SE Asia grid (SE)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_3dMstNe
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dMstNe(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeA3MstE %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg3_3d_mst_Ne_data )

    ! Return the list of fields and number of fields to process
    ! from each of the GeosFp raw met data files
    CALL GetNFields( tavg3_3d_mst_Ne_data, nFields_3dMstNe, fields_3dMstNe )
    CALL GetNFields( allFieldsList,        nAllFields,      allFields      )
    
    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_mst_Ne_file ), nFields_3dMstNe
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested CH output file
    IF ( doNestCh ) THEN
       fName = TRIM( tempDirTmplNestCh ) // TRIM( dataTmplNestCh )
       gName = 'nested CH'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                     )    
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstE'                   )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  L025x03125+1, TIMES_A3,     &
                          xMid_025x03125(I0_ch:I1_ch),                      &
                          yMid_025x03125(J0_ch:J1_ch),                      &
                          zEdge_025x03125,                       a3Mins,    &
                          gName,     fName,        fOutNestCh              )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu ) THEN
       fName = TRIM( tempDirTmplNestEu ) // TRIM( dataTmplNestEu )
       gName = 'nested EU'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                     )   
       CALL StrRepl     ( fname,     '%%%%%%',  'A3mstE'                   )
       CALL NcOutFileDef( I_NestEu,  J_NestEu,  L025x03125+1,    TIMES_A3,  &
                          xMid_025x03125(I0_eu:I1_eu),                      &
                          yMid_025x03125(J0_eu:J1_eu),                      &
                          zEdge_025x03125,                       a3Mins,    &
                          gName,     fName,        fOutNestEu              )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa ) THEN
       fName = TRIM( tempDirTmplNestNa ) // TRIM( dataTmplNestNa )
       gName = 'nested NA'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                     )   
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstE'                   )
       CALL NcOutFileDef( I_NestNa,  J_NestNa,  L025x03125+1,    TIMES_A3,  &
                          xMid_025x03125(I0_na:I1_na),                      &
                          yMid_025x03125(J0_na:J1_na),                      &
                          zEdge_025x03125,                       a3Mins,    &
                          gName,     fName,        fOutNestNa              )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe ) THEN
       fName = TRIM( tempDirTmplNestSe ) // TRIM( dataTmplNestSe )
       gName = 'nested SE'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                     )   
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstE'                   )
       CALL NcOutFileDef( I_NestSe,  J_NestSe,  L025x03125+1,    TIMES_A3,  &
                          xMid_025x03125(I0_se:I1_se),                      &
                          yMid_025x03125(J0_se:J1_se),                      &
                          zEdge_025x03125,                       a3Mins,    &
                          gName,     fName,        fOutNestSe              )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                  )   
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'                )
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25+1,      TIMES_A3,  &
                          xMid_2x25, nc_yMid_2x25, zEdge_2x25,   a3Mins,    &
                          gName,     fName,        fOut2x25                )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                  )   
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'                )
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5+1,       TIMES_A3,  &
                          xMid_4x5,  nc_yMid_4x5,  zEdge_4x5,    a3Mins,    &
                          gName,     fName,        fOut4x5                 )
    ENDIF
    
    !----- (lzh,06/20/2014)------------
    ! Open nested 0625 CH output file
    IF ( doNestCh05 ) THEN
       fName = TRIM( tempDirTmplNestCh05 ) // TRIM( dataTmplNestCh05 )
       gName = 'nested CH 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'    )
       CALL NcOutFileDef( I_NestCh05,  J_NestCh05, L05x0625,  TIMES_A3,  &
                          xMid_05x0625(I0_ch05:I1_ch05),          &
                          yMid_05x0625(J0_ch05:J1_ch05),          &
                          zEdge_05x0625,                a3Mins,    &
                          gName,    fName,        fOut05NestCh          )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'    )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          zEdge_05x0625,                a3Mins,    &
                          gName,    fName,       fOut05NestEu          )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'    )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          zEdge_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestNa          )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe05 ) THEN
       fName = TRIM( tempDirTmplNestSe05 ) // TRIM( dataTmplNestSe05 )
       gName = 'nested SE 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstE'    )
       CALL NcOutFileDef( I_NestSe05,  J_NestSe05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_se05:I1_se05),          &
                          yMid_05x0625(J0_se05:J1_se05),          &
                          zEdge_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestSe         )
    ENDIF
    !------(finish edit)---------------
    
    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process3dMstNe( nFields_3dMstNe, fields_3dMstNe ) ! tavg3_3d_mst_Ne
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A3mstE output files'
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
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeA3mstE %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'   
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeA3MstE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dMstNe
!
! !DESCRIPTION: Subroutine Process3dMstNe regrids the GeosFp met fields from 
!  the "tavg3\_3d\_mst\_Ne" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dMstNe( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  11 Jan 2012 - R. Yantosca - Initial version
!  17 Jan 2012 - R. Yantosca - Bug fix: flip data in vertical immediately
!                              after reading.  Use pointers for efficiency
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  08 Oct 2013 - R. Yantosca - Now save output to nested EU and SE grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: F,        H,        L,       LR
    INTEGER                 :: hhmmss 

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        Z,       T
    INTEGER                 :: XNestCh,  YNestCh,  ZNestCh, TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  ZNestEu, TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  ZNestNa, TNestNa
    INTEGER                 :: XNestSe,  YNestSe,  ZNestSe, TNestSe
    INTEGER                 :: X2x25,    Y2x25,    Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,     Z4x5,    T4x5
    INTEGER                 :: st4d(4),  ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, L025x03125+1 )
    REAL*4                  :: Q2x25( I2x25,      J2x25,      L2x25+1      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5,       L4x5+1       )

    ! (lzh,06/20/2014) 0.5x0.625
    INTEGER                 :: XNestCh05,  YNestCh05, ZNestCh05, TNestCh05
    INTEGER                 :: XNestEu05,  YNestEu05, ZNestEu05, TNestEu05
    INTEGER                 :: XNestNa05,  YNestNa05, ZNestNa05, TNestNa05
    INTEGER                 :: XNestSe05,  YNestSe05, ZNestSe05, TNestSe05
    REAL*4,  TARGET         :: Q05    ( I05x0625, J05x0625, L05x0625+1 )    

    ! Pointer arrays
    REAL*4,  POINTER        :: Qflip(:,:,:)
    REAL*4,  POINTER        :: Ptr  (:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

     ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dMstNe %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested CH grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'lev',  ZNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested EU grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'lev',  ZNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested NA grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'lev',  ZNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
    ENDIF

    ! Nested SE grid
    IF ( doNestSe ) THEN
       CALL NcGet_DimLen( fOutNestSe, 'lon',  XNestSe )
       CALL NcGet_DimLen( fOutNestSe, 'lat',  YNestSe ) 
       CALL NcGet_DimLen( fOutNestSe, 'lev',  ZNestSe ) 
       CALL NcGet_DimLen( fOutNestSe, 'time', TNestSe )
    ENDIF

    ! 2 x 2.5 global grid       
    IF ( do2x25 ) THEN
       CALL NcGet_DimLen( fOut2x25,   'lon',  X2x25   )
       CALL NcGet_DimLen( fOut2x25,   'lat',  Y2x25   ) 
       CALL NcGet_DimLen( fOut2x25,   'lev',  Z2x25   ) 
       CALL NcGet_DimLen( fOut2x25,   'time', T2x25   )
    ENDIF

    ! 4x5 global grid
    IF ( do4x5 ) THEN
       CALL NcGet_DimLen( fOut4x5,    'lon',  X4x5    )
       CALL NcGet_DimLen( fOut4x5,    'lat',  Y4x5    )   
       CALL NcGet_DimLen( fOut4x5,    'lev',  Z4x5    )   
       CALL NcGet_DimLen( fOut4x5,    'time', T4x5    )
    ENDIF

    ! (lzh, 06/21/2014) 0.5x0.625
    ! Nested CH grid
    IF ( doNestCh05 ) THEN
       CALL NcGet_DimLen( fOut05NestCh, 'lon',  XNestCh05 )
       CALL NcGet_DimLen( fOut05NestCh, 'lat',  YNestCh05 )
       CALL NcGet_DimLen( fOut05NestCh, 'lev',  ZNestCh05 )
       CALL NcGet_DimLen( fOut05NestCh, 'time', TNestCh05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestEu05 ) THEN
       CALL NcGet_DimLen( fOut05NestEu, 'lon',  XNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lat',  YNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lev',  ZNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'time', TNestEu05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestNa05 ) THEN
       CALL NcGet_DimLen( fOut05NestNa, 'lon',  XNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lat',  YNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lev',  ZNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'time', TNestNa05 )
    ENDIF

    ! Nested SE grid
    IF ( doNestSe05 ) THEN
       CALL NcGet_DimLen( fOut05NestSe, 'lon',  XNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lat',  YNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lev',  ZNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'time', TNestSe05 )
    ENDIF

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_mst_Ne_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )

       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'lev',  Z ) 
       CALL NcGet_DimLen( fIn, 'time', T )

       !====================================================================
       ! Process data
       !====================================================================
       
       ! Loop over data fields
       DO F = 1, nFields
          
          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0
          Q05   = 0e0        ! (lzh,06/21/2014)

          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading    ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Start and count index arrays for netCDF
          ! (There is only one data block per file)
          st4d = (/ 1, 1, 1, 1 /)
          ct4d = (/ X, Y, Z, 1 /)
          
          ! Read data
          CALL NcRd( Q, fIn, TRIM( name ), st4d, ct4d )

          ! Strip fill values
          WHERE( Q == FILL_VALUE ) Q = 0e0    

          ! Flip levels in vertical
          Qflip => Q( :, :, Z:1:-1 )

          !-----------------------------------------------------------------
          ! Regrid data fields
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Loop over the A-3 times and vertical levels
          DO L = 1, Z

             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridGeosFpTo2x25( 0, Qflip(:,:,L), Q2x25(:,:,L) )
             ENDIF
                   
             ! Regrid to 4x5 
             IF ( do4x5 ) THEN
                CALL RegridGeosFpTo4x5 ( 0, Qflip(:,:,L), Q4x5(:,:,L)  )
             ENDIF
             ! Regrid to 0.5 x 0.625 (lzh, 06/21/2014)
             IF ( do05x0625 ) THEN
                CALL RegridGeosFpTo05x0625( 0, Qflip(:,:,L), Q05(:,:,L) )
             ENDIF
             
          ENDDO
                
          !-----------------------------------------------------------------
          ! Write netCDF output (all except QI, QL)
          !-----------------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
          ! Nested CH (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Qflip( I0_ch:I1_ch, J0_ch:J1_ch, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestCh, YNestCh, ZNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Qflip( I0_eu:I1_eu, J0_eu:J1_eu, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestEu, YNestEu, ZNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Qflip( I0_na:I1_na, J0_na:J1_na, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestNa, YNestNa, ZNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested SE (point to proper slice of global data)
          IF ( doNestSe ) THEN
             Ptr  => Qflip( I0_se:I1_se, J0_se:J1_se, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestSe, YNestSe, ZNestSe, 1 /)
             CALL NcWr( Ptr, fOutNestSe, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st4d = (/ 1,     1,     1,     H  /)
             ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st4d, ct4d )
          ENDIF
          
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st4d  = (/ 1,    1,    1,    H /)
             ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st4d, ct4d )
          ENDIF
          
          !-------(lzh, 06/21/2014)----------
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh05 ) THEN
             Ptr  => Q05( I0_ch05:I1_ch05, J0_ch05:J1_ch05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestCh05, YNestCh05, ZNestCh05, 1 /)
             CALL NcWr( Ptr, fOut05NestCh, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu05 ) THEN
             Ptr  => Q05( I0_eu05:I1_eu05, J0_eu05:J1_eu05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestEu05, YNestEu05, ZNestEu05, 1 /)
             CALL NcWr( Ptr, fOut05NestEu, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa05 ) THEN
             Ptr  => Q05( I0_na05:I1_na05, J0_na05:J1_na05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestNa05, YNestNa05, ZNestNa05, 1 /)
             CALL NcWr( Ptr, fOut05NestNa, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested SE (point to proper slice of global data)
          IF ( doNestSe05 ) THEN
             Ptr  => Q05( I0_se05:I1_se05, J0_se05:J1_se05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestSe05, YNestSe05, ZNestSe05, 1 /)
             CALL NcWr( Ptr, fOut05NestSe, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF
          !-------(finish edit)--------------

          ! Free pointer memory
          NULLIFY( Qflip )
       ENDDO

       !--------------------------------------------------------------------
       ! Close input file
       !--------------------------------------------------------------------
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )       
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% LEAVING ROUTINE Process3dMstNe %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dMstNe
!EOC
END MODULE GeosFpA3MstEModule

