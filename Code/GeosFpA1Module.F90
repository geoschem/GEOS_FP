!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpA1Module
!
! !DESCRIPTION: Module GeosFpA1Module contains routines to create the 
!  GEOS-Chem average 1-hr data files from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE: 

MODULE GeosFpA1Module
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
  PUBLIC  :: GeosFpMakeA1
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process2dFlxNx
  PRIVATE :: Process2dLndNx
  PRIVATE :: Process2dRadNx
  PRIVATE :: Process2dSlvNx
  PRIVATE :: Process2dAlbedo
  PRIVATE :: GeosFpSeaIceBins
  PRIVATE :: GeosFpCreateLwi
  PRIVATE :: GeosFpRegridLwi
  PRIVATE :: GeosFpAdjustSnomas
  PRIVATE :: GeosFpProcessAlbedo
  PRIVATE :: GeosFpProcessTropp
!
! !DEFINED_PARAMETERS:
!
  INTEGER, PARAMETER :: N_ICE   = 10         ! # of sea ice bins
  REAL*4,  PARAMETER :: BINSIZE = 1e0/N_ICE  ! Sea ice bin size
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  05 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  09 Jan 2012 - R. Yantosca - Add driver routine Process2dAlbedo
!  09 Jan 2012 - R. Yantosca - Updated comments, cosmetic changes
!  11 Jan 2012 - R. Yantosca - Now put debugging kludges in #if blocks
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Renamed to GeosFpA1Module; adjusted for COARDS
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
! !REVISION HISTORY: 
!  05 Jan 2012 - R. Yantosca - Initial version, based on GeosFpCnModule
!  01 Feb 2012 - R. Yantosca - Make all global attribute names lowercase
!  20 Sep 2013 - R. Yantosca - Change and/or add attribute names for COARDS
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
    lName = 'GEOS-FP 1-hour time-averaged (A1) fields, processed for GEOS-Chem input'
    CALL NcDef_Glob_Attributes( fOut, 'Title',              TRIM( lName ) )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',            TRIM( lName ) )

    ! References
    lName = "www.geos-chem.org; wiki.geos-chem.org"
    CALL NcDef_Glob_Attributes( fOut, 'References',         TRIM( lName ) )

    ! Filename
    lName = NotDir( outFileName )
    CALL NcDef_Glob_Attributes( fOut, 'Filename',           TRIM( lName ) )
    
    ! History
    sysTime = SystemTimeStamp()
    lName = 'File generated on: ' // TRIM( sysTime )
    CALL NcDef_Glob_Attributes( fOut, 'History' ,           TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ProductionDateTime', TRIM( lName ) )

    ! Format
    lName = "NetCDF-3" ;
    CALL NcDef_Glob_Attributes( fOut, 'Format' ,            TRIM( lName ) )
                                                            
    ! Format                                                
    lName = "global" ;                                      
    CALL NcDef_Glob_Attributes( fOut, 'SpatialCoverage',    TRIM( lName ) )
                                                            
    ! Conventions                                           
    lName = 'COARDS'                                        
    CALL NcDef_Glob_Attributes( fOut, 'Conventions',        TRIM( lName ) )
                                                            
    ! Version                                               
    lName = 'GEOS5-FP'                                      
    CALL NcDef_Glob_Attributes( fOut, 'Version',            TRIM( lName ) )
                                                            
    ! Model                                                 
    lName = 'GEOS-FP'                                       
    CALL NcDef_Glob_Attributes( fOut, 'Model',              TRIM( lName ) )
                                                            
    ! NLayers                                               
    lName = '72'                                            
    CALL NcDef_Glob_Attributes( fOut, 'Nlayers',            TRIM( lName ) )
                                                            
    ! Start Date (hardwire to 2011/01/01)                   
    lName = '20110101'                                      
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',         TRIM( lName ) )
                                                            
    ! Start Time                                            
    lName = '00:00:00.0'                                    
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',         TRIM( lName ) )
                                                            
    ! End Date (hardwire to 2011/01/01)                     
    lName = '20110101'                                      
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',           TRIM( lName ) )
                                                            
    ! End Time                                              
    lName = '00:00:00.0'                                    
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',           TRIM( lName ) )
                                                            
    ! Delta-time                                            
    lName = '000000'                                        
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Time',         TRIM( lName ) )

    ! Pick DI and DJ attributes based on the grid
    SELECT CASE ( TRIM( gridName ) )
       CASE( 'native', 'nested CH', 'nested NA', 'nested EU' )
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
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lon',          TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lat',          TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'lon',  X, idLon  )
    CALL NcDef_Dimension( fOut, 'lat',  Y, idLat  )
    CALL NcDef_Dimension( fOut, 'time', T, idTime )

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
    units   = UnitsForTime( yyyymmdd )
    delta_t = '0000-00-00 01:00:00'
    begin_d = yyyymmdd_string
    begin_t = '000000'
    incr    = '010000'
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

    ! ALBEDO
    IF ( StrPos( 'ALBEDO', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface albedo' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'ALBEDO', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! CLDTOT
    IF ( StrPos( 'CLDTOT', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Total cloud fraction' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'CLDTOT', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! EFLUX
    IF ( StrPos( 'EFLUX', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Latent heat flux (positive upward)'
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'EFLUX', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! EVAP
    IF ( StrPos( 'EVAP', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface evaporation' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'EVAP', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRSEAICE
    IF ( StrPos( 'FRSEAICE', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of sea ice on surface' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRSEAICE', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FRSNO
    IF ( StrPos( 'FRSNO', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fractional snow-covered area' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRSNO', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! GRN
    IF ( StrPos( 'GRN', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Vegetation greenness fraction' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GRN', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! GWETROOT
    IF ( StrPos( 'GWETROOT', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Root zone soil wetness' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GWETROOT', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! GWETTOP
    IF ( StrPos( 'GWETTOP', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Top soil wetness' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GWETTOP', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! HFLUX
    IF ( StrPos( 'HFLUX', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Sensible heat flux (positive upward)' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'HFLUX', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! LAI
    IF ( StrPos( 'LAI', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Leaf area index' 
       units = 'm2 m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LAI', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! LWI (derived from FRLANDICE + other fields)
    IF ( StrPos( 'FRSEAICE', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Land/water/ice flags' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LWI', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! LWGNT
    IF ( StrPos( 'LWGNT', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Net longwave flux at the ground' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LWGNT', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF


    ! LWTUP
    IF ( StrPos( 'LWTUP', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Upward longwave flux at top of atmosphere (TOA)' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LWTUP', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PARDF
    IF ( StrPos( 'PARDF', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface downward PAR diffuse flux' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PARDF', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PARDR
    IF ( StrPos( 'PARDR', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface downward PAR beam flux' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PARDR', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PBLH
    IF ( StrPos( 'PBLH', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Planetary boundary layer height above surface' 
       units = 'm'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PBLH', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PRECANV
    IF ( StrPos( 'PRECANV', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface precipitation flux from anvils' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECANV', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PRECCON
    IF ( StrPos( 'PRECCON', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface precipitation flux from convection' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECCON', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PRECLSC
    IF ( StrPos( 'PRECLSC', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface precipitation flux from large-scale' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECLSC', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PRECSNO
    IF ( StrPos( 'PRECLSC', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface precipitation flux from snow' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECSNO', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PRECTOT
    IF ( StrPos( 'PRECTOT', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Total surface precipitation flux' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECTOT', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! QV2M
    IF ( StrPos( 'QV2M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Specific humidity at 2m above the displacement height' 
       units = 'kg kg-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'QV2M', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    IF ( StrPos( 'FRSEAICE', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN

       ! SEAICE00
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 0-10% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE00', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )

       ! SEAICE10
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 10-20% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE10', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
       
       ! SEAICE20
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 20-30% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE20', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )

       ! SEAICE30
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 30-40% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE30', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )

       ! SEAICE40
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 40-50% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE40', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
       
       ! SEAICE50
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 50-60% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE50', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
       
       ! SEAICE60
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 60-70% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE60', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
       
       ! SEAICE70
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 70-80% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE70', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
       
       ! SEAICE80
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 80-90% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE80', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )

       ! SEAICE90
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 90-100% sea ice coverage' 
       units = '1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE90', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! SLP
    IF ( StrPos( 'SLP', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Sea level pressure' 
       units = 'hPa'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SLP', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! SNODP
    IF ( StrPos( 'SNODP', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Snow depth' 
       units = 'm'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SNODP', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! SNOMAS
    IF ( StrPos( 'SNOMAS', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Snow mass' 
       units = 'kg m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SNOMAS', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! SWGDN
    IF ( StrPos( 'SWGDN', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface incident shortwave flux' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SWGDN', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! SWGNT
    IF ( StrPos( 'SWGNT', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Net surface downward shortwave flux' 
       units = 'W m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SWGNT', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! TROPPT
    IF ( StrPos( 'TROPPT', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Temperature-based tropopause pressure' 
       units = 'hPa'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'TROPPT', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! TS
    IF ( StrPos( 'TS', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface skin temperature' 
       units = 'K'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'TS', NF_FLOAT, 3, var3, vId         )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! T2M
    IF ( StrPos( 'T2M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Temperature 2m above displacement height' 
       units = 'K'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'T2M', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! U10M
    IF ( StrPos( 'U10M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Eastward wind 10m above displacement height' 
       units = 'm s-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'U10M', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! USTAR
    IF ( StrPos( 'USTAR', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Friction velocity' 
       units = 'm s-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'USTAR', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! V10M
    IF ( StrPos( 'V10M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Northward wind 10m above displacement height' 
       units = 'm s-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'V10M', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! Z0M
    IF ( StrPos( 'Z0M', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Roughness length, momentum' 
       units = 'm'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'Z0M', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! saj
    ! T10M 
    IF ( StrPos( 'T10M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       vId   = vId + 1
       lName = 'Temperature at 10 m above the displacement height'
       units = 'K'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'T10M', NF_FLOAT, 3, var3, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

   ! Q850
   IF ( StrPos( 'Q850', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)
       vId   = vId + 1
       lName = 'Specific humidity at 850 hPa'
       units = 'kg kg-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'Q850', NF_FLOAT, 3, var3, vId       )
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
! !IROUTINE: GeosFpMakeA1
!
! !DESCRIPTION: Routine GeosFpMakeA1
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (surface values) from 
!       the GEOS-FP raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to disk in netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeA1
!
! !REVISION HISTORY:
!  05 Jan 2012 - R. Yantosca - Initial version, based on GeosFpA3Module.F90
!  11 Jan 2012 - R. Yantosca - Now call StrCompress to remove white space
!                              in the input file name.
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_2dFlxNx
    INTEGER                 :: nFields_2dLndNx
    INTEGER                 :: nFields_2dRadNx
    INTEGER                 :: nFields_2dSlvNx
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_2dFlxNx(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_2dLndNx(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_2dRadNx(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_2dSlvNx(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeA1 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg1_2d_flx_Nx_data ) // ',' // &
                    TRIM( tavg1_2d_lnd_Nx_data ) // ',' // &
                    TRIM( tavg1_2d_rad_Nx_data ) // ',' // &
                    TRIM( tavg1_2d_slv_Nx_data )

    ! Return the list of fields and number of fields to process
    ! from each of the MERRA raw met data files
    CALL GetNFields( tavg1_2d_flx_Nx_data, nFields_2dFlxNx, fields_2dFlxNx )
    CALL GetNFields( tavg1_2d_lnd_Nx_data, nFields_2dLndNx, fields_2dLndNx )
    CALL GetNFields( tavg1_2d_rad_Nx_data, nFields_2dRadNx, fields_2dRadNx )
    CALL GetNFields( tavg1_2d_slv_Nx_data, nFields_2dSlvNx, fields_2dSlvNx )
    CALL GetNFields( allFieldsList,        nAllFields,      allFields      )
    
    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( tavg1_2d_flx_Nx_file ), nFields_2dFlxNx
    WRITE( IU_LOG, 100 ) TRIM( tavg1_2d_lnd_Nx_file ), nFields_2dLndNx
    WRITE( IU_LOG, 100 ) TRIM( tavg1_2d_rad_Nx_file ), nFields_2dRadNx
    WRITE( IU_LOG, 100 ) TRIM( tavg1_2d_slv_Nx_file ), nFields_2dSlvNx
    WRITE( IU_LOG, 110 ) N_ICE
    WRITE( IU_LOG, 120 ) 1
    WRITE( IU_LOG, 130 ) nAllFields + N_ICE + 1

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% # of fractional sea ice bins      : ', i5 )
120 FORMAT( '%%% # of land/water/ice flags fields  : ', i5 )
130 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested China output file
    IF ( doNestCh ) THEN
       fName = TRIM( tempDirTmplNestCh ) // TRIM( dataTmplNestCh )
       gName = 'nested CH'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.       )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  TIMES_A1,    &
                          xMid_025x03125(I0_ch:I1_ch),       &
                          yMid_025x03125(J0_ch:J1_ch),       &
                          a1Mins,    gName,     fName,       &
                          fOutNestCh                        )
    ENDIF

    ! Open nested Europe output file
    IF ( doNestEu ) THEN
       fName = TRIM( tempDirTmplNestEu ) // TRIM( dataTmplNestEu )
       gName = 'nested EU'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.       )
       CALL NcOutFileDef( I_NestEu,  J_NestEu,  TIMES_A1,    &
                          xMid_025x03125(I0_eu:I1_eu),       &
                          yMid_025x03125(J0_eu:J1_eu),       &
                          a1Mins,    gName,     fName,       &
                          fOutNestEu                        )
    ENDIF

    ! Open nested North America output file
    IF ( doNestNa ) THEN
       fName = TRIM( tempDirTmplNestNa ) // TRIM( dataTmplNestNa )
       gName = 'nested NA'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.       )
       CALL NcOutFileDef( I_NestNa,  J_NestNa,  TIMES_A1,    &
                          xMid_025x03125(I0_na:I1_na),       &
                          yMid_025x03125(J0_na:J1_na),       &
                          a1Mins,    gName,     fName,       &
                          fOutNestNa                        )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.       )
       CALL NcOutFileDef( I2x25,     J2x25,     TIMES_A1,    &
                          xMid_2x25, yMid_2x25, a1Mins,      &
                          gName,     fName,     fOut2x25    )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A1    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.       )
       CALL NcOutFileDef( I4x5,      J4x5,      TIMES_A1,    &
                          xMid_4x5,  yMid_4x5,  a1Mins,      &
                          gName,     fName,     fOut4x5     )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process2dFlxNx ( nFields_2dFlxNx, fields_2dFlxNx )  ! tavg1_2d_flx_Nx
    CALL Process2dLndNx ( nFields_2dLndNx, fields_2dLndNx )  ! tavg1_2d_lnd_Nx 
    CALL Process2dRadNx ( nFields_2dRadNx, fields_2dRadNx )  ! tavg1_2d_rad_Nx
    CALL Process2dSlvNx ( nFields_2dSlvNx, fields_2dSlvNx )  ! tavg1_2d_slv_Nx
    CALL Process2dAlbedo(                                 )  ! Sfc albedo
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A1 output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( doNestEu ) CALL NcCl( fOutNestEu )
    IF ( doNestNa ) CALL NcCl( fOutNestNa )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeA1 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeA1
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process2dFlxNx
!
! !DESCRIPTION: Subroutine Process2dFlxNx regrids the GEOS-FP met fields 
!  from the "tavg1\_2d\_flx\_Nx" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process2dFlxNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  05 Jan 2012 - R. Yantosca - Initial version, based on GeosFpCnModule.F90
!  06 Jan 2012 - R. Yantosca - Now call GeosFpCreateLwi to make the LWI field
!  06 Jan 2012 - R. Yantosca - Now call GeosFpSeaIceBins to compute the
!                              fractional sea ice bins SEAICE{00..90}
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  11 Jan 2012 - R. Yantosca - Now put debugging kludge in an #if block
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  15 Feb 2012 - R. Yantosca - Bug fix: save out proper nested grid window
!                              for sea ice data.  Now use pointer.
!  19 Sep 2013 - R. Yantosca - Also save out nested Europe grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: H,        F,        S,      hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: ct3d(3),  st3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q      ( I025x03125, J025x03125        )
    REAL*4, TARGET          :: lwi    ( I025x03125, J025x03125        )
    REAL*4, TARGET          :: ice    ( I025x03125, J025x03125, N_ICE )
    REAL*4                  :: Q2x25  ( I2x25,      J2x25             )
    REAL*4                  :: lwi2x25( I2x25,      J2x25             )
    REAL*4                  :: ice2x25( I2x25,      J2x25,      N_ICE )
    REAL*4                  :: Q4x5   ( I4x5,       J4x5              )
    REAL*4                  :: lwi4x5 ( I4x5,       J4x5              )
    REAL*4                  :: ice4x5 ( I4x5,       J4x5,       N_ICE )

    ! Pointers
    REAL*4, POINTER         :: ptr(:,:)

     ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=8       ) :: name2
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process2dFlxNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested Europe grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested North America grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
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
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A1

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a1Mins(H) / 60 ) * 10000 + 3000
       
       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg1_2d_flx_Nx_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )

       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'time', T )

       !====================================================================
       ! Process data
       !====================================================================

       ! Loop over data fields
       DO F = 1, nFields

          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Skip if the fieldname is empty
          IF ( name == '' .or. name == 'PS' ) CYCLE

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0
                    
          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading     ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Start and count index arrays for netCDF
          ! (There is only one data block per file)
          st3d  = (/ 1, 1, 1 /)
          ct3d  = (/ X, Y, 1 /)

          ! Read data from file
          CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )
          
          ! Replace missing values with zeroes
          WHERE( Q == FILL_VALUE ) Q = 0e0
          
          !-----------------------------------------------------------------
          ! Pre-regrid special handling: create derived fields
          !-----------------------------------------------------------------
          IF ( name == 'FRSEAICE' ) THEN

             ! Echo info
             msg = '%%% Computing fractional sea ice coverage and LWI flags'
             WRITE( IU_LOG, '(a)' ) TRIM( msg )

             IF ( doNative ) THEN

                ! Create the LWI fields on the global grid (if necessary)
                CALL GeosFpCreateLwi( Q,          mapNative,              &
                                      I025x03125, J025x03125, lwi        )

                ! Bin sea ice for native grid output 
                CALL GeosFpSeaIceBins( Q,         BINSIZE,    mapNative,  &
                                       ice,       I025x03125, J025x03125 )
             ENDIF

             IF ( doNestCh ) THEN

                !----------------------------------------------------------
                ! NESTED CHINA GRID: land/water/ice flags
                !----------------------------------------------------------
                Ptr  => lwi( I0_ch:I1_ch, J0_ch:J1_ch )
                st3d = (/ 1,       1,       H /)
                ct3d = (/ XNestCh, YNestCh, 1 /)
                CALL NcWr( Ptr, fOutNestCh, 'LWI', st3d, ct3d )
                NULLIFY( Ptr )

                !----------------------------------------------------------
                ! NESTED CHINA GRID: sea ice bins
                !----------------------------------------------------------
                DO S = 1, N_ICE
                   WRITE( name2, 200 ) S-1
 200               FORMAT( 'SEAICE', i1, '0' )
                   Ptr  => ice( I0_ch:I1_ch, J0_ch:J1_ch, S )
                   st3d = (/ 1,       1,       H  /)
                   ct3d = (/ XNestCh, YNestCh, 1  /)
                   CALL NcWr( Ptr, fOutNestCh, name2, st3d, ct3d )
                   NULLIFY( Ptr )
                ENDDO

             ENDIF

             IF ( doNestEu ) THEN

                !----------------------------------------------------------
                ! NESTED EU GRID: land/water/ice flags
                !----------------------------------------------------------
                Ptr  => lwi( I0_eu:I1_eu, J0_eu:J1_eu )
                st3d = (/ 1,       1,       H /)
                ct3d = (/ XNestEu, YNestEu, 1 /)
                CALL NcWr( Ptr, fOutNestEu, 'LWI', st3d, ct3d )
                NULLIFY( Ptr )

                !----------------------------------------------------------
                ! NESTED EU GRID: sea ice bins
                !----------------------------------------------------------
                DO S = 1, N_ICE
                   WRITE( name2, 200 ) S-1
                   Ptr  => ice( I0_eu:I1_eu, J0_eu:J1_eu, S )
                   st3d = (/ 1,       1,       H  /)
                   ct3d = (/ XNestEu, YNestEu, 1  /)
                   CALL NcWr( Ptr, fOutNestEu, name2, st3d, ct3d )
                   NULLIFY( Ptr )
                ENDDO

             ENDIF

             IF ( doNestNa ) THEN

                !----------------------------------------------------------
                ! NESTED NA GRID: land/water/ice flags
                !----------------------------------------------------------
                Ptr  => lwi( I0_na:I1_na, J0_na:J1_na )
                st3d = (/ 1,       1,       H /)
                ct3d = (/ XNestNa, YNestNa, 1 /)
                CALL NcWr( Ptr, fOutNestNa, 'LWI', st3d, ct3d )
                NULLIFY( Ptr )

                !----------------------------------------------------------
                ! NESTED NA GRID: sea ice bins
                !----------------------------------------------------------
                DO S = 1, N_ICE
                   WRITE( name2, 200 ) S-1
                   Ptr  => ice( I0_na:I1_na, J0_na:J1_na, S )
                   st3d = (/ 1,       1,       H  /)
                   ct3d = (/ XNestNa, YNestNa, 1  /)
                   CALL NcWr( Ptr, fOutNestNa, name2, st3d, ct3d )
                   NULLIFY( Ptr )
                ENDDO

             ENDIF

             IF ( do2x25 ) THEN 

                !-----------------------------------------------------------
                ! 2 x 2.5 GRID: Land/water/ice flags
                !-----------------------------------------------------------

                ! Create the 2 x 2.5 LWI field
                CALL GeosFpCreateLwi( Q, mapTo2x25, I2x25, J2x25, lwi2x25 )

                ! Write LWI to disk
                st3d = (/ 1,     1,     H  /)
                ct3d = (/ X2x25, Y2x25, 1  /)
                CALL NcWr( lwi2x25, fOut2x25, 'LWI', st3d, ct3d )

                !-----------------------------------------------------------
                ! 2 x 2.5 GRID: Sea ice bins
                !-----------------------------------------------------------
                
                ! Bin sea ice for 2 x 2.5 output 
                CALL GeosFpSeaIceBins( Q,       BINSIZE, mapTo2x25,  &
                                       ice2x25, I2x25,   J2x25        )

                ! Write each sea ice field to disk
                DO S = 1, N_ICE
                   WRITE( name2, 200 ) S-1
                   st3d = (/ 1,     1,     H  /)
                   ct3d = (/ X2x25, Y2x25, 1  /)
                   CALL NcWr( ice2x25(:,:,S), fOut2x25, name2, st3d, ct3d )
                ENDDO

             ENDIF

             IF ( do4x5 ) THEN

                !-----------------------------------------------------------
                ! 4 x 5 GRID: Land/water/ice flags
                !-----------------------------------------------------------

                ! Create the 4 x 5 LWI field
                CALL GeosFpCreateLwi( Q, mapTo4x5, I4x5, J4x5, lwi4x5 )

                ! Write LWI to disk
                st3d = (/ 1,    1,    H  /)
                ct3d = (/ X4x5, Y4x5, 1  /)
                CALL NcWr( lwi4x5, fOut4x5, 'LWI', st3d, ct3d )

                !-----------------------------------------------------------
                ! 4 x 5 GRID: Sea ice bins
                !-----------------------------------------------------------

                ! Bin sea ice for 4x5 output
                CALL GeosFpSeaIceBins( Q,       BINSIZE, mapTo4x5,   &
                                       ice4x5,  I4x5,    J4x5         )

                ! Write each sea ice bin to disk
                DO S = 1, N_ICE
                   WRITE( name2, 200 ) S-1
                   st3d = (/ 1,    1,    H  /)
                   ct3d = (/ X4x5, Y4x5, 1  /)
                   CALL NcWr( ice4x5(:,:,S), fOut4x5, name2, st3d, ct3d )
                ENDDO
             ENDIF
          ENDIF
             
          !-----------------------------------------------------------------
          ! Regrid to 2 x 2.5 &  4 x 5
          !-----------------------------------------------------------------
          msg = '%%% Regridding  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Do the regridding
          IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, Q, Q2x25 )
          IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, Q, Q4x5  )

          !-----------------------------------------------------------------
          ! Post-regrid special handling
          !-----------------------------------------------------------------
          SELECT CASE( name )
             CASE( 'PRECANV', 'PRECCON', 'PRECLSC', 'PRECTOT', 'USTAR' )
                ! These fields are always positive-definite
                IF ( do2x25 ) WHERE( Q2x25 < 0e0 ) Q2x25 = 0e0
                IF ( do4x5  ) WHERE( Q4x5  < 0e0 ) Q4x5  = 0e0
             CASE DEFAULT
                ! Do Nothing
          END SELECT

          !----------------------------------------------------------------
          ! Write netCDF output
          !----------------------------------------------------------------
          msg = '%%% Archiving   ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestCh, YNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Q( I0_eu:I1_eu, J0_eu:J1_eu )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestEu, YNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Q( I0_na:I1_na, J0_na:J1_na )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestNa, YNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st3d = (/ 1,     1,     H  /)
             ct3d = (/ X2x25, Y2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )
          ENDIF
          
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st3d = (/ 1,    1,    H /)
             ct3d = (/ X4x5, Y4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )
          ENDIF
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
    msg = '%%%%%% LEAVING ROUTINE Process2dFlxNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process2dFlxNx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process2dLndNx
!
! !DESCRIPTION:  Subroutine Process2dLndNx regrids the GEOS-FP met fields 
!  from the "tavg1\_2d\_lnd\_Nx" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process2dLndNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  05 Jan 2012 - R. Yantosca - Initial version, based on GeosFpCnModule.F90
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  11 Jan 2012 - R. Yantosca - Now put debugging kludge in an #if block
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Also save out nested EU grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: H,        F,        hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: ct3d(3),  st3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q    ( I025x03125, J025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5       )

    ! Pointers
    REAL*4, POINTER         :: ptr(:,:)

     ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process2dFlxNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested Europe grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested North American grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
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
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A1

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a1Mins(H) / 60 ) * 10000 + 3000
       
       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg1_2d_lnd_Nx_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )

       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'time', T )

       !====================================================================
       ! Process data
       !====================================================================

       ! Loop over data fields
       DO F = 1, nFields

          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Skip if the fieldname is empty
          IF ( name == '' ) CYCLE

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0

          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading     ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Start and count index arrays for netCDF
          ! (There is only one data block per file)
          st3d  = (/ 1, 1, 1 /)
          ct3d  = (/ X, Y, 1 /)

          ! Read data from file
          CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )
          
          ! Replace missing values with zeroes
          WHERE( Q == FILL_VALUE ) Q = 0e0

          !-----------------------------------------------------------------
          ! Pre-regrid handling
          !-----------------------------------------------------------------

          ! Adjust SNOMAS to be consistent w/ GEOS-Chem usage
          IF ( name == 'SNOMAS' ) CALL GeosFpAdjustSnomas( Q )  

          !-----------------------------------------------------------------
          ! Do the regridding!
          !-----------------------------------------------------------------
          msg = '%%% Regridding  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Regrid to global grids
          IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, Q, Q2x25 )
          IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, Q, Q4x5  )
 
          !-----------------------------------------------------------------
          ! Post-regrid handling
          !-----------------------------------------------------------------
          SELECT CASE( name )
             CASE( 'FRSNO', 'GRN',   'GWETROOT', 'GWETTOP', 'LAI',  &
                   'PARDF', 'PARDR', 'SNODP',    'SNOMAS'          )
                ! These fields are always positive-definite
                IF ( do2x25 ) WHERE( Q2x25 < 0e0 ) Q2x25 = 0e0
                IF ( do4x5  ) WHERE( Q4x5  < 0e0 ) Q4x5  = 0e0
             CASE DEFAULT
                ! Nothing
          END SELECT

          !-----------------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------------
          msg = '%%% Archiving   ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestCh, YNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Q( I0_eu:I1_eu, J0_eu:J1_eu )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestEu, YNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Q( I0_na:I1_na, J0_na:J1_na )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestNa, YNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF
          
          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st3d = (/ 1,     1,     H  /)
             ct3d = (/ X2x25, Y2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )
          ENDIF
          
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st3d = (/ 1,    1,    H /)
             ct3d = (/ X4x5, Y4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )
          ENDIF
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
    msg = '%%%%%% LEAVING ROUTINE Process2dLndNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process2dLndNx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process2dRadNx
!
! !DESCRIPTION:  Subroutine Process2dRadNx regrids the GEOS-FP met fields 
!  from the "tavg1\_2d\_rad\_Nx" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process2dRadNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  05 Jan 2012 - R. Yantosca - Initial version, based on GeosFpA3Module.F90
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  11 Jan 2012 - R. Yantosca - Now put debugging kludge in an #if block
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Now save nested Europe grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: H,        F,        hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: ct3d(3),  st3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q    ( I025x03125, J025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5       )

    ! Pointers
    REAL*4, POINTER         :: ptr(:,:)

     ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=8       ) :: name2
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process2dRadNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested Europe grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested North America grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
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
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A1

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a1Mins(H) / 60 ) * 10000 + 3000
       
       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg1_2d_rad_Nx_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )
       
       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'time', T )
       
       !====================================================================
       ! Process data
       !====================================================================

       ! Loop over data fields
       DO F = 1, nFields

          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Skip certain fields
          IF ( name == '' .or. name == 'PS' .or. name == 'ALBEDO' ) CYCLE
          
          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0
          
          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )
          
          ! Skip if the field is empty
          IF ( name == '' ) CYCLE
            
          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading    ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Start and count index arrays for netCDF
          ! (There is only one data block per file)
          st3d  = (/ 1, 1, 1 /)
          ct3d  = (/ X, Y, 1 /)

          ! Read data from file
          CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )
          
          ! Replace missing values with zeroes
          WHERE( Q == FILL_VALUE ) Q = 0e0

          !-----------------------------------------------------------------
          ! Regrid to coarse resolution
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Regrid
          IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, Q, Q2x25 )
          IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, Q, Q4x5  )

          !-----------------------------------------------------------------
          ! Post-regrid handling
          !-----------------------------------------------------------------
          SELECT CASE( name )

             ! These fields should be positive-definite
             CASE( 'CLDTOT', 'LWGNT', 'LWTUP', 'SWGDN', 'SWTUP' )
                IF ( do2x25 ) WHERE( Q2x25 < 0e0 ) Q2x25 = 0e0
                IF ( do4x5  ) WHERE( Q4x5  < 0e0 ) Q4x5  = 0e0

             CASE DEFAULT
                ! Do Nothing

          END SELECT

          !-----------------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestCh, YNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Q( I0_eu:I1_eu, J0_eu:J1_eu )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestEu, YNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF
          
          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Q( I0_na:I1_na, J0_na:J1_na )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestNa, YNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st3d = (/ 1,     1,     H  /)
             ct3d = (/ X2x25, Y2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )
          ENDIF
          
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st3d = (/ 1,    1,    H /)
             ct3d = (/ X4x5, Y4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )
          ENDIF
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
    msg = '%%%%%% LEAVING ROUTINE Process2dRadNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process2dRadNx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process2dSlvNx
!
! !DESCRIPTION:  Subroutine Process2dSlvNx regrids the GEOS-FP met fields 
!  from the "tavg1\_2d\_slv\_Nx" file and saves output to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process2dSlvNx( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version, based on GeosFpA3Module.F90
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  11 Jan 2012 - R. Yantosca - Now put debugging kludge in an #if block
!  17 Jan 2012 - R. Yantosca - Divide native-grid winds by pressures
!                              after regridding
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Now save out nested Europe grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: H,        F,        hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: ct3d(3),  st3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q    ( I025x03125, J025x03125 )
    REAL*4, TARGET          :: P    ( I025x03125, J025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25      )
    REAL*4                  :: P2x25( I2x25,      J2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5       )
    REAL*4                  :: P4x5 ( I4x5,       J4x5       )

    ! Pointers
    REAL*4, POINTER         :: ptr(:,:)

     ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process2dSlvNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested Europe grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested China grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
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
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A1

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a1Mins(H) / 60 ) * 10000 + 3000
       
       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg1_2d_slv_Nx_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )

       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'time', T )

       !====================================================================
       ! Process surface pressure (need for use below)
       !====================================================================

       msg = '%%% Reading     PS'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st3d  = (/ 1, 1, 1 /)
       ct3d  = (/ X, Y, 1 /)

       ! Read data from file
       CALL NcRd( P, fIn, 'PS', st3d, ct3d )
       
       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0

       ! Convert from [hPa] to [Pa]
       Q = Q / 100e0

       ! Regrid to 2 x 2.5
       IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, P, P2x25 )
       IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, P, P4x5  )

       !====================================================================
       ! Process all other data fields
       !====================================================================

       ! Loop over data fields
       DO F = 1, nFields

          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Skip if the fieldname is empty
          IF ( name == '' ) CYCLE

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0

          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading    ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )

          ! Start and count index arrays for netCDF
          ! (There is only one data block per file)
          st3d  = (/ 1, 1, 1 /)
          ct3d  = (/ X, Y, 1 /)

          ! Read data from file
          CALL NcRd( Q, fIn, TRIM( name ), st3d, ct3d )

          ! Fill missing values in TROPP and other fields
          SELECT CASE( name )
             CASE( 'TROPPT', 'TROPPV', 'TROPPB' )
                CALL GeosFpProcessTropp( Q )
             CASE DEFAULT
                WHERE( Q == FILL_VALUE ) Q = 0e0
          END SELECT
          
          !-----------------------------------------------------------------
          ! Pre-regrid handling
          !-----------------------------------------------------------------
          SELECT CASE( name )
             CASE(  'SLP', 'TROPPT', 'TROPPV', 'TROPPB' )
                Q = Q / 100e0                      ! Pa -> hPa 
             CASE( 'U10M', 'V10M' )               
                Q = Q * P                          ! Multiply winds by pressure
             CASE DEFAULT
                ! Nothing
          END SELECT

          !-----------------------------------------------------------------
          ! Regrid data 
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Regrid
          IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, Q, Q2x25 )
          IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, Q, Q4x5  )

          !-----------------------------------------------------------------
          ! Post-regrid handling
          !-----------------------------------------------------------------
          SELECT CASE( name )

             ! These fields are always positive-definite
             CASE( 'QV2M',  'T2M', 'TS' )
                IF ( do2x25  ) WHERE( Q2x25 < 0e0 ) Q2x25 = 0e0
                IF ( do4x5   ) WHERE( Q4x5  < 0e0 ) Q4x5  = 0e0
             
             ! Divide winds by pressures
             CASE( 'U10M', 'V10M' )
                IF ( doNative ) Q     = Q     / P
                IF ( do2x25   ) Q2x25 = Q2x25 / P2x25
                IF ( do4x5    ) Q4x5  = Q4x5  / P4x5

             CASE DEFAULT
                ! Nothing

          END SELECT

          !-----------------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestCh, YNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Q( I0_eu:I1_eu, J0_eu:J1_eu )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestEu, YNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Q( I0_na:I1_na, J0_na:J1_na )
             st3d = (/ 1,       1,       H /)
             ct3d = (/ XNestNa, YNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name ), st3d, ct3d )
             NULLIFY( Ptr )
          ENDIF
          
          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st3d = (/ 1,     1,     H  /)
             ct3d = (/ X2x25, Y2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st3d, ct3d )
          ENDIF
          
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st3d = (/ 1,    1,    H /)
             ct3d = (/ X4x5, Y4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st3d, ct3d )
          ENDIF
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
    msg = '%%%%%% LEAVING ROUTINE Process2dSlvNx %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process2dSlvNx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process2dAlbedo
!
! !DESCRIPTION: Subroutine Process2dAlbedo creates the daily average albedo
!  field.  This routine is a wrapper for GeosFpProcessAlbedo.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process2dAlbedo()
!
! !REMARKS:
!   Rationale for doing this: 
!   ----------------------------------------------------------------------
!   The GEOS-FP ALBEDO field is only defined where it is daylight.
!   Some places in GEOS-Chem require an ALBEDO field even at night (i.e.
!   as a proxy for determining land surface).  Therefore compute the
!   daily average albedo and return to the data processing routine above.
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  19 Sep 2013 - R. Yantosca - Now save out nested Europe grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: H,        F,        hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,        T
    INTEGER                 :: XNestCh,  YNestCh,  TNestCh
    INTEGER                 :: XNestEu,  YNestEu,  TNestEu
    INTEGER                 :: XNestNa,  YNestNa,  TNestNa
    INTEGER                 :: X2x25,    Y2x25,    T2x25
    INTEGER                 :: X4x5,     Y4x5,     T4x5
    INTEGER                 :: ct3d(3),  st3d(3)

    ! Data arrays
    REAL*4, TARGET          :: Q    ( I025x03125, J025x03125, TIMES_A1 )
    REAL*4                  :: Q2x25( I2x25,      J2x25                )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5                 )

    ! Pointers
    REAL*4, POINTER         :: ptr(:,:)

     ! Character strings and arrays
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process2dAlbedo %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
    ENDIF

    ! Nested Europe grid
    IF ( doNestEu ) THEN
       CALL NcGet_DimLen( fOutNestEu, 'lon',  XNestEu )
       CALL NcGet_DimLen( fOutNestEu, 'lat',  YNestEu ) 
       CALL NcGet_DimLen( fOutNestEu, 'time', TNestEu )
    ENDIF

    ! Nested China grid
    IF ( doNestNa ) THEN
       CALL NcGet_DimLen( fOutNestNa, 'lon',  XNestNa )
       CALL NcGet_DimLen( fOutNestNa, 'lat',  YNestNa ) 
       CALL NcGet_DimLen( fOutNestNa, 'time', TNestNa )
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
    ! Read each hour of ALBEDO data and store into the Q array
    !=======================================================================

    ! Zero data arrays
    Q     = 0e0

    ! Loop over the number of files per day
    DO H = 1, TIMES_A1

       !-----------------------------------------------------------------
       ! Prepare for file input
       !-----------------------------------------------------------------

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a1Mins(H) / 60 ) * 10000 + 3000
       
       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg1_2d_rad_Nx_file )
       CALL expandDate( fNameInput, yyyymmdd, hhmmss )
       
       ! Echo info
       msg = '%%% Opening ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Open the netCDF4 file for input
       CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
       ! Get the dimensions from the netCDF file
       CALL NcGet_DimLen( fIn, 'lon',  X )
       CALL NcGet_DimLen( fIn, 'lat',  Y ) 
       CALL NcGet_DimLen( fIn, 'time', T )
              
       !-----------------------------------------------------------------
       ! Read data
       !-----------------------------------------------------------------
       msg = '%%% Reading ALBEDO'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       
       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st3d  = (/ 1, 1, 1 /)
       ct3d  = (/ X, Y, 1 /)

       ! Read data from file
       CALL NcRd( Q(:,:,H), fIn, 'ALBEDO', st3d, ct3d )
       
       ! Replace missing values with zeroes
       WHERE( Q == FILL_VALUE ) Q = 0e0

       !-----------------------------------------------------------------
       ! Close input file
       !-----------------------------------------------------------------
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )          

    ENDDO

    !=======================================================================
    ! Process ALBEDO data
    !=======================================================================

    ! Compute daily average albedo
    CALL GeosFpProcessAlbedo( Q )

    !-----------------------------------------------------------------
    ! Regrid to coarse resolution
    !-----------------------------------------------------------------
    msg = '%%% Regridding ALBEDO'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Regrid
    IF ( do2x25 ) CALL RegridGeosFpTo2x25( 0, Q(:,:,1), Q2x25 )
    IF ( do4x5  ) CALL RegridGeosFpTo4x5 ( 0, Q(:,:,1), Q4x5  )

    ! Make sure ALBEDO is positive-definite
    IF ( do2x25 ) WHERE( Q2x25 < 0e0 ) Q2x25 = 0e0
    IF ( do4x5  ) WHERE( Q4x5  < 0e0 ) Q4x5  = 0e0
   
    !-----------------------------------------------------------------
    ! Write average daily albedo to netCDF files
    !-----------------------------------------------------------------

    ! Save to disk
    msg = '%%% Archiving  ALBEDO'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Write the daily average albedo to disk
    DO H = 1, TIMES_A1
          
       ! Nested China (point to proper slice of global data)
       IF ( doNestCh ) THEN
          Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch, 1 )
          st3d = (/ 1,       1,       H /)
          ct3d = (/ XNestCh, YNestCh, 1 /)
          CALL NcWr( Ptr, fOutNestCh, 'ALBEDO', st3d, ct3d )
          NULLIFY( Ptr )
       ENDIF

       ! Nested EU (point to proper slice of global data)
       IF ( doNestEu ) THEN
          Ptr  => Q( I0_eu:I1_eu, J0_eu:J1_eu, 1 )
          st3d = (/ 1,       1,       H /)
          ct3d = (/ XNestEu, YNestEu, 1 /)
          CALL NcWr( Ptr, fOutNestEu, 'ALBEDO', st3d, ct3d )
          NULLIFY( Ptr )
       ENDIF

       ! Nested NA (point to proper slice of global data)
       IF ( doNestNa ) THEN
          Ptr  => Q( I0_na:I1_na, J0_na:J1_na, 1 )
          st3d = (/ 1,       1,       H /)
          ct3d = (/ XNestNa, YNestNa, 1 /)
          CALL NcWr( Ptr, fOutNestNa, 'ALBEDO', st3d, ct3d )
          NULLIFY( Ptr )
       ENDIF

       ! Write 2 x 2.5 data
       IF ( do2x25 ) THEN
          st3d = (/ 1,     1,     H  /)
          ct3d = (/ X2x25, Y2x25, 1  /)
          CALL NcWr( Q2x25, fOut2x25, 'ALBEDO', st3d, ct3d )
       ENDIF
       
       ! Write 4x5 data
       IF ( do4x5 ) THEN
          st3d = (/ 1,    1,    H /)
          ct3d = (/ X4x5, Y4x5, 1 /)
          CALL NcWr( Q4x5, fOut4x5, 'ALBEDO', st3d, ct3d )
       ENDIF
       
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================
    msg = '%%%%%% LEAVING ROUTINE Process2dAlbedo %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process2dAlbedo
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpSeaIceBins
!
! !DESCRIPTION: Subroutine GeosFpSeaIceBins bins the FRSEAICE field into 
!  bins for 2 x 2.5 and 4 x 5 output.  For each coarse grid box, the number of 
!  fine grid boxes having a sea ice fraction within a particular bin is
!  computed.  Typically the bins will be percentage decades (e.g. 0-10%, 
!  10-20%, 20-30%, etc.)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpSeaIceBins( Ice, BinSize, map, IceOut, IMX, JMX )
!
! !INPUT PARAMETERS: 
!
    ! Sea ice fraction (Nx grid)
    REAL*4,       INTENT(IN)        :: Ice(I025x03125,J025x03125)

    ! Size of each fractional sea ice bin
    REAL*4,       INTENT(IN)        :: BinSize

    ! Mapping weight object
    TYPE(MapObj), POINTER, OPTIONAL :: map(:,:)

    ! Dimensions of coarse grid
    INTEGER,      INTENT(IN)        :: IMX, JMX  
!
! !OUTPUT PARAMETERS:
!
    ! Binned sea ice fraction, output grid
    REAL*4,       INTENT(OUT)       :: IceOut(IMX,JMX,N_ICE)
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  17 Aug 2010 - R. Yantosca - Initial version, based on RegridTau
!  25 Aug 2010 - R. Yantosca - Renamed to "GeosFpSeaIceBins"
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    LOGICAL :: isNative
    INTEGER :: B, I, J, T, nPoints, Nx, Ny, X, Y 
    REAL*4  :: sum_Wn, weight

    ! Zero output variable
    IceOut = 0e0

    ! If we are processing at native resolution, then 
    ! we don't need to pass the mapping weight object.
    isNative = ( .not. PRESENT( map ) )

    ! Loop over coarse grid boxes
    DO J = 1, JMX
    DO I = 1, IMX

       ! Number of "fine" grid boxes in each dimension
       ! that comprise a "coarse" grid box
       nPoints = map(I,J)%nPoints

       !---------------------------------------------------------------
       ! Place fractional sea ice data into N_ICE bins
       !---------------------------------------------------------------

       ! Zero mapping variables
       sum_Wn = 0e0

       ! Loop over "fine" grid boxes
       DO Ny = 1, nPoints
       DO Nx = 1, nPoints
          
          ! Mapping weight (i.e. the fraction of each fine
          ! grid box that fits into the coarse grid box)
          weight = map(I,J)%weight(Nx,Ny)

          ! Avoid useless clock cycles if the mapping weight is zero
          IF ( weight > 0d0 ) THEN

             ! Indices of each "fine" grid box that makes up the "coarse" box
             X               = map(I,J)%xInd(Nx)
             Y               = map(I,J)%yInd(Ny)

             ! Sum of the mapping weights over all of the "fine" grid
             ! boxes (X,Y) that make up the "coarse" grid box (I,J)
             sum_Wn          = sum_Wn + map(I,J)%weight(Nx,Ny)

             ! Compute the bin number, based on the value of the 
             ! sea ice fraction on the fine "Nx" grid
             B               = INT( Ice(X,Y) / BinSize ) + 1

             ! Make sure B lies in the range 1...N_ICE
             B               = MAX( MIN( B, N_ICE ), 1 )
       
             ! Add the number of fine boxes having the particular sea ice 
             ! fraction to each bin B.  We just need to add the mapping 
             ! weight, which accounts for the fraction of the fine box 
             ! (X,Y) that is located inside the coarse box (I,J).
             IceOut(I,J,B) = IceOut(I,J,B) + map(I,J)%weight(Nx,Ny)
          ENDIF
       ENDDO
       ENDDO

       !---------------------------------------------------------------
       ! Compute fractional sea ice coverage in each bin
       !---------------------------------------------------------------

       ! Normalize each fractional sea ice bin by the total
       ! # of fine boxes that fit into the coarse box
       DO B = 1, N_ICE
          IceOut(I,J,B) = IceOut(I,J,B) / sum_Wn 
       ENDDO

       ! Safety check!  The sum of all bins should add up to 1, 
       ! within roundoff tolerance of 1e-4.
       IF ( ABS( 1e0 - SUM( IceOut(I,J,:) ) ) >= 1e-4 ) THEN
          WRITE( 6, '(a)' ) 'SEA ICE BINS DO NOT ADD UP TO 1!'
          WRITE( 6, 100 ) I, J, T, SUM( IceOut(I,J,:) )
100       FORMAT( 'I, J, T, SUM: ', 3i4, 1x, f13.7 )
          CALL EXIT(1)
       ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE GeosFpSeaIceBins
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpCreateLwi
!
! !DESCRIPTION: Subroutine GeosFpCreateLwi creates the GEOS-5 style 
!  land/water/ice (LWI) flags field and then regrids it to coarser resolution. 
!  LWI is used for backwards compatibility w/ existing GEOS-Chem routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpCreateLwi( frSeaIce, map, IMX, JMX, lwiOut )
!
! !INPUT PARAMETERS:
!
    ! Sea ice fraction, from the tavg1_2d_flx_Nx file
    REAL*4,       INTENT(IN)  :: frSeaIce(I025x03125,J025x03125)

    ! Object containing mapping weights
    TYPE(MapObj), POINTER     :: map(:,:)

    ! Dimensions of output array LWIOUT
    INTEGER,      INTENT(IN)  :: IMX, JMX                     
!
! !OUTPUT PARAMETERS:
!
    ! Regridded land-water indices on the output grid
    REAL*4,       INTENT(OUT) :: lwiOut(IMX,JMX)
!
! !REMARKS:
!  LWI = 0 are ocean boxes
!  LWI = 1 are land or land-ice boxes
!  LWI = 2 are sea ice boxes
!
! !REVISION HISTORY: 
!  06 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  06 Jan 2012 - R. Yantosca - Skip regridding for native resolution
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J
    REAL*4  :: lwiIn(I025x03125,J025x03125)

    ! Loop over # of A1 times
    ! Intitialize the LWI array w/ the time-invariant part
    lwiIn = lwiMask

    ! Also factor in the time-varying sea ice
    DO J = 1, J025x03125
    DO I = 1, I025x03125
       IF ( frSeaIce(I,J) > 0.5e0 .and. lwiIn(I,J) < 1e0 ) THEN
          lwiIn(I,J) = 2e0
       ENDIF
    ENDDO
    ENDDO

    IF ( IMX == I025x03125 .and. JMX == J025x03125 ) THEN

       ! Skip regridding if the output grid is the native grid.
       ! Simply assign lwiOut = lwiIn and return.
       lwiOut = lwiIn

    ELSE

       ! If the MAP object is passed, then 
       ! regrid the LWI field to coarse resolution
       CALL GeosFpRegridLwi( lwiIn, lwiOut, map, IMX, JMX )
 
    ENDIF
       
  END SUBROUTINE GeosFpCreateLwi
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpRegridLwi
!
! !DESCRIPTION: This routine regrids the land-water indices (LWI) field.  
!  Instead of an actual area regridding, we pick the mode of the LWI values
!  of each 0.5 x 0.667 box that fits into a coarse grid box.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpRegridLwi( lwiIn, lwiOut, map, IMX, JMX )
!
! !INPUT PARAMETERS:
!
    ! Land-water indices on the 0.5. x 0.666 grid
    REAL*4,       INTENT(IN)  :: lwiIn(I025x03125,J025x03125)

    ! Object that contains the mapping weights to the output grid
    TYPE(MapObj), POINTER     :: map(:,:)

    ! Dimensions of the output grid
    INTEGER,      INTENT(IN)  :: IMX, JMX
!
! !OUTPUT PARAMETERS:
!
    ! Regridded land-water indices on the output grid
    REAL*4,       INTENT(OUT) :: lwiOut(IMX,JMX)
!
! !REVISION HISTORY:
!  06 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    INTEGER :: I, J, X, Y, Nx, Ny, nPoints, index
    INTEGER :: hist(0:2)
    REAL*4  :: mode(1)

    ! Loop over grid boxes
    DO J = 1, JMX
    DO I = 1, IMX

       ! Number of "fine" grid boxes in each dimension
       ! that comprise a "coarse" grid box
       nPoints = map(I,J)%nPoints

       ! Zero the histogram array
       hist = 0
       
       ! Loop over "fine" grid boxes
       DO Ny = 1, nPoints
       DO Nx = 1, nPoints
          
          ! Avoid useless clock cycles if the mapping weight is zero
          IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN

             ! Indices of each "fine" grid box that makes up the "coarse" box
             X             = map(I,J)%xInd(Nx)
             Y             = map(I,J)%yInd(Ny)

             ! Sort each LWI value on the "fine" grid into a histogram
             ! Possible values of LWI are 0, 1, 2
             index       = INT( lwiIn(X,Y) ) 
             hist(index) = hist(index) + 1 
          ENDIF
       ENDDO
       ENDDO

       ! The bin in the histogram w/ the most counts is the MODE,
       ! so we'll use that as the regridded value of LWI.
       !
       ! NOTE: The result from MAXLOC is indexed starting from 1 and not
       ! from zero...so we have to subtract 1 and then save to LWIOUT.
       !
       ! ALSO NOTE: If two or more elements of HIST have the same value,
       ! then MAXLOC will pick the one that comes first in array order. 
       ! So if a coarse box is 50% water and 50% ice, this regridding scheme 
       ! will label the box as water. (LWIOUT(I,J)=0).  For 50% land and 50% 
       ! ice, the coarse box will be labeled as land (LWIOUT(I,J)=1). 
       mode        = MAXLOC( hist )
       lwiOut(I,J) = mode(1) - 1
    ENDDO
    ENDDO

  END SUBROUTINE GeosFpRegridLwi
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpAdjustSnomas
!
! !DESCRIPTION: Routine GeosFpAdjustSnomas will adjust the GEOS-FP SNOMAS 
!  field to make it more similar to the GEOS-5 SNOMAS field.  This is necessary
!  for backward compatibility with existing GEOS-Chem routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpAdjustSnomas( Q )
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4, INTENT(INOUT) :: Q(I025x03125,J025x03125)  ! Raw SNOMAS data
!
! !REMARKS:
!  NOTE: This routine was originally developed for the MERRA data processing
!  code, but the same algorithm also can be used for GEOS-FP data.
!                                                                             .
!  Original comments from MERRA code:
!  ----------------------------------
!  The SNOMAS field in MERRA differs from that in the GEOS-5 ops data:
!                                                                             .
!  From the GEOS-5 File Specification Document:
!                                                                             .
!     SNOMAS: The mass of snow in per unit of land area in meters of 
!     liquid-water-equivalent depth (i.e., 10^3 kg/m2). In grid boxes 
!     with no land (FRLAND+FRLANDICE=0) it is set to _FillValue (= 1e15). 
!     Where FRLANDICE>0.9 it is arbitrarily set to 4 meters. Over other 
!     land areas it represents an average over the non-glaciated part.
!                                                                             .
!  From the MERRA File Specification Document:
!                                                                             .
!     SNOMAS: The mass of snow per unit of ice-free land area (FRLAND), 
!     in kg/m2. In grid boxes with no land it is set to _FillValue (=1e15). 
!     Over other land areas it represents an average over the nonglaciated
!     part.
!                                                                             .
!  Max Suarez (Max.J.Suarez@nasa.gov) clarifies this difference:
!                                                                             .
!     Early versions of GEOS had been writing SNOMAS in meters.  This was 
!     changed to mm (or kg/m^2) in  all recent versions, including 5_2, 
!     MERRA, and the current development tags.  But the forward processing 
!     spec was not updated. The MERRA spec, however, is correct.
!                                                                             .
!     To further complicate matters, the variable called SNOMAS in MERRA 
!     comes from a very different part of the code that in 5_x.  It is 
!     in a land collection intended to have representative values over 
!     the ice-free land portion of the grid box.  This applies to all 
!     variables in that collection.  SNOMAS in particular makes this 
!     clear in the glossary definition in the spec.
!                                                                             .
!     Forward processing (FP) puts out a grid averaged SNOMAS, including
!     ice-covered areas, where the "SNOMAS" was arbitrarily set to 4000 mm. 
!     Neither FP nor MERRA includes ocean or freshwater regions with snow 
!     over ice.
!                                                                             .
!     To answer your question, to convert MERRA SNOMAS to 5.2 SNOMAS,
!     use this equation:
!                                                                             .
!     SNOMAS_5.X = ( SNOMAS_merra * FRLAND_merra + 
!                    4000         * FRLANDICE_merra ) /
!                  ( FRLAND_merra + FRLANDICE_merra )
!                                                                             .
!     Sorry about the confusion, but it seemed silly in MERRA to continue
!     writing an invented value over glaciers. The next FP system will be
!     like MERRA in this regard.  In the longer term, we plan to have a 
!     better snow/ice parameterization over permanent glaciers.   
!                                                                             .
!  Therefore, we shall implement the algorithm that Max Suarez described
!  above.  This will make the output SNOMAS field similar to GEOS-5,
!  which will allow better backward compatibility w/ existing code.
!                                                                             .
!  Also note: Liquid water equivalent height is defined as such:
!     1 m H2O = 10^3 kg/m2   ==>  10^-3 m H2O = 1 mm H2O = 1 kg/m2. 
!
! !REVISION HISTORY: 
!  06 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!    
    INTEGER :: I, J
    REAL*4  :: den

    ! Loop over surface grid boxes
    DO J = 1, J025x03125
    DO I = 1, I025x03125

       ! The denominator is the sum of the land ice and land fractions in 
       ! he grid box.  Nonzero denotes that we are over land and not ocean
       den = frLand(I,J) + frLandIce(I,J)

       ! Test if the division is possilble
       IF ( den > 0e0 ) THEN

          ! If so, then compute the SNOMAS value according to the
          ! algorithm described above
          Q(I,J) = (  Q(I,J) * frLand(I,J) + 4.0e3 * frLandIce(I,J) ) / den

       ELSE

          ! Otherwise, then we are over the ocean, so SNOMAS = 0
          Q(I,J) = 0e0
          
       ENDIF

   ENDDO
   ENDDO

  END SUBROUTINE GeosFpAdjustSnomas
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpProcessAlbedo
!
! !DESCRIPTION: Subroutine GeosFpProcessAlbedo computes the daily average
!  albedo from the MERRA raw data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpProcessAlbedo( Q )
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Input:  surface albedo [unitless] for each hour of the day
    ! Output: daily average surface albedo [unitless]
    REAL*4, INTENT(INOUT) :: Q( I025x03125, J025x03125, TIMES_A1 ) 
!
! !REMARKS:
!   Rationale for doing this: 
!   ----------------------------------------------------------------------
!   The GEOS-FP ALBEDO field is only defined where it is daylight.
!   Some places in GEOS-Chem require an ALBEDO field even at night (i.e.
!   as a proxy for determining land surface).  Therefore compute the
!   daily average albedo and return to the data processing routine above.
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  09 Jan 2012 - R. Yantosca - Now return avg albedo in 1st slot of Q
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: I, J, T

    ! Arrays
    ! For albedo processing
    REAL*4  :: A ( I025x03125, J025x03125 )
    INTEGER :: Ct( I025x03125, J025x03125 )

    ! Initialization
    A  = 0e0
    Ct = 0

    ! Sum up albedo over the entire day
    DO T = 1, TIMES_A1
    DO J = 1, J025x03125
    DO I = 1, I025x03125
       IF ( Q(I,J,T) > 0e0 ) THEN
          A (I,J) = A (I,J) + Q(I,J,T)
          Ct(I,J) = Ct(I,J) + 1
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    
    ! Zero data in Q array
    Q = 0e0

    ! Compute average albedo, store in 1st slot of Q
    ! Assign 0.85 for snow/ice over poles 
    DO J = 1, J025x03125
    DO I = 1, I025x03125
       IF ( ct(I,J) > 0 ) THEN
          Q(I,J,1) = A(I,J) / REAL( ct(I,J) )
       ELSE
          Q(I,J,1) = 0.85e0
       ENDIF
    ENDDO
    ENDDO

  END SUBROUTINE GeosFpProcessAlbedo
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpProcessTropp
!
! !DESCRIPTION: Subroutine "GeosProcessTropp" replaces any missing values in 
!  the TROPP field with a zonal mean average.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpProcessTropp( Q )
!
! !INPUT/OUTPUT PARAMETERS: 
!
    REAL*4, INTENT(INOUT) :: Q(:,:)  ! TROPP [hPa]
!
! !REMARKS:
!   Rationale for doing this: 
!   ----------------------------------------------------------------------
!   Sometimes the TROPP field has missing values, so we need to replace 
!   those with the average of the other boxes in the same latitude.  
!   Otherwise those missing values will get reset to zeroes (in routine
!   Geos5MakeA1Files), and those zeroes will propagate through the 
!   regridding process.  
!                                                                             .
!   If a box has a zero TROPP value, then when it is regridded to a 
!   coarser resolution, the resultant TROPP values will not be realistic.  
!   This will cause GEOS-Chem to diagnose the tropopause as much higher 
!   than it should.  Resetting the TROPP values according to the algorithm
!   below will avoid this problem.
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  09 Jan 2012 - R. Yantosca - Now make Q a 2-D array
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    INTEGER                 :: IX, JX, TX
    INTEGER                 :: I,  I2, J,  J2, T
    REAL*4                  :: tot, ct
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Geos5ProcessTropp begins here!
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpProcessTropp %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Test if missing values are found
    IF ( ANY( Q == FILL_VALUE ) ) THEN 

       ! If yes, then echo a message
       msg = '%%% Missing data values found in TROPP!  Removing these ...'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ELSE

       ! If no, echo a message, then exit
       msg = '%%% No missing data values found in TROPP!  Continuing on ...'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpProcessTropp %%%%%%%%%%'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       WRITE( IU_LOG, '(a)' ) '%%%'
       RETURN

    ENDIF

    ! Dimensions of the array
    IX = SIZE( Q, 1 )
    JX = SIZE( Q, 2 )
    
    ! Loop over grid boxes
    DO J = 1, JX
    DO I = 1, IX

       ! Replace "missing" values with a zonal average pressure
       IF ( Q(I,J) == FILL_VALUE ) THEN
             
          ! Zero summing variables
          tot = 0e0
          ct  = 0e0
             
          ! Sum up "good" boxes at this latitude
          DO I2 = 1, IX
             IF ( Q(I2,J) < FILL_VALUE ) THEN
                tot = tot + Q(I2,J)
                ct  = ct  + 1e0
             ENDIF
          ENDDO
          
          ! Avoid div by zero
          IF ( ct > 0e0 ) THEN 

             ! Replace "bad" value with zonal mean of "good" values
             Q(I,J) = tot / ct 
             
          ELSE
                
             ! If ct==0 then we have no good data at this latitude
             IF ( J > JX/2 ) THEN
                
                ! Northern hemisphere
                ! Then search down until the next good value
                DO J2 = J, 1, -1 
                   IF ( Q(I,J2) < FILL_VALUE ) THEN
                      Q(I,J) = Q(I,J2)
                   ENDIF
                ENDDO
                
             ELSE

                ! Southern hemisphere
                ! Then search up until the next good value
                DO J2 = 1, J
                   IF ( Q(I,J2) < FILL_VALUE ) THEN
                      Q(I,J) = Q(I,J2)
                   ENDIF
                ENDDO
                
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    ENDDO

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpProcessTropp %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    WRITE( IU_LOG, '(a)' ) '%%%'

  END SUBROUTINE GeosFpProcessTropp
!EOC
END MODULE GeosFpA1Module

