!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57A1Module
!
! !DESCRIPTION: Module Geos57A1Module contains routines to create the 
!  GEOS-Chem average 1-hr data files from the MERRA raw data.
!\\
!\\
! !INTERFACE: 

MODULE Geos57A1Module
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
  PUBLIC  :: Geos57MakeA1
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
!  PRIVATE :: Process2dFlxNx
!  PRIVATE :: Process2dLndNx
!  PRIVATE :: Process2dRadNx
!  PRIVATE :: Process2dSlvNx
!  PRIVATE :: Geos57SeaIceBins
!  PRIVATE :: Geos57CreateLwi
!  PRIVATE :: Geos57RegridLwi
!  PRIVATE :: Geos57AdjustSnomas
!  PRIVATE :: Geos57ProcessAlbedo
!  PRIVATE :: Geos57ProcessTropp
!
! !PRIVATE TYPES:
!
  REAL*4, TARGET, ALLOCATABLE :: Q2x25  (:,:,:,:)     ! 2 x 2.5 output array
  REAL*4, TARGET, ALLOCATABLE :: Q4x5   (:,:,:,:)     ! 4 x 5   output array
  REAL*4, TARGET, ALLOCATABLE :: Ice2x25(:,:,:,:)     ! 2 x 2.5 sea ice bins
  REAL*4, TARGET, ALLOCATABLE :: Ice4x5 (:,:,:,:)     ! 4 x 5   sea ice bins
  REAL*4, TARGET, ALLOCATABLE :: Lwi2x25(:,:,:  )     ! 2 x 2.5 land/water/ice
  REAL*4, TARGET, ALLOCATABLE :: Lwi4x5 (:,:,:  )     ! 4 x 5   land/water/ice
!
! !DEFINED_PARAMETERS:
!
  INTEGER, PARAMETER          :: N_ICE   = 10         ! # of sea ice bins
  REAL*4,  PARAMETER          :: BINSIZE = 1e0/N_ICE  ! Sea ice bin size
!
! !REVISION HISTORY:
!  05 Jan 2012 - R. Yantosca - Initial version, based on MERRA
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
!  05 Jan 2012 - R. Yantosca - Initial version, based on Geos57CnModule
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
    lName = 'GEOS-5.7.2 1-hour time-averaged (A1) fields for GEOS-Chem'
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

    ! Start Date
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',  yyyymmdd_string )

    ! Start Time
    lName = '00:00:00.00000'
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',  TRIM( lName )   )

    ! End Date
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',    yyyymmdd_string )

    ! End Time
    lName = '23:59:59.99999'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',    TRIM( lName )   )

    ! Delta-time
    lName = '10000'
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
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'ALBEDO', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! CLDTOT
    IF ( StrPos( 'CLDTOT', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Total cloud fraction' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'CLDTOT', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
    ENDIF

    ! FRSEAICE
    IF ( StrPos( 'FRSEAICE', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of sea ice on surface' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRSEAICE', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FRSNO
    IF ( StrPos( 'FRSNO', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fractional snow-covered area' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'FRSNO', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! GRN
    IF ( StrPos( 'GRN', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Vegetation greenness fraction' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GRN', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! GWETROOT
    IF ( StrPos( 'GWETROOT', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Root zone soil wetness' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GWETROOT', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! GWETTOP
    IF ( StrPos( 'GWETTOP', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Top soil wetness' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'GWETTOP', NF_FLOAT, 3, var3, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
    ENDIF

    ! LWI
    !IF ( StrPos( 'LWI', tavg1_2d__Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Land/water/ice flags' 
       units = 'unitless'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LWI', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    !ENDIF

    ! LWTUP
    IF ( StrPos( 'LWTUP', tavg1_2d_rad_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Upward longwave flux at top of atmosphere (TOA)' 
       units = ''
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'LWTUP', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
    ENDIF

    ! PRECTOT
    IF ( StrPos( 'PRECTOT', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Total surface precipitation flux' 
       units = 'kg m-2 s-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PRECTOT', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! QV2M
    IF ( StrPos( 'QV2M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Specific humidity at 2m above the displacement height' 
       units = 'kg kg-1'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'QV2M', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    IF ( StrPos( 'FRSEAICE', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN

       ! SEAICE00
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 0-10% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE00', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )

       ! SEAICE10
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 10-20% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE10', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       
       ! SEAICE20
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 20-30% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE20', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )

       ! SEAICE30
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 30-40% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE30', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )

       ! SEAICE40
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 40-50% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE40', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       
       ! SEAICE50
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 50-60% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE50', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       
       ! SEAICE60
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 60-70% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE60', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       
       ! SEAICE70
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 70-80% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE70', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       
       ! SEAICE80
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 80-90% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE80', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )

       ! SEAICE90
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of grid box that has 90-100% sea ice coverage' 
       units = 'fraction'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SEAICE90', NF_FLOAT, 3, var3, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )

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
    ENDIF

    ! SNOMAS
    IF ( StrPos( '', tavg1_2d_lnd_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Snow mass' 
       units = 'kg m-2'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SNOMAS', NF_FLOAT, 3, var3, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
    ENDIF

    ! T2M
    IF ( StrPos( 'T2M', tavg1_2d_slv_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Temperature 2m above displacement height' 
       units = 'K'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'SLV', NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
    ENDIF

    ! USTAR
    IF ( StrPos( 'USTAR', tavg1_2d_flx_Nx_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Friction velocity' 
       units = 'm -s'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'USTAR', NF_FLOAT, 3, var3, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
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
! !IROUTINE: Geos57MakeA1
!
! !DESCRIPTION: Routine Geos57MakeA1
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (surface values) from 
!       the MERRA raw data files (HDF4-EOS format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data in a format that GEOS-Chem can read.
! \end{enumerate}
! This routine is called directly from the main program Geos57Driver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57MakeA1
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!  17 Aug 2010 - R. Yantosca - Added Ice2x25, Ice4x5 arrays for computing
!                              binned fractional sea ice fields
!  24 Aug 2010 - R. Yantosca - Now also construct land/water/ice flags field
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
    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57MakeA1 %%%%%%%%%%'
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

    ! Allocate module arrays on 2 x 2.5 grid
    IF ( do2x25 ) THEN 
       ALLOCATE( Q2x25  ( I2x25, J2x25, TIMES_A1, nAllFields ) )
       ALLOCATE( Ice2x25( I2x25, J2x25, TIMES_A1, N_ICE      ) )
       ALLOCATE( Lwi2x25( I2x25, J2x25, TIMES_A1             ) )
       Q2x25   = 0e0
       Ice2x25 = 0e0
       Lwi2x25 = 0e0
    ENDIF

    ! Allocate 4x5 module arrays
    IF ( do4x5 ) THEN
       ALLOCATE( Q4x5  ( I4x5, J4x5, TIMES_A1, nAllFields ) )
       ALLOCATE( Ice4x5( I4x5, J4x5, TIMES_A1, N_ICE      ) ) 
       ALLOCATE( Lwi4x5( I4x5, J4x5, TIMES_A1             ) )
       Q4x5   = 0e0
       Ice4x5 = 0e0
       Lwi4x5 = 0e0
    ENDIF

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested China output file
    IF ( doNestCh ) THEN
       fName = dataTmplNestCh
       gName = 'SEA4CRS'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%',     'a1'         )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  TIMES_A1,    &
                          xMid_025x03125(I0_ch:I1_ch),       &
                          yMid_025x03125(J0_ch:J1_ch),       &
                          a1Mins,    gName,     fName,       &
                          fOutNestCh                        )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = dataTmpl2x25
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%',      'a1'        )
       CALL NcOutFileDef( I2x25,     J2x25,     TIMES_A1,    &
                          xMid_2x25, yMid_2x25, a1Mins,      &
                          gName,     fName,     fOut2x25    )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = dataTmpl4x5
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000      )      
       CALL StrRepl     ( fName,     '%%',     'a1'         )
       CALL NcOutFileDef( I4x5,      J4x5,      TIMES_A1,    &
                          xMid_4x5,  yMid_4x5,  a1Mins,      &
                          gName,     fName,     fOut4x5     )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
!    CALL Process2dFlxNx( nFields_2dFlxNx, fields_2dFlxNx, &
!                         foffset )
!    CALL Process2dLndNx( nFields_2dLndNx, fields_2dLndNx, offset )
!    CALL Process2dRadNx( nFields_2dRadNx, fields_2dRadNx, offset )
!    CALL Process2dSlvNx( nFields_2dSlvNx, fields_2dSlvNx, offset )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A1 output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57MakeA1 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Geos57MakeA1
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process2dFlxNx
!!
!! !DESCRIPTION: Subroutine Process2dFlxFx regrids the MERRA met fields from 
!!  the "tavg1\_2d\_flx\_Nx" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process2dFlxNx( nFields, fields, offset )
!!
!! !INPUT PARAMETERS:
!!
!    INTEGER,          INTENT(IN)    :: nFields     ! # of fields to process
!    CHARACTER(LEN=*), INTENT(IN)    :: fields(:)   ! List of field names
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,          INTENT(INOUT) :: offset      ! Offset for output arrays
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!!  12 Aug 2010 - R. Yantosca - Now use separate arrays for archiving PRECTOT
!!  17 Aug 2010 - R. Yantosca - Now call FracSeaIceBins to compute binned
!!                              fractional sea ice fields
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER                 :: T, F
!
!    ! Arrays
!    REAL*4, TARGET          :: Q     (I05x0666,J05x0666,TIMES_A1)
!    REAL*4                  :: PrcAnv(I05x0666,J05x0666,TIMES_A1)
!    REAL*4                  :: PrcCon(I05x0666,J05x0666,TIMES_A1)
!    REAL*4                  :: PrcLsc(I05x0666,J05x0666,TIMES_A1)
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_05x0666(:,:)
!    REAL*4, POINTER         :: ptr_2x25   (:,:)
!    REAL*4, POINTER         :: ptr_4x5    (:,:)
!
!     ! Character strings and arrays
!    CHARACTER(LEN=8       ) :: name
!    CHARACTER(LEN=MAX_CHAR) :: fileHDF
!    CHARACTER(LEN=MAX_CHAR) :: fileName
!    CHARACTER(LEN=MAX_CHAR) :: file2x25
!    CHARACTER(LEN=MAX_CHAR) :: file4x5
!    CHARACTER(LEN=MAX_CHAR) :: msg
!
!    !=======================================================================
!    ! Initialization
!    !=======================================================================
!
!    ! Zero all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process2dFlxNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg1_2d_flx_Nx_file )
!    CALL expandDate( fileHDF, yyyymmdd, 000000 )
!
!    ! Echo info
!    msg = '%%% Reading ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Open the HDF4-EOS file for input
!    CALL He4SetVerbose( VERBOSE )
!    CALL He4GridOpen( fileHDF )
!    CALL He4GridGetDimInfo
!    CALL He4GridReadX
!    CALL He4GridReadY
!    CALL He4GetNymdNhms
!
!    !=======================================================================
!    ! Process data
!    !=======================================================================
!
!    ! Loop over data fields
!    DO F = 1, nFields
!
!       ! Save field name into an 8-char variable. 
!       ! This will truncate field names longer than 8 chars.
!       name = TRIM( fields(F) )
!
!       ! Skip if the fieldname is empty
!       IF ( name == '' .or. name == 'PS' ) CYCLE
!
!       !--------------------------------
!       ! Read data from HDF-EOS file
!       !--------------------------------
!       IF ( name == 'PRECTOT' ) THEN
!
!          ! PRECTOT isn't read from disk, but instead is the sum of the 
!          ! PRECANV + PRECCON + PRECLSC fields, which are archived below
!          ! We assume that PRECANV, PRECCON, PRECLSC are all read from
!          ! disk first (i.e. put PRECTOT in the data list after all three
!          ! of these fields).
!          msg = '%%% Regridding   ' // name
!          WRITE( IU_LOG, '(a)' ) TRIM( msg )
!          Q  = PrcAnv + PrcCon + PrcLsc
!
!       ELSE 
!
!          ! Zero data aray
!          Q = 0e0
!          
!          ! Read data from HDF file
!          msg = '%%% Reading      ' // name
!          WRITE( IU_LOG, '(a)' ) TRIM( msg )
!          CALL He4GridReadData( name, Q )
!
!          ! Replace fill values with zeroes
!          WHERE( Q == FILL_VALUE ) Q = 0e0
!    
!          !-----------------------------
!          ! Pre-regrid handling
!          !-----------------------------
!          SELECT CASE( name )
!
!             CASE( 'PRECANV' )
!                msg = '%%% Constructing PRECTOT from ' // name
!                WRITE( IU_LOG, '(a)' ) TRIM( msg )
!                PrcAnv = Q
!
!             CASE( 'PRECCON' )
!                msg = '%%% Constructing PRECTOT from ' // name
!                WRITE( IU_LOG, '(a)' ) TRIM( msg )
!                PrcCon = Q
!
!             CASE( 'PRECLSC' )
!                msg = '%%% Constructing PRECTOT from ' // name
!                WRITE( IU_LOG, '(a)' ) TRIM( msg )
!                PrcLsc = Q
!
!             CASE( 'FRSEAICE' )
!                msg = '%%% Computing fractional sea ice coverage and LWI flags'
!                WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!                IF ( do2x25 ) THEN 
!
!                   ! Create the 2 x 2.5 land/water/ice flags field
!                   CALL Geos57CreateLwi( Q, mapNxTo2x25, I2x25, J2x25, Lwi2x25 )
!
!                   ! Bin sea ice for 2 x 2.5 output 
!                   CALL Geos57SeaIceBins( Q,       BINSIZE, mapNxTo2x25,  &
!                                         Ice2x25, I2x25,   J2x25        )
!                ENDIF
!
!                IF ( do4x5 ) THEN
!
!                   ! Create the 4 x 5 land/water/ice flags field
!                   CALL Geos57CreateLwi( Q, mapNxTo4x5, I4x5, J4x5, Lwi4x5 )
!
!                   ! Bin sea ice for 4x5 output
!                   CALL Geos57SeaIceBins( Q,       BINSIZE, mapNxTo4x5,   &
!                                         Ice4x5,  I4x5,    J4x5         )
!                ENDIF
!
!             CASE DEFAULT
!                ! Do Nothing
!          END SELECT
!       ENDIF
!
!       !--------------------------------
!       ! Regrid to 2 x 2.5 &  4 x 5
!       !--------------------------------
!       msg = '%%% Regridding   ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over the A1 times
!       !$OMP PARALLEL DO                                 &
!       !$OMP DEFAULT( SHARED )                           &
!       !$OMP PRIVATE( T, ptr_05x0666, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A1
!
!          !-----------------------------
!          ! Do the regridding!
!          !-----------------------------
!
!          ! Point to time-slice of native data array
!          ptr_05x0666 => Q(:,:,T)
!
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,T,F+offset)
!             CALL RegridGeos57NTo2x25( 0, ptr_05x0666, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,T,F+offset)
!             CALL RegridGeos57NTo4x5( 0, ptr_05x0666, ptr_4x5 )
!          ENDIF
!
!          !-----------------------------
!          ! Post-regrid handling
!          !-----------------------------
!          SELECT CASE( name )
!
!             ! These fields are always positive-definite
!             CASE( 'PRECANV', 'PRECCON', 'PRECLSC', 'PRECTOT', 'USTAR' )
!                IF ( do2x25 ) WHERE( ptr_2x25 < 0e0 ) ptr_2x25 = 0e0
!                IF ( do4x5  ) WHERE( ptr_4x5  < 0e0 ) ptr_4x5  = 0e0
!             CASE DEFAULT
!                ! Do Nothing
!          END SELECT
!       ENDDO
!       !!$OMP END PARALLEL DO
!    ENDDO
!          
!    !=======================================================================
!    ! Cleanup & quit
!    !=======================================================================
!
!    ! Increment offset for next routine
!    offset = offset + nFields
!
!    ! Nullify all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!  
!    msg = '%%%%%% LEAVING ROUTINE Process2dFlxNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process2dFlxNx
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process2dLndNx
!!
!! !DESCRIPTION: Subroutine Process2dLndFx regrids the MERRA met fields from 
!!  the "tavg1\_2d\_lnd\_Nx" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process2dLndNx( nFields, fields, offset )
!!
!! !INPUT PARAMETERS:
!!
!    INTEGER,          INTENT(IN)    :: nFields     ! # of fields to process
!    CHARACTER(LEN=*), INTENT(IN)    :: fields(:)   ! List of field names
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,          INTENT(INOUT) :: offset      ! Offset for output arrays
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!!  25 Aug 2010 - R. Yantosca - Now call Geos57AdjustSnomas to restore missing
!!                              SNOMAS data over the poles.  Use FRLANDICE
!!                              as a proxy for where to add the default values.
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER                 :: T, F
!
!    ! Arrays
!    REAL*4, TARGET          :: Q(I05x0666,J05x0666,TIMES_A1)
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_05x0666(:,:)
!    REAL*4, POINTER         :: ptr_2x25   (:,:)
!    REAL*4, POINTER         :: ptr_4x5    (:,:)
!
!     ! Character strings and arrays
!    CHARACTER(LEN=8       ) :: name
!    CHARACTER(LEN=MAX_CHAR) :: fileHDF
!    CHARACTER(LEN=MAX_CHAR) :: fileName
!    CHARACTER(LEN=MAX_CHAR) :: file2x25
!    CHARACTER(LEN=MAX_CHAR) :: file4x5
!    CHARACTER(LEN=MAX_CHAR) :: msg
!
!    !=======================================================================
!    ! Initialization
!    !=======================================================================
!
!    ! Zero all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process2dLndNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg1_2d_lnd_Nx_file )
!    CALL expandDate( fileHDF, yyyymmdd, 000000 )
!
!    ! Echo info
!    msg = '%%% Reading ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Open the HDF4-EOS file for input
!    CALL He4SetVerbose( VERBOSE )
!    CALL He4GridOpen( fileHDF )
!    CALL He4GridGetDimInfo
!    CALL He4GridReadX
!    CALL He4GridReadY
!    CALL He4GetNymdNhms
!
!    !=======================================================================
!    ! Process data
!    !=======================================================================
!
!    ! Loop over data fields
!    DO F = 1, nFields
!
!       ! Save field name into an 8-char variable. 
!       ! This will truncate field names longer than 8 chars.
!       name = TRIM( fields(F) )
!
!       ! Skip if the fieldname is empty
!       IF ( name == '' ) CYCLE
!
!       !-----------------------------------
!       ! Read data 
!       !-----------------------------------
!
!       ! Zero data array
!       Q = 0e0
!
!       ! Read data from HDF-EOS
!       msg = '%%% Reading    ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       CALL He4GridReadData( name, Q )
!
!       ! Remove fill values
!       WHERE( Q == FILL_VALUE ) Q = 0e0
!
!       !-----------------------------------
!       ! Regrid data to 2 x 2.5 & 4 x 5
!       !-----------------------------------
!       msg = '%%% Regridding ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over # of A-3 times
!       !$OMP PARALLEL DO                                 &
!       !$OMP DEFAULT( SHARED )                           &
!       !$OMP PRIVATE( T, ptr_05x0666, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A1
!
!          ! Point to time-slice of native data array
!          ptr_05x0666 => Q(:,:,T)
!
!          !-----------------------------
!          ! Pre-regrid handling
!          !-----------------------------
!          IF ( name == 'SNOMAS' ) THEN
!             CALL Geos57AdjustSnomas( ptr_05x0666 )  ! Restore data over poles
!          ENDIF
!         
!          !-----------------------------
!          ! Do the regridding!
!          !-----------------------------
!
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,T,F+offset)
!             CALL RegridGeos57NTo2x25( 0, ptr_05x0666, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,T,F+offset)
!             CALL RegridGeos57NTo4x5( 0, ptr_05x0666, ptr_4x5 )
!          ENDIF
! 
!          !-----------------------------
!          ! Post-regrid handling
!          !-----------------------------
!          SELECT CASE( name )
!             ! These fields are always positive-definite
!             CASE( 'FRSNO', 'GRN',   'GWETROOT', 'GWETTOP', 'LAI',   &
!                   'PARDF', 'PARDR', 'PRECSNO',  'SNOMAS',  'SNODP' )
!                IF ( do2x25 ) WHERE( ptr_2x25 < 0e0 ) ptr_2x25 = 0e0
!                IF ( do4x5  ) WHERE( ptr_4x5  < 0e0 ) ptr_4x5  = 0e0
!             CASE DEFAULT
!                ! Nothing
!          END SELECT
!       ENDDO
!       !$OMP END PARALLEL DO
!    ENDDO
!
!    !=======================================================================
!    ! Cleanup & quit
!    !=======================================================================
!
!    ! Increment offset
!    offset = offset + nFields
!
!    ! Nullify all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!
!    ! Echo info    
!    msg = '%%%%%% LEAVING ROUTINE Process2dLndNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process2dLndNx
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process2dRadNx
!!
!! !DESCRIPTION: Subroutine Process2dRadFx regrids the MERRA met fields from 
!!  the "tavg1\_2d\_rad\_Nx" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process2dRadNx( nFields, fields, offset )
!!
!! !INPUT PARAMETERS:
!!
!    INTEGER,          INTENT(IN)    :: nFields     ! # of fields to process
!    CHARACTER(LEN=*), INTENT(IN)    :: fields(:)   ! List of field names
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,          INTENT(INOUT) :: offset      ! Offset for output arrays
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER                 :: F, T
!
!    ! Arrays
!    REAL*4, TARGET          :: Q(I05x0666,J05x0666,TIMES_A1)
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_05x0666(:,:)
!    REAL*4, POINTER         :: ptr_2x25   (:,:)
!    REAL*4, POINTER         :: ptr_4x5    (:,:)
!
!     ! Character strings and arrays
!    CHARACTER(LEN=8       ) :: name
!    CHARACTER(LEN=MAX_CHAR) :: fileHDF
!    CHARACTER(LEN=MAX_CHAR) :: fileName
!    CHARACTER(LEN=MAX_CHAR) :: file2x25
!    CHARACTER(LEN=MAX_CHAR) :: file4x5
!    CHARACTER(LEN=MAX_CHAR) :: msg
!
!    !=======================================================================
!    ! Initialization
!    !=======================================================================
!
!    ! Zero all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process2dRadNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg1_2d_rad_Nx_file )
!    CALL expandDate( fileHDF, yyyymmdd, 000000 )
!
!    ! Echo info
!    msg = '%%% Reading ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Open the HDF4-EOS file for input
!    CALL He4SetVerbose( VERBOSE )
!    CALL He4GridOpen( fileHDF )
!    CALL He4GridGetDimInfo
!    CALL He4GridReadX
!    CALL He4GridReadY
!    CALL He4GetNymdNhms
!
!    !=======================================================================
!    ! Process data
!    !=======================================================================
!
!    ! Loop over data fields
!    DO F = 1, nFields
!
!       ! Save field name into an 8-char variable. 
!       ! This will truncate field names longer than 8 chars.
!       name = TRIM( fields(F) )
!
!       ! Skip if the field is empty
!       IF ( name == '' ) CYCLE
!
!       !--------------------------------
!       ! Read data 
!       !--------------------------------
!       Q = 0e0
!       msg = '%%% Reading    ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       CALL He4GridReadData( name, Q )
!          
!       ! Remove fill values
!       WHERE( Q == FILL_VALUE ) Q = 0e0
!
!       !--------------------------------
!       ! Special ALBEDO handling
!       !--------------------------------
!       IF ( name == 'ALBEDO' ) THEN
!          CALL Geos57ProcessAlbedo( Q ) 
!       ENDIF
!            
!       !--------------------------------
!       ! Regrid to 2 x 2.5 & 4 x 5
!       !--------------------------------
!       msg = '%%% Regridding ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over times
!       !$OMP PARALLEL DO                                 &
!       !$OMP DEFAULT( SHARED )                           &
!       !$OMP PRIVATE( T, ptr_05x0666, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A1
!
!          !-----------------------------
!          ! Do the regridding!
!          !-----------------------------
!
!          ! Point to time-slice of native data array
!          ptr_05x0666 => Q(:,:,T)
!
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,T,F+offset)
!             CALL RegridGeos57NTo2x25( 0, ptr_05x0666, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,T,F+offset)
!             CALL RegridGeos57NTo4x5( 0, ptr_05x0666, ptr_4x5 )
!          ENDIF
!
!          !-----------------------------
!          ! Post-regrid handling
!          !-----------------------------
!          SELECT CASE( name )
!             ! These fields should be positive-definite
!             CASE( 'ALBEDO', 'CLDTOT', 'LWTUP', 'SWGDN' )
!                IF ( do2x25 ) WHERE( ptr_2x25 < 0e0 ) ptr_2x25 = 0e0
!                IF ( do4x5  ) WHERE( ptr_4x5  < 0e0 ) ptr_4x5  = 0e0
!             CASE DEFAULT
!                ! Do Nothing
!          END SELECT
!       ENDDO
!       !$OMP END PARALLEL DO
!    ENDDO
!          
!    !=======================================================================
!    ! Cleanup & quit
!    !=======================================================================
!
!    ! Increment offset for next routine
!    offset = offset + nFields
!
!    ! Nullify all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!
!    ! Echo info    
!    msg = '%%%%%% LEAVING ROUTINE Process2dRadNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process2dRadNx
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process2dSlvNx
!!
!! !DESCRIPTION: Subroutine Process2dSlvFx regrids the MERRA met fields from 
!!  the "tavg1\_2d\_slv\_Nx" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process2dSlvNx( nFields, fields, offset )
!!
!! !INPUT PARAMETERS:
!!
!    INTEGER,          INTENT(IN)    :: nFields     ! # of fields to process
!    CHARACTER(LEN=*), INTENT(IN)    :: fields(:)   ! List of field names
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,          INTENT(INOUT) :: offset      ! Offset for output arrays
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER                 :: T, F
!
!    ! Arrays
!    REAL*4, TARGET          :: Q    (I05x0666,J05x0666,TIMES_A1)
!    REAL*4, TARGET          :: P    (I05x0666,J05x0666,TIMES_A1)
!    REAL*4, TARGET          :: P2x25(I2x25,   J2x25,   TIMES_A1)
!    REAL*4, TARGET          :: P4x5 (I4x5,    J4x5,    TIMES_A1)
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_05x0666(:,:)
!    REAL*4, POINTER         :: ptr_2x25   (:,:)
!    REAL*4, POINTER         :: ptr_4x5    (:,:)
!
!     ! Character strings and arrays
!    CHARACTER(LEN=8       ) :: name
!    CHARACTER(LEN=MAX_CHAR) :: fileHDF
!    CHARACTER(LEN=MAX_CHAR) :: fileName
!    CHARACTER(LEN=MAX_CHAR) :: file2x25
!    CHARACTER(LEN=MAX_CHAR) :: file4x5
!    CHARACTER(LEN=MAX_CHAR) :: msg
!
!    !=======================================================================
!    ! Initialization
!    !=======================================================================
!
!    ! Zero all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process2dSlvNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg1_2d_slv_Nx_file )
!    CALL expandDate( fileHDF, yyyymmdd, 000000 )
!
!    ! Echo info
!    msg = '%%% Reading ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Open the HDF4-EOS file for input
!    CALL He4SetVerbose( VERBOSE )
!    CALL He4GridOpen( fileHDF )
!    CALL He4GridGetDimInfo
!    CALL He4GridReadX
!    CALL He4GridReadY
!    CALL He4GetNymdNhms
!
!    !=======================================================================
!    ! Process surface pressure
!    !=======================================================================
!
!    ! Read surface pressure (all times)
!    msg = '%%% Reading      PS'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridReadData( 'PS', P )
!
!    ! Loop over times
!    !$OMP PARALLEL DO                               &
!    !$OMP DEFAULT( SHARED )                         &
!    !$OMP PRIVATE( T, ptr_05x0666, ptr_2x25, ptr_4x5 )
!    DO T = 1, TIMES_A1
!
!       ! Point to time-slice of native data array
!       ptr_05x0666 => P(:,:,T)
!
!       ! Regrid to 2 x 2.5
!       IF ( do2x25 ) THEN
!          ptr_2x25 => P2x25(:,:,T)
!          CALL RegridGeos57NTo2x25( 0, ptr_05x0666, ptr_2x25 )
!       ENDIF
!
!       ! Regrid to 4x5 
!       IF ( do4x5 ) THEN
!          ptr_4x5 => P4x5(:,:,T)
!          CALL RegridGeos57NTo4x5( 0, ptr_05x0666, ptr_4x5  )
!       ENDIF
!    ENDDO
!    !$OMP END PARALLEL DO
!
!    !=======================================================================
!    ! Process data
!    !=======================================================================
!
!    ! Loop over data fields
!    DO F = 1, nFields
!
!       ! Save field name into an 8-char variable. 
!       ! This will truncate field names longer than 8 chars.
!       name = TRIM( fields(F) )
!
!       !-----------------------------
!       ! Read data 
!       !-----------------------------
!
!       ! Zero data array
!       Q = 0e0
!
!       ! Read data from file
!       msg = '%%% Reading    ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       CALL He4GridReadData( name, Q )
!
!       ! Fill missing values in TROPP and other fields
!       SELECT CASE( name )
!          CASE( 'TROPPT', 'TROPPV', 'TROPPB' )
!             CALL Geos57ProcessTropp( Q )
!          CASE DEFAULT
!             WHERE( Q == FILL_VALUE ) Q = 0e0
!       END SELECT
!
!       !-----------------------------
!       ! Pre-regrid handling
!       !-----------------------------
!       SELECT CASE( name )
!          CASE( 'PS',     'SLP',              &
!                'TROPPT', 'TROPPV', 'TROPPB' )
!             Q = Q / 100e0                      ! Pa -> hPa 
!          CASE( 'U10M', 'V10M' )               
!             Q = Q * P                          ! Multiply winds by pressure
!          CASE DEFAULT
!             ! Nothing
!       END SELECT
!
!       !-----------------------------
!       ! Regrid data to 2x25, 4x5
!       !-----------------------------
!       msg = '%%% Regridding ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over A1 times
!       !$OMP PARALLEL DO                                 &
!       !$OMP DEFAULT( SHARED )                           &
!       !$OMP PRIVATE( T, ptr_05x0666, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A1
!
!          ! Point to time-slice of native data array
!          ptr_05x0666 => Q(:,:,T)
!
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,T,F+offset)
!             CALL RegridGeos57NTo2x25( 0, ptr_05x0666, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,T,F+offset)
!             CALL RegridGeos57NTo4x5( 0, ptr_05x0666, ptr_4x5 )
!          ENDIF
!
!          !---------------------------
!          ! Post-regrid handling
!          !---------------------------
!          SELECT CASE( name )
!
!             ! These fields are always positive-definite
!             CASE( 'QV2M',  'T2M', 'TS' )
!                IF ( do2x25  ) WHERE( ptr_2x25 < 0e0 ) ptr_2x25 = 0e0
!                IF ( do4x5   ) WHERE( ptr_4x5  < 0e0 ) ptr_4x5  = 0e0
!             
!             ! Divide winds by pressures
!             CASE( 'U10M', 'V10M' )
!                IF ( do2x25  ) ptr_2x25 = ptr_2x25 / P2x25(:,:,T)
!                IF ( do4x5   ) ptr_4x5  = ptr_4x5  / P4x5 (:,:,T)
!
!             CASE DEFAULT
!                ! Nothing
!          END SELECT
!
!       ENDDO
!       !$OMP END PARALLEL DO
!    ENDDO
!
!    !=======================================================================
!    ! Cleanup & quit
!    !=======================================================================
!
!    ! Increment offset for next routine
!    offset = offset + nFields
!
!    ! Nullify all pointers
!    NULLIFY( ptr_05x0666, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!
!    ! Echo info    
!    msg = '%%%%%% LEAVING ROUTINE Process2dSlvNx %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process2dSlvNx
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57SeaIceBins
!!
!! !DESCRIPTION: Subroutine Geos57SeaIceBins bins the FRSEAICE field into 
!!  bins for 2 x 2.5 and 4 x 5 output.  For each coarse grid box, the number of 
!!  fine grid boxes having a sea ice fraction within a particular bin is
!!  computed.  Typically the bins will be percentage decades (e.g. 0-10%, 
!!  10-20%, 20-30%, etc.)
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57SeaIceBins( IceNx, BinSize, map, IceOut, IMX, JMX )
!!
!! !INPUT PARAMETERS: 
!!
!    ! Sea ice fraction (Nx grid)
!    REAL*4,       INTENT(IN)  :: IceNx(I05x0666,J05x0666,TIMES_A1) 
!
!    ! Size of each fractional sea ice bin
!    REAL*4,       INTENT(IN)  :: BinSize
!
!    ! Mapping weight object
!    TYPE(MapObj), POINTER     :: map(:,:)
!
!    ! Dimensions of coarse grid
!    INTEGER,      INTENT(IN)  :: IMX, JMX  
!!
!! !OUTPUT PARAMETERS:
!!
!    ! Binned sea ice fraction, output grid
!    REAL*4,       INTENT(OUT) :: IceOut(IMX,JMX,TIMES_A1,N_ICE)
!!
!! !REMARKS:
!!
!! !REVISION HISTORY: 
!!  17 Aug 2010 - R. Yantosca - Initial version, based on RegridTau
!!  25 Aug 2010 - R. Yantosca - Renamed to "Geos57SeaIceBins"
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!
!    ! Local variables
!    INTEGER :: B, I, J, T, nPoints, Nx, Ny, X, Y 
!    REAL*4  :: sum_Wn
!
!    ! Zero output variable
!    IceOut = 0e0
!
!    ! Loop over coarse grid boxes
!    !$OMP PARALLEL DO  &
!    !$OMP DEFAULT( SHARED ) &
!    !$OMP PRIVATE( I, J, T, nPoints, Nx, Ny, X, Y, sum_Wn, B )
!    DO T = 1, TIMES_A1
!    DO J = 1, JMX
!    DO I = 1, IMX
!
!       ! Number of "fine" grid boxes in each dimension
!       ! that comprise a "coarse" grid box
!       nPoints = map(I,J)%nPoints
!
!       !---------------------------------------------------------------
!       ! Place fractional sea ice data into N_ICE bins
!       !---------------------------------------------------------------
!
!       ! Zero mapping variables
!       sum_Wn = 0e0
!
!       ! Loop over "fine" grid boxes
!       DO Ny = 1, nPoints
!       DO Nx = 1, nPoints
!          
!          ! Avoid useless clock cycles if the mapping weight is zero
!          IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN
!
!             ! Indices of each "fine" grid box that makes up the "coarse" box
!             X               = map(I,J)%xInd(Nx)
!             Y               = map(I,J)%yInd(Ny)
!
!             ! Sum of the mapping weights over all of the "fine" grid
!             ! boxes (X,Y) that make up the "coarse" grid box (I,J)
!             sum_Wn          = sum_Wn + map(I,J)%weight(Nx,Ny)
!
!             ! Compute the bin number, based on the value of the 
!             ! sea ice fraction on the fine "Nx" grid
!             B               = INT( IceNx(X,Y,T) / BinSize ) + 1
!
!             ! Make sure B lies in the range 1...N_ICE
!             B               = MAX( MIN( B, N_ICE ), 1 )
!       
!             ! Add the number of fine boxes having the particular sea ice 
!             ! fraction to each bin B.  We just need to add the mapping 
!             ! weight, which accounts for the fraction of the fine box 
!             ! (X,Y) that is located inside the coarse box (I,J).
!             IceOut(I,J,T,B) = IceOut(I,J,T,B) + map(I,J)%weight(Nx,Ny)
!          ENDIF
!       ENDDO
!       ENDDO
!
!       !---------------------------------------------------------------
!       ! Compute fractional sea ice coverage in each bin
!       !---------------------------------------------------------------
!
!       ! Normalize each fractional sea ice bin by the total
!       ! # of fine boxes that fit into the coarse box
!       DO B = 1, N_ICE
!          IceOut(I,J,T,B) = IceOut(I,J,T,B) / sum_Wn 
!       ENDDO
!
!       ! Safety check!  The sum of all bins should add up to 1, 
!       ! within roundoff tolerance of 1e-4.
!       IF ( ABS( 1e0 - SUM( IceOut(I,J,T,:) ) ) >= 1e-4 ) THEN
!          WRITE( 6, '(a)' ) 'SEA ICE BINS DO NOT ADD UP TO 1!'
!          WRITE( 6, 100 ) I, J, T, SUM( IceOut(I,J,T,:) )
!100       FORMAT( 'I, J, T, SUM: ', 3i4, 1x, f13.7 )
!          CALL EXIT(1)
!       ENDIF
!
!    ENDDO
!    ENDDO
!    ENDDO
!    !$OMP END PARALLEL DO
!
!  END SUBROUTINE Geos57SeaIceBins
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57CreateLwi
!!
!! !DESCRIPTION: Subroutine Geos57CreateLwi creates the GEOS-5 style 
!!  land/water/ice (LWI) flags field and then regrids it to coarser resolution. 
!!  LWI is used for backwards compatibility w/ existing GEOS-Chem routines.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57CreateLwi( frSeaIce, map, IMX, JMX, lwiOut )
!!
!! !INPUT PARAMETERS:
!!
!    ! Sea ice fraction, from the tavg1_2d_flx_Nx file
!    REAL*4,       INTENT(IN)  :: frSeaIce(I05x0666,J05x0666,TIMES_A1)  
!
!    ! Object containing mapping weights
!    TYPE(MapObj), POINTER     :: map(:,:)
!
!    ! Dimensions of output array LWIOUT
!    INTEGER,      INTENT(IN)  :: IMX, JMX                     
!!
!! !OUTPUT PARAMETERS:
!!
!    ! Regridded land-water indices on the output grid
!    REAL*4,       INTENT(OUT) :: lwiOut(IMX,JMX,TIMES_A1)
!!
!! !REMARKS:
!!  LWI = 0 are ocean boxes
!!  LWI = 1 are land or land-ice boxes
!!  LWI = 2 are sea ice boxes
!!
!! !REVISION HISTORY: 
!!  25 Aug 2010 - R. Yantosca - Initial version
!!  31 Aug 2010 - R. Yantosca - Now only assign sea ice (LWI=2) to boxes where
!!                              FRSEAICE > 0.5.  This will eliminate spurious
!!                              ice coverage (e.g. over Hudson Bay in summer).
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER        :: I, J, T
!    REAL*4, TARGET :: lwiIn(I05x0666,J05x0666,TIMES_A1)
!
!    ! Loop over # of A1 times
!    !$OMP PARALLEL DO        &
!    !$OMP DEFAULT( SHARED )  &
!    !$OMP PRIVATE( I, J, T )
!    DO T = 1, TIMES_A1
!
!       ! Intitialize the LWI array w/ the time-invariant part
!       lwiIn(:,:,T) = lwiMask
!
!       ! Also factor in the time-varying sea ice
!       DO J = 1, J05x0666
!       DO I = 1, I05x0666
!          IF ( frSeaIce(I,J,T) > 0.5e0 .and. lwiIn(I,J,T) < 1e0 ) THEN
!             lwiIn(I,J,T) = 2e0
!          ENDIF
!       ENDDO
!       ENDDO
!
!       ! Regrid the LWI field to coarse resolution
!       CALL Geos57RegridLwi( lwiIn(:,:,T), lwiOut(:,:,T), map, IMX, JMX )
!
!    ENDDO
!    !$OMP END PARALLEL DO
!    
!  END SUBROUTINE Geos57CreateLwi
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57RegridLwi
!!
!! !DESCRIPTION: This routine regrids the land-water indices (LWI) field.  
!!  Instead of an actual area regridding, we pick the mode of the LWI values
!!  of each 0.5 x 0.667 box that fits into a coarse grid box.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57RegridLwi( lwiIn, lwiOut, map, IMX, JMX )
!!
!! !INPUT PARAMETERS:
!!
!    ! Land-water indices on the 0.5. x 0.666 grid
!    REAL*4,       INTENT(IN)  :: lwiIn(I05x0666,J05x0666)
!
!    ! Object that contains the mapping weights to the output grid
!    TYPE(MapObj), POINTER     :: map(:,:)
!
!    ! Dimensions of the output grid
!    INTEGER,      INTENT(IN)  :: IMX, JMX
!!
!! !OUTPUT PARAMETERS:
!!
!    ! Regridded land-water indices on the output grid
!    REAL*4,       INTENT(OUT) :: lwiOut(IMX,JMX)
!!
!! !REVISION HISTORY: 
!!  25 Aug 2010 - R. Yantosca - Initial version, based on GEOS-5
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!    ! Local variables
!    INTEGER :: I, J, X, Y, Nx, Ny, nPoints, index
!    INTEGER :: hist(0:2)
!    REAL*4  :: mode(1)
!
!    ! Loop over grid boxes
!    DO J = 1, JMX
!    DO I = 1, IMX
!
!       ! Number of "fine" grid boxes in each dimension
!       ! that comprise a "coarse" grid box
!       nPoints = map(I,J)%nPoints
!
!       ! Zero the histogram array
!       hist = 0
!       
!       ! Loop over "fine" grid boxes
!       DO Ny = 1, nPoints
!       DO Nx = 1, nPoints
!          
!          ! Avoid useless clock cycles if the mapping weight is zero
!          IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN
!
!             ! Indices of each "fine" grid box that makes up the "coarse" box
!             X             = map(I,J)%xInd(Nx)
!             Y             = map(I,J)%yInd(Ny)
!
!             ! Sort each LWI value on the "fine" grid into a histogram
!             ! Possible values of LWI are 0, 1, 2
!             index       = INT( lwiIn(X,Y) ) 
!             hist(index) = hist(index) + 1 
!          ENDIF
!       ENDDO
!       ENDDO
!
!       ! The bin in the histogram w/ the most counts is the MODE,
!       ! so we'll use that as the regridded value of LWI.
!       !
!       ! NOTE: The result from MAXLOC is indexed starting from 1 and not
!       ! from zero...so we have to subtract 1 and then save to LWIOUT.
!       !
!       ! ALSO NOTE: If two or more elements of HIST have the same value,
!       ! then MAXLOC will pick the one that comes first in array order. 
!       ! So if a coarse box is 50% water and 50% ice, this regridding scheme 
!       ! will label the box as water. (LWIOUT(I,J)=0).  For 50% land and 50% 
!       ! ice, the coarse box will be labeled as land (LWIOUT(I,J)=1). 
!       mode        = MAXLOC( hist )
!       lwiOut(I,J) = mode(1) - 1
!    ENDDO
!    ENDDO
!
!  END SUBROUTINE Geos57RegridLwi
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57AdjustSnomas
!!
!! !DESCRIPTION: Routine Geos57AdjustSnomas will adjust the MERRA SNOMAS field
!!  to make it more similar to the GEOS-5 SNOMAS field.  This is necessary for
!!  backward compatibility with existing GEOS-Chem routines.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57AdjustSnomas( Q )
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    REAL*4, INTENT(INOUT) :: Q(I05x0666,J05x0666)  ! MERRA SNOMAS data [kg/m2]
!!
!! !REMARKS:
!!  The SNOMAS field in MERRA differs from that in the GEOS-5 ops data:
!!                                                                             .
!!  From the GEOS-5 File Specification Document:
!!                                                                             .
!!     SNOMAS: The mass of snow in per unit of land area in meters of 
!!     liquid-water-equivalent depth (i.e., 10^3 kg/m2). In grid boxes 
!!     with no land (FRLAND+FRLANDICE=0) it is set to _FillValue (= 1e15). 
!!     Where FRLANDICE>0.9 it is arbitrarily set to 4 meters. Over other 
!!     land areas it represents an average over the non-glaciated part.
!!                                                                             .
!!  From the MERRA File Specification Document:
!!                                                                             .
!!     SNOMAS: The mass of snow per unit of ice-free land area (FRLAND), 
!!     in kg/m2. In grid boxes with no land it is set to _FillValue (=1e15). 
!!     Over other land areas it represents an average over the nonglaciated
!!     part.
!!                                                                             .
!!  Max Suarez (Max.J.Suarez@nasa.gov) clarifies this difference:
!!                                                                             .
!!     Early versions of GEOS had been writing SNOMAS in meters.  This was 
!!     changed to mm (or kg/m^2) in  all recent versions, including 5_2, 
!!     MERRA, and the current development tags.  But the forward processing 
!!     spec was not updated. The MERRA spec, however, is correct.
!!                                                                             .
!!     To further complicate matters, the variable called SNOMAS in MERRA 
!!     comes from a very different part of the code that in 5_x.  It is 
!!     in a land collection intended to have representative values over 
!!     the ice-free land portion of the grid box.  This applies to all 
!!     variables in that collection.  SNOMAS in particular makes this 
!!     clear in the glossary definition in the spec.
!!                                                                             .
!!     Forward processing (FP) puts out a grid averaged SNOMAS, including
!!     ice-covered areas, where the "SNOMAS" was arbitrarily set to 4000 mm. 
!!     Neither FP nor MERRA includes ocean or freshwater regions with snow 
!!     over ice.
!!                                                                             .
!!     To answer your question, to convert MERRA SNOMAS to 5.2 SNOMAS,
!!     use this equation:
!!                                                                             .
!!     SNOMAS_5.X = ( SNOMAS_merra * FRLAND_merra + 
!!                    4000         * FRLANDICE_merra ) /
!!                  ( FRLAND_merra + FRLANDICE_merra )
!!                                                                             .
!!     Sorry about the confusion, but it seemed silly in MERRA to continue
!!     writing an invented value over glaciers. The next FP system will be
!!     like MERRA in this regard.  In the longer term, we plan to have a 
!!     better snow/ice parameterization over permanent glaciers.   
!!                                                                             .
!!  Therefore, we shall implement the algorithm that Max Suarez described
!!  above.  This will make the output SNOMAS field similar to GEOS-5,
!!  which will allow better backward compatibility w/ existing code.
!!                                                                             .
!!  Also note: Liquid water equivalent height is defined as such:
!!     1 m H2O = 10^3 kg/m2   ==>  10^-3 m H2O = 1 mm H2O = 1 kg/m2. 
!!
!! !REVISION HISTORY: 
!!  25 Aug 2010 - R. Yantosca - Initial version, based on GEOS-5
!!  27 Aug 2010 - R. Yantosca - Updated to use Max Suarez algorithm
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES
!!    
!    INTEGER :: I, J
!    REAL*4  :: den
!
!    ! Loop over surface grid boxes
!    DO J = 1, J05x0666
!    DO I = 1, I05x0666
!
!       ! The denominator is the sum of the land ice and land fractions in 
!       ! he grid box.  Nonzero denotes that we are over land and not ocean
!       den = frLand(I,J) + frLandIce(I,J)
!
!       ! Test if the division is possilble
!       IF ( den > 0e0 ) THEN
!
!          ! If so, then compute the SNOMAS value according to the
!          ! algorithm described above
!          Q(I,J) = (  Q(I,J) * frLand(I,J) + 4.0e3 * frLandIce(I,J) ) / den
!
!       ELSE
!
!          ! Otherwise, then we are over the ocean, so SNOMAS = 0
!          Q(I,J) = 0e0
!          
!       ENDIF
!
!   ENDDO
!   ENDDO
!
!  END SUBROUTINE Geos57AdjustSnomas
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57ProcessAlbedo
!!
!! !DESCRIPTION: Subroutine Geos57ProcessAlbedo computes the daily average
!!  albedo from the MERRA raw data.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57ProcessAlbedo( Q )
!!
!! !INPUT/OUTPUT PARAMETERS: 
!!
!    REAL*4, INTENT(INOUT) :: Q(:,:,:)  ! TROPP [hPa]
!!
!! !REVISION HISTORY: 
!!  23 Jul 2010 - R. Yantosca - Initial version, based on Geos5RegridModule.f90
!!
!! !REMARKS:
!!   Rationale for doing this: 
!!   ----------------------------------------------------------------------
!!   The MERRA ALBEDO field is only defined where it is daylight.
!!   Some places in GEOS-Chem require an ALBEDO field even at night (i.e.
!!   as a proxy for determining land surface).  Therefore compute the
!!   daily average albedo and return to the data processing routine above.
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER :: I, J, T
!
!    ! Arrays
!    ! For albedo processing
!    REAL*4  :: A ( SIZE( Q, 1 ), SIZE( Q, 2 ) )
!    INTEGER :: Ct( SIZE( Q, 1 ), SIZE( Q, 2 ) )
!
!    ! Initialization
!    A  = 0e0
!    Ct = 0
!
!    ! Sum up albedo over the entire day
!    DO T = 1, SIZE( Q, 3 )
!    DO J = 1, SIZE( Q, 2 )
!    DO I = 1, SIZE( Q, 1 )
!       IF ( Q(I,J,T) > 0e0 ) THEN
!          A (I,J) = A (I,J) + Q(I,J,T)
!          Ct(I,J) = Ct(I,J) + 1
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
!
!    ! Compute average albedo
!    ! Assign 0.85 for snow/ice over poles 
!    DO T = 1, SIZE( Q, 3 )
!    DO J = 1, SIZE( Q, 2 )
!    DO I = 1, SIZE( Q, 1 )
!       IF ( ct(I,J) > 0 ) THEN
!          Q(I,J,T) = A(I,J) / REAL( ct(I,J) )
!       ELSE
!          Q(I,J,T) = 0.85e0
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
!
!  END SUBROUTINE Geos57ProcessAlbedo
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Geos57ProcessTropp
!!
!! !DESCRIPTION: Subroutine "GeosProcessTropp" replaces any missing values in 
!!  the TROPP field with a zonal mean average.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Geos57ProcessTropp( Q )
!!
!! !INPUT/OUTPUT PARAMETERS: 
!!
!    REAL*4, INTENT(INOUT) :: Q(:,:,:)  ! TROPP [hPa]
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A3Module.F90
!!
!! !REMARKS:
!!   Rationale for doing this: 
!!   ----------------------------------------------------------------------
!!   Sometimes the TROPP field has missing values, so we need to replace 
!!   those with the average of the other boxes in the same latitude.  
!!   Otherwise those missing values will get reset to zeroes (in routine
!!   Geos5MakeA1Files), and those zeroes will propagate through the 
!!   regridding process.  
!!
!!   If a box has a zero TROPP value, then when it is regridded to a 
!!   coarser resolution, the resultant TROPP values will not be realistic.  
!!   This will cause GEOS-Chem to diagnose the tropopause as much higher 
!!   than it should.  Resetting the TROPP values according to the algorithm
!!   below will avoid this problem.
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!    ! Local variables
!    INTEGER                 :: IX, JX, TX
!    INTEGER                 :: I,  I2, J,  J2, T
!    REAL*4                  :: tot, ct
!    CHARACTER(LEN=MAX_CHAR) :: msg
!
!    !=======================================================================
!    ! Geos5ProcessTropp begins here!
!    !=======================================================================
!
!    ! Echo info
!    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57ProcessTropp %%%%%%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Test if missing values are found
!    IF ( ANY( Q == FILL_VALUE ) ) THEN 
!
!       ! If yes, then echo a message
!       msg = '%%% Missing data values found in TROPP!  Removing these ...'
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ELSE
!
!       ! If no, echo a message, then exit
!       msg = '%%% No missing data values found in TROPP!  Continuing on ...'
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       msg = '%%%%%%%%%% LEAVING ROUTINE Geos57ProcessTropp %%%%%%%%%%'
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       WRITE( IU_LOG, '(a)' ) '%%%'
!       RETURN
!
!    ENDIF
!
!    ! Dimensions of the array
!    IX = SIZE( Q, 1 )
!    JX = SIZE( Q, 2 )
!    TX = SIZE( Q, 3 )
!
!    ! Loop over grid boxes
!    DO T = 1, TX
!    DO J = 1, JX
!    DO I = 1, IX
!
!       ! Replace "missing" values with a zonal average pressure
!       IF ( Q(I,J,T) == FILL_VALUE ) THEN
!             
!          ! Zero summing variables
!          tot = 0e0
!          ct  = 0e0
!             
!          ! Sum up "good" boxes at this latitude
!          DO I2 = 1, IX
!             IF ( Q(I2,J,T) < FILL_VALUE ) THEN
!                tot = tot + Q(I2,J,T)
!                ct  = ct  + 1e0
!             ENDIF
!          ENDDO
!          
!          ! Avoid div by zero
!          IF ( ct > 0e0 ) THEN 
!
!             ! Replace "bad" value with zonal mean of "good" values
!             Q(I,J,T) = tot / ct 
!             
!          ELSE
!                
!             ! If ct==0 then we have no good data at this latitude
!             IF ( J > JX/2 ) THEN
!                
!                ! Northern hemisphere
!                ! Then search down until the next good value
!                DO J2 = J, 1, -1 
!                   IF ( Q(I,J2,T) < FILL_VALUE ) THEN
!                      Q(I,J,T) = Q(I,J2,T)
!                   ENDIF
!                ENDDO
!                
!             ELSE
!
!                ! Southern hemisphere
!                ! Then search up until the next good value
!                DO J2 = 1, J
!                   IF ( Q(I,J2,T) < FILL_VALUE ) THEN
!                      Q(I,J,T) = Q(I,J2,T)
!                   ENDIF
!                ENDDO
!                
!             ENDIF
!          ENDIF
!       ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
!
!    ! Echo info
!    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57ProcessTropp %%%%%%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    WRITE( IU_LOG, '(a)' ) '%%%'
!
!
!  END SUBROUTINE Geos57ProcessTropp
!EOC
END MODULE Geos57A1Module

