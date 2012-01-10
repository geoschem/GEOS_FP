!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57A3Module
!
! !DESCRIPTION: Module Geos57A3Module contains routines to create the 
!  GEOS-Chem average 3-hr data files from the Geos57 raw data.
!\\
!\\
! !INTERFACE: 

MODULE Geos57A3Module
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

# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Geos57MakeA3
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dCldNv
  PRIVATE :: Process3dUdtNv
  PRIVATE :: Process3dOptDep
  !PRIVATE :: Process3dMstNv
  !PRIVATE :: Process3dMstNe
  !PRIVATE :: Process3dQdtNv
  !PRIVATE :: Process3dRadNv
  PRIVATE :: RegridTau
!
! !REVISION HISTORY:
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
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
  SUBROUTINE NcOutFileDef( X,        Y,           Z,    T,    &
                           xMid,     yMid,        zMid, time, &
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
    REAL*4,           INTENT(IN)    :: zMid(Z)       ! Array of alt centers
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,    DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr,  msg
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
    lName = 'GEOS-5.7.2 time-averaged 3-hour (A3) fields for GEOS-Chem'
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
    lName = '00:00:00.000000'
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',  TRIM( lName )   )

    ! End Date
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',    yyyymmdd_string )

    ! End Time
    lName = '23:59:59.99999'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',    TRIM( lName )   )

    ! Delta-time
    lName = '30000'
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
    CALL NcDef_Glob_Attributes( fOut, 'Delta_lon',   TRIM( DI    )   )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_lat',   TRIM( DJ    )   )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'lon',  X,   idLon  )
    CALL NcDef_Dimension( fOut, 'lat',  Y,   idLat  )
    CALL NcDef_Dimension( fOut, 'lev',  Z,   idLev  )
    CALL NcDef_Dimension( fOut, 'time', T,   idTime )
    CALL NcDef_Dimension( fOut, 'ap',   Z+1, idAp   )
    CALL NcDef_Dimension( fOut, 'bp',   Z+1, idBp   )

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
  
    ! Level index array
    var1    = (/ idLev /)
    vId     = vId + 1
    lName   = 'levels'
    units   = 'unitless'
    CALL NcDef_Variable      ( fOut, 'lev', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    ! Time index array
    var1    = (/ idTime /)
    vId     = vId + 1
    lName   = 'time'
    units   = UnitsForTime( yyyymmdd )
    delta_t = '0000-00-00 03:00:00'
    begin_d = yyyymmdd_string
    begin_t = '000000'
    incr    = '030000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Hybrid grid Ap index array
    var1    = (/ idAp /)
    vId     = vId + 1
    lName   = 'Hybrid grid A parameter [ P = ap + ( bp * Psurface ) ]'
    units   = 'hPa'
    CALL NcDef_Variable      ( fOut, 'ap', NF_FLOAT, 1, var1, vId            )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 
 
    ! Hybrid grid Bp index array
    var1    = (/ idBp /)
    vId     = vId + 1
    lName   = 'Hybrid grid B parameter [ P = ap + ( bp * Psurface ) ]'
    units   = 'unitless'
    CALL NcDef_Variable      ( fOut, 'bp', NF_FLOAT, 1, var1, vId            )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! CLOUD
    IF ( StrPos( 'CLOUD', tavg3_3d_rad_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Total cloud fraction in grid box'
       units = 'unitless'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'CLOUD', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! DTRAIN
    IF ( StrPos( 'DTRAIN', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Detrainment cloud mass flux'
       units = 'kg m-2 s-2'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'DTRAIN', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! OMEGA
    IF ( StrPos( 'OMEGA', tavg3_3d_udt_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Vertical pressure velocity'
       units = 'Pa s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'OMEGA', NF_FLOAT, 4, var4, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! OPTDEPTH
    IF ( StrPos( 'OPTDEPTH', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Total in-cloud optical thickness (visible band)'
       units = 'unitless'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'OPTDEPTH', NF_FLOAT, 4, var4, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! QI
    IF ( StrPos( 'QI', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Cloud ice water mixing ratio'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QI', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! QL
    IF ( StrPos( 'QL', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Cloud liquid water mixing ratio'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QL', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! RH
    IF ( StrPos( 'RH', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Relative humidity'
       units = 'fraction'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'RH', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! TAUCLI
    IF ( StrPos( 'TAUCLI', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'In-cloud ice optical thickness (visible band)'
       units = 'unitless'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'TAUCLI', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! TAUCLW
    IF ( StrPos( 'TAUCLW', tavg3_3d_cld_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'In-cloud water optical thickness (visible band)'
       units = 'unitless'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'TAUCLW', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! U
    IF ( StrPos( 'U', tavg3_3d_udt_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Eastward component of wind'
       units = 'm s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'U', NF_FLOAT, 4, var4, vId          )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! V
    IF ( StrPos( 'V', tavg3_3d_udt_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Northward component of wind'
       units = 'm s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'V', NF_FLOAT, 4, var4, vId          )
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
    CALL NcWr( xMid, fOut, 'lon',  (/ 1 /), (/ X   /) )
    CALL NcWr( yMid, fOut, 'lat',  (/ 1 /), (/ Y   /) )
    CALL NcWr( zMid, fOut, 'lev',  (/ 1 /), (/ Z   /) )
    CALL NcWr( time, fOut, 'time', (/ 1 /), (/ T   /) )
    CALL NcWr( Ap,   fOut, 'ap',   (/ 1 /), (/ Z+1 /) )
    CALL NcWr( Bp,   fOut, 'bp',   (/ 1 /), (/ Z+1 /) )

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
! !IROUTINE: Geos57MakeA3
!
! !DESCRIPTION: Routine Geos57MakeA3 is the the driver routine for 
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (surface values) from 
!       the Geos57 raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data in a format that GEOS-Chem can read.
! \end{enumerate}
! This routine is called directly from the main program Geos57Driver.F90
!\\
!\\
! !INTERFACE:
  SUBROUTINE Geos57MakeA3
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A6Module.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_3dCldNv
    INTEGER                 :: nFields_3dMstNv
    INTEGER                 :: nFields_3dMstNe
    INTEGER                 :: nFields_3dQdtNv
    INTEGER                 :: nFields_3dRadNv
    INTEGER                 :: nFields_3dUdtNv
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dCldNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dMstNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dMstNe(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dQdtNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dRadNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dUdtNv(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57MakeA3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg3_3d_cld_Nv_data ) // ',' // &
                    TRIM( tavg3_3d_mst_Nv_data ) // ',' // &
                    TRIM( tavg3_3d_mst_Ne_data ) // ',' // &
                    TRIM( tavg3_3d_qdt_Nv_data ) // ',' // &
                    TRIM( tavg3_3d_rad_Nv_data ) // ',' // &
                    TRIM( tavg3_3d_udt_Nv_data )

    ! Return the list of fields and number of fields to process
    ! from each of the Geos57 raw met data files
    CALL GetNFields( tavg3_3d_cld_Nv_data, nFields_3dCldNv, fields_3dCldNv )
    CALL GetNFields( tavg3_3d_mst_Nv_data, nFields_3dMstNv, fields_3dMstNv )
    CALL GetNFields( tavg3_3d_mst_Ne_data, nFields_3dMstNe, fields_3dMstNe )
    CALL GetNFields( tavg3_3d_qdt_Nv_data, nFields_3dQdtNv, fields_3dQdtNv )
    CALL GetNFields( tavg3_3d_rad_Nv_data, nFields_3dRadNv, fields_3dRadNv )
    CALL GetNFields( tavg3_3d_udt_Nv_data, nFields_3dUdtNv, fields_3dUdtNv )
    CALL GetNFields( allFieldsList,        nAllFields,      allFields      )
    
    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_cld_Nv_file ), nFields_3dCldNv
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_mst_Nv_file ), nFields_3dMstNv
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_mst_Ne_file ), nFields_3dMstNe
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_qdt_Nv_file ), nFields_3dQdtNv
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_rad_Nv_file ), nFields_3dRadNv
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_udt_Nv_file ), nFields_3dUdtNv
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )


    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested China output file
    IF ( doNestCh ) THEN
       fName = dataTmplNestCh
       gName = 'SEA4CRS'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                 )      
       CALL StrRepl     ( fName,     '%%',     'a3'                    )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  L025x03125, TIMES_A3,   &
                          xMid_025x03125(I0_ch:I1_ch),                  &
                          yMid_025x03125(J0_ch:J1_ch),                  &
                          zMid_025x03125,                   a3MinsI,    &
                          gName,    fName,      fOutNestCh             )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = dataTmpl2x25
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                 )      
       CALL StrRepl     ( fName,     '%%',      'a3'                   )
       CALL NcOutFileDef( I2x25,     J2x25,     L2x25,      TIMES_A3,   &
                          xMid_2x25, yMid_2x25, zMid_2x25,  a3MinsI,    &
                          gName,     fName,     fOut2x25               )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = dataTmpl4x5
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                 )      
       CALL StrRepl     ( fName,     '%%',     'a3'                    )
       CALL NcOutFileDef( I4x5,      J4x5,      L4x5,       TIMES_A3,   &
                          xMid_4x5,  yMid_4x5,  zMid_4x5,   a3MinsI,    &
                          gName,     fName,     fOut4x5                )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process3dCldNv ( nFields_3dCldNv, fields_3dCldNv ) ! tavg3_3d_cld_Nv
!    CALL Process3dMstNv( nFields_3dMstNv, fields_3dMstNv ) ! tavg3_3d_mst_Nv
!    CALL Process3dMstNe( nFields_3dMstNe, fields_3dMstNe ) ! tavg3_3d_mst_Ne
!    CALL Process3dQdtNv( nFields_3dQdtNv, fields_3dQdtNv ) ! tavg3_3d_qdt_Nv
!    CALL Process3dRadNv( nFields_3dRadNv, fields_3dRadNv ) ! tavg3_3d_rad_Nv
    CALL Process3dUdtNv ( nFields_3dUdtNv, fields_3dUdtNv ) ! tavg3_3d_udt_Nv
    CALL Process3dOptDep(                                 ) ! optical depths
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A3 output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57MakeA3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'   
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Geos57MakeA3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dCldNv
!
! !DESCRIPTION: Subroutine Process3dCldNv regrids the Geos57 met fields from 
!  the "tavg3\_3d\_cld\_Nv" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dCldNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REMARKS:
!  The cloud fraction field CLOUD and cloud optical depth fields TAUCLI, 
!  TAUCLW, and OPTDEPTH are processed separately in routine Process3dOptDep.  
!  This is because these fields must all be regridded together using the
!  algorithm developed by Hongyu Liu (in routine RegridTau).!
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version
!  10 Jan 2012 - R. Yantosca - Activate parallel loop over vertical levels
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
    INTEGER                 :: X2x25,    Y2x25,    Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,     Z4x5,    T4x5
    INTEGER                 :: st4d(4),  ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Q      ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: QI     ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: QL     ( I025x03125, J025x03125, L025x03125 )
    REAL*4                  :: Q2x25  ( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: QI_2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: QL_2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: Q4x5   ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: QI_4x5 ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: QL_4x5 ( I4x5,       J4x5,       L4x5       )

    ! Pointer arrays
    REAL*4, POINTER         :: Ptr(:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Initialization
    !=======================================================================

     ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dChmFv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'lev',  ZNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
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

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR DEBUGGING -- Sample data is only up to hour 19:30:00
       if ( hhmmss > 193000 ) CYCLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_cld_Nv_file )
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

       ! Zero QI, QL arrays
       QI      = 0e0
       QI_2x25 = 0e0
       QI_4x5  = 0e0
       QL      = 0e0
       QL_2x25 = 0e0
       QL_4x5  = 0e0

       !====================================================================
       ! Process data
       !====================================================================
       
       ! Loop over data fields
       DO F = 1, nFields
          
          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )
          
          ! Skip certain fieldnames
          SELECT CASE ( name ) 
             CASE( '' )                               ! Null string
                CYCLE
             CASE( 'QI', 'QL' )                       ! These fields are
                CYCLE                                 !  derived, not read
             CASE( 'TAUCLI', 'TAUCLW', 'OPTDEPTH' )   ! These fields are
                CYCLE                                 !  procesed elsewhere
             CASE DEFAULT
                ! Nothing
          END SELECT

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0

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

          !-----------------------------------------------------------------
          ! Process data (or save for later special handling)
          !-----------------------------------------------------------------
          SELECT CASE( name )

             CASE ( 'QIAN', 'QILS' )
                QI   = QI  + Q        ! QI = QIAN + QILS; regrid below

             CASE ( 'QLAN', 'QLLS' )
                QL   = QL  + Q        ! QL = QLAN + QLLS; regrid below

             CASE DEFAULT 

                !-----------------------------------------------------------
                ! Regrid data fields (all except QI, QL)
                !-----------------------------------------------------------
                msg = '%%% Regridding ' // name
                WRITE( IU_LOG, '(a)' ) TRIM( msg )

                ! Loop over the A-3 times and vertical levels
!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) & 
!$OMP PRIVATE( L, LR )
                DO L = 1, L025x03125

                   ! Reverse level index for output arrays
                   LR = L025x03125 - L + 1

                   ! Regrid to 2 x 2.5
                   IF ( do2x25 ) THEN
                      CALL RegridGeos57To2x25( 0, Q(:,:,LR), Q2x25(:,:,L) )
                   ENDIF

                    ! Regrid to 4x5 
                   IF ( do4x5 ) THEN
                      CALL RegridGeos57To4x5 ( 0, Q(:,:,LR), Q4x5(:,:,L)  )
                   ENDIF

                ENDDO
!$OMP END PARALLEL DO
                
                !-----------------------------------------------------------
                ! Write netCDF output (all except QI, QL)
                !-----------------------------------------------------------
                msg = '%%% Archiving  ' // name
                WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
                ! Nested China (point to proper slice of global data)
                IF ( doNestCh ) THEN
                   Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch, : )
                   st4d = (/ 1,       1,       1,          H /)
                   ct4d = (/ XNestCh, YNestCh, L025x03125, 1 /)
                   CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st4d, ct4d )
                ENDIF
                
                ! Write 2 x 2.5 data
                IF ( do2x25 ) THEN
                   st4d = (/ 1,     1,     1,     H  /)
                   ct4d = (/ X2x25, Y2x25, L2x25, 1  /)
                   CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st4d, ct4d )
                ENDIF
       
                ! Write 4x5 data
                IF ( do4x5 ) THEN
                   st4d  = (/ 1,    1,    1,    H /)
                   ct4d  = (/ X4x5, Y4x5, L4x5, 1 /)
                   CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st4d, ct4d )
                ENDIF
                
          END SELECT

       ENDDO

       !====================================================================
       ! Regrid QI and QL
       !====================================================================   
       msg = '%%% Regridding QI and QL'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       ! Loop over the A-3 times and vertical levels
!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) & 
!$OMP PRIVATE( L, LR )
       DO L = 1, L025x03125
          
          ! Reverse level index for output arrays
          LR = L025x03125 - L + 1

          ! Regrid to 2 x 2.5
          IF ( do2x25 ) THEN
             CALL RegridGeos57To2x25( 0, QI(:,:,LR), QI_2x25(:,:,L) )
             CALL RegridGeos57To2x25( 0, QL(:,:,LR), QL_2x25(:,:,L) )
          ENDIF
          
          ! Regrid to 4x5 
          IF ( do4x5 ) THEN
             CALL RegridGeos57To4x5 ( 0, QI(:,:,LR), QI_4x5(:,:,L)  )
             CALL RegridGeos57To4x5 ( 0, QL(:,:,LR), QL_4x5(:,:,L)  )
          ENDIF

       ENDDO
!$OMP END PARALLEL DO

       !====================================================================
       ! Write QI and QL to netCDF
       !==================================================================== 

       msg = '%%% Archiving  QI and QL'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )

       !-----------------------------
       ! SEA4CRS NESTED CHINA GRID
       !-----------------------------
       IF ( doNestCh ) THEN

          ! netCDF indices
          st4d = (/ 1,       1,       1,          H /)
          ct4d = (/ XNestCh, YNestCh, L025x03125, 1 /)

          ! QI
          Ptr  => QI( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'QI', st4d, ct4d )

          ! QL
          Ptr  => QL( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'QL', st4d, ct4d )

       ENDIF
       
       !-----------------------------
       ! 2 x 2.5 GLOBAL GRID
       !-----------------------------
       IF ( do2x25 ) THEN

          ! netCDF indices
          st4d = (/ 1,     1,     1,     H  /)
          ct4d = (/ X2x25, Y2x25, L2x25, 1  /)

          ! Write data
          CALL NcWr( QI_2x25, fOut2x25, 'QI', st4d, ct4d )
          CALL NcWr( QL_2x25, fOut2x25, 'QL', st4d, ct4d )

       ENDIF

       !-----------------------------
       ! 4 x 5 GLOBAL GRID
       !-----------------------------
       IF ( do4x5 ) THEN

          ! netCDF indices
          st4d  = (/ 1,    1,    1,    H /)
          ct4d  = (/ X4x5, Y4x5, L4x5, 1 /)

          ! Write data
          CALL NcWr( QI_4x5, fOut4x5, 'QI', st4d, ct4d )
          CALL NcWr( QL_4x5, fOut4x5, 'QL', st4d, ct4d )

       ENDIF

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
    msg = '%%%%%% LEAVING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dCldNv
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dUdtNv
!
! !DESCRIPTION: Subroutine Process3dUdtNv regrids the Geos57 met fields from 
!  the "tavg3\_3d\_udt\_Nv" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dUdtNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version
!  10 Jan 2012 - R. Yantosca - Activate parallel loop over vertical levels
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: F,        H,       L,       LR
    INTEGER                 :: hhmmss 

    ! Variables for netCDF I/O
    INTEGER                 :: X,        Y,       Z,       T
    INTEGER                 :: XNestCh,  YNestCh, ZNestCh, TNestCh
    INTEGER                 :: X2x25,    Y2x25,   Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,    Z4x5,    T4x5
    INTEGER                 :: st3d(3),  ct3d(3), st4d(4), ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: P    ( I025x03125, J025x03125             )
    REAL*4                  :: Q2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: P2x25( I2x25,      J2x25                  )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: P4x5 ( I4x5,       J4x5                   )

    ! Pointer arrays
    REAL*4, POINTER         :: Ptr(:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

     ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dUdtNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'lev',  ZNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
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

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR DEBUGGING -- Sample data is only up to hour 19:30:00
       if ( hhmmss > 193000 ) CYCLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_udt_Nv_file )
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
       ! Process surface pressure data
       !====================================================================

       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st3d = (/ 1, 1, 1 /)
       ct3d = (/ X, Y, 1 /)
       
       ! Read surface pressure
       CALL NcRd( P, fIn, 'PS', st3d, ct3d )

       ! Regrid surface pressure
       IF ( do2x25 ) CALL RegridGeos57To2x25( 0, P, P2x25 )
       IF ( do4x5  ) CALL RegridGeos57To4x5 ( 0, P, P4x5  )

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

          !-----------------------------------------------------------------
          ! Regrid data fields 
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Loop over the A-3 times and vertical levels
!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE ( L, LR )
          DO L = 1, L025x03125
             
             ! Reverse level index for output arrays
             LR = L025x03125 - L + 1

             !%%% Pre-regrid special handling %%%
             !%%% Multiply U, V by pressure   %%%
             IF ( name == 'U' .or. name == 'V' ) THEN
                Q(:,:,LR) = Q(:,:,LR) * P
             ENDIF

             !%%% Regrid to 2 x 2.5 %%%
             IF ( do2x25 ) THEN
                CALL RegridGeos57To2x25( 0, Q(:,:,LR), Q2x25(:,:,L) )

                ! Post-regrid handling, divide winds by pressure
                IF ( name == 'U' .or. name == 'V' ) THEN
                   Q2x25(:,:,L) = Q2x25(:,:,L) / P2x25
                ENDIF
             ENDIF

             ! %%% Regrid to 4 x 5 %%%
             IF ( do4x5 ) THEN
                CALL RegridGeos57To4x5 ( 0, Q(:,:,LR), Q4x5(:,:,L)  )

                ! Post-regrid handling, divide winds by pressure
                IF ( name == 'U' .or. name == 'V' ) THEN
                   Q4x5(:,:,L) = Q4x5(:,:,L) / P4x5
                ENDIF
             ENDIF             
          ENDDO
!$OMP END PARALLEL DO
                
          !-----------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          IF ( doNestCh ) THEN
             Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch, : )
             st4d = (/ 1,       1,       1,          H /)
             ct4d = (/ XNestCh, YNestCh, L025x03125, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name ), st4d, ct4d )
          ENDIF
                
          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st4d = (/ 1,     1,     1,     H  /)
             ct4d = (/ X2x25, Y2x25, L2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st4d, ct4d )
          ENDIF
       
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st4d  = (/ 1,    1,    1,    H /)
             ct4d  = (/ X4x5, Y4x5, L4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st4d, ct4d )
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
    msg = '%%%%%% LEAVING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dUdtNv
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dOptDep
!
! !DESCRIPTION: Subroutine Process3dOptDep regrids the CLOUD, TAUCLI,
!  TAUCLW, and OPTDEPTH fields using Hongyu Liu's regridding algorithm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dOptDep()
!
! !REMARKS:
!  The OPTDEPTH and CLOUD fields are regridded following the algorithm of
!  Hongyu Liu (cf. routine RegridTau).  OPTDEPTH = TAUCLI + TAUCLW.
!
! !REVISION HISTORY: 
!  09 Jan 2012 - R. Yantosca - Initial version
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
    INTEGER                 :: X2x25,    Y2x25,    Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,     Z4x5,    T4x5
    INTEGER                 :: st3d(3),  st4d(4)
    INTEGER                 :: ct3d(3),  ct4d(4)

    ! Data arrays
    REAL*4, TARGET          :: Cld      ( I025x03125, J025x03125, L025x03125 )
    REAL*4, TARGET          :: OptD     ( I025x03125, J025x03125, L025x03125 )
    REAL*4, TARGET          :: TauI     ( I025x03125, J025x03125, L025x03125 )
    REAL*4, TARGET          :: TauW     ( I025x03125, J025x03125, L025x03125 )
    REAL*4                  :: Cld_2x25 ( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: OptD_2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: TauI_2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: TauW_2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: Cld_4x5  ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: OptD_4x5 ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: TauI_4x5 ( I4x5,       J4x5,       L4x5       )
    REAL*4                  :: TauW_4x5 ( I4x5,       J4x5,       L4x5       )

    ! Pointer arrays
    REAL*4, POINTER         :: ptr(:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Initialization
    !=======================================================================

    Cld   = 0e0
    OptD  = 0e0
    TauI  = 0e0
    Tauw  = 0e0

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dOptDep %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Nested China grid
    IF ( doNestCh ) THEN
       CALL NcGet_DimLen( fOutNestCh, 'lon',  XNestCh )
       CALL NcGet_DimLen( fOutNestCh, 'lat',  YNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'lev',  ZNestCh ) 
       CALL NcGet_DimLen( fOutNestCh, 'time', TNestCh )
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

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! Zero cloud arrays for each hour
       Cld       = 0e0
       Cld_2x25  = 0e0
       Cld_4x5   = 0e0
       TauI      = 0e0
       TauI_2x25 = 0e0
       TauI_4x5  = 0e0
       TauW      = 0e0
       TauW_2x25 = 0e0
       TauW_4x5  = 0e0
       OptD      = 0e0
       OptD_2x25 = 0e0
       OptD_4x5  = 0e0

       ! GMT time of day (hh:mm:ss)
       hhmmss    = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR DEBUGGING -- Sample data is only up to hour 19:30:00
       if ( hhmmss > 193000 ) CYCLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       !==================================================================
       ! Read CLOUD data from tavg3_3d_rad_Nv
       !==================================================================

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_rad_Nv_file )
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

       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st4d = (/ 1, 1, 1, 1 /)
       ct4d = (/ X, Y, Z, 1 /)

       ! CLOUD
       msg = '%%% Reading    CLOUD' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( Cld, fIn, 'CLOUD', st4d, ct4d )
       WHERE( Cld == FILL_VALUE ) Cld = 0e0    

       ! Close input netCDF file
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )

       !==================================================================
       ! Read TAUCLI, TAUCLW from tavg3_3d_cld_Nv
       !==================================================================

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_cld_Nv_file )
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

       ! Start and count index arrays for netCDF
       ! (There is only one data block per file)
       st4d = (/ 1, 1, 1, 1 /)
       ct4d = (/ X, Y, Z, 1 /)

       ! TAUCLI
       msg = '%%% Reading    TAUCLI' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( TauI, fIn, 'TAUCLI', st4d, ct4d )
       WHERE( TauI == FILL_VALUE ) TauI = 0e0    

       ! TAUCLW
       msg = '%%% Reading    TAUCLW' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( TauW, fIn, 'TAUCLW', st4d, ct4d )
       WHERE( TauW == FILL_VALUE ) TauI = 0e0    

       ! OPTDEPTH (construct from TAUCLI and TAUCLW)
       OptD = TauI + TauW

       ! Close input netCDF file
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )

       !====================================================================
       ! Regrid cloud & optical depth fields w/ Hongyu Liu's algorithm
       !====================================================================

       msg = '%%% Regridding CLOUD, TAUCLI, TAUCLW, OPTDEPTH'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
    
       ! Regrid to 2 x 2.5 (reverse levels of output array Q2x25)
       IF ( do2x25 ) THEN
          CALL RegridTau( OptD,      OptD_2x25, TauI,  TauI_2x25,  &
                          TauW,      TauW_2x25, Cld,   Cld_2x25,   &
                          mapTo2x25, I2x25,     J2x25, L2x25      )
       ENDIF
       
       ! Regrid to 4x5 (reverse levels of output array Q4x5)
       IF ( do4x5 ) THEN
          CALL RegridTau( OptD,      OptD_4x5,  TauI,  TauI_4x5,   &
                          TauW,      TauW_4x5,  Cld,   Cld_4x5,    &
                          mapTo4x5,  I4x5,      J4x5,  L4x5       )
       ENDIF
              
       !====================================================================
       ! Write to netCDF
       !====================================================================
 
       msg = '%%% Archiving  CLOUD, TAUCLI, TAUCLW, OPTDEPTH'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
       !%%%%% SEA4CRS NESTED CHINA GRID %%%%%
       IF ( doNestCh ) THEN

          ! netCDF index arrasy
          st4d = (/ 1,       1,       1,          H /)
          ct4d = (/ XNestCh, YNestCh, L025x03125, 1 /)

          ! CLOUD
          Ptr  => Cld( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'CLOUD',    st4d, ct4d )

          ! TAUCLI
          Ptr  => TauI( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW
          Ptr  => TauW( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD( I0_ch:I1_ch, J0_ch:J1_ch, : )
          CALL NcWr( Ptr, fOutNestCh, 'OPTDEPTH', st4d, ct4d )

       ENDIF

       !%%%%% 2 x 2.5 GLOBAL GRID %%%%%
       IF ( do2x25 ) THEN

          ! netCDF indices
          st4d = (/ 1,     1,     1,     H  /)
          ct4d = (/ X2x25, Y2x25, L2x25, 1  /)

          ! Write to disk
          CALL NcWr( Cld_2x25,  fOut2x25, 'CLOUD',    st4d, ct4d )
          CALL NcWr( TauI_2x25, fOut2x25, 'TAUCLI',   st4d, ct4d )
          CALL NcWr( TauW_2x25, fOut2x25, 'TAUCLW',   st4d, ct4d )
          CALL NcWr( OptD_2x25, fOut2x25, 'OPTDEPTH', st4d, ct4d )

       ENDIF

       !%%%%% 4 x 5 GLOBAL GRID %%%%%
       IF ( do4x5 ) THEN

          ! netCDF indices
          st4d = (/ 1,    1,    1,    H  /)
          ct4d = (/ X4x5, Y4x5, L4x5, 1  /)

          ! Write to disk
          CALL NcWr( Cld_4x5,  fOut4x5, 'CLOUD',    st4d, ct4d )
          CALL NcWr( TauI_4x5, fOut4x5, 'TAUCLI',   st4d, ct4d )
          CALL NcWr( TauW_4x5, fOut4x5, 'TAUCLW',   st4d, ct4d )
          CALL NcWr( OptD_4x5, fOut4x5, 'OPTDEPTH', st4d, ct4d )

       ENDIF
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================
    msg = '%%%%%% LEAVING ROUTINE Process3dOptDep %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dOptDep

!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process3dChmFe
!!
!! !DESCRIPTION: Subroutine Process3dChmFe regrids the Geos57 met fields from 
!!  the "tavg3\_3d\_chm\_Fe" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process3dChmFe( nFields, fields )
!!
!! !INPUT PARAMETERS:
!!
!    INTEGER,          INTENT(IN)    :: nFields     ! # of fields to process
!    CHARACTER(LEN=*), INTENT(IN)    :: fields(:)   ! List of field names
!!
!! !REVISION HISTORY: 
!!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A6Module.F90
!!  11 Aug 2010 - R. Yantosca - Now do not average to 6hr time resolution
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!
!    ! Scalars
!    INTEGER                 :: F, L, LR, T
!
!    ! Arrays
!    REAL*4, TARGET          :: Q( I1x125, J1x125, L1x125+1, TIMES_A3 )
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_1x125(:,:)
!    REAL*4, POINTER         :: ptr_2x25 (:,:)
!    REAL*4, POINTER         :: ptr_4x5  (:,:)
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
!    NULLIFY( ptr_1x125, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process3dChmFe %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg3_3d_chm_Fe_file )
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
!    CALL He4GridReadZ
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
!       !-----------------------------
!       ! Read data 
!       !-----------------------------
!
!       ! Zero data array
!       Q = 0e0
!
!       ! Read from file
!       msg = '%%% Reading    ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       CALL He4GridReadData( name, Q )
!
!       ! Remove fill values
!       WHERE( Q == FILL_VALUE ) Q = 0e0
!
!       !--------------------------------
!       ! Regrid to 2 x 2.5 & 4 x 5 
!       !--------------------------------
!       msg = '%%% Regridding ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over the A3 times and # of levels
!       !%%% For some reason the parallel loop doesn't work, disable it %%%
!       !!$OMP PARALLEL DO                                      &
!       !!$OMP DEFAULT( SHARED )                                &
!       !!$OMP PRIVATE( T, L, LR, ptr_1x125, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A3
!       DO L = 1, L1x125+1
!            
!          ! Reverse level index for output arrays
!          LR = ( L1x125+1 ) - L + 1
!
!          !-----------------------------
!          ! Do the regridding!
!          !-----------------------------
!          
!          ! Point to time-slice of native data array
!          ptr_1x125 => Q(:,:,L,T)
!          
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25e(:,:,LR,T)
!             CALL RegridGeos57FTo2x25( 0, ptr_1x125, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5e(:,:,LR,T)
!             CALL RegridGeos57FTo4x5( 0, ptr_1x125, ptr_4x5 )
!          ENDIF
!
!       ENDDO
!       ENDDO
!       !!$OMP END PARALLEL DO
!
!    ENDDO
!          
!    !=======================================================================
!    ! Cleanup & quit
!    !=======================================================================
!
!    ! Nullify all pointers
!    NULLIFY( ptr_1x125, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!  
!    msg = '%%%%%% LEAVING ROUTINE Process3dChmFe %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process3dChmFe
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process3dMstCp
!!
!! !DESCRIPTION: Subroutine Process3dMstCp regrids the Geos57 met fields from 
!!  the "tavg3\_3d\_mst\_Cp" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process3dMstCp( nFields, fields, offset )
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
!!  03 Aug 2010 - R. Yantosca - Initial version
!!  11 Aug 2010 - R. Yantosca - Now do not average to 6hr time resolution
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!
!    ! Scalars
!    INTEGER                 :: F, I, J, L, T
!
!    ! Arrays for 3hr data
!    REAL*4                  :: P_P(I125x125,J125x125,           TIMES_A3)
!    REAL*4                  :: PLE(I125x125,J125x125,L125x125+1,TIMES_A3)
!    REAL*4                  :: Q_P(I125x125,J125x125,L125x125_P,TIMES_A3)
!    REAL*4, TARGET          :: Q  (I125x125,J125x125,L125x125,  TIMES_A3)
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_125x125(:,:)
!    REAL*4, POINTER         :: ptr_2x25   (:,:)
!    REAL*4, POINTER         :: ptr_4x5    (:,:)
!
!     ! Character strings and arrays
!    CHARACTER(LEN=10      ) :: name
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
!    ! Zero pointers and data arrays
!    NULLIFY( ptr_125x125, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process3dMstCp %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    !=======================================================================
!    ! Read instantaneous surface pressure from the inst3_3d_asm_Cp
!    ! at 1.25 x 1.25 and compute the array of pressure at level edges 
!    !=======================================================================
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( inst3_3d_asm_Cp_file )
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
!    CALL He4GridReadZ
!    CALL He4GetNymdNhms
!
!    ! Read surface pressure data
!    P_P  = 0e0
!    msg  = '%%% Reading PS on Cp grid'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridReadData( 'PS', P_P )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!
!    ! Make PLE data [hPa] for the hybrid grid from the surface pressure
!    CALL MakePLE( P_P, PLE )
!    
!    !=======================================================================
!    ! Open the tavg3_3d_qdt_Cp HDF-EOS file for input
!    !=======================================================================
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg3_3d_mst_Cp_file )
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
!    CALL He4GridReadZ
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
!       !-----------------------------
!       ! Read data from HDF file
!       !-----------------------------
!
!       ! Zero data arrays
!       Q_P = 0e0
!       Q   = 0e0
!
!       ! Read data 
!       msg  = '%%% Reading                ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       CALL He4GridReadData( name, Q_P )
!          
!       ! Remove fill values
!       WHERE( Q_P == FILL_VALUE ) Q_P = 0e0
!       
!       !-----------------------------
!       ! P -> Hybrid regridding
!       !-----------------------------
!       msg  = '%%% P -> Hybrid interp for ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!       
!       ! Interpolate P -> hybrid
!       !$OMP PARALLEL DO         & 
!       !$OMP DEFAULT( SHARED )   &
!       !$OMP PRIVATE( I, J, T )
!       DO T = 1, TIMES_A3
!       DO J = 1, J125x125
!       DO I = 1, I125x125
!          CALL P2Hybrid72L( L125x125_P, PEDGE_P,      Q_P(I,J,:,T), &
!                            L125x125,   PLE(I,J,:,T), Q  (I,J,:,T) )
!       ENDDO
!       ENDDO
!       ENDDO
!       !$OMP END PARALLEL DO
!
!       !--------------------------------
!       ! Regrid data horizontally
!       ! (no level flipping needed)
!       !--------------------------------
!       msg  = '%%% Regridding             ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over the A3 times and # of levels
!       !$OMP PARALLEL DO                                    &
!       !$OMP DEFAULT( SHARED )                              &
!       !$OMP PRIVATE( T, L, ptr_125x125, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A3
!       DO L = 1, L1x125
!            
!          ! Point to time-slice of native data array
!          ptr_125x125 => Q(:,:,L,T)
!           
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,L,T,F+offset)
!             CALL RegridGeos57CTo2x25( 0, ptr_125x125, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,L,T,F+offset)
!             CALL RegridGeos57CTo4x5( 0, ptr_125x125, ptr_4x5 )
!          ENDIF
!       ENDDO
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
!    NULLIFY( ptr_125x125, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!  
!    msg = '%%%%%% LEAVING ROUTINE Process3dMstCp %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process3dMstCp
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Process3dQdtCp
!!
!! !DESCRIPTION: Subroutine Process3dChmFe regrids the Geos57 met fields from 
!!  the "tavg3\_3d\_qdt\_Cp" file and saves to the GEOS-Chem file format.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Process3dQdtCp( nFields, fields, offset )
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
!!  30 Jul 2010 - R. Yantosca - Initial version
!!  11 Aug 2010 - R. Yantosca - Now do not average to 6hr time resolution
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!
!    ! Scalars
!    INTEGER                 :: F, I, J, L, T
!
!    ! Arrays
!    REAL*4                  :: DqI( I125x125, J125x125, L125x125,   TIMES_A3 )
!    REAL*4                  :: DqL( I125x125, J125x125, L125x125,   TIMES_A3 )
!    REAL*4                  :: DqV( I125x125, J125x125, L125x125,   TIMES_A3 )
!    REAL*4                  :: P_P( I125x125, J125x125,             TIMES_A3 )
!    REAL*4                  :: PLE( I125x125, J125x125, L125x125+1, TIMES_A3 )
!    REAL*4                  :: Q_P( I125x125, J125x125, L125x125_P, TIMES_A3 )
!    REAL*4, TARGET          :: Q  ( I125x125, J125x125, L125x125,   TIMES_A3 )
!
!    ! Pointer arrays
!    REAL*4, POINTER         :: ptr_125x125(:,:)
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
!    ! Zero pointers and data arrays
!    NULLIFY( ptr_125x125, ptr_2x25, ptr_4x5 )
!
!    ! Echo info    
!    msg = '%%%%%% ENTERING ROUTINE Process3dQdtCp %%%%%%'
!    WRITE( IU_LOG, '(a)' ) '%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!    !=======================================================================
!    ! Read instantaneous surface pressure from the inst3_3d_asm_Cp
!    ! at 1.25 x 1.25 and compute the array of pressure at level edges 
!    !=======================================================================
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( inst3_3d_asm_Cp_file )
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
!    CALL He4GridReadZ
!    CALL He4GetNymdNhms
!
!    ! Read surface pressure data
!    P_P = 0e0
!    msg  = '%%% Reading PS on Cp grid'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridReadData( 'PS', P_P )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!
!    ! Make PLE data [hPa] for the hybrid grid from the surface pressure
!    CALL MakePLE( P_P, PLE )
!    
!    !=======================================================================
!    ! Open the tavg3_3d_qdt_Cp HDF-EOS file for input
!    !=======================================================================
!
!    ! Create filename from the template
!    fileHDF = TRIM( dataDirHDF ) // TRIM( tavg3_3d_qdt_Cp_file )
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
!    CALL He4GridReadZ
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
!       !----------------------------------------------------
!       ! MOISTQ = DQIDTMST + DQLDTMST + DQVDTMST
!       !
!       ! We assume that MOISTQ is at the end of the data
!       ! list, and that DQIDTMST, DQLDTMST, DQVDTMST will
!       ! be processed first.  We shall archive the sum of
!       ! the DQ*DTMST fields into the MSTQ array, which
!       ! will then be assigned to the Q array, which 
!       ! will be regridded to 2 x 2.5 and 4 x .5
!       !----------------------------------------------------
!       IF ( name == 'MOISTQ' ) THEN
!          
!          ! Echo info                        
!          msg  = '%%% Constructing           ' // name
!          WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!          ! Use previously vertically interpolated values for MOISTQ
!          ! Also convert MOISTQ from kg/kg/s to g/kg/day
!          ! for backwards compatibility with GEOS-Chem
!          ! Unit conversion = 1000 g/kg * 86400 s/day 
!          Q      = ( DqI + DqL + DqV ) * 1000e0 * 86400e0
!
!       ELSE
!
!          !-----------------------------
!          ! Read data from HDF file
!          !-----------------------------
!
!          ! Zero data arrays
!          Q_P = 0e0
!          Q   = 0e0
!
!          ! Read data 
!          msg  = '%%% Reading                ' // name
!          WRITE( IU_LOG, '(a)' ) TRIM( msg )
!          CALL He4GridReadData( name, Q_P )
!          
!          ! Remove fill values
!          WHERE( Q_P == FILL_VALUE ) Q_P = 0e0
!                 
!          !-----------------------------
!          ! P -> Hybrid regridding
!          !-----------------------------
!          msg  = '%%% P -> Hybrid interp for ' // name
!          WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!          ! Process 3-hr avg data
!          !$OMP PARALLEL DO         & 
!          !$OMP DEFAULT( SHARED )   &
!          !$OMP PRIVATE( I, J, T )
!          DO T = 1, TIMES_A3
!          DO J = 1, J125x125
!          DO I = 1, I125x125
!             CALL P2Hybrid72L( L125x125_P, PEDGE_P,      Q_P(I,J,:,T), &
!                               L125x125,   PLE(I,J,:,T), Q  (I,J,:,T) )
!          ENDDO
!          ENDDO
!          ENDDO
!          !$OMP END PARALLEL DO
!
!          !-----------------------------
!          ! Archive for MOISTQ
!          !-----------------------------
!          SELECT CASE( name )
!             CASE( 'DQIDTMST' )
!                DqI = Q
!             CASE( 'DQLDTMST' )
!                DqL = Q
!             CASE( 'DQVDTMST' ) 
!                DqV = Q
!             CASE DEFAULT
!                ! Do nothing
!          END SELECT
!
!       ENDIF
!
!       !--------------------------------
!       ! Regrid data horizontally
!       ! (no level flipping needed)
!       !--------------------------------
!       msg  = '%%% Regridding             ' // name
!       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!       ! Loop over the A3 times and # of levels
!       !$OMP PARALLEL DO                                      &
!       !$OMP DEFAULT( SHARED )                                &
!       !$OMP PRIVATE( T, L, ptr_125x125, ptr_2x25, ptr_4x5 )
!       DO T = 1, TIMES_A3
!       DO L = 1, L1x125
!            
!          ! Point to time-slice of native data array
!          ptr_125x125 => Q(:,:,L,T)
!           
!          ! Regrid to 2 x 2.5
!          IF ( do2x25 ) THEN
!             ptr_2x25 => Q2x25(:,:,L,T,F+offset)
!             CALL RegridGeos57CTo2x25( 0, ptr_125x125, ptr_2x25 )
!          ENDIF
!
!          ! Regrid to 4x5 
!          IF ( do4x5 ) THEN
!             ptr_4x5 => Q4x5(:,:,L,T,F+offset)
!             CALL RegridGeos57CTo4x5( 0, ptr_125x125, ptr_4x5 )
!          ENDIF
!       ENDDO
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
!    NULLIFY( ptr_125x125, ptr_2x25, ptr_4x5 )
!
!    ! Detach from grid and close HDF file
!    msg = '%%% Closing ' // TRIM( fileHDF )
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!    CALL He4GridClose( fileHDF )
!    CALL He4CleanUpIndexFields
!  
!    msg = '%%%%%% LEAVING ROUTINE Process3dQdtCp %%%%%%'
!    WRITE( IU_LOG, '(a)' ) TRIM( msg )
!
!  END SUBROUTINE Process3dQdtCp
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegridTau
!
! !DESCRIPTION: This routine regrids the GEOS-5 optical depth and cloud 
!  fraction fields from the 0.5 x 0.666 native resolution grid to a GEOS-Chem 
!  "coarse" grid (e.g. 1 x 1.25, 2 x 2.5, 4 x 5).
!
! !REMARKS:
!  Regridding formulae from Hongyu Liu:
!
!  Eq. 1: Optical depth:
!  ---------------------
!                                                                             .
!           TAU                   {     Fn * Wn            TAUn      }
!      -------------  =  A  =  SUM{ --------------- * -------------  }
!      ( TAU + 7.7 )              {  SUM( Fn * Wn )   ( TAUn + 7.7 ) }
!
!                                   |_____________|   |____________|
!     in the code these are:          Fn_Wn_Norm         optRatio
!                                                        tauiRatio
!                                                        tauwRAtio
!                                                                             .
!  and to solve for TAU from this equation, we manipulate the right hand 
!  side of the equation (called A), as follows:
!                                                                             .
!                                A
!                 TAU = 7.7 * --------
!                              1 - A
!                                                                             .
!  Eq. 2. Cloud fraction:
!  ----------------------
!                                                                             .
!                        SUM( Fn * Wn )
!                  F  =  ----------------   
!                           SUM( Wn )
!                                                                             .
!  In both above equations:
!                                                                             .
!     TAU    = In-cloud optical depth regridded to "coarse" grid box
!     TAUn   = In-cloud optical depth in a 0.5 x 0.666 grid box
!                                                                             .
!     F      = Cloud fraction regriddded to "coarse" grid box
!     Fn     = Cloud fraction in a 0.5 x 0.666 grid box
!                                                                             .
!     Wn     = Fraction of 0.5 x 0.666 box in the "coarse" grid box
!
! !INTERFACE:
!
  SUBROUTINE RegridTau( optIn,  optOut,  tauiIn, tauiOut, &
                        tauwIn, tauwOut, fIn,    fOut,    &
                        map,    IMX,     JMX,    LMX )
!
! !INPUT PARAMETERS: 
!
    ! Dimensions of coarse grid
    INTEGER,      INTENT(IN)  :: IMX, JMX, LMX    

    ! Mapping weight object
    TYPE(MapObj), POINTER     :: map(:,:)

    ! Input total in-cloud optical depth (= water OD + ice OD)
    REAL*4,       INTENT(IN)  :: optIn  ( I025x03125, J025x03125, L025x03125 ) 

    ! Input in-cloud ice optical depth
    REAL*4,       INTENT(IN)  :: tauiIn ( I025x03125, J025x03125, L025x03125 )

    ! Input in-cloud water optical depth
    REAL*4,       INTENT(IN)  :: tauwIn ( I025x03125, J025x03125, L025x03125 )
 
    ! Input cloud fraction
    REAL*4,       INTENT(IN)  :: fIn    ( I025x03125, J025x03125, L025x03125 )
!
! !OUTPUT PARAMETERS:
!
    ! Output total in-cloud optical depth
    REAL*4,       INTENT(OUT) :: optOut ( IMX, JMX, LMX )

    ! Output in-cloud ice optical depth
    REAL*4,       INTENT(OUT) :: tauiOut( IMX, JMX, LMX )

    ! Output in-cloud water path optical depth
    REAL*4,       INTENT(OUT) :: tauwOut( IMX, JMX, LMX )

    ! Output cloud fraction
    REAL*4,       INTENT(OUT) :: fOut   ( IMX, JMX, LMX )
!
! !REVISION HISTORY: 
!  29 Jul 2010 - R. Yantosca - Initial version, based on GEOS-5
!  02 Aug 2010 - R. Yantosca - Now flip output arrays in vertical
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: I,         J,         L,        LR,   nPoints
    INTEGER :: Nx,        Ny,        X,        Y
    REAL*4  :: sum_Fn_Wn, sum_Wn,    Fn_Wn
    REAL*4  :: optRatio,  tauiRatio, tauwRatio
    REAL*4  :: optRHS,    tauiRHS,   tauwRHS
    
    ! Loop over coarse grid boxes
    !$OMP PARALLEL DO                                                  &
    !$OMP DEFAULT( SHARED )                                            &
    !$OMP PRIVATE( I,         J,         L,       nPoints, sum_Fn_Wn ) &
    !$OMP PRIVATE( sum_Wn,    optRHS,    tauiRHS, tauwRHS, Nx        ) &
    !$OMP PRIVATE( Ny,        X,         Y,       Fn_Wn,   optRatio  ) &
    !$OMP PRIVATE( tauiRatio, tauwRatio, LR                          )
    DO L = 1, LMX

       ! Reverse level index for output arrays
       LR = LMX - L + 1

    DO J = 1, JMX
    DO I = 1, IMX

       ! Number of "fine" grid boxes in each dimension
       ! that comprise a "coarse" grid box
       nPoints = map(I,J)%nPoints

       !---------------------------------
       ! Regrid cloud fraction & OD
       !---------------------------------

       ! Zero summing variables
       sum_Fn_Wn = 0e0
       sum_Wn    = 0e0
       optRHS    = 0e0
       tauiRHS   = 0e0
       tauwRHS   = 0e0
       
       ! Loop over "fine" grid boxes
       DO Ny = 1, nPoints
       DO Nx = 1, nPoints
          
          ! Avoid useless clock cycles if the mapping weight is zero
          IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN

             ! Indices of each "fine" grid box that makes up the "coarse" box
             X          = map(I,J)%xInd(Nx)
             Y          = map(I,J)%yInd(Ny)

             ! Sum of the mapping weights over all of the "fine" grid
             ! boxes (X,Y) that make up the "coarse" grid box (I,J)
             sum_Wn     = sum_Wn    + map(I,J)%weight(Nx,Ny)

             ! Sum of the cloud fraction * mapping weights over all of the 
             ! "fine" grid boxes (X,Y) that make up the "coarse" grid box (I,J)
             sum_Fn_Wn  = sum_Fn_Wn + ( fIn(X,Y,L) * map(I,J)%weight(Nx,Ny) )

             ! Compute the cloud fraction * mapping weight
             Fn_Wn      = fIn(X,Y,L) * map(I,J)%weight(Nx,Ny)
        
             ! Compute the term TAU / ( TAU + 7.7 ) for the total,
             ! ice-path, and water-path opbical depths.  NOTE: We don't have 
             ! to worry about div by zero due to the +7.7 in the denominator.
             optRatio   = optIn(X,Y,L)  / ( optIn(X,Y,L)  + 7.7e0 )
             tauiRatio  = tauiIn(X,Y,L) / ( tauiIn(X,Y,L) + 7.7e0 )
             tauwRatio  = tauwIn(X,Y,L) / ( tauwIn(X,Y,L) + 7.7e0 )
          
             ! Compute the right hand side of the equation for the total OD, 
             ! ice-path OD, and water-path OD.  For computational expediency, 
             ! we'll divide by SUM( Fn * Wn ) outside this DO loop.
             optRHS     = optRHS  + ( Fn_Wn * optRatio  )
             tauiRHS    = tauiRHS + ( Fn_Wn * tauiRatio )
             tauwRHS    = tauwRHS + ( Fn_Wn * tauwRatio )
          ENDIF
       ENDDO
       ENDDO

       ! Output cloud fraction.  NOTE, we don't have to worry about div by
       ! zero since SUM( Wn ) will always be greater than zero (there is 
       ! always at least 1 "fine" small box in the "coarse" box).
       fOut(I,J,LR) = sum_Fn_Wn / sum_Wn

       ! Total optical depth on the coarse grid
       IF ( IsSafeDiv( optRHS, sum_Fn_Wn-optRHS ) ) THEN
          optOut(I,J,LR) = 7.7e0 * optRHS / ( sum_Fn_Wn - optRHS )
       ELSE
          optOut(I,J,LR) = 0e0
       ENDIF

       ! Ice-path optical depth on the coarse grid
       IF ( IsSafeDiv( tauiRHS, sum_Fn_Wn-tauiRHS ) ) THEN
          tauiOut(I,J,LR) = 7.7e0 * tauiRHS / ( sum_Fn_Wn - tauiRHS )
       ELSE
          tauiOut(I,J,LR) = 0e0
       ENDIF

       ! Water-path optical depth on the coarse grid
       IF ( IsSafeDiv( tauwRHS, sum_Fn_Wn-tauwRHS ) ) THEN
          tauwOut(I,J,LR) = 7.7e0 * tauwRHS / ( sum_Fn_Wn - tauwRHS )
       ELSE
          tauwOut(I,J,LR) = 0e0
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE RegridTau
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsSafeDiv
!
! !DESCRIPTION: Function IsSafeDiv returns TRUE if the numerator N and 
!  denominator D may be divided safely (i.e. without resulting in a 
!  division-by-zero, Not-a-Number (NaN), or Infinity), or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsSafeDiv( N, D ) RESULT( isSafe )
!
! !INPUT PARAMETERS: 
!
    REAL*4, INTENT(IN) :: N        ! Numerator
    REAL*4, INTENT(IN) :: D        ! Denominator
!
! !RETURN VALUE:
!    
    LOGICAL            :: isSafe   ! Returns TRUE if it's safe to divide N/D 
!
! !REVISION HISTORY: 
!   29 Sep 2008 - R. Yantosca - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: MaxExp, MinExp

    ! Maxinum 
    MaxExp = MAXEXPONENT( N )
    MinExp = MINEXPONENT( N )
    
    ! Test if it's safe to divide 
    IF ( ( D                         == 0      )  .or. &
         ( EXPONENT(N) - EXPONENT(D) >= MaxExp )  .or. &
         ( EXPONENT(N) - EXPONENT(D) <= MinExp ) ) THEN
       isSafe = .FALSE.
    ELSE
       isSafe = .TRUE.
    ENDIF

  END FUNCTION IsSafeDiv
!EOC
END MODULE Geos57A3Module

