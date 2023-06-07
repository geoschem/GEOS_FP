!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpA3MstCModule
!
! !DESCRIPTION: Module GeosFpA3MstCModule contains routines to create the
!  GEOS-Chem average 3-hr data files (moist parameters on level centers)
!  from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE:

MODULE GeosFpA3MstCModule
!
! !USES:
!
  ! GEOS-5.7.x data modules
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
  PUBLIC  :: GeosFpMakeA3MstC
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dMstNv
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  12 Jan 2012 - R. Yantosca - Now just save out fields on level centers
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Renamed to GeosFpA3MstCModule.F90
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
!  01 Feb 2012 - R. Yantosca - Make all global attribute names lowercase
!  20 Sep 2013 - R. Yantosca - Change and/or add attributes for COARDS
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
    LOGICAL            :: is_nc4
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

    !Output nc4 now to accomodate large file size for global 0.25x0.3125 data
    is_nc4 = .TRUE.

    ! Open netCDF file for writing
    CALL NcCr_Wr( fOut, TRIM( outFileName ), WRITE_NC4=is_nc4 )

    ! Turn filling off
    CALL NcSetFill( fOut, NF_NOFILL, omode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------

    ! Title string
    lName = 'GEOS-FP time-averaged 3-hour moist parameters on level centers (A3mstC), processed for GEOS-Chem input'
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
    lName = "NetCDF-4" ; !(jxu, 2015/09/004, convert nc3 to nc4)
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
       CASE( 'native', 'nested AS', 'nested EU', 'nested NA', '0.25 x 0.3125 global')
          DI = '0.3125'
          DJ = '0.25'
       CASE( 'nested AS 05', 'nested EU 05', 'nested NA 05')
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
    cal     = 'gregorian'
    lName   = 'time'
    units   = UnitsForTime( yyyymmdd )
    delta_t = '0000-00-00 03:00:00'
    begin_d = yyyymmdd_string
    begin_t = '000000'
    incr    = '030000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',       TRIM( cal     )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Level index array
    var1    = (/ idLev /)
    lName   = 'levels'
    units   = '1'
    CALL NcDef_Variable      ( fOut, 'lev', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    )

    ! Latitude index array
    var1    = (/ idLat /)
    lName   = 'latitude'
    units   = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    )

    ! Longitude index array
    var1    = (/ idLon /)
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! DQRCU
    IF ( StrPos( 'DQRCU', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Precipitation production rate -- convective'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'DQRCU', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! DQRLSAN
    IF ( StrPos( 'DQRLSAN', tavg3_3d_mst_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Precipitation production rate -- large scale + anvil'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'DQRLSAN', NF_FLOAT, 4, var4, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FPCU
    IF ( StrPos( 'FPCU', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Fraction of box undergoing convective precipitation'
       units = 'fraction'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'FPCU', NF_FLOAT, 4, var4, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! FPLSAN
    IF ( StrPos( 'FPLSAN', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Fraction of box undergoing large scale + anvil precipitation'
       units = 'fraction'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'FPLSAN', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! REEVAPCN
    IF ( StrPos( 'REEVAPCN', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Evaporation of precipitating convective condensate'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'REEVAPCN', NF_FLOAT, 4, var4, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! REEVAPLS(AN)
    IF ( StrPos( 'REEVAPLS', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)
       
       lName = 'Evaporation of precipitating large-scale & anvil condensate'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'REEVAPLS', NF_FLOAT, 4, var4, vId   )
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
    CALL NcWr( xMid, fOut, 'lon',  (/ 1 /), (/ X   /) )
    CALL NcWr( yMid, fOut, 'lat',  (/ 1 /), (/ Y   /) )
    CALL NcWr( zMid, fOut, 'lev',  (/ 1 /), (/ Z   /) )
    CALL NcWr( time, fOut, 'time', (/ 1 /), (/ T   /) )

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
! !IROUTINE: GeosFpMakeA3MstC
!
! !DESCRIPTION: Routine GeosFpMakeA3MstC is the the driver routine for
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (moist parameters on level
!       centers) from the GEOS-FP raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeA3MstC()
!
! !REVISION HISTORY:
!  11 Aug 2010 - R. Yantosca - Initial version, based on GeosFpA6Module.F90
!  12 Jan 2012 - R. Yantosca - Now just process fields on level centers
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Now save output to nested Europe grid
!  23 Sep 2013 - R. Yantosca - Now define netCDF latitude such that the poles
!                              are at -90/+90.  This facilitates the GIGC
!                              using ESMF/MAPL.
!  08 Oct 2013 - R. Yantosca - Now save output to nested SE grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_3dMstNv
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dMstNv(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeA3MstC %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg3_3d_mst_Nv_data )

    ! Return the list of fields and number of fields to process
    ! from each of the GeosFp raw met data files
    CALL GetNFields( tavg3_3d_mst_Nv_data,  nFields_3dMstNv, fields_3dMstNv )
    CALL GetNFields( allFieldsList,         nAllFields,      allFields      )

    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_mst_Nv_file ), nFields_3dMstNv
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested AS output file
    IF ( doNestAs ) THEN
       fName = TRIM( tempDirTmplNestAs ) // TRIM( dataTmplNestAs )
       gName = 'nested IN'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I_NestAs,  J_NestAs,     L025x03125, TIMES_A3,  &
                          xMid_025x03125(I0_as:I1_as),                    &
                          yMid_025x03125(J0_as:J1_as),                    &
                          zMid_025x03125,                      a3Mins,    &
                          gName,     fName,        fOutNestAs            )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu ) THEN
       fName = TRIM( tempDirTmplNestEu ) // TRIM( dataTmplNestEu )
       gName = 'nested EU'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I_NestEu,  J_NestEu,     L025x03125, TIMES_A3,  &
                          xMid_025x03125(I0_eu:I1_eu),                    &
                          yMid_025x03125(J0_eu:J1_eu),                    &
                          zMid_025x03125,                      a3Mins,    &
                          gName,     fName,        fOutNestEu            )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa ) THEN
       fName = TRIM( tempDirTmplNestNa ) // TRIM( dataTmplNestNa )
       gName = 'nested NA'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I_NestNa,  J_NestNa,     L025x03125, TIMES_A3,  &
                          xMid_025x03125(I0_na:I1_na),                    &
                          yMid_025x03125(J0_na:J1_na),                    &
                          zMid_025x03125,                      a3Mins,    &
                          gName,     fName,        fOutNestNa            )
    ENDIF

    ! Open 0.5x0.625 output file
    IF ( doGlobal05 ) THEN
      fName = TRIM( tempDirTmplGlobal05 ) // TRIM( dataTmplGlobal05 )
      gName = '0.5x0.625 global'
      CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
      CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
      CALL NcOutFileDef( I05x0625,     J05x0625,        L05x0625,      TIMES_A3,  &
                         xMid_05x0625, yMid_05x0625, zMid_05x0625,  a3Mins,    &
                         gName,     fName,        fOutGlobal05             )
   ENDIF

    ! Open 0.25x0.3125 output file
    IF ( do025x03125 ) THEN
       fName = TRIM( tempDirTmpl025x03125 ) // TRIM( dataTmpl025x03125 )
       gName = '0.25x0.3125 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I025x03125,     J025x03125,        L025x03125,      TIMES_A3,  &
                          xMid_025x03125, yMid_025x03125, zMid_025x03125,  a3Mins,    &
                          gName,     fName,        fOut025x03125              )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25,      TIMES_A3,  &
                          xMid_2x25, nc_yMid_2x25, zMid_2x25,  a3Mins,    &
                          gName,     fName,        fOut2x25              )
    ENDIF

    ! Open 4 x 5 output file
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'              )
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5,       TIMES_A3,  &
                          xMid_4x5,  nc_yMid_4x5,  zMid_4x5,   a3Mins,    &
                          gName,     fName,        fOut4x5               )
    ENDIF

    ! Open nested IN output file
    IF ( doNestAs05 ) THEN
       fName = TRIM( tempDirTmplNestAs05 ) // TRIM( dataTmplNestAs05 )
       gName = 'nested IN 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'    )
       CALL NcOutFileDef( I_NestAs05,  J_NestAs05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_as05:I1_as05),          &
                          yMid_05x0625(J0_as05:J1_as05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestAs         )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'    )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,       fOut05NestEu          )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3mstC'    )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestNa          )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process3dMstNv( nFields_3dMstNv, fields_3dMstNv ) ! tavg3_3d_mst_Nv

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A3mstC output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestAs    ) CALL NcCl( fOutNestAs    )
    IF ( doNestEu    ) CALL NcCl( fOutNestEu    )
    IF ( doNestNa    ) CALL NcCl( fOutNestNa    )
    IF ( do025x03125 ) CALL NcCl( fOut025x03125 )
    IF ( doNestAs05  ) CALL NcCl( fOut05NestAs  )
    IF ( doNestEu05  ) CALL NcCl( fOut05NestEu  )
    IF ( doNestNa05  ) CALL NcCl( fOut05NestNa  )
    IF ( doGlobal05  ) CALL NcCl( fOutGlobal05  )
    IF ( do2x25      ) CALL NcCl( fOut2x25      )
    IF ( do4x5       ) CALL NcCl( fOut4x5       )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeA3MstC %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeA3MstC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dMstNv
!
! !DESCRIPTION: Subroutine Process3dMstNv regrids the GeosFp met fields from
!  the "tavg3\_3d\_mst\_Nv" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dMstNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY:
!  09 Jan 2012 - R. Yantosca - Initial version
!  10 Jan 2012 - R. Yantosca - Activate parallel loop over vertical levels
!  17 Jan 2012 - R. Yantosca - Bug fix: flip data in vertical immediately
!                              after reading.
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Now save output to nested Europe grid (EU)
!  08 Oct 2013 - R. Yantosca - Now save output to nested SE Asia grid (SE)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: F,       H,       L,       LR
    INTEGER                 :: hhmmss

    ! Variables for netCDF I/O
    INTEGER                 :: X,       Y,       Z,       T
    INTEGER                 :: XNestAs, YNestAs, ZNestAs, TNestAs
    INTEGER                 :: XNestEu, YNestEu, ZNestEu, TNestEu
    INTEGER                 :: XNestNa, YNestNa, ZNestNa, TNestNa
    INTEGER                 :: X025x03125, Y025x03125, Z025x03125, T025x03125
    INTEGER                 :: X05x0625, Y05x0625, Z05x0625, T05x0625
    INTEGER                 :: X2x25,   Y2x25,   Z2x25,   T2x25
    INTEGER                 :: X4x5,    Y4x5,    Z4x5,    T4x5
    INTEGER                 :: st4d(4), ct4d(4)
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: FC   ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: FL   ( I025x03125, J025x03125, L025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5,       L4x5       )

    INTEGER                 :: XNestAs05,  YNestAs05, ZNestAs05, TNestAs05
    INTEGER                 :: XNestEu05,  YNestEu05, ZNestEu05, TNestEu05
    INTEGER                 :: XNestNa05,  YNestNa05, ZNestNa05, TNestNa05
    REAL*4,  TARGET         :: Q05    ( I05x0625, J05x0625, L05x0625 )

    ! Pointer arrays
    REAL*4, POINTER         :: Qflip(:,:,:)
    REAL*4, POINTER         :: Ptr  (:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=10      ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

     ! Echo info
    msg = '%%%%%% ENTERING ROUTINE Process3dMstNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Nested AS grid
    IF ( doNestAs ) THEN
       CALL NcGet_DimLen( fOutNestAs, 'lon',  XNestAs )
       CALL NcGet_DimLen( fOutNestAs, 'lat',  YNestAs )
       CALL NcGet_DimLen( fOutNestAs, 'lev',  ZNestAs )
       CALL NcGet_DimLen( fOutNestAs, 'time', TNestAs )
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

    ! 0.25x0.3125 global grid
    IF ( do025x03125 ) THEN
       CALL NcGet_DimLen( fOut025x03125,   'lon',  X025x03125   )
       CALL NcGet_DimLen( fOut025x03125,   'lat',  Y025x03125   )
       CALL NcGet_DimLen( fOut025x03125,   'lev',  Z025x03125   )
       CALL NcGet_DimLen( fOut025x03125,   'time', T025x03125   )
    ENDIF

    ! 0.5x0.625 global grid
    IF ( doGlobal05 ) THEN
      CALL NcGet_DimLen( fOutGlobal05,   'lon',  X05x0625   )
      CALL NcGet_DimLen( fOutGlobal05,   'lat',  Y05x0625   )
      CALL NcGet_DimLen( fOutGlobal05,   'lev',  Z05x0625   )
      CALL NcGet_DimLen( fOutGlobal05,   'time', T05x0625   )
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

    ! Nested AS grid
    IF ( doNestAs05 ) THEN
       CALL NcGet_DimLen( fOut05NestAs, 'lon',  XNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lat',  YNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lev',  ZNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'time', TNestAs05 )
    ENDIF

    ! Nested EU grid
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

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_mst_Nv_file )
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

          ! Skip certain field names
          SELECT CASE( name )
             CASE( 'FPCU', 'FPLSAN' )    ! These are derived, not read
                CYCLE
             CASE DEFAULT
                ! Nothing
          END SELECT

          ! Save 1st 8 characters of NAME for netCDF output
          name8 = TRIM( name )

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

          ! Flip levels in the vertical
          QFlip => Q( :, :, Z:1:-1 )

          !-----------------------------------------------------------------
          ! Regrid data fields
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name8
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Loop over the A-3 times and vertical levels
          DO L = 1, Z

             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridGeosFpTo2x25( 0, Qflip(:,:,L), Q2x25(:,:,L) )
             ENDIF

             ! Regrid to 4 x 5
             IF ( do4x5 ) THEN
                CALL RegridGeosFpTo4x5 ( 0, Qflip(:,:,L), Q4x5(:,:,L)  )
             ENDIF
             ! Regrid to 0.5 x 0.625
             IF ( do05x0625 ) THEN
                CALL RegridGeosFpTo05x0625( 0, Qflip(:,:,L), Q05(:,:,L) )
             ENDIF

          ENDDO

          !-----------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------
          msg = '%%% Archiving  ' // name8
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Nested AS (point to proper slice of global data)
          IF ( doNestAs ) THEN
             Ptr  => Qflip( I0_as:I1_as, J0_as:J1_as, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestAs, YNestAs, ZNestAs, 1 /)
             CALL NcWr( Ptr, fOutNestAs, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu ) THEN
             Ptr  => Qflip( I0_eu:I1_eu, J0_eu:J1_eu, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestEu, YNestEu, ZNestEu, 1 /)
             CALL NcWr( Ptr, fOutNestEu, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa ) THEN
             Ptr  => Qflip( I0_na:I1_na, J0_na:J1_na, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestNa, YNestNa, ZNestNa, 1 /)
             CALL NcWr( Ptr, fOutNestNa, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Global 0.25x0.3125
          IF ( do025x03125 ) THEN
             Ptr  => Qflip
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ X025x03125, Y025x03125, Z025x03125, 1 /)
             CALL NcWr( Ptr, fOut025x03125, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Global 0.5x0.625
          IF ( doGlobal05 ) THEN
            Ptr  => Q05
            st4d = (/ 1,       1,       1,       H /)
            ct4d = (/ X05x0625, Y05x0625, Z05x0625, 1 /)
            CALL NcWr( Ptr, fOutGlobal05, TRIM( name8 ), st4d, ct4d )
            NULLIFY( Ptr )
         ENDIF

          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st4d = (/ 1,     1,     1,     H  /)
             ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name8 ), st4d, ct4d )
          ENDIF

          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st4d  = (/ 1,    1,    1,    H /)
             ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name8 ), st4d, ct4d )
          ENDIF

          ! Nested AS (point to proper slice of global data)
          IF ( doNestAs05 ) THEN
             Ptr  => Q05( I0_as05:I1_as05, J0_as05:J1_as05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestAs05, YNestAs05, ZNestAs05, 1 /)
             CALL NcWr( Ptr, fOut05NestAs, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu05 ) THEN
             Ptr  => Q05( I0_eu05:I1_eu05, J0_eu05:J1_eu05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestEu05, YNestEu05, ZNestEu05, 1 /)
             CALL NcWr( Ptr, fOut05NestEu, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa05 ) THEN
             Ptr  => Q05( I0_na05:I1_na05, J0_na05:J1_na05, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestNa05, YNestNa05, ZNestNa05, 1 /)
             CALL NcWr( Ptr, fOut05NestNa, TRIM( name8 ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Free pointer memory
          NULLIFY( Qflip )
       ENDDO

       !====================================================================
       ! Close input file
       !====================================================================
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================

    ! Echo info
    msg = '%%%%%% LEAVING ROUTINE Process3dMstNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dMstNv
!EOC
END MODULE GeosFpA3MstCModule
