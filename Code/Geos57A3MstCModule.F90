!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57A3MstCModule
!
! !DESCRIPTION: Module Geos57A3MstCModule contains routines to create the 
!  GEOS-Chem average 3-hr data files (moist parameters on level centers) 
!  from the GEOS-5.7.2 raw data.
!\\
!\\
! !INTERFACE: 

MODULE Geos57A3MstCModule
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
  PUBLIC  :: Geos57MakeA3MstC
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dMstNv
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Prior to 1/19/12: 
!% Abandon this line of development, until further notice
!%  PRIVATE :: ProcessFracPrecip
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  09 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  12 Jan 2012 - R. Yantosca - Now just save out fields on level centers
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
    lName = 'GEOS-5.7.2 time-averaged 3-hour moist parameters on level centers (A3mstC) fields for GEOS-Chem'
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

    ! DQRCU
    IF ( StrPos( 'DQRCU', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Precipitation production rate -- convective'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'DQRCU', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! DQRLSAN
    IF ( StrPos( 'DQRLSAN', tavg3_3d_mst_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Precipitation production rate -- large scale + anvil'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'DQRLSAN', NF_FLOAT, 4, var4, vId    )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FPCU
    IF ( StrPos( 'FPCU', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of box undergoing convective precipitation'
       units = 'fraction'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'FPCU', NF_FLOAT, 4, var4, vId       )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! FPLSAN
    IF ( StrPos( 'FPLSAN', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Fraction of box undergoing large scale + anvil precipitation'
       units = 'fraction'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'FPLSAN', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! REEVAPCN
    IF ( StrPos( 'REEVAPCN', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Evaporation of precipitating convective condensate'
       units = 'kg kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'REEVAPCN', NF_FLOAT, 4, var4, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
    ENDIF

    ! REEVAPLS(AN)
    IF ( StrPos( 'REEVAPLS', tavg3_3d_mst_Nv_data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Evaporation of precipitating large-scale & anvil condensate'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'REEVAPLS', NF_FLOAT, 4, var4, vId   )
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
! !IROUTINE: Geos57MakeA3MstC
!
! !DESCRIPTION: Routine Geos57MakeA3MstC is the the driver routine for 
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (moist parameters on level
!       centers) from the GEOS-5.7.2 raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program Geos57Driver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57MakeA3MstC()
!
! !REVISION HISTORY: 
!  11 Aug 2010 - R. Yantosca - Initial version, based on Geos57A6Module.F90
!  12 Jan 2012 - R. Yantosca - Now just process fields on level centers
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
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
    msg = '%%%%%%%%%% ENTERING ROUTINE Geos57MakeA3MstC %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg3_3d_mst_Nv_data ) 

    ! Return the list of fields and number of fields to process
    ! from each of the Geos57 raw met data files
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

    ! Open nested China output file
    IF ( doNestCh ) THEN
       fName = TRIM( tempDirTmplNestCh ) // TRIM( dataTmplNestCh )
       gName = 'SEA4CRS'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstC'              )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,  L025x03125, TIMES_A3,  &
                          xMid_025x03125(I0_ch:I1_ch),                 &
                          yMid_025x03125(J0_ch:J1_ch),                 &
                          zMid_025x03125,                   a3Mins,    &
                          gName,     fName,     fOutNestCh            )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstC'              )
       CALL NcOutFileDef( I2x25,     J2x25,     L2x25,      TIMES_A3,  &
                          xMid_2x25, yMid_2x25, zMid_2x25,  a3Mins,    &
                          gName,     fName,     fOut2x25              )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                )      
       CALL StrRepl     ( fName,     '%%%%%%',  'A3mstC'              )
       CALL NcOutFileDef( I4x5,      J4x5,      L4x5,       TIMES_A3,  &
                          xMid_4x5,  yMid_4x5,  zMid_4x5,   a3Mins,    &
                          gName,     fName,     fOut4x5               )
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
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Geos57MakeA3MstC %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'   
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Geos57MakeA3MstC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dMstNv
!
! !DESCRIPTION: Subroutine Process3dMstNv regrids the Geos57 met fields from 
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
    INTEGER                 :: XNestCh, YNestCh, ZNestCh, TNestCh
    INTEGER                 :: X2x25,   Y2x25,   Z2x25,   T2x25
    INTEGER                 :: X4x5,    Y4x5,    Z4x5,    T4x5
    INTEGER                 :: st4d(4), ct4d(4)
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: FC   ( I025x03125, J025x03125, L025x03125 )
    REAL*4,  TARGET         :: FL   ( I025x03125, J025x03125, L025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5,       L4x5       )

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
#if defined( USE_SAMPLE_DATA )
       if ( hhmmss > 193000 ) CYCLE
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Prior to 1/19/12:
!% Abandon this line of development, until further notice (bmy, 1/19/12)
!%          !-----------------------------------------------------------------
!%          ! Pre-regrid handling (flip levels in the vertical)
!%          !-----------------------------------------------------------------
!%          SELECT CASE( name )
!%             CASE( 'DQRCU'   )
!%                FC = QFlip
!%             CASE( 'DQRLSAN' )
!%                FL = QFlip
!%             CASE DEFAULT
!%                ! Nothing
!%          END SELECT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
          !-----------------------------------------------------------------
          ! Regrid data fields 
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name8
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Loop over the A-3 times and vertical levels
          DO L = 1, Z

             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridGeos57To2x25( 0, Qflip(:,:,L), Q2x25(:,:,L) )
             ENDIF

             ! Regrid to 4 x 5
             IF ( do4x5 ) THEN
                CALL RegridGeos57To4x5 ( 0, Qflip(:,:,L), Q4x5(:,:,L)  )
             ENDIF       

          ENDDO
                
          !-----------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------
          msg = '%%% Archiving  ' // name8
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Nested China (point to proper slice of global data)
          ! Flip levels in the vertical
          IF ( doNestCh ) THEN
             Ptr  => Qflip( I0_ch:I1_ch, J0_ch:J1_ch, : )
             st4d = (/ 1,       1,       1,       H /)
             ct4d = (/ XNestCh, YNestCh, ZNestCh, 1 /)
             CALL NcWr( Ptr, fOutNestCh, TRIM( name8 ), st4d, ct4d )
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

          ! Free pointer memory
          NULLIFY( Qflip )
       ENDDO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Prior to 1/19/12:
!% Abandon this line of development, until further notice (bmy, 1/19/12)
!%       !====================================================================
!%       ! Process FPCU (fraction of box undergoing conv precip)
!%       !====================================================================
!%
!%       !---------------------------------
!%       ! Create FPCU from DQRCU
!%       !---------------------------------
!%       msg = '%%% Creating   FPCU'
!%       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!%
!%       ! Native grid
!%       IF ( doNative ) THEN 
!%          CALL ProcessFracPrecip( FC,         Q,          mapNative, &
!%                                  I025x03125, J025x03125, L025x03125 )
!%       ENDIF
!%
!%       ! 2 x 2.5 grid
!%       IF ( do2x25 ) THEN
!%          CALL ProcessFracPrecip( FC, Q2x25, mapTo2x25, I2x25, J2x25, L2x25 )
!%       ENDIF
!%
!%       ! 4 x 5 grid
!%       IF ( do4x5 ) THEN 
!%          CALL ProcessFracPrecip( FC, Q4x5,  mapTo4x5,  I4x5,  J4x5,  L4x5  )
!%       ENDIF
!%
!%       !---------------------------------
!%       ! Write FPCU to netCDF
!%       !---------------------------------
!%       msg = '%%% Archiving  FPCU'
!%       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!%
!%       ! Nested China (point to proper slice of global data)
!%       IF ( doNestCh ) THEN
!%          Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch, : )
!%          st4d = (/ 1,       1,       1,       H /)
!%          ct4d = (/ XNestCh, YNestCh, ZNestCh, 1 /)
!%          CALL NcWr( Ptr, fOutNestCh, 'FPCU', st4d, ct4d )
!%          NULLIFY( Ptr )
!%       ENDIF
!%       
!%       ! Write 2 x 2.5 data
!%       IF ( do2x25 ) THEN 
!%          st4d = (/ 1,     1,     1,     H  /)
!%          ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
!%          CALL NcWr( Q2x25, fOut2x25, 'FPCU', st4d, ct4d )
!%       ENDIF
!%       
!%       ! Write 4x5 data
!%       IF ( do4x5 ) THEN
!%          st4d  = (/ 1,    1,    1,    H /)
!%          ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
!%          CALL NcWr( Q4x5, fOut4x5, 'FPCU', st4d, ct4d )
!%       ENDIF
!%                
!%       !====================================================================
!%       ! Process FPLSAN (fraction of box undergoiong LS+anvil precip)
!%       !====================================================================
!%
!%       !---------------------------------
!%       ! Create FPLSAN from DQRLSAN
!%       !---------------------------------
!%       msg = '%%% Creating   FPLSAN'
!%       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!%
!%       ! Native grid
!%       IF ( doNative ) THEN 
!%          CALL ProcessFracPrecip( FL,         Q,          mapNative, &
!%                                  I025x03125, J025x03125, L025x03125 )
!%       ENDIF
!%
!%       ! 2 x 2.5 grid
!%       IF ( do2x25 ) THEN
!%          CALL ProcessFracPrecip( FL, Q2x25, mapTo2x25, I2x25, J2x25, L2x25 )
!%       ENDIF
!%
!%       ! 4 x 5 grid
!%       IF ( do4x5 ) THEN 
!%          CALL ProcessFracPrecip( FL, Q4x5,  mapTo4x5,  I4x5,  J4x5,  L4x5  )
!%       ENDIF
!%
!%       !---------------------------------
!%       ! Write FPLSAN to netCDF
!%       !---------------------------------
!%       msg = '%%% Archiving  FPCU'
!%       WRITE( IU_LOG, '(a)' ) TRIM( msg )
!%
!%       ! Nested China (point to proper slice of global data)
!%       IF ( doNestCh ) THEN
!%          Ptr  => Q( I0_ch:I1_ch, J0_ch:J1_ch, : )
!%          st4d = (/ 1,       1,       1,       H /)
!%          ct4d = (/ XNestCh, YNestCh, ZNestCh, 1 /)
!%          CALL NcWr( Ptr, fOutNestCh, 'FPLSAN', st4d, ct4d )
!%          NULLIFY( Ptr )
!%       ENDIF
!%       
!%       ! Write 2 x 2.5 data
!%       IF ( do2x25 ) THEN
!%          st4d = (/ 1,     1,     1,     H  /)
!%          ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
!%          CALL NcWr( Q2x25, fOut2x25, 'FPLSAN', st4d, ct4d )
!%       ENDIF
!%       
!%       ! Write 4x5 data
!%       IF ( do4x5 ) THEN
!%          st4d  = (/ 1,    1,    1,    H /)
!%          ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
!%          CALL NcWr( Q4x5, fOut4x5, 'FPLSAN', st4d, ct4d )
!%       ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Prior to 1/19/12:
!% Abandon this line of development, until further notice (bmy, 1/19/12)
!%!------------------------------------------------------------------------------
!%!          Harvard University Atmospheric Chemistry Modeling Group            !
!%!------------------------------------------------------------------------------
!%!BOP
!%!
!%! !IROUTINE: RegridTau
!%!
!%! !DESCRIPTION: This routine regrids the GEOS-5 optical depth and cloud 
!%!  fraction fields from the 0.5 x 0.666 native resolution grid to a GEOS-Chem 
!%!  "coarse" grid (e.g. 1 x 1.25, 2 x 2.5, 4 x 5).
!%!
!%! !INTERFACE:
!%!
!%  SUBROUTINE ProcessFracPrecip( fpIn, fpOut, map, IMX, JMX, LMX )
!%!
!%! !INPUT PARAMETERS: 
!%!
!%    ! Dimensions of coarse grid
!%    INTEGER,      INTENT(IN)  :: IMX, JMX, LMX    
!%
!%    ! Mapping weight object
!%    TYPE(MapObj), POINTER     :: map(:,:)
!%
!%    ! Input total in-cloud optical depth (= water OD + ice OD)
!%    REAL*4,       INTENT(IN)  :: fpIn  ( I025x03125, J025x03125, L025x03125 ) 
!%!
!%! !OUTPUT PARAMETERS:
!%!
!%    ! Output total in-cloud optical depth
!%    REAL*4,       INTENT(OUT) :: fpOut ( IMX, JMX, LMX )
!%!
!%! !REVISION HISTORY: 
!%!  12 Jan 2012 - R. Yantosca - Initial version, based on RegridTau
!%!  17 Jan 2012 - R. Yantosca - Now flip levels in the vertical in the 
!%!                              calling routine.  Omit level flipping here.
!%!EOP
!%!------------------------------------------------------------------------------
!%!BOC
!%!
!%! !LOCAL VARIABLES:
!%!
!%    INTEGER           :: I,         J,     L,  LR, nPoints
!%    INTEGER           :: Nx,        Ny,    X,  Y
!%    REAL*4            :: sum_Fn_Wn, sum_Wn   
!%!
!%! !DEFINED PARAMETERS:
!%!
!%    REAL*4, PARAMETER :: maxPrecipFrac = 0.3e0
!%
!%    IF ( IMX == I025x03125 .and. JMX == J025x03125 ) THEN
!%
!%       !====================================================================
!%       ! For the native grid: Assume that if the 0.25 x 0.3125 box has
!%       ! nonzero precipitation, then assume that 100% of the box is 
!%       ! precipitationg.  Therefore, fpOut will either be 0 or 1.
!%       !====================================================================
!%!$OMP PARALLEL DO        &
!%!$OMP DEFAULT( SHARED )  &
!%!$OMP PRIVATE( I, J, L )
!%       DO L = 1, LMX
!%       DO J = 1, JMX
!%       DO I = 1, IMX
!%          IF ( fpIn(I,J,L) > 0e0 ) THEN
!%             fpOut(I,J,L) = maxPrecipFrac   ! Box is precipitating
!%          ELSE
!%             fpOut(I,J,L) = 0e0             ! Box is not precipitating
!%          ENDIF
!%       ENDDO
!%       ENDDO
!%       ENDDO
!%!$OMP END PARALLEL DO
!%
!%    ELSE
!%
!%       !====================================================================
!%       ! For coarser grids: Use the mapping weights (fraction of each fine
!%       ! box that fits into the coarse box) to determine the overall
!%       ! precipitating fraction of the box.  In this case, fpOut will
!%       ! be fractional (varying from 0.000 to 1.000)
!%       !====================================================================
!%
!%       ! Loop over coarse grid boxes
!%!$OMP PARALLEL DO                                                       &
!%!$OMP DEFAULT( SHARED )                                                 &
!%!$OMP PRIVATE( I,  J,  L, nPoints, sum_Fn_Wn, sum_Wn, Nx, Ny, X, Y )
!%       DO L = 1, LMX
!%       DO J = 1, JMX
!%       DO I = 1, IMX
!%
!%          ! Zero output array
!%          fpOut(I,J,L) = 0e0
!%
!%          ! Number of "fine" grid boxes in each dimension
!%          ! that comprise a "coarse" grid box
!%          nPoints = map(I,J)%nPoints
!%             
!%          !---------------------------------
!%          ! Regrid cloud fraction & OD
!%          !---------------------------------
!%          
!%          ! Zero summing variables
!%          sum_Fn_Wn = 0e0
!%          sum_Wn    = 0e0
!%          
!%          ! Loop over "fine" grid boxes
!%          DO Ny = 1, nPoints
!%          DO Nx = 1, nPoints
!%          
!%             ! Avoid useless clock cycles if the mapping weight is zero
!%             IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN
!%
!%                ! Indices of each "fine" grid box that makes up  
!%                ! the "coarse" box
!%                X          = map(I,J)%xInd(Nx)
!%                Y          = map(I,J)%yInd(Ny)
!%
!%                ! Sum of the mapping weights over all of the "fine" grid
!%                ! boxes (X,Y) that make up the "coarse" grid box (I,J)
!%                sum_Wn     = sum_Wn    + map(I,J)%weight(Nx,Ny)
!%
!%                ! If the grid box has nonzero precpitation, then consider
!%                ! that the fraction MAXPRECIPFRAC of the grid box is 
!%                ! precipitating.  Sum the mapping weight of this grid box 
!%                ! into SUM_FN_WN.
!%                IF ( fpIn(X,Y,L) > 0e0 ) THEN
!%                   sum_Fn_Wn = sum_Fn_Wn + ( map(I,J)%weight(Nx,Ny)  &
!%                                         *   maxPrecipFrac          )
!%                ENDIF
!%             ENDIF
!%          ENDDO
!%          ENDDO
!%
!%          ! Divide SUM_FN_WN by SUM_WN to get the overall fraction of the 
!%          ! box that is precipitating.  NOTE, we don't have to worry about 
!%          ! div by zero since SUM( Wn ) will always be greater than zero 
!%          ! (there is always at least 1 "fine" small box in the "coarse" box).
!%          fpOut(I,J,L) = sum_Fn_Wn / sum_Wn
!%
!%       ENDDO
!%       ENDDO
!%       ENDDO
!%!$OMP END PARALLEL DO
!%
!%    ENDIF
!%
!%  END SUBROUTINE ProcessFracPrecip
!%EOC
!%!------------------------------------------------------------------------------
!%!          Harvard University Atmospheric Chemistry Modeling Group            !
!%!------------------------------------------------------------------------------
!%!BOP
!%!
!%! !IROUTINE: IsSafeDiv
!%!
!%! !DESCRIPTION: Function IsSafeDiv returns TRUE if the numerator N and 
!%!  denominator D may be divided safely (i.e. without resulting in a 
!%!  division-by-zero, Not-a-Number (NaN), or Infinity), or FALSE otherwise.
!%!\\
!%!\\
!%! !INTERFACE:
!%!
!%  FUNCTION IsSafeDiv( N, D ) RESULT( isSafe )
!%!
!%! !INPUT PARAMETERS: 
!%!
!%    REAL*4, INTENT(IN) :: N        ! Numerator
!%    REAL*4, INTENT(IN) :: D        ! Denominator
!%!
!%! !RETURN VALUE:
!%!    
!%    LOGICAL            :: isSafe   ! Returns TRUE if it's safe to divide N/D 
!%!
!%! !REVISION HISTORY: 
!%!   29 Sep 2008 - R. Yantosca - Initial Version
!%!EOP
!%!------------------------------------------------------------------------------
!%!BOC
!%
!%    ! Local variables
!%    INTEGER :: MaxExp, MinExp
!%
!%    ! Maxinum 
!%    MaxExp = MAXEXPONENT( N )
!%    MinExp = MINEXPONENT( N )
!%    
!%    ! Test if it's safe to divide 
!%    IF ( ( D                         == 0      )  .or. &
!%         ( EXPONENT(N) - EXPONENT(D) >= MaxExp )  .or. &
!%         ( EXPONENT(N) - EXPONENT(D) <= MinExp ) ) THEN
!%       isSafe = .FALSE.
!%    ELSE
!%       isSafe = .TRUE.
!%    ENDIF
!%
!%  END FUNCTION IsSafeDiv
!EOC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE Geos57A3MstCModule

