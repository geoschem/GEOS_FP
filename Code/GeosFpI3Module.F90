!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpI3Module
!
! !DESCRIPTION: Module GeosFpI3Module contains routines to create the 
!  GEOS-Chem instantaneous 3-hr data files from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE: 
!
MODULE GeosFpI3Module
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
  PUBLIC  :: GeosFpMakeI3
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: ProcessI33dAsmNv
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  03 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  04 Jan 2012 - R. Yantosca - Add extra global attributes
!  09 Jan 2012 - R. Yantosca - Now close input file w/in the hourly do loop
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Renamed to GeosFpI3Module
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
!  03 Jan 2012 - R. Yantosca - Initial version, based on GeosFpCnModule
!  04 Jan 2012 - R. Yantosca - Add extra global attributes
!  01 Feb 2012 - R. Yantosca - Make all global attribute names lowercase
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Change and/or add attributes for COARDS
!  23 Sep 2013 - R. Yantosca - Add calendar attribute to time
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
    INTEGER            :: var1(1), var3(3), var4(4)

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
    lName = 'GEOS-FP instantaneous 3-hour (I3) fields, processed for GEOS-Chem input'
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
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lon',            TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lat',            TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'lon',  X,   idLon  )
    CALL NcDef_Dimension( fOut, 'lat',  Y,   idLat  )
    CALL NcDef_Dimension( fOut, 'lev',  Z,   idLev  )
    CALL NcDef_Dimension( fOut, 'time', T,   idTime )

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
    units   = '1'
    CALL NcDef_Variable      ( fOut, 'lev', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    ! Time index array
    var1    = (/ idTime /)
    vId     = vId + 1
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

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! PS
    IF ( StrPos( 'PS', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'Surface pressure' 
       units = 'hPa'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PS',  NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PV (aka EPV)
    IF ( StrPos( 'PV', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Ertel potential vorticity' 
       units = 'K m-2 kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PV', NF_FLOAT, 4, var4, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! QV
    IF ( StrPos( 'QV', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Specific humidity' 
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QV', NF_FLOAT, 4, var4, vId         )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! T
    IF ( StrPos( 'T', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Temperature' 
       units = 'K'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'T', NF_FLOAT, 4, var4, vId    )
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
! !IROUTINE: GeosFpMakeI3
!
! !DESCRIPTION: Routine GeosFpMakeI3 is the the driver routine for 
! \begin{enumerate}
! \item Extracting instantaneous 3-hr data fields (surface values) from 
!       the GEOS-FP raw data files (netCDF4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeI3
!
! !REVISION HISTORY: 
!  03 Jan 2012 - Initial version, based on MERRA
!  11 Jan 2012 - R. Yantosca - Now call StrCompress to remove white space
!                              in the input file name.
!  19 Jan 2012 - R. Yantosca - Now write output to temporary data directories
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file
!  20 Sep 2013 - R. Yantosca - Now save output to nested Europe grid
!  23 Sep 2013 - R. Yantosca - Now define netCDF latitude such that the poles
!                              are at -90/+90.  This facilitates the GIGC
!                              using ESMF/MAPL.
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
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeI3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Return the list of fields and number of fields to process
    ! from each of the MERRA raw met data files
    CALL GetNFields( inst3_3d_asm_Nv_data, nFields, fields )

    ! Total number of fields that we will process
    nAllFields = nFields

    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( inst3_3d_asm_Nv_file ), nFields
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested China output file
    IF ( doNestCh ) THEN
       fname = TRIM( tempDirTmplNestCh ) // TRIM( dataTmplNestCh )
       gName = 'SEAC4RS'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I_NestCh,  J_NestCh,     L025x03125, TIMES_A3,   &
                          xMid_025x03125(I0_ch:I1_ch),                     &
                          yMid_025x03125(J0_ch:J1_ch),                     &
                          zMid_025x03125,                      a3MinsI,    &
                          gName,     fName,        fOutNestCh             )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu ) THEN
       fname = TRIM( tempDirTmplNestEu ) // TRIM( dataTmplNestEu )
       gName = 'nested EU'
       CALL ExpandDate  ( fName,     yyyymmdd,  000000                    )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I_NestEu,  J_NestEu,     L025x03125, TIMES_A3,   &
                          xMid_025x03125(I0_eu:I1_eu),                     &
                          yMid_025x03125(J0_eu:J1_eu),                     &
                          zMid_025x03125,                      a3MinsI,    &
                          gName,     fName,        fOutNestEu             )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa ) THEN
       fname = TRIM( tempDirTmplNestNa ) // TRIM( dataTmplNestNa )
       gName = 'nested NA'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )    
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I_NestNa,  J_NestNa,     L025x03125, TIMES_A3,   &
                          xMid_025x03125(I0_na:I1_na),                     &
                          yMid_025x03125(J0_na:J1_na),                     &
                          zMid_025x03125,                      a3MinsI,    &
                          gName,     fName,        fOutNestNa             )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fname = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25,      TIMES_A3,   &
                          xMid_2x25, nc_yMid_2x25, zMid_2x25,  a3MinsI,    &
                          gName,     fName   ,     fOut2x25               )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fname = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5,       TIMES_A3,   &
                          xMid_4x5,  nc_yMid_4x5,  zMid_4x5,   a3MinsI,    &
                          gName,     fName,        fOut4x5                )
    ENDIF

    ! Regrid fields from the various raw data files
    CALL ProcessI33dAsmNv( nFields, fields )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing I3 output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( doNestCh ) CALL NcCl( fOutNestCh )
    IF ( doNestEu ) CALL NcCl( fOutNestEu )
    IF ( doNestNa ) CALL NcCl( fOutNestNa )
    IF ( do2x25   ) CALL NcCl( fOut2x25   )
    IF ( do4x5    ) CALL NcCl( fOut4x5    )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeI3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeI3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessI33dAsmNv
!
! !DESCRIPTION: Subroutine  ProcessI33dAsmNv regrids the GEOS-FP met fields 
!  from the "inst3t\_3d\_asm\_Nv" file and saves to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ProcessI33dAsmNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REVISION HISTORY: 
!  04 Jan 2012 - R. Yantosca - Initial version, based on MERRA
!  09 Jan 2012 - R. Yantosca - Now close input file w/in the hourly do loop
!  09 Jan 2012 - R. Yantosca - Remove fOut* arguments, they are passed via
!                              the module GeosFpInputsModule.F90
!  10 Jan 2012 - R. Yantosca - Activate parallel loop over vertical levels
!  17 Jan 2012 - R. Yantosca - Bug fix: flip data in vertical immediately
!                              after reading.
!  17 Jan 2012 - R. Yantosca - Nullify pointers after using them  
!  15 Feb 2012 - R. Yantosca - Now save output to nested NA grid netCDF file 
!  20 Sep 2013 - R. Yantosca - Now save output to nested Europe grid
 
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
    INTEGER                 :: X2x25,    Y2x25,    Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,     Z4x5,    T4x5
    INTEGER                 :: st3d(3),  st4d(4)
    INTEGER                 :: ct3d(3),  ct4d(4)

    ! Data arrays (NOTE: 2d or 3d refers to spatial dimensions)
    REAL*4,  TARGET         :: Q2d     ( I025x03125, J025x03125              )
    REAL*4                  :: Q2d_2x25( I2x25,      J2x25                   )
    REAL*4                  :: Q2d_4x5 ( I4x5,       J4x5                    )
    REAL*4,  TARGET         :: Q3d     ( I025x03125, J025x03125, L025x03125  )
    REAL*4                  :: Q3d_2x25( I2x25,      J2x25,      L2x25       )
    REAL*4                  :: Q3d_4x5 ( I4x5,       J4x5,       L4x5        )

    ! Pointers
    REAL*4,  POINTER        :: Ptr_2d(:,:)
    REAL*4,  POINTER        :: Ptr_3d(:,:,:)
    REAL*4,  POINTER        :: Qflip (:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: fNameNested
    CHARACTER(LEN=MAX_CHAR) :: fName2x25
    CHARACTER(LEN=MAX_CHAR) :: fName4x5
    CHARACTER(LEN=MAX_CHAR) :: msg

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
    ! NOTE: For constant file, hardwire date to 2011/01/01
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE ProcessI33dAsmNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( a3minsI(H) / 60 ) * 10000

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( inst3_3d_asm_Nv_file )
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
          
          ! Save field name into an 9-char variable. 
          ! This will truncate field names longer than 8 chars.
          name  = TRIM( fields(F) )
          
          ! Skip if fieldname is empty
          IF ( name == '' ) CYCLE

          ! Zero data arrays
          Q2d      = 0e0
          Q2d_2x25 = 0e0
          Q2d_4x5  = 0e0
          Q3d      = 0e0
          Q3d_2x25 = 0e0
          Q3d_4x5  = 0e0

          ! Test field name
          IF ( TRIM( name ) == 'PS' ) THEN
          
             !==============================================================
             ! Special handling for surface pressure data
             ! since this field is defined at the surface only
             !==============================================================

             ! Start and count index arrays for netCDF
             ! (There is only one data block per file)
             st3d = (/ 1, 1, 1 /)
             ct3d = (/ X, Y, 1 /)

             !--------------------------------------------------------------
             ! Read data
             !--------------------------------------------------------------
             msg = '%%% Reading     ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             CALL NcRd( Q2d, fIn, TRIM( name ), st3d, ct3d )
             
             ! Replace missing values with zeroes
             WHERE( Q2d == FILL_VALUE ) Q2d = 0e0

             ! Convert from [Pa] to [hPa]
             Q2d = Q2d / 100e0

             !--------------------------------------------------------------
             ! Regrid data
             !--------------------------------------------------------------
             msg = '%%% Regridding  ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridGeosFpto2x25( 0, Q2d, Q2d_2x25 )
             ENDIF
             
             ! Regrid to 4x5 
             IF ( do4x5 ) THEN
                CALL RegridGeosFpTo4x5 ( 0, Q2d, Q2d_4x5  )
             ENDIF
 
             !--------------------------------------------------------------
             ! Write netCDF output
             !--------------------------------------------------------------
             msg = '%%% Archiving   ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Nested China (point to proper slice of global data)
             IF ( doNestCh ) THEN
                Ptr_2d => Q2d( I0_ch:I1_ch, J0_ch:J1_ch )
                st3d   = (/ 1,       1,       H /)
                ct3d   = (/ XNestCh, YNestCh, 1 /)
                CALL NcWr( Ptr_2d, fOutNestCh, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

             ! Nested EU (point to proper slice of global data)
             IF ( doNestEu ) THEN
                Ptr_2d => Q2d( I0_eu:I1_eu, J0_eu:J1_eu )
                st3d   = (/ 1,       1,       H /)
                ct3d   = (/ XNestEu, YNestEu, 1 /)
                CALL NcWr( Ptr_2d, fOutNestEu, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

             ! Nested NA (point to proper slice of global data)
             IF ( doNestNa ) THEN
                Ptr_2d => Q2d( I0_na:I1_na, J0_na:J1_na )
                st3d   = (/ 1,       1,       H /)
                ct3d   = (/ XNestNa, YNestNa, 1 /)
                CALL NcWr( Ptr_2d, fOutNestNa, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF
             
             ! Write 2 x 2.5 data
             IF ( do2x25 ) THEN
                st3d  = (/ 1,     1,     H  /)
                ct3d  = (/ X2x25, Y2x25, 1  /)
                CALL NcWr( Q2d_2x25, fOut2x25, TRIM( name ), st3d, ct3d )
             ENDIF
       
             ! Write 4x5 data
             IF ( do4x5 ) THEN
                st3d  = (/ 1,    1,    H /)
                ct3d  = (/ X4x5, Y4x5, 1 /)
                CALL NcWr( Q2d_4x5, fOut4x5, TRIM( name ), st3d, ct3d )
             ENDIF

          ELSE
          
             !==============================================================
             ! Process all other fields
             !==============================================================

             ! Start and count index arrays for netCDF
             ! (There is only one data block per file)
             st4d = (/ 1, 1, 1, 1 /)
             ct4d = (/ X, Y, Z, 1 /)

             ! Special handling: "EPV" is known as "PV" 
             ! for backwards compatibility
             name8 = name
             IF ( TRIM( name8 ) == 'PV' ) name8 = 'EPV'

             !--------------------------------------------------------------
             ! Read data
             !--------------------------------------------------------------
             msg = '%%% Reading     ' // name8
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             CALL NcRd( Q3d, fIn, TRIM( name8 ), st4d, ct4d )
             
             ! Replace missing values with zeroes
             WHERE( Q3d == FILL_VALUE ) Q3d = 0e0

             ! Flip data in the vertical
             Qflip => Q3d( :, :, Z:1:-1 )

             !--------------------------------------------------------------
             ! Regrid data
             !--------------------------------------------------------------
             msg = '%%% Regridding  ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Loop over vertical levels
             DO L = 1, Z 

                ! Regrid to 2 x 2.5
                IF ( do2x25 ) THEN
                   CALL RegridGeosFpto2x25( 0, Qflip(:,:,L), Q3d_2x25(:,:,L) )
                ENDIF
             
                ! Regrid to 4x5 
                IF ( do4x5 ) THEN
                   CALL RegridGeosFpTo4x5 ( 0, Qflip(:,:,L), Q3d_4x5(:,:,L)  )
                ENDIF
             
             ENDDO

             !--------------------------------------------------------------
             ! Write netCDF output
             !--------------------------------------------------------------
             msg = '%%% Archiving   ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Nested China (point to proper slice of global data)
             IF ( doNestCh ) THEN
                Ptr_3d => Qflip( I0_ch:I1_ch, J0_ch:J1_ch, : )
                st4d   = (/ 1,       1,       1,       H /)
                ct4d   = (/ XNestCh, YNestCh, ZNestCh, 1 /)
                CALL NcWr( Ptr_3d, fOutNestCh, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested Europe (point to proper slice of global data)
             IF ( doNestEu ) THEN
                Ptr_3d => Qflip( I0_eu:I1_eu, J0_eu:J1_eu, : )
                st4d   = (/ 1,       1,       1,       H /)
                ct4d   = (/ XNestEu, YNestEu, ZNestEu, 1 /)
                CALL NcWr( Ptr_3d, fOutNestEu, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested North America (point to proper slice of global data)
             IF ( doNestNa ) THEN
                Ptr_3d => Qflip( I0_na:I1_na, J0_na:J1_na, : )
                st4d   = (/ 1,       1,       1,       H /)
                ct4d   = (/ XNestNa, YNestNa, ZNestNa, 1 /)
                CALL NcWr( Ptr_3d, fOutNestNa, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF
             
             ! Write 2 x 2.5 data
             IF ( do2x25 ) THEN
                st4d  = (/ 1,     1,     1,     H  /)
                ct4d  = (/ X2x25, Y2x25, Z2x25, 1  /)
                CALL NcWr( Q3d_2x25, fOut2x25, TRIM( name ), st4d, ct4d )
             ENDIF
       
             ! Write 4x5 data
             IF ( do4x5 ) THEN
                st4d  = (/ 1,    1,    1,    H /)
                ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
                CALL NcWr( Q3d_4x5, fOut4x5, TRIM( name ), st4d, ct4d )
             ENDIF

             ! Free pointer memory
             NULLIFY( Qflip )
          ENDIF
       ENDDO

       !-----------------------------------------------------------------
       ! Close input file
       !-----------------------------------------------------------------
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcCl( fIn )
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================
    msg = '%%%%%% EXITING ROUTINE ProcessI33dAsmNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE ProcessI33dAsmNv
!EOC
END MODULE GeosFpI3Module

