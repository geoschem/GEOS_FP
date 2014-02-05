!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpVOmgModule
!
! !DESCRIPTION: Module GeosFpVOmgModule contains routines to create the 
!  GEOS-Chem average 3-hr data files (winds etc.) from the GEOS-FP raw data.
!\\
!\\
! !INTERFACE: 

MODULE GeosFpA3VOmgModule
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
  PUBLIC  :: GeosFpMakeA3VOmg
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dAsmNv
  PRIVATE :: Process3dCldNv
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  09 Jan 2014 - R. Yantosca - Initial version
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
!  09 Jan 2014 - R. Yantosca - Initial version
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
    lName = 'GEOS-FP time-averaged 3-hour variance of OMEGA (A3vOmg), processed for GEOS-Chem input'
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
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',       TRIM( cal     )  )
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
    vId     =  vId + 1
    var1    = (/ idLon /)
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! VAR_OMEGA
    IF ( StrPos( 'VAR_OMEGA', tavg3_3d_asm_Nv_data ) >= 0 ) THEN
       var4     = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'Variance of OMEGA (Vertical pressure velocity)'
       units = '1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'VAR_OMEGA', NF_FLOAT, 4, var4, vId  )
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
! !IROUTINE: GeosFpMakeA3VOmg
!
! !DESCRIPTION: Routine GeosFpMakeA3VOmg is the the driver routine for 
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (dynamical fields) from 
!       the GEOS-FP raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program GeosFpDriver3.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpMakeA3Dyn()
!
! !REVISION HISTORY: 
!  09 Jan 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_3dCldNv
    INTEGER                 :: nFields_3dAsmNv
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dCldNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dAsmNv(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE GeosFpMakeA3Dyn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )      
       CALL StrRepl     ( fName,     '%%%%%%',    'A3vOmg'               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                    )      
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25,      TIMES_A3,  &
                          xMid_2x25, nc_yMid_2x25, zMid_2x25,  a3Mins,    &
                          gName,     fName,        fOut2x25              )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                ) 
       CALL StrRepl     ( fName,     '%%%%%%',    'A3vOmg'               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                    )      
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5,       TIMES_A3,  &
                          xMid_4x5,  nc_yMid_4x5,  zMid_4x5,   a3Mins,    &
                          gName,     fName   ,     fOut4x5               )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process3dAsmNv()
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A3dyn output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( do2x25   ) CALL NcCl( fOut2x25 )
    IF ( do4x5    ) CALL NcCl( fOut4x5  )

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE GeosFpMakeA3Dyn %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'   
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE GeosFpMakeA3Dyn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dAsmNv
!
! !DESCRIPTION: Subroutine Process3dAsmNv regrids the GEOS-FP met fields 
!  from the "tavg3\_3d\_udt\_Nv" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dAsmNv()
!
! !INPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  09 Jan 2014 - R. Yantosca - Initial versions
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
    INTEGER                 :: X2x25,    Y2x25,   Z2x25,   T2x25
    INTEGER                 :: X4x5,     Y4x5,    Z4x5,    T4x5
    INTEGER                 :: st3d(3),  ct3d(3), st4d(4), ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Q    ( I025x03125, J025x03125, L025x03125 )
    REAL*4                  :: Q2x25( I2x25,      J2x25,      L2x25      )
    REAL*4                  :: Q4x5 ( I4x5,       J4x5,       L4x5       )

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
    msg = '%%%%%% ENTERING ROUTINE Process3dAsmNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

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

       ! Create input filename from the template
       fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_asm_Nv_file )
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
          CALL NcRd( Q, fIn, 'OMEGA', st4d, ct4d )

          ! Strip fill values
          WHERE( Q == FILL_VALUE ) Q = 0e0    

          ! Flip data in vertical
          Qflip => Q( :, :, Z:1:-1 )

          !-----------------------------------------------------------------
          ! Regrid data fields 
          !-----------------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          !%%% Compute variance of OMEGA in each 2 x 2.5 box %%%
          !IF ( do2x25 ) THEN
          !   CALL ComputeVarOmega2x25( Qflip, Q2x25, GRID='2x25' )
          !ENDIF

          !%%% Compute variance of OMEGA in each 4 x 5 box %%%
          IF ( do4x5 ) THEN
             CALL ComputeVarOmega4x5( Qflip, Q4x5,  GRID='4x5'  )
          ENDIF
                
          !-----------------------------------------------------------
          ! Write netCDF output
          !-----------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st4d = (/ 1,     1,     1,     H  /)
             ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, 'VAR_OMEGA', st4d, ct4d )
          ENDIF
       
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st4d  = (/ 1,    1,    1,    H /)
             ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, 'VAR_OMEGA', st4d, ct4d )
          ENDIF
          
          ! Free pointer memory
          NULLIFY( Qflip )
       ENDDO
         
       !--------------------------------------------------------------------
       ! Close input file
       !--------------------------------------------------------------------
       msg = '%%% Closing ' // TRIM( fNameInput )
       WRITE( IU_LOG, '(a)' ) TsRIM( msg )
       CALL NcCl( fIn )       
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% LEAVING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dAsmNv
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ComputeVarOmega2x25( Omega, VarOmega )
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: Omega(I025x03125,J025x03125,L025x03125)
!
! !OUTPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: VarOmega(I2x25,J2x25,L2x25)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  09 Jan 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    

  END SUBROUTINE ComputeVarOmega2x25
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ComputeVarOmega4x5( Omega, VarOmega )
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: Omega(I025x03125,J025x03125,L025x03125)
!
! !OUTPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: VarOmega(I4x5,J4x5,L4x5)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  09 Jan 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER :: X = ( I025x03125 / I4x5 ) * ( J025x03125 / J4x5 )
!
! !LOCAL VARIABLES:
!
    REAL*4             :: DataIn (X)
    REAL*4             :: DataOut(X)
 
    PRINT*, X


  END SUBROUTINE ComputeVarOmega4x5
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: variance
!
! !DESCRIPTION: Computes the variance of a data array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Variance( DATA, N, VAR )
!
! !USES:
!  
    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: N
    REAL*4,  INTENT(IN)  :: DATA(N)
!
! !OUTPUT PARAMETERS:
!
    REAL*4,  INTENT(OUT) :: VAR
!
! !REMARKS:
!  Algorithm from Numerical Recipes.
! 
! !REVISION HISTORY: 
!  09 Jan 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: J
    REAL*4  :: AVE
    REAL*4  :: S

    ! Error check
    IF ( N .le. 1 ) THEN
       PRINT*, 'N must be at least 2'
       STOP
    ENDIF
    
    ! Initialize 
    VAR    = 0d0

    ! Sum and average of the data
    S      = SUM( DATA )
    AVE    = S / DBLE( N )

    ! Compute intermediate quantities
    DO J = 1, N
       S   = DATA(J) - AVE
       VAR = VAR     + ( S * S )
    ENDDO
  
    ! Variance
    VAR    = VAR / DBLE( N - 1 )
  
  END SUBROUTINE Variance
!EOC
END MODULE GeosFpvOmgModule

