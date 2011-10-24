! $Id: TestNcdfUtilities.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: TestNcdfUtilities.F90
!
! !DESCRIPTION: Program TestNcdfUtilities.F90 is the standalone driver that
!  tests if the libNcUtils.a file was built correctly.
!\\
!\\
! !INTERFACE:
!
PROGRAM TestNcdfUtilities
!
! !USES: 
!
  IMPLICIT NONE
!
! !BUGS:  
!  None known at this time
!
! !SEE ALSO: 
!  m_do_err_out.F90
!  m_netcdf_io_checks.F90
!  m_netcdf_io_close.F90
!  m_netcdf_io_create.F90
!  m_netcdf_io_define.F90
!  m_netcdf_io_get_dimlen.F90
!  m_netcdf_io_handle_err.F90
!  m_netcdf_io_open.F90
!  m_netcdf_io_read.F90
!  m_netcdf_io_write.F90
!
! !SYSTEM ROUTINES: 
!  None
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Global private variables
  INTEGER, PARAMETER :: ILONG = 72               ! # of longitude grid points
  INTEGER, PARAMETER :: ILAT  = 46               ! # of latitude  grid points
  INTEGER, PARAMETER :: IVERT = 55               ! # of pressure  levels 
  INTEGER            :: pCt                      ! # of passed tests
  INTEGER            :: tCt                      ! # of total tests
  INTEGER            :: I                        ! Loop index
  INTEGER            :: longdeg, latdeg          ! For longdat, latdat
  REAL*8             :: longdat(ILONG)           ! Longitude data
  REAL*8             :: latdat(ILONG)            ! Latitude data
  REAL*8             :: prsdat(IVERT)            ! Altitude data

!  ! Initialize
!  pCt = 0
!  tCt = 0
!
!  ! Longitude data
!  longdeg = 360.0 / REAL( ILONG )
!  if ( mod( 360, ILONG) /= 0 ) longdeg = longdeg + 1
!  do i = 1, ILONG
!     longdat(i) = i*longdeg
!  enddo
!
!  ! Writing latitude data point
!  latdeg  = 180.0 / REAL( ILAT )
!  if ( mod( 180, ILAT ) /= 0 ) latdeg = latdeg + 1
!  do i = 1, ilong
!     latdat(i) = -90 + (i-0.5)*latdeg
!  enddo
!
!  ! Pressure
!  do i = 1, IVERT
!     prsdat(i) = 1000.00 - (i-1)*(920.00/IVERT) 
!  enddo

  ! Echo info
  WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE( 6, '(a)' ) '%%%  Testing libNcdfUtilities.a  %%%'
  WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

!  ! Create a netCDF file
!  CALL TestNcdfCreate

  ! And try to read it back
  CALL TestNcdfRead
!BOC

CONTAINS

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TestNcdfCreate
!
! !DESCRIPTION: Subroutine TestNcdfCreate creates a netCDF file
!  named \texttt{my\_filename.nc} with the following variables:
!
!  \begin{description}
!  \item[PSF] Surface pressure (2D variable)
!  \item[KEL] Temperature (3D variable)
!  \end{description}
!
!  Fake values are used for the data.  An unlimited dimension is employed
!  to write out several records of kel.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TestNcdfCreate
!
! !USES:
!
    ! Modules 
    USE m_netcdf_io_create
    USE m_netcdf_io_define
    USE m_netcdf_io_write
    USE m_netcdf_io_close
    
    ! Include files
    INCLUDE "netcdf.inc"
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: lon_id,    lat_id,    prs_id
    INTEGER             :: ct1d(1),   ct2d(2),   ct3d(3)
    INTEGER             :: st1d(1),   st2d(2),   st3d(3)
    INTEGER             :: fId, vId
    INTEGER             :: var1(1),   var2(2),   var3(3)
    INTEGER             :: omode, i
    REAL*8              :: psf(ilong, ilat)        ! surface pressure
    REAL*8              :: kel(ilong, ilat, ivert) ! temperature
    
    !=========================================================================
    ! Create the netCDF file
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Begin netCDF file creation test ==='

    CALL NcCr_Wr( fId, 'my_filename.nc' )
    
    !=========================================================================
    ! Define the dimensions
    !=========================================================================
    
    ! Longitude dimension                      
    WRITE( 6, '(a)' ) 'Writing XDim (# lons)   to netCDF file'    
    CALL NcDef_Dimension( fId, 'XDim', ilong, lon_id )
    
    ! Latitude dimension
    WRITE( 6, '(a)' ) 'Writing YDim (# lats)   to netCDF file'    
    CALL NcDef_Dimension( fId, 'YDim', ilat , lat_id )
    
    ! Altitude dimension
    WRITE( 6, '(a)' ) 'Writing ZDim (# alts)   to netCDF file'    
    CALL NcDef_Dimension( fId, 'ZDim', ivert, prs_id )
    
    !=========================================================================
    ! Define the variables and variable attributes
    !=========================================================================

    ! Define longitude variable
    var1 = (/ lon_id /)
    CALL NcDef_variable( fId, 'LON', NF_FLOAT, 1, var1, vId )
    CALL NcDef_var_attributes( fId, vId,  'long_name', 'Longitude'   )
    CALL NcDef_var_attributes( fId, vId,  'units',     'degree_east' )
  
    ! Define latitude variable
    var1 = (/ lat_id /)
    CALL NcDef_variable( fId, 'LAT', NF_FLOAT, 1, var1, vId )
    CALL NcDef_var_attributes( fId, vId, 'long_name', 'Latitude'     )
    CALL NcDef_var_attributes( fId, vId, 'units',     'degree_north' )
    
    ! Define vertical (pressure) variable
    var1 = (/ prs_id /)
    CALL NcDef_variable( fId, 'PLEV', NF_FLOAT, 1, var1, vId )
    CALL NcDef_var_attributes( fId, vId, 'long_name', 'Pressure' )
    CALL NcDef_var_attributes( fId, vId, 'units',     'hPa'      )
    
    ! Define surface pressure variable
    var2 = (/ lon_id, lat_id /)
    CALL NcDef_variable( fId, 'PS', NF_FLOAT, 2, var2, vId )
    CALL NcDef_var_attributes(  fId, vId, 'long_name', 'Surface Pressure' )
    CALL NcDef_var_attributes ( fId, vId, 'units',     'hPa'              )
    
    ! Define 
    var3 = (/lon_id, lat_id, prs_id /)
    CALL NcDef_variable( fId, 'T', NF_FLOAT, 3, var3, vId )
    CALL NcDef_var_attributes( fId, vId, 'long_name', 'Temperature' )
    CALL NcDef_var_attributes( fId, vId, 'units',     'K')
    
    !=========================================================================
    ! Define the variables and variable attributes
    !=========================================================================
    CALL NcDef_Glob_Attributes( fId,'title','Surf. Pres. & Temp. output' )
    CALL NcDef_Glob_Attributes( fId,'history', &
                                    'Initial file for testing - 051005' )
    CALL NcDef_Glob_Attributes( fId,'Conventions','COARDS' )
    
    !=========================================================================
    ! Set fill mode
    !=========================================================================

    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    !=========================================================================
    ! %%% END OF DEFINITION SECTION %%%
    !=========================================================================
    CALL NcEnd_def( fId )
    
    !=========================================================================
    ! Write longitude
    !=========================================================================
    WRITE( 6, '(a)' ) 'Writing LON  (1D array) to netCDF file'    
    st1d = (/ 1     /)
    ct1d = (/ ilong /)
    call NcWr( longdat, fId, 'LON', st1d, ct1d )
    
    !=========================================================================
    ! Write latitude
    !=========================================================================
    WRITE( 6, '(a)' ) 'Writing LAT  (1D array) to netCDF file'    
    st1d = (/ 1    /)
    ct1d = (/ ilat /)
    call NcWr( latdat, fId, 'LAT', st1d, ct1d )
    
    !=========================================================================
    ! Write pressure levels
    !=========================================================================
    WRITE( 6, '(a)' ) 'Writing PLEV (1D array) to netCDF file'    
    st1d = (/ 1     /)
    ct1d = (/ ivert /)
    call NcWr( prsdat, fId, 'PLEV', st1d, ct1d )
    
    !=========================================================================
    ! Write surface pressure
    !=========================================================================
    WRITE( 6, '(a)' ) 'Writing PS   (2D array) to netCDF file'  
    psf  = 1.0
    ct2d = (/ ilong, ilat /)
    st2d = (/ 1,     1    /)
    CALL NcWr( psf, fId, 'PS', st2d, ct2d )
    
    !=========================================================================
    ! Write temperature
    !=========================================================================
    WRITE( 6, '(a)' ) 'Writing T    (3D array) to netCDF file'      
    kel  = 1.0
    ct3d = (/ ilong, ilat, ivert /)
    st3d = (/ 1,     1,    1     /)
    CALL NcWr( kel, fId, 'T', st3d, ct3d )
    
    !=========================================================================
    ! Close the netCDF file
    !=========================================================================

    CALL NcCl( fId )

    ! Echo info
    WRITE( 6, '(a)' ) '=== End netCDF file creation test ==='

  END SUBROUTINE TestNcdfCreate
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TestNcdfRead
!
! !DESCRIPTION: Routine TestNcdfRead extracts the following fields from
!  the netCDF file \texttt{my\_filename.nc}:
!
!  \begin{description}
!  \item[PSF] Surface pressure (2D variable)
!  \item[KEL] Temperature (3D variable).
!  \end{description}
!
!  Note that the file \texttt{my\_filename.nc} was created with fake data
!  values by subroutine TestNcdfCreate.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TestNcdfRead
!
! !USES:
!
    ! Modules
    USE m_netcdf_io_open
    USE m_netcdf_io_close     
    USE m_netcdf_io_get_dimlen
    USE m_netcdf_io_read

    ! Include files
    INCLUDE "netcdf.inc"   ! netCDF include file
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
 
    INTEGER             :: fId,     rc,      I
    INTEGER             :: XDim,    YDim,    ZDim
    INTEGER             :: ct1d(1), ct2d(2), ct3d(3)
    INTEGER             :: st1d(1), st2d(2), st3d(3)
    REAL*8, ALLOCATABLE :: lon(:),  lat(:),  plev(:)
    REAL*8, ALLOCATABLE :: ps(:,:), t(:,:,:)
    CHARACTER(LEN=255)  :: fileName

    !=========================================================================
    ! Open the netCDF file
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Begin netCDF file reading test ==='
    
    filename = '/as/priv/ftp/gmao/gmao/DAS.ops.asm.tavg1_2d_flx_Nx.GEOS572.20110815_0030.V01.nc4'

    CALL Ncop_Rd( fId, TRIM( fileName ) )
    
    !=========================================================================
    ! Get the dimensions
    !=========================================================================
    CALL Ncget_Dimlen( fId, 'lon',  XDim )
    CALL Ncget_Dimlen( fId, 'lat',  YDim )
    CALL Ncget_Dimlen( fId, 'time', ZDim )

    !=========================================================================
    ! Read the T variable
    !=========================================================================
!
    ! Read data    
    ALLOCATE( t( XDim, YDim, ZDim ) )
    st3d = (/ 1,    1,    1    /)
    ct3d = (/ XDim, YDim, ZDim /)
    CALL NcRd( t, fId, 'USTAR', st3d, ct3d )

    print*, 'USTAR(1,1,1): ', T(1,1,1)

    ! Close netCDF file
    CALL NcCl( fId )

    ! Cleanup
    IF ( ALLOCATED( lon  ) ) DEALLOCATE( lon  )
    IF ( ALLOCATED( lat  ) ) DEALLOCATE( lat  )
    IF ( ALLOCATED( plev ) ) DEALLOCATE( plev )
    IF ( ALLOCATED( ps   ) ) DEALLOCATE( ps   )
    IF ( ALLOCATED( t    ) ) DEALLOCATE( t    )

    ! Echo info
    WRITE( 6, '(a)' ) '=== End of netCDF file read test! ==='
    
  END SUBROUTINE TestNcdfRead
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check
!
! !DESCRIPTION: Subroutine that prints "PASSED" or "FAILED" after each test.  
!  Also increments the various counters of passed or failed tests.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check( msg, rc, passCt, totCt )
!
! !INPUT PARAMETERS:
!    
    CHARACTER(LEN=*), INTENT(IN)    :: msg     ! message to print
    INTEGER,          INTENT(IN)    :: rc      ! Return code
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: passCt  ! # of passed tests
    INTEGER,          INTENT(INOUT) :: totCt   ! # of total tests
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC  
!
! !LOCAL VARIABLES:
!
    INTEGER :: s

    ! length of message
    s = LEN( msg )

    IF ( rc == 0 ) THEN
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 45-s ) // 'PASSED'
       passCt = passCt + 1
    ELSE
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 45-s ) // 'FAILED'
    ENDIF

    totCt = totCt + 1

  END SUBROUTINE Check
!EOC

END PROGRAM TestNcdfUtilities

