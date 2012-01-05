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
PROGRAM CreateTemplateFile
!
! !USES: 
!
  ! Modules for netCDF write
  USE m_netcdf_io_create
  USE m_netcdf_io_define
  USE m_netcdf_io_write
  USE m_netcdf_io_close

  ! Modules for netCDF read
  USE m_netcdf_io_open
  USE m_netcdf_io_close     
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read

  IMPLICIT NONE

  ! Include files
  INCLUDE "netcdf.inc"   ! netCDF include file
!
! !REVISION HISTORY: 
!  05 Jan 2012 - R. Yantosca - Initial version
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Lon & lat
  INTEGER, PARAMETER :: IMX = 1152
  INTEGER, PARAMETER :: JMX = 721
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER            :: fId,      vId,    oMode
  INTEGER            :: idLon,    idLat,  idTime
  CHARACTER(LEN=255) :: fileName, fileId, str
  
  ! Index arrays 
  INTEGER            :: st3d(3),  ct3d(3)
  INTEGER            :: var1(1),  var3(3)
  
  ! Data arrays
  REAL*8             :: lon      (IMX      )
  REAL*8             :: lat      (    JMX  )
  INTEGER            :: time     (        1)
  REAL*4             :: frLake   (IMX,JMX  )
  REAL*4             :: frLand   (IMX,JMX  )
  REAL*4             :: frLandIce(IMX,JMX  )
  REAL*4             :: frOcean  (IMX,JMX  )
  REAL*4             :: lwiMask  (IMX,JMX  )
  
  !=========================================================================
  ! Open the netCDF file
  !=========================================================================
  
  ! Echo info
  WRITE( 6, '(a)' ) '=== Begin readinag from constant data file ==='

  ! File name w/ constant data
  fileName = '/as/cache-old/2011-08/bmy/DAS.ops.asm.const_2d_asm_Nx.GEOS572.00000000_0000.V01.nc4'
  
  ! Open input file name
  CALL Ncop_Rd( fId, TRIM( fileName ) )
  
  ! Read index fields
  CALL NcRd( lon,       fId, 'lon',       (/1/), (/IMX/) )
  CALL NcRd( lat,       fId, 'lat',       (/1/), (/JMX/) )
  CALL NcRd( time,      fId, 'time',      (/1/), (/1  /) )

  ! netCDF indices
  st3d = (/ 1,   1,   1 /)
  ct3d = (/ IMX, JMX, 1 /)
  
  ! Read data fields
  CALL NcRd( frLake,    fId, 'FRLAKE',    st3d,  ct3d    )
  CALL NcRd( frLand,    fId, 'FRLAND',    st3d,  ct3d    )
  CALL NcRd( frLandIce, fId, 'FRLANDICE', st3d,  ct3d    )
  CALL NcRd( frOcean,   fId, 'FROCEAN',   st3d,  ct3d    )
    
  ! Close the file
  CALL NcCl( fId )
  
  WRITE( 6, '(a)' ) '=== End reading from constant data file ==='
  
  !=========================================================================
  ! Create LWI mask: this time, only land vs. ocean
  !=========================================================================
  
  ! Default value
  lwiMask = 1e0
  
  ! Ocean has a value of 0
  WHERE( frOcean > 0 ) lwiMask = 0e0

  ! Land has a value of 1
  WHERE ( frLand + frLandIce + frLake > 0 ) lwiMask = 1e0

  !=========================================================================
  ! Create the netCDF file; define dimensions and global attributes
  !=========================================================================
  
  ! Echo info
  WRITE( 6, '(a)' ) '=== Begin writing to template file ==='

  ! Create file (netCDF-3)
  fileName = 'Geos57TemplateFile.nc'
  CALL NcCr_Wr( fId, TRIM( fileName ) )
  
  ! Turn filling off
  CALL NcSetFill( fId, NF_NOFILL, oMode )
  
  ! Longitude dimension                      
  CALL NcDef_Dimension( fId, 'lon',  IMX, idLon  )
  CALL NcDef_Dimension( fId, 'lat',  JMX, idLat  )
  CALL NcDef_Dimension( fId, 'time', 1,   idTime )

  ! Global attributes
  CALL NcDef_Glob_Attributes( fId, 'Title',       'Template file'   )
  CALL NcDef_Glob_Attributes( fId, 'Filename',     TRIM( fileName ) )
  CALL NcDef_Glob_Attributes( fId, 'Conventions', 'COARDS'          )
  CALL NcDef_Glob_Attributes( fId, 'Format' ,     'netCDF-3'        )
  CALL NcDef_Glob_Attributes( fId, 'Model',       'GEOS-5'          )
  CALL NcDef_Glob_Attributes( fId, 'NLayers',     '1'               )
  CALL NcDef_Glob_Attributes( fId, 'Start_Date',  '20110101'        )
  CALL NcDef_Glob_Attributes( fId, 'Start_Time',  '00:00:00.000000' )
  CALL NcDef_Glob_Attributes( fId, 'End_Date',    '20110101'        )
  CALL NcDef_Glob_Attributes( fId, 'End_Time',    '00:00:00.000000' )
  CALL NcDef_Glob_Attributes( fId, 'Model',       'GEOS5'           )
  CALL NcDef_Glob_Attributes( fId, 'Delta_lon',   '0.3125'          )
  CALL NcDef_Glob_Attributes( fId, 'Delta_lat',   '0.25'            )
  CALL NcDef_Glob_Attributes( fId, 'Delta_Time',  '000000'          )
  
  !=========================================================================
  ! Define the variables and variable attributes
  !=========================================================================

  ! Define longitude variable
  var1 = (/ idLon /)
  vId  = 0
  CALL NcDef_variable( fId, 'lon', NF_DOUBLE, 1, var1, vId )
  CALL NcDef_var_attributes( fId, vId,  'long_name', 'Longitude'   )
  CALL NcDef_var_attributes( fId, vId,  'units',     'degree_east' )
  
  ! Define latitude variable
  var1 = (/ idLat /)
  vId  = vId + 1
  CALL NcDef_variable( fId, 'lat', NF_DOUBLE, 1, var1, vId )
  CALL NcDef_var_attributes( fId, vId, 'long_name', 'Latitude'     )
  CALL NcDef_var_attributes( fId, vId, 'units',     'degree_north' )
  
  ! Define vertical (pressure) variable
  var1 = (/ idTime /)
  vId  = vId + 1
  CALL NcDef_Variable      ( fId, 'time', NF_INT,  1, var1, vId      )
  CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'time'      )
  str  = 'minutes since 2011-01-01 00:00:00'
  CALL NcDef_Var_Attributes( fId, vId, 'units',          TRIM( str ) ) 
  str  = '0000-00-00 00:00:00'
  CALL NcDef_Var_Attributes( fId, vId, 'delta_t',        TRIM( str ) ) 
  CALL NcDef_Var_Attributes( fId, vId, 'begin_date',     '20110101'  )
  CALL NcDef_Var_Attributes( fId, vId, 'begin_date',     '20110101'  )
  CALL NcDef_Var_Attributes( fId, vId, 'time_increment', '000000'    )
  
  ! LWI mask
  var3 = (/ idLon, idLat, idTime /)
  vId  = vId + 1
  CALL NcDef_variable( fId, 'LWI', NF_FLOAT, 3, var3, vId )
  CALL NcDef_var_attributes( fId, vId, 'long_name',      'LWI'       )
  CALL NcDef_var_attributes( fId, vId, 'units',          'unitless'  )
  CALL NcDef_var_attributes( fId, vId, 'gamap_category', 'GMAO-2D'   )
  
  ! FRLAND
  var3 = (/ idLon, idLat, idTime /)
  vId  = vId + 1
  CALL NcDef_variable( fId, 'FRLAND', NF_FLOAT, 3, var3, vId )
  CALL NcDef_var_attributes( fId, vId, 'long_name',      'FRLAND'    )
  CALL NcDef_var_attributes( fId, vId, 'units',          'fraction'  )
  CALL NcDef_var_attributes( fId, vId, 'gamap_category', 'GMAO-2D'   )
  
  ! FRLANDICE
  var3 = (/ idLon, idLat, idTime /)
  vId  = vId + 1
  CALL NcDef_variable( fId, 'FRLANDIC', NF_FLOAT, 3, var3, vId )
  CALL NcDef_var_attributes( fId, vId, 'long_name',      'FRLANDIC'  )
  CALL NcDef_var_attributes( fId, vId, 'units',          'fraction'  )
  CALL NcDef_var_attributes( fId, vId, 'gamap_category', 'GMAO-2D'   )
  
  !=========================================================================
  ! %%% END OF DEFINITION SECTION %%%
  !
  ! Write output to netCDF file
  !=========================================================================
  CALL NcEnd_def( fId )
  
  ! Write index arrays
  CALL NcWr( lon,  fId, 'lon',  (/1/), (/IMX/) )
  CALL NcWr( lat,  fId, 'lat',  (/1/), (/JMX/) )
  CALL NcWr( time, fId, 'time', (/1/), (/1  /) )
  
  ! netCDF indices
  
  ! Write LWI mask to file
  st3d = (/ 1,   1,   1 /)
  ct3d = (/ IMX, JMX, 1 /)
  CALL NcWr( lwiMask,   fId, 'LWI',      st3d, ct3d )
  
  ! Write FRLAND to file
  st3d = (/ 1,   1,   1 /)
  ct3d = (/ IMX, JMX, 1 /)
  CALL NcWr( frLand,    fId, 'FRLAND',   st3d, ct3d )
  
  ! Write FRLANDIC to file
  st3d = (/ 1,   1,   1 /)
  ct3d = (/ IMX, JMX, 1 /)
  CALL NcWr( frLandIce, fId, 'FRLANDIC', st3d, ct3d )
  
  ! Close the fle
  CALL NcCl( fId )

  ! Echo info
  WRITE( 6, '(a)' ) '=== End writing to template file ==='
!EOC
END PROGRAM CreateTemplateFile

