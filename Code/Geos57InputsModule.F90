!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57InputsModule
!
! !DESCRIPTION: Geos57InputsModule contains variables that are used by
!  the various regridding routines.  These are initialized from the
!  \texttt{Geos57Driver.input} file.
!\\
!\\
! !INTERFACE: 

MODULE Geos57InputsModule
! 
! !USES:
!
  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!

  ! Limits
  INTEGER,      PARAMETER :: MAX_FLDS   = 200
  INTEGER,      PARAMETER :: MAX_CHAR   = 1024
                                        
  ! Number of times per file            
  INTEGER,      PARAMETER :: TIMES_A1   = 24          ! # of A1 data times
  INTEGER,      PARAMETER :: TIMES_A3   = 8           ! # of A3 data times
  INTEGER,      PARAMETER :: TIMES_A6   = 4           ! # of A6 data times
  INTEGER,      PARAMETER :: TIMES_I6   = 4           ! # of I6 data times
                                        
  ! File units                          
  INTEGER,      PARAMETER :: IU_LOG     = 6           ! Log file output
  INTEGER,      PARAMETER :: IU_TXT     = 70          ! Text file input
  INTEGER,      PARAMETER :: IU_BIN     = 71
  INTEGER,      PARAMETER :: A1_2x25    = 20          ! 2x25  A3 file
  INTEGER,      PARAMETER :: A1_4x5     = 21          ! 4x5   A3 file 
  INTEGER,      PARAMETER :: A3_2x25    = 31          ! 2x25  A3 file
  INTEGER,      PARAMETER :: A3_4x5     = 32          ! 4x5   A3 file 
  INTEGER,      PARAMETER :: A6_2x25    = 41          ! 2x25  A6 file 
  INTEGER,      PARAMETER :: A6_4x5     = 42          ! 4x5   A6 file
  INTEGER,      PARAMETER :: CN_2x25    = 51          ! 2x25  A6 file 
  INTEGER,      PARAMETER :: CN_4x5     = 52          ! 4x5   A6 file
  INTEGER,      PARAMETER :: I6_2x25    = 61          ! 2x25  I6 file 
  INTEGER,      PARAMETER :: I6_4x5     = 62          ! 4x5   I6 file
                                        
  ! Grid sizes                          
  INTEGER,      PARAMETER :: I05x0666   = 540         ! 0.5  x 0.666 lon dim
  INTEGER,      PARAMETER :: J05x0666   = 361         ! 0.5  x 0.666 lat dim
  INTEGER,      PARAMETER :: L05x0666   = 72          ! 0.5  x 0.666 alt dim
  INTEGER,      PARAMETER :: I125x125   = 288         ! 1.25 x 1.25  lon dim
  INTEGER,      PARAMETER :: J125x125   = 144         ! 1.25 x 1.25  lat dim
  INTEGER,      PARAMETER :: L125x125   = 72          ! 1.25 x 1.25  alt dim
  INTEGER,      PARAMETER :: L125x125_P = 42          ! 1.25 x 1.25  alt dim
  INTEGER,      PARAMETER :: I1x125     = 288         ! 1.0  x 1.25  lon dim
  INTEGER,      PARAMETER :: J1x125     = 181         ! 1.0  x 1.25  lat dim
  INTEGER,      PARAMETER :: L1x125     = 72          ! 1.0  x 1.25  alt dim
  INTEGER,      PARAMETER :: I2x25      = 144         ! 2.0  x 2.5   lon dim
  INTEGER,      PARAMETER :: J2x25      = 91          ! 2.0  x 2.5   lat dim
  INTEGER,      PARAMETER :: L2x25      = 72          ! 2.0  x 2.5   alt dim
  INTEGER,      PARAMETER :: I4x5       = 72          ! 4.0  x 5.0   lon dim
  INTEGER,      PARAMETER :: J4x5       = 46          ! 4.0  x 5.0   lat dim
  INTEGER,      PARAMETER :: L4x5       = 72          ! 4.0  x 5.0   alt dim
!
! !PUBLIC TYPES:
!
  TYPE MapObj
     INTEGER              :: I, J                     ! Lon & lat
     INTEGER              :: nPoints                  ! Size for weight array
     INTEGER,     POINTER :: xInd(:)                  ! Lon indices
     INTEGER,     POINTER :: yInd(:)                  ! Lat indices
     REAL*4,      POINTER :: weight(:,:)              ! Array of mapping wts
  END TYPE MapObj
!
! !PUBLIC DATA MEMBERS:
!
  ! Objects
  TYPE(MapObj),   POINTER :: mapNxTo2x25(:,:)         ! Mapping Nx -> 2 x 2.5
  TYPE(MapObj),   POINTER :: mapNxTo4x5(:,:)          ! Mapping Nx -> 4 x 5
  TYPE(MapObj),   POINTER :: mapFxTo2x25(:,:)         ! Mapping Fx -> 2 x 2.5
  TYPE(MapObj),   POINTER :: mapFxTo4x5(:,:)          ! Mapping Fx -> 4 x 5

  ! Scalars
  LOGICAL                 :: do2x25                   ! Save out 2 x 2.25
  LOGICAL                 :: do4x5                    ! Save out 4 x 5?
  LOGICAL                 :: doMakeCn
  LOGICAL                 :: VERBOSE                  ! Do debug printout?
  INTEGER                 :: yyyymmdd                 ! Today's date
  REAL*4                  :: FILL_VALUE = 1e15        ! Fill value in HDF file
  CHARACTER(LEN=MAX_CHAR) :: dataDirNcdf              ! netCDF data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmpl2x25             ! 2x25  file template
  CHARACTER(LEN=MAX_CHAR) :: dataTmpl4x5              ! 4x5   file template
  CHARACTER(LEN=MAX_CHAR) :: const_2d_asm_Nx_file     ! const_2d_chm_Nx file
  CHARACTER(LEN=MAX_CHAR) :: const_2d_asm_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: inst3_3d_asm_Nv_file     ! inst3_3d_asm_Nv file
  CHARACTER(LEN=MAX_CHAR) :: inst3_3d_asm_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_cld_Nv_file     ! tavg3_3d_cld_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_cld_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Nv_file     ! tavg3_3d_mst_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Ne_file     ! tavg3_3d_mst_Ne file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Ne_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_rad_Nv_file     ! tavg3_3d_rad_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_rad_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_udt_Nv_file     ! tavg3_3d_udt_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_udt_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_lsf_Nv_file     ! tavg3_3d_lsf_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_lsv_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_flx_Nx_file     ! tavg1_2d_flx_Nx 
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_flx_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_lnd_Nx_file     ! tavg1_2d_lnd_Nx 
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_lnd_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_rad_Nx_file     ! tavg1_2d_rad_Nx 
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_rad_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_slv_Nx_file     ! tavg1_2d_slv_Nx 
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_slv_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_2d_ocn_Nx_file     ! tavg3_2d_chm_Fx 
  CHARACTER(LEN=MAX_CHAR) :: tavg3_2d_ocn_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: weightFileNxTo2x25       ! Mapping weights for
  CHARACTER(LEN=MAX_CHAR) :: weightFileNxTo4x5        !  Nx grid and Fx grid
  CHARACTER(LEN=MAX_CHAR) :: weightFileFxTo2x25       !  to GEOS-Chem 2 x 2.5 
  CHARACTER(LEN=MAX_CHAR) :: weightFileFxTo4x5        !  and 4 x 5 grids
  CHARACTER(LEN=MAX_CHAR) :: lwiMaskFile              ! Mask file for LWI
  CHARACTER(LEN=MAX_CHAR) :: frLandIceFile            ! File for FRLANDICE
  CHARACTER(LEN=MAX_CHAR) :: frLandFile               ! File for FRLAND

  ! Arrays
  INTEGER                 :: a1Hours  (TIMES_A1)           ! A3 data times
  INTEGER                 :: a3Hours  (TIMES_A3)           ! A3 data times
  INTEGER                 :: a6Hours  (TIMES_A6)           ! A6 data times
  INTEGER                 :: i6Hours  (TIMES_I6)           ! A6 data times
  REAL*4                  :: lwiMask  (I05x0666,J05x0666)  ! LWI mask
  REAL*4                  :: frLandIce(I05x0666,J05x0666)  ! FRLANDICE data 
  REAL*4                  :: frLand   (I05x0666,J05x0666)  ! FRLAND data
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Geos57Initialize 
  PUBLIC  :: Geos57Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ReadMappingWeights 
!
! !REVISION HISTORY:
!  30 Aug 2011 - R. Yantosca - Initial version, based on MerraInputsModule
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
! !IROUTINE: Geos57Initialize

!
! !DESCRIPTION: This routine deallocates all previously-allocated 
!  module arrays and pointer objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57Initialize
!
! !REVISION HISTORY: 
!  30 Aug 2011 - R. Yantosca - Initial version, based on MerraInputsModule
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    INTEGER                 :: nPts
    CHARACTER(LEN=MAX_CHAR) :: line, temp

    !-----------------------------------------------------------------------
    ! Read the file with the date (passed from the Perl script "doGeos5")
    !-----------------------------------------------------------------------
    
    ! Get day of year
    READ( 5, '(i8)', ERR=990 ) yyyymmdd

    ! Get times of day (these are listed in the GEOS-5 file spec)
    a1Hours = (/ 003000, 013000, 023000, 033000, 043000, 053000, &
                 063000, 073000, 083000, 093000, 103000, 113000, &
                 123000, 133000, 143000, 153000, 163000, 173000, &
                 183000, 193000, 203000, 213000, 223000, 233000 /)

    a3Hours = (/ 013000, 043000, 073000, 103000,  &
                 133000, 163000, 193000, 223000 /)

    a6Hours = (/ 030000, 090000, 150000, 210000 /)

    i6Hours = (/ 000000, 060000, 120000, 180000 /)
    
    !-----------------------------------------------------------------------
    ! Read the file with the filename templates and fields to pull
    !-----------------------------------------------------------------------

    ! Open the file
    OPEN( IU_TXT, FILE='./Geos57Driver.input', STATUS='old' )

    ! Read each line
    DO 

       ! Read a line from the file
       READ( IU_TXT, '(a)', END=100 ) line

       ! Read the various options from the file
       SELECT CASE( TRIM( line ) )

          CASE( '==> Turn on debug print output?' )
             READ( IU_TXT,   *,   ERR=999 ) VERBOSE

          CASE( '==> Local Raw Data Path' )
             READ( IU_TXT, '(a)', ERR=999 ) dataDirNcdf

          CASE( '==> 2 x 2.5 Output' )
             READ( IU_TXT, '(a)', ERR=999 ) dataTmpl2x25
             READ( IU_TXT,   *,   ERR=999 ) do2x25

          CASE( '==> 4 x 5 Output' )
             READ( IU_TXT, '(a)', ERR=999 ) dataTmpl4x5
             READ( IU_TXT,   *,   ERR=999 ) do4x5

          CASE( '==> const_2d_asm_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) const_2d_asm_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) const_2d_asm_Nx_data
             READ( IU_TXT,   *,   ERR=999 ) doMakeCn

          CASE( '==> tavg1_2d_flx_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_flx_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_flx_Nx_data

          CASE( '==> tavg1_2d_lnd_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_lnd_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_lnd_Nx_data

          CASE( '==> tavg1_2d_rad_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_rad_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_rad_Nx_data

          CASE( '==> tavg1_2d_slv_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_slv_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_slv_Nx_data

          CASE( '==> tavg1_2d_ocn_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_ocn_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_ocn_Nx_data

          CASE( '==> tavg1_2d_ocn_Nx' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_ocn_Nx_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg1_2d_ocn_Nx_data

          CASE( '==> tavg3_3d_cld_Nv' ) 
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_cld_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_cld_Nv_data
             
          CASE( '==> tavg3_3d_mst_Nv' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_mst_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_mst_Nv_data

          CASE( '==> tavg3_3d_mst_Ne' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_mst_Ne_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_mst_Ne_data

          CASE( '==> tavg3_3d_rad_Nv' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_rad_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) tavg3_3d_rad_Nv_data

          CASE( '==> tavg3_3d_udt_Nv' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) inst3_3d_udt_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) inst3_3d_udt_Nv_data

          CASE( '==> tavg3_3d_lsf_Nv' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) inst6_3d_lsf_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) inst6_3d_lsf_Nv_data

          CASE( '==> inst3_3d_asm_Nv' )
             READ( IU_TXT, '(a)', ERR=999 ) temp
             READ( IU_TXT, '(a)', ERR=999 ) inst3_3d_asm_Nv_file
             READ( IU_TXT, '(a)', ERR=999 ) inst3_3d_asm_Nv_data

          CASE( '==> Mapping Weight Files' ) 
             READ( IU_TXT, '(a)', ERR=999 ) weightFileNxTo2x25
             READ( IU_TXT, '(a)', ERR=999 ) weightFileNxTo4x5
             READ( IU_TXT, '(a)', ERR=999 ) weightFileFxTo2x25
             READ( IU_TXT, '(a)', ERR=999 ) weightFileFxTo4x5

          CASE( '==> Template Files' ) 
             READ( IU_TXT, '(a)', ERR=999 ) lwiMaskFile
             READ( IU_TXT, '(a)', ERR=999 ) frLandIceFile
             READ( IU_TXT, '(a)', ERR=999 ) frLandFile

          CASE DEFAULT
             ! Nothing

       END SELECT
    ENDDO

    !-----------------------------------------------------------------------
    ! Define the mapping weight objects for cloud regridding 
    ! after reading the file 
    !-----------------------------------------------------------------------
100 CONTINUE

    ! Close file
    CLOSE( IU_TXT )
    
    ! Mapping weights
    IF ( do2x25 ) THEN

       ! Nx grid to 2 x 2.5 grid
       nPts = ( I05x0666 / I2x25 ) + 2
       CALL ReadMappingWeights( weightFileNxTo2x25,              &
                                I2x25, J2x25, nPts, mapNxTo2x25 )

       ! Fx grid to 2 x 2.5 grid
       nPts = ( I1x125 / I2x25 ) + 2
       CALL ReadMappingWeights( weightFileFxTo2x25,                &
                                I2x25, J2x25, nPts, mapFxTo2x25 )
    ENDIF

    IF ( do4x5 ) THEN

       ! Nx grid to 4 x 5 grid
       nPts = ( I05x0666 / I4x5 ) + 2
       CALL ReadMappingWeights( weightFileNxTo4x5,               &
                                I4x5,  J4x5,  nPts, mapNxTo4x5  )

       ! Fx grid to 4 x 5 grid
       nPts = ( I1x125 / I4x5 ) + 2
       CALL ReadMappingWeights( weightFileFxTo4x5,               &
                                I4x5,  J4x5,  nPts, mapFxTo4x5  )
    ENDIF

    !-----------------------------------------------------------------------
    ! Read data from template files
    !-----------------------------------------------------------------------

    ! Default values for LWI regridding
    CALL ReadTemplateFile( lwiMaskFile, lwiMask )
    
    ! FRLANDICE data (for SNOMAS regridding)
    CALL ReadTemplateFile( frLandIceFile, frLandIce )

    ! FRLANDICE data (for SNOMAS regridding)
    CALL ReadTemplateFile( frLandFile, frLand )

    !-----------------------------------------------------------------------
    ! Verbose output for debugging
    !-----------------------------------------------------------------------
    IF ( VERBOSE ) THEN
       PRINT*, 'YYYYMMDD        : ', yyyymmdd
       PRINT*, 'i6Hours         : ', a1Hours
       PRINT*, 'a3Hours         : ', a3Hours
       PRINT*, 'a6Hours         : ', a6Hours
       PRINT*, 'do2x25          : ', do2x25
       PRINT*, 'do4x5           : ', do4x5
       PRINT*, 'doMakeCn        : ', doMakeCn
       PRINT*, 'dataDirHDF      : ', TRIM( dataDirHDF            )
       PRINT*, 'dataFile2x25    : ', TRIM( dataTmpl2x25          )
       PRINT*, 'dataFile4x5     : ', TRIM( dataTmpl4x5           )
       PRINT*, 'const_2d_asm_Nx : ', TRIM( const_2d_asm_Nx_file  )
       PRINT*, '                  ', TRIM( const_2d_asm_Nx_data  )
       PRINT*, 'tavg1_2d_flx_Nx : ', TRIM( tavg1_2d_flx_Nx_file  )
       PRINT*, '                  ', TRIM( tavg1_2d_flx_Nx_data  )
       PRINT*, 'tavg1_2d_lnd_Nx : ', TRIM( tavg1_2d_lnd_Nx_file  )
       PRINT*, '                  ', TRIM( tavg1_2d_lnd_Nx_data  )
       PRINT*, 'tavg1_2d_ocn_Nx : ', TRIM( tavg1_2d_ocn_Nx_file  )
       PRINT*, '                  ', TRIM( tavg1_2d_ocn_Nx_data  )
       PRINT*, 'tavg1_2d_rad_Nx : ', TRIM( tavg1_2d_rad_Nx_file  )
       PRINT*, '                  ', TRIM( tavg1_2d_rad_Nx_data  )
       PRINT*, 'tavg1_2d_slv_Nx : ', TRIM( tavg1_2d_slv_Nx_file  )
       PRINT*, '                  ', TRIM( tavg1_2d_slv_Nx_data  )
       PRINT*, 'tavg3_3d_cld_Nv : ', TRIM( tavg3_2d_cld_Nv_file  )
       PRINT*, '                  ', TRIM( tavg3_2d_chm_Fx_data  )
       PRINT*, 'tavg3_3d_mst_Nv : ', TRIM( tavg3_3d_mst_Nv_file  )
       PRINT*, '                  ', TRIM( tavg3_3d_mst_Nv_data  )
       PRINT*, 'tavg3_3d_mst_Ne : ', TRIM( tavg3_3d_mst_Ne_file  )
       PRINT*, '                  ', TRIM( tavg3_3d_mst_Ne_data  )
       PRINT*, 'tavg3_3d_rad_Nv : ', TRIM( tavg3_3d_rad_Nv_file  )
       PRINT*, '                  ', TRIM( tavg3_3d_rad_Nv_data  )
       PRINT*, 'tavg3_3d_udt_Nv : ', TRIM( tavg3_3d_udt_Nv_file  )
       PRINT*, '                  ', TRIM( tavg3_3d_udt_Cp_data  )
       PRINT*, 'tavg3_3d_lsv_Nv : ', TRIM( tavg3_3d_lsf_Nv_file  )
       PRINT*, '                  ', TRIM( tavg3_3d_lsf_Cp_data  )
       PRINT*, 'inst3_3d_asm_Nv : ', TRIM( inst3_3d_asm_Nv_file  )
       PRINT*, '                  ', TRIM( inst3_3d_asm_Nv_data  )
       PRINT*, 'WeightsNxTo2x25 : ', TRIM( weightFileNxTo2x25    )
       PRINT*, 'WeightsNxTo4x5  : ', TRIM( weightFileNxTo4x5     )
       PRINT*, 'WeightsFxTo2x25 : ', TRIM( weightFileFxTo2x25    )
       PRINT*, 'WeightsFxTo4x5  : ', TRIM( weightFileFxTo4x5     )
       PRINT*, 'lwiMaskFile     : ', TRIM( lwiMaskFile           )
       PRINT*, 'frLandIceFile   : ', TRIM( frLandIceFile         )
       PRINT*, 'frLandFile      : ', TRIM( frLandFile            )
    ENDIF

    ! Write a message to denote if we are using pressure-weighting
    ! when regridding the U, V, U10M, V10M winds (bmy, 10/11/08)
    WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    WRITE( 6, '(a)' ) '%%%  Regridding U, V winds weighted by pressure  %%%'
    WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    ! Return to calling program
    RETURN

    !-----------------------------------------------------------------------
    ! The date file had an I/O error; exit w/ error message
    !-----------------------------------------------------------------------
990 CONTINUE
    PRINT*, '%%% Error reading DATE!'
    CALL EXIT(1)

    !-----------------------------------------------------------------------
    ! The "filename.input" file had an I/O error; exit w/ error message
    !-----------------------------------------------------------------------
999 CONTINUE
    PRINT*, '%%% Error reading "filename.input"!'
    CALL EXIT(1)

  END SUBROUTINE Geos57Initialize
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadMappingWeights
!
! !DESCRIPTION: This routine reads the mapping weights from the GEOS-5 
!  0.5 x 0.667 grid to coarser resolution grids (e.g. 2 x 2.5, 4 x 5).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadMappingWeights( fileName, IMX, JMX, nPts, map )
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN) :: fileName   ! Name of file w/ weight info
    INTEGER,          INTENT(IN) :: IMX, JMX   ! Lon & lat dims of coarse grid
    INTEGER,          INTENT(IN) :: nPts       ! # of points to read in
!
! !OUTPUT PARAMETERS:
!
    TYPE(MapObj),     POINTER    :: map(:,:)
!
! !REVISION HISTORY: 
!   25 Sep 2008 - R. Yantosca - Initial Version
!
! !REMARKS:
!   If the MAP object is not defined, ReadMappingWeights will allocate
!   it and initialize it here.  The user is responsible for deallocating
!   it elsewhere.
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    INTEGER                      :: I, J, Nx, Ny, rc
    CHARACTER(LEN=11)            :: fmtStr

    !========================================================================
    ! Initialize the pointer object if necessary
    !========================================================================
    IF ( .not. ASSOCIATED ( map ) ) THEN

       ! Allocate space for the MAP variable
       ALLOCATE( map( IMX, JMX ), STAT=rc )
  
       ! Loop over lats & lons
       DO J = 1, JMX
       DO I = 1, IMX

          ! Allocate pointer fields of MAP
          ALLOCATE( map(I,J)%xInd( nPts ),         STAT=rc )
          ALLOCATE( map(I,J)%yInd( nPts ),         STAT=rc )
          ALLOCATE( map(I,J)%weight( nPts, nPts ), STAT=rc )

          ! Initialize fields
          map(I,J)%I       = I
          map(I,J)%I       = J
          map(I,J)%nPoints = nPts
          map(I,J)%xInd    = 0
          map(I,J)%yInd    = 0
          map(I,J)%weight  = 0e0
       ENDDO
       ENDDO
    ENDIF

    !========================================================================
    ! Read data
    !========================================================================


    ! Pick the format string for the various resolutions (as created by
    ! IDL program ctm_getweight.pro).  The format string length varies so that
    ! we don't introduce any extraneous values at 
    ! introdu
    IF ( IMX == 72 .and. JMX == 46 ) THEN
       fmtStr = '(3x,12f7.3)'                           ! 4 x 5
    ELSE IF ( IMX == 144 .and. JMX == 91 ) THEN
       fmtStr = '(3x,12f8.4)'                           ! 2 x 2.5
    ELSE IF ( IMX == 288 .and. JMX == 181 ) THEN
       fmtStr = '(3x,12f9.5)'                           ! 1 x 1.25
    ELSE IF ( IMX == 360 .and. JMX == 181 ) THEN
       fmtStr = '(3x,12f8.4)'                           ! 1 x 1           
    ENDIF

    ! Open file with weights
    OPEN( 10, FILE=TRIM( fileName ), IOSTAT=rc )
    IF ( rc /= 0 ) THEN 
       WRITE( 6, '(a)' ) 'Cannot open ' // TRIM( fileName )
    ENDIF

    ! Read data
    DO 

       ! Read "coarse" grid box indices to file
       READ( 10, '(2i4)', IOSTAT=rc ) I, J

       ! Test for end-of-file
       IF ( rc < 0 ) EXIT

       ! Convert from IDL notation to F90 notation
       I = I + 1
       J = J + 1
     
       ! Read lon & lat indices of "fine" boxes that comprise a "coarse" box
       READ( 10, '(3x,12i4)'   ) ( map(I,J)%xInd(Nx), Nx=1,nPts )
       READ( 10, '(3x,12i4)'   ) ( map(I,J)%yInd(Ny), Ny=1,nPts )

       ! Convert from IDL notation to F90 notation
       map(I,J)%xInd(:) =  map(I,J)%xInd(:) + 1
       map(I,J)%yInd(:) =  map(I,J)%yInd(:) + 1

       ! Read mapping weights (fraction of each "fine" box that 
       ! lies within each "coarse" box
       READ( 10, fmtStr ) (( map(I,J)%weight(Nx,Ny), Nx=1,nPts ), Ny=1,nPts )

    ENDDO

    ! Close file
    CLOSE( 10 )

  END SUBROUTINE ReadMappingWeights
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadTemplateFile
!
! !DESCRIPTION: This routine deallocates all previously-allocated 
!  module arrays and pointer objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadTemplateFile( fileName, dataArray )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: fileName
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: dataArray(:,:)
!
! !REVISION HISTORY: 
!  25 Aug 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IOS

    ! Open the LWI mask file
    OPEN( IU_BIN,     FILE=TRIM( fileName ), STATUS='OLD',        &
          IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

    ! Exit w/ error msg if 
    IF ( IOS /= 0 ) THEN 
       PRINT*, '### ERROR OPENING ', TRIM( fileName )
       CALL EXIT(1)
    ENDIF

    ! Read the data
    READ( IU_BIN, IOSTAT=IOS ) dataArray

    ! Exit w/ error msg if failure
    IF ( IOS > 0 ) THEN 
       PRINT*, '### ERROR READING DATA IN ', TRIM( fileName )
       CALL EXIT(1)
    ENDIF

    ! Close the file
    CLOSE( IU_BIN )

  END SUBROUTINE ReadTemplateFile
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Geos57Cleanup
!
! !DESCRIPTION: This routine deallocates all previously-allocated 
!  module arrays and pointer objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57Cleanup
!
! !REVISION HISTORY: 
!  25 Sep 2008 - R. Yantosca - Initial Version
!  17 Aug 2010 - R. Yantosca - Now deallocate mapNxTo2x25, mapNxTo4x5,
!                              mapFxTo2x25, mapFxTo4x5 objects
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: I, J

    !======================================================================
    ! Deallocate 2x25 mapping weight objects
    !======================================================================
    IF ( do2x25 ) THEN

       ! Echo info
       IF ( VERBOSE ) WRITE( 6, 100 ) 
100    FORMAT( 'Deallocating mapping weight objects for 2 x 2.5 grid' )

       ! Loop over 2 x 2.5 boxes
       DO J = 1, J2x25
       DO I = 1, I2x25

          !-------------------------------------------------
          ! Deallocate Nx grid to 2 x 2.5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapNxTo2x25(I,J)%xInd  ) ) THEN 
             DEALLOCATE( mapNxTo2x25(I,J)%xInd )
          ENDIF
 
          IF ( ASSOCIATED( mapNxTo2x25(I,J)%yInd  ) ) THEN 
             DEALLOCATE( mapNxTo2x25(I,J)%yInd )
          ENDIF
          
          IF ( ASSOCIATED( mapNxTo2x25(I,J)%weight) ) THEN
             DEALLOCATE( mapNxTo2x25(I,J)%weight )
          ENDIF

          !-------------------------------------------------
          ! Deallocate Fx grid to 2 x 2.5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapFxTo2x25(I,J)%xInd  ) ) THEN 
             DEALLOCATE( mapFxTo2x25(I,J)%xInd  )
          ENDIF

          IF ( ASSOCIATED( mapFxTo2x25(I,J)%yInd  ) ) THEN 
             DEALLOCATE( mapFxTo2x25(I,J)%yInd  )
          ENDIF
          
          IF ( ASSOCIATED( mapFxTo2x25(I,J)%weight) ) THEN
             DEALLOCATE( mapFxTo2x25(I,J)%weight)
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapNxTo2x25 ) ) DEALLOCATE( mapNxTo2x25 )
       IF ( ASSOCIATED( mapFxTo2x25 ) ) DEALLOCATE( mapFxTo2x25 )
    ENDIF

    !======================================================================
    ! Deallocate 4x5 mapping weight objects
    !======================================================================
    IF ( do4x5 ) THEN

       IF ( VERBOSE ) WRITE( 6, 110 ) 
110    FORMAT( 'Deallocating mapping weight objects for 4 x 5 grid' )

       ! Loop over 4 x 5 boxes
       DO J = 1, J4x5
       DO I = 1, I4x5

          !-------------------------------------------------
          ! Deallocate Nx grid to 4x5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapNxTo4x5(I,J)%xInd  ) ) THEN 
             DEALLOCATE( mapNxTo4x5(I,J)%xInd )
          ENDIF
 
          IF ( ASSOCIATED( mapNxTo4x5(I,J)%yInd  ) ) THEN 
             DEALLOCATE( mapNxTo4x5(I,J)%yInd )
          ENDIF
          
          IF ( ASSOCIATED( mapNxTo4x5(I,J)%weight) ) THEN
             DEALLOCATE( mapNxTo4x5(I,J)%weight )
          ENDIF

          !-------------------------------------------------
          ! Deallocate Fx grid to 2 x 2.5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapFxTo4x5(I,J)%xInd  ) ) THEN 
             DEALLOCATE( mapFxTo4x5(I,J)%xInd  )
          ENDIF

          IF ( ASSOCIATED( mapFxTo4x5(I,J)%yInd  ) ) THEN 
             DEALLOCATE( mapFxTo4x5(I,J)%yInd  )
          ENDIF
          
          IF ( ASSOCIATED( mapFxTo4x5(I,J)%weight) ) THEN
             DEALLOCATE( mapFxTo4x5(I,J)%weight)
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapNxTo4x5 ) ) DEALLOCATE( mapNxTo4x5 )
       IF ( ASSOCIATED( mapFxTo4x5 ) ) DEALLOCATE( mapFxTo4x5 )
    ENDIF

  END SUBROUTINE Geos57Cleanup

END MODULE Geos57InputsModule
!EOC
