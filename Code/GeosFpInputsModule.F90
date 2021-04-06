!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GeosFpInputsModule
!
! !DESCRIPTION: GeosFpInputsModule contains variables that are used by
!  the various regridding routines.  These are initialized from the
!  \texttt{GeosFpDriver.input} file.
!\\
!\\
! !INTERFACE:

MODULE GeosFpInputsModule
!
! !USES:
!
  ! Modules for reading netCDF
  USE m_netcdf_io_open
  USE m_netcdf_io_close
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read

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

  ! File units
  INTEGER,      PARAMETER :: IU_LOG     = 6           ! Log file output
  INTEGER,      PARAMETER :: IU_TXT     = 70          ! Text file input
  INTEGER,      PARAMETER :: IU_BIN     = 71

  ! Grid size dimensions
  INTEGER,      PARAMETER :: I025x03125 = 1152        ! 0.25 x 0.3125 lon dim
  INTEGER,      PARAMETER :: J025x03125 = 721         ! 0.25 x 0.3125 lat dim
  INTEGER,      PARAMETER :: L025x03125 = 72          ! 0.25 x 0.3125 alt dim
  INTEGER,      PARAMETER :: I05x0666   = 540         ! 0.5  x 0.666  lon dim
  INTEGER,      PARAMETER :: J05x0666   = 361         ! 0.5  x 0.666  lat dim
  INTEGER,      PARAMETER :: L05x0666   = 72          ! 0.5  x 0.666  alt dim
  INTEGER,      PARAMETER :: I125x125   = 288         ! 1.25 x 1.25   lon dim
  INTEGER,      PARAMETER :: J125x125   = 144         ! 1.25 x 1.25   lat dim
  INTEGER,      PARAMETER :: L125x125   = 72          ! 1.25 x 1.25   alt dim
  INTEGER,      PARAMETER :: I1x125     = 288         ! 1.0  x 1.25   lon dim
  INTEGER,      PARAMETER :: J1x125     = 181         ! 1.0  x 1.25   lat dim
  INTEGER,      PARAMETER :: L1x125     = 72          ! 1.0  x 1.25   alt dim
  INTEGER,      PARAMETER :: I2x25      = 144         ! 2.0  x 2.5    lon dim
  INTEGER,      PARAMETER :: J2x25      = 91          ! 2.0  x 2.5    lat dim
  INTEGER,      PARAMETER :: L2x25      = 72          ! 2.0  x 2.5    alt dim
  INTEGER,      PARAMETER :: I4x5       = 72          ! 4.0  x 5.0    lon dim
  INTEGER,      PARAMETER :: J4x5       = 46          ! 4.0  x 5.0    lat dim
  INTEGER,      PARAMETER :: L4x5       = 72          ! 4.0  x 5.0    alt dim
! -----(lzh,06/20/2014)-----
  INTEGER,      PARAMETER :: I05x0625   = 576         ! 0.5  x 0.625  lon dim
  INTEGER,      PARAMETER :: J05x0625   = 361         ! 0.5  x 0.625  lat dim
  INTEGER,      PARAMETER :: L05x0625   = 72          ! 0.5  x 0.625  alt dim
!-----(finish edit)---------

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
  TYPE(MapObj),   POINTER :: mapNative(:,:)           ! Placeholder map object
  TYPE(MapObj),   POINTER :: mapTo2x25(:,:)           ! Map native -> 2 x 2.5
  TYPE(MapObj),   POINTER :: mapTo4x5(:,:)            ! Map native -> 4 x 5
! (lzh, 06/20/2014)
  TYPE(MapObj),   POINTER :: mapTo05x0625(:,:)        ! Map native -> 0.5 x 0.625

  ! NetCDF file Handles

  ! Scalars
  LOGICAL                 :: doNative                 ! Process native grid?
  LOGICAL                 :: doNestCh                 ! Save nested CH grid?
  INTEGER                 :: I0_ch,    J0_ch          ! LL corner of CH grid
  INTEGER                 :: I1_ch,    J1_ch          ! UR corner of CH grid
  INTEGER                 :: I_NestCh, J_NestCh       ! NestCh dimensions
  LOGICAL                 :: doNestEu                 ! Save nested EU grid?
  INTEGER                 :: I0_eu,    J0_eu          ! LL corner of EU grid
  INTEGER                 :: I1_eu,    J1_eu          ! UR corner of EU grid
  INTEGER                 :: I_NestEu, J_NestEu       ! NestNa dimensions
  LOGICAL                 :: doNestNa                 ! Save nested NA grid?
  INTEGER                 :: I0_na,    J0_na          ! LL corner of NA grid
  INTEGER                 :: I1_na,    J1_na          ! UR corner of NA grid
  INTEGER                 :: I_NestNa, J_NestNa       ! NestNa dimensions
  LOGICAL                 :: doNestSe                 ! Save nested SE grid?
  INTEGER                 :: I0_se,    J0_se          ! LL corner of SE grid
  INTEGER                 :: I1_se,    J1_se          ! UR corner of SE grid
  INTEGER                 :: I_NestSe, J_NestSe       ! NestSe dimensions
! (jxu, 2015/12/08, add nested Asia)
  LOGICAL                 :: doNestAs                 ! Save nested Asia grid?
  INTEGER                 :: I0_as,    J0_as          ! LL corner of Asia grid
  INTEGER                 :: I1_as,    J1_as          ! UR corner of Dia grid
  INTEGER                 :: I_NestAs, J_NestAs       ! NestAs dimensions
! (jxu, end)
  !(jxu, 2016/02/13, add 0.25 global)
  LOGICAL                 :: do025x03125              ! Save out 0.25 x 0.3125
  !(jxu, end)
  LOGICAL                 :: do2x25                   ! Save out 2 x 2.25
  LOGICAL                 :: do4x5                    ! Save out 4 x 5?
  LOGICAL                 :: doMakeCn
  LOGICAL                 :: VERBOSE                  ! Do debug printout?
  INTEGER                 :: yyyymmdd                 ! Today's date
  INTEGER                 :: fIn                      ! NC fId; input
  INTEGER                 :: fOutNestCh               ! NC fId; output CH grid
  INTEGER                 :: fOutNestEu               ! NC fId; output EU grid
  INTEGER                 :: fOutNestNa               ! NC fId; output NA grid
  INTEGER                 :: fOutNestSe               ! NC fId; output SE grid
! (jxu, 2015/12/08, add nested Asia)
  INTEGER                 :: fOutNestAs               ! NC fId; output IN grid
! (jxu, end)
!(jxu, 2016/02/13, add 0.25 global)
 INTEGER                 :: fOut025x03125             ! NC fId; output 025x03125
!(jxu, end)
  INTEGER                 :: fOut2x25                 ! NC fId; output 2x25
  INTEGER                 :: fOut4x5                  ! NC fId; output 4x5
  REAL*4                  :: FILL_VALUE = 1e15        ! Fill value in HDF file
  REAL*4                  :: Ap(L025x03125+1)         ! Hybrid grid "A" array
  REAL*4                  :: Bp(L025x03125+1)         ! Hybrid grid "B" array
  CHARACTER(LEN=8)        :: yyyymmdd_string          ! String for YYYYMMDD
  CHARACTER(LEN=MAX_CHAR) :: inputDataDir             ! netCDF data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestCh           ! NstCh file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestCh        ! NstCh temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestCh        ! NstCh data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestNa           ! NstNa file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestNa        ! NstNa temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestNa        ! NstNa data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestEu           ! NstEu file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestEu        ! NstEu temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestEu        ! NstEu data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestSe           ! NstSe file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestSe        ! NstSe temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestSe        ! NstSe data dir
! (jxu, 2015/12/08, add nested Asia)
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestAs           ! NstIn file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestAs        ! NstIn temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestAs        ! NstIn data dir
! (jxu, end)
!(jxu, 2016/02/13, add 0.25 global)
 CHARACTER(LEN=MAX_CHAR) :: dataTmpl025x03125         ! 025x03125  file template
 CHARACTER(LEN=MAX_CHAR) :: tempDirTmpl025x03125      ! 025x03125  temp dir
 CHARACTER(LEN=MAX_CHAR) :: dataDirTmpl025x03125      ! 025x03125  data dir
!(jxu, end)
  CHARACTER(LEN=MAX_CHAR) :: dataTmpl2x25             ! 2x25  file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmpl2x25          ! 2x25  temp dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmpl2x25          ! 2x25  data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmpl4x5              ! 4x5   file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmpl4x5           ! 4x5   temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmpl4x5           ! 4x5   data dir
  CHARACTER(LEN=MAX_CHAR) :: const_2d_asm_Nx_file     ! const_2d_chm_Nx file
  CHARACTER(LEN=MAX_CHAR) :: const_2d_asm_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: inst3_3d_asm_Nv_file     ! inst3_3d_asm_Nv file
  CHARACTER(LEN=MAX_CHAR) :: inst3_3d_asm_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_asm_Nv_file     ! tavg3_3d_udt_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_asm_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_cld_Nv_file     ! tavg3_3d_cld_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_cld_Nv_data_c   !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_cld_Nv_data_d   !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Nv_file     ! tavg3_3d_mst_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Ne_file     ! tavg3_3d_mst_Ne file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_mst_Ne_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_rad_Nv_file     ! tavg3_3d_rad_Nv file
  CHARACTER(LEN=MAX_CHAR) :: tavg3_3d_rad_Nv_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_flx_Nx_file     ! tavg1_2d_flx_Nx
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_flx_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_lnd_Nx_file     ! tavg1_2d_lnd_Nx
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_lnd_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_rad_Nx_file     ! tavg1_2d_rad_Nx
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_rad_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_slv_Nx_file     ! tavg1_2d_slv_Nx
  CHARACTER(LEN=MAX_CHAR) :: tavg1_2d_slv_Nx_data     !  and list of data flds
  CHARACTER(LEN=MAX_CHAR) :: weightFileTo2x25       ! Mapping weights for
  CHARACTER(LEN=MAX_CHAR) :: weightFileTo4x5        !  Nx grid and Fx grid
  CHARACTER(LEN=MAX_CHAR) :: templateFile             ! Mask file for LWI

  ! Arrays
  INTEGER                 :: a1Mins   (TIMES_A1)               ! A1 data times
  INTEGER                 :: a3MinsI  (TIMES_A3)               ! Inst A1 times
  INTEGER                 :: a3Mins   (TIMES_A3)               ! A3 data times
  REAL*4                  :: lwiMask  (I025x03125,J025x03125)  ! LWI mask
  REAL*4                  :: frLandIce(I025x03125,J025x03125)  ! FRLANDICE data
  REAL*4                  :: frLand   (I025x03125,J025x03125)  ! FRLAND data

! ----- (lzh,06/20/2014) -----------
  LOGICAL                 :: doNestCh05               ! Save nested CH grid?
  INTEGER                 :: I0_ch05,    J0_ch05      ! LL corner of CH grid
  INTEGER                 :: I1_ch05,    J1_ch05      ! UR corner of CH grid
  INTEGER                 :: I_NestCh05, J_NestCh05   ! NestCh dimensions
  LOGICAL                 :: doNestEu05               ! Save nested EU grid?
  INTEGER                 :: I0_eu05,    J0_eu05      ! LL corner of EU grid
  INTEGER                 :: I1_eu05,    J1_eu05      ! UR corner of EU grid
  INTEGER                 :: I_NestEu05, J_NestEu05   ! NestNa dimensions
  LOGICAL                 :: doNestNa05               ! Save nested NA grid?
  INTEGER                 :: I0_na05,    J0_na05      ! LL corner of NA grid
  INTEGER                 :: I1_na05,    J1_na05      ! UR corner of NA grid
  INTEGER                 :: I_NestNa05, J_NestNa05   ! NestNa dimensions
  LOGICAL                 :: doNestSe05               ! Save nested SE grid?
  INTEGER                 :: I0_se05,    J0_se05      ! LL corner of SE grid
  INTEGER                 :: I1_se05,    J1_se05      ! UR corner of SE grid
  INTEGER                 :: I_NestSe05, J_NestSe05   ! NestSe dimensions
!(jxu, 2015/12/08, add nested Asia)
  LOGICAL                 :: doNestAs05               ! Save nested Asia grid?
  INTEGER                 :: I0_as05,    J0_as05      ! LL corner of IN grid
  INTEGER                 :: I1_as05,    J1_as05      ! UR corner of IN grid
  INTEGER                 :: I_NestAs05, J_NestAs05   ! NestAs dimensions
!(jxu, end)
  LOGICAL                 :: do05x0625                ! Save out 0.5 x 0.625

  INTEGER                 :: fOut05NestCh             ! NC fId; output CH grid
  INTEGER                 :: fOut05NestEu             ! NC fId; output EU grid
  INTEGER                 :: fOut05NestNa             ! NC fId; output NA grid
  INTEGER                 :: fOut05NestSe             ! NC fId; output SE grid
!(jxu, 2015/12/08, add nested Asia)
  INTEGER                 :: fOut05NestAs             ! NC fId; output In grid
!(jxu, end)
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestCh05           ! NstCh file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestCh05        ! NstCh temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestCh05        ! NstCh data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestNa05           ! NstNa file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestNa05        ! NstNa temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestNa05        ! NstNa data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestEu05           ! NstEu file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestEu05        ! NstEu temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestEu05        ! NstEu data dir
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestSe05           ! NstSe file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestSe05        ! NstSe temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestSe05        ! NstSe data dir
!(jxu, 2015/12/08, add nested Asia)
  CHARACTER(LEN=MAX_CHAR) :: dataTmplNestAs05           ! NstIn file template
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplNestAs05        ! NstIn temporary dir
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplNestAs05        ! NstIn data dir
!(lb, 2021/03/22, added global 05)  
  LOGICAL                 :: doGlobal05 
  CHARACTER(LEN=MAX_CHAR) :: dataTmplGlobal05
  CHARACTER(LEN=MAX_CHAR) :: tempDirTmplGlobal05
  CHARACTER(LEN=MAX_CHAR) :: dataDirTmplGlobal05
  INTEGER                 :: fOutGlobal05
!(jxu, end)
  CHARACTER(LEN=MAX_CHAR) :: weightFileTo05x0625      ! Mapping weights
! ----- (finish edit) -----------
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GeosFpInitialize
  PUBLIC  :: GeosFpCleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ReadMappingWeights
!
! !REVISION HISTORY:
!  30 Aug 2011 - R. Yantosca - Initial version, based on MerraInputsModule
!  21 Dec 2011 - R. Yantosca - Now add a3Mins, a3MinsI, a1Mins variables
!  03 Jan 2012 - R. Yantosca - Add Ap, Bp arrays for hybrid grid definition
!  05 Jan 2012 - R. Yantosca - ReadTemplateFile now reads netCDF data
!  06 Jan 2012 - R. Yantosca - Define "doNative" logical as a convenience
!                              variable to denote when to read native data
!  11 Jan 2012 - R. Yantosca - Now split fields from tavg3_3d_cld_Nv and
!                              tavg3_3d_mst_Ne into output files
!  20 Jan 2012 - R. Yantosca - Now use lowercase "output" string for all grids
!  15 Feb 2012 - R. Yantosca - Add variables for nested NA grid
!  19 Sep 2013 - R. Yantosca - Renamed to GeosFpInputsModule
!  19 Sep 2013 - R. Yantosca - Added variables for nested Europe grid  (EU)
!  19 Sep 2013 - R. Yantosca - Remove reference to tavg3_3d_rad_Nv collection
!  19 Sep 2013 - R. Yantosca - Remove variables for collections that we don't
!                              use any longer
!  08 Oct 2013 - R. Yantosca - Added variables for nested SE Asia Grid (SE)
!  05 Dec 2013 - R. Yantosca - Add back the tavg3_3d_rad_Nv variables
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
! !IROUTINE: GeosFpInitialize

!
! !DESCRIPTION: This routine deallocates all previously-allocated
!  module arrays and pointer objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpInitialize
!
! !REVISION HISTORY:
!  30 Aug 2011 - R. Yantosca - Initial version, based on MerraInputsModule
!  21 Dec 2011 - R. Yantosca - Now also initialize a3Mins, a3MinsI, a1Mins
!  11 Jan 2012 - R. Yantosca - Split fields from the tavg3_3d_cld_Nv and
!                              tavg3_3d_Mst_Ne collections into multiple
!                              netCDF output files.
!  15 Feb 2012 - R. Yantosca - Read information about nested NA grid
!  16 May 2013 - R. Yantosca - Bug fix: set doNative=T if doNestNa=T
!  19 Sep 2013 - R. Yantosca - Added nested Europe grid (EU)
!  19 Sep 2013 - R. Yantosca - Adjust reading of data fields to GEOS-FP vars
!  08 Oct 2013 - R. Yantosca - Added nested SE asia grid (SE)
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    INTEGER                 :: nPts, ios, T
    CHARACTER(LEN=MAX_CHAR) :: line, temp

    !-----------------------------------------------------------------------
    ! Read the file with the date (passed from the Perl script "doGeos5")
    !-----------------------------------------------------------------------

    ! Get day of year
    READ( 5, '(i8)', ERR=990 ) yyyymmdd

    ! Save day of year in a string
    WRITE( yyyymmdd_string, '(i8)' ) yyyymmdd

    ! 1-hourly data timestamps
    DO T = 1, TIMES_A1
       a1Mins(T)  = ( ( T - 1 ) * 60 ) + 30
    ENDDO

    ! 3-hourly timestamps
    DO T = 1, TIMES_A3
       a3MinsI(T) = ( ( T - 1 ) * 180 )
       a3Mins (T) = ( ( T - 1 ) * 180 ) + 90
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the file with the filename templates and fields to pull
    !-----------------------------------------------------------------------

    ! Open the file
    OPEN( IU_TXT, FILE='./GeosFpDriver.input', STATUS='old', ERR=999 )

    ! Read each line
    DO

       ! Read a line from the file
       READ( IU_TXT, '(a)', END=100 ) line

       ! Read the various options from the file
       SELECT CASE( TRIM( line ) )

          CASE( '==> Turn on debug print output?' )
             READ( IU_TXT,   *,      ERR=999 ) VERBOSE

          CASE( '==> Local Raw Data Path' )
             READ( IU_TXT, '(a)',    ERR=999 ) inputDataDir

          CASE( '==> Nested CH output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestCh
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestCh
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestCh
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestCh
             READ( IU_TXT,   *,      ERR=999 ) I0_ch, J0_ch, I1_ch, J1_ch
             I_NestCh = I1_ch - I0_ch + 1
             J_NestCh = J1_ch - J0_ch + 1

          CASE( '==> Nested EU output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestEu
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestEu
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestEu
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestEu
             READ( IU_TXT,   *,      ERR=999 ) I0_eu, J0_eu, I1_eu, J1_eu
             I_NestEu = I1_eu - I0_eu + 1
             J_NestEu = J1_eu - J0_eu + 1

          CASE( '==> Nested NA output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestNa
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestNa
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestNa
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestNa
             READ( IU_TXT,   *,      ERR=999 ) I0_na, J0_na, I1_na, J1_na
             I_NestNa = I1_na - I0_na + 1
             J_NestNa = J1_na - J0_na + 1

          CASE( '==> Nested SE output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestSe
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestSe
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestSe
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestSe
             READ( IU_TXT,   *,      ERR=999 ) I0_se, J0_se, I1_se, J1_se
             I_NestSe = I1_se - I0_se + 1
             J_NestSe = J1_se - J0_se + 1

          !(jxu, 2015/12/08, add nested Asia)
          CASE( '==> Nested AS output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestAs
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestAs
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestAs
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestAs
             READ( IU_TXT,   *,      ERR=999 ) I0_as, J0_as, I1_as, J1_as
             I_NestAs = I1_as - I0_as + 1
             J_NestAs = J1_as - J0_as + 1
         !(jxu, end)

          ! ----- (lzh, 06/20/2014) ---------------
          CASE( '==> Nested 0625 CH output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestCh05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestCh05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestCh05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestCh05
             READ( IU_TXT,   *,      ERR=999 ) I0_ch05, J0_ch05, I1_ch05, J1_ch05
             I_NestCh05 = I1_ch05 - I0_ch05 + 1
             J_NestCh05 = J1_ch05 - J0_ch05 + 1

          CASE( '==> Nested 0625 EU output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestEu05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestEu05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestEu05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestEu05
             READ( IU_TXT,   *,      ERR=999 ) I0_eu05, J0_eu05, I1_eu05, J1_eu05
             I_NestEu05 = I1_eu05 - I0_eu05 + 1
             J_NestEu05 = J1_eu05 - J0_eu05 + 1

          CASE( '==> Nested 0625 NA output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestNa05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestNa05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestNa05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestNa05
             READ( IU_TXT,   *,      ERR=999 ) I0_na05, J0_na05, I1_na05, J1_na05
             I_NestNa05 = I1_na05 - I0_na05 + 1
             J_NestNa05 = J1_na05 - J0_na05 + 1

          CASE( '==> Nested 0625 SE output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestSe05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestSe05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestSe05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestSe05
             READ( IU_TXT,   *,      ERR=999 ) I0_se05, J0_se05, I1_se05, J1_se05
             I_NestSe05 = I1_se05 - I0_se05 + 1
             J_NestSe05 = J1_se05 - J0_se05 + 1

          !(jxu, 2015/12/08, add nested Asia)
          CASE( '==> Nested 0625 AS output' )
             READ( IU_TXT,   *,      ERR=999 ) doNestAs05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplNestAs05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplNestAs05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplNestAs05
             READ( IU_TXT,   *,      ERR=999 ) I0_as05, J0_as05, I1_as05, J1_as05
             I_NestAs05 = I1_as05 - I0_as05 + 1
             J_NestAs05 = J1_as05 - J0_as05 + 1
          !(jxu, end)
          ! ----- (finish edit) ----------------------
          
          !(lb, 2021/03/22)
          CASE( '==> 0.5 x 0.625 output' )
             READ( IU_TXT,   *,      ERR=999 ) doGlobal05
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmplGlobal05
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmplGlobal05
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmplGlobal05

          !(jxu, 2016/02/13, add 0.25 global)
          CASE( '==> 0.25 x 0.3125 output' )
             READ( IU_TXT,   *,      ERR=999 ) do025x03125
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmpl025x03125
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmpl025x03125
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmpl025x03125
         !(jxu, end)

          CASE( '==> 2 x 2.5 output' )
             READ( IU_TXT,   *,      ERR=999 ) do2x25
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmpl2x25
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmpl2x25
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmpl2x25

          CASE( '==> 4 x 5 output' )
             READ( IU_TXT,   *,      ERR=999 ) do4x5
             READ( IU_TXT, '(a)',    ERR=999 ) dataTmpl4x5
             READ( IU_TXT, '(a)',    ERR=999 ) tempDirTmpl4x5
             READ( IU_TXT, '(a)',    ERR=999 ) dataDirTmpl4x5

          CASE( '==> const_2d_asm_Nx' )
             READ( IU_TXT, '(a)',    ERR=999 ) const_2d_asm_Nx_file
             READ( IU_TXT, '(a)',    ERR=999 ) const_2d_asm_Nx_data
             READ( IU_TXT,   *,      ERR=999 ) doMakeCn

          CASE( '==> tavg1_2d_flx_Nx' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_flx_Nx_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_flx_Nx_data

          CASE( '==> tavg1_2d_lnd_Nx' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_lnd_Nx_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_lnd_Nx_data

          CASE( '==> tavg1_2d_rad_Nx' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_rad_Nx_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_rad_Nx_data

          CASE( '==> tavg1_2d_slv_Nx' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_slv_Nx_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg1_2d_slv_Nx_data

          CASE( '==> tavg3_3d_asm_Nv' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_asm_Nv_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_asm_Nv_data

          CASE( '==> tavg3_3d_cld_Nv' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_cld_Nv_file
             READ( IU_TXT, '(7x,a)', ERR=999 ) tavg3_3d_cld_Nv_data_c
             READ( IU_TXT, '(7x,a)', ERR=999 ) tavg3_3d_cld_Nv_data_d

          CASE( '==> tavg3_3d_mst_Ne' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_mst_Ne_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_mst_Ne_data

          CASE( '==> tavg3_3d_mst_Nv' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_mst_Nv_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_mst_Nv_data

          CASE( '==> tavg3_3d_rad_Nv' )
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_rad_Nv_file
             READ( IU_TXT, '(a)',    ERR=999 ) tavg3_3d_rad_Nv_data

          CASE( '==> inst3_3d_asm_Nv' )
             READ( IU_TXT, '(a)',    ERR=999 ) inst3_3d_asm_Nv_file
             READ( IU_TXT, '(a)',    ERR=999 ) inst3_3d_asm_Nv_data

          CASE( '==> Mapping Weight Files' )
             READ( IU_TXT, '(a)',    ERR=999 ) weightFileTo05x0625   !(lzh)
             READ( IU_TXT, '(a)',    ERR=999 ) weightFileTo2x25
             READ( IU_TXT, '(a)',    ERR=999 ) weightFileTo4x5

          CASE( '==> Template Files' )
             READ( IU_TXT, '(a)',    ERR=999 ) templateFile

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

    ! Define a convenience switch for the native grid
    ! (jxu, 2015/12/08, add nested Asia by adding .or. do NestAs)
    doNative = ( doNestCh .or. doNestEu .or. doNestNa .or. doNestSe .or. doNestAs .or. do025x03125)
    !(jxu,end)  We add do 025x03125 to this native grid (Chi Li)

    ! (lzh, 06/20/2014)
    !(jxu, 2015/12/08, add nested Asia by adding .or. do NestAs)
    do05x0625 = ( doNestCh05 .or. doNestEu05 .or. doNestNa05 .or. doNestSe05 .or. doNestAs05 .or. doGlobal05)
    !(jxu,end)

    ! Mapping weights: native grid (use as placeholder for routines below)
    IF ( doNative ) THEN
       nPts = 1
       CALL ReadMappingWeights( '',                                       &
                                I025x03125,  J025x03125, nPts, mapNative )
    ENDIF

    ! (lzh,06/20/2014) Mapping weights to 0.5 x 0.625 grid
    IF ( do05x0625 ) THEN
       nPts = ( I025x03125 / I05x0625 ) + 2
       CALL ReadMappingWeights( weightFileTo05x0625,                         &
                       I05x0625,      J05x0625,       nPts, mapTo05x0625 )
    ENDIF

    ! Mapping weights to 2 x 2.5 grid
    IF ( do2x25 ) THEN
       nPts = ( I025x03125 / I2x25 ) + 2
       CALL ReadMappingWeights( weightFileTo2x25,                         &
                                I2x25,      J2x25,       nPts, mapTo2x25 )
    ENDIF

    ! (lzh,06/20/2014) Mapping weights to 0.5 x 0.625 grid
    IF ( do05x0625 ) THEN
       nPts = ( I025x03125 / I05x0625 ) + 2
       CALL ReadMappingWeights( weightFileTo05x0625,                         &
                       I05x0625,      J05x0625,       nPts, mapTo05x0625 )
    ENDIF

    ! Mapping weights to 4 x 5 grid
    IF ( do4x5 ) THEN
       nPts = ( I025x03125 / I4x5 ) + 2
       CALL ReadMappingWeights( weightFileTo4x5,                          &
                                I4x5,       J4x5,        nPts, mapTo4x5  )
    ENDIF

    !-----------------------------------------------------------------------
    ! Read data from template files
    !-----------------------------------------------------------------------

    ! FRLANDICE data (for SNOMAS regridding)
    CALL ReadTemplateFile( templateFile, lwiMask, frLand, frLandIce )

    !-----------------------------------------------------------------------
    ! Define hybrid-grid index arrays
    !-----------------------------------------------------------------------

    ! Ap [hPa] for 72 levels (73 edges)
    Ap =  (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
             1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
             4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
             7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
             1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
             1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
             2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
             2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02, &
             1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
             7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01, &
             4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01, &
             1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01, &
             9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00, &
             4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00, &
             1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01, &
             6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01, &
              2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02, &
             6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02, &
             1.000000d-02 /)

    ! Bp [unitless] for 72 levels (73 edges)
    Bp =  (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01, &
             9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01, &
             8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01, &
             7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01, &
             6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01, &
             4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01, &
             2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01, &
             6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
             0.000000d+00 /)

    !-----------------------------------------------------------------------
    ! Verbose output for debugging
    !-----------------------------------------------------------------------
    IF ( VERBOSE ) THEN
       PRINT*, 'YYYYMMDD        : ', yyyymmdd
       PRINT*, 'aMins           : ', a1Mins
       PRINT*, 'a3MinsI         : ', a3MinsI
       PRINT*, 'a3Mins          : ', a3Mins
       PRINT*, 'doNative        : ', doNative
       PRINT*, 'doNestCh        : ', doNestCh
       PRINT*, ' I0, J0, I1, J1 : ', I0_ch, J0_ch, I1_ch, J1_ch
       PRINT*, ' ICH, JCH       : ', I_NestCh, J_NestCh
       PRINT*, 'doNestEu        : ', doNestEu
       PRINT*, ' I0, J0, I1, J1 : ', I0_eu, J0_eu, I1_eu, J1_eu
       PRINT*, ' INA, JNA       : ', I_NestEu, J_NestEu
       PRINT*, 'doNestNa        : ', doNestNa
       PRINT*, ' I0, J0, I1, J1 : ', I0_na, J0_na, I1_na, J1_na
       PRINT*, ' INA, JNA       : ', I_NestNa, J_NestNa
       PRINT*, 'doNestSe        : ', doNestSe
       PRINT*, ' I0, J0, I1, J1 : ', I0_se, J0_se, I1_se, J1_se
       PRINT*, ' ISE, JSE       : ', I_NestSe, J_NestSe
!(jxu, 2015/12/08, add nested Asia)
       PRINT*, 'doNestAs        : ', doNestAs
       PRINT*, ' I0, J0, I1, J1 : ', I0_as, J0_as, I1_as, J1_as
       PRINT*, ' IIN, JIN       : ', I_NestAs, J_NestAs
!(jxu, end)
!! (lzh, 05/25/2014)
       PRINT*, 'doNestCh05        : ', doNestCh05
       PRINT*, ' I0, J0, I1, J1 : ', I0_ch05, J0_ch05, I1_ch05, J1_ch05
       PRINT*, ' ICH, JCH       : ', I_NestCh05, J_NestCh05
       PRINT*, 'doNestEu05        : ', doNestEu05
       PRINT*, ' I0, J0, I1, J1 : ', I0_eu05, J0_eu05, I1_eu05, J1_eu05
       PRINT*, ' INA, JNA       : ', I_NestEu05, J_NestEu05
       PRINT*, 'doNestNa05        : ', doNestNa05
       PRINT*, ' I0, J0, I1, J1 : ', I0_na05, J0_na05, I1_na05, J1_na05
       PRINT*, ' INA, JNA       : ', I_NestNa05, J_NestNa05
       PRINT*, 'doNestSe05        : ', doNestSe05
       PRINT*, ' I0, J0, I1, J1 : ', I0_se05, J0_se05, I1_se05, J1_se05
       PRINT*, ' ISE, JSE       : ', I_NestSe05, J_NestSe05
!(jxu, 2015/12/08, add nested Asia)
       PRINT*, 'doNestAs05        : ', doNestAs05
       PRINT*, ' I0, J0, I1, J1 : ', I0_as05, J0_as05, I1_as05, J1_as05
       PRINT*, ' IIN, JIN       : ', I_NestAs05, J_NestAs05
!(jxu, end)
!! --end edit---------
       !(jxu, 2016/02/13, add 0.25 global)
       PRINT*, 'do025x03125     : ', do025x03125
       PRINT*, 'doGlobal05      : ', doGlobal05
       !(jxu, end)
       PRINT*, 'do2x25          : ', do2x25
       PRINT*, 'do4x5           : ', do4x5
       PRINT*, 'doMakeCn        : ', doMakeCn
       PRINT*, 'inputDataDir     : ', TRIM( inputDataDir           )
       PRINT*, 'dataTmplNestCh  : ', TRIM( dataTmplNestCh          )
       PRINT*, 'tempDirNestCh   : ', TRIM( tempDirTmplNestCh       )
       PRINT*, 'dataDirNestCh   : ', TRIM( dataDirTmplNestCh       )
       PRINT*, 'dataTmplNestEu  : ', TRIM( dataTmplNestEu          )
       PRINT*, 'tempDirNestEu   : ', TRIM( tempDirTmplNestEu       )
       PRINT*, 'dataDirNestEu   : ', TRIM( dataDirTmplNestEu       )
       PRINT*, 'dataTmplNestNa  : ', TRIM( dataTmplNestNa          )
       PRINT*, 'tempDirNestNa   : ', TRIM( tempDirTmplNestNa       )
       PRINT*, 'dataDirNestNa   : ', TRIM( dataDirTmplNestNa       )
       PRINT*, 'dataTmplNestSe  : ', TRIM( dataTmplNestSe          )
       PRINT*, 'tempDirNestSe   : ', TRIM( tempDirTmplNestSe       )
       PRINT*, 'dataDirNestSe   : ', TRIM( dataDirTmplNestSe       )
!(jxu, 2015/12/08, add nested Asia)
       PRINT*, 'dataTmplNestAs  : ', TRIM( dataTmplNestAs          )
       PRINT*, 'tempDirNestAs   : ', TRIM( tempDirTmplNestAs       )
       PRINT*, 'dataDirNestAs   : ', TRIM( dataDirTmplNestAs       )
!(jxu, end)
       PRINT*, 'dataTmplGlobal05    : ', TRIM( dataTmplGlobal05          )
       PRINT*, 'tempDirTmplGlobal05 : ', TRIM( tempDirTmplGlobal05       )
       PRINT*, 'dataDirTmplGlobal05 : ', TRIM( dataDirTmplGlobal05       )
!(jxu, 2016/02/13, add 0.25 global)
       PRINT*, 'dataTmpl025x03125    : ', TRIM( dataTmpl025x03125            )
       PRINT*, 'tempDirTmpl025x03125 : ', TRIM( tempDirTmpl025x03125         )
       PRINT*, 'dataDirTmpl025x03125 : ', TRIM( dataDirTmpl025x03125         )
!(jxu, end)
       PRINT*, 'dataTmpl2x25    : ', TRIM( dataTmpl2x25            )
       PRINT*, 'tempDirTmpl2x25 : ', TRIM( tempDirTmpl2x25         )
       PRINT*, 'dataDirTmpl2x25 : ', TRIM( dataDirTmpl2x25         )
       PRINT*, 'dataFile4x5     : ', TRIM( dataTmpl4x5             )
       PRINT*, 'tempDirTmpl4x5  : ', TRIM( tempDirTmpl4x5          )
       PRINT*, 'dataDirFile4x5  : ', TRIM( dataDirTmpl4x5          )
       PRINT*, 'const_2d_asm_Nx : ', TRIM( const_2d_asm_Nx_file    )
       PRINT*, '                  ', TRIM( const_2d_asm_Nx_data    )
       PRINT*, 'tavg1_2d_flx_Nx : ', TRIM( tavg1_2d_flx_Nx_file    )
       PRINT*, '                  ', TRIM( tavg1_2d_flx_Nx_data    )
       PRINT*, 'tavg1_2d_lnd_Nx : ', TRIM( tavg1_2d_lnd_Nx_file    )
       PRINT*, '                  ', TRIM( tavg1_2d_lnd_Nx_data    )
       PRINT*, 'tavg1_2d_rad_Nx : ', TRIM( tavg1_2d_rad_Nx_file    )
       PRINT*, '                  ', TRIM( tavg1_2d_rad_Nx_data    )
       PRINT*, 'tavg1_2d_slv_Nx : ', TRIM( tavg1_2d_slv_Nx_file    )
       PRINT*, '                  ', TRIM( tavg1_2d_slv_Nx_data    )
       PRINT*, 'tavg3_3d_asm_Nv : ', TRIM( tavg3_3d_asm_Nv_file    )
       PRINT*, '                  ', TRIM( tavg3_3d_asm_Nv_data    )
       PRINT*, 'tavg3_3d_cld_Nv : ', TRIM( tavg3_3d_cld_Nv_file    )
       PRINT*, '                  ', TRIM( tavg3_3d_cld_Nv_data_c  )
       PRINT*, '                  ', TRIM( tavg3_3d_cld_Nv_data_d  )
       PRINT*, 'tavg3_3d_mst_Nv : ', TRIM( tavg3_3d_mst_Ne_file    )
       PRINT*, '                  ', TRIM( tavg3_3d_mst_Ne_data    )
       PRINT*, 'tavg3_3d_mst_Ne : ', TRIM( tavg3_3d_mst_Nv_file    )
       PRINT*, '                  ', TRIM( tavg3_3d_mst_Nv_data    )
       PRINT*, 'inst3_3d_asm_Nv : ', TRIM( inst3_3d_asm_Nv_file    )
       PRINT*, '                  ', TRIM( inst3_3d_asm_Nv_data    )
       PRINT*, 'WeightsNxTo2x25 : ', TRIM( weightFileTo2x25        )
       PRINT*, 'WeightsNxTo4x5  : ', TRIM( weightFileTo4x5         )
       PRINT*, 'lwiMaskFile     : ', TRIM( templateFile            )
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

  END SUBROUTINE GeosFpInitialize
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
! !REMARKS:
!   If the MAP object is not defined, ReadMappingWeights will allocate
!   it and initialize it here.  The user is responsible for deallocating
!   it elsewhere.
!
! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version, based on MERRA
!  06 Jan 2012 - R. Yantosca - Now can define a placeholder map object
!                              for the native grid
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
    ! Special handling: fill placeholder map object for native grid
    !========================================================================
    IF ( IMX == I025x03125 .and. JMX == J025x03125 ) THEN

       ! Here the fine grid is equal to the coarse grid,
       ! so xInd=1, yInd=1, and weight=1.0
       DO J = 1, JMX
       DO I = 1, IMX
          map(I,J)%xInd    = 1
          map(I,J)%yInd    = 1
          map(I,J)%weight  = 1e0
       ENDDO
       ENDDO

       ! Return w/o reading data
       RETURN
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
    ELSE IF ( IMX == 576 .and. JMX == 361 ) THEN
       fmtStr = '(3x,12f8.4)'                    ! 0.5 x 0.625 !(lzh,06/20/2014)
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
! !DESCRIPTION: This routin reads template data (e.g. LWI, FRLAND, FRLANDICE)
!  from a netCDF file into a data array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadTemplateFile( fileName, lwiMask, frLand, frLandIce )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: fileName
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: lwiMask  (:,:)
    REAL*4,           INTENT(OUT) :: frLand   (:,:)
    REAL*4,           INTENT(OUT) :: frLandIce(:,:)
!
! !REVISION HISTORY:
!  05 Jan 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: fId

    ! Arrays
    INTEGER :: ct3d(3), st3d(3)

    ! Open netCDF file for input
    CALL NcOp_Rd( fId, TRIM( fileName ) )

    ! netCDF indices (we know that the file is 0.25 x 0.3125)
    st3d = (/ 1,          1,          1 /)
    ct3d = (/ I025x03125, J025x03125, 1 /)

    ! Read data
    CALL NcRd( lwiMask,   fId, 'LWI',      st3d, ct3d )
    CALL NcRd( frLand,    fId, 'FRLAND',   st3d, ct3d )
    CALL NcRd( frLandIce, fId, 'FRLANDIC', st3d, ct3d )

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE ReadTemplateFile
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GeosFpCleanup
!
! !DESCRIPTION: This routine deallocates all previously-allocated
!  module arrays and pointer objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFpCleanup
!
! !REVISION HISTORY:
!  25 Oct 2011 - R. Yantosca - Initial version, based on MERRA
!  06 Jan 2012 - R. Yantosca - Now deallocate the native-grid map object
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: I, J

    !======================================================================
    ! Deallocate native-grid mapping weight objects
    !======================================================================
    IF ( doNative ) THEN

       ! Echo info
       IF ( VERBOSE ) WRITE( 6, 100 )
100    FORMAT( 'Deallocating mapping weight objects for native grid' )

       ! Loop over 2 x 2.5 boxes
       DO J = 1, J025x03125
       DO I = 1, I025x03125

          !-------------------------------------------------
          ! Deallocate Nx grid to 2 x 2.5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapNative(I,J)%xInd  ) ) THEN
             DEALLOCATE( mapNative(I,J)%xInd )
          ENDIF

          IF ( ASSOCIATED( mapNative(I,J)%yInd  ) ) THEN
             DEALLOCATE( mapNative(I,J)%yInd )
          ENDIF

          IF ( ASSOCIATED( mapNative(I,J)%weight) ) THEN
             DEALLOCATE( mapNative(I,J)%weight )
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapNative ) ) DEALLOCATE( mapNative )
    ENDIF

!----- (lzh, 06/20/2014)-----
    !======================================================================
    ! Deallocate 0.5x0.625 mapping weight objects
    !======================================================================
    IF ( do05x0625 ) THEN

       ! Echo info
       IF ( VERBOSE ) WRITE( 6, 101 )
101    FORMAT( 'Deallocating mapping weight objects for 0.5 x 0.625 grid' )

       ! Loop over 0.5 x 0.625 boxes
       DO J = 1, J05x0625
       DO I = 1, I05x0625

          !-------------------------------------------------
          ! Deallocate Nx grid to 0.5 x 0.625 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapTo05x0625(I,J)%xInd  ) ) THEN
             DEALLOCATE( mapTo05x0625(I,J)%xInd )
          ENDIF

          IF ( ASSOCIATED( mapTo05x0625(I,J)%yInd  ) ) THEN
             DEALLOCATE( mapTo05x0625(I,J)%yInd )
          ENDIF

          IF ( ASSOCIATED( mapTo05x0625(I,J)%weight) ) THEN
             DEALLOCATE( mapTo05x0625(I,J)%weight )
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapTo05x0625 ) ) DEALLOCATE( mapTo05x0625 )
    ENDIF
!-----(finish edit)----------

    !======================================================================
    ! Deallocate 2x25 mapping weight objects
    !======================================================================
    IF ( do2x25 ) THEN

       ! Echo info
       IF ( VERBOSE ) WRITE( 6, 110 )
110    FORMAT( 'Deallocating mapping weight objects for 2 x 2.5 grid' )

       ! Loop over 2 x 2.5 boxes
       DO J = 1, J2x25
       DO I = 1, I2x25

          !-------------------------------------------------
          ! Deallocate Nx grid to 2 x 2.5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapTo2x25(I,J)%xInd  ) ) THEN
             DEALLOCATE( mapTo2x25(I,J)%xInd )
          ENDIF

          IF ( ASSOCIATED( mapTo2x25(I,J)%yInd  ) ) THEN
             DEALLOCATE( mapTo2x25(I,J)%yInd )
          ENDIF

          IF ( ASSOCIATED( mapTo2x25(I,J)%weight) ) THEN
             DEALLOCATE( mapTo2x25(I,J)%weight )
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapTo2x25 ) ) DEALLOCATE( mapTo2x25 )
    ENDIF

    !======================================================================
    ! Deallocate 4x5 mapping weight objects
    !======================================================================
    IF ( do4x5 ) THEN

       IF ( VERBOSE ) WRITE( 6, 120 )
120    FORMAT( 'Deallocating mapping weight objects for 4 x 5 grid' )

       ! Loop over 4 x 5 boxes
       DO J = 1, J4x5
       DO I = 1, I4x5

          !-------------------------------------------------
          ! Deallocate Nx grid to 4x5 object fields
          !-------------------------------------------------
          IF ( ASSOCIATED( mapTo4x5(I,J)%xInd  ) ) THEN
             DEALLOCATE( mapTo4x5(I,J)%xInd )
          ENDIF

          IF ( ASSOCIATED( mapTo4x5(I,J)%yInd  ) ) THEN
             DEALLOCATE( mapTo4x5(I,J)%yInd )
          ENDIF

          IF ( ASSOCIATED( mapTo4x5(I,J)%weight) ) THEN
             DEALLOCATE( mapTo4x5(I,J)%weight )
          ENDIF

       ENDDO
       ENDDO

       ! Free the objects themselves
       IF ( ASSOCIATED( mapTo4x5 ) ) DEALLOCATE( mapTo4x5 )
    ENDIF

  END SUBROUTINE GeosFpCleanup
!EOC
END MODULE GeosFpInputsModule
