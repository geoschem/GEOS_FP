!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Geos57RegridModule
!
! !DESCRIPTION: Module Geos57RegridModule contains arrays and variables used 
!  to regrid the GEOS-5 data from 0.5 x 0.666 to 1 x 1, 1 x 1.25, 2 x 2.5, 
!  and 4 x 5 resolution. (bmy, 3/3/04)
!
!\subsection*{Overview}
!  GeosRegridModule uses the regridding software MAP\_A2A from S-J Lin.  
!  This is area-preserving mapping.  For example, if you have a quantity 
!  such as kg/m2/s or W/m2, MAP\_A2A will multiply by the area on the
!  input grid, then regrid, and divide by the area on the output grid,
!  such that the total quantity is preserved.  
!
!\subsection*{Regridding of the winds}
!
!  Comments from Claire Carouge (ccarouge@seas.harvard.edu)
!
! \begin{quote}
!  There are two types of physical quantities:
!  \begin{description}
!  \item[Intensive] A quantity whose value doesn't 
!                   change with the grid cell size.
!  \item[Extensive] A quantity whose value changes with the
!                   grid cell size. 
!  \end{description} 
!
!  For example, if you have a uniformly distributed atmosphere and a
!  grid and you increase the resolution by 2 (grid cells are 4 times
!  smaller), then temperature, velocities, concentrations, pressure won't
!  change at each grid cell but the mass will be smaller for each grid
!  cell. So mass is an extensive variable and the others are intensives.
!
!  The MAP\_A2A algorithm is set up to regrid extensive quantities.  So
!  in order for us to use it to regrid winds, we must first multiply the 
!  winds by the pressure on the input grid in order to create an extensive 
!  quantity (e.g. a "mass flux").  Then we must divide by the pressure on
!  the output grid to convert back to a wind.
!
!  In other regridding routines, (e.g. NCREGRID), you can specify if 
!  the quantity you are regridding is intensive or extensive, and it
!  will do the regridding accordingly.
!
!  I've been looking at the MAP\_A2A regridding algorithm to answer two 
!  different questions:
!
! \begin{enumerate}
!  \item \textbf{Why does it look so complicated?}  Apparently we are using a 
!  more elaborate regridding algorithm than NCREGRID. It probably has some 
!  nice qualities that NCREGRID hasn't. Considering the small differences 
!  seen on Bob's graphs, I don't think it is essential to use this algorithm 
!  but it can't hurt.
!
!  \item \textbf{Is the treatment of the poles coherent with TPCORE?}
!  In the regridding algorithm, there is some special treatment on the
!  poles. This is due to the fact, the algorithm is based on the
!  calculation of the slopes between neighbor grid cells. So we need
!  different values from neighbor grid cells. The problem with the poles is
!  to define the grid cells $j-1$ (resp. $j+1$) at the South Pole (resp. North
!  Pole). The method used in the regridding algorithm is the same than the
!  one used in tpcore: to access the grid cell $j-1$ at South Pole you need
!  to go southward. If you start from the cell, $i=1$ then you arrive in the
!  cell $i=im/2+1$ (sphere). So in the code the values are "crossed" at the
!  poles. Also, they don't take the $j=1$ grid cells but the $j=2$ because the
!  poles are supposed to be only one circular grid cell, so the points
!  $(i=1, j=1)$ and $(i=im/2+1, j=1)$ are supposed to be the same grid cell. 
!  So the cell $j-1$ is the cell $(i=im/2+1, j=2)$ ...
!
!  In addition, the wind values for the pole grid cells are averaged at the
!  end. This is also coherent with tpcore as we consider the poles are one
!  circular grid cell.
! \end{enumerate}
! 
!  So with the explanation of the multiplication by P, I think everything
!  is done correctly in the regridding. The differences we see are
!  explained by the differences in algorithms used but we can't say one is
!  better than the other.
! \end{quote}
!
! !INTERFACE: 
!
MODULE Geos57RegridModule
! 
! !USES:
!
  USE Geos57InputsModule, ONLY : I025x03125, J025x03125, L025x03125
  USE Geos57InputsModule, ONLY : I05x0666,   J05x0666,   L05x0666
  USE Geos57InputsModule, ONLY : I1x125,     J1x125,     L1x125  
  USE Geos57InputsModule, ONLY : I2x25,      J2x25,      L2x25   
  USE Geos57InputsModule, ONLY : I4x5,       J4x5,       L4x5    

  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE        :: map_a2a
  PRIVATE        :: YMAP
  PRIVATE        :: ppm_lat
  PRIVATE        :: xmap
  PRIVATE        :: ppm_cycle
  PRIVATE        :: lmppm
  PRIVATE        :: huynh
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC         :: RegridGeos57to1x125
  PUBLIC         :: RegridGeos57to2x25
  PUBLIC         :: RegridGeos57to4x5
  PUBLIC         :: Geos57RegridInit
!
! !PUBLIC DATA MEMBERS:
! 
  !--------------------------
  ! 0.25 x 0.3125 resolution
  !--------------------------
  REAL*4, PUBLIC :: xedge_025x03125( I025x03125 + 1 )  ! Lon edges
  REAL*4, PUBLIC :: yedge_025x03125( J025x03125 + 1 )  ! Lat edges
  REAL*4, PUBLIC :: sine_025x03125 ( J025x03125 + 1 )  ! SIN( lat edges ) 
  REAL*4, PUBLIC :: xmid_025x03125 ( I025x03125     )  ! Lon centers
  REAL*4, PUBLIC :: ymid_025x03125 ( J025x03125     )  ! Lat centers
  REAL*4, PUBLIC :: zmid_025x03125 ( L025x03125     )  ! Vertical levels

  !--------------------------
  ! 0.5 x 0.666 resolution
  !--------------------------
  REAL*4, PUBLIC :: xedge_05x0666( I05x0666 + 1 )      ! Lon edges
  REAL*4, PUBLIC :: yedge_05x0666( J05x0666 + 1 )      ! Lat edges
  REAL*4, PUBLIC :: sine_05x0666 ( J05x0666 + 1 )      ! SIN( lat edges ) 
  REAL*4, PUBLIC :: xmid_05x0666 ( I05x0666     )      ! Lon centers
  REAL*4, PUBLIC :: ymid_05x0666 ( J05x0666     )      ! Lat centers
  REAL*4, PUBLIC :: zmid_05x0666 ( L05x0666     )      ! Vertical levels

  !--------------------------
  ! 1 x 1.25 resolution
  !--------------------------
  REAL*4, PUBLIC :: xedge_1x125( I1x125 + 1 )          ! Lon edges
  REAL*4, PUBLIC :: yedge_1x125( J1x125 + 1 )          ! Lat edges
  REAL*4, PUBLIC :: sine_1x125 ( J1x125 + 1 )          ! SIN( lat edges )
  REAL*4, PUBLIC :: xmid_1x125 ( I1x125     )          ! Lon centers
  REAL*4, PUBLIC :: ymid_1x125 ( J1x125     )          ! Lat centers
  REAL*4, PUBLIC :: zmid_1x125 ( L1x125     )          ! Vertical levels


  !--------------------------
  ! 2 x 2.5 resolution 
  !--------------------------
  REAL*4, PUBLIC :: xedge_2x25( I2x25 + 1 )            ! Lon edges
  REAL*4, PUBLIC :: yedge_2x25( J2x25 + 1 )            ! Lat edges
  REAL*4, PUBLIC :: sine_2x25 ( J2x25 + 1 )            ! SIN( lat edges )
  REAL*4, PUBLIC :: xmid_2x25 ( I2x25     )            ! Lon centers
  REAL*4, PUBLIC :: ymid_2x25 ( J2x25     )            ! Lat centers
  REAL*4, PUBLIC :: zmid_2x25 ( L2x25     )            ! Vertical levels

  !--------------------------
  ! 4 x 5 resolution
  !--------------------------
  REAL*4, PUBLIC :: xedge_4x5( I4x5 + 1 )              ! Lon edges  
  REAL*4, PUBLIC :: yedge_4x5( J4x5 + 1 )              ! Lat edges
  REAL*4, PUBLIC :: sine_4x5 ( J4x5 + 1 )              ! SIN( lat edges )
  REAL*4, PUBLIC :: xmid_4x5 ( I4x5     )              ! Lon centers
  REAL*4, PUBLIC :: ymid_4x5 ( J4x5     )              ! Lat centers
  REAL*4, PUBLIC :: zmid_4x5 ( L4x5     )              ! Vertical levels
!
! !AUTHOR:
! Original MAP_A2A code from S-J Lin
! Modified by Bob Yantosca and placed into F90 module format
!
! !REMARKS:
!  The GEOS-5.7.x 0.25 x 0.3125 grid is centered on (-180,-90).
!  The GEOS-Chem  1.0  x 1.25   grid is centered on (-180,-90).
!  The GEOS-Chem  2.0  x 2.5    grid is centered on (-180,-90).
!  The GEOS-Chem  4.0  x 5.0    grid is centered on (-180,-90).
!
! !REVISION HISTORY:
!  23 Jul 2010 - R. Yantosca - Initial version, based on MerraRegridModule.f90
!  26 Oct 2011 - R. Yantosca - Now make all lon & lat arrays public
!  03 Jan 2012 - R. Yantosca - Now define index arrays for vertical levels
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Degrees to Radians
  REAL*4, PARAMETER  :: D2R = 3.141592658979323d0 / 180d0

  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegridGeos57to1x125
!
! !DESCRIPTION: Subroutine RegridGeos57to1x125 is a wrapper for MAP\_A2A.  
!  It is called to regrid from the GEOS-5.7.x chemistry forcing ("F") 
!  grid (1 x 1.25) to the GEOS-Chem 4 x 5 grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RegridGeos57to1x125( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    ! IV = 0 is scalar field; IV = 1 is vector field
    INTEGER, INTENT(IN)  :: iv

    ! Input data on 1 x 1.25 grid
    REAL*4,  INTENT(IN)  :: q1(I025x03125,J025x03125)
!
! !OUTPUT PARAMETERS:
!
    ! Output data on 4 x 5 grid
    REAL*4,  INTENT(OUT) :: q2(I1x125,J1x125)
!
! !REVISION HISTORY: 
!  26 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER :: T

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0e0 ) ) THEN
       q2 = 0e0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I025x03125,     J025x03125,  xedge_025x03125,  &
                  sine_025x03125, q1,          I1x125,           &
                  J1x125,         xedge_1x125, sine_1x125,       &
                  q2,             0,           iv               )

  END SUBROUTINE RegridGeos57to1x125
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegridGeos57To2x25
!
! !DESCRIPTION: Subroutine RegridGeos57To2x25 is a wrapper for MAP\_A2A.
!  It is called to regrid from the GEOS-5.7.x native grid (0.25 x0 x 0.3125)
!  to the GEOS-Chem 2 x 2.5 grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RegridGeos57To2x25( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    ! IV = 0 is scalar field; IV = 1 is vector field
    INTEGER, INTENT(IN)  :: iv
    
    ! Input data on 0.5 x 0.666 grid
    REAL*4,  INTENT(IN)  :: q1(I025x03125,J025x03125)  
!
! !OUTPUT PARAMETERS:
!
    ! Output data on 2 x 2.5 grid
    REAL*4,  INTENT(OUT) :: q2(I2x25,J2x25)
!
! !REVISION HISTORY: 
!  26 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: T

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0e0 ) ) THEN
       q2 = 0e0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I025x03125,     J025x03125,  xedge_025x03125,  &
                  sine_025x03125, q1,          I2x25,            &
                  J2x25,          xedge_2x25,  sine_2x25,        &
                  q2,             0,           iv               )

  END SUBROUTINE RegridGeos57To2x25
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegridGeos57NTo4x5
!
! !DESCRIPTION: Subroutine RegridGeos57NTo2x25 is a wrapper for MAP\_A2A.
!  It is called to regrid from the GEOS-5.7.x native ("N") grid (0.5 x 0.666)
!  to the GEOS-Chem 4 x 5 grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RegridGeos57to4x5( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    ! IV = 0 is scalar field; IV = 1 is vector field
    INTEGER, INTENT(IN)  :: iv

    ! Input data on 0.5 x 0.666 grid
    REAL*4,  INTENT(IN)  :: q1(I025x03125,J025x03125)  
!
! !OUTPUT PARAMETERS:
!
    ! Output data on 4 x 5 grid
    REAL*4,  INTENT(OUT) :: q2(I4x5,J4x5)
!
! !REVISION HISTORY: 
!  26 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER :: T

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0e0 ) ) THEN
       q2 = 0e0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I025x03125,     J025x03125,  xedge_025x03125,  &
                  sine_025x03125, q1,          I4x5,             &
                  J4x5,           xedge_4x5,   sine_4x5,         &
                  q2,             0,           iv               )

  END SUBROUTINE RegridGeos57to4x5
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A
!
! !DESCRIPTION: Subroutine MAP\_A2A is a orizontal arbitrary grid to arbitrary 
!  grid conservative high-order mapping regridding routine by S-J Lin.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_A2A( im, jm, lon1, sin1, q1, &
                      in, jn, lon2, sin2, q2, ig, iv )
!
! !INPUT PARAMETERS:
!

    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole; 
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig 

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    !  Quantity on INPUT grid
    REAL*4,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:

    ! Regridded quantity on OUTPUT grid
    REAL*4,  INTENT(OUT) :: q2(in,jn)
!
! !AUTHOR:
!   Original subroutine by S-J Lin (GSFC)
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i,j,k
    REAL*4               :: qtmp(in,jm)

    !===================================================================
    ! MAP_A2A begins here!
    !
    ! Mapping in the E-W direction
    ! If both grids have the same longitude dimension, don't call XMAP
    !===================================================================    
    IF ( im .eq. in ) THEN
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
    ELSE
       CALL xmap(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )
    ENDIF
    
    !===================================================================
    ! Mapping in the N-S direction
    ! If both grids have the same latitude dimension, don't call YMAP 
    !===================================================================    
    IF ( jm .eq. jn ) THEN
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
    ELSE
       CALL ymap(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)
    ENDIF

  END SUBROUTINE Map_A2A
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ymap
!
! !DESCRIPTION: Subroutine Ymap performs area preserving mapping in N-S from 
!  an arbitrary resolution to another.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Ymap( im, jm, sin1, q1, jn, sin2, q2, ig, iv )

!
! !INPUT PARAMETERS:
!

    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im            

    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm            

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn           

    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig            

    ! IV=0: scalar; 
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv            

    ! Original southern edge of the cell sin(lat1)  
    REAL*4,  INTENT(IN)  :: sin1(jm+1-ig) 
                                            
    ! Original data at center of the cell
    REAL*4,  INTENT(IN)  :: q1(im,jm)      
                                          
    ! Target cell's southern edge sin(lat2)
    REAL*4,  INTENT(IN)  :: sin2(jn+1-ig) 
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT) :: q2(im,jn)     
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: S.-J. Lin
!   First version: piece-wise constant mapping
!   Apr 1, 2000
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*4               :: al(im,jm), ar(im,jm), a6(im,jm), dy1(jm)
    REAL*4,  PARAMETER   :: r3 = 1./3., r23 = 2./3. 
    REAL*4               :: pl, pr, qsum, esl, dy, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    ! Construct subgrid PP distribution
    call ppm_lat(im, jm, ig, q1, al, ar, a6, 3, iv)
    
    do 1000 i=1,im
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             pl = (sin2(j)-sin1(m)) / dy1(m)
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                pr = (sin2(j+1)-sin1(m)) / dy1(m)
                q2(i,j) = al(i,m) + 0.5*(a6(i,m)+ar(i,m)-al(i,m)) &
               &                    *(pr+pl)-a6(i,m)*r3*(pr*(pr+pl)+pl**2)
                j0 = m
                goto 555
             else

                ! South most fractional area
                qsum = (sin1(m+1)-sin2(j))*(al(i,m)+0.5*(a6(i,m)+ &
                &              ar(i,m)-al(i,m))*(1.+pl)-a6(i,m)*  &
                &               (r3*(1.+pl*(1.+pl))))

                do mm=m+1,jm-ig

                   ! locate the eastern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else

                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      esl = dy / dy1(mm)
                      qsum = qsum + dy*(al(i,mm)+0.5*esl* &
                     &       (ar(i,mm)-al(i,mm)+a6(i,mm)*(1.-r23*esl)))
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue

    !===================================================================
    ! Final processing for poles
    !===================================================================
    if ( ig .eq. 0 .and. iv .eq. 0 ) then
        
       ! South pole
       sum = 0.
       do i=1,im
          sum = sum + q2(i,1)
       enddo
        
       sum = sum / float(im)
       do i=1,im
          q2(i,1) = sum
       enddo
        
       ! North pole:
       sum = 0.
       do i=1,im
          sum = sum + q2(i,jn)
       enddo
        
       sum = sum / float(im)
       do i=1,im
          q2(i,jn) = sum
       enddo

    endif

  END SUBROUTINE Ymap
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ppm_Lat
!
! !DESCRIPTION: Subroutine Ppm\_Lat is called by Ymap.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Ppm_Lat( im, jm, ig, q, al, ar, a6, jord, iv )
!
! !INPUT PARAMETERS:
!
    ! Arguments
    INTEGER           :: im, jm          ! Dimensions
    INTEGER           :: ig              ! ig=0: scalar pole to pole
                                         ! ig=1: D-grid u-wind; 
                                         !   not defined at poles because 
                                         !   of staggering
    REAL*4            :: q(im,jm-ig)
    INTEGER           :: jord
    INTEGER           :: iv              ! iv=0 scalar
                                         ! iv=1 vector
!
! !OUTPUT PARAMETERS:
!
    REAL*4            :: al(im,jm-ig)
    REAL*4            :: ar(im,jm-ig)
    REAL*4            :: a6(im,jm-ig)
!
! !AUTHOR:
!  Written by S-J Lin
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*4            :: dm(im,jm-ig)
    REAL*4, PARAMETER :: r3 = 1./3. 
    INTEGER           :: i, j, im2, iop, jm1
    REAL*4            :: tmp, qmax, qmin, qop
    
    ! PPM_LAT begins here
    ! Compute dm: linear slope
    do j=2,jm-1-ig
       do i=1,im
          dm(i,j) = 0.25*(q(i,j+1) - q(i,j-1))
          qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
          qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
          dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
       enddo
    enddo

    im2 = im/2
    jm1 = jm - 1

    ! Poles:
    if (iv .eq. 1 ) then

       !===============================================================
       ! u-wind (ig=1)
       ! v-wind (ig=0)
       !===============================================================

       ! SP
       do i=1,im
          if( i .le. im2) then
             qop = -q(i+im2,2-ig)
          else
             qop = -q(i-im2,2-ig)
          endif
          tmp = 0.25*(q(i,2) - qop)
          qmax = max(q(i,2),q(i,1), qop) - q(i,1)
          qmin = q(i,1) - min(q(i,2),q(i,1), qop)
          dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
       
       ! NP
       do i=1,im
          if( i .le. im2) then
             qop = -q(i+im2,jm1)
          else
             qop = -q(i-im2,jm1)
          endif
          tmp = 0.25*(qop - q(i,jm1-ig))
          qmax = max(qop,q(i,jm-ig), q(i,jm1-ig)) - q(i,jm-ig)
          qmin = q(i,jm-ig) - min(qop,q(i,jm-ig), q(i,jm1-ig))
          dm(i,jm-ig) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
    else
        
       !===============================================================
       ! Scalar:
       ! This code segment currently works only if ig=0
       !===============================================================

       ! SP
       do i=1,im2
          tmp = 0.25*(q(i,2)-q(i+im2,2))
          qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
          qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
          dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
        
       do i=im2+1,im
          dm(i, 1) =  - dm(i-im2, 1)
       enddo

       ! NP
       do i=1,im2
          tmp = 0.25*(q(i+im2,jm1)-q(i,jm1))
          qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
          qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
          dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo

       do i=im2+1,im
          dm(i,jm) =  - dm(i-im2,jm)
       enddo
     endif
      
     do j=2,jm-ig
        do i=1,im
           al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
     enddo
      
     do j=1,jm-1-ig
        do i=1,im
           ar(i,j) = al(i,j+1)
        enddo
     enddo
     
     if ( iv .eq. 1 ) then
        
        if ( ig .eq. 0 ) then

           !============================================================
           ! Vector: ig=0
           !============================================================
           do i=1,im2
              al(i,    1) = -al(i+im2,2)
              al(i+im2,1) = -al(i,    2)
           enddo
           
           do i=1,im2
              ar(i,    jm) = -ar(i+im2,jm1)
              ar(i+im2,jm) = -ar(i,    jm1)
           enddo
        else

           !============================================================
           ! ig=1 : SP
           !============================================================
           do i=1,im
              if( i .le. im2) then
                 iop = i+im2
              else
                 iop = i-im2
              endif
              al(i,1) = 0.5*(q(i,1)-q(iop,1)) - r3*(dm(iop,1) + dm(i,1))
           enddo

           !============================================================
           ! NP
           !============================================================
           do i=1,im
              if( i .le. im2) then
                 iop = i+im2
              else
                 iop = i-im2
              endif
              ar(i,jm1) = 0.5*(q(i,jm1)-q(iop,jm1)) - &
             &                 r3*(dm(iop,jm1) + dm(i,jm1))
            enddo
        endif
     else

        ! Scalar (works for ig=0 only):
        do i=1,im2
           al(i,    1) = al(i+im2,2)
           al(i+im2,1) = al(i,    2)
        enddo
        
        do i=1,im2
           ar(i,    jm) = ar(i+im2,jm1)
           ar(i+im2,jm) = ar(i,    jm1)
        enddo
     endif
      
     do j=1,jm-ig
        do i=1,im
           a6(i,j) = 3.*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
        enddo
        call lmppm(dm(1,j), a6(1,j), ar(1,j), al(1,j),  q(1,j), im, jord-3)
     enddo

  END SUBROUTINE Ppm_Lat
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xmap
!
! !DESCRIPTION: Subroutine Xmap performs area preserving mapping in E-W 
!  from an arbitrary resolution to another.  Periodic domain will be assumed, 
!  i.e., the eastern wall bounding cell im is $lon1(im+1) = lon1(1)$; 
!  Note the equal sign is true geophysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Xmap( im, jm, lon1, q1, in, lon2, q2 )
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im           

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: in           

    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm           

    ! Original western edge of the cell
    REAL*4,  INTENT(IN)  :: lon1(im+1)   

    ! Original data at center of the cell
    REAL*4,  INTENT(IN)  :: q1(im,jm)    

    ! Target cell's western edge
    REAL*4,  INTENT(IN)  :: lon2(in+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT) :: q2(in,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: S.-J. Lin
!   First version: piece-wise constant mapping
!   Apr 1, 2000
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*4               :: qtmp(-im:im+im)
    REAL*4               :: al(-im:im+im)
    REAL*4               :: ar(-im:im+im)
    REAL*4               :: a6(-im:im+im)
    REAL*4               :: x1(-im:im+im+1)
    REAL*4               :: dx1(-im:im+im)
    REAL*4               :: pl, pr, qsum, esl, dx
    INTEGER              :: iord = 3
    LOGICAL              :: found
!
! !DEFINED PARAMETERS:
!
    REAL*4,  PARAMETER   :: r3 = 1./3., r23 = 2./3. 

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo

    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo
    
    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'failed in xmap'
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo
    
    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'failed in xmap'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo
    
    !write(6,*) 'i1,i2=',i1,i2

    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       ! Construct subgrid PP distribution
       call ppm_cycle(im, q1(1,j), al(1), ar(1), a6(1), qtmp(0), iord)
       
       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
             al(i) = al(im+i)
             ar(i) = ar(im+i)
             a6(i) = a6(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
             al(i) =   al(i-im)
             ar(i) =   ar(i-im)
             a6(i) =   a6(i-im)
          enddo
       endif
        
       i0 = i1
        
       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             pl = (lon2(i)-x1(m)) / dx1(m)
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                pr = (lon2(i+1)-x1(m)) / dx1(m)
                q2(i,j) = al(m) + 0.5*(a6(m)+ar(m)-al(m)) &
               &                  *(pr+pl)-a6(m)*r3*(pr*(pr+pl)+pl**2)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                qsum = (x1(m+1)-lon2(i))*(al(m)+0.5*(a6(m)+ &
               &              ar(m)-al(m))*(1.+pl)-a6(m)*   &
               &               (r3*(1.+pl*(1.+pl))))

                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)

                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      esl = dx / dx1(mm)
                      qsum = qsum + dx*(al(mm)+0.5*esl* &
                     &              (ar(mm)-al(mm)+a6(mm)*(1.-r23*esl)))
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100       continue
123       q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue

   END SUBROUTINE Xmap
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ppm_Cycle
!
! !DESCRIPTION: Subroutine Ppm\_Cycle is called by Xmap.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Ppm_Cycle( im, q, al, ar, a6, p, iord )
!
! !INPUT PARAMETERS:
!
     INTEGER, INTENT(IN)  :: im, iord
     REAL*4,  INTENT(IN)  :: q(1)
!
! !OUTPUT PARAMETERS:
!
     REAL*4,  INTENT(OUT) :: al(1), ar(1), a6(1), p(0:im+1)
!
! !AUTHOR:
!   Originally written by S-J Lin 
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     REAL*4               :: dm(0:im), tmp, qmax, qmin
     INTEGER              :: i, lmt
!
! !DEFINED PARAMETERS:
! 
     REAL*4,  PARAMETER   :: r3 = 1./3. 

     ! PPM_CYCLE begins here!
     p(0) = q(im)
     do i=1,im
        p(i) = q(i)
     enddo
     p(im+1) = q(1)

     ! 2nd order slope
     do i=1,im
        tmp = 0.25*(p(i+1) - p(i-1))
        qmax = max(p(i-1), p(i), p(i+1)) - p(i)
        qmin = p(i) - min(p(i-1), p(i), p(i+1))
        dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
     enddo
     dm(0) = dm(im)

     do i=1,im
        al(i) = 0.5*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
     enddo

     do i=1,im-1
        ar(i) = al(i+1)
     enddo
     ar(im) = al(1)

     if(iord .le. 6) then
        do i=1,im
           a6(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
        enddo
        lmt = iord - 3
        if(lmt.le.2) call lmppm(dm(1),a6(1),ar(1),al(1),p(1),im,lmt)
     else
        call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
     endif

   END SUBROUTINE Ppm_cycle
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lmppm
!
! !DESCRIPTION: Subroutine Lmppm is called by Ppm\_Cycle.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Lmppm( dm, a6, ar, al, p, im, lmt )
!
! !INPUT PARAMETERS:
!
     ! Longitude dimension
     INTEGER, INTENT(IN) :: im

     ! LMT = 0: full monotonicity
     ! LMT = 1: semi-monotonic constraint (no undershoot)
     ! LMT = 2: positive-definite constraint
     INTEGER, INTENT(IN) :: lmt
!
! !OUTPUT PARAMETERS:
!
     REAL*4              :: a6(im),ar(im),al(im),p(im),dm(im)
!
! !AUTHOR:
!   Originally written by S-J Lin 
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     INTEGER           :: i
     REAL*4            :: da1, da2, fmin, a6da
!
! !DEFINED PARAMETERS:
!
     REAL*4, PARAMETER :: r12 = 1./12. 

     ! LMPPM begins here!
     if(lmt.eq.0) then

        ! Full constraint
        do 100 i=1,im
           if(dm(i) .eq. 0.) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0.
           else
              da1  = ar(i) - al(i)
              da2  = da1**2
              a6da = a6(i)*da1
              if(a6da .lt. -da2) then
                 a6(i) = 3.*(al(i)-p(i))
                 ar(i) = al(i) - a6(i)
              elseif(a6da .gt. da2) then
                 a6(i) = 3.*(ar(i)-p(i))
                 al(i) = ar(i) - a6(i)
              endif
           endif
100     continue

     elseif(lmt.eq.1) then

        ! Semi-monotonic constraint
        do 150 i=1,im
           if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 150
           if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0.
           elseif(ar(i) .gt. al(i)) then
              a6(i) = 3.*(al(i)-p(i))
              ar(i) = al(i) - a6(i)
           else
              a6(i) = 3.*(ar(i)-p(i))
              al(i) = ar(i) - a6(i)
           endif
150     continue
           
     elseif(lmt.eq.2) then

        ! Positive definite constraint
        do 250 i=1,im
           if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 250
           fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
           if(fmin.ge.0.) go to 250
           if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0.
           elseif(ar(i) .gt. al(i)) then
              a6(i) = 3.*(al(i)-p(i))
              ar(i) = al(i) - a6(i)
           else
              a6(i) = 3.*(ar(i)-p(i))
              al(i) = ar(i) - a6(i)
           endif
250     continue
     endif

   END SUBROUTINE Lmppm
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Huynh
!
! !DESCRIPTION: Subroutine Huynh enforces Huynh's 2nd constraint in 
!  the 1D periodic domain.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Huynh( im, ar, al, p, d2, d1 )
!
! !INPUT PARAMETERS:
!
    INTEGER :: im
!
! !OUTPUT PARAMETERS:
!
    REAL*4  :: ar(im), al(im), p(im), d2(im), d1(im)
!
! !AUTHOR:
!   Originally written by S-J Lin 
!
! !REVISION HISTORY: 
!   21 Sep 2000 - R. Yantosca - Converted to F90 freeform format and inserted
!                               into "Geos3RegridModule".  Added F90 type 
!                               declarations to be consistent with  
!                               TypeModule.f90.  Also updated comments. 
!   08 Nov 2006 - R. Yantosca - Inserted into Geos57RegridModule.f90
!   12 Dec 2008 - R. Yantosca - Added ProTeX documentation headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i
    REAL*4  :: pmp, lac, pmin, pmax

    !===================================================================
    ! HUYNH begins here!
    ! Compute d1 and d2
    !===================================================================
    d1(1) = p(1) - p(im)
    do i=2,im
       d1(i) = p(i) - p(i-1)
    enddo
    
    do i=1,im-1
       d2(i) = d1(i+1) - d1(i)
    enddo
    d2(im) = d1(1) - d1(im)

    !===================================================================
    ! Constraint for AR
    ! i = 1
    !===================================================================
    pmp   = p(1) + 2.0 * d1(1)
    lac   = p(1) + 0.5 * (d1(1)+d2(im)) + d2(im) 
    pmin  = min(p(1), pmp, lac)
    pmax  = max(p(1), pmp, lac)
    ar(1) = min(pmax, max(ar(1), pmin))
    
    do i=2, im
       pmp   = p(i) + 2.0*d1(i)
       lac   = p(i) + 0.5*(d1(i)+d2(i-1)) + d2(i-1)
       pmin  = min(p(i), pmp, lac)
       pmax  = max(p(i), pmp, lac)
       ar(i) = min(pmax, max(ar(i), pmin))
    enddo
     
    !==================================================================
    ! Constraint for AL
    !==================================================================
    do i=1, im-1
       pmp   = p(i) - 2.0*d1(i+1)
       lac   = p(i) + 0.5*(d2(i+1)-d1(i+1)) + d2(i+1)
       pmin  = min(p(i), pmp, lac)
       pmax  = max(p(i), pmp, lac)
       al(i) = min(pmax, max(al(i), pmin))
    enddo

    !==================================================================
    ! i=im
    !==================================================================
    i = im
    pmp    = p(im) - 2.0*d1(1)
    lac    = p(im) + 0.5*(d2(1)-d1(1)) + d2(1)
    pmin   = min(p(im), pmp, lac)
    pmax   = max(p(im), pmp, lac)
    al(im) = min(pmax, max(al(im), pmin))

    !==================================================================
    ! compute A6 (d2)
    !==================================================================
    do i=1, im
       d2(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
    enddo

  END SUBROUTINE Huynh
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Geos57RegridInit
!
! !DESCRIPTION: Subroutine Geos57RegridInit initializes the longitude and 
!  latitude edge arrays for 0.5 x 0.666, 1 x 1.25, 2 x 2.5, and 4 x 5 grids.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57RegridInit
!
! !REMARKS:
!  Computation is done in REAL*8 and then casted to REAL*4 in order
!  to get correct values for the high-resolution grids. 
!
! !REVISION HISTORY: 
!  25 Oct 2011 - R. Yantosca - InitialVersion, based on MerraRegridModule
!  26 Oct 2011 - R. Yantosca - Now initialize lon & lat center arrays
!  03 Jan 2012 - R. Yantosca - Now define index arrays for vertical levels
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L
    REAL*8  :: DI, DJ

    !======================================================================
    ! GEOS-5.7.x NATIVE GRID
    ! 0.25 x 0.3125 resolution; centered on (-180,-90)
    !======================================================================

    ! Size of box
    DI = 0.3125d0
    DJ = 0.25d0

    ! Lon edges
    DO I = 0, I025x03125
       xedge_025x03125(I+1)       = -180d0 - DI/2d0 + ( DI * I )
    ENDDO

    ! Lat edges
    DO J = 0, J025x03125 
       yedge_025x03125(J+1)       =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO

    ! Lon centers
    DO I = 0, I025x03125-1
       xmid_025x03125(I+1)        = -180d0          + ( DI * I )
    ENDDO                         
                                  
    ! Lat centers                 
    DO J = 0, J025x03125-1              
       ymid_025x03125(J+1)        =  -90d0          + ( DJ * J ) 
    ENDDO  

    ! Reset poles
    yedge_025x03125(1           ) = -90e0
    yedge_025x03125(J025x03125+1) = +90e0
    ymid_025x03125 (1           ) = -89.9375d0
    ymid_025x03125 (J025x03125  ) = +89.9375d0

    ! Sine of latitude edges
    DO J = 1, J025x03125+1
       sine_025x03125(J)          = SIN( yedge_025x03125(J) * D2R )
    ENDDO

    ! Vertical levels
    DO L = 1, L025x03125
       zmid_025x03125(L) = L
    ENDDO

    !======================================================================
    ! 0.5 x 0.666 resolution; centered on (-180,-90)
    !======================================================================

    ! Size of box
    DI = 2d0 / 3d0
    DJ = 0.5d0

    ! Lon edges
    DO I = 0, I05x0666
       xedge_05x0666(I+1)         = -180d0 - DI/2d0 + ( DI * I )
    ENDDO                         
                                  
    ! Lat edges                   
    DO J = 0, J05x0666            
       yedge_05x0666(J+1)         =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO                         
        
    ! Lon centers
    DO I = 0, I05x0666-1
       xmid_05x0666(I+1)          = -180d0          + ( DI * I )
    ENDDO                         
                                  
    ! Lat centers                 
    DO J = 0, J05x0666-1              
       ymid_05x0666(J+1)          =  -90d0          + ( DJ * J ) 
    ENDDO  
                          
    ! Reset poles                 
    yedge_05x0666(1         )     = -90e0
    yedge_05x0666(J05x0666+1)     = +90e0
    ymid_05x0666 (1         )     = -89.875d0
    ymid_05x0666 (J05x0666  )     = +89.87500
                                  
    ! Sine of latitude edges      
    DO J = 1, J05x0666+1          
       sine_05x0666(J)            = SIN( yedge_05x0666(J) * D2R )
    ENDDO

    ! Vertical levels
    DO L = 1, L05x0666
       zmid_05x0666(L) = L
    ENDDO

    !======================================================================
    ! 1 x 1.25 resolution; centered on (-180,-90)
    !======================================================================
    
    ! Size of box
    DI = 1.25d0
    DJ = 1.00d0

    ! Lon edges
    DO I = 0, I1x125
       xedge_1x125(I+1)           = -180d0 - DI/2d0 + ( DI * I )
    ENDDO                         
                                  
    ! Lat edges                   
    DO J = 0, J1x125              
       yedge_1x125(J+1)           =  -90d0 - DJ/2d0 + ( DJ * J ) 
    ENDDO   

    ! Lon centers
    DO I = 0, I1x125-1
       xmid_1x125(I+1)            = -180d0          + ( DI * I )
    ENDDO                         
                                  
    ! Lat centers                 
    DO J = 0, J1x125-1              
       ymid_1x125(J+1)            =  -90d0          + ( DJ * J ) 
    ENDDO                         
                                  
    ! Reset poles                 
    yedge_1x125(1       )         = -90e0
    yedge_1x125(J1x125+1)         = +90e0
    ymid_1x125 (1       )         = -89.75e0
    ymid_1x125 (J1x125  )         = +89.75e0
                                  
    ! Sine of latitude edges      
    DO J = 1, J1x125+1            
       sine_1x125(J)              = SIN( yedge_1x125(J) * D2R )
    ENDDO

    ! Vertical levels
    DO L = 1, L1x125
       zmid_1x125(L) = L
    ENDDO

    !======================================================================
    ! 2 x 2.5 resolution; centered on (-180,-90)
    !======================================================================

    ! Size of box
    DI = 2.5d0
    DJ = 2.0d0

    ! Lon edges
    DO I = 0, I2x25
       xedge_2x25(I+1)            = -180d0 - DI/2d0 + ( DI * I )
    ENDDO                         
                                  
    ! Lat edges                   
    DO J = 0, J2x25               
       yedge_2x25(J+1)            =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO                         
            
    ! Lon centers
    DO I = 0, I2x25-1
       xmid_2x25(I+1)             = -180d0          + ( DI * I )
    ENDDO
   
    ! Lat centers         
    DO J = 0, J2x25-1              
       ymid_2x25(J+1)             =  -90d0          + ( DJ * J )
    ENDDO 
         
    ! Reset poles                 
    yedge_2x25(1      )           = -90e0
    yedge_2x25(J2x25+1)           = +90e0
    ymid_2x25 (1      )           = -89.75e0
    ymid_2x25 (J2x25  )           = +89.75e0
                                  
    ! Sine of latitude edges      
    DO J = 1, J2x25+1             
       sine_2x25(J)               = SIN( yedge_2x25(J) * D2R )
    ENDDO

    ! Vertical levels
    DO L = 1, L2x25
       zmid_2x25(L) = L
    ENDDO

    !======================================================================
    ! 4 x 5 Grid; centered on (-180,-90)
    !======================================================================

    ! Size of box
    DI = 5d0
    DJ = 4d0

    ! Lon edges
    DO I = 0, I4x5
       xedge_4x5(I+1)            = -180d0 - DI/2d0 + ( DI * I )
    ENDDO                        
                         
    ! Lat edges         
    DO J = 0, J4x5               
       yedge_4x5(J+1)            =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO                        

    ! Lon centers
    DO I = 0, I4x5-1
       xmid_4x5(I+1)             = -180d0          + ( DI * I )
    ENDDO
   
    ! Lat centers         
    DO J = 0, J4x5-1              
       ymid_4x5(J+1)             =  -90d0          + ( DJ * J )
    ENDDO 
              
    ! Reset poles                
    yedge_4x5(1     )            = -90e0   
    yedge_4x5(J4x5+1)            = +90e0
    ymid_4x5 (1     )            = -89e0
    ymid_4x5 (J4x5  )            = +89e0
                                 
    ! Sine of latitude edges     
    DO J = 1, J4x5+1             
       sine_4x5(J)               = SIN( yedge_4x5(J) * D2R )
    ENDDO

    ! Vertical levels
    DO L = 1, L4x5
       zmid_4x5(L) = L
    ENDDO

  END SUBROUTINE Geos57RegridInit
!EOC
END MODULE Geos57RegridModule

