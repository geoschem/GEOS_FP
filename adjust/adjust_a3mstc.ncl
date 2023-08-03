;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_a3mstc
;
; !DESCRIPTION: Routine to adjust the existing GEOS-FP "A3mstC" (moist
;  variables, 3-hour time averaged) data files for better COARDS compliance.
;\\
;\\
; !USES:
;
 ; NCL routines
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

 ; Local routines
 load "./add_coards_var_atts.ncl"
 load "./add_coards_global_atts.ncl"
;
; !INTERFACE:
;
 undef( "adjust_a3mstc" )
 procedure adjust_a3mstc( inFile  : string, 
                          outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-FP "A3mstC" file w/ input data
;  outFile  : Name of GEOS-FP "A3mstC" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn, time, lev, lat, lon, DQRCU, DQRLSAN, REEVAPCN, REEVAPLS, fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Ensures that the poles have latitudes -90 and +90 degrees.
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Read data from the input file
  ;=========================================================================
  fIn                    = addfile( inFile, "r" )
  time                   = fIn->time
  lev                    = fIn->lev
  lat                    = fIn->lat
  lon                    = fIn->lon
  DQRCU                  = fIn->DQRCU
  DQRLSAN                = fIn->DQRLSAN
  REEVAPCN               = fIn->REEVAPCN
  REEVAPLS               = fIn->REEVAPLS

  ;=========================================================================
  ; For ESMF/MAPL, we need to make sure that the latitudes at the poles
  ; are forced to -90/+90 degrees.  This is a requirement of MAPL.
  ; We don't have to worry about this for nested grids.
  ;=========================================================================
  if ( ( fIn@delta_lat .eq. 2 ) .or. ( fIn@delta_lat .eq. 4 ) ) then

    ; Force the poles to be -90/+90
    lat(0)               = -90.0
    lat(dimsizes(lat)-1) =  90.0  

    ; Redefine the latitude dimension for all relevant variables
    lat&lat              = lat
    DQRCU&lat            = lat
    DQRLSAN&lat          = lat
    REEVAPCN&lat         = lat
    REEVAPLS&lat         = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  DQRCU&time             = time
  DQRLSAN&time           = time
  REEVAPCN&time          = time
  REEVAPLS&time          = time

  ; The units of the lev index array has to be 1
  ; Redefine lev for all relevant variables
  lev@units              = "1"
  lev&lev                = lev
  DQRCU&lev              = lev
  DQRLSAN&lev            = lev
  REEVAPCN&lev           = lev
  REEVAPLS&lev           = lev

  ; Add common COARDS variable attributes
  add_coards_var_atts( DQRCU    )
  add_coards_var_atts( DQRLSAN  )
  add_coards_var_atts( REEVAPCN )
  add_coards_var_atts( REEVAPLS )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lev      )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( DQRCU    )
  ;printVarSummary( DQRLSAN  )
  ;printVarSummary( REEVAPCN )
  ;printVarSummary( REEVAPLS )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->DQRCU            = DQRCU
  fOut->DQRLSAN          = DQRLSAN
  fOut->REEVAPCN         = REEVAPCN
  fOut->REEVAPLS         = REEVAPLS

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
  ; Redefine the title string
  fOut@Title             = "GEOS-FP time-averaged 3-hour moist parameters on level centers (A3mstC), processed for GEOS-Chem input"

end
;EOC