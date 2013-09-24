;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_a3dyn
;
; !DESCRIPTION: Routine to adjust the existing GEOS_FP "A3dyn" (dynamics,
;  3-hour time averaged) data files for better COARDS compliance.
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
 undef( "adjust_a3dyn" )
 procedure adjust_a3dyn( inFile  : string, 
                         outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-5.7.2 "A3dyn" file w/ input data
;  outFile  : Name of GEOSFP "A3dyn" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn, time, lev, lat, lon, DTRAIN, OMEGA, RH, U, V, fOut
;
; !REMARKS:
;  Also add a calendar attribute onto the "time" index variable.
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
  DTRAIN                 = fIn->DTRAIN
  OMEGA                  = fIn->OMEGA
  RH                     = fIn->RH
  U                      = fIn->U
  V                      = fIn->V

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
    DTRAIN&lat           = lat
    OMEGA&lat            = lat
    RH&lat               = lat
    U&lat                = lat
    V&lat                = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  DTRAIN&time            = time
  OMEGA&time             = time
  RH&time                = time
  U&time                 = time
  V&time                 = time

  ; The units of the lev index array has to be 1
  ; Redefine lev for all relevant variables
  lev@units              = "1"
  lev&lev                = lev
  DTRAIN&lev             = lev
  OMEGA&lev              = lev
  RH&lev                 = lev
  U&lev                  = lev
  V&lev                  = lev

  ; Quantities that are unitless need to have a unit string of "1"
  RH@units               = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( DTRAIN )
  add_coards_var_atts( OMEGA  )
  add_coards_var_atts( RH     )
  add_coards_var_atts( U      )
  add_coards_var_atts( V      )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time   )
  ;printVarSummary( lev    )
  ;printVarSummary( lon    )
  ;printVarSummary( lat    )
  ;printVarSummary( DTRAIN )
  ;printVarSummary( OMEGA  )
  ;printVarSummary( RH     )
  ;printVarSummary( U      )
  ;printVarSummary( V      )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->DTRAIN           = DTRAIN
  fOut->OMEGA            = OMEGA
  fOut->RH               = RH
  fOut->U                = U
  fOut->V                = V

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )

  ; Redefine the title string
  fOut@Title             = "GEOS-FP time-averaged 3-hour dynamical parameters (A3dyn), processed for GEOS-Chem input"

end
;EOC