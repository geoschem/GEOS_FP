;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_cn
;
; !DESCRIPTION: Routine to adjust the existing GEOS-FP "CN" (constant)
;  data files for better COARDS compliance.
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
 undef( "adjust_cn" )
 procedure adjust_cn( inFile  : string, 
                      outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-FP "CN" file w/ input data
;  outFile  : Name of GEOS-FP "CN" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn,     time, lat, lon, FRLAKE, FRLAND, FRLANDIC, 
       FROCEAN, PHIS, ap,  bp,  lev,    fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Ensures that the poles have latitudes -90 and +90 degrees.
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;  24 Sep 2013 - R. Yantosca - Split off Ap, Bp handling to adjust_apbp.ncl
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Read data from the input file
  ;=========================================================================
  fIn                    = addfile( inFile, "r" )
  time                   = fIn->time
  lat                    = fIn->lat
  lon                    = fIn->lon
  FRLAKE                 = fIn->FRLAKE
  FRLAND                 = fIn->FRLAND
  FRLANDIC               = fIn->FRLANDIC
  FROCEAN                = fIn->FROCEAN
  PHIS                   = fIn->PHIS

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
    FRLAKE&lat           = lat
    FRLAND&lat           = lat
    FRLANDIC&lat         = lat
    FROCEAN&lat          = lat
    PHIS&lat             = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  FRLAKE&time            = time
  FRLAND&time            = time
  FRLANDIC&time          = time
  FROCEAN&time           = time
  PHIS&time              = time

  ; Unitless quantities should have units = "1"
  FRLAKE@units           = "1"
  FRLAND@units           = "1"
  FRLANDIC@units         = "1"
  FROCEAN@units          = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( FRLAKE   )
  add_coards_var_atts( FRLAND   )
  add_coards_var_atts( FRLANDIC )
  add_coards_var_atts( FROCEAN  )
  add_coards_var_atts( PHIS     )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lev      )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( FRLAKE   )
  ;printVarSummary( FRLAND   )
  ;printVarSummary( FRLANDIC )
  ;printVarSummary( FROCEAN  )
  ;printVarSummary( PHIS     )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->FRLAKE           = FRLAKE
  fOut->FRLAND           = FRLAND
  fOut->FRLANDIC         = FRLANDIC
  fOut->FROCEAN          = FROCEAN
  fOut->PHIS             = PHIS

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
 
  ; Redefine the title string
  fOut@Title             = "GEOS-FP constant parameters (CN), processed for GEOS-Chem input"

end
;EOC