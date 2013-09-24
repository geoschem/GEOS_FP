;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_a3mste
;
; !DESCRIPTION: Routine to adjust the existing GEOS-FP "A3mstE" (moist 
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
 undef( "adjust_a3mste" )
 procedure adjust_a3mste( inFile  : string, 
                          outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-5.7.2 "A3mstE" file w/ input data
;  outFile  : Name of GEOSFP "A3mstE" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fileDyn, fDyn,  time,    lev,   lat,     lon, 
       CMFMC,   PFICU, PFILSAN, PFLCU, PFLLSAN, fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Ensures that the poles have latitudes -90 and +90 degrees.
;  (3) Moves the CMFMC field from A3dyn to A3mstE, so that we can keep
;       all variables on level edges in the same file.  This is necessary
;       when reading netCDF files with ESMF/MAPL.
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Read data from the A3mstE input file
  ;=========================================================================
  fIn                    = addfile( inFile, "r" )
  time                   = fIn->time
  lev                    = fIn->lev
  lat                    = fIn->lat
  lon                    = fIn->lon
  PFICU                  = fIn->PFICU
  PFILSAN                = fIn->PFILSAN
  PFLCU                  = fIn->PFLCU
  PFLLSAN                = fIn->PFLLSAN

  ;=========================================================================
  ; Read CMFMC data from the A3dyn file
  ; Change vertical dimension from "ap" to "lev"
  ;=========================================================================
  fileDyn                = str_sub_str( inFile, "A3mstE", "A3dyn" )
  fDyn                   = addfile( fileDyn, "r" )
  CMFMC                  = fDyn->CMFMC
  CMFMC!1                = "lev"
  CMFMC&lev              = lev

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
    CMFMC&lat            = lat 
    PFICU&lat            = lat
    PFILSAN&lat          = lat
    PFLCU&lat            = lat
    PFLLSAN&lat          = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  CMFMC&time             = time 
  PFICU&time             = time
  PFILSAN&time           = time
  PFLCU&time             = time
  PFLLSAN&time           = time

  ; The units of the lev index array has to be 1
  ; Redefine lev for all relevant variables
  lev@units              = "1"
  lev&lev                = lev
  CMFMC&lev              = lev 
  PFICU&lev              = lev
  PFILSAN&lev            = lev
  PFLCU&lev              = lev
  PFLLSAN&lev            = lev

  ; Add common COARDS variable attributes
  add_coards_var_atts( CMFMC   )
  add_coards_var_atts( PFICU   )
  add_coards_var_atts( PFILSAN )
  add_coards_var_atts( PFLCU   )
  add_coards_var_atts( PFLLSAN )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lev      )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( CMFMC    )
  ;printVarSummary( PFICU    )
  ;printVarSummary( PFILSAN  )
  ;printVarSummary( PFLCU    )
  ;printVarSummary( PFLLSAN  )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->CMFMC            = CMFMC
  fOut->PFICU            = PFICU
  fOut->PFILSAN          = PFILSAN
  fOut->PFLCU            = PFLCU
  fOut->PFLLSAN          = PFLLSAN

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )

  ; Redefine the title string
  fOut@Title             = "GEOS-FP time-averaged 3-hour moist parameters on level edges (A3mstE), processed for GEOS-Chem input"
end
;EOC