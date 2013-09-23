;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_cn
;
; !DESCRIPTION: Routine to adjust the existing GEOS_FP "CN" (constant)
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
 procedure adjust_cn( inFile   : string, 
                      apBpFile : string, 
                      outFile  : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-5.7.2 "CN" file w/ input data
;  apBpFile ; Name of file w/ Ap and Bp data
;  outFile  : Name of GEOSFP "CN" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn,     time, lat, lon, FRLAKE, FRLAND, FRLANDIC, 
       FROCEAN, PHIS, ap,  bp,  lev,    fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Also merges the ap and bp fields into the CN file.  We need to take 
;       these out of the 3-D files because ESMF/MAPL cannot read netCDF 
;       files in which there are variables with more than one vertical 
;       dimension.
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
  lat                    = fIn->lat
  lon                    = fIn->lon
  FRLAKE                 = fIn->FRLAKE
  FRLAND                 = fIn->FRLAND
  FRLANDIC               = fIn->FRLANDIC
  FROCEAN                = fIn->FROCEAN
  PHIS                   = fIn->PHIS

  ;=========================================================================
  ; Read the AP and BP fields
  ;=========================================================================
  fApBp                  = addfile( apBpFile, "r" )
  ap                     = fApBp->ap
  bp                     = fApBp->bp

  ; Rename the dimensions of ap and bp to "lev"
  ap!0                   = "lev"
  bp!0                   = "lev"

  ; Create a new index variable for levels
  ; Redefine "lev" for all relevant variables
  lev                    = fspan( 1.0, 73.0, 73 )
  lev!0                  = "lev"
  lev@long_name          = "levels"
  lev@units              = "1"
  lev&lev                = lev
  ap&lev                 = lev
  bp&lev                 = lev

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
  bp@units               = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( FRLAKE   )
  add_coards_var_atts( FRLAND   )
  add_coards_var_atts( FRLANDIC )
  add_coards_var_atts( FROCEAN  )
  add_coards_var_atts( PHIS     )
  add_coards_var_atts( ap       )
  add_coards_var_atts( bp       )

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
  ;printVarSummary( ap       )
  ;printVarSummary( bp       )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->FRLAKE           = FRLAKE
  fOut->FRLAND           = FRLAND
  fOut->FRLANDIC         = FRLANDIC
  fOut->FROCEAN          = FROCEAN
  fOut->PHIS             = PHIS
  fOut->ap               = ap
  fOut->bp               = bp

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
end
;EOC