;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_cn
;
; !DESCRIPTION: Routine to adjust the existing GEOS_FP data files for 
;  COARDS compliance.  Also merges the ap and bp fields into the CN file.
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
 local fIn,     time, lat, lon, frLake, frLand, frLandIc, 
       frOcean, phis, ap,  bp,  lev,    fOut
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
  frLake                 = fIn->FRLAKE
  frLand                 = fIn->FRLAND
  frLandIc               = fIn->FRLANDIC
  frOcean                = fIn->FROCEAN
  phis                   = fIn->PHIS

  ;=========================================================================
  ; Read the AP and BP fields
  ;=========================================================================
  fApBp                  = addfile( apBpFile, "r" )
  ap                     = fApBp->ap
  bp                     = fApBp->bp

  ; Create a new index variable for levels
  lev                    = fspan( 1.0, 73.0, 73 )
  lev!0                  = "lev"
  lev@long_name          = "levels"
  lev@units              = "1"
  lev&lev                = lev

  ; Rename the dimensions of ap and bp to "lev"
  ap!0                   = "lev"
  bp!0                   = "lev"
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
    frLake&lat           = lat
    frLand&lat           = lat
    frLandIc&lat         = lat
    frOcean&lat          = lat
    phis&lat             = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  time@calendar          = "gregorian"
  time&time              = time

  ; Unitless quantities should have units = "1"
  frLake@units           = "1"
  frLand@units           = "1"
  frLandIc@units         = "1"
  frOcean@units          = "1"
  bp@units               = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( frLake   )
  add_coards_var_atts( frLand   )
  add_coards_var_atts( frLandIc )
  add_coards_var_atts( frOcean  )
  add_coards_var_atts( phis     )
  add_coards_var_atts( ap       )
  add_coards_var_atts( bp       )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->FRLAKE           = frLake
  fOut->FRLAND           = frLand
  fOut->FRLANDIC         = frLandIc
  fOut->FROCEAN          = frOcean
  fOut->PHIS             = phis
  fOut->ap               = ap
  fOut->bp               = bp

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lev      )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( frLake   )
  ;printVarSummary( frLand   )
  ;printVarSummary( frLandIc )
  ;printVarSummary( frOcean  )
  ;printVarSummary( phis     )
  ;printVarSummary( ap       )
  ;printVarSummary( bp       )
  
end
;EOC