;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_i3
;
; !DESCRIPTION: Routine to adjust the existing GEOS_FP "I3" (instantaneous
;  3-hour) data files for better COARDS compliance.
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
 undef( "adjust_i3" )
 procedure adjust_i3( inFile  : string, 
                      outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-5.7.2 "CN" file w/ input data
;  apBpFile ; Name of file w/ Ap and Bp data
;  outFile  : Name of GEOSFP "CN" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn, time, lev, lat, lon, PS, PV, QV, T, fOut
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
  PS                     = fIn->PS
  PV                     = fIn->PV
  QV                     = fIn->QV
  T                      = fIn->T

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
    PS&lat               = lat
    PV&lat               = lat
    QV&lat               = lat
    T&lat                = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine "time" for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  PS&time                = time
  PV&time                = time
  QV&time                = time
  T&time                 = time

  ; The units of the lev index array has to be 1
  ; Redefine "lev" for all relevant variables
  lev@units              = "1"
  lev&lev                = lev
  PV&lev                 = lev
  QV&lev                 = lev
  T&lev                  = lev

  ; Add common COARDS variable attributes
  add_coards_var_atts( PS )
  add_coards_var_atts( PV )
  add_coards_var_atts( QV )
  add_coards_var_atts( T  )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time )
  ;printVarSummary( lev  )
  ;printVarSummary( lon  )
  ;printVarSummary( lat  )
  ;printVarSummary( PS   )
  ;printVarSummary( PV   )
  ;printVarSummary( QV   )
  ;printVarSummary( T    )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->PS               = PS
  fOut->PV               = PV
  fOut->QV               = QV
  fOut->T                = T

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
end
;EOC