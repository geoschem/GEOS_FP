;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_a3cld
;
; !DESCRIPTION: Routine to adjust the existing GEOS-FP "A3cld" (cloud 
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
 undef( "adjust_a3cld" )
 procedure adjust_a3cld( inFile  : string, 
                         outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-FP "A3cld" file w/ input data
;  outFile  : Name of GEOS-FP "A3cld" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn,      time, lev, lat, lon,    CLOUD, 
       OPTDEPTH, QCCU, QI,  QL,  TAUCLI, TAUCLW, fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Ensures that the poles have latitudes -90 and +90 degrees.
;  (3) Omits saving out CFAN, CFCU.  We really only need CLOUD.
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
  CLOUD                  = fIn->CLOUD
  OPTDEPTH               = fIn->OPTDEPTH
  QCCU                   = fIn->QCCU
  QI                     = fIn->QI
  QL                     = fIn->QL
  TAUCLI                 = fIn->TAUCLI
  TAUCLW                 = fIn->TAUCLW

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
    CLOUD&lat            = lat
    OPTDEPTH&lat         = lat
    QCCU&lat             = lat
    QI&lat               = lat
    QL&lat               = lat
    TAUCLI&lat           = lat
    TAUCLW&lat           = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  CLOUD&time             = time
  OPTDEPTH&time          = time
  QCCU&time              = time
  QI&time                = time
  QL&time                = time
  TAUCLI&time            = time
  TAUCLW&time            = time

  ; The units of the lev index array has to be 1
  ; Redefine lev for all relevant variables
  lev@units              = "1"
  lev&lev                = lev
  CLOUD&lev              = lev
  OPTDEPTH&lev           = lev
  QCCU&lev               = lev
  QI&lev                 = lev
  QL&lev                 = lev
  TAUCLI&lev             = lev
  TAUCLW&lev             = lev

  ; Quantities that are unitless need to have a unit string of "1"
  CLOUD@units            = "1"
  OPTDEPTH@units         = "1"
  TAUCLI@units           = "1"
  TAUCLW@units           = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( CLOUD    )
  add_coards_var_atts( OPTDEPTH )
  add_coards_var_atts( QCCU     )
  add_coards_var_atts( QI       )
  add_coards_var_atts( QL       )
  add_coards_var_atts( TAUCLI   )
  add_coards_var_atts( TAUCLW   )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lev      )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( CLOUD    )
  ;printVarSummary( OPTDEPTH )
  ;printVarSummary( QCCU     )
  ;printVarSummary( QI       )
  ;printVarSummary( QL       )
  ;printVarSummary( TAUCLI   )
  ;printVarSummary( TAUCLW   )
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lev              = lev
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->CLOUD            = CLOUD
  fOut->OPTDEPTH         = OPTDEPTH
  fOut->QCCU             = QCCU
  fOut->QI               = QI
  fOut->QL               = QL
  fOut->TAUCLI           = TAUCLI
  fOut->TAUCLW           = TAUCLW

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
  ; Redefine the title string
  fOut@Title             = "GEOS-FP time-averaged 3-hour cloud parameters (A3cld), processed for GEOS-Chem input"

end
;EOC