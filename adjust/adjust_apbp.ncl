;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !MODULE: adjust_apbp
;
; !DESCRIPTION: Creates a COARDS-compliant file for the GEOS-FP Ap and Bp
;  hybrid-grid parameters.
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
 undef( "adjust_apbp" )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-5.7.2 "Ap/Bp" file w/ input data
;  outFile  : Name of GEOSFP "Ap/Bp" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;  local fIn, time, lat, lon, lev, ap,  bp, fOut
;
; !REMARKS:
;  (1) You only have to run this file once, as the GEOS-FP Ap and Bp 
;       parameters are the same for all horizontal resolutions. 
;  (2) Looks for the input file GEOS572.Input.ApBp.nc in this directory.
;
; !REVISION HISTORY:
;  24 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Read the AP and BP fields
  ;=========================================================================
  inFile                   = "./GEOS572.Input.ApBp.nc"
  fIn                      = addfile( inFile, "r" )
  ap                       = fIn->ap
  bp                       = fIn->bp

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; Rename the dimensions of ap and bp to "lev"
  ap!0                     = "lev"
  bp!0                     = "lev"

  ; Create a new index variable for levels
  ; Redefine "lev" for all relevant variables
  lev                     = fspan( 1.0, 73.0, 73 )
  lev!0                   = "lev"
  lev@long_name           = "levels"
  lev@units               = "1"
  lev&lev                 = lev
  ap&lev                  = lev
  bp&lev                  = lev

  ; Unitless quantities should have units = "1"
  bp@units                 = "1"

  ; Add common COARDS variable attributes
  add_coards_var_atts( ap )
  add_coards_var_atts( bp )

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( lev )
  ;printVarSummary( ap  )
  ;printVarSummary( bp  )
  ;exit

  ; Open output file (remove prior version)
  outFile                  = "./GEOSFP.ApBp.nc"
  system( "rm -f " + outFile )
  fOut                     = addfile( outFile, "c" )
  fOut->lev                = lev
  fOut->ap                 = ap
  fOut->bp                 = bp

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
  ; Redefine some global parameters.  
  ; Ap and Bp do not depend on time, lat, or lon.
  date                      = systemfunc( "date" )
  fOut@Title                = "GEOS-FP hybrid grid Ap and Bp parameters"
  fOut@Filename             = str_sub_str( outFile, "./", "" )
  fOut@History              = "File created on: " + date
  fOut@ProductionDateTime   = "File created on: " + date
  fOut@ModificationDateTime = "File created on: " + date
  delete( fOut@Start_Date )
  delete( fOut@Start_Time )
  delete( fOut@End_Date   )
  delete( fOut@End_Time   )
  delete( fOut@Delta_Time )
  delete( fOut@Delta_Lon  )
  delete( fOut@Delta_Lat  )
end
;EOC