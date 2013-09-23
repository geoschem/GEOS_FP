;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_all
;
; !DESCRIPTION: Driver routine for the GEOS-FP adjustment code.  Calls 
;  other scripts to adjust GEOS-FP data files for better COARDS compliance.
;\\
;\\
; !USES:
;
 ; NCL routines
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

 ; Local routines
 load "./adjust_a3cld.ncl"
 load "./adjust_a3dyn.ncl"
 load "./adjust_cn.ncl"
 load "./adjust_i3.ncl"
;
; !INTERFACE:
;
 undef( "adjust_all" )
;
; !REMARKS:
;  Test code -- will make more general later.  Need to loop over all
;  files of a given type in a directory.
;
; !REVISION HISTORY:
;  11 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  dir      = "/as/tmp/all/tmp/bmy/GEOS_FP/adjust/"

  ;=========================================================================
  ; Adjust the CN file
  ;=========================================================================
  inFile   = dir + "GEOS572.20110101.CN.4x5.nc"
  apBpFile = dir + "GEOS572.20110101.ApBp.nc"
  outFile  = dir + "GEOSFP.20110101.CN.4x5.nc"
  adjust_cn( inFile, apBpFile, outFile )

  ;=========================================================================
  ; Adjust the I3 file
  ;=========================================================================
  inFile   = dir + "GEOS572.20130921.I3.4x5.nc"
  outFile  = dir + "GEOSFP.20130921.I3.4x5.nc"
  adjust_i3( inFile, outFile )

  ;=========================================================================
  ; Adjust the A3cld file
  ;=========================================================================
  inFile   = dir + "GEOS572.20130921.A3cld.4x5.nc"
  outFile  = dir + "GEOSFP.20130921.A3cld.4x5.nc"
  adjust_a3cld( inFile, outFile )

  ;=========================================================================
  ; Adjust the A3dyn file
  ;=========================================================================
  inFile   = dir + "GEOS572.20130921.A3dyn.4x5.nc"
  outFile  = dir + "GEOSFP.20130921.A3dyn.4x5.nc"
  adjust_a3dyn( inFile, outFile )
end
;EOC