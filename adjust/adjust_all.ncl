;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_all
;
; !DESCRIPTION: NCL template script
;\\
;\\
; !USES:
;
 ; NCL code
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

 ; Our code
 load "./adjust_cn.ncl"
;
; !INTERFACE:
;
 undef( "adjust_all" )
;
; !LOCAL VARIABLES:
;
; local
;
; !REMARKS:
;  Test code
;
; !REVISION HISTORY:
;  11 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  dir      = "/as/tmp/all/tmp/bmy/GEOS_FP/adjust/"
  inFile   = dir + "GEOS572.20110101.CN.4x5.nc"
  apBpFile = dir + "GEOS572.20110101.ApBp.nc"
  outFile  = dir + "GEOSFP.20110101.CN.4x5.nc"

  tmp = notdir( outFile )

  ; Adjust the CN file
  system( "rm -f " + outFile )
  adjust_cn( inFile, apBpFile, outFile )

end
;EOC