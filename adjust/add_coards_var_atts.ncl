;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: add_coards_var_atts
;
; !DESCRIPTION: Simple routine to add several common attributes for 
;  COARDS compliance to a variable.
;\\
;\\
; !USES:
;
 ; NCL code
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;
; !INTERFACE:
;
 undef( "add_coards_var_atts" )
 procedure add_coards_var_atts( var )
;
; !LOCAL VARIABLES:
;  none
;
; !REMARKS:
;  Adds the following attributes for better COARDS compliance:
;  (1) missing_value
;  (2) _FillValue
;  (3) scale_factor
;  (4) add_offset
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ; Add common attributes for COARDS compliance
  var@missing_value = 1e+15
  var@_FillValue    = 1e+15
  var@scale_factor  = 1.0
  var@add_offset    = 0.0

end
;EOC