;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: notdir
;
; !DESCRIPTION: Given a full directory path, returns just the filename part.
;  This is the analog to the GNU Make "notdir" function.
;\\
;\\
; !USES:
;
 ; NCL code
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;
; !INTERFACE:
;
 undef( "notdir" )
 function notdir( path:string )
;
; !LOCAL VARIABLES:
;
 local result
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ; Split by the Unix path separator "/"
  result = str_split( path, "/" )

  ; Return the last element of the RESULT array
  return( result( dimsizes( result ) - 1  ) )

end
;EOC