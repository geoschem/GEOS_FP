;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: add_coards_global_atts
;
; !DESCRIPTION: Simple routine to (1) copy global attribute values from
;  an input file to an output file, and (2) rename the global attributes
;  for COARDS compliance.
;\\
;\\
; !USES:
;
 ; NCL routines
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

 ; Local routines
 load "./notdir.ncl"
;
; !INTERFACE:
;
 undef( "add_coards_global_atts" )
 procedure add_coards_global_atts( fIn  : file, 
                                   fOut : file  )
;
; !INPUT PARAMETERS:
;  inFile  : File handle for input  file (opened with addfile)
;  oufFile : File handle for output file (opened with addfile)
;
; !LOCAL VARIABLES:
;
 local fileName
;
; !REMARKS:
;  NOTE: Due to a bug, some files have the global attribute "Delta_time" 
;  of "delta_time".  We test for both here.
;
; !REVISION HISTORY:
;  23 Sep 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Copy global attributes from fIn to fOut
  ; and also make sure attributes begin w/ upper case
  ;=========================================================================

  ; Basic contact info
  fOut@Title                = str_sub_str( fIn@title, "GEOS-5.7.2", "GEOS-FP" )
  fOut@Contact              = fIn@contact
  fOut@References           = fIn@references

  ; Filename (disregard full directory path)
  fileName                  = notdir( fIn@filename )
  fileName                  = str_sub_str( fileName, "GEOS572", "GEOSFP" )
  fOut@Filename             = fileName

  ; Production & modification dates
  fOut@History              = fIn@history
  fOut@ProductionDateTime   = fIn@history
  fOut@ModificationDateTime = "Last modified on: " + systemfunc( "date" )

  ; Format
  fOut@Format               = fIn@format

  ; Spatial coverage
  if ( ( fIn@delta_lat .eq. 2 ) .or. ( fIn@delta_lat .eq. 4 ) ) then
    fOut@SpatialCoverage    = "global"
  else
    fOut@SpatialCoverage    = "regional"
  end if

  ; File info
  fOut@Conventions          = "COARDS"
  fOut@Version              = "GEOS-FP"
  fOut@Model                = "GEOS-5"
  fOut@Nlayers              = fIn@nlayers

  ; Start & end dates
  fOut@Start_Date           = fIn@start_date
  fOut@Start_Time           = fIn@start_time
  fOut@End_Date             = fIn@end_date
  fOut@End_Time             = fIn@end_time

  ; Time increment
  if      ( isatt( fIn, "delta_time" ) ) then
    fOut@Delta_Time         = fIn@delta_time
  else if ( isatt( fIn, "Delta_time" ) ) 
    fOut@Delta_Time         = fIn@Delta_time
  end if 
  end if

  ; Lon, lat increments
  fOut@Delta_Lon            = fIn@delta_lon
  fOut@Delta_Lat            = fIn@delta_lat

end
;EOC