;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !MODULE: make_native_sfc_wind
;
; !DESCRIPTION: Creates the native surface wind file used for the GEOS-Chem
;  dust emissions module (by Dave Ridley et al).
;\\
;\\
; !USES:
;
 ; NCL routines
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

 ; Local routines
 load "../adjust/notdir.ncl"
;
; !INTERFACE:
;
 undef( "make_native_sfc_wind" )
;
; !CALLING SEQUENCE:
;  ncl 'date="YYYYMMDD"'       \
;      'inDir="/input/dir/"'   \
;      'outDir="/output/dir/"' make_native_sfc_wind.ncl

; !REMARKS:
;  (1) Make sure the directory paths end in a "/".
;  (2) Assumes all data files for a given day have the same dimensions,
;      (which is true for GEOS-FP raw data).
;
; !REVISION HISTORY:
;  09 Oct 2013 - R. Yantosca - Initial version
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Error check arguments
  ;=========================================================================

  ; Make sure date is passed
  if ( .not. isvar( "date" ) ) then
    print( "Date not specified!  Exiting ..." )
    exit
  end if

  ; Make sure inDir is passed
  if ( .not. isvar( "inDir" ) ) then
    print( "Input directory not specified!  Exiting ..." )
    exit
  end if

  ; Make sure outDir is passed
  if ( .not. isvar( "outDir" ) ) then
    print( "Output directory not specified!  Exiting ..." )
    exit
  end if
 
  ;=========================================================================
  ; Compute wind speed from U10M and V10M
  ;
  ; NOTE: addfiles doesn't work on the GEOS-FP met field files (netCDF-4), 
  ; probably due  to the chunking settings.  Read files individually.
  ;=========================================================================

  ; Find all GEOS572*A1 files in this directory
  files = systemfunc( "ls -1 " + inDir + "GEOS.fp.*slv_Nx." + date + "*.nc4" )

  ; Loop over files
  do n = 0, dimsizes( files )-1 

    ; Read data from each file
    fIn  = addfile( files(n), "r" )

    ; Is it the first iteration?
    if ( n .eq. 0 ) then

      ;=====================================================================
      ; Create index arrays
      ;=====================================================================

      ; Create time index array (minutes since start of today)
      time                 = (/   30,   90,  150,  210,  270,  330, \
                                 390,  450,  510,  570,  630,  690, \
                                 750,  810,  870,  930,  990, 1050, \
                                1110, 1170, 1230, 1290, 1350, 1410  /)
      time!0               = "time"
      time&time            = time
      y                    = str_get_cols( date, 0, 3 )
      m                    = str_get_cols( date, 4, 5 )
      d                    = str_get_cols( date, 6, 7 )
      time@time_increment  = "010000"
      time@begin_time      = "000000"
      time@begin_date      = date
      time@delta_t         = "0000-00-00 01:00:00"
      time@units           = "minutes since "+y+"-"+m+"-"+d+" 00:00:00.0"
      time@long_name       = "time"

      ; Create lat index array
      lat                  = doubletofloat( fIn->lat )
      lat!0                = "lat"
      lat&lat              = lat
      lat@units            = "degrees_north"
      lat@long_name        = "latitude"
 
      ; Create lon index array
      lon                  = doubletofloat( fIn->lon )
      lon!0                = "lon"
      lon&lon              = lon
      lon@units            = "degrees_east"
      lon@long_name        = "longitude"

      ;=====================================================================
      ; Create data array for surface wind speed
      ;=====================================================================

      ; Get dimensions
      nTime                = dimsizes( time )
      nLat                 = dimsizes( lat  )
      nLon                 = dimsizes( lon  )

      ; Define the speed array
      speed                = new( (/ nTime, nLat, nLon /), float )
      speed!0              = "time"
      speed!1              = "lat"
      speed!2              = "lon"
      speed&time           = time
      speed&lat            = lat
      speed&lon            = lon
      speed@add_offset     = 0.0
      speed@scale_factor   = 1.0
      speed@missing_value  = 1e15
      speed@_FillValue     = 1e15
      speed@gamap_category = "GMAO-2D"
      speed@units          = "m s-1"
      speed@long_name      = "Wind speed at 10m above the surface"

      ; Print info about speed (comment out if needed)
      ;printVarSummary( speed )
    end if

    ;=======================================================================
    ; On each iteration, compute the surface wind speed
    ;=======================================================================
    u                     = fIn->U10M
    v                     = fIn->V10M
    speed(n,:,:)          = sqrt( u^2 + v^2 ) 
    
  end do

  ;=========================================================================
  ; Write data to the output file
  ;=========================================================================

  ; Open file
  outFile                   = outDir + "GEOSFP." + date + ".SPEED.025x03125.nc"
  print( "Writing " + outFile + " to disk" )
  system( "rm -f " + outFile )
  fOut                      = addfile( outFile, "c" )

  ; Write global attributes
  fOut@Title                = "GEOS-FP 1-hour native resolution wind speed"
  fOut@Contact              = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)" ;
  fOut@References           = "www.geos-chem.org; wiki.geos-chem.org"
  fOut@Filename             = notdir( outFile )
  modDate                   = systemfunc( "date" )
  fOut@History              = "File generated on: " + modDate
  fOut@ProductionDateTime   = "File generated on: " + modDate
  fOut@ModificationDateTime = "Last modified on: "  + modDate
  fOut@Format               = "netCDF-3"
  fOut@SpatialCoverage      = "global"
  fOut@Conventions          = "COARDS"
  fOut@Version              = "GEOS-FP"
  fOut@Model                = "GEOS-5"
  fOut@Nlayers              = 72
  fOut@Start_Date           = date
  fOut@Start_Time           = "00:00:00.00000"
  fOut@End_Date             = date
  fOut@End_Time             = "23:59:59.99999"
  fOut@Delta_Time           = "010000"
  fOut@Delta_Lon            = "0.3125"
  fOut@Delta_Lat            = "0.25"

  ; Write index & data arrays
  fOut->time                = time
  fOut->lat                 = lat
  fOut->lon                 = lon
  fOut->SPEED               = speed

end
;EOC