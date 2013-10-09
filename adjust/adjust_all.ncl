;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !MODULE: adjust_all
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
 load "./adjust_a1.ncl"
 load "./adjust_a3cld.ncl"
 load "./adjust_a3dyn.ncl"
 load "./adjust_a3mstc.ncl"
 load "./adjust_a3mste.ncl"
 load "./adjust_cn.ncl"
 load "./adjust_i3.ncl"
;
; !INTERFACE:
;
 undef( "adjust_all" )
;
; !CALLING SEQUENCE:
;  ncl 'dir="/directory/path/to/files/"' adjust_all.ncl

; !REMARKS:
;  (1) Make sure the directory path ends in a "/".
;
; !REVISION HISTORY:
;  11 Sep 2013 - R. Yantosca - Initial version
;  09 Oct 2013 - R. Yantosca - Rewrote to handle missing files 
;  09 Oct 2013 - R. Yantosca - Need to redefine the "files" array with :=
;                              because its size can change
;EOP
;------------------------------------------------------------------------------
;BOC
begin

  ;=========================================================================
  ; Make sure the user specifes an input directory
  ;=========================================================================
  if ( .not. isvar( "dir" ) ) then
    print( "Input directory not specified!  Exiting ..." )
    exit
  end if

  ;=========================================================================
  ; Adjust the A1 files
  ;=========================================================================

  ; Find all GEOS572*A1 files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*A1*.nc" )
  
  ; Loop over each A1 file
  do n = 0, dimsizes( files )-1 

    ; Set input file name
    inFile  = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then 

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_a1( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the A3cld files
  ;=========================================================================

  ; Find all GEOS572*A3cld files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*A3cld*.nc" )

  ; Loop over each A3cld file
  do n = 0, dimsizes( files )-1 
    
    ; Set input file name
    inFile  = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_a3cld( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the A3dyn files
  ;=========================================================================

  ; Find all GEOS572*A3dyn files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*A3dyn*.nc" )

  ; Loop over each A3dyn file
  do n = 0, dimsizes( files )-1
    
    ; Set input file name
    inFile  = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) 

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_a3dyn( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the A3mstC files
  ;=========================================================================

  ; Find all GEOS572*A3mstC files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*A3mstC*.nc" )

  ; Loop over each A3mstC file
  do n = 0, dimsizes( files )-1
    
    ; Set input file name
    inFile  = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_a3mstc( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the A3mstE files
  ;=========================================================================

  ; Find all GEOS572*A3mstE files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*A3mstE*.nc" )

  ; Loop over each CN file
  do n = 0, dimsizes( files )-1 
    
    ; Set file names
    inFile = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then 

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_a3mste( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the CN files
  ;=========================================================================
  
  ; Find all GEOS572*CN files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*CN*.nc" )

  ; Loop over each CN file
  do n = 0, dimsizes( files )-1 
    
    ; Set input file name
    inFile  = files(n)
    
    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust files
      adjust_cn( inFile, outFile )
    end if
  end do

  ;=========================================================================
  ; Adjust the I3 file
  ;=========================================================================

  ; Find all GEOS572*I3 files in this directory
  files := systemfunc( "ls -1 " + dir + "GEOS572*I3*.nc" )

  ; Loop over each I3 file
  do n = 0, dimsizes( files )-1   
    
    ; Set file names
    inFile  = files(n)

    ; Make sure input file exists
    if ( .not. ismissing( inFile ) ) then 

      ; Set output file name
      outFile = str_sub_str( inFile, "GEOS572", "GEOSFP" )
      print( "Now processing: "+ notdir( outFile ) )

      ; Readjust the A3cld files
      adjust_i3( inFile, outFile )
    end if
  end do

end
;EOC
