;-----------------------------------------------------------------------
;+
; NAME:
;        LWI_MASK
;
; PURPOSE:
;        Creates the LWI mask template used by the MERRA data
;        processing code.  Also saves out the FRLANDICE and FRLAND
;        template files required by MerraA1Module.F90 for regridding 
;        the MERRA SNOMAS field.
;
; CATEGORY:
;        MERRA regrdiding utilities
;
; CALLING SEQUENCE:
;        LWI_MASK [, Keywords ]
;
; INPUTS:
;        None
;
; KEYWORD PARAMETERS:
;        LWIMASKFILE -> Name for the output LWI mask file.  
;             Default is lwi_mask.05x0666.bin.
;
;        FRLANDICEFILE -> Name for the output FRLANDICE file.
;             Default is frlandice.05x0666.bin
;
;        FRLANDFILE -> Name for the output FRLAND file.
;             Default is frland.05x0666.bin
;
;        _EXTRA=e -> Picks up extra keywords.
;
; OUTPUTS:
;        None
;
; SUBROUTINES:
;        External Subroutines Required:
;        ===========================================
;        CTM_TYPE  (function)   CTM_GRID (function) 
;        HDF_GETSD (function)   OPEN_FILE
;
; REQUIREMENTS:
;        Requires routines from the GAMAP package.
;
; NOTES:
;        None
;
; EXAMPLE:
;        LWI_MASK
;
;             ; Creates template files for MERRA regridding.
;
; MODIFICATION HISTORY:
;        bmy, 25 Aug 2010: INITIAL VERSION
;        bmy, 25 Aug 2010: - Now save out FRLANDICE itself
;        bmy, 25 Aug 2010: - Now save out FRLAND 
;
;-
; Copyright (C) 2005-2010, Bob Yantosca, Harvard University
; This software is provided as is without any warranty whatsoever. 
; It may be freely used, copied or distributed for non-commercial 
; purposes. This copyright notice must be kept with any copy of 
; this software. If this software shall be used commercially or 
; sold as part of a larger package, please contact the author.
; Bugs and comments should be directed to yantosca@seas.harvard.edu
; with subject "IDL routine lwi_mask.pro"
;-----------------------------------------------------------------------


pro lwi_mask, LwiMaskFile=LwiMaskFile,     $
              FrLandIceFile=FrLandIceFile, $
              FrLandFile=FrLandFile,       $
              _EXTRA=e

   ;-------------------------
   ; Initialization
   ;-------------------------

   ; External functions
   FORWARD_FUNCTION CTM_Type, CTM_Grid, HDF_GetSd

   ; File name (change if necessary)
   Dir     = '/as/scratch/bmy/MERRA/'
   File    =  Dir + 'MERRA300.prod.assim.const_2d_asm_Nx.00000000.hdf'

   ; MERRA 0.5 x 0.666 grid
   InType  = CTM_Type( 'MERRA', Res=[2D/3D, 1D/2D] )
   InGrid  = CTM_Grid( InType )

   ;-------------------------
   ; Open land frac file
   ;-------------------------

   ; Open the HDF file and get the file ID # (FID)
   fId = HDF_SD_START( File, /Read )
   IF ( fId lt 0 ) then MESSAGE, 'Error opening file!'

   ; Get data fields
   Time    = HDF_GetSd( fId, 'Time'      )    ; Hours sine 0GMT 1/1/1993
   Land    = HDF_GetSd( fId, 'FRLAND'    )    ; Land fraction
   LandIce = HDF_GetSd( fId, 'FRLANDICE' )    ; Land ice fraction
   Lake    = HDF_GetSd( fId, 'FRLAKE'    )    ; Lake fraction
   Ocean   = HDF_GetSd( fId, 'FROCEAN'   )    ; Ocean fraction

   ; Close file
   HDF_SD_END, fId

   ;-------------------------
   ; Construct data array
   ;-------------------------

   ; Output array
   Data = FltArr( InGrid.IMX, InGrid.JMX ) + 1e0

   ; Assign ocean boxes a value of 0
   Ind     = Where( Ocean gt 0e0 )
   if ( Ind[0] ge 0 ) then Data[Ind] = 0e0

   ; Assign land boxes a value of 1
   Ind     = Where( Land + LandIce + Lake gt 0.0 )
   if ( Ind[0] ge 0 ) then Data[Ind] = 1e0

   ;### Debug
   ;tvmap, data, ingrid.xmid, ingrid.ymid, /sample
   ;print, Min( data, max=m ), mWhere

   ;-------------------------
   ; Write LWI mask file
   ;-------------------------

   ; Output file name
   if ( N_Elements( LwiMaskFile ) eq 0 ) $
      then LwiMaskFile = 'lwi_mask.05x0666.bin'

   ; Open file for output
   Open_File, LwiMaskFile, Ilun,                     $
              /Get_Lun, /F77_Unformatted,            $
              /Write,   Swap_Endian=Little_Endian(), $
              _EXTRA=e

   ; Save data array
   WriteU, Ilun, Data
  
   ; Close file
   Close,    Ilun
   Free_Lun, Ilun

   ;-------------------------
   ; Write FRLANDICE file
   ;-------------------------
;-----------------------------------------------------------------------
; NOTE: This is no longer needed, we use the new algorithm
; Preserve the code here (bmy, 8/26/10)
;   ; Output file name
;   if ( N_Elements( FrLandIceFile ) eq 0 ) $
;      then FrLandIceFile = 'frlandice_mmH2O.05x0666.bin'
;
;   ; Create output array
;   Data = FltArr( InGrid.IMX, InGrid.JMX )
;
;   ; Set default snow mass [mm H2O = kg/m2] as 4000 where
;   ; FRLANDICE > 0.9.  This is what they did in GEOS-5.
;   Ind = Where( LandIce gt 0.9 )
;   if ( Ind[0] ge 0 ) then Data[Ind] = 4000e0
;-----------------------------------------------------------------------

   ; Output file name
   if ( N_Elements( FrLandIceFile ) eq 0 ) $
      then FrLandIceFile = 'frlandice.05x0666.bin'

   ; Open file for output
   Open_File, FrLandIceFile, Ilun,                   $
              /Get_Lun, /F77_Unformatted,            $
              /Write,   Swap_Endian=Little_Endian(), $
              _EXTRA=e

   ; Save data array
   WriteU, Ilun, LandIce
  
   ; Close file
   Close,    Ilun
   Free_Lun, Ilun

   ;-------------------------
   ; Write FRLAND file
   ;-------------------------

   ; Output file name
   if ( N_Elements( FrLandFile ) eq 0 ) $
      then FrLandFile = 'frland.05x0666.bin'

   ; Open file for output
   Open_File, FrLandFile, Ilun,                      $
              /Get_Lun, /F77_Unformatted,            $
              /Write,   Swap_Endian=Little_Endian(), $
              _EXTRA=e

   ; Save data array
   WriteU, Ilun, Land
  
   ; Close file
   Close,    Ilun
   Free_Lun, Ilun

end
