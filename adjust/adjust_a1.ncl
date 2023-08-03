;------------------------------------------------------------------------------
;          Harvard University Atmospheric Chemistry Modeling Group            !
;------------------------------------------------------------------------------
;BOP
;
; !ROUTINE: adjust_a1
;
; !DESCRIPTION: Routine to adjust the existing GEOS-FP "A1" (surface
;  variables, 1-hour time averaged) data files for better COARDS compliance.
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
 undef( "adjust_a1" )
 procedure adjust_a1( inFile  : string, 
                      outFile : string  )
;
; !INPUT PARAMETERS:
;  inFile   : Name of GEOS-FP "A1" file w/ input data
;  outFile  : Name of GEOS-FP "A1" file w/ vars & atts adjusted for COARDS
;
; !LOCAL VARIABLES:
;
 local fIn,      time,     lat,      lon,      ALBEDO,   CLDTOT,  
       EFLUX,    EVAP,     FRSEAICE, FRSNO,    GRN,      GWETROOT, GWETTOP, 
       HFLUX,    LAI,      LWI,      LWGNT,    LWTUP,    PARDF,    PARDR,    
       PBLH,     PRECANV,  PRECCON,  PRECLSC,  PRECSNO,  PRECTOT,  QV2M,
       SEAICE00, SEAICE10, SEAICE20, SEAICE30, SEAICE40, SEAICE50, SEAICE60,
       SEAICE70, SEAICE80, SEAICE90, SLP,      SNODP,    SNOMAS,   SWGDN,
       TROPPT,   TS,       T2M,      U10M,     USTAR,    V10M,     Z0M,
       T10M,     Q850,     fOut
;
; !REMARKS:
;  (1) Adds a calendar attribute onto the "time" index variable.
;  (2) Ensures that the poles have latitudes -90 and +90 degrees.
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
  lat                    = fIn->lat
  lon                    = fIn->lon
  ALBEDO                 = fIn->ALBEDO       
  CLDTOT                 = fIn->CLDTOT
  EFLUX                  = fIn->EFLUX
  EVAP                   = fIn->EVAP
  FRSEAICE               = fIn->FRSEAICE
  FRSNO                  = fIn->FRSNO
  GRN                    = fIn->GRN
  GWETROOT               = fIn->GWETROOT
  GWETTOP                = fIn->GWETTOP
  HFLUX                  = fIn->HFLUX
  LAI                    = fIn->LAI
  LWI                    = fIn->LWI
  LWGNT                  = fIn->LWGNT
  LWTUP                  = fIn->LWTUP
  PARDF                  = fIn->PARDF
  PARDR                  = fIn->PARDR
  PBLH                   = fIn->PBLH
  PRECANV                = fIn->PRECANV
  PRECCON                = fIn->PRECCON
  PRECLSC                = fIn->PRECLSC
  PRECSNO                = fIn->PRECSNO
  PRECTOT                = fIn->PRECTOT
  QV2M                   = fIn->QV2M
  SEAICE00               = fIn->SEAICE00
  SEAICE10               = fIn->SEAICE10
  SEAICE20               = fIn->SEAICE20
  SEAICE30               = fIn->SEAICE30
  SEAICE40               = fIn->SEAICE40
  SEAICE50               = fIn->SEAICE50
  SEAICE60               = fIn->SEAICE60
  SEAICE70               = fIn->SEAICE70
  SEAICE80               = fIn->SEAICE80
  SEAICE90               = fIn->SEAICE90
  SLP                    = fIn->SLP
  SNODP                  = fIn->SNODP
  SNOMAS                 = fIn->SNOMAS
  SWGDN                  = fIn->SWGDN
  TROPPT                 = fIn->TROPPT
  TS                     = fIn->TS
  T2M                    = fIn->T2M
  U10M                   = fIn->U10M
  USTAR                  = fIn->USTAR
  V10M                   = fIn->V10M
  Z0M                    = fIn->Z0M
  T10M                   = fIn->T10M
  Q850                   = fIn->Q850 

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
    ALBEDO&lat           = lat
    CLDTOT&lat           = lat
    EFLUX&lat            = lat
    EVAP&lat             = lat
    FRSEAICE&lat         = lat
    FRSNO&lat            = lat
    GRN&lat              = lat
    GWETROOT&lat         = lat
    GWETTOP&lat          = lat
    HFLUX&lat            = lat
    LAI&lat              = lat
    LWI&lat              = lat
    LWGNT&lat            = lat
    LWTUP&lat            = lat
    PARDF&lat            = lat
    PARDR&lat            = lat
    PBLH&lat             = lat
    PRECANV&lat          = lat     
    PRECCON&lat          = lat
    PRECLSC&lat          = lat
    PRECSNO&lat          = lat
    PRECTOT&lat          = lat
    QV2M&lat             = lat
    SEAICE00&lat         = lat
    SEAICE10&lat         = lat
    SEAICE20&lat         = lat
    SEAICE30&lat         = lat
    SEAICE40&lat         = lat
    SEAICE50&lat         = lat
    SEAICE60&lat         = lat
    SEAICE70&lat         = lat
    SEAICE80&lat         = lat
    SEAICE90&lat         = lat
    SLP&lat              = lat
    SNODP&lat            = lat
    SNOMAS&lat           = lat
    SWGDN&lat            = lat
    TROPPT&lat           = lat
    TS&lat               = lat
    T2M&lat              = lat
    U10M&lat             = lat
    USTAR&lat            = lat
    V10M&lat             = lat
    Z0M&lat              = lat
    T10M&lat             = lat
    Q850&lat             = lat
  end if

  ;=========================================================================
  ; Modify variable attributes for COARDS compliance
  ;=========================================================================

  ; The time index array should have a calendar specified
  ; Redefine time for all relevant variables
  time@calendar          = "gregorian"
  time&time              = time
  ALBEDO&time            = time
  CLDTOT&time            = time
  EFLUX&time             = time
  EVAP&time              = time
  FRSEAICE&time          = time
  FRSNO&time             = time
  GRN&time               = time
  GWETROOT&time          = time
  GWETTOP&time           = time
  HFLUX&time             = time
  LAI&time               = time
  LWI&time               = time
  LWGNT&time             = time
  LWTUP&time             = time
  PARDF&time             = time
  PARDR&time             = time
  PBLH&time              = time
  PRECANV&time           = time
  PRECCON&time           = time
  PRECLSC&time           = time
  PRECSNO&time           = time
  PRECTOT&time           = time
  QV2M&time              = time
  SEAICE00&time          = time
  SEAICE10&time          = time
  SEAICE20&time          = time
  SEAICE30&time          = time
  SEAICE40&time          = time
  SEAICE50&time          = time
  SEAICE60&time          = time
  SEAICE70&time          = time
  SEAICE80&time          = time
  SEAICE90&time          = time
  SLP&time               = time
  SNODP&time             = time
  SNOMAS&time            = time
  SWGDN&time             = time
  TROPPT&time            = time
  TS&time                = time
  T2M&time               = time
  U10M&time              = time
  USTAR&time             = time
  V10M&time              = time
  Z0M&time               = time
  T10M&time              = time
  Q850&time              = time

  ; Quantities that are unitless need to have a unit string of "1"
  ALBEDO@units           = "1" ;
  CLDTOT@units           = "1" ;
  FRSEAICE@units         = "1" ;
  FRSNO@units            = "1" ;
  GRN@units              = "1" ;
  GWETROOT@units         = "1" ;
  GWETTOP@units          = "1" ;
  LWI@units              = "1" ;
  SEAICE00@units         = "1" ;
  SEAICE10@units         = "1" ;
  SEAICE20@units         = "1" ;
  SEAICE30@units         = "1" ;
  SEAICE40@units         = "1" ;
  SEAICE50@units         = "1" ;
  SEAICE60@units         = "1" ;
  SEAICE70@units         = "1" ;
  SEAICE80@units         = "1" ;
  SEAICE90@units         = "1" ;

  ; Add missing units
  LWGNT@units            = "W m-2"
  LWTUP@units            = "W m-2"

  ; Add common COARDS variable attributes
  add_coards_var_atts( ALBEDO   ) 
  add_coards_var_atts( CLDTOT   )
  add_coards_var_atts( EFLUX    )
  add_coards_var_atts( EVAP     )
  add_coards_var_atts( FRSEAICE )
  add_coards_var_atts( FRSNO    )
  add_coards_var_atts( GRN      )
  add_coards_var_atts( GWETROOT )
  add_coards_var_atts( GWETTOP  )
  add_coards_var_atts( HFLUX    )
  add_coards_var_atts( LAI      )
  add_coards_var_atts( LWI      )
  add_coards_var_atts( LWGNT    )
  add_coards_var_atts( LWTUP    )
  add_coards_var_atts( PARDF    )
  add_coards_var_atts( PARDR    )
  add_coards_var_atts( PBLH     )
  add_coards_var_atts( PRECANV  )
  add_coards_var_atts( PRECCON  )
  add_coards_var_atts( PRECLSC  )
  add_coards_var_atts( PRECSNO  )
  add_coards_var_atts( PRECTOT  )
  add_coards_var_atts( QV2M     )
  add_coards_var_atts( SEAICE00 )
  add_coards_var_atts( SEAICE10 )
  add_coards_var_atts( SEAICE20 )
  add_coards_var_atts( SEAICE30 )
  add_coards_var_atts( SEAICE40 )
  add_coards_var_atts( SEAICE50 )
  add_coards_var_atts( SEAICE60 )
  add_coards_var_atts( SEAICE70 )
  add_coards_var_atts( SEAICE80 )
  add_coards_var_atts( SEAICE90 )
  add_coards_var_atts( SLP      )
  add_coards_var_atts( SNODP    )
  add_coards_var_atts( SNOMAS   )
  add_coards_var_atts( SWGDN    )
  add_coards_var_atts( TROPPT   )
  add_coards_var_atts( TS       )
  add_coards_var_atts( T2M      )
  add_coards_var_atts( U10M     )
  add_coards_var_atts( USTAR    )
  add_coards_var_atts( V10M     )
  add_coards_var_atts( Z0M      )
  add_coards_var_atts( T10M     )
  add_coards_var_atts( Q850     )  

  ;=========================================================================
  ; Save to output file
  ;=========================================================================

  ; Print data (uncomment for debugging)
  ;printVarSummary( time     )
  ;printVarSummary( lon      )
  ;printVarSummary( lat      )
  ;printVarSummary( ALBEDO   ) 
  ;printVarSummary( CLDTOT   )
  ;printVarSummary( EFLUX    )
  ;printVarSummary( EVAP     )
  ;printVarSummary( FRSEAICE )
  ;printVarSummary( FRSNO    )
  ;printVarSummary( GRN      ) 
  ;printVarSummary( GWETROOT )
  ;printVarSummary( GWETTOP  )
  ;printVarSummary( HFLUX    )
  ;printVarSummary( LAI      )
  ;printVarSummary( LWI      )
  ;printVarSummary( LWGNT    )
  ;printVarSummary( LWTUP    )
  ;printVarSummary( PARDF    )
  ;printVarSummary( PARDR    )
  ;printVarSummary( PBLH     )
  ;printVarSummary( PRECANV  )
  ;printVarSummary( PRECCON  )
  ;printVarSummary( PRECLSC  )
  ;printVarSummary( PRECSNO  )
  ;printVarSummary( PRECTOT  )
  ;printVarSummary( QV2M     )
  ;printVarSummary( SEAICE00 )
  ;printVarSummary( SEAICE10 )
  ;printVarSummary( SEAICE20 )
  ;printVarSummary( SEAICE30 )
  ;printVarSummary( SEAICE40 )
  ;printVarSummary( SEAICE50 )
  ;printVarSummary( SEAICE60 )
  ;printVarSummary( SEAICE70 )
  ;printVarSummary( SEAICE80 )
  ;printVarSummary( SEAICE90 )
  ;printVarSummary( SLP      )
  ;printVarSummary( SNODP    )
  ;printVarSummary( SNOMAS   )
  ;printVarSummary( SWGDN    )
  ;printVarSummary( TROPPT   )
  ;printVarSummary( TS       )
  ;printVarSummary( T2M      )
  ;printVarSummary( U10M     )
  ;printVarSummary( USTAR    )
  ;printVarSummary( V10M     )
  ;printVarSummary( Z0M      )
  ;printVarSummary( T10M     )
  ;printVarSummary( Q850     )  
  ;exit

  ; Open output file (remove prior version)
  system( "rm -f " + outFile )
  fOut                   = addfile( outFile, "c" )
  fOut->time             = time
  fOut->lat              = lat
  fOut->lon              = lon
  fOut->ALBEDO           = ALBEDO
  fOut->CLDTOT           = CLDTOT
  fOut->EFLUX            = EFLUX
  fOut->EVAP             = EVAP
  fOut->FRSEAICE         = FRSEAICE
  fOut->FRSNO            = FRSNO
  fOut->GRN              = GRN
  fOut->GWETROOT         = GWETROOT
  fOut->GWETTOP          = GWETTOP
  fOut->HFLUX            = HFLUX
  fOut->LAI              = LAI
  fOut->LWI              = LWI
  fOut->LWGNT            = LWGNT
  fOut->LWTUP            = LWTUP
  fOut->PARDF            = PARDF
  fOut->PARDR            = PARDR
  fOut->PBLH             = PBLH
  fOut->PRECANV          = PRECANV
  fOut->PRECCON          = PRECCON
  fOut->PRECLSC          = PRECLSC
  fOut->PRECSNO          = PRECSNO
  fOut->PRECTOT          = PRECTOT
  fOut->QV2M             = QV2M
  fOut->SEAICE00         = SEAICE00 
  fOut->SEAICE10         = SEAICE10
  fOut->SEAICE20         = SEAICE20
  fOut->SEAICE30         = SEAICE30
  fOut->SEAICE40         = SEAICE40
  fOut->SEAICE50         = SEAICE50
  fOut->SEAICE60         = SEAICE60
  fOut->SEAICE70         = SEAICE70
  fOut->SEAICE80         = SEAICE80
  fOut->SEAICE90         = SEAICE90
  fOut->SLP              = SLP
  fOut->SNODP            = SNODP
  fOut->SNOMAS           = SNOMAS 
  fOut->SWGDN            = SWGDN
  fOut->TROPPT           = TROPPT
  fOut->TS               = TS
  fOut->T2M              = T2M 
  fOut->U10M             = U10M
  fOut->USTAR            = USTAR 
  fOut->V10M             = V10M
  fOut->Z0M              = Z0M
  fOut->T10M             = T10M
  fOut->Q850             = Q850

  ; Copy global attributes from the input file and
  ; adjust the attribute names to be COARDS compliant
  add_coards_global_atts( fIn, fOut )
  
  ; Redefine the title string
  fOut@Title             = "GEOS-FP 1-hour time-averaged parameters (A1), processed for GEOS-Chem input"

end
;EOC