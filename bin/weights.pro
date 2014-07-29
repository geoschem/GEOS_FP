pro weights, res

   ; Old grid
   OldType = CTM_Type( 'GEOS5', res=[ 5D/16D, 1D/4D ] )
   OldGrid = CTM_Grid( OldType )
   
   ; New grid
   NewType = CTM_Type( 'GEOS5', res=Res )
   NewGrid = CTM_Grid( NewType )

   ; Format string
   case ( Res[0] ) of
      0.625  : fmtStr = '(3x,12f8.4)'
      2    : fmtStr = '(3x,12f8.4)'
      2.5  : fmtStr = '(3x,12f8.4)'
      4    : fmtStr = '(3x,12f7.3)'  
      5    : fmtStr = '(3x,12f7.3)'  
      else : fmtStr = '(3x,12f8.4)'  
   endcase

   ; Filename
   case ( Res[0] ) of
      0.625   : WeightFile = 'weights_025x03125_to_05x0625.txt'
      2    : WeightFile = 'weights_025x03125_to_2x25.txt'
      2.5  : WeightFile = 'weights_025x03125_to_2x25.txt'
      4    : WeightFile = 'weights_025x03125_to_4x5.txt'
      5    : WeightFile = 'weights_025x03125_to_4x5.txt'
      else : WeightFile = ''
   endcase

   ; Make the weight file
   CTM_GetWeight, OldGrid, NewGrid, $
                  WeightFile=WeightFile, WeightFormat=fmtStr

end


;    IF ( IMX == 72 .and. JMX == 46 ) THEN
;       fmtStr = '(3x,12f7.3)'                           ! 4 x 5
;    ELSE IF ( IMX == 144 .and. JMX == 91 ) THEN
;       fmtStr = '(3x,12f8.4)'                           ! 2 x 2.5
;    ELSE IF ( IMX == 288 .and. JMX == 181 ) THEN
;       fmtStr = '(3x,12f9.5)'                           ! 1 x 1.25
;    ELSE IF ( IMX == 360 .and. JMX == 181 ) THEN
;       fmtStr = '(3x,12f8.4)'                           ! 1 x 1           
;    ENDIF
