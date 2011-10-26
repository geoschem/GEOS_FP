! $Id: CharpakModule.f,v 1.3 2009/07/29 21:33:41 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CharpakModule
!
! !DESCRIPTION: Module CharpakModule (bmy, 12/3/98, 11/8/06) contains the
!  CHARPAK routines by R. D. Stewart as well as some wrapper subroutines 
!  by Martin Schultz.
!\\
!\\
! !INTERFACE: 
! 
      MODULE CharpakModule
!
! !PUBLIC MEMBER FUNCTIONS:
!     
      ! Subroutines
      PUBLIC :: addtxt
      PUBLIC :: AlgnTxt
      PUBLIC :: BldTxt
      PUBLIC :: CntMat
      PUBLIC :: CntRep
      PUBLIC :: CopyTxt
      PUBLIC :: CSWAP
      PUBLIC :: DExtNum
      PUBLIC :: DFMAT
      PUBLIC :: DNum2Txt
      PUBLIC :: DRndDig
      PUBLIC :: DTabNum
      PUBLIC :: EWWDD
      PUBLIC :: FillStr
      PUBLIC :: FWWDD
      PUBLIC :: GetNum
      PUBLIC :: IExtNum
      PUBLIC :: IFMAT
      PUBLIC :: InsTxt
      PUBLIC :: INum2Txt
      PUBLIC :: KeyExt
      PUBLIC :: LetSub
      PUBLIC :: LINE
      PUBLIC :: MakFmt
      PUBLIC :: PGSTR
      PUBLIC :: PigLat
      PUBLIC :: PLINE
      PUBLIC :: REPLICATE
      PUBLIC :: RmvRep
      PUBLIC :: RmvTxt
      PUBLIC :: SCANSTR
      PUBLIC :: SExtNum
      PUBLIC :: ShtTxt
      PUBLIC :: SNum2Txt
      PUBLIC :: SQUEEZE
      PUBLIC :: SRndDig
      PUBLIC :: STRCOMPRESS
      PUBLIC :: STRIP
      PUBLIC :: STRLOWCASE
      PUBLIC :: STRPUT
      PUBLIC :: STRTRIM
      PUBLIC :: STRUPCASE
      PUBLIC :: TRANLC
      PUBLIC :: TRANUC
      PUBLIC :: TrimStr
      PUBLIC :: Txt2Dnum
      PUBLIC :: Txt2Inum
      PUBLIC :: Txt2Snum
      PUBLIC :: TxtExt
      PUBLIC :: TxtLoc
      PUBLIC :: TxtMac
      PUBLIC :: TxtSub
      PUBLIC :: StrRepl
      PUBLIC :: YMD_EXTRACT
      PUBLIC :: expandDate
      PUBLIC :: MakeCharArrayFromCharList
      
      ! Functions
      PUBLIC :: IMagD
      PUBLIC :: IMagI
      PUBLIC :: IMagS
      PUBLIC :: IsAlpCh
      PUBLIC :: IsChar
      PUBLIC :: IsNum
      PUBLIC :: IsPrnCh
      PUBLIC :: LENTRIM
      PUBLIC :: STRLEN
      PUBLIC :: STRPOS
!
! !REMARKS:
!  CHARPAK by R.D. Stewart, 1992-1998
!  F90 module by Bob Yantosca (12/3/98)
!
! !REVISION HISTORY:
!  (1) Bob Yantosca added "IMPLICIT NONE" statements.
!  (2) EXTERNAL attributes are commented out since all of the subroutines
!      and functions listed in CHARPAK_MOD are have module scope, and
!      thus can be "seen" by all of the other subroutines. 
!  (3) Functions have been converted to the F90 style, e.g
!         FUNCTION myfunction( myarg )  result ( myresult ) 
!  (4) Function LENTRIM has been replaced by a wrapper to the
!      the F90 intrinsic function, LEN_TRIM.
!  (5) Added STRREPL, YMD_EXTRACT, EXPAND_DATE from GEOS-Chem (bmy, 11/8/06)
!  (6) Added ProTeX headers (for some routines, for now) (bmy, 7/28/09)
!  26 Oct 2011 - R. Yantosca - Added MakeCharArrayFromCharList
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS

      SUBROUTINE AddTxt(Tadd,Text)
!
!     Add (append) text contained in Tadd to the end of Text
!     NOTE: Any blanks trailing Text are ignored.  However,
!           blanks preceding non-blank characters in Tadd
!           are attached. (e.g., Tadd= '  xx', Text = 'yyz '
!           will return text = 'yyz  xx')
!
!     Author: R.D. Stewart
!     Date:   Sept. 27, 1994
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                   CHARPAK.FOR
!
      CHARACTER*(*) Tadd,Text
      INTEGER ilen,i1,i2,iadd
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

      ilen = LEN(text)
      i1 = LENTRIM(text)
      i2 = LENTRIM(Tadd)
      IF (i2.GT.0) THEN
        iadd = MIN0(ilen,i1+i2)
        text(i1+1:iadd) = Tadd(1:i2)
      ENDIF

      RETURN
      END SUBROUTINE addtxt
! 
!============================================================================
!
      SUBROUTINE AlgnTxt(col,ijust,ins,Tins,Tout)
!
!     PURPOSE: This routine transfer characters from one string
!              to another.  The text is aligned at a specified
!              column and can be left, right, or centered
!              justified.  In addition, the text can be inserted
!              or just over-write the existing text.
!
!     INPUT:   col: reference column
!              IF  ijust < 0 left justify text at col
!                  ijust = 0 center justify text at col
!                  ijust > 0 right justify text at col
!              IF  ins = true, insert text
!                  ins = false, overwrite any existing text
!              Tins: characters to insert into Tout
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                   CHARPAK.FOR
!        InsTxt                    CHARPAKA.FOR
!        CopyTxt                   CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) Tins,Tout
      INTEGER col,ijust
      LOGICAL ins,add
      INTEGER ic,LI,LIT,LT,T1,T2
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

!     Modified 5-DEC-1995 so that when text is right justified at col,
!     trailing blanks are ignored.
      IF (ijust.GT.0) THEN
        LI=LENTRIM(Tins)
      ELSE
        LI = LEN(Tins)
      ENDIF

      IF (LI.GT.0) THEN
        LT = LEN(Tout)
        add = .true.
        T1 = 1
        T2 = LI
        IF (ijust.LT.0) THEN
!         Left justify text at column, ic.
          ic = col
          IF (ic.LE.0) THEN
!           Not all of the text will be inserted/copied
            ic = 1
            T1 = 2-ic
            T2 = LI
            IF (T1.GT.LI) add = .false.
          ENDIF
        ELSEIF (ijust.GT.0) THEN
!         Right justify text at column, ic.
          ic = col - LI
          IF (ic.LE.0) THEN
!           Not all of the text will be inserted/copied
            T1 = LI+1-col
            T2 = LI
            ic = 1
            IF (T1.GT.LI) add = .false.
          ENDIF
        ELSEIF (ijust.EQ.0) THEN
!         Center justify text at column, ic.
          LIT = LENTRIM(Tins)
          ic = col - LIT/2
          IF (ic.LE.0) THEN
!           Not all of the text will be inserted/copied
            T1 = 2-ic
            T2 = LI
            ic = 1
            IF (T1.GT.LI) add = .false.
          ENDIF
        ENDIF

        IF (add) THEN
!         Transfer Text
          IF (ins) THEN
            CALL InsTxt(ic,Tins(T1:T2),Tout)
          ELSE
            CALL CopyTxt(ic,Tins(T1:T2),Tout)
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE AlgnTxt
! 
!============================================================================
!
      SUBROUTINE BldTxt(Add,Bld,nspc)
!
!     Add characters in ADD to Bld with nspc spaces between them.
!     AUTHOR: RD STEWART
!     DATE:   DEC 6, 1995
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) Add,Bld
      INTEGER ilen,nspc
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

      ILEN=LENTRIM(BLD)+MAX(0,nspc+1)
      CALL AddTxt(Add,Bld(ilen:))

      RETURN
      END SUBROUTINE BldTxt
! 
!============================================================================
!
      SUBROUTINE CntMat(str1,str2,imat)
!
!     Count the number of characters in str1 that match
!     a character in str2.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!          LENTRIM                CharPak
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!     COMMENTS: Revised slightly (2-5-1996) so that trailing
!               blanks in str1 are ignored.  Revised again
!               on 3-6-1996.
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) str1,str2
      INTEGER imat
      INTEGER L1,L2,i,j !,LENTRIM
      LOGICAL again
!      EXTERNAL LENTRIM

      L1 = MAX(1,LENTRIM(str1))
      L2 = LEN(str2)
      imat = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            imat = imat+1
            again = .false.
          ELSEIF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CntMat
! 
!============================================================================
!
      SUBROUTINE CntRep(str1,str2,imat)
!
!     Count the number of times characters in str1 occur in
!     str2.  NOTE: CntRep is not quite the same as CntMat,
!     but can be used interchangably in some cases.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JAN. 17, 1995
!     AUTHOR: R.D. STEWART
!
!     COMMENTS: Revised slightly (2-5-1996) so that trailing
!               blanks in str1 are ignored. Revised again
!               on 3-6-1996.
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) str1,str2
      INTEGER imat
      INTEGER L1,L2,i,j !,LENTRIM
      LOGICAL again
!      EXTERNAL LENTRIM

      L1 = MAX(1,LENTRIM(str1))
      L2 = LEN(str2)
      imat = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            imat = imat+1
          ENDIF
          IF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CntRep
! 
!============================================================================
!
      SUBROUTINE CopyTxt(col,str1,str2)
!
!     PURPOSE: Write all of the characters in str1 into variable
!              str2 beginning at column, col.  If the length of str1
!              + col is longer than the number of characters str2
!              can store, some characters will not be transfered to
!              str2.  Any characters already existing in str2 will
!              will be overwritten.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   DEC. 24, 1993
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) str2,str1
      INTEGER col,ilt1,i1,i,j,ic

      i1 = LEN(str2)
      IF (i1.GT.0) THEN
        ilt1 = LEN(str1)
        IF (ilt1.GT.0) THEN
          ic = MAX0(col,1)
          i = 1
          j = ic
          DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
            str2(j:j) = str1(i:i)
            i = i + 1
            j = ic + (i-1)
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CopyTxt
! 
!============================================================================
!
      SUBROUTINE CSWAP(i,j,str)
!
!     Swap substring characters i and j in string, str.
!
!      AUTHOR: Robert D. Stewart
!        DATE: FEB 22, 1996
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER i,j,ilen
      CHARACTER*(*) str
      CHARACTER*1 swap

      ilen=LEN(str)
      IF (ilen.GT.1) THEN
        IF ((i.GE.1).AND.(i.LE.ilen)) THEN
          IF ((j.GE.1).AND.(j.LE.ilen)) THEN
            swap = str(j:j)
            str(j:j) = str(i:i)
            str(i:i) = swap
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CSWAP
! 
!============================================================================
!
      SUBROUTINE DExtNum(col,txt,lst,dnum,iflg)
!
!     This routine searches through txt beginning at column,
!     col, and extracts a section of text delimited by
!     the string of characters in lst.  After a string of
!     characters is selected, an attempt is made to convert
!     them to a double precision number, dnum.  A status flag
!     is returned with the following meaning(s)
!
!        iflg <= 0 no errors
!        iflg > 0 an error occured
!
!     COMMENTS: The work space only needs to be as large as
!               the largest string of characters representing
!               a number.  50 characters should be sufficient
!               for almost any practical situation.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JAN 23, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                    CHARPAK.FOR
!        STRIP                      CHARPAK.FOR
!        TxtExt                     CHARPAKB.FOR
!        Txt2DNum                   CHARPAK.FOR
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col,iflg
      CHARACTER*(*) txt,lst
      DOUBLE PRECISION dnum
      CHARACTER*50 work,fmat*15
      INTEGER ilen !,LENTRIM
!      EXTERNAL LENTRIM

      IF (lst.EQ.'\D') THEN
!       Use default delimiters for numeric data.
        CALL TxtExt(' ,=()',txt,col,work,iflg)
      ELSE
!       User specified list.
        CALL TxtExt(lst,txt,col,work,iflg)
      ENDIF

      IF (iflg.LE.0) THEN
!       Now try to convert the text to a double precision
!       number.
        CALL STRIP(work)
        ilen=LENTRIM(work)
        CALL MakFmt('F',ilen,0,fmat,iflg)
        IF (iflg.LE.0) THEN
          CALL Txt2Dnum(fmat,work,dnum,iflg)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE DExtNum
! 
!============================================================================
!
      SUBROUTINE DFMAT(col,dnum,fspec,text)
!
!     PURPOSE: Format a double precision number (dnum) and
!              store it as a character string beginning at
!              column COL in TEXT.  FSPEC is an adaptive
!              (or "smart") format specifier.
!
!      INPUTS:  col - a valid integer in the range 1 to LEN(text)
!             dnum - a double precision number
!
!       FSPEC: The specification of fspec is adaptive in the
!              since that the width of the field is automatically
!              adjusted to the necessary length so that the
!              number can be displayed in the appropriate format
!              (i.e., no more ******** in the output).  Also,
!              integer (I), floating point (F), and exponential (E)
!              format are intelligently applied to all data types
!              (e.g., you can apply an "I" format to a real number-
!              this just means you will not get any fractional
!              part).  If a negative sign is included in fspec,
!              the number is rounded.  Perhaps the examples below
!              will make it clearer
!                fspec = 'F3' = Fwww.xxx
!                        'F5' = Fwww.xxxxx
!                        '-F3' = round to the third significant
!                                digit and display as appropriate.
!                        'I' = display integer part of number
!                        'E3' = 0.xxxE+00
!                        'E4' = 0.xxxxE+00
!                        '-E1' = x.xE+00
!
!              Using the short-hand specifiers makes it much
!              easier to handle the output of numeric data, but
!              you will probably just have to "play" with them
!              a bit to get a feel for what the output looks like.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 23, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        FillStr                     CHARPAK.FOR
!        TRANU!                      CHARPAK.FOR
!        CopyTxt                     CHARPAK.FOR
!        LENTRIM                     CHARPAK.FOR
!        FWWDD                       CHARPAKA.FOR
!        EWWDD                       CHARPAKA.FOR
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col
      DOUBLE PRECISION dnum
      CHARACTER*(*) fspec,text
      INTEGER iw
      CHARACTER*35 number
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

!     Ensure number is blank
      CALL FillStr(' ',number)

!     Ensure that fspec is all uppercase letters with no extra spaces
!     between characters.  NOTE: fspec may be modified upon return.
      CALL TRANUC(fspec)

!     Call routine to translate number into an appropriate
!     character string.
      IF (INDEX(fspec,'I').GT.0) THEN
        CALL FWWDD(dnum,fspec,number)
      ELSEIF (INDEX(fspec,'E').GT.0) THEN
        CALL EWWDD(dnum,fspec,number)
      ELSE
        CALL FWWDD(dnum,fspec,number)
      ENDIF

      iw = LENTRIM(number)
      CALL CopyTxt(col,number(1:iw),text)

      RETURN
      END SUBROUTINE DFMAT
! 
!============================================================================
!
      SUBROUTINE DNum2Txt(dnum,icol,fmat,txt,iflg)
!
!     <DNum2Txt> converts a double precision number into
!     a string of equivalent characters and stores them
!     in the character variable Txt beginning at column,
!     icol.  The number is formated according to the
!     the standard FORTRAN specifier passed in fmat.
!     iflg is returned as a status flag indicating the success
!     or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      DOUBLE PRECISION dnum
      INTEGER icol,iflg
      CHARACTER*(*) txt,fmat

      WRITE(txt(icol:),fmt=fmat,iostat=iflg) dnum

      RETURN
      END SUBROUTINE DNum2Txt
! 
!============================================================================
!
      SUBROUTINE DRndDig(din,nd,dout)
!
!     <DRndDig> rounds a double precision number din to the
!     number of significant digits indicated by integer paramter
!     nd.  If nd < 0, rounding will be performed to the right
!     of the decimal place.  For nd > 0, rounding is performed to
!     the left of the decimal place.  The rounded number is
!     returned in dout.
!
!      AUTHOR: Robert D. Stewart
!        DATE: FEB 13, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      DOUBLE PRECISION din,dout
      INTEGER nd
      DOUBLE PRECISION xmult

!     Create multiplier for rounding
      xmult=10.0D+00**(nd)

!     Round number.
      dout = DNINT(din*xmult)/xmult

      RETURN
      END SUBROUTINE DRndDig
! 
!============================================================================
!
      SUBROUTINE DTabNum(col,dnum,fspec,txt)
!
!     PURPOSE: This routine is intended to help facilitate
!              the creation of tables.  Frequently, tables
!              look better when numbers are aligned on a
!              decimal point.  DTabTxt automates
!              the process of converting a double precision
!              number (dnum) into a string of characters
!              and then centering the decimal place about
!              column, col, in Txt.  The number is formated
!              according to the "smart" specifier fspec
!              (see DFMAT).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col,idec,ilen,icol
      DOUBLE PRECISION dnum
      CHARACTER*(*) fspec,txt,num*35
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

      num = ' '
      CALL DFMAT(1,dnum,fspec,num)
      idec=INDEX(num,'.')
      ilen=LENTRIM(num)
      IF (idec.GT.0) THEN
        icol=MAX(1,col-(idec-1))
      ELSE
!       No decimal place (integer), so right justify number at col
        icol=MAX(1,col-(ilen-1))
      ENDIF
      CALL AlgnTxt(icol,-1,.FALSE.,num(:ilen),txt)


      RETURN
      END SUBROUTINE DTabNum
! 
!============================================================================
!
      SUBROUTINE EWWDD(xnum,fspec,number)
!
!     EWWDD takes are D.P. number, xnum, and dynamically
!     creates a format specifier of type 'E', then,
!     the number is transfered to the charaacter variable
!     number.  As an added "benefit" this routine puts a 0
!     in front of the decimal place (when appropriate), e.g.,
!     instead of returning .03E+00, 0.03E+00 will be returned.
!
!     This routine only handles the fspec = 'E' format.  The
!     routine FWWDD handles the specs 'I' and 'F'.
!     Typically, the SFMAT, DFMAT, and IFMAT driver
!     routines should be used to format numbers.
!
!     IMPORTANT NOTE: I highly recommend that the
!     SFMAT, DFMAT, and IFMAT driver routines be used
!     instead of directly calling this routine.  Also,
!     I may be changing the format/algorithm in the
!     EWWDD and FWWDD routines to take advantage of the newer
!     routines in CHARAPAK.  Using the higher level
!     xFMAT routines should less the need to change
!     code in the future if these routines are modified.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 22, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        TRANU!                      CHARPAK.FOR
!        CopyTxt                     CHARPAK.FOR
!        InsTxt                      CHARPAK.FOR
!        IsNum                       CHARPAK.FOR
!        LENTRIM                     CHARPAK.FOR
!        IMagD                       CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER flen,ic,i,iexp,id,iochk
      LOGICAL lsd
      DOUBLE PRECISION xnum,xtmp
      CHARACTER*(*) number,fspec
      CHARACTER*12 fmat
      CHARACTER*2 tnum
      DOUBLE PRECISION ZERO,UNITY,HALF
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,UNITY=1.0D+00)
      INTEGER ww,itmp ! LENTRIM,IMagD
!      LOGICAL ISNUM
!      EXTERNAL ISNUM,LENTRIM,IMagD

!     Ensure character variables are blank
      tnum = '  '

!     Ensure that fspec is all uppercase letters with no extra spaces
!     between characters.
      call tranuc(fspec)

!     length of fspec minus trailing blanks
      flen = LENTRIM(fspec)
      ww = LEN(number)

!     Decode short-hand specifier and build a valid FORTRAN
!     format specifier.

!     Locate the three numeric digit in fspec.  These will be used to
!     decide how many decimal places will be generated.
      ic = 0
      DO i=1,flen
        IF (ISNUM(fspec(i:i))) THEN
          IF (ic.LT.3) THEN
            ic = ic + 1
            tnum(ic:ic) = fspec(i:i)
          ENDIF
        ENDIF
      ENDDO

!     Here iexp is the nearest power of ten to xnum.
      iexp = IMagD(xnum)

      IF (ic.GT.0) THEN
        read(tnum,fmt='(i2)',iostat=iochk) id
        id = MIN0(id,ww-8)
      ELSE
        id = 0
      ENDIF

!     Store xnum in xtmp
      xtmp = xnum

      lsd = (INDEX(fspec,'-').GT.0)
      IF (lsd) THEN
        id = MAX0(id-1,0)
        fmat = '(1PExxx.xxx)'
        WRITE(fmat(5:11),fmt='(I3,A1,I3)',iostat=iochk) ww,'.',id
      ELSE
!       Need at least one decimal place
        id = MAX0(1,id)
        fmat = '(Exxx.xxx)'
        WRITE(fmat(3:9),fmt='(I3,A1,I3)',iostat=iochk) ww,'.',id
      ENDIF
      CALL STRIP(fmat)
      WRITE(number,fmt=fmat,iostat=iochk) ABS(xtmp)
      CALL STRIP(number)
      IF (number(1:1).EQ.'.') THEN
        IF (xtmp.LT.zero) THEN
!          Add -zero before decimal place
           CALL InsTxt(1,'-0',number)
        ELSE
!         Add -zero before decimal place
          CALL InsTxt(1,'0',number)
        ENDIF
      ELSEIF (xtmp.LT.zero) THEN
!       Add - before the number
        CALL InsTxt(1,'-',number)
      ENDIF

      itmp = LENTRIM(number)

!     Remove any excess blanks
      CALL STRIP(number)
      CALL TRANUC(number)

      RETURN
      END SUBROUTINE EWWDD
! 
!============================================================================
!
      SUBROUTINE FillStr(ch,text)
!
!     PURPOSE: Fill text with the (single) character in ch.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      CHARACTER*1 ch
      INTEGER ilen,i

      ilen = LEN(text)
      DO i=1,ilen
        text(i:i) = ch
      ENDDO

      RETURN
      END SUBROUTINE FillStr
! 
!============================================================================
!
      SUBROUTINE FWWDD(xnum,fspec,number)
!
!     FWWDD takes a D.P. number, xnum, and dynamically
!     creates a format specifier.  Then, the number is
!     transfered to the character variable number.  As
!     an added "benefit" this routine puts a 0 in front
!     of the decimal place (when appropriate), e.g.,
!     instead of returning .03, 0.03 will be returned.
!
!     This routine can handle fspec = 'F' or 'I', the
!     routine EWWDD handles the exponential specifiers
!     (fspec='E').
!
!     IMPORTANT NOTE: I highly recommend that the
!     SFMAT, DFMAT, and IFMAT driver routines be used
!     instead of directly calling this routine.  Also,
!     I may be changing the format/algorithm in the
!     EWWDD and FWWDD routines to take advantage of some of
!     the newer routines in CHARPAK.  Using the higher level
!     xFMAT routines should less the need to change
!     code in the future if these routines are modified.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 22, 1993
!        REVISIONS:
!           Nov. 30, 1995 (rds)
!             Incorrectly handled '-Fx' specifiers when ABS(xnum) < 1.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        TRANUC                      CHARPAK.FOR
!        CopyTxt                     CHARPAK.FOR
!        InsTxt                      CHARPAK.FOR
!        IsNum                       CHARPAK.FOR
!        LENTRIM                     CHARPAK.FOR
!        IMagD                       CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER flen,ic,i,iasc,iexp,id,iochk
      LOGICAL lsd
      DOUBLE PRECISION xnum,xtmp
      CHARACTER*(*) number,fspec
      CHARACTER*12 fmat
      CHARACTER*2 tnum
      DOUBLE PRECISION ZERO,UNITY,HALF
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,UNITY=1.0D+00)
      INTEGER ww,itmp !LENTRIM,IMagD
!      LOGICAL IsNum
!      EXTERNAL IsNum,LENTRIM,IMagD

!     Ensure character variables are blank
      tnum = '  '

!     Ensure that fspec is all uppercase letters with no extra spaces
!     between characters.
      call tranuc(fspec)

!     length of fspec minus trailing blanks
      flen = LENTRIM(fspec)
      ww = LEN(number)

!     Decode short-hand specifier and build a valid FORTRAN
!     format specifier.

!     Locate the three numeric digit in fspec.  These will be used to
!     decide how many decimal places will be generated.
      ic = 0
      DO i=1,flen
        IF (ISNUM(fspec(i:i))) THEN
          IF (ic.LT.3) THEN
            ic = ic + 1
            tnum(ic:ic) = fspec(i:i)
          ENDIF
        ENDIF
      ENDDO

!     Here iexp is the nearest power of ten to xnum.
      iexp = IMagD(xnum)

      IF (ic.GT.0) THEN
        read(tnum,fmt='(i2)',iostat=iochk) iasc
        iasc = MIN0(iasc,ww-3-iexp)
      ELSE
        iasc = 0
      ENDIF

!     Store xnum in xtmp
      xtmp = xnum

      lsd = (INDEX(fspec,'-').GT.0)
      IF (lsd) THEN
!       Always keep at least one significant digit.
        iasc = MAX0(iasc,1)
      ENDIF

      IF (iexp.LT.0) THEN
        IF (lsd) THEN
          id = MAX0(ABS(iexp) + iasc - 1,0)
        ELSE
          id = iasc
        ENDIF
      ELSE
        IF (lsd) THEN
!         Need to round xtmp.  This section modified by RDS 30-Nov-1995
!         to fix a bug in the -Fx algorith for ABS(xnum) < 1.
          IF (ABS(xnum).LT.unity) THEN
            itmp = iasc+iexp-1
            CALL DRndDig(xnum,itmp,xtmp)
            id = MAX0(itmp,0)
          ELSE
            itmp = MAX0(iasc-iexp-1,-iexp)
            CALL DRndDig(xnum,itmp,xtmp)
            id = MAX0(itmp,0)
          ENDIF
        ELSE
          id = iasc
        ENDIF
      ENDIF

      CALL MakFmt('F',ww,id,fmat,iochk)
      IF (iochk.NE.0) THEN
        PAUSE 'BAD TROUBLE IN <FWWDD> ROUTINE'
        fmat='(F5.0)'
      ENDIF

      WRITE(number,fmt=fmat,iostat=iochk) ABS(xtmp)
      CALL STRIP(number)
      IF (number(1:1).EQ.'.') THEN
        IF (xtmp.LT.zero) THEN
!          Add -zero before decimal place
           CALL InsTxt(1,'-0',number)
        ELSE
!         Add -zero before decimal place
          CALL InsTxt(1,'0',number)
        ENDIF
      ELSEIF (xtmp.LT.zero) THEN
!       Add - before the number
        CALL InsTxt(1,'-',number)
      ENDIF

      itmp = LENTRIM(number)

      IF (INDEX(fspec,'I').GT.0) THEN
!       Remove all characters to the left (and including)
!       the decimal point
        IF (ABS(xtmp).LT.unity) THEN
          number = '0'
        ELSE
          id = INDEX(number,'.')
          DO i=id,ww
            number(i:i) = ' '
          ENDDO
        ENDIF
      ELSEIF (iasc.EQ.0) THEN
!       Remove all characters to the left (and including)
!       the decimal point
        IF (ABS(xtmp).LT.unity) THEN
          number = '0'
        ELSE
          id = INDEX(number,'.')
          DO i=id,ww
            number(i:i) = ' '
          ENDDO
        ENDIF
      ELSEIF (number(itmp:itmp).EQ.'.') THEN
        number(itmp:itmp) = ' '
      ENDIF

!     Remove any excess blanks
      CALL STRIP(number)

      RETURN
      END SUBROUTINE FWWDD
! 
!============================================================================
!
      SUBROUTINE GetNum(STR,PC, IC, DNUM, IFLG)
!     Get the first DOUBLE PRECISION number in str begining at position
!     STR(ic:) delimited by the characters in PC.  A reasonable choice
!     for PC = PC = '= ,'.  After reading the number into DNUM, IC is
!     incremented to the next position in str (i.e. where the next number
!     may start).  IFLG is returned with the following values
!
!      -2: A number could not be read into dnum and no more data expected
!          (line doesn't end with ampersand).
!      -1: A number could not be read into dnum and additional data
!          for this variable is expected on the next line.
!       0: A number was read into dnum and str may contain additional
!          data.
!       > 0 (error)
!
!     NOTE: if iflg < 0, then IC is returned equal to 1.
!
!     DATE:   MAR. 6, 1996
!     AUTHOR: R.D. STEWART
!     REVISIONS: MAR. 10, 1996 (bug fix)
!     COMMENTS:
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER IC,IFLG,ilen
      DOUBLE PRECISION DNUM
      CHARACTER*(*) PC,STR
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

!     Use the CHARPAK routine DEXTNUM
      CALL TrimStr('!',str)
      ilen=LENTRIM(str)
      IF (ilen.GT.0) THEN
        CALL DEXTNUM(IC,STR,PC,DNUM,IFLG)
        IF (iflg.NE.0) THEN
!         No more data on this line.
          IF (str(ilen:ilen).EQ.'&') THEN
!           Continue on next line.
            iflg=-1
          ELSE
            iflg=-2
          ENDIF
        ELSE
          IFLG=0
        ENDIF
      ELSE
        IFLG=-2
      ENDIF

      RETURN
      END SUBROUTINE GetNum
! 
!============================================================================
!
      SUBROUTINE IExtNum(col,txt,lst,inum,iflg)
!
!     This routine searches through txt beginning at column,
!     col, and extracts a section of text delimited by
!     the string of characters in lst.  After a string of
!     characters is selected, an attempt is made to convert
!     them to an integer number, inum.  A status flag
!     is returned with the following meaning(s)
!
!        iflg <= 0 no errors
!        iflg > 0 an error occured
!
!     COMMENTS: The work space only needs to be as large as
!               the largest string of characters representing
!               a number.  50 characters should be sufficient
!               for almost any practical situation.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JAN 23, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                    CHARPAK.FOR
!        STRIP                      CHARPAK.FOR
!        TxtExt                     CHARPAKB.FOR
!        Txt2DNum                   CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col,iflg
      CHARACTER*(*) txt,lst
      INTEGER inum
      CHARACTER*50 work,fmat*15
      INTEGER ilen !,LENTRIM
!      EXTERNAL LENTRIM

      IF (lst.EQ.'\D') THEN
!       Use default delimiters for numeric data.
        CALL TxtExt(' ,=()',txt,col,work,iflg)
      ELSE
!       User specified list.
        CALL TxtExt(lst,txt,col,work,iflg)
      ENDIF

      IF (iflg.LE.0) THEN
!       Now try to convert the text to a double precision
!       number.
        CALL STRIP(work)
        ilen=LENTRIM(work)
        CALL MakFmt('I',ilen,0,fmat,iflg)
        IF (iflg.LE.0) THEN
          CALL Txt2Inum(fmat,work,inum,iflg)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE IExtNum
! 
!============================================================================
!
      SUBROUTINE IFMAT(col,inum,fspec,text)
!
!     PURPOSE: Format an integer number (inum) and
!              store it as a character string beginning at
!              column COL in TEXT.  FSPEC is an adaptive
!              (or "smart") format specifier.
!
!      INPUTS:  col - a valid integer in the range 1 to LEN(text)
!             inum - an integer number
!
!       FSPEC: see comments in DFMAT
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 23, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        TRANU!                     CHARPAK.FOR
!        LENTRIM                    CHARPAK.FOR
!        CopyTxt                    CHARPAK.FOR
!        FWWDD                      CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col
      INTEGER inum,ic,iw  !,LENTRIM
      CHARACTER*(*) fspec,text
      CHARACTER*1 achr,number*35
      DOUBLE PRECISION dnum
!      EXTERNAL LENTRIM

      dnum = inum
      call tranuc(fspec)

!     For integers, a "sensible" or limited set of format specifiers
!     are allowed.
      ic = INDEX(fspec,'E')
      achr = 'I'
      IF (ic.GT.0) THEN
        fspec(ic:ic) = 'I'
        achr = 'E'
      ELSE
        ic = INDEX(fspec,'F')
        IF (ic.GT.0) THEN
          fspec(ic:ic) = 'I'
          achr = 'F'
        ENDIF
      ENDIF
      CALL FWWDD(dnum,fspec,number)

!     Restore the orignial format specifier.
      ic = INDEX(fspec,'I')
      fspec(ic:ic) = achr

      iw = LENTRIM(number)
      CALL CopyTxt(col,number(1:iw),text)

      RETURN
      END SUBROUTINE IFMAT
! 
!============================================================================
!
      FUNCTION IMagD( dnum ) RESULT( ix )
!
!    <IMagD> computes the integer magnitude of a
!     double precision number, dnum.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JUN 3, 1994
!      REVISION: MAR. 21, 1996 bug-fix for very small (or large) 
!                double precision numbers (e.g., 1D-209).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                  CHARPAK.FOR
!        SCANSTR                  CHARPAK.FOR
!        PGSTR                    CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      DOUBLE PRECISION dnum
      CHARACTER*8 text
      CHARACTER*5 fmat
      INTEGER ic,ix,iochk,ilen,ipos !,LENTRIM,ipos  !, imagd
!      EXTERNAL LENTRIM

      WRITE(text,'(1PE8.1)',iostat=iochk) ABS(dnum)
      CALL SCANSTR('DEde',text,ipos)
      IF (ipos.EQ.0) THEN
        CALL SCANSTR('+-',text,ipos)
        IF (ipos.EQ.0) THEN
          WRITE(*,'(1X,A,1X,G12.6)')
     &       '<IMagD> WARNING: Unable to compute exponent for ',dnum
! For F90 style function (bmy, 12/4/98)
!          IMAGD=0
	   ix = 0
          RETURN
        ENDIF
      ELSE
        ipos=ipos+1
      ENDIF

      ILEN=LENTRIM(text)
      IC=ILEN-IPOS+1
      CALL MakFmt('I',ic,0,fmat,iochk)
      IF (iochk.EQ.0) THEN
        READ(text(ipos:ilen),fmt=fmat,iostat=iochk) ix
        IF (iochk.NE.0) THEN
          WRITE(*,'(1x,A,1X,I6)') '<IMagD> ERROR CODE: ',iochk
!          WRITE(*,*) text
!          WRITE(*,*) ' ilen: ',ilen
!          WRITE(*,*) ' ipos: ',ipos
!          WRITE(*,*) '   ic: ',ic
!          WRITE(*,*) ' dnum: ',dnum
        ELSE
          !IMagD=ix
        ENDIF
      ELSE
        WRITE(*,'(1X,A,1X,I6)') '<IMagD> ERROR CODE: ',iochk
        WRITE(*,'(1X,A,1X,G12.6)')
     &       'WARNING: Cannot create format specifier for ',dnum
! For F90-style function (bmy, 12/4/98)
        !IMAGD=0
        ix = 0
      ENDIF

      RETURN
      END FUNCTION IMagD
! 
!============================================================================
!
      FUNCTION IMagI(inum) RESULT( ix )
!
!    <IMagI> computes the integer magnitude of an
!     integer number, inum.
!
!     COMMENTS: Just a stub for IMagD
!
!      AUTHOR: Robert D. Stewart
!        DATE: JUN 3, 1994
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        IMagD                      CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER inum, ix
!      INTEGER IMagD
!      EXTERNAL IMagD

! For the F90-style function (bmy, 12/7/98)
!      IMagI = IMagD(DBLE(inum))
      Ix    = IMagD(DBLE(inum))

      RETURN
      END FUNCTION IMagI
! 
!============================================================================
!
      FUNCTION IMagS(snum) RESULT( ix )
!
!    <IMagS> computes the integer magnitude of a
!     single precision number, snum.
!
!     COMMENTS: Just a stub for IMagD
!
!      AUTHOR: Robert D. Stewart
!        DATE: JUN 3, 1994
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        IMagD                      CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      REAL snum
      INTEGER ix
!      INTEGER IMagD
!      EXTERNAL IMagD

! For the F90-style function statement (bmy, 12/7/98)
!      IMagS = IMagD(DBLE(snum))
      ix    = IMagD(DBLE(snum))

      RETURN
      END FUNCTION IMagS
! 
!============================================================================
!
      SUBROUTINE InsTxt(col,Tins,Tout)
!
!     PURPOSE: Insert the characters in Tins into Tout beginning at
!              column, col, in Tout.  Characters already existing
!              are shifted to the right in TEXT.  They may also
!              be shifted out the end and lost if the inserted
!              character string is too large.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        ShtTxt                    CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) Tins,Tout
      INTEGER col,LI,ic,LT

      LI = LEN(Tins)
      LT = LEN(Tout)
      ic = MAX0(col,1)
      CALL ShtTxt(LI,Tout(ic:LT))
      CALL CopyTxt(ic,Tins,Tout)

      RETURN
      END SUBROUTINE InsTxt
! 
!============================================================================
!
      SUBROUTINE INum2Txt(inum,icol,fmat,txt,iflg)
!
!     <INum2Txt> converts an integer number into
!     a string of equivalent characters and stores them
!     in the character variable Txt beginning at column,
!     icol.  The number is formated according to the
!     the standard FORTRAN specifier passed in fmat.
!     iflg is returned as a status flag indicating the success
!     or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER inum
      INTEGER icol,iflg
      CHARACTER*(*) txt,fmat

      WRITE(txt(icol:),fmt=fmat,iostat=iflg) inum

      RETURN
      END SUBROUTINE INum2Txt
! 
!============================================================================
!
      LOGICAL FUNCTION IsAlpCh(ch)
!
!     Is ch an alphanumeric character (i.e., an upper or lowercase
!     letter from A to Z or a numger from 0 to 9).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 ch
      LOGICAL Ok
      INTEGER iasc

      iasc = ICHAR(ch)
      ok = .false.
      IF ((iasc.GE.65).AND.(iasc.LE.90)) THEN
!       Uppercase Character
        ok = .true.
      ELSEIF ((iasc.GE.97).AND.(iasc.LE.122)) THEN
!       Lowercase Character
        ok = .true.
      ELSEIF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
!       Number
        ok = .true.
      ENDIF
      IsAlpCh = Ok

      RETURN
      END FUNCTION IsAlpCh
! 
!============================================================================
!
      LOGICAL FUNCTION IsChar(ch)
!
!     Is ch a character (i.e., an upper or lowercase
!     letter from A to Z).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 ch
      LOGICAL Ok
      INTEGER iasc

      iasc = ICHAR(ch)
      ok = .false.
      IF ((iasc.GE.65).AND.(iasc.LE.90)) THEN
!       Uppercase Character
        ok = .true.
      ELSEIF ((iasc.GE.97).AND.(iasc.LE.122)) THEN
!       Lowercase Character
        ok = .true.
      ENDIF
      IsChar = Ok

      RETURN
      END FUNCTION IsChar
! 
!============================================================================
!
      FUNCTION IsNum( ch ) RESULT ( lnum )
!
!     Returned as true if ch is a numeric character (i.e., one of
!     the numbers from 0 to 9).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   NOV. 11, 1993
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 ch
      INTEGER iasc  !, isnum
      LOGICAL lnum

      iasc = ICHAR(ch)
      lnum = .FALSE.
      IF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
        lnum = .TRUE.
      ENDIF
      !IsNum = lnum  !(bmy, 12/4/98)

      RETURN
      END FUNCTION IsNum
! 
!============================================================================
!
      LOGICAL FUNCTION IsPrnCh(ch)
!
!     Is the character in ch printable (i.e., ASCII code
!     from 33 to 255).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 ch
      LOGICAL Ok
      INTEGER iasc

      iasc = ICHAR(ch)
      ok = .true.
      IF (iasc.LT.33) THEN
        ok = .false.
      ENDIF
      IsPrnCh = Ok

      RETURN
      END FUNCTION IsPrnCh
! 
!============================================================================
!
      SUBROUTINE KeyExt(Nblk,key,ch,text,word,iflg)
!
!     PURPOSE: KeyExt provides an extremely flexible way to extract
!              parameters associated with a keyword in a
!              NAMELIST-like manner.  For example, if text='INP=1.00',
!              the number 1.00 can be extracted by the following
!
!                 CALL KeyExt(1,'INP','=',text,word,iflg)
!
!              "Word" would be returned with '1.00' in it.  ch can
!              can contain one or more delimiting characters.
!              (see TxtExt for more information).
!
!              Nblk is the number of "blocks" (i.e., sections of
!              text separated by characters in ch) to extract
!              from text.  If Nblk>1, each block of text in
!              word will be seperated by the first character
!              in ch.
!
!              Although they must be unique, partial keywords can
!              also be used.  In the above example, identical
!              results would be obtained from KeyExt if
!              text='INPUT = 1.00' instead of text='INP=1.00'
!              KeyExt IS keyword case sensitive (i.e., INPUT<>Input).
!
!              A status flag is also returned with the following
!              meaning(s)
!
!                 IF iflg = 0, task completed sucessfully
!                 IF iflg = 1, keyword not found
!                    iflg = 2, Insufficient data in TEXT
!                    iflg = 3, out of word storage-space
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CntMat                    CHARPAKA.FOR
!        TxtExt                    CHARPAKB.FOR
!        FillStr                   CHARPAKA.FOR
!        other routines (indirectly) at a lower level
!
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 1st, 1995
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------

      CHARACTER*(*) key,ch,text,word
      INTEGER Nblk,ikey,iflg,col,ic1
      INTEGER Kmax,Wmax,ilen,i,Tmax,imat
      LOGICAL again,more

!     Check for keyword
      ikey=INDEX(text,key)
      IF (ikey.GT.0) THEN
        Kmax = LEN(key)
        Tmax = LEN(text)

!       Extract the next sequence of characters delimited
!       by characters passed in ch
        ic1 = ikey+Kmax

!       Move ic1 to the location of the next character matching
!       a character in ch.
        again = (ic1.LE.Tmax)
        DO WHILE (again)
          CALL CntMat(text(ic1:ic1),ch,imat)
          IF (imat.GT.0) THEN
            again=.false.
            more = .true.
          ELSEIF (ic1.LT.Tmax) THEN
            ic1=ic1+1
          ELSE
!           No matching delimiters between keyword and associated
!           text.  Set ic1 to key+kmax
            more = .false.
            again = .false.
            ic1=ikey+Kmax
          ENDIF
        ENDDO

!       Now extract blocks of data
        col = ic1+1
        Wmax=LEN(Word)
        ilen = 1
        i = 1
        again = (more.AND.(i.LE.Nblk))
        DO WHILE (again)
          CALL TxtExt(ch,text,col,word(ilen:Wmax),iflg)
          IF (iflg.EQ.1) THEN
!           Insufficient Data in text to continue to read current
!           block.
            again=.false.
            iflg = 2
          ELSEIF (i.LT.Nblk) THEN
!           Increment counters and word-length
            i = i+1
            ilen = ilen + (col-ic1)
            ic1=col
            IF (iflg.LT.0) THEN
!             Insufficient Data in text to read next block
              again=.false.
              iflg = 2
            ELSEIF (ilen.GE.Wmax) THEN
!             Insufficient word-space to continume
              again=.false.
              iflg=3
            ELSE
!             Insert Space charcter between blocks selected
!             from the first character in ch.  NOTE: iflg
!             is used as a temporary variable here.
              iflg=ilen-1
              Word(iflg:iflg) = ch(1:1)
            ENDIF
          ELSE
!           Done
            again=.false.
            iflg=0
          ENDIF
        ENDDO
      ELSE
!       Keyword not found
        CALL FillStr(' ',Word)
        iflg = 1
      ENDIF

      RETURN
      END SUBROUTINE KeyExt
! 
!============================================================================
!
      FUNCTION LENTRIM( text ) RESULT ( nchr )
!
! Wrapper for the F90 intrinsic function "LEN_TRIM"
! Used for compatibility with the rest of the charpak functions
! (bmy, 12/8/98)

      IMPLICIT NONE   

      INTEGER           nchr
      CHARACTER (LEN=*) text

      nchr = LEN_TRIM( text )

      RETURN
      END FUNCTION LENTRIM
! 
!============================================================================
!
      SUBROUTINE LetSub(ch,sub,Text)
!
!     PURPOSE: LetSub can be used to replace a set of characters
!              in text with a the specific character passed in
!              sub.  The set of characters to be replaced are
!              listed in ch.  For example, if ch=' =' and
!              sub='~', all th equal signs and blanks in text
!              will be replaced upon return with a ~ (squiggle).
!
!     Author: R.D. Stewart
!     Date:   Jan. 6, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CntMat                    CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 sub
      CHARACTER*(*) ch,Text
      INTEGER Tmax,i,imat

      Tmax=LEN(Text)
      DO i=1,Tmax
        CALL CntMat(Text(i:i),ch,imat)
        IF (imat.GT.0) THEN
          Text(i:i) = sub
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE LetSub
! 
!============================================================================
!
      SUBROUTINE LINE(ic1,ic2,style,text)
!
!     Generate a line from column ic1 to ic2 in the character
!     variable text.  The line style is determined
!     the the character in style.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   JUL 1, 1993
!     AUTHOR: R.D. STEWART
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER ic1,ic2
      CHARACTER*(*) text
      CHARACTER*1 style
      INTEGER i,j1,j2,ilen

      ilen = LEN(text)
      j1 = MAX0(1,ic1)
      j2 = MIN0(ilen,ic2)

      DO i=j1,j2
        text(i:i) = style
      ENDDO

      RETURN
      END SUBROUTINE LINE
! 
!============================================================================
!
      SUBROUTINE MakFmt(type,iw,id,fmat,iochk)
!
!     Make a valid format specifier of type, type, with a field
!     width, iw and fractional part id digits wide. iochk is
!     returned <= zero, if the format specifier was created
!     sucessfully.
!
!     The work space is dimensioned so that it can handle
!     virtually any practical (conceivable?) iw and id.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JAN 23, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                    CHARPAK.FOR
!        STRIP                      CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) fmat,type
      INTEGER iw,id
      INTEGER ilen,itmp,iochk,i
      CHARACTER*25 work
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

      itmp=LENTRIM(type)
      IF (itmp.GT.0) THEN
        CALL TRANUC(type)
        IF (INDEX(type,'I').GT.0) THEN
          WRITE(work,fmt='(a,a,I5,a)',iostat=iochk)
     @          '(',type(1:itmp),iw,')'
        ELSE
          WRITE(work,fmt='(a,a,I5,a,I5,a)',iostat=iochk)
     @          '(',type(1:itmp),iw,'.',id,')'
        ENDIF
      ELSE
!       Invalid type-- default to F
        type='F'
        WRITE(work,fmt='(a,a,I5,a,I5,a)',iostat=iochk)
     @        '(','F',iw,'.',id,')'
      ENDIF

      IF (iochk.EQ.0) THEN
!       Remove extra blanks and attempt to transfer to fmat
        itmp=LEN(fmat)
        CALL STRIP(work)
        ilen=LENTRIM(work)
        IF (ilen.GT.itmp) THEN
          iochk=2
        ELSE
          fmat(1:ilen) = work(1:ilen)
          DO i=ilen+1,itmp
            fmat(i:i) = ' '
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE MakFmt
! 
!============================================================================
!
      SUBROUTINE PGSTR(str)
!     PGSTR (PigLat String) swaps charcters in str so that, if str is a 
!     word, str
!     is returned as the "Pig Latin" equivalent of the original str.
!     Trailing blanks in str are ignored.
!     NOTE: the original string is destroyed on exit, namely it is
!           replaced with the PigLat equivalent of the original string.
!
!      AUTHOR: Robert D. Stewart
!        DATE: FEB. 22, 1996
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!---------------------------------------------------------------------------- 
      CHARACTER*(*) str
      INTEGER ilen,I,J,imid !,LENTRIM
!      EXTERNAL LENTRIM

      ilen= LENTRIM(str)
      IF (ilen.GT.0) THEN
        imid = ilen/2
        DO i=1,imid
          j = ilen-i + 1
          CALL CSWAP(i,j,str)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE PGSTR
! 
!============================================================================
!
      SUBROUTINE PigLat(str)
!
!     Convert a string of words (groups of letters separated by
!     blanks or various puncutation marks) into "Pig Latin"
!     Words should contain less than 80 characters.
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) str
      CHARACTER*80 word,pg
      INTEGER col,iflg,i1,i2,ilen !,LENTRIM
!      EXTERNAL LENTRIM

      iflg=0
      col = 1
      DO WHILE (iflg.EQ.0)
        word = ' '
        CALL TxtExt('":;,. !?''''',str,col,word,iflg)
        IF (iflg.LE.0) THEN
          pg = ' '
          pg = word
          CALL PGSTR(pg)
          ilen=LENTRIM(word)
          i1 = INDEX(str,word(:ilen))
          i2 = i1 + ilen-1
          str(i1:i2) = pg(:ilen)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE PigLat
! 
!============================================================================
!
      SUBROUTINE PLINE(io,ic1,ic2,style)
!
!     Print a line from column ic1 to ic2 on device io.
!     NOTE: the LINE routine can also be used to generate
!           lines.  The current maximum line length is
!           132 columns.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JAN 23, 1993

!     CODE DEPENDENCIES:
!      Routine Name                  File
!        FillStr                   CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER io,ic1,ic2,itmp
      CHARACTER*132 text
      CHARACTER*1 style

      itmp = MIN0(ic2,132)
      IF (itmp.GT.ic1) THEN
        IF (ic1.Gt.0) THEN
          CALL FillStr(' ',text)
          CALL FillStr(style,text(ic1:itmp))
          WRITE(io,'(A)') text(1:itmp)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PLINE
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! replicate: create a string with n occurences of ch
! difference to IDL: not a function, string must exist before
      SUBROUTINE REPLICATE(ch,count,text)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1   ch
      INTEGER       count
      CHARACTER*(*) text
!      EXTERNAL      FILLSTR, LINE

      CALL FILLSTR(' ',text)
      CALL LINE(1,count,ch,text)

      RETURN
      END SUBROUTINE REPLICATE
! 
!============================================================================
!
      SUBROUTINE RmvRep(ch,sub,Text,ilen)
!
!     RmvRep replaces multiple occurances of the character, ch,
!     in Text with a single occurance of the character in sub.
!     The modified string is returned in Text left justified.
!     ilen is returned as the length of the modified string
!     minus trailing blanks.
!
!     COMMENTS: ilen will always be less than or equal to the
!               original length of text, and the original
!               string will be destroyed/modified.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!        AUTHOR: Robert D. Stewart
!          DATE: Jan. 2, 1995
!     REVISIONS: Mar. 6, 1996 bug fix (ilen not returned as length of
!                text minus trailing blanks).
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 sub,ch,buf
      CHARACTER*(*) Text
      INTEGER Tmax,ilen,i !,LENTRIM
      LOGICAL srt
!      EXTERNAL LENTRIM

!     Compute length of string
      Tmax=LEN(Text)

!     Srt indicates the start of a replacement sequence.
      srt=.true.
      ilen = 0
      DO i=1,Tmax
        IF (Text(i:i).EQ.ch) THEN
          IF (srt) THEN
!           Perform (macro) substitution.
            srt=.false.
            ilen = ilen+1
            Text(ilen:ilen) = sub
          ENDIF
        ELSE
!         Add current character to Text at position ilen.
!         srt is set  true so that the next occurance of a
!         character in ch will begin the removal procedure.
!         NOTE: buf is used to avoid self-referencing a
!         character variable a no-no in F77 (it can also
!         improve performance on some computers.
          ilen = ilen+1
          buf = Text(i:i)
          Text(ilen:ilen) = buf
          srt = .true.
        ENDIF
      ENDDO

!     Fill remainder of Text with blanks.
      DO i=ilen+1,Tmax
        Text(i:i) = ' '
      ENDDO
      ilen=LENTRIM(text)

      RETURN
      END SUBROUTINE RmvRep
! 
!============================================================================
!
      SUBROUTINE RmvTxt(ch,Text,ilen)
!
!     RmvTxt removes all occurances of the character, ch,
!     in Text.
!
!     The modified string is returned in Text left justified.
!     ilen is returned as the length of the modified string
!     minus trailing blanks.
!
!     COMMENTS: ilen will always be less than or equal to the
!               original length of text, and the original
!               string will be destroyed/modified.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!        AUTHOR: Robert D. Stewart
!          DATE: Jan. 2, 1995
!     REVISIONS: Mar. 6, 1996 bug fix (ilen not returned as length of
!                text minus trailing blanks).
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*1 ch,buf
      CHARACTER*(*) Text
      INTEGER Tmax,ilen,i !,LENTRIM
!      EXTERNAL LENTRIM

!     Compute length of string
      Tmax=LEN(Text)

      ilen = 0
      DO i=1,Tmax
        IF (Text(i:i).NE.ch) THEN
!         Add current character to Text at position ilen.
!         NOTE: buf is used to avoid self-referencing a
!         character variable (a no-no in FORTRAN77)
          ilen = ilen+1
          buf = Text(i:i)
          Text(ilen:ilen) = buf
        ENDIF
      ENDDO

!     Fill remainder of Text with blanks.
      DO i=ilen+1,Tmax
        Text(i:i) = ' '
      ENDDO
      ilen=LENTRIM(text)

      RETURN
      END SUBROUTINE RmvTxt
! 
!============================================================================
!
      SUBROUTINE SCANSTR(ch,str,ipos)
!
!     Return the position (as ipos) of the first character in ch
!     that matches a character in str.  If no matching characters
!     are found, ipos is returned equal to zero.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!          CntMat                   CharPak
!
!     DATE:   MAR. 6, 1996
!     AUTHOR: R.D. STEWART
!     COMMENTS:
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) ch,str
      INTEGER ipos,imax,imat
      LOGICAL again

      ipos=0
      imax=LEN(str)
      again=(ipos.LE.imax)
      IF (again) ipos=1
      DO WHILE (again)
        CALL CntMat(ch,str(ipos:ipos),imat)
        IF (imat.GT.0) THEN
          again=.FALSE.
        ELSEIF (ipos.LT.imax) THEN
          ipos=ipos+1
        ELSE
          again=.FALSE.
          ipos=0
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SCANSTR
! 
!============================================================================
!
      SUBROUTINE SExtNum(col,txt,lst,snum,iflg)
!
!     This routine searches through txt beginning at column,
!     col, and extracts a section of text delimited by
!     the string of characters in lst.  After a string of
!     characters is selected, an attempt is made to convert
!     them to a single precision number, snum.  A status flag
!     is returned with the following meaning(s)
!
!        iflg <= 0 no errors
!        iflg > 0 an error occured
!
!     COMMENTS: The work space only needs to be as large as
!               the largest string of characters representing
!               a number.  50 characters should be sufficient
!               for almost any practical situation.
!
!      AUTHOR: Robert D. Stewart
!        DATE: JAN 23, 1995
!
!     REVISIONS: 04-FEB-1996 snum was incorrectly declared as
!                an integer.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                    CHARPAK.FOR
!        STRIP                      CHARPAK.FOR
!        TxtExt                     CHARPAK.FOR
!        Txt2SNum                   CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col,iflg
      CHARACTER*(*) txt,lst
      REAL snum
      CHARACTER*50 work,fmat*15
      INTEGER ilen !,LENTRIM
!      EXTERNAL LENTRIM

      IF (lst.EQ.'\D') THEN
!       Use default delimiters for numeric data.
        CALL TxtExt(' ,=()',txt,col,work,iflg)
      ELSE
!       User specified list.
        CALL TxtExt(lst,txt,col,work,iflg)
      ENDIF

      IF (iflg.LE.0) THEN
!       Now try to convert the text to a double precision
!       number.
        CALL STRIP(work)
        ilen=LENTRIM(work)
        CALL MakFmt('F',ilen,0,fmat,iflg)
        IF (iflg.LE.0) THEN
          CALL Txt2Snum(fmat,work,snum,iflg)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE SExtNum
! 
!============================================================================
!
      SUBROUTINE SFMAT(col,snum,fspec,text)
!
!     PURPOSE: Format a single precision number (snum) and
!              store it as a character string beginning at
!              column COL in TEXT.  FSPEC is an adaptive
!              (or "smart") format specifier.
!
!      INPUTS:  col - a valid integer in the range 1 to LEN(text)
!             dnum - a single precision number
!
!       FSPEC: see comments in DFMAT
!
!      COMMENTS: This routine is just a "stub" or wrapper
!                function for the DFMAT routine.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 23, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        DFMAT                      CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER col
      REAL snum
      CHARACTER*(*) fspec,text
      DOUBLE PRECISION dnum

      dnum = snum
      CALL DFMAT(col,dnum,fspec,text)

      RETURN
      END SUBROUTINE SFMAT
! 
!============================================================================
!
      SUBROUTINE ShtTxt(ispc,Text)
!
!     PURPOSE: Shift the characters in TEXT ABS(ispc) spaces
!              to the left (ispc<0) or right (ispc>0).  Blanks
!              are added at the left (or right) as padding for
!              the newly "exposed" location in TEXT.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 2, 1995
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      INTEGER ispc
      INTEGER num,ilen,i,j
      CHARACTER*1 ch

      ilen = LEN(text)
      num = ABS(ispc)
      IF (ispc.LT.0) THEN
!       Shift characters to the left
        DO i=1,num
          DO j=2,ilen
!           NOTE: ch is used as to directly avoid self-referencing
!           a character variable.
            ch = text(j:j)
            text(j-1:j-1) = ch
          ENDDO
          text(ilen:ilen) = ' '
        ENDDO
      ELSE
!       Shift characters to the right
        DO i=1,num
          DO j=ilen-1,1,-1
            ch = text(j:j)
            text(j+1:j+1) = ch
          ENDDO
          text(1:1) = ' '
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE ShtTxt
! 
!============================================================================
!
      SUBROUTINE SNum2Txt(snum,icol,fmat,txt,iflg)
!
!     <SNum2Txt> converts a single precision number into
!     a string of equivalent characters and stores them
!     in the character variable Txt beginning at column,
!     icol.  The number is formated according to the
!     the standard FORTRAN specifier passed in fmat.
!     iflg is returned as a status flag indicating the success
!     or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      REAL snum
      INTEGER icol,iflg
      CHARACTER*(*) txt,fmat

      WRITE(txt(icol:),fmt=fmat,iostat=iflg) snum

      RETURN
      END SUBROUTINE SNum2Txt
! 
!============================================================================
!
      SUBROUTINE SQUEEZE(text)
!
!     PURPOSE: Remove any extra and leading blanks from character
!              variable text.
!
!    COMMENTS: The original "text" is destroyed on exit.  This
!              routine is superseded by the RmvRep routines, but
!              is included for backwards compatibility.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      CHARACTER*1 ch
      INTEGER iasc,icnt,ilen,i
      LOGICAL space

      ilen = LEN(text)
      IF (ilen.GT.1) THEN
        icnt = 1
        space = .false.
        DO I=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
            SPACE = .true.
          ELSEIF (space) THEN
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
            space = .false.
          ENDIF
        ENDDO
        DO I=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SQUEEZE
!
!============================================================================
!
      SUBROUTINE SRndDig(sin,nd,sout)
!
!     <SRndDig> rounds a single precision number sin to the
!     number of significant digits indicated by integer paramter
!     nd.  If nd < 0, rounding will be performed to the right
!     of the decimal place.  For nd > 0, rounding is performed to
!     the left of the decimal place.  The rounded number is
!     returned in sout.
!
!      AUTHOR: Robert D. Stewart
!        DATE: FEB 13, 1993
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      REAL sin,sout
      INTEGER nd
      DOUBLE PRECISION xmult

!     Create multiplier for rounding
      xmult=10.0D+00**(nd)

!     Round number.
      sout = DNINT(sin*xmult)/xmult

      RETURN
      END SUBROUTINE SRndDig
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strcompress: remove white spaces (blank or non-printable) from string
! difference to IDL: not a function, removeall flag must always be passed
      SUBROUTINE STRCOMPRESS(text,removeall)
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      LOGICAL       removeall

!      EXTERNAL      SQUEEZE,STRIP


      IF(removeall) THEN
         CALL strip(text) 
      ELSE 
         CALL squeeze(text)
      ENDIF

      RETURN
      END SUBROUTINE STRCOMPRESS
! 
!============================================================================
!
      SUBROUTINE STRIP(text)
!
!     PURPOSE: Strip blanks and null characters for the variable TEXT.
!
!     COMMENTS: The original "text" is destroyed upon exit.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) TEXT
      INTEGER ilen,iasc,icnt,i
      CHARACTER*1 ch

      ilen = LEN(text)
      IF (ilen.GT.1) THEN
        icnt = 1
        DO i=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
!           Keep character
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
!       Fill remainder of text with blanks
        DO i=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE STRIP
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strlen: length of a string 
      FUNCTION STRLEN( text ) RESULT ( ILEN )
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      INTEGER       ilen
!      INTEGER LENTRIM
!      EXTERNAL LENTRIM

! For F90-style function declaration
      !STRLEN = LENTRIM(text)
      ILEN  = LENTRIM( text )

      RETURN
      END FUNCTION STRLEN
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strlowcase: convert to lower case (identical to tranlc)
      SUBROUTINE STRLOWCASE(text)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
!      EXTERNAL      TRANLC

      CALL TRANLC(text)


      RETURN
      END SUBROUTINE STRLOWCASE
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strpos: return (first!) position of substring in string
! As in IDL, -1 is returned if expression not found
! difference to IDL: optional pos parameter not supported
      INTEGER FUNCTION STRPOS(expr,searchstr)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*)  expr,searchstr
      CHARACTER*4096 buf1, buf2
      INTEGER        pos, dum
      LOGICAL        ok
   
!      INTEGER        STRLEN 
!      EXTERNAL       TXTLOC, STRLEN


      strpos = -1   ! fail safe


      if (STRLEN(expr) .EQ. 0) RETURN
      if (STRLEN(searchstr) .EQ. 0) RETURN
      if (LEN(expr).GT.4096 .OR. LEN(searchstr).GT.4096) RETURN 
     
      CALL TxtLoc(expr,searchstr,buf1,buf2,.FALSE.,.TRUE.,pos,dum,ok)

      IF(pos.LE.0)pos=-1

      strpos = pos 

      RETURN
      END FUNCTION STRPOS
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strput: overwrite substring in dest with source
! (identical to COPYTXT, but different arrangement of parmeters)
      SUBROUTINE STRPUT(dest,source,pos)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) dest, source
      INTEGER       pos,ilen,slen

!      INTEGER       STRLEN
!      EXTERNAL      STRLEN

      ilen = LEN(dest)

      if (pos.GT.ilen .OR. pos.LE.0) RETURN

      slen = MIN(ilen-pos,STRLEN(source))
      dest(pos:pos+slen-1) = source(1:slen)


      RETURN
      END SUBROUTINE STRPUT
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strtrim: remove leading blanks and convert all trailing blanks to
! true spaces
! difference to IDL: not a function, no flag parameter, since there is
! no way to actually shorten the string in FORTRAN
! (adapted from SQUEEZE)
      SUBROUTINE STRTRIM(text)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*)    text
      CHARACTER*1      ch
      INTEGER          iasc,icnt,ilen,i
      LOGICAL          space


      ilen = LEN(text)

      IF (ilen.GT.1) THEN
        icnt = 1
        space = .false.
        DO I=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
            SPACE = .true.
          ELSEIF (space) THEN
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
        DO I=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE STRTRIM
! 
!============================================================================
!
!  mgs interface: a small set of wrapper routines with names and 
!     parameters as close to IDL as possible
!
!  17 Jun 1998
!
! strupcase: convert to upper case characters (identical to tranuc)
      SUBROUTINE STRUPCASE(text)

!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
!      EXTERNAL      TRANUC

      CALL TRANUC(text)

      RETURN
      END SUBROUTINE STRUPCASE
! 
!============================================================================
!
      SUBROUTINE TRANLC(text)
!
!     PURPOSE: Tranlate a character variable to all lowercase letters.
!              Non-alphabetic characters are not affected.
!
!    COMMENTS: The original "text" is destroyed.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO I=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
          text(i:i) = CHAR(iasc+32)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE TRANLC
! 
!============================================================================
!
!
      SUBROUTINE TRANUC(text)
!
!     PURPOSE: Tranlate a character variable to all upper case letters.
!              Non-alphabetic characters are not affected.
!
!    COMMENTS: The original "text" is destroyed.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO i=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE TRANUC
! 
!============================================================================
!
      SUBROUTINE TrimStr(ch,str)
!
!     Trim (delete) all characters to the right of a specified
!     character, ch.
!
!      AUTHOR: Robert D. Stewart
!        DATE: MAR. 6, 1996
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) str,ch*1
      INTEGER ix

      ix=INDEX(str,ch)
      IF (ix.GT.0) THEN
        CALL FillStr(' ',str(ix:))
      ENDIF

      RETURN
      END SUBROUTINE TrimStr
! 
!============================================================================
!
      SUBROUTINE Txt2Dnum(fmat,txt,dnum,iflg)
!     <Txt2Dnum> attempts to convert the string of characters
!     in txt into a double precision number.  fmat is the
!     VALID format specifier to use in the internal read
!     statement.  iflg is returned as a status flag indicating
!     the success or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!     COMMENTS: Generally, the Fxx.0 format should be used to convert
!               string of characters to a number.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) txt,fmat
      DOUBLE PRECISION dnum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) dnum

      RETURN
      END SUBROUTINE Txt2Dnum
! 
!============================================================================
!
      SUBROUTINE Txt2Inum(fmat,txt,Inum,iflg)
!     <Txt2Inum> attempts to convert the string of characters
!     in txt into a integer number.  fmat is the
!     VALID format specifier to use in the internal read
!     statement.  iflg is returned as a status flag indicating
!     the success or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!     COMMENTS: Generally, the Fxx.0 format should be used to convert
!               string of characters to a number.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) txt,fmat
      INTEGER inum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) inum

      RETURN
      END SUBROUTINE Txt2Inum
! 
!============================================================================
!
      SUBROUTINE Txt2Snum(fmat,txt,snum,iflg)
!     <Txt2Snum> attempts to convert the string of characters
!     in txt into a single precision number.  fmat is the
!     VALID format specifier to use in the internal read
!     statement.  iflg is returned as a status flag indicating
!     the success or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!     COMMENTS: Generally, the Fxx.0 format should be used to convert
!               string of characters to a number.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) txt,fmat
      REAL snum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) snum

      RETURN
      END SUBROUTINE Txt2Snum
! 
!============================================================================
!
      SUBROUTINE TxtExt(ch,text,col,word,iflg)
!
!     PURPOSE: TxtExt extracts a sequence of characters from
!              text and transfers them to word.  The extraction
!              procedure uses a set of character "delimiters"
!              to denote the desired sequence of characters.
!              For example if ch=' ', the first character sequence
!              bracketed by blank spaces will be returned in word.
!              The extraction procedure begins in column, col,
!              of TEXT.  If text(col:col) = ch (any character in
!              the character string), the text is returned beginning
!              with col+1 in text (i.e., the first match with ch
!              is ignored).
!
!              After completing the extraction, col is incremented to
!              the location of the first character following the
!              end of the extracted text.
!
!              A status flag is also returned with the following
!              meaning(s)
!
!                 IF iflg = -1, found a text block, but no more characters
!                               are available in TEXT
!                    iflg = 0, task completed sucessfully (normal term)
!                    iflg = 1, ran out of text before finding a block of
!                              text.
!
!       COMMENTS: TxtExt is short for Text Extraction.  This routine
!                 provides a set of powerful line-by-line
!                 text search and extraction capabilities in
!                 standard FORTRAN.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CntMat                    CHARPAK.FOR
!        TxtExt                    CHARPAK.FOR
!        FillStr                   CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR
!
!        other routines are indirectly called.
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 1st, 1995
!
!      REVISIONS: FEB 22, 1996.  Slight bug fix (introduced by a
!        (recent = FLIB 1.04) change in the CntMat routine)
!        so that TxtExt correctlyhandles groups of characters
!        delimited by blanks).
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) ch,text,word
      INTEGER col,iflg
      INTEGER Tmax,T1,T2,imat
      LOGICAL again,prev

!     Length of text
      Tmax = LEN(text)

!     Fill Word with blanks
      CALL FillStr(' ',Word)

      IF (col.GT.Tmax) THEN
!       Text does not contain any characters past Tmax.
!       Reset col to one and return flag = {error condition}
        iflg = 1
        col = 1
      ELSEIF (col.EQ.Tmax) THEN
!       End of TEXT reached
        CALL CntMat(ch,text(Tmax:Tmax),imat)
        IF (imat.EQ.0) THEN
!         Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
        ELSE
!         Same error condition as if col.GT.Tmax
          iflg = 1
        ENDIF
      ELSE
!       Make sure column is not less than 1
        IF (col.LT.1) col=1
        CALL CntMat(ch,text(col:col),imat)
        IF (imat.GT.0) THEN
          prev=.true.
        ELSE
          prev=.false.
        ENDIF
        T1=col
        T2 = T1

        again = .true.
        DO WHILE (again)
!         Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
!           Current character in TEXT matches one (or more) of the
!           characters in ch.
            IF (prev) THEN
              IF (T2.LT.Tmax) THEN
!               Keep searching for a block of text
                T2=T2+1
                T1=T2
              ELSE
!               Did not find any text blocks before running
!               out of characters in TEXT.
                again=.false.
                iflg=1
              ENDIF
            ELSE
!             Previous character did not match ch, so terminate.
!             NOTE: This is "NORMAL" termination of the loop
              again=.false.
              T2=T2-1
              iflg = 0
            ENDIF
          ELSEIF (T2.LT.Tmax) THEN
!           Add a letter to the current block of text
            prev = .false.
            T2=T2+1
          ELSE
!           Reached the end of the characters in TEXT before reaching
!           another delimiting character.  A text block was identified
!           however.
            again=.false.
            iflg=-1
          ENDIF
        ENDDO

        IF (iflg.EQ.0) THEN
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
        ELSE
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE TxtExt
! 
!============================================================================
!
      SUBROUTINE TxtLoc(Pat,Str,Pbuf,Sbuf,CS,RTB,c1,c2,ok)
!
!     TxtLoc (Text Location) search for a pattern of characters (Pat)
!     and returns the start/stop column of the pattern (c1:c2) occurring
!     in Str.  Pbuf and Sbuf are temporary workspaces for Pat and
!     Str, respectively and should be of greater or equal length.
!     If CS is true, the search is Case Sensitive, if CS
!     is false, the search is case insensitive (e.g., a
!     case insensitive search would mean Input=INPUT).
!     If RTB (remove trailing blanks) is true, any trailing
!     blanks in Pat are ignored.
!
!     Author: R.D. Stewart
!     Date:   Jan. 18, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        LENTRIM                   CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR
!        TRANL!                    CHARPAK.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) Pat,Str,Pbuf,Sbuf
      LOGICAL CS,RTB,ok
      INTEGER c1,c2,iplen !,LENTRIM
!      EXTERNAL LENTRIM

      IF (rtb) THEN
        iplen=LENTRIM(Pat)
      ELSE
        iplen=LEN(Pat)
      ENDIF

      IF (cs) THEN
!       Case Sensitive Case (easy)
        c1=INDEX(Str,Pat(1:iplen))
      ELSE
!       ***** Case Insensitive Search *****
!           (Quick and Dirty Algorithm)
        CALL CopyTxt(1,Pat(1:iplen),Pbuf)
        CALL CopyTxt(1,Str,Sbuf)
        CALL TRANLC(Pbuf(1:iplen))
        CALL TRANLC(Sbuf)
        c1=INDEX(Sbuf,Pbuf(1:iplen))
      ENDIF

      IF (c1.GT.0) THEN
        c2 = c1 + (iplen-1)
        ok = .true.
      ELSE
        ok=.false.
      ENDIF

      RETURN
      END SUBROUTINE TxtLoc
! 
!============================================================================
!
      SUBROUTINE TxtMac(cnop,CS,RTB,mac,rpl,Text,ilen)
!
!     TxtMac replaces cnop occurances of mac within text by
!     the string, rpl.
!     ilen is returned as the length of text minus trailing
!     blanks.  All string dimensions are dynamically determined
!     (i.e., by the calling program).  If cnop='*', all
!     occurances will be replaced.
!
!       CS = Case sensitive (true or false)
!      RTB = Remmove trailing blanks from macro replacement.
!         (See TxtLoc for more info on CS and RTB)
!
!     COMMENTS: This routine can be used as a primitive form of
!               macro substitution. NOTE: The algorithm could
!               use some additional work to improve performance.
!
!     Author: R.D. Stewart
!     Date:   Jan. 6, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CopyTxt                   CHARPAK.FOR   (indirect)
!        TRANLC                    CHARPAK.FOR   (indirect)
!        InsTxt                    CHARPAKA.FOR  (indirect)
!        ShtTxt                    CHARPAK.FOR   (indirect)
!        LENTRIM                   CHARPAK.FOR
!        TxtLoc                    CHARPAKA.FOR
!        TxtSub                    CHARPAKA.FOR
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      INTEGER nop,ilen
      CHARACTER*(*) cnop,mac,rpl,Text
      LOGICAL again,ok,CS,RTB
      INTEGER Tmax,Rmax,Mmax,ic
      INTEGER icnt,iochk,c1,c2 !,LENTRIM
!      EXTERNAL LENTRIM

!     Temporary WorkSpace.  A workspace of 255 characters should
!     cover most circumstances.
      CHARACTER*256 Tbuf,Mbuf
      INTEGER Blen
      PARAMETER (Blen=256)

!     Compute the length of the character strings
      Tmax=LEN(Text)
      Mmax=LEN(mac)
      IF (RTB) THEN
        Rmax=LENTRIM(rpl)
      ELSE
        Rmax=LEN(rpl)
      ENDIF
      ilen=MAX0(Tmax,Mmax)
      IF (ilen.GT.Blen) THEN
!       Print Warning Message about internal
!       buffer space.
        Tbuf='INTERNAL BUFFER WARNING: Insufficient workspace.'
        WRITE(*,'(1x,a79)') Tbuf
        Tbuf='Workspace available:            '
        WRITE(Tbuf(22:27),'(I6)') Blen
        WRITE(*,'(1x,a79)') Tbuf
        Tbuf='Workspace needed:            '
        WRITE(Tbuf(19:24),'(I6)') ilen
        WRITE(*,'(1x,a79)') Tbuf
        WRITE(*,'(1x,a)') 'You will need to re-compile <TxtMac>'
        PAUSE
      ENDIF

!     Compute Number of replacement ops
      IF (cnop.EQ.'*') THEN
!       All occurances.  NOTE: can't have more macro substitutions
!       than characters in Text.
        nop=Tmax
      ELSE
        READ(cnop,fmt='(I9)',iostat=iochk) nop
        IF (iochk.NE.0) THEN
!         Default to all occurances
          nop=Tmax
        ENDIF
      ENDIF

      CALL TxtLoc(mac,text,Mbuf,Tbuf,CS,RTB,c1,c2,ok)
      IF (ok) THEN
        icnt = 1
        ic=1
        again = (icnt.LE.nop)
      ELSE
        again=.false.
      ENDIF
      DO WHILE (again)
        CALL TxtSub(rpl,Rmax,Text(ic:Tmax),c1,c2)
        IF (icnt.LT.nop) THEN
          ic=ic+c1+Rmax-1
          CALL TxtLoc(mac,Text(ic:Tmax),Mbuf,Tbuf,CS,RTB,c1,c2,ok)
          IF (.not.ok) THEN
!           Nothing left to replace, so quit.
            again=.false.
          ELSE
            icnt=icnt+1
          ENDIF
        ELSE
!         Completed requested number of replacements. Done.
          again=.false.
        ENDIF
      ENDDO
      ilen=LENTRIM(Text)

      RETURN
      END SUBROUTINE TxtMac
! 
!============================================================================
!
      SUBROUTINE TxtSub(Tnew,Nlen,Told,c1,c2)
!
!      TxtSub (Text Substitute) replaces the string of
!      characters Told(c1:c2) with the string of characters
!      in Tnew(1:nlen).  If nlen is larger than (c1:c2),
!      characters in Tnew are shifted to the right.
!      If nlen is smaller than (c1:c2), characters in Tnew
!      are shifted to the left.
!
!     Author: R.D. Stewart
!     Date:   Jan. 18, 1995
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        InsTxt                    CHARPAKA.FOR
!        ShtTxt                    CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR (indirectly)
!
!----------------------------------------------------------------------------
! Also add the IMPLICIT NONE statement, for safety's sake (bmy, 12/3/98)
      IMPLICIT NONE   
!----------------------------------------------------------------------------
      CHARACTER*(*) Tnew,Told
      INTEGER Nlen,c1,c2,ix,diff
      INTEGER Tlen

      Tlen=c2-c1+1
      diff= Tlen - Nlen
      IF (diff.EQ.0) THEN
        Told(c1:c2) = Tnew(1:nlen)
      ELSEIF (diff.LT.0) THEN
!       Tnew is bigger than Told(c1:c2)
        Told(c1:c2)=Tnew(1:Tlen)
        CALL InsTxt(c2+1,Tnew(Tlen+1:nlen),Told)
      ELSEIF (diff.GT.0) THEN
!       Tnew is smaller than Told(c1:c2)
        ix=nlen+c1-1
        Told(c1:ix)=Tnew(1:nlen)
        CALL ShtTxt(-diff,Told(ix+1:))
      ENDIF

      RETURN
      END SUBROUTINE TxtSub
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: strrepl
!
! !DESCRIPTION: Subroutine STRREPL replaces all instances of PATTERN within
!  a string STR with replacement text REPLTXT. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE StrRepl( STR, PATTERN, REPLTXT )
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)    :: PATTERN   ! Search pattern
      CHARACTER(LEN=*), INTENT(IN)    :: REPLTXT   ! Replacement text
!
! !INPUT/OUTPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(INOUT) :: STR       ! String to be modified
!
! !REMARKS:
!     REPLTXT must have the same # of characters as PATTERN.
!
! !REVISION HISTORY: 
!     25 Jun 2002 - R. Yantosca - Initial version
!     20 Jul 2004 - R. Yantosca - Replace LEN_TRIM with LEN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                         :: I1, I2

      !=================================================================
      ! STRREPL begins here!
      !=================================================================

      ! Error check: make sure PATTERN and REPLTXT have the same # of chars
      IF ( LEN( PATTERN ) /= LEN( REPLTXT ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 
     &    'STRREPL: PATTERN and REPLTXT must have same # of characters!'
         WRITE( 6, '(a)' ) 'STOP in STRREPL (charpak_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      ! Loop over all instances of PATTERN in STR
      DO 

         ! I1 is the starting location of PATTERN w/in STR  
         I1 = INDEX( STR, PATTERN )

         ! If pattern is not found, then return to calling program
         IF ( I1 < 1 ) RETURN

         ! I2 is the ending location of PATTERN w/in STR
         I2 = I1 + LEN_TRIM( PATTERN ) - 1
      
         ! Replace text
         STR(I1:I2) = REPLTXT

      ENDDO
         
      ! Return to calling program
      END SUBROUTINE StrRepl
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymd_extract
!
! !DESCRIPTION: Subroutine YMD\_EXTRACT extracts the year, month, and date 
!  from an integer variable in YYYYMMDD format.  It can also extract the 
!  hours, minutes, and seconds from a variable in HHMMSS format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE YMD_EXTRACT( NYMD, Y, M, D )
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER, INTENT(IN)  :: NYMD   ! YYYY/MM/DD date or hh:mm:ss time
!
! !OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(OUT) :: Y      ! Year  (YYYY) or hour    (hh)
      INTEGER, INTENT(OUT) :: M      ! Month (MM)   or minutes (mm)
      INTEGER, INTENT(OUT) :: D      ! Day   (DD)   or seconds (ss)
!
! !REVISION HISTORY: 
!  21 Nov 2001 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      REAL*8               :: REM

      !=================================================================
      ! YMD_EXTRACT begins here!
      !=================================================================

      ! Extract YYYY from YYYYMMDD 
      Y = INT( DBLE( NYMD ) / 1d4 )

      ! Extract MM from YYYYMMDD
      REM = DBLE( NYMD ) - ( DBLE( Y ) * 1d4 )
      M   = INT( REM / 1d2 )

      ! Extract DD from YYYYMMDD
      REM = REM - ( DBLE( M ) * 1d2 )
      D   = INT( REM )

      ! Return to calling program
      END SUBROUTINE YMD_EXTRACT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: expandDate
!
! !DESCRIPTION: Subroutine expandDate replaces YYYYMMDD and hhmmss 
!  tokens within a filename string with the actual values. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE expandDate( FILENAME, YYYYMMDD, HHMMSS )
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)    :: YYYYMMDD   ! Date
      INTEGER,          INTENT(IN)    :: HHMMSS     ! Time
!
! !INPUT/OUTPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME   ! File name to be modified
!
! !REVISION HISTORY: 
!  (1 ) Bug fix for Linux: use ENCODE statement to convert number to string 
!        instead of F90 internal read. (bmy, 9/29/03)
!  (2 ) Now replace 2 and 4 digit year strings for all models (bmy, 10/23/03)
!  (3 ) Renamed LINUX to LINUX_PGI (bmy, 12/2/03)
!  (4 ) Now do not replace "ss" with seconds, as the smallest GEOS-CHEM
!        timestep is in minutes. (bmy, 7/20/04)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                         :: YYYY, YY, MM, DD, HH, II, SS
      CHARACTER(LEN=2)                :: MM_STR, DD_STR
      CHARACTER(LEN=2)                :: HH_STR, II_STR, SS_STR
      CHARACTER(LEN=2)                :: YY_STR
      CHARACTER(LEN=4)                :: YYYY_STR

      !=================================================================
      ! EXPAND_DATE begins here!
      !=================================================================

      ! Extract today's date into year, month, and day sections
      CALL YMD_EXTRACT( YYYYMMDD, YYYY, MM, DD )

      ! Extract today's time into HH, MM, and SS sections
      ! (rename minutes to II so as not to overwrite MM)
      CALL YMD_EXTRACT( HHMMSS, HH, II, SS )

      ! 2-digit year number (e.g. "97" instead of "1997")
      YY = YYYY - 1900
      IF ( YY >= 100 ) YY = YY - 100

      ! For other platforms, use an F90 internal write (bmy, 9/29/03)
      WRITE( YYYY_STR, '(i4.4)' ) YYYY
      WRITE( YY_STR,   '(i2.2)' ) YY
      WRITE( MM_STR,   '(i2.2)' ) MM
      WRITE( DD_STR,   '(i2.2)' ) DD
      WRITE( HH_STR,   '(i2.2)' ) HH
      WRITE( II_STR,   '(i2.2)' ) II

      ! Replace YYYY, MM, DD, HH tokens w/ actual values 
      CALL STRREPL( FILENAME, 'YYYY', YYYY_STR )
      CALL STRREPL( FILENAME, 'YY',   YY_STR   )
      CALL STRREPL( FILENAME, 'MM',   MM_STR   )
      CALL STRREPL( FILENAME, 'DD',   DD_STR   )
      CALL STRREPL( FILENAME, 'hh',   HH_STR   )
      CALL STRREPL( FILENAME, 'mm',   II_STR   )

      ! Return to calling program
      END SUBROUTINE expandDate
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MakeCharArrayFromCharList
!
! !DESCRIPTION: Subroutine makeCharArrayFromCharList takes a separated 
!  word list, and places each word into a separate element of a character 
!  array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MakeCharArrayFromCharList( list, separator, array )
!
! !INPUT PARAMETERS: 
!

      ! String with comma-separated words
      CHARACTER(LEN=*), INTENT(IN)  :: list

      ! String for separator text
      CHARACTER(LEN=1), INTENT(IN)  :: separator
!
! !OUTPUT PARAMETERS:
!
      ! Array of substrings
      CHARACTER(LEN=*), INTENT(OUT) :: array(:)
!
! !REVISION HISTORY: 
!  17 Jan 2006 - R. Yantosca - Initial version
!  08 Nov 2006 - R. Yantosca - Now set the output "array" argument to '' 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: P, N, ind(2048)
      CHARACTER(LEN=1) :: C

      ! Initialize
      N     = 1
      ind   = 0
      array = ''

      ! Find the positions of all the commas in LIST
      DO P = 1, LEN( list )

         ! Look at each character individually
         C = list(P:P)

         ! If a comma... 
         IF ( C == separator ) THEN 

            ! Increment comma
            N      = N + 1
            ind(N) = P 
         ENDIF
      ENDDO

      ! Add the position of the end of the string into IND
      ind(N+1) = LEN( list )

      ! Save text between the commas into ARRAY
      DO P = 1, N
         IF ( P == N ) THEN 
            array(P) = list( ind(P)+1:ind(P+1)   )
         ELSE
            array(P) = list( ind(P)+1:ind(P+1)-1 )
         ENDIF
      ENDDO

      END SUBROUTINE MakeCharArrayFromCharList
!EOC     
      END MODULE charpakModule

