#!/usr/bin/perl -w

# $Id: Dates.pm,v 1.2 2009/07/30 20:00:02 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Dates
#
# !DESCRIPTION: This Perl package contains handy algorithms for date and
#  time conversions, as well as getting the date and time from the system.
#\\
#\\
# !INTERFACE:
#
package Dates;
#
# !USES:
#
  require 5.003;   # need this version of Perl or newer
  use English;     # Use English language
  use Carp;        # Get detailed error messages
  use strict;      # Force explicit variable declarations (like IMPLICIT NONE)
#
#
# !PUBLIC MEMBER FUNCTIONS:
#  &julDay($$$)       
#  &mint($)       
#  &calDate($)      
#  &addDate($$)      
#  &getDayOfYear($) 
#  &getDayOfWeek($)
#  &getLocalTime()
#  &getUtcTime()  
#  &ymdCombine($$$)
#  &ymdExtract($)
#
# !CALLING SEQUENCE:
#  use Dates qw( function-name1, function-name2, ... );
#
# !REVISION HISTORY: 
#  bmy, 01 May 2002 - INITIAL VERSION
#  bmy, 10 Feb 2006 - added getDayOfYear routine
#  bmy, 22 Feb 2006 - bug fix in getDayOfYear routine
#  bmy, 07 Feb 2008 - added getDayOfWeek routine
#  bmy, 14 Feb 2008 - added getLocalTime, getUtcTime routines
#  bmy, 20 Feb 2008 - Bug fix: add 1 to month in getLocalTime, getUtcTime
#  bmy, 29 Jul 2009 - Added ProTeX headers
#  bmy, 29 Jul 2009 - Added routines &ymdCombine, &ymdExtract
#EOP
#------------------------------------------------------------------------------
#BOC

BEGIN {

  #=========================================================================
  # The BEGIN method lists the names to export to the calling routine
  #=========================================================================
  use Exporter ();
  use vars     qw( $VERSION @ISA @EXPORT_OK );

  $VERSION   = 1.00;                                   # version number
  @ISA       = qw( Exporter );                         # export method
  @EXPORT_OK = qw( &julDay       &mint         
                   &calDate      &addDate      
                   &getDayOfYear &getDayOfWeek 
                   &getLocalTime &getUtcTime 
                   &ymdCombine   &ymdExtract   );      # export on request
}

#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: julDay
#
# !DESCRIPTION: Subroutine julDay returns the astronomical Julian day.
#\\
#\\
# !INTERFACE:
#
sub julDay($$$) {
#
# !INPUT PARAMETERS:
#
  # $year  : Current year
  # $month : Current month
  # $day   : Current day (can be fractional, e.g. 17.25)
  my( $year, $month, $day ) = @_;
#
# !CALLING SEQUENCE:
#  $jd = julDay( $year, $month, $day );  
#
# !REMARKS:
#  (1) Algorithm taken from "Practical Astronomy With Your Calculator",
#      Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
#  (2) Requires function mint(x), defined in this package.
#  (3) Hardwired for Gregorian dates only.  Let's not worry about the
#       Julian Calendar.
#
# !REVISION HISTORY:
#  01 May 2002 - R. Yantosca - Initial version.
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#    
  my( $year1, $month1, $x1, $a  ) = ( 0,   0,   0.0, 0.0 );
  my( $b,     $c,      $d,  $jd ) = ( 0.0, 0.0, 0.0, 0.0 );
  
  #=========================================================================
  # julDay begins here!
  #=========================================================================

  # Compute YEAR1 and MONTH1
  if ( ( $month == 1 ) || ( $month == 2 ) ) {
    $year1  = $year  - 1;
    $month1 = $month + 12;
  } else {
    $year1  = $year;
    $month1 = $month;
  }
  
  # Compute the A term
  $x1 = $year / 100.0;
  $a  = &mint( $x1 );
    
  # Compute the "B" term according to Gregorian or Julian calendar
  $b = 2.0 - $a + mint( $a / 4.0 );
  
  # Compute the "C" term for BC dates (YEAR1 <= 0 ) 
  # or AD dates (YEAR1 > 0)
  if ( $year1 < 0 ) {
    $x1 = ( 365.25 * $year1 ) - 0.75;
    $c  = &mint( $x1 );
  } else {
    $x1 = 365.25 * $year1;
    $c  = &mint( $x1 );
  # Return variables

  }

  # Compute the D term
  $x1 = 30.6001 * ( $month1 + 1 );
  $d  = &mint( $x1 );
  
  # Compute the Julian day
  $jd = $b + $c + $d + $day + 1720994.5;
  
  # Return to calling program
  return $jd;
  
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: mint
#
# !DESCRIPTION: Subroutine mint computes the modified integer function (which
#  is required by subroutine julDay).
#\\
#\\
# !INTERFACE:
#
sub mint($) {
#
# !INPUT PARAMETERS:
#
  # $x: Value to be tested
  my( $x ) = @_;
#
# !CALLING SEQUENCE:
#  $value = mint( $x );
#
# !REMARKS:
#  The modified integer function is defined as follows
#
#            { -INT( ABS( X ) );  X <  0
#     MINT = { 
#            {  INT( ABS( X ) );  X >= 0
#
# !REVISION HISTORY:
#  01 May 2002 - R. Yantosca - Initial version
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#  
  my $value = 0.0;
  
  #=========================================================================
  # mint begins here!
  #=========================================================================
  if ( $x < 0 ) { 
    $value = -int( abs( $x ) );
  } else {
    $value = int( abs( $x ) );
  }
    
  # Return to calling program
  return $value;
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: calDate
#
# !DESCRIPTION:  Subroutine calDate takes an astronomical Julian day and
#  returns the corresponding year, month, and day.
#\\
#\\
# !INTERFACE:
#
sub calDate($) {
#
# !INPUT PARAMETERS:
#
  # $jdIn: Astronomical Julian date for which you want to return Y/M/D
  my( $jdIn ) = @_;
#
# !CALLING SEQUENCE:
#  ( $year, $month, $day ) = calDate( $jdIn );
#
# !REMARKS:
#  Algorithm taken from "Practical Astronomy With Your Calculator",
#  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
#
# !REVISION HISTORY:
#  01 May 2002 - R. Yantosca - Initial version
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#  
  my( $a,   $b,  $c,  $d,   ) = ( 0.0, 0.0, 0.0, 0.0 );
  my( $day, $e,  $f,  $fday ) = ( 0.0, 0.0, 0.0, 0.0 );
  my( $g,   $i,  $j,  $jd,  ) = ( 0.0, 0.0, 0.0, 0.0 );
  my( $m,   $y              ) = ( 0.0, 0.0,          );

  #=========================================================================
  # calDate begins here!
  #=========================================================================
  $jd = $jdIn + 0.5;
  $i  = int( $jd );
  $f  = $jd - int( $i );
  
  if ( $i > 2299160.0 ) {
    $a = int( ( $i - 1867216.25 ) / 36524.25 );
    $b = $i + 1 + $a - int( $a / 4 );
  } else {
    $b = $i
  }

  # Compute intermediate quantities
  $c = $b + 1524.0;
  $d = int( ( $c - 122.1 ) / 365.25 );
  $e = int( 365.25 * $d );
  $g = int( ( $c - $e ) / 30.6001 );
    
  # day is the day number
  $day = $c - $e + $f - int( 30.6001 * $g );
  
  # fday is the fractional day number
  $fday = $day - int( $day );
  
  # m is the month number
  if ( $g < 13.5 ) {
    $m = $g - 1;
  } else {
    $m = $g - 13;
  }
  
  # y is the year number
  if ( $m > 2.5 ) {
    $y = $d - 4716.0;
  } else {
    $y = $d - 4715.0;
  }
  
  # Return values: year, month, day (incl. fractional part)
  return ( $y, $m, $day );
  
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: addDate
#
# !DESCRIPTION: Subroutine addDate adds a specified number of days to a 
#  calendar date and returns the result.  This is useful for straddling the 
#  end of a month or the end of a year.
#\\
#\\
# !INTERFACE:
#
sub addDate($$) {
#
# !INPUT PARAMETERS:
#
  # $nymd0   : Starting date in YYYYMMDD format
  # $addDays : Number of days to add to $nymd0
  my( $nymd0, $addDays ) = @_;
#
# !CALLING SEQUENCE:
#  $nymd1 = addDate( $nymd0, $addDays );
#
# !REVISION HISTORY:
#  20 Aug 2001 - R. Yantosca - Initial version
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  my( $y0,    $m0,    $d0 ) = ( 0,   0,   0.0 );
  my( $y1,    $m1,    $d1 ) = ( 0,   0,   0.0 );
  my( $nymd1, $nhms1, $jd ) = ( 0.0, 0.0, 0.0 );
  
  #=========================================================================
  # addDate begins here!
  #========================================================================= 
  
  # Translate $nymd0 into $y0, $m0, $d0
  $y0 = int( $nymd0 / 10000 );
  $m0 = int( ( $nymd0 - int( $y0 * 10000 ) ) / 100 );
  $d0 = $nymd0 - int( $y0 * 10000 ) - int( $m0 * 100 );
  
  # Compute the astronomical julian day for the starting date
  $jd = &julDay( $y0, $m0, $d0 );
  
  # Add the offset to jd
  $jd = $jd + $addDays;
  
  # Convert the new Julian day back to a calendar date
  ( $y1, $m1, $d1 ) = &calDate( $jd );
    
  # Convert the new calendar date into YYYYMMDD format
  $nymd1 = int( $y1 * 10000 ) + int( $m1 * 100 ) + int( $d1 );

  # Return to calling program
  return $nymd1;
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getDayOfYear
#
# !DESCRIPTION: Subroutine getDayOfYear takes in a date in YYYYMMDD format
#  and returns the day of year (0-365) or (0-366 in leap years).
#\\
#\\
# !INTERFACE:
#
sub getDayOfYear($) {
#
# !INPUT PARAMETERS:
#
  # $nymd : Date in YYYYMMDD format
  my ( $nymd ) = @_;
#
# !CALLING SEQUENCE:
#  $doy = getDayOfYear( $nymd );
#
# !REVISION HISTORY:
#  10 Feb 2006 - R. Yantosca - Initial version
#  22 Feb 2006 - R. Yantosca - Bug fix: old version returned 1st day of month
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Convert YMD to year-month-date
  my $year     = int( $nymd / 10000 );
  my $month    = int( ( $nymd - int( $year * 10000 ) ) / 100 );
  my $day      = $nymd - int( $year * 10000 ) - int( $month * 100 );

  # Get Julian day for 1st day of the year
  my $jd0      = julDay( $year, 1, 1 );

  # Get Julian day for start of the month
  my $jd1      = julDay( $year, $month, $day );

  # Day of year DDD value for 1st & last days of month
  my $doy      = $jd1 - $jd0 + 1;

  # Return to the calling program
  return( $doy );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getDayOfWeek
#
# !DESCRIPTION: Subroutine getDayOfWeek returns the day of week (0=Sun, 
#  1=Mon, ..., 6=Sat) corresponding to a given YYYY/MM/DD date.
#\\
#\\
# !INTERFACE:
#
sub getDayOfWeek($) {
#
# !INPUT PARAMETERS:
#
  # $nymd : Date in YYYYMMDD format
  my ( $nymd ) = @_;
#
# !CALLING SEQUENCE:
#  $doy = getDayOfWeek( $nymd );
#
# !REMARKS:
#  Algorithm taken from "Practical Astronomy With Your Calculator",
#  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
#
# !REVISION HISTORY:
#  07 Feb 2008 - R. Yantosca - Initial version
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Convert YMD to year-month-date
  my $year     = int( $nymd / 10000 );
  my $month    = int( ( $nymd - int( $year * 10000 ) ) / 100 );
  my $day      = $nymd - int( $year * 10000 ) - int( $month * 100 );

  # Get day of week -- Algorithm from Peter Duffett-Smith
  my $jd       = julDay( $year, $month, $day );
  my $tmp      = ( $jd + 1.5 ) / 7.0;
  my $dow      = int( ( ( $tmp - int( $tmp ) ) * 7.0 ) + 0.5 );

  # Return to calling program
  return( $dow );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getLocalTime
#
# !DESCRIPTION: Subroutine getLocalTime returns the local time as a string 
#  in the format "YYYY/MM/DD hh:mm:ss".  Output is from the Perl localtime 
#  function.
#\\
#\\
# !INTERFACE:
#
sub getLocalTime() {
#
# !CALLING SEQUENCE:
#  $timeStr = getLocalTime();
#
# !REVISION HISTORY:
#  14 Feb 2008 - R. Yantosca - Initial version
#  20 Feb 2008 - R. Yantosca - Bug fix: Need to add 1 to the month 
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Call the Perl localtime function
  my @r = localtime(); 

  # Format the output string
  my $str = sprintf( "%04d", $r[5] + 1900 ) . "/" . 
            sprintf( "%02d", $r[4] + 1    ) . "/" .
            sprintf( "%02d", $r[3]        ) . " " .
            sprintf( "%02d", $r[2]        ) . ":" .
            sprintf( "%02d", $r[1]        ) . ":" .
            sprintf( "%02d", $r[0]        );     

  # Return
  return( $str );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getUtcTime
#
# !DESCRIPTION: Subroutine getUtcTime returns the UTC (aka GMT) time as a
#  string in the format "YYYY/MM/DD hh:mm:ss".  Output is from the Perl 
#  gmtime function.
#\\
#\\
# !INTERFACE:
#
sub getUtcTime() {
#
# !CALLING SEQUENCE:
#  $timeStr = getUtcTime();
#
# !REVISION HISTORY:
#  14 Feb 2008 - R. Yantosca - Initial version
#  20 Feb 2008 - R. Yantosca - Bug fix: Need to add 1 to the month 
#  29 Jul 2009 - R. Yantosca - Added ProTeX headers
#EOP
#------------------------------------------------------------------------------
#BOC

  # Call the Perl gmtime function
  my @r = gmtime(); 

  # Format the output string
  my $str = sprintf( "%04d", $r[5] + 1900 ) . "/" . 
            sprintf( "%02d", $r[4] + 1    ) . "/" .
            sprintf( "%02d", $r[3]        ) . " " .
            sprintf( "%02d", $r[2]        ) . ":" .
            sprintf( "%02d", $r[1]        ) . ":" .
            sprintf( "%02d", $r[0]        );     

  # Return
  return( $str );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: ymdCombine
#
# !DESCRIPTION: Subroutine ymdCombine combines separate year, month, and date
#  values into a single YYYYMMDD value.  (Can also be used to combine separate
#  values for hours, minutes, seconds into a single hhmmss time value.)
#\\
#\\
# !INTERFACE:
#
sub ymdCombine($$$) {
#
# !INPUT PARAMETERS:
#
  # $year  : Year  (YYYY) or hours   (hh)
  # $month : Month (MM)   or minutes (mm)
  # $day   : Day   (DD)   or seconds (ss)
  my ( $year, $month, $day ) = @_;
#
# !CALLING SEQUENCE:
#  $date = &ymdCombine( $year, $month, $day );
#
# !REVISION HISTORY:
#  29 Jul 2009 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
  my $date = ( $year * 10000 ) + ( $month * 100 ) + $day;
  return( $date );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: ymdExtract
#
# !DESCRIPTION: Subroutine ymdExtract splits a YYYYMMDD date into separate 
#  year, month, and day values.  Can also be used to split a hhmmss time 
#  into hours, minutes, and seconds 
#\\
#\\
# !INTERFACE:
#
sub ymdExtract($) {
#
# !INPUT PARAMETERS:
#
  # $date: Date in YYYYMMDD format (or time in hhmmss format)
  my ( $date ) = @_;
#
# !CALLING SEQUENCE:
#  ( $year, $month, $day ) = &ymdExtract( $date );
#
# !REVISION HISTORY:
#  29 Jul 2009 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
  my $y  = int( $date / 10000 );
  my $m  = int( ( $date - int( $y * 10000 ) ) / 100 );
  my $d  = $date - int( $y * 10000 ) - int( $m * 100 );
  return ( $y, $m, $d );
}
#EOC

END {}
