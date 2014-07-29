#!/usr/bin/perl -w

#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: NrtUtils.pm
#
# !DESCRIPTION: Contains functions that are used to generate NRT simulations
#  with GEOS-Chem, esp. in support of aircraft campaigns.
#\\
#\\
# !INTERFACE:
#
package NrtUtils;
#
# !USES:
#
require 5.003;      # need this version of Perl or newer
use English;        # Use English language
use strict;         # Force explicit variable declarations (like IMPLICIT NONE)
#
# !PUBLIC MEMBER FUNCTIONS:
#  &checkDir      : Ensures that a directory exists
#  &fmtStr        : Pads a date string w/ leading zeroes if necessary
#  &makeInputGeos : Creates a new input.geos file for each day of simulation
#  &makeRunScript : Creates a GEOS-Chem run script for NRT simulations
#  &replaceDate   : Replaces YYYY, MM, DD tokens in a string w/ date values
#
# !CALLING SEQUENCE:
#  use NrtUtils qw( function-name1, function-name2, ... );
#
# !REVISION HISTORY:
#  20 Jun 2013 - R. Yantosca - Initial version, moved other routines here
#  30 Jul 2013 - R. Yantosca - Added function &checkDir
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
  @ISA       = qw( Exporter       );                   # export method
  @EXPORT_OK = qw( &fmtStr        
                   &makeInputGeos          
                   &makeRunScript
                   &replaceDate   
                   &checkDir      );                   # export on request
}
#EOP
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: fmtStr
#
# !DESCRIPTION: Routine fmtStr returns a date/time string in either YYYYMMDD 
#  or HHMMSS format.  The string is padded with leading zeroes if necessary.
#\\
#\\
# !INTERFACE:
#
sub fmtStr($) {
#
# !INPUT PARAMETERS:
#
  my ( $num ) = @_;   # Value to pad (if necessary)
#
# !RETURN VALUE:
#
  my $str     = "";   # Modified string
#
# !CALLING SEQUENCE:
#  $dateStr = &fmtStr( 20040101 );
#  $dateStr = &fmtStr( 0        );
#
# !REVISION HISTORY:
#  23 May 2013 - R. Yantosca - Initial version, based on NRT-ARCTAS
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  my $tmp = int( $num );

  # Pad w/ proper # of leading zeroes (if necessary)
  if    ( $tmp == 0                      ) { $str = "000000";    }
  elsif ( $tmp >= 100   && $tmp < 1000   ) { $str = "000$tmp";   }
  elsif ( $tmp >= 1000  && $tmp < 10000  ) { $str = "00$tmp";    }
  elsif ( $tmp >= 10000 && $tmp < 100000 ) { $str = "0$tmp";     }
  else                                     { $str = "$tmp";      }

  # Return to calling program
  return( $str );
}
#EOP
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: makeInputGeos
#
# !DESCRIPTION: Routine makeInputGeos constructs the "input.geos" file for 
#  GEOS-Chem.  It reads a pre-defined template file and then just replaces 
#  tokens # with the start and end date and time.
#\\
#\\
# !INTERFACE:
#
sub makeInputGeos($$$$$$) {
#
# !INPUT PARAMETERS:
#
  # $date1    : Starting date for GEOS-Chem model run (e.g. 20040101) 
  # $time1    : Starting time for GEOS-Chem model run (e.g. 000000  ) 
  # $date2    : Ending   date for GEOS-Chem model run (e.g. 20040102)
  # $time2    : Ending   time for GEOS-Chem model run (e.g. 000000  ) 
  # $template : Path for input.geos "template" file
  # $fileName : Path for input.geos file (w/ dates replaced)
  my ( $date1, $time1, $date2, $time2, $inFile, $outFile ) = @_;
#
# !CALLING SEQUENCE:
# &makeInputGeos( 20130101,             000000, 
#                 20130102,             000000, 
#                "input.geos.template", "input.geos" );
#
# !REVISION HISTORY:
#  23 May 2013 - R. Yantosca - Initial version, adapted from NRT-ARCTAS
#  31 Jul 2013 - R. Yantosca - Change permission of input.geos to chmod 777
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  my @lines  = "";
  my $line   = "";
  my $dStr1  = &fmtStr( $date1 );
  my $dStr2  = &fmtStr( $date2 );
  my $tStr1  = &fmtStr( $time1 );
  my $tStr2  = &fmtStr( $time2 );

  #------------------------------  
  # Read template file
  #------------------------------ 

  # Read template "input.geos" file into an array
  open( I, "$inFile" ) or croak( "Cannot open $inFile!\n" );
  @lines = <I>;
  close( I );

  #------------------------------  
  # Create "input.geos" file
  #------------------------------ 

  # Open file
  open( O, ">$outFile") or die "Can't open $outFile\n";

  # Loop thru each line
  foreach $line ( @lines ) {
    
    # Remove newline character
    chomp( $line );

    # Replace start & end dates
    $line =~ s/{DATE1}/$dStr1/g;
    $line =~ s/{TIME1}/$tStr1/g;
    $line =~ s/{DATE2}/$dStr2/g; 
    $line =~ s/{TIME2}/$tStr2/g;

    # Write to output file
    print O "$line\n";
  }

  # Close output file
  close( O );

  # Make sure the output file has chmod 777
  chmod( 0777, $outFile );

  # Exit
  return(0);
}

#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: makeRunScript
#
# !DESCRIPTION: Routine makeRunScript creates the script which will run the 
#  GEOS-Chem model.
#\\
#\\
# !INTERFACE:
#
sub makeRunScript($$$$$$$$$$$$) {
#
# !INPUT PARAMETERS:
#
  # $date_0    : Today's YYYYMMDD date
  # $date_1    : Tomorrow's YYYYMMDD date
  # $jobDir    : Directory where run scripts will be sent
  # $globDir   : GEOS-Chem run directory, for global simulation
  # $globExe   : GEOS-Chem executable, for global simulation
  # $nestDir   : GEOS-Chem run directory, for nested simulation
  # $nestExe   : GEOS-Chem executable, for nested simulation
  # $tpbcDir   : TPCORE boundary conditions directory
  # $tpbcFile  : TPCORE boundary conditions filename (e.g. BC.YYYYMMDD)
  # $logFile   : Log file for today's NRT simulation
  # $web       : Web update command
  # $runScript : Name of script file to be created (in $jobDir)
  #
  my( $date_0,   $date_1,  $jobDir,  $globDir,  
      $globExe,  $nestDir, $nestExe, $tpBcDir,  
      $tpBcFile, $logFile, $web,     $runScript ) = @_;
#
# !REMARKS:
#  Creates a tcsh script that is used to run the GEOS-Chem global and
#  nested simulations.  This needs to be submitted to a queue.
#
# !REVISION HISTORY:
#  21 Jun 2013 - R. Yantosca - Initial version, based on ARCTAS
#  01 Jul 2013 - R. Yantosca - Corrected a few typos
#  31 Jul 2013 - R. Yantosca - Make everything 
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
#-----begin "HERE document" ---------------------------------------------------
  # Define the script text
  my $txt = <<EOF;
#!/bin/tcsh -f

#-------------------------------------------------------------
# sleepNrt.$date_0: Run script for SEAC4RS NRT for $date_0
# Generated by sleepNrt
#-------------------------------------------------------------

# Define local variables
set DATE     = "$date_0"
set BCDIR    = "$tpBcDir"
set BCFILE   = "$tpBcFile"
set BCFILE_0 = `echo \${BCFILE} | sed 's/YYYYMMDD/'$date_0'/g'`
set JOBDIR   = "$jobDir"
set GLOBDIR  = "$globDir"
set GLOBEXE  = "$globExe"
set NESTDIR  = "$nestDir"
set NESTEXE  = "$nestExe"
set LOGFILE  = "$logFile"
set WEB      = "$web"

#------------------------------------------------------------
# Do some error checking before starting simulations
#------------------------------------------------------------

# Make everything chmod 777 by default
umask 000

# Check if the directory for TPCORE BC's exists
if ( !( -d \${BCDIR} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find directory:"    >> \${LOGFILE}
   echo "\${BCDIR}"                            >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}  
   exit(1)
endif

# Check if the global run dir exists
if ( !( -d \${GLOBDIR} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find directory:"    >> \${LOGFILE}
   echo "\$GLOBDIR"                            >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(2)
endif

# Check if the nested run dir exists
if ( !( -d \${NESTDIR} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find directory:"    >> \${LOGFILE} 
   echo "\$NESTDIR"                            >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(3)
endif

# Check if the global run dir exists
if ( !( -f \${GLOBDIR}\/\${GLOBEXE} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find executable:"   >> \${LOGFILE}
   echo "\${GLOBDIR}\/\${GLOBEXE}"              >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(4)
endif

# Check if the nested run dir exists
if ( !( -f \${NESTDIR}\/\${NESTEXE} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find executable:"   >> \${LOGFILE}
   echo "\${NESTDIR}\/\${NESTEXE}"              >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(5)
endif

#------------------------------------------------------------
# Run the global simulation first (to save TPCORE BC's)
#------------------------------------------------------------

# Switch to global simulation directory
cd \${GLOBDIR}
echo "Switched to directory: \$PWD" >> \${LOGFILE}

# Echo info
echo "===> $date_0: Global run began at' `date`" >> \${LOGFILE}

# Run the global simulation
\${GLOBEXE} >> \${LOGFILE}
if ( \$status != 0 ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Run exited abnormally!"       >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(6)
endif

# Echo info
echo "===> $date_0: Global run ended at' `date`" >> \${LOGFILE}
echo ""                                       >> \${LOGFILE}

#------------------------------------------------------------
# Run the nested simulation (using TPCORE BC's)
#------------------------------------------------------------

# Switch to global simulation directory
cd \${NESTDIR}
echo "Switched to directory: \$PWD" >> \${LOGFILE}
echo ""

# Make sure the TPCORE boundary file exists
if ( !( -f \${BCDIR}\/\${BCFILE_0} ) ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Could not find BC file:"      >> \${LOGFILE}
   echo "\${BCDIR}\/\${BCFILE_0}"               >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(7)
endif

# Echo info
echo "===> $date_0: Nested run began at' `date`" >> \${LOGFILE}
echo ""                                       >> \${LOGFILE}

# Run the nested simulation
\${NESTEXE} >> \${LOGFILE}
if ( \$status != 0 ) then
   echo '-----------------------------------' >> \${LOGFILE}
   echo "ERROR: Run exited abnormally!"       >> \${LOGFILE}
   echo '-----------------------------------' >> \${LOGFILE}
   exit(8)
endif

# Echo info
echo "===> $date_0: Nested run ended at' `date`" >> \${LOGFILE}
echo ""     

#------------------------------------------------------------
# Copy output files to the web directory
#------------------------------------------------------------

# Echo info
echo "===> $date_0: Copying output files to web directory"

# Copy bpch file
cp -f \${NESTDIR}/bpch/ctm.bpch.\${DATE} \${WEB}

# Copy timeseries files
cp -f \${NESTDIR}/timeseries/ts*\${DATE}* \${WEB}

# Copy planeflight files
cp -f \${NESTDIR}/plane/plane.log.\${DATE} \${WEB}

#------------------------------------------------------------
# Cleanup and quit
#------------------------------------------------------------
unset BCDIR
unset BCFILE
unset BCFILE_0
unset DATE
unset JOBDIR
unset GLOBDIR
unset GLOBEXE
unset NESTDIR
unset NESTEXE
unset LOGFILE
unset WEB

exit(0)
EOF
#-----end "HERE document" ---------------------------------------------------

  # Write text to script 
  open( J, ">$runScript" ) or die "Can't open $runScript!\n"; 
  print J "$txt\n";
  close( J ); 

  # Make run script executable (chmod 777)
  chmod( 0777, $runScript );

  # Return normally
  return(0);
}
#EOP
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: replaceDate
#
# !DESCRIPTION: Routine replaceDate replaces date tokens (YYYY, MM, DD) in
#  a string with the actual year, month, and date values.
#\\
#\\
# !INTERFACE:
#
sub replaceDate($$) {
#
# !INPUT PARAMETERS:
#
  my ( $str, $date ) = @_;  # $str: String w/ tokens; 
                            # $date: YYYYMMDD date
#
# !RETURN VALUE:
#
  my $newStr = "";          # Updated string 
#
# !CALLING SEQUENCE:
#  $newStr = &replaceDate( "file.YYYYMMDD", 20130101 );
#
# !REVISION HISTORY:
#  23 May 2013 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  my $yyyy = substr( $date, 0, 4 );    # Extract year  from $date
  my $mm   = substr( $date, 4, 2 );    # Extract month from $date
  my $dd   = substr( $date, 6, 2 );    # Extract day   from $date

  # Replace tokens
  $newStr =  $str;          
  $newStr =~ s/YYYY/$yyyy/g;           # Replace year 
  $newStr =~ s/MM/$mm/g;               # Replace month
  $newStr =~ s/DD/$dd/g;               # Replace day

  # Return modified string
  return( $newStr );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: checkDir
#
# !DESCRIPTION: Checks to make sure a directory exists.  If not, then it
#  will exit with an error code.
#\\
#\\
# !INTERFACE:
#
sub checkDir($) { 
#
# !INPUT PARAMETERS:
#
  # Directory to be tested
  my ( $dir ) =  @_;
#
# !CALLING SEQUENCE:
#  &checkDir( $dir );
#
# !REVISION HISTORY:
#  30 Jul 2013 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

  # Halt execution if directory is not found
  if ( !( -d $dir ) ) {    
    print "Directory $dir does not exist!  Exiting.";
    exit(999)
  }
  
  # Otherwise return w/ error status
  return( $? );
}
#EOC

END {}
