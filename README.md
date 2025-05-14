[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/GEOS_FP/blob/master/LICENSE.txt)

# GEOS_FP

Development of the code and scripts used to process the "raw" GEOS-FP met data for input into GEOS-Chem

## Instructions for processing GEOS-FP met fields

**by Bob Yantosca (@yantosca), 13 Feb 2014**

**Updated obsolete info and fixed broken links: 14 May 2025**

### 1. The GEOS-FP “raw” data files

The GEOS-FP “raw” data files (native resolution 0.25&deg; x 0.3125&deg; x 72 vertical levels) are grouped into several “collections”, or file types.  The collections of “raw” data files that you will need to download are described in more detail on the [GEOS-FP wiki page](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-FP#GEOS-FP_data_file_collections).

Within each GEOS-FP data collection, several different “raw” data files are archived for each day.  Collections containing surface data contain 24 files per day, while collections containing 3D data contain 8 individual files per day.  The data file names contain the collection name plus a timestamp, such as:

```console
GEOS.fp.asm.tavg3_3d_rad_Nv.YYYYMMDD_hhmm.V01.nc4
GEOS.fp.asm.tavg3_3d_mst_Ne.YYYYMMDD_hhmm.V01.nc4
etc.
```

File names contain the following elements:

1. `GEOS.fp.asm`, which indicates that the data is the GEOS-FP assimilated met product  
2. The name of each [data collection](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-FP#GEOS-FP_data_file_collections) (i.e. tavg3_3d_rad_Nv, etc.)  
3. A date and time stamp
4. A version stamp (i.e. `V01`)
5. the `.nc4` suffix which indicates that each file is stored in the netCDF-4 data format

The GEOS-FP “raw” data files may be downloaded from the server ftp.nccs.nasa.gov.  In the following sections we will describe the scripts that you can use to download these data files to your computer system.

## 2. Regridding the “raw” GEOS-FP data files into files that GEOS-Chem will read

The software described in this manual will perform the following functions:

* Download the GEOS-FP “raw” data files to your computer system  
* Extract and horizontally regrid fields that you select to files that GEOS-Chem can read.  Supported horizontal grids are:
    * 4 x 5 global
    * 2 x 2.5 global
    * 0.25 x 0.3125 nested grids (China, Europe, N. America, and SE Asia)
  
For each day of “raw” GEOS-FP data, the GEOS-FP data processing software will extract and regrid a subset of fields to the horizontal grids mentioned above. It will create a new set of output files (in netCDF-3 format) that GEOS-Chem can read directly.  These files use the following naming convention:

```console
GEOSFP.20110101.CN.*.nc
GEOSFP.YYYYMMDD.A1.*.nc
GEOSFP.YYYYMMDD.A3cld.*.nc
GEOSFP.YYYYMMDD.A3dyn.*.nc
GEOSFP.YYYYMMDD.A3mstC.*.nc
GEOSFP.YYYYMMDD.A3mstE.*.nc
GEOSFP.YYYYMMDD.I3.*.nc
```
The date string `YYYYMMDD` represents the day of data.  The wildcard character * will be replaced by a string representing the horizontal resolution of the data (e.g. `4x5`, `2x25`, `025x03125.CH`, `025x03125.EU`, `025x03125.NA`, `025x03125.SE`).

In the sections below, you will learn how to download and run the GEOS-FP data processing software.

## 3. Requirements

The GEOS-FP “raw” met data files are stored in the netCDF4 file format.  You will need to install the netCDF4 library, which also requires HDF5 and zlib libraries to be built.  Someone at your group may have done this already; if so, then you don’t need to re-install netCDF4.

## 4. Downloading the GEOS-FP data processing software for the first time

The GEOS-FP data processing software is distributed via the Git version control system.  If you are downloading this software for the very first time, then type at the Unix prompt:

```console
$ git clone https://github.com/geoschem/GEOS_FP
```

This will cause a directory called `GEOS_FP`, containing several subdirectories, to be downloaded to your disk space.  We will look at the contents of the GEOS_FP directory in more detail in the following sections.  

The `GEOS_FP` directory is a complete clone of the repository stored on GitHub. You can use the gitk browser to examine revision history of the GEOS_FP directory.  You can also use the git gui to commit additional changes.  See the GEOS-Chem wiki page [*Using Git with GEOS-Chem*](http://wiki.seas.harvard.edu/geos-chem/index.php/Using_Git_with_GEOS-Chem) for more information about these commands.  (This page describes how to use Git with the GEOS-Chem code, but you would use the same Git commands when working with the GEOS-FP data processing code.)

## 5. Obtaining updated versions of the GEOS-FP data processing software

The `git clone` command only has to be done once.  To apply further updates to the software into your already-downloaded `GEOS_FP` directory structure, you can use the `git pull` command:

```console
$ git pull https://github.com/geoschem/GEOS_FP main
```

Here, `main` is the branch containing the updates that you want to apply to your copy of the code.  Most of the time you will be pulling updates into the `main` branch.  

## 6. Subdirectories of the GEOS_FP root directory

The GEOS-FP data processing software is stored in the following subdirectories of GEOS_FP: 

1. `Code/`  : Contains Fortran code  
2. `adjust/`  : Contains NCL code (not used so much anymore)  
3. `bin/`  : Contains the executable and the input files which specify data directories & options for the data reprocessing  
4. `doc/` : Manual pages will be built here when you type make doc  
5. `jobs/` : Contains job scripts created by the doGeosFpMulti driver script  
6. `lib/` : Fortran library files (*.a) will be built here during compilation  
7. `logs/` : log files from the GEOS-FP data processing jobs will be sent here.  
8. `mod/` : Fortran module files (*.mod) will be built here during compilation  
9. `perl/` : Contains scripts (written in Perl) to download and process GEOS-FP data

Of these directories, `Code/`, `bin/`, `logs/`, and `perl/` are the most important.  You can more or less neglect the others.

## 7. Setting the proper environment variables

Before you build the executables, you will need to define some environment variables that specify the location of your netCDF4 library installation.

### 7a. With csh or tcsh

If you use the `csh` or `tcsh` Unix shells, then include these lines in your `.cshrc` file:

```csh
# Root path of your GEOS-Chem-Libraries installation  
# (this may be different on your system, so ask where it is!)  
setenv ROOT_LIB_DIR  /opt/GEOS-Chem-Libraries/ifort

# Links to the netCDF bin, include, lib directories   
setenv GC_BIN      $ROOT_LIB_DIR/nc4/bin  
setenv GC_INCLUDE  $ROOT_LIB_DIR/nc4/include    
setenv GC_LIB      $ROOT_LIB_DIR/nc4/lib  
setenv BIN_NETCDF  $GC_BIN  
setenv INC_NETCDF  $GC_INCLUDE  
setenv LIB_NETCDF  $GC_LIB
```

Then type `source ~/.cshrc` at the Unix prompt to apply the changes.  

### 7b. With bash

If, on the other hand, you use the bash Unix shell, then add these lines:

```bash
# Root path of your GEOS-Chem-Libraries installation  
# (this may be different on your system, so ask where it is!)  
export ROOT_LIB_DIR=/opt/GEOS-Chem-Libraries/ifort  

# Links to the netCDF bin, include, lib directories   
export GC_BIN=$ROOT_LIB_DIR/nc4/bin  
export GC_INCLUDE=$ROOT_LIB_DIR/nc4/include    
export GC_LIB=$ROOT_LIB_DIR/nc4/lib  
export BIN_NETCDF=$GC_BIN  
export INC_NETCDF=$GC_INCLUDE  
export LIB_NETCDF=$GC_LIB
```

Then type `source ~/.bashrc` at the Unix prompt to apply the changes.  

## 8. Compiling the Fortran code

To compile the Fortran source code into executables, type:

```console
$ cd GEOS_FP/Code  
$ make all
```

This will create the following executables:

```console
GeosFpDriver0.x  
GeosFpDriver1.x  
GeosFpDriver2.x  
GeosFpFixA3Cld.x
```

Copies of these executable files will be created in both the Code/ and bin/ subdirectories of GEOS_FP.  Of these, you can ignore the `GeosFpFixA3Cld.x`, as that was created for a specific reprocessing of certain erroneous fields.  This was done in late 2013 by Sajeev Philip.

## 9. Editing the GeosFpDriver.input file

In the bin/ subdirectory, there is an input file named `GeosFpDriver.input` that lets you specify the following settings:

* Filename structure of each GEOS_FP “raw” data file  
* The fields in each GEOS_FP “raw” data file that you want to extract and/or regrid  
* The data directory to which the GEOS_FP “raw” data files have been downloaded  
* The types of output you want to create, including  
* 4x5 global output  
* 2x25 global output  
* 0.25&deg; x 0.3125&deg; nested grids for  
  * CH: China  
  * EU: Europe  
  * NA: North America  
  * SE: Southeast Asia  
* Your “scratch” directory (i.e. temporary directory used for creating the files)  
* The active data directory into which you want the files to be stored in

We shall examine each of the relevant sections of `GeosFpDriver.input` below.  Lines beginning with the `#` character are comments and will be ignored.

### 9a
```console
==> Turn on debug print output?  
F
```
If this value is set to `T`, then extra debugging output will be printed to the log files.  Setting this to `F` will not print the debugging output.  Most of the time you can leave this set to `F`.

### 9b
```console
==> const_2d_asm_Nx  
GEOS.fp.asm.const_2d_asm_Nx.00000000_0000.V01.nc4  
FRLAKE,FRLAND,FRLANDICE,FROCEAN,PHIS  
F
```

This section of `GeosFpDriver.input` specifies the fields that you want to extract and regrid from the GEOS-FP time-invariant data file ([data collection](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-FP#GEOS-FP_data_file_collections) const_2d_asm_nx).  Here we specify the filename structure (line 2) and the fields that you want to extract/regrid (line 3).  Line 4 lets you toggle the creation of the `GEOSFP.20110101.CN.*.nc` files on or off.  Normally you will leave this set to `F`, since you only have to create these files once (and this has already been done).

### 9c
```console
==> tavg1_2d_flx_Nx  
GEOS.fp.asm.tavg1_2d_flx_Nx.YYYYMMDD_hhmm.V01.nc4  
EFLUX,EVAP,FRSEAICE,HFLUX,PBLH,PRECANV,PRECCON,PRECLSC,PRECSNO,PRECTOT,USTAR,Z0M

==> tavg1_2d_lnd_Nx  
GEOS.fp.asm.tavg1_2d_lnd_Nx.YYYYMMDD_hhmm.V01.nc4  
FRSNO,GRN,GWETROOT,GWETTOP,LAI,PARDF,PARDR,SNODP,SNOMAS

==> tavg1_2d_rad_Nx  
GEOS.fp.asm.tavg1_2d_rad_Nx.YYYYMMDD_hhmm.V01.nc4  
ALBEDO,CLDTOT,LWGNT,LWTUP,SWGDN

==> tavg1_2d_slv_Nx  
GEOS.fp.asm.tavg1_2d_slv_Nx.YYYYMMDD_hhmm.V01.nc4  
QV2M,SLP,TROPPT,TS,T2M,U10M,V10M,T10M,Q850,TO3

==> tavg3_3d_asm_Nv  
GEOS.fp.asm.tavg3_3d_asm_Nv.YYYYMMDD_hhmm.V01.nc4  
OMEGA,U,V

==> tavg3_3d_cld_Nv  
GEOS.fp.asm.tavg3_3d_cld_Nv.YYYYMMDD_hhmm.V01.nc4  
A3cld: QCCU,QIAN,QILS,QI,QLAN,QLLS,QL,TAUCLI,TAUCLW,OPTDEPTH  
A3dyn: DTRAIN,RH

==> tavg3_3d_mst_Ne  
GEOS.fp.asm.tavg3_3d_mst_Ne.YYYYMMDD_hhmm.V01.nc4  
CMFMC,PFICU,PFILSAN,PFLCU,PFLLSAN

==> tavg3_3d_mst_Nv  
GEOS.fp.asm.tavg3_3d_mst_Nv.YYYYMMDD_hhmm.V01.nc4  
DQRCU,DQRLSAN,REEVAPCN,REEVAPLSAN

==> tavg3_3d_rad_Nv  
GEOS.fp.asm.tavg3_3d_rad_Nv.YYYYMMDD_hhmm.V01.nc4  
CLOUD

==> inst3_3d_asm_Nv  
GEOS.fp.asm.inst3_3d_asm_Nv.YYYYMMDD_hhmm.V01.nc4  
PS,PV,QV,T
```
This section lets you specify the file name structure for each of the GEOS_FP “raw” data file types.  The software will replace the date and time tokens `YYYYMMDD_hhmm` with the data and time of each data file.  Furthermore, for each file type, you can specify a list of individual met fields that you want to be extracted and regridded into the GEOS-Chem netCDF output files.

Normally you will not have to edit this section, as the list of GEOS-FP fields that we will archive for GEOS-Chem will not change.

**SPECIAL NOTE:** Fields from the `tavg3_3d_cld_Nv` collection are used to create two different GEOS-Chem output data files, the `A3cld` and `A3dyn` files.  The first line following the file name structure lets you specify the fields that will be included into `A3cld`, and the line after that lets you specify fields that will be included into `A3dyn`.

### 9d
```console
==> Local Raw Data Path  
/as/scratch/bmy/GEOS_FP/
```

The above section lets you specify the data directory into which the GEOS-FP “raw” data files are stored.  The above directory is what we use at Harvard; consult with your sysadmin as to where this directory is located on your computer system.

NOTE: Make sure the directory path you specify ends with a `/` character.

### 9e
```console
==> Nested CH output  
F  
GEOSFP.YYYYMMDD.%%%%%%.025x03125.CH.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/as/data/geos-rw/GEOS_0.25x0.3125_CH.d/GEOS_FP/YYYY/MM/  
  801  421  1025 581   # China grid
```

This section lets you extract GEOS-FP met data for the CH (China) nested grid at 0.25&deg; x 0.3125&deg; resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the CH nested grid, or `F` to skip creating files for the CH nested grid.   
* **Line 3:** Specify the file name structure of the output files for the `CH` nested grid.    
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, `I3`).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path without using a temporary directory, then make sure that both Line 4 and Line 5 are identical.  
* **Line 6:** Specify the lon & lat indices of the lower left corner `(I0,J0)` and upper right corner `(I1,J1)` that define the `CH` nested grid region.

### 9f
```console
==> Nested EU output  
F  
GEOSFP.YYYYMMDD.%%%%%%.025x03125.EU.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_0.25x0.3125/GEOS_0.25x0.3125_EU.d/GEOS_FP/YYYY/MM/  
  161 400 385 601  # NOTE: Need correct indices
```
This section lets you extract the GEOS-FP met data for the EU (European) nested grid at 0.25o x 0.3125o resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the EU nested grid, or `F` to skip creating files for the EU nested grid.   
* **Line 3:** Specify the file name structure of the output files for the EU nested grid.   
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, `I3`).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path without using a temporary directory, then make sure that both Line 4 and Line 5 are identical.  
* **Line 6:** Specify the lon & lat indices of the lower left corner `(I0,J0)` and upper right corner `(I1,J1)` that define the `EU` nested grid region.

### 9g
```console
==> Nested NA output  
T  
GEOSFP.YYYYMMDD.%%%%%%.025x03125.NA.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_0.25x0.3125/GEOS_0.25x0.3125_NA.d/GEOS_FP/YYYY/MM/  
  161 400 385 601
```
This section lets you extract the GEOS-FP met data for the NA (North American) nested grid at 0.25&deg; x 0.3125&deg; resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the NA nested grid, or `F` to skip creating files for the `NA` nested grid.   
* **Line 3:** Specify the file name structure of the output files for the NA nested grid.  
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, `I3`).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path, then make sure that both Line 4 and Line 5 are identical.  
* **Line 6:** Specify the lon & lat indices indices of the lower left corner `(I0,J0)` and upper right corner `(I1,J1)` that define the NA nested grid region.

### 9h
```console
==> Nested SE output  
F  
GEOSFP.YYYYMMDD.%%%%%%.025x03125.SE.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_0.25x0.3125/GEOS_0.25x0.3125_SE.d/GEOS_FP/YYYY/MM/  
  817 321 993 481
```

This section lets you extract the GEOS-FP met data for the SE (SE Asian) nested grid at 0.25&deg; x 0.3125&deg; resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the SE nested grid, or `F` to skip creating files for the `SE` global grid.  
* **Line 3:** Specify the file name structure of the output files for the SE nested grid.  
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, I3).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path, then make sure that both Line 4 and Line 5 are identical.  
* **Line 6:** Specify the lon & lat indices indices of the lower left corner `(I0,J0)` and upper right corner `(I1,J1)` that define the `SE` nested grid region.

### 9i
```console
==> 2 x 2.5 output  
T  
GEOSFP.YYYYMMDD.%%%%%%.2x25.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_2x2.5.d/GEOS_FP/YYYY/MM/
```
This section lets you extract the GEOS-FP met data for the global grid at 2o x 2.5o resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the 2&deg; x 2.5&deg; global grid, or `F` to skip creating files for the 2&deg; x 2.5&deg; global grid.  
* **Line 3:** Specify the file name structure of the output files for the 2o x 2.5o global grid nested grid.  
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, `I3`).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path, then make sure that both Line 4 and Line 5 are identical.

### 9j 
==> 4 x 5 output
```console
T  
GEOSFP.YYYYMMDD.%%%%%%.4x5.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_4x5/GEOS_FP/YYYY/MM/
```

This section lets you extract the GEOS-FP met data for the global grid at 4&deg;x 5&deg resolution:

* **Line 2:** Specify `T` to create GEOS-FP data files for the 4&deg; x 5&deg; global grid, or `F` to skip creating files for the 4&deg; x 5&deg; global grid.   
* **Line 3:** Specify the file name structure of the output files for the 4&deg; x 5&deg; global grid.  
  * The software will replace the `%%%%%%` with one of several strings that denote different files (e.g. `A1`, `A3cld`, `A3dyn`, `A3mstC`, `A3mstE`, `I3`).  
* **Line 4:** Specify a temporary directory where the output files will be created.  
* **Line 5:** Specify the data directory path into which the output files will be moved once they have been created.  For example, this will be the directory path on your FTP server.  
  * If you wish to write the output data files directly into your FTP server path, then make sure that both Line 4 and Line 5 are identical.

### 9k
```console
==> Native-resolution wind speed output  
F  
GEOSFP.YYYYMMDD.SPEED.025x03125.nc  
/as/tmp/all/tmp/bmy/GEOS_FP/  
/mnt/gcgrid-rw/GEOS_0.25x0.3125.d/GEOS_FP/YYYY/MM/
```

This section lets you extract the wind speed fields at the native 0.25&deg; x 0.3125&deg; resolution, which are needed for David Ridley’s new dust mobilization algorithm.  NOTE: This feature is not fully implemented.

### 9l
```console
==> Mapping Weight Files  
weights_025x03125_to_2x25.txt  
weights_025x03125_to_4x5.txt
```

This section lets you specify the files that contain the mapping weights for the `MAP_A2A` algorithm.  The mapping weights determine the fraction of each 0.25&deg; x 0.3125&deg; grid box that fits inside each 2&deg; x 2.5&deg; or 4&deg; x 5&deg; grid box.  It is more computationally efficient to have computed these mapping weights once and then to read them in each time the software is called.  

By default, the mapping weights files are stored in the `bin/` subdirectory.  Therefore, you should not have to modify this section.

### 9m
```console
==> Template Files  
GeosFpTemplateFile.n
```
This section lets you specify the template file, which contains the mask for creating the land-water indices field (LWI).  By default, this file is stored in the `bin/` subdirectory.  Therefore, you should not have to modify this section.

## 10. Scripts that control the regridding process

In the perl/subdirectory, there are several scripts (mostly written in the Perl language) that control the entire data download and regridding process.  These are:

Scripts for downloading data

* `getGeosFp`: Downloads 1 day of GEOS-FP met data from the NASA FTP site  
* `checkGeosFp`: Checks to see if the GEOS-FP data files were downloaded

Scripts for file management

* `cleanJobs`: Removes job scripts in the jobs/ subdirectory  
* `cleanLogs`: Removes log files in the logs/ subdirectory  
* `moveGeosFp:` Moves files from temp dir to data dir (called by doGeosFp)  
* `delGeosFp`: Removes GeosFp “raw” met data files (for manual use)  
* `purgeGeosFp`: Removes old GEOS-FP “raw” data (called by doGeosFp)

Scripts for regridding data 

* `Dates.pm`: Perl module containing common subroutines  
* `doGeosFpMulti`: Main driver script; extracts/regrids 1 day of GEOS-FP data  
* `doGeosFp.input`: Input file with settings for doGeosFpMulti  
* `runMet`: Called by doGeosFp

Scripts for automating the data download and regrididng process

* `schedGeosFp`: Calls purgeGeosFp, getGeosFp, doGeosFp  
* `schedGeosFpInAdvance`: Schedules regridding jobs in days advance

Scripts you can ignore for now

* `doGeosFpFixA3Cld`  
* `makeGeosFpWindSpeed`  
* `make_native_sfc_wind.ncl`

Of these, you will probably use `getGeosFp`, `doGeosFpMulti`, `schedGeosFp`, and `schedGeosFpInAdvance` most frequently.  We shall look at these in more detail in each the following sections.

## 11. Editing the doGeosFp.input file

The `doGeosFp.input` file lets you specify several directory data paths that are used by the various scripts in the `perl/` subdirectory.

### 11a
```console
==> Raw Met Data Directory  
/as/scratch/bmy/GEOS_FP
```
In this section, you specify the directory path where the “raw” GEOS-FP data files have been downloaded into (or are stored) on your system.  This should be the same directory that you specified in the `==> Local Raw Data Path` line in the `bin/GeosFpDriver.input` file.

NOTE: It is OK to omit the trailing `/` character in the directory path that you specify here.

### 11b
```console
==> Code Directory  
../bin

==> Job Directory  
../jobs

==> Log Directory  
../logs

==> Temporary Directory  
../jobs

==> Program Executable  
GeosFpDriver{THREAD}.x

==> Defaults for Executable  
../bin/GeosFpDriver.input

==> Submit Statement  s
{JOB} &
```

This section tells the Perl scripts where to find the various directory paths for the Fortran Code, job scripts, log files.  It also indicates the naming convention used by the Fortran executable files, as well as where the executables can find the `GeosFpDriver.input` file.

Under normal circumstances, you will not have to modify this section of doGeosFp.input.

### 11c
```console
==> Sleep Time [s]  
300
```

The scripts will keep testing if all of the “raw” GEOS-FP data files have been downloaded to disk.  If not, the scripts will go to sleep for a specified period of time.  You can set that period of time (in seconds) here.

### 11d
```console
==> Emails for Notification  
user\@email.com
```

If you choose, the scripts can send an email to one or more persons to denote that it has completed extracting and/or regridding GEOS-FP met field data for a given date.  Specify a list of comma-separated email addresses here.

NOTE: Put a slash before the `@` sign in all email addresses.  This is necessary, because the `@`symbol is used in Perl to declare array variables.  The `@` tells Perl to interpret the `@` sign as a literal string instead of indicationg the start of an array.

### 11e
```console
#-----------------------------------------------------  
# File sizes for processed GEOS-FP met data:  
# Used by purgeGeosFp to remove old "raw" files          
# before getting new files  
#-----------------------------------------------------

==> CH Nested-Grid Met Fields  
/as/data/geos/GEOS_0.25x0.3125_NA/GEOS_FP/YYYY/MM  
 159985884  GEOSFP.YYYYMMDD.A1.025x03125.CH.nc  
 584242324  GEOSFP.YYYYMMDD.A3cld.025x03125.CH.nc  
 417317376  GEOSFP.YYYYMMDD.A3dyn.025x03125.CH.nc  
 333854836  GEOSFP.YYYYMMDD.A3mstC.025x03125.CH.nc  
 423113044  GEOSFP.YYYYMMDD.A3mstE.025x03125.CH.nc  
 251550928  GEOSFP.YYYYMMDD.I3.025x03125.CH.nc

==> EU Nested-Grid Met Fields  
  91857936  GEOSFP.YYYYMMDD.A1.025x03125.EU.nc  
 328290588  GEOSFP.YYYYMMDD.A3cld.025x03125.EU.nc  
 234494100  GEOSFP.YYYYMMDD.A3dyn.025x03125.EU.nc  
 187596040  GEOSFP.YYYYMMDD.A3mstC.025x03125.EU.nc  
 237751068  GEOSFP.YYYYMMDD.A3mstE.025x03125.EU.nc  
 141349284  GEOSFP.YYYYMMDD.I3.025x03125.EU.nc

==> NA Nested-Grid Met Fields  
/as/data/geos/GEOS_0.25x0.3125_NA/GEOS_FP/YYYY/MM  
 205087116  GEOSFP.YYYYMMDD.A1.025x03125.NA.nc  
 733023288  GEOSFP.YYYYMMDD.A3cld.025x03125.NA.nc  
 523589040  GEOSFP.YYYYMMDD.A3dyn.025x03125.NA.nc  
 418872100  GEOSFP.YYYYMMDD.A3mstC.025x03125.NA.nc  
 530861208  GEOSFP.YYYYMMDD.A3mstE.025x03125.NA.nc  
 315609504  GEOSFP.YYYYMMDD.I3.025x03125.NA.nc

==> SE Nested-Grid Met Fields  
 128594824  GEOSFP.YYYYMMDD.A1.025x03125.SE.nc  
 459604948  GEOSFP.YYYYMMDD.A3cld.025x03125.SE.nc  
 328290124  GEOSFP.YYYYMMDD.A3dyn.025x03125.SE.nc  
 262632896  GEOSFP.YYYYMMDD.A3mstC.025x03125.SE.nc  
 332849812  GEOSFP.YYYYMMDD.A3mstE.025x03125.SE.nc  
 197887516  GEOSFP.YYYYMMDD.I3.025x03125.SE.nc

==> 2x25 Global Met Fields  
/as/data/geos/GEOS_2x2.5/GEOS_FP/YYYY/MM  
  59141184  GEOSFP.YYYYMMDD.A1.2x25.nc  
 211346220  GEOSFP.YYYYMMDD.A3cld.2x25.nc  
 150962340  GEOSFP.YYYYMMDD.A3dyn.2x25.nc  
 120770584  GEOSFP.YYYYMMDD.A3mstC.2x25.nc  
 153059148  GEOSFP.YYYYMMDD.A3mstE.2x25.nc  
  90998100  GEOSFP.YYYYMMDD.I3.2x25.nc

==> 4x5 Global Met Fields  
/as/data/geos/GEOS_4x5/GEOS_FP/YYYY/MM  
  14959212  GEOSFP.YYYYMMDD.A1.4x5.nc  
  53420372  GEOSFP.YYYYMMDD.A3cld.4x5.nc  
  38158028  GEOSFP.YYYYMMDD.A3dyn.4x5.nc  
  30527044  GEOSFP.YYYYMMDD.A3mstC.4x5.nc  
  38688120  GEOSFP.YYYYMMDD.A3mstE.4x5.nc  
  23001984  GEOSFP.YYYYMMDD.I3.4x5.nc
```
If you so choose, the data processing software can delete the previous day’s GEOS-FP “raw” met field data files before starting to extract and/or regrid the current day’s “raw” data.  The script purgeGeosFp will check the output data files (matching them against the file names and sizes listed here) in order to determine if the output are of the proper size.  You should not have to change any of the file sizes above.

**NOTE:** You should only use the purgeGeosFp script if you do not have sufficient disk space to store the GEOS-FP “raw” data files.  You may want to store the “raw” data on your server for a period of time so that you can reprocess days if the need arises.

## 12. Downloading GEOS-FP “raw” data with the getGeosFp script

You can use the getGeosFp script to download 1 day of GEOS-FP “raw” met data at a time.  The data will be sent to the Raw Met Data Directory that you specified in doGeosFp.input (which is the same as the `Local Raw Met Path` in `bin/GeosFpDriver.input`).

To start the data download process, type at the Unix prompt:

```console
$ cd GEOS_FP/perl  
$ getGeosFp YYYYMMDD
```

**NOTE:** We recommend that you use the `schedGeosFp` or `schedGeosFpInAdvance` scripts, which will call `getGeosFp` for you.  This will let you start the data download process at a time when the file transfer speeds are greater (i.e. overnight).  We will discuss these scripts in the sections below.

## 13. Extracting and regridding GEOS-FP data with the doGeosFpMulti script

You can use the `getGeosFp` script to download 1 day of GEOS-FP “raw” met data at a time.  The data will be sent to the Raw Met Data Directory that you specified in `doGeosFp.input` (which is the same as the `Local Raw Met Path` in `bin/GeosFpDriver.input`).

To start the data download process, type at the Unix prompt:

```console
$ cd GEOS_FP/perl  
$ doGeosFpMulti YYYYMMDD
```

**NOTE:** We recommend that you use the `schedGeosFp` or `schedGeosFpInAdvance` scripts, which will call `doGeosFp` for you.  This will let you start the data downloading and regridding at a time when the file transfer speeds are greater (i.e. overnight).  We will discuss these scripts in the sections below.

## 14. Automating the downloading and regridding process with schedGeosFp

In order to make the GEOS-FP data extraction and regridding process more convenient for you, we have created a couple of scripts that you can use to (1) download several days day of “raw” GEOS-FP data from the NASA FTP site to your disk and (2) to start the extraction and regridding process.  

The first script is called schedGeosFp,which takes the following arguments

```console
$ schedGeosFp YYYYMMDD nDays when
```

where

* `YYYYMMDD` is the starting date of the GEOS-FP data that you want to download and regrid,

* `nDays` is the number of days of GEOS-FP data (including `YYYYMMDD`) that you want to download and regrid

* when lets you specify a time when you want the data downloading and processing to begin.  Some allowable values for when are:  
  * `now`   
  * `“2am tomorrow"`
  * `"8am Sunday"`
  * `"8am 2014/07/01"`  
  * NOTE: Using the `“”` characters will tell `schedGeosFp` to treat a space-separated string (e.g. `“2am tomorrow”`) as a single entry.

For example, to download and regrid the GEOS-FP data (starting immediately) for the dates 2014/01/01, 2014/01/02, 2014/01/03, 2014/01/04, and 2014/01/05, type at the Unix prompt:

```console
$ cd GEOS_FP/perl  
$ schedGeosFp 20140101 5 now
```

This will produce the following output:
```console
at now <<EOF  
./getGeosFp 20140101  
./doGeosFpMulti 20140101  
./getGeosFp 20140102  
./doGeosFpMulti 20140102  
./getGeosFp 20140103  
./doGeosFpMulti 20140103  
./getGeosFp 20140104  
./doGeosFpMulti 20140104  
./getGeosFp 20140105  
./doGeosFpMulti 20140105  
EOF  
job 72 at Tue Feb 11 12:18:00 2014
```
As you can see, the `schedGeosFp` script calls `getGeosFp` and `doGeosFpMulti` for each day of data that you want to process.  

The Unix at command schedules the job to run at the time you specified via the when argument.   You can get a list of all running jobs by typing:

`at -l`

or you can delete a running job with

`at -r JOBNUMBER`

where `JOBNUMBER` is the number of the job (shown via `at -l)`.  Depending on the flavor of your operating system (i.e. Linux, Ubuntu, Fedora, CentOS, MacOS), your at command may use slightly different options for these commands.  Check your manual pages to be sure.

NOTE: At present we have disabled the call to `purgeGeosFp` in order to prevent inadvertent deletion of “raw” GEOS-FP met data.  We developed `purgeGeosFp` for use with the SEAC4RS mission, as we did not have enough storage space at Harvard to store the “raw” GEOS-FP met data files indefinitely.

If you would like to restore the call to `purgeGeosFp`, then uncomment the **green line of code** in the schedGeosFp script, as shown below:

```perl
# Add commands to the string  
#-----------------------------------------------------------------------------  
# Prior to 2/11/14:  
# For now, disable the call to "purgeGeosFp" script. (bmy, 2/11/14)  
#$cmd .= "./purgeGeosFp $daten";    # <== Uncomment this line
#-----------------------------------------------------------------------------  
$cmd .= "./getGeosFp $daten";  
$cmd .= "./doGeosFpMulti $daten";
```

## 15. Using schedGeosFpInAdvance to schedule several data processing jobs

The GEOS-FP met data product is the active data product produced by NASA/GMAO.  The GEOS-FP “raw” data files for a given date are usually posted to the NASA FTP server by about 7AM ET the next calendar day.  For example, the entire set of GEOS-FP “raw” data files for 2014/02/14 should be available by 7AM ET on 2014/02/15.  Normally, you will want to download and process each day of GEOS_FP data as soon as it is ready.

One way to do this would be to invoke the schedGeosFp script several times.  For example, typing the following at the Unix prompt will schedule GEOS-FP data processing jobs for data dates 2014/02/15 thru 2014/02/20:

```console
$ cd GEOS_FP/perl  
$ schedGeosFp 20140215 1 “7am 2014/02/16”  
$ schedGeosFp 20140216 1 “7am 2014/02/17”  
$ schedGeosFp 20140217 1 “7am 2014/02/18”  
$ schedGeosFp 20140218 1 “7am 2014/02/19”  
$ schedGeosFp 20140219 1 “7am 2014/02/20”  
$ schedGeosFp 20140220 1 “7am 2014/02/21”
```

This is perfectly fine, albeit a little cumbersome.  Also having to manually specify the start dates is prone to error. In order to simplify this process, we have developed the schedGeosFpInAdvance script.  This script will call schedGeosFp repeatedly, as demonstrated above, so that you don’t have to do it yourself.

The schedGeosFpInAdvance script will assume that the GEOS-FP data processing jobs will start at 7am.  If you wish to change this, you can modify the this line of source code in schedGeosFpInAdvance:

```perl
      my $when = "7am"; # <===== You can change the start time here!
```

Also, if your system has weekly scheduled downtime, you can tell `schedGeosFpInAdvance` to start a little later by uncommenting and modifying the red lines of code:

```perl
#-------------------------------------------------------------------  
# NOTE: At Harvard, there is a weekly maintenance period   
# every Monday.  This code below allows you to start the data   
# processing job a little bit# later on Mondays.  &getDayOfWeek  
# returns 0=Sun 1=Mon 2=Tue 3=Wed 4=Thu 5=Fri 6=Sat.  
# (bmy, 2/11/14)  
## Start the job later on Mondays due to the maintenance window  
#if ( &getDayOfWeek( $tomorrow ) == 1 ) { $when = "10am"; }    # <== uncomment this line
#else                                   { $when = "7am";  }    # <== uncomment this line
#--------------------------------------------------------------------

Once you have set up `schedGeosFpInAdvance` with the proper start time, you can call it to start submitting GEOS-FP data processing jobs.  To download the and regrid the GEOS_FP met data for the same dates as in the above example (2014/02/15 thru 2014/02/20), we now type:

```console
$ cd GEOS_FP/perl  
$ schedGeosFpInAdvance 20140215 6
```

where in this instance 2014021 indicates the starting date, and 6 indicates the total number of days of GEOS_FP met data (including the starting date) to process.

## 16. Log file output from the data processing jobs

All log files from the GEOS-FP data processing software will get set to the Log Directory that you specified in the `perl/doGeosFp.input` file.  By default, this directory is `GEOS_FP/logs`.  Looking at the log files will help you determine if the data processing job finished properly.  

In order to save execution time time, the GEOS-FP data processing software will split the workload in parallel over 3 different Unix threads.  For each met field date that is being processed, 3 log files will be created. Typing the following at the Unix prompt:

```console
$ cd GEOS_FP/logs  
$ ls log*
```
will show output similar to this:

```console
log.doGeosFpMulti.20140215.22786-0
log.doGeosFpMulti.20140215.22786-1
log.doGeosFpMulti.20140215.22786-2
```
The number following the date (`22786`) is the job id # of the Unix process.  This number is unique, and ensures that log file output from previous data processing jobs for the same met field date will not overwrite each other.  

The last number in the filename (`0`, `1`, `2`) indicates the thread number.:

* The log file ending in `0` is created by the executable `bin/GeosFpDriver0.x`.    
  This executable creates the `CN` (if necessary), `A1`, and `I3` output files.

* The log file ending in `1` is created by the executable `bin/GeosFpDriver1.x`.    
  This executable creates the `A3cld` and `A3dyn` output files.

* The log file ending in `2` is created by the executable `bin/GeosFpDriver2.x`.    
  This executable creates the `A3mstC` and `A3mstE` output files.

After some time you may want to delete old log files in `GEOS_FP/logs`.  You can use the the following script to do this for you:

```console
$ cd GEOS_FP/perl  
$ cleanLogs
```

## 17. Job scripts

The `doGeosFpMulti` driver script will create several job scripts, which will then call the corresponding executable files to start the GEOS_FP data processing.  These scripts are created in the Job Directory that you specified in the `perl/doGeosFp.input` file.  By default, this directory is `GEOS_FP/jobs`.  Normally you won’t have to work with these job scripts directly.  

After some time you may want to delete old job fils in `GEOS_FP/jobs`.  You can use the following script to do this for you:

```console
$ cd GEOS_FP/perl  
$ cleanJobs
```

**18.  Removing “raw” data files once you no longer need them

As mentioned above, we recommend keeping the GEOS-FP “raw” met field files for a period of time so that you can reprocess them scratch if there is a problem in the regridded data files.  But when you are finally ready to delete these from your Raw Met Data Directory (i.e. the same path you specified in `doGeosFp.input`), you can use the `delGeosFp` script.  Type:

To delete “raw” met field files for a specific date, type:a specific date, 

```console
$ cd GEOS_FP/perl  
$ delGeosFp 20140215
```

Or, to delete an entire month, type:

```console
$ cd GEOS_FP/perl  
$ delGeosFp 201402
```

Or an entire year:

```console
$ cd GEOS_FP/perl  
$ delGeosFp 2014
```

The script will ask you three times if you REALLY want to delete the files.  You must answer `YES` (or `Y`) each time or the script will exit without doing anything.

**USE WITH CAUTION!!!!**
