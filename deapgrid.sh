#!/bin/bash

########################################################
## SET THESE VARIABLES
########################################################
#PROCESSOFFSET=4000	# XXX use this to offset PROCESS (job) number, OR SET TO ZERO OTHERWISE - IT MUST BE SET
PROCESSOFFSET=0
let HISTOGRAMS_PER_JOB=10

SOURCEFILEDIR="/pnfs/annie/persistent/users/moflaher/deapfit/DEAPfit"
SOURCEFILE="DEAP_fitting.tar.gz"
INPUTDIR="/pnfs/annie/persistent/users/moflaher/deapfit_inputs"
OUTPUTDIR="/pnfs/annie/persistent/users/moflaher/deapfit_results"
OUTPUTDIRFALLBACK="/pnfs/annie/persistent/users/moflaher/deapfit_results" # fallback if ifdh mkdir fails to work*
OUTFILEBASE="R1214_to_R1231_PMTStability_mk2";

# NOTE: ifdh mkdir doesn't accept the -p argument. XXX actually it doesn't seem to work at all!!
# ENSURE OUTPUT DIRECTORY EXISTS
#OUTPUTDIR="."

#INPUTFILE="LEDRun1178S0LEDs5And6And10And28_PulseWindowOnly_PMTStability_Run0.root"
#INPUTFILE="R1214S0_1000V_PMTStability_Run0.root"
#INPUTFILE="R1214S1_1100V_PMTStability_Run32676.root"
#INPUTFILE="R1215S0_1200V_PMTStability_Run0.root"
#INPUTFILE="R1216S0_1300V_PMTStability_Run0.root"
#INPUTFILE="R1217S0_1400V_PMTStability_Run0.root"
#INPUTFILE="R1218S0_1500V_PMTStability_Run0.root"
#INPUTFILE="R1219S0_1600V_PMTStability_Run32664.root"
#INPUTFILE="R1220S0_1700V_PMTStability_Run0.root"
#INPUTFILE="R1227S0_1050V_PMTStability_Run0.root"
#INPUTFILE="R1228S0_1150V_PMTStability_Run0.root"
#INPUTFILE="R1229S0_1250V_PMTStability_Run0.root"
#INPUTFILE="R1230S0_1350V_PMTStability_Run32668.root"
#INPUTFILE="R1231S0_1450V_PMTStability_Run0.root"
INPUTFILES=("R1227S0_1050V_PMTStability_Run0.root" "R1214S1_1100V_PMTStability_Run32676.root" "R1228S0_1150V_PMTStability_Run0.root" "R1215S0_1200V_PMTStability_Run0.root" "R1229S0_1250V_PMTStability_Run0.root" "R1216S0_1300V_PMTStability_Run0.root" "R1230S0_1350V_PMTStability_Run32668.root" "R1217S0_1400V_PMTStability_Run0.root" "R1231S0_1450V_PMTStability_Run0.root" "R1218S0_1500V_PMTStability_Run0.root" "R1219S0_1600V_PMTStability_Run32664.root" "R1220S0_1700V_PMTStability_Run0.root")
echo "input file list is: ${INPUTFILES[@]}"

# these are suitable with one input file, but not many
# generate output filename
#OUTFILEBASE=${INPUTFILE%.*} # same as input filename, just strip off extention and add '_PROCESSNUM'
# derive output directory
#RUN=${INPUTFILE%%_*}
#REM=${INPUTFILE##${RUN}_}
#VOLTS=${REM%%_*}
#OUTFNAME=${RUN}_${VOLTS}
#OUTPUTDIRBASE=${OUTPUTDIR}
#OUTPUTDIR=${OUTPUTDIR}/${OUTFNAME}
#echo "Input file is ${INPUTDIR}/${INPUTFILE}, Output dir is ${OUTPUTDIR}, Output basename is ${OUTFNAME}, Run ${RUN}, PMTs at ${VOLTS}"

## make the output directory if it doesn't exist
## =============================================
echo "checking if output directory already exists"
if [ -d ${OUTPUTDIR} ]; then
    echo "output directory already exists"
else
    echo "making the output directory ${OUTPUTDIR}"
    ifdh mkdir ${OUTPUTDIR}
fi
if [ ! -d ${OUTPUTDIR} ]; then
    echo "ifdh mkdir is useless, switching output directory to ${OUTPUTDIRFALLBACK}" # BUT WHY?
    OUTPUTDIR=${OUTPUTDIRFALLBACK}
fi

########################################################
## SET THESE VARIABLES
########################################################


echo "setting up software base"
#export CODE_BASE=/grid/fermiapp/products
export CODE_BASE=/cvmfs/fermilab.opensciencegrid.org/products
source ${CODE_BASE}/common/etc/setup
export PRODUCTS=${PRODUCTS}:${CODE_BASE}/larsoft  # only artdaq has new enough ROOT!

echo "setting up products"
setup ifdhc   # for copying geometry & flux files
export IFDH_CP_MAXRETRIES=2  # default 8 tries is silly
setup fife_utils

setup root v6_18_04b -q debug:e17  # has to be >6.12 or there's a memory leak in TF1Convolution
source ${ROOTSYS}/bin/thisroot.sh

# calculate which histograms this job will be fitting
let PROCESSNUM=${PROCESS}  # conversion to number for arithmetic
let THECOUNTER=${PROCESSNUM}+${PROCESSOFFSET}
let STARTINDEX=$((${THECOUNTER}*${HISTOGRAMS_PER_JOB}))
echo "PROCESSNUM=${PROCESSNUM}, PROCESSOFFSET=${PROCESSOFFSET}, THECOUNTER=PROCESSNUM+PROCESSOFFSET=${THECOUNTER}"
echo "HISTOGRAMS_PER_JOB=${HISTOGRAMS_PER_JOB}"
echo "STARTINDEX=THECOUNTER*HISTOGRAMS_PER_JOB=${STARTINDEX}"
let ENDINDEX=${STARTINDEX}+${HISTOGRAMS_PER_JOB}

echo "this is job PROCESS=${PROCESS}, will process histogram ${ENDINDEX} histograms from offset ${STARTINDEX}"

# build the other filenames
OUTFILE=${OUTFILEBASE}_${PROCESSNUM}.root
OUTLOG=${OUTFILEBASE}_${PROCESSNUM}.log
echo "output files will be ${OUTPUTDIR}/${OUTFILE}, ${OUTPUTDIR}/${OUTLOG}"

# Skip job if output file already exists
echo "checking if output file already exists"
if [ -f ${OUTPUTDIR}/${OUTFILE} ]; then
    echo "output file already exists, skipping this job"
    exit 0
else
    echo "it doesn't"
fi
# TODO: add something to create a new temp outputdir based on job num and run the job anyway?

# copy the source files
echo "searching for source files in ${SOURCEFILEDIR}/${SOURCEFILE}"
echo "ifdh ls ${SOURCEFILEDIR}"
ifdh ls ${SOURCEFILEDIR}
ifdh ls ${SOURCEFILEDIR}/${SOURCEFILE} 1>/dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "copying source files"
  ifdh cp -D ${SOURCEFILEDIR}/${SOURCEFILE} .
else
  echo "source file not found in ${SOURCEFILEDIR}!"
fi

# copy the input files
#echo "copying the input file ${INPUTDIR}/${INPUTFILE}"
for AFILE in ${INPUTFILES[@]}; do
	echo "copying the input file ${INPUTDIR}/${AFILE}"
	ifdh cp -D ${INPUTDIR}/${AFILE} .
	if [ ! -f ${AFILE} ]; then echo "input file not found!!!"; exit 15; fi
done

echo "compiling application"
tar -xzf ${SOURCEFILE}
make

# run executable here, rename the output file
NOWS=`date "+%s"`
DATES=`date "+%Y-%m-%d %H:%M:%S"`
echo "checkpoint start @ ${DATES} s=${NOWS}"
echo " "

./deapfit --start ${STARTINDEX} --numhistos ${HISTOGRAMS_PER_JOB} ${INPUTFILES[@]} 2>&1 | tee ${OUTLOG}
RESULTSFILE="DEAPfitter_outfile.root"     # this is the file generated by the deapfit executable

echo " "
NOWF=`date "+%s"`
DATEF=`date "+%Y-%m-%d %H:%M:%S"`
let DS=${NOWF}-${NOWS}
echo "checkpoint finish @ ${DATEF} s=${NOWF}  ds=${DS}"
echo " "

echo "copying the output files to ${OUTPUTDIR}"
# copy back the output files
echo "copying ${OUTLOG} to ${OUTPUTDIR}"
ifdh cp -D ${OUTLOG} ${OUTPUTDIR}
if [ $? -ne 0 ]; then echo "something went wrong with the copy of ${OUTLOG}?!"; fi
# only attempt this if there was a results file, there may have been no good histos
# if we attempt to copy an output file that doesn't exist, the whole job hangs!
if [ -f ${RESULTSFILE} ]; then
	echo "renaming results file from ${RESULTSFILE} to ${OUTFILE}"
	mv ${RESULTSFILE} ${OUTFILE}
	echo "copying ${OUTFILE} to ${OUTPUTDIR}"
	ifdh cp -D ${OUTFILE} ${OUTPUTDIR}
	if [ $? -ne 0 ]; then echo "something went wrong with the copy of ${OUTFILE}?!"; fi
else
	echo "This job produced no results file!"
fi

# clean things up
#rm -rf ./*   <<< DON'T nuke everything, this will remove condor logs
#                 and when it fails to transfer those out and they're not present,
#                 that will cause the job to be put on HOLD
