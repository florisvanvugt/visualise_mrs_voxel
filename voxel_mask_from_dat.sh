#!/bin/bash
#This script creates a voxel mask image that can be used to 
#determine the tissue components within an MRS voxel.  

#This code is based on jn_voxel_mask_from_rda.sh.  This code retrieves
#the voxel coordinates from the .dat file, rather than the .rda file. 
#This is important because it removes the need to save the .rda files
#from the scans. 

# Jamie Near 02/08/2010
# Developed with a great deal of help of Mark Jenkinson
# Using some bits of code from Matt Taylor (get_voxel_from_rda.sh)

# Matt Taylor's original code was not designed to handle MRS voxels that 
# are tilted with respect to the scanner coordinate space.
# This code has been modified to handle MRS voxels that are tilted
# in any orientation

# This script now accepts structural images with any resolution (not 
# only 1x1x1mm).  Modification by Mark Jenkinson 

#NOTE:  THIS SCRIPT IS CURRENTLY FULL OF BUGS.  MANY VOXEL ORIENTATIONS WILL FAIL


# FVV Note Nov 2017: This script should produce struc/anat_mask.nii.gz which is a mask corresponding
# to the MRS voxel placement.


STRUCFILE=$1
MRS_DAT=$2
OUTFILE=$3


USAGE="Usage: jn_voxel_mask_from_dat.sh <STRUC_NII> <MRS_DAT> [<OUT_FILE>]"

## FVV Nov 2017; added a few tests to ensure the user gives the correct input
if [ -z $STRUCFILE ]; then
    echo $USAGE
    exit
fi

if [ -z $MRS_DAT ]; then
    echo $USAGE
    exit
fi



if [ ! -e $STRUCFILE ]; then
    echo "Structural file does not exist."
    exit
fi

if [ ! -e $MRS_DAT ]; then
    echo "MRS data file does not exist."
    exit
fi




#SCANDIR=${HOME}/documents/data/studies/${STUDYID}
#MRS_DAT=${SCANDIR}/${INDIV}/${SUBDIR}/*.dat
#STRUCFILE=${SCANDIR}/${INDIV}/images/${INDIV}/Anatomy*nii.gz
TEMPDIR=./struc/temp
MATFILE=${TEMPDIR}/mat

if [ -z $OUTFILE ]; then # if the user did not supply an output file
    OUTFILE=./struc/vox_anat_mask
    echo "Output to $OUTFILE"
    exit
fi


echo "Voxel information ${MRS_DAT}; anatomical ${STRUCFILE}; output to ${OUTFILE}"


mkdir -p ${TEMPDIR}

echo "Getting voxel information from .dat file"
POS1=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sPosition.dSag | cut -d '=' -f2 | awk 'NR<2')
POS2=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sPosition.dCor | cut -d '=' -f2 | awk 'NR<2')
POS3=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sPosition.dTra | cut -d '=' -f2 | awk 'NR<2')
VOX1=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.dReadoutFOV | cut -d '=' -f2 | awk 'NR<2')
VOX2=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.dPhaseFOV | cut -d '=' -f2 | awk 'NR<2')
VOX3=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.dThickness | cut -d '=' -f2 | awk 'NR<2')
ZED1=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sNormal.dSag | cut -d '=' -f2 | awk 'NR<2')
ZED2=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sNormal.dCor | cut -d '=' -f2 | awk 'NR<2')
ZED3=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.sNormal.dTra | cut -d '=' -f2 | awk 'NR<2')
ROT=$(strings ${MRS_DAT} | egrep ^sSpecPara.sVoI.dInPlaneRot | cut -d '=' -f2 | awk 'NR<2')



# Test whether ZED3 is 1 (ZED3 may be a float)
ZED3one=$(echo $ZED3'==1' | bc -l)

# FVV Nov 2017: Adjusted the line below which failed for me when $ZED3 is not a float
# I manually defined a variable that checked ZED3 against 1 using bc and then add that into the mix here.
# I think the test was supposed to check whether $ZED3 equals 1 (using float comparison) AND $ROT equals 0 or is not set.
if [ "$ZED3one" -eq 1 -a -z "$ROT" ];
then
    echo "There is *no* voxel rotation information"
    ROW1=-1.0
    ROW2=0.0
    ROW3=0.0
    COL1=0.0
    COL2=1.0
    COL3=0.0
    ZED1=0.0
    ZED2=0.0
    ZED3=1.0
    
else

    echo "There *is* voxel rotation information"

    # FVV Added these three lines to make sure that ROT is not just empty but set to 0 explicitly.
    if [ -z "$ROT" ]; then
	ROT=0
    fi
    
    # Change the sign of the ZED vector:
    ZED1=`echo "$ZED1 * -1.0" | bc -l`
    ZED2=`echo "$ZED2 * -1.0" | bc -l`
    ZED3=`echo "$ZED3 * -1.0" | bc -l`

    echo $ZED1
    echo $ZED2
    echo $ZED3
    
    # Now make a rotation matrix that rotates about the ZED vector by an angle ROT:
    R11=`echo "c(${ROT})+(${ZED1}^2) * (1-c(${ROT}))" | bc -l`
    R12=`echo "(${ZED1}*${ZED2}*(1-c(${ROT})))-(${ZED3}*s(${ROT}))" | bc -l`
    R13=`echo "(${ZED1}*${ZED3}*(1-c(${ROT})))+(${ZED2}*s(${ROT}))" | bc -l`
    R21=`echo "(${ZED2}*${ZED1}*(1-c(${ROT})))+(${ZED3}*s(${ROT}))" | bc -l`
    R22=`echo "c(${ROT})+(${ZED2}^2) * (1-c(${ROT}))" | bc -l`
    R23=`echo "(${ZED2}*${ZED3}*(1-c(${ROT})))-(${ZED1}*s(${ROT}))" | bc -l`
    R31=`echo "(${ZED3}*${ZED1}*(1-c(${ROT})))-(${ZED2}*s(${ROT}))" | bc -l`
    R32=`echo "(${ZED3}*${ZED2}*(1-c(${ROT})))+(${ZED1}*s(${ROT}))" | bc -l`
    R33=`echo "c(${ROT})+(${ZED3}^2) * (1-c(${ROT}))" | bc -l`

    
    echo "${R11} ${R12} ${R13} 0" > ${MATFILE}_Raxisangle
    echo "${R21} ${R22} ${R23} 0" >> ${MATFILE}_Raxisangle
    echo "${R31} ${R32} ${R33} 0" >> ${MATFILE}_Raxisangle
    echo "0 0 0 1" >> ${MATFILE}_Raxisangle

    # Now, find a column vector that is perpindicular to ZED, but whose first
    # term is zero:
    colStart1=0
    colStart3=`echo "1/sqrt(((-1.0 * ${ZED3} / ${ZED2})^2)+1)" | bc -l`
    echo $colStart3
    colStart2=`echo "-1.0 * ${colStart3} * ${ZED3} / ${ZED2}" | bc -l `
    echo $colStart2

    echo "${colStart1} 0 0 0" > ${MATFILE}_colStart
    echo "${colStart2} 0 0 0" >> ${MATFILE}_colStart
    echo "${colStart3} 0 0 0" >> ${MATFILE}_colStart
    echo "0 0 0 1" >> ${MATFILE}_colStart

    echo Calling convert_xfm
    # Now rotate the colStart vector using the Raxisangle matrix to find the new colVector:
    convert_xfm -omat ${MATFILE}_colVector -concat ${MATFILE}_Raxisangle ${MATFILE}_colStart

    # Now extract the values of COL1 COL2 and COl3:
    COL1=$(strings ${MATFILE}_colVector | cut -d ' ' -f1 | awk 'NR==1')
    COL2=$(strings ${MATFILE}_colVector | cut -d ' ' -f1 | awk 'NR==2')
    COL3=$(strings ${MATFILE}_colVector | cut -d ' ' -f1 | awk 'NR==3')

    # Now, to get RowVector, do the cross product of COL and ZED:
    ROW1=`echo "($COL2 * $ZED3) - ($COL3 * $ZED2)" |bc`
    ROW2=`echo "($COL3 * $ZED1) - ($COL1 * $ZED3)" |bc`
    ROW3=`echo "($COL1 * $ZED2) - ($COL2 * $ZED1)" |bc`

    # Change the sign of the third element of the zed vector (Not sure why, but this is required):
    ZED3=`echo "$ZED3 * -1.0" | bc -l`


fi


# Print the values that were found
#echo subject and voxel = ${INDIV} ${SUBDIR}
echo POS = $POS1 $POS2 $POS3
echo VOX = $VOX1 $VOX2 $VOX3
echo ROW = $ROW1 $ROW2 $ROW3
echo COL = $COL1 $COL2 $COL3
echo ZED = $ZED1 $ZED2 $ZED3

# Get the voxel dimensions of the structural image, as we may need these 
# in the future
VST1=$(fslinfo ${STRUCFILE} | egrep ^pixdim1 | cut -d ' ' -f9)
VST2=$(fslinfo ${STRUCFILE} | egrep ^pixdim2 | cut -d ' ' -f9)
VST3=$(fslinfo ${STRUCFILE} | egrep ^pixdim3 | cut -d ' ' -f9)

# Print the values that were found
echo ${VST1} ${VST2} ${VST3}

# Change the sign of the x- and y-coordinates on the **ASSUMPTION** that 
# the nifti version of the scanner coordinates and the spectroscopy version
# of the scanner coordinates are related by x -> -x, y -> -y, z -> z

POS1=`echo "$POS1 * -1.0" | bc -l`
POS2=`echo "$POS2 * -1.0" | bc -l`
ROW1=`echo "$ROW1 * -1.0" | bc -l`
ROW2=`echo "$ROW2 * -1.0" | bc -l`
COL1=`echo "$COL1 * -1.0" | bc -l`
COL2=`echo "$COL2 * -1.0" | bc -l`

# Get the transformation from the structural image coordinate space (st)
# into the scanner coordinate space (sc)
fslval ${STRUCFILE} sto_xyz:1 > ${MATFILE}_st2sc
fslval ${STRUCFILE} sto_xyz:2 >> ${MATFILE}_st2sc
fslval ${STRUCFILE} sto_xyz:3 >> ${MATFILE}_st2sc
fslval ${STRUCFILE} sto_xyz:4 >> ${MATFILE}_st2sc

# Convert the _st2sc matrix so that the input coords (st) are in (FLIRT) mm, not voxels
# deal with potential handedness issues for NEUROLOGICAL vs RADIOLOGICAL
if [ `fslorient ${STRUCFILE}` = NEUROLOGICAL ] ; then
    echo NEUROLOGICAL COORDINATES
    xsize=`fslval ${STRUCFILE} dim1`;
    xsize=`echo "($xsize - 1) * $VST1" | bc -l`;
    xmul=`echo "$VST1 * -1" | bc -l`;
    echo "$xmul 0 0 $xsize" > ${MATFILE}_nvox2fmm
else
    echo NOT NEUROLOGICAL COORDINATES
    echo "$VST1 0 0 0" > ${MATFILE}_nvox2fmm
fi
echo "0 $VST2 0 0" >> ${MATFILE}_nvox2fmm
echo "0 0 $VST3 0" >> ${MATFILE}_nvox2fmm
echo "0 0 0 1" >> ${MATFILE}_nvox2fmm
convert_xfm -omat ${MATFILE}_fmm2nvox -inverse ${MATFILE}_nvox2fmm
convert_xfm -omat ${MATFILE}_st2sc -concat ${MATFILE}_st2sc ${MATFILE}_fmm2nvox

# We will need the inverse of this matrix:
convert_xfm -omat ${MATFILE}_sc2st -inverse ${MATFILE}_st2sc

# Now we are ready to make the initial mask.  Start by puting it in the 
# bottom corner of the image.  Only the size will be correct.
# convert VOX1 from mm to structural voxel units
VOX1vox=`printf %.0f $(echo $VOX1 / $VST1 | bc -l)`;
VOX2vox=`printf %.0f $(echo $VOX2 / $VST2 | bc -l)`;
VOX3vox=`printf %.0f $(echo $VOX3 / $VST3 | bc -l)`;
# FVV Added `printf %.0f` to the call above, because when you give non-rounded numbers to the scdript below (e.g. 19.9998) it will floor them, causing too small voxel locations
#echo fslmaths ${STRUCFILE} -mul 0 -add 1 -roi 0 ${VOX1vox} 0 ${VOX2vox} 0 ${VOX3vox} 0 1 ${TEMPDIR}_vox_start
fslmaths ${STRUCFILE} -mul 0 -add 1 -roi 0 ${VOX1vox} 0 ${VOX2vox} 0 ${VOX3vox} 0 1 ${TEMPDIR}_vox_start


# Okay, now we need to make a rotation matrix to rotate from the spectroscopy
# voxel coordinate space (sp), into the scanner coordinate space (sc)
echo $(echo ${ROW1} ${COL1} ${ZED1}) 0 > ${MATFILE}_sp2sc_R
echo $(echo ${ROW2} ${COL2} ${ZED2}) 0 >> ${MATFILE}_sp2sc_R
echo $(echo ${ROW3} ${COL3} ${ZED3}) 0 >> ${MATFILE}_sp2sc_R
echo 0 0 0 1 >> ${MATFILE}_sp2sc_R

# Now we need to make two translation matrices to translate from the spectroscopy
# voxel coordinate space, into the scanner coordinate space:
# The first translates by the voxel centre position in scanner space
echo 0 0 0 $(echo ${POS1})> ${MATFILE}_sc_Tc
echo 0 0 0 $(echo ${POS2})>> ${MATFILE}_sc_Tc
echo 0 0 0 $(echo ${POS3})>> ${MATFILE}_sc_Tc
echo 0 0 0 0 >> ${MATFILE}_sc_Tc

# And the second translates back by half of the voxel dimensions.
echo 1 0 0 $(echo -.5*${VOX1} | bc)> ${MATFILE}_sc_Tv
echo 0 1 0 $(echo -.5*${VOX2} | bc)>> ${MATFILE}_sc_Tv
echo 0 0 1 $(echo -.5*${VOX3} | bc)>> ${MATFILE}_sc_Tv
echo 0 0 0 1 >> ${MATFILE}_sc_Tv

# Next, we must concatenate sc2st with sc_Tc to make a new matrix, 
# which we will call:  st_Tc

convert_xfm -omat ${MATFILE}_st_Tc -concat ${MATFILE}_sc2st ${MATFILE}_sc_Tc

# Now we have all of the matrices that we need.  All we have to do now is 
# put them together by concatenating (as well as a bit of tricky manipulation
# as you will soon see).

# We need to make a matrix that contains the 4x4 identity matrix combined with 
# the last column of the st_Tc matrix.  We will call the resulting matrix 
# Id_st_Tc.  This is where the tricky manipulation happens:
col4=`cat ${MATFILE}_st_Tc | awk '{ print $4 }'`
t1=`echo $col4 | cut -d' ' -f1`
t2=`echo $col4 | cut -d' ' -f2`
t3=`echo $col4 | cut -d' ' -f3`

echo 1 0 0 $t1 > ${MATFILE}_Id_st_Tc
echo 0 1 0 $t2 >> ${MATFILE}_Id_st_Tc
echo 0 0 1 $t3 >> ${MATFILE}_Id_st_Tc
echo 0 0 0 1 >> ${MATFILE}_Id_st_Tc


# Now we just need to do some contatenations:
convert_xfm -omat ${MATFILE}_sp2st_R -concat ${MATFILE}_sc2st ${MATFILE}_sp2sc_R
convert_xfm -omat ${MATFILE}_sp2st_Tv_R -concat ${MATFILE}_sp2st_R ${MATFILE}_sc_Tv
convert_xfm -omat ${MATFILE}_final -concat ${MATFILE}_Id_st_Tc ${MATFILE}_sp2st_Tv_R

## FVV: no longer using flirt for the final transformation, instead using ANTS, which seemed slightly more accurate.
# If you don't have ANTS you can uncomment the following line and comment out anything that follows.
#flirt -in ${TEMPDIR}_vox_start -ref ${STRUCFILE} -out ${OUTFILE} -applyxfm -init ${MATFILE}_final -interp trilinear # nearestneighbour  # trilinear

#fslview ${STRUCFILE} ${OUTFILE}

#clean up
#rm -r ${TEMPDIR}


# Convert the fsl-style transformation matrix into something ants can use
cp ${MATFILE}_final ${MATFILE}_final.txt

AFFINETRANSF=${TEMPDIR}_vox_start2final.txt
c3d_affine_tool ${MATFILE}_final.txt -ref ${STRUCFILE} -src ${TEMPDIR}_vox_start -fsl2ras -oitk $AFFINETRANSF


antsApplyTransforms --default-value 0 --dimensionality 3 --interpolation Linear --input ${TEMPDIR}_vox_start.nii.gz --output ${TEMPDIR}_vox_final.nii.gz --reference-image ${STRUCFILE} --transform [ $AFFINETRANSF, 0 ]
cp ${TEMPDIR}_vox_final.nii.gz ${OUTFILE}

