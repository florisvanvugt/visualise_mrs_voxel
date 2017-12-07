# MRS Voxel Visualisation

This script is based heavily on work by Near and Taylor. All I (Floris) did was some tweaking.

This little script takes an MRS `.dat` file (from the scanner) and a subject anatomical T1 image, and produces a mask file that can be overlaid on the subject anatomical T1 image that shows the coverage of the voxel.

## Requirements

* FSL
* ANTS (you could also do the transformation without ANTS but the latter seems more accurate)


## Usage

```
voxel_mask_from_dat.sh <STRUCTURAL> <MRS> <OUTFILE>
```

Where `<STRUCTURAL>` is the filename of the structural (anatomical image), `<MRS>` is the name of the MRS `.dat` file, and `<OUTFILE>` is the desired output filename for the mask.



## Inner workings
Essentially what the script does is it first creates a mask file where the voxel is in a corner of the image (but has the correct size already) and then it's translated and rotated into the correct position.



