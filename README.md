# LarmorCa: Protein Ca-based Chemical Shift Predictor
 
- Predicts Backbone Protein (HN, N, C, HA, CA, CB) chemical shifts
- Predictors were generated using the Random Forest machine learning technique

## Install
```shell
$ cd /path/to/LARMORCA/
$ make clean
$ make 
```

## Usage manual
```shell
$ bin/larmorca -h
Usage: larmorca [-options] <PDBfile>

Options:
         -csfile : chemical shifts list file (string)
         -trj : trajectory file (string)
         -skip : skip rate for reading trajectory (integer)
         -start : frame at which to start reading trajectory (integer)
         -stop : frame at which to stop reading trajectoryframe (integer)
         -identification : ID tag used in output (string)
         -predictorType: which predictors to use; training or full (string: "train" or "full")
         -printError : output errors, rather than chemical shifts (flag)
         -errorType: if -printError, this specifies type of error to compute (string: "MAE", "RMSE", "wMAE", or "wRMSE")

```

## Examples
```shell
$ # predict chemical shifts from a coordinate file 
$ bin/larmorca -csfile tests/cs.dat -identification ID tests/file.pdb (PDB format)
$ # predict chemical shifts from a trajectory file (DCD format) 
$ bin/larmorca -csfile tests/cs.dat -identification ID -trj tests/file.dcd tests/file.pdb
```

## output
### format
_trajectory-frame, residue-number, residue-name, nucleus, reference-shifts, predicted-shifts, measured-shifts, id-tag_

### example
```shell
$ bin/larmorca -csfile tests/cs.dat -identification 1SCL tests/struct.pdb
  
1 1 HA MET 4.48 4.53644 3.73 ID
1 1 CA MET 55.4 54.94 53.3 ID
1 1 C MET 173.6 174.172 170.2 ID
1 2 H ASN 8.4 8.47179 6.88 ID
1 2 HA ASN 4.74 4.73207 4.28 ID
1 2 CA ASN 53.1 52.0674 51.7 ID
1 2 C ASN 175.2 175.192 174.3 ID
1 2 N ASN 118.7 122.102 115.7 ID
1 3 HA ILE 4.17 3.80757 4 ID
1 3 CA ILE 61.1 63.1395 64.1 ID
1 3 C ILE 176.4 177.449 174.9 ID
1 3 N ILE 119.9 119.898 117.4 ID
1 4 HA PHE 4.62 4.22713 3.6 ID
1 4 CA PHE 57.7 60.0242 61 ID
1 4 C PHE 175.8 177.382 176.3 ID
1 4 N PHE 120.3 118.762 122.9 ID
1 5 HA GLU 4.35 4.0578 3.63 ID
1 5 CA GLU 56.6 58.7379 58.1 ID
1 5 C GLU 176.6 178.694 178.3 ID
1 5 N GLU 120.2 118.591 117.7 ID
1 6 HA MET 4.48 4.02625 3.24 ID
1 6 CA MET 55.4 57.9186 58.7 ID
1 6 C MET 173.6 177.884 177 ID
1 6 N MET 119.6 119.895 118.2 ID
1 7 HA LEU 4.34 3.9542 3.95 ID
1 7 CA LEU 55.1 57.3532 56 ID
1 7 C LEU 177.6 178.538 178.9 ID
1 7 N LEU 121.8 118.766 118.5 ID
1 8 HA ARG 4.34 3.91707 3.3 ID
1 8 CA ARG 56 59.4871 60.2 ID
1 8 C ARG 176.3 178.189 178.8 ID
1 8 N ARG 120.5 119.398 123.8 ID
1 9 HA ILE 4.17 3.69591 3.49 ID
1 9 CA ILE 61.1 64.2022 65 ID
1 9 C ILE 176.4 178.337 177.7 ID
1 9 N ILE 119.9 118.095 119.7 ID
1 10 HA ASP 4.64 4.45631 4.49 ID
1 10 CA ASP 54.2 56.5941 57.8 ID
1 10 C ASP 176.3 177.661 177.3 ID
1 10 N ASP 120.4 117.88 117.7 ID
  ...
```

## Chemical shift list file
### file format
_residue-name, residue-number, nucleus, measured-shifts, error_

### example
```shell
$ head -10 tests/cs.dat
  
  CYT 6 C1' 91.438 0.00
  URA 7 C1' 92.388 0.00
  CYT 8 C1' 89.297 0.00
  ADE 9 C1' 85.504 0.00
  GUA 10 C1' 81.3 0.00
  GUA 1 C1' 89.246 0.00
  GUA 2 C1' 90.258 0.00
  GUA 3 C1' 90.657 0.00
  URA 4 C1' 90.835 0.00
  GUA 5 C1' 89.925 0.00
  ...
```

## LarmorD parameter file
### file format
_nucleus, neighbor-residue-name, neighbor-atom-name, alpha_

### example for H5''
```shell
$ head -68 data/parameters.txt

  H5'' GUA C1' 3.67815480447
  H5'' GUA C2' 1.63100734262
  H5'' GUA C3' 2.89281770448
  H5'' GUA C4' 2.3548429102
  H5'' GUA C5' 0.0321335101989
  H5'' GUA P -0.171617457218
  H5'' GUA O5' 1.57702037737
  H5'' GUA O3' 0.326168827858
  H5'' GUA C2 -6.64414148375
  H5'' GUA C4 -14.2514277865
  H5'' GUA C5 -18.4752206521
  H5'' GUA C6 -23.3660033205
  H5'' GUA C8 11.2256815971
  H5'' GUA N1 -15.0357687755
  H5'' GUA N2 18.8071249973
  H5'' GUA N3 7.49438758708
  H5'' GUA N7 6.47382380489
  H5'' GUA N9 -1.33876233081
  H5'' GUA O6 13.3592994749
  H5'' ADE C1' -6.87501637071
  H5'' ADE C2' -0.125715293229
  H5'' ADE C3' 0.934354550638
  H5'' ADE C4' 0.256912097145
  H5'' ADE C5' 0.571155032105
  H5'' ADE P -0.433798841887
  H5'' ADE O5' 1.5556382895
  H5'' ADE O3' 0.852536186953
  H5'' ADE C2 19.465835991
  H5'' ADE C4 -2.63457000853
  H5'' ADE C5 -6.20431462114
  H5'' ADE C6 -5.52564589385
  H5'' ADE C8 3.72779586352
  H5'' ADE N1 8.97415774656
  H5'' ADE N3 -3.65500648528
  H5'' ADE N6 -21.6733704635
  H5'' ADE N7 0.39668878159
  H5'' ADE N9 1.89834133117
  H5'' URA C1' -7.8361401694
  H5'' URA C2' -1.05475191362
  H5'' URA C3' -1.23351633568
  H5'' URA C4' 4.8561240664
  H5'' URA C5' 0.0
  H5'' URA P 1.86617611577
  H5'' URA O5' 1.01063675972
  H5'' URA O3' 0.949438715561
  H5'' URA C2 3.44987727435
  H5'' URA C4 -7.50237209988
  H5'' URA C5 -14.857005414
  H5'' URA C6 3.44656331225
  H5'' URA N1 2.71399930928
  H5'' URA N3 19.9073925365
  H5'' URA O4 23.4150612805
  H5'' CYT C1' 0.706640452452
  H5'' CYT C2' 0.166171606825
  H5'' CYT C3' -3.1035073294
  H5'' CYT C4' -1.57037444807
  H5'' CYT C5' 1.32626298398
  H5'' CYT P -0.440649244862
  H5'' CYT O5' -2.63138895134
  H5'' CYT O3' 1.92397549643
  H5'' CYT C2 8.20045630532
  H5'' CYT C4 1.29679400155
  H5'' CYT C5 4.77736668512
  H5'' CYT C6 1.05351000601
  H5'' CYT N1 2.78643038779
  H5'' CYT N3 6.38091065052
  H5'' CYT N4 -12.0144974986
  H5'' CYT O2 -2.71778957123
  ...
```

## Limitations
Currently, LarmorD does not predict chemical shifts for unprotonated 13C and 15N nuclei.

## Licence
```
...
```
