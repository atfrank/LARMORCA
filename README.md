# LarmorCα: Protein Cα-based Chemical Shift Predictor
 
- Predicts Backbone Protein (HN, N, C, Hα, Cα, Cβ) chemical shifts
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

## Chemical shift list file
### file format
_residue-name, residue-number, nucleus, measured-shifts, error_

### example
```shell
$ head -10 tests/cs.dat
  
  MET 1 C 170.2 0
  MET 1 CA 53.3 0
  MET 1 HA 3.73 0
  ASN 2 C 174.3 0
  ASN 2 CA 51.7 0
  ASN 2 HA 4.28 0
  ASN 2 H 6.88 0
  ASN 2 N 115.7 0
  ILE 3 C 174.9 0
  ILE 3 CA 64.1 0
  ...
```

## output (chemical shifts)
### format
_trajectory-frame, residue-number, residue-name, nucleus, reference-shifts, predicted-shifts, measured-shifts, id-tag_

### example
```shell
$ bin/larmorca -csfile tests/cs.dat -identification ID tests/file.pdb
  
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

## output (errors)
### format
_trajectory-frame, error-HN, error-Hα, error-Cα, error-C, error-Cβ, error-N, error-total, id-tag_

### example
```shell
$ bin/larmorca -csfile tests/cs.dat -identification ID -printError -errorType MAE -trj tests/file.dcd tests/file.pdb 
  
  1 1.41216 0.22443 0.863814 0.923743 0 2.46342 5.88757 ID
  2 1.49718 0.298723 1.68938 1.25538 0 2.67265 7.41332 ID
  3 1.47061 0.279231 1.69115 1.24907 0 2.90211 7.59217 ID
  4 1.50804 0.295767 1.59486 1.24795 0 2.7518 7.39842 ID
  5 1.44317 0.287217 1.59268 1.19073 0 2.66489 7.17869 ID
  6 1.57143 0.280846 1.61944 1.23007 0 2.68759 7.38938 ID
  7 1.53636 0.286255 1.51561 1.18623 0 2.68117 7.20563 ID
  8 1.40482 0.306484 1.6368 1.25667 0 2.77738 7.38215 ID
  9 1.49456 0.296307 1.60117 1.19395 0 2.7616 7.34758 ID
  10 1.65376 0.303159 1.82448 1.31903 0 2.78786 7.8883 ID
  11 1.38815 0.278989 1.58481 1.23251 0 2.75756 7.24201 ID
  12 1.57407 0.261466 1.32346 1.08008 0 2.52171 6.7608 ID
  ...
```

