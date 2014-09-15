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

