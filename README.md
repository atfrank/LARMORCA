# CASA: Cα-based Solvent Accessibilities
  
- Predicts residue solvent accessibilities based on Cα coordinates using random decision trees.

## Install
```shell
$ cd /path/to/CASA/
$ make clean
$ make 
```

## Usage manual
```shell
$ bin/casa -h
Usage: casa [-options] <PDBfile>

Options:
         -trj : trajectory file (string)
         -skip : skip rate for reading trajectory (integer)
         -start : frame at which to start reading trajectory (integer)
         -stop : frame at which to stop reading trajectoryframe (integer)
         -identification : ID tag used in output (string)

```

## Examples
```shell
$ # predict solvent accessibilities from a coordinate file 
$ bin/casa tests/file.pdb (PDB format)
$
$ # predict chemical shifts from a trajectory file (DCD format) 
$ bin/casa -trj tests/file.dcd tests/file.pdb
```

## Output
### format
_trajectory-frame, residue-number, residue-name, sasa, id-tag_

### example
```shell
$ bin/casa -identification IDtag tests/file.pdb
1 1 GLY 122.846  IDtag
1 2 ARG 149.313  IDtag
1 3 GLY 72.4997  IDtag
1 4 LEU 136.496  IDtag
1 5 GLY 62.3456  IDtag
1 6 PRO 44.3154  IDtag
1 7 LEU 44.8277  IDtag
1 8 GLN 71.9106  IDtag
1 9 ILE 21.4488  IDtag
1 10 TRP 17.9812  IDtag
1 11 GLN 55.4005  IDtag
1 12 THR 36.2837  IDtag
1 13 ASP 35.8978  IDtag
1 14 PHE 39.6585  IDtag
1 15 THR 34.9928  IDtag
1 16 LEU 32.4873  IDtag
1 17 GLU 58.6507  IDtag
1 18 PRO 26.4722  IDtag
1 19 ARG 85.2815  IDtag
1 20 MET 16.2541  IDtag
1 21 ALA 24.2885  IDtag
1 22 PRO 51.0041  IDtag
1 23 ARG 62.8826  IDtag
1 24 SER 10.2901  IDtag
1 25 TRP 14.1995  IDtag
1 26 LEU 4.68346  IDtag
1 27 ALA 4.64075  IDtag
1 28 VAL 5.02193  IDtag
1 29 THR 15.2532  IDtag
1 30 VAL 5.16215  IDtag
1 31 ASP 17.4655  IDtag
1 32 THR 20.1837  IDtag
1 33 ALA 29.7804  IDtag
1 34 SER 73.6202  IDtag
1 35 SER 64.7972  IDtag
1 36 ALA 18.918  IDtag
1 37 ILE 18.4492  IDtag
1 38 VAL 13.0663  IDtag
1 39 VAL 17.029  IDtag
1 40 THR 19.3027  IDtag
1 41 GLN 56.7396  IDtag
1 42 HIS 41.1943  IDtag
1 43 GLY 27.3733  IDtag
1 44 ARG 92.6463  IDtag
1 45 VAL 26.2965  IDtag
1 46 THR 30.85  IDtag
1 47 SER 25.1629  IDtag
1 48 VAL 49.4112  IDtag
1 49 ALA 10.9413  IDtag
1 50 ALA 4.48318  IDtag
1 51 GLN 46.5995  IDtag
1 52 HIS 80.4594  IDtag
1 53 HIS 41.4464  IDtag
1 54 TRP 29.2906  IDtag
1 55 ALA 51.4916  IDtag
1 56 THR 47.8386  IDtag
1 57 ALA 6.27355  IDtag
1 58 ILE 20.5302  IDtag
1 59 ALA 78.8612  IDtag
1 60 VAL 56.856  IDtag
1 61 LEU 47.8338  IDtag
1 62 GLY 48.2455  IDtag
1 63 ARG 123.272  IDtag
1 64 PRO 18.3451  IDtag
1 65 LYS 82.1203  IDtag
1 66 ALA 5.84372  IDtag
1 67 ILE 4.05915  IDtag
1 68 LYS 48.4715  IDtag
1 69 THR 12.7724  IDtag
1 70 ASP 38.6929  IDtag
1 71 ASN 56.2789  IDtag
1 72 GLY 45.2089  IDtag
1 73 SER 56.1783  IDtag
1 74 CYS 15.419  IDtag
1 75 PHE 16.0736  IDtag
1 76 THR 63.5141  IDtag
1 77 SER 24.4761  IDtag
1 78 LYS 151.046  IDtag
1 79 SER 43.9767  IDtag
1 80 THR 19.8881  IDtag
1 81 ARG 61.1255  IDtag
1 82 GLU 88.7863  IDtag
1 83 TRP 19.6943  IDtag
1 84 LEU 9.76135  IDtag
1 85 ALA 48.5693  IDtag
1 86 ARG 155.987  IDtag
1 87 TRP 68.4943  IDtag
1 88 GLY 50.7  IDtag
1 89 ILE 26.8913  IDtag
1 90 ALA 27.7666  IDtag
1 91 HIS 67.5432  IDtag
1 92 THR 40.624  IDtag
1 93 THR 40.6842  IDtag
1 94 GLY 54.7335  IDtag
1 95 ILE 110.908  IDtag
1 96 PRO 121.547  IDtag
1 97 GLY 71.9464  IDtag
1 98 GLN 69.2684  IDtag
1 99 ALA 71.2321  IDtag
1 100 MET 40.4118  IDtag
1 101 VAL 12.5065  IDtag
1 102 GLU 79.3639  IDtag
1 103 ARG 106.515  IDtag
1 104 ALA 11.4834  IDtag
1 105 ASN 21.7986  IDtag
1 106 ARG 114.12  IDtag
1 107 LEU 23.4228  IDtag
1 108 LEU 8.6411  IDtag
1 109 LYS 61.658  IDtag
1 110 ASP 64.6166  IDtag
1 111 LYS 56.622  IDtag
1 112 ILE 6.25301  IDtag
1 113 ARG 75.7197  IDtag
1 114 VAL 58.756  IDtag
1 115 LEU 8.8543  IDtag
1 116 ALA 8.39724  IDtag
1 117 GLU 84.7471  IDtag
1 118 GLY 62.5702  IDtag
1 119 ASP 91.8669  IDtag
1 120 GLY 60.2122  IDtag
1 121 PHE 52.5697  IDtag
1 122 MET 57.7558  IDtag
1 123 LYS 109.824  IDtag
1 124 ARG 70.1557  IDtag
1 125 ILE 20.3269  IDtag
1 126 PRO 58.378  IDtag
1 127 THR 70.9881  IDtag
1 128 SER 109.852  IDtag
1 129 LYS 103.318  IDtag
1 130 GLN 48.5196  IDtag
1 131 GLY 48.0237  IDtag
1 132 GLU 92.7902  IDtag
1 133 LEU 15.1281  IDtag
1 134 LEU 9.87831  IDtag
1 135 ALA 27.3036  IDtag
1 136 LYS 84.6839  IDtag
1 137 ALA 8.77687  IDtag
1 138 MET 11.9091  IDtag
1 139 TYR 67.4171  IDtag
1 140 ALA 62.0968  IDtag
1 141 LEU 41.4506  IDtag
1 142 ASN 40.7067  IDtag
1 143 HIS 144.485  IDtag
```
