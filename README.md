# X-ray diffractogram

Computes the X-ray scattering intensity from a molecular-dynamics trajectory. Atomic form factors (including anomalous dispersion at the chosen
photon energy) are supplied by xraydb through a small Python script. This requires python3 and xraydb (https://xraypy.github.io/XrayDB/installation.html) to be installed.


## Compilation:

gfortran compute_intensity.F90 -o compute_intensity

## Preparation

The program reads two input files (the user has to supply them) and writes one or more result files.

1. movie.xyz: Multi-frame XYZ trajectory. Coordinates are expected to be in Angstroms (converted to nm internally). The box is treated as infinite (no periodic wrapping) unless the box-reading line in subroutine read_movie is enabled.

2. q_vals.dat: List of q-values for which the profile is calculated, one value per line (unit nm^-1). An example file is provided.

## Run

./xray_intensity xxx

where xxx is x-ray photon energy in eV.


## Methods

Three ways of calculations are available, and can be switched by the logical parameter flags near the top of compute_intensity.F90:

- precise: Is default. Direct double sum over all atom pairs and frames. Slow, but the reference result. Output: intensity_precise.dat

- fast: Bins pair distances into a histogram (up to d_max=20nm) and uses a tabulated sinc. Faster in case of many atoms. Output: intensity_fast.dat

- use_debye_waller: Averages over the trajectory with a Debye–Waller factor. Faster for long movies, but assumes isotropic harmonic motion, a constant atom count and ordering, and a fixed box. Output: intensity_dbf.dat

Each output file contains two columns: q and I(q).

temporary files of "f_vals_xxx.dat" with xxx=pid are created and deleted in the current working directory.
