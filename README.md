# streaminginstability
Research code for upcoming paper on lunar formation via the streaming instability. 
All code is meant to work on the output of Athena, which can be found at 
https://princetonuniversity.github.io/Athena-Cversion/.  
Nothing here is guaranteed to work 100%: keep your bug swatter handy.

## Included files
All runnable scripts can be run with `-h` or `--help` to print usage

plotutils.py: general function utility file, not runnable

calcsurfdens.py: calculates the surface density at all points for 2D (parallelized) simulation data

carrera.py: outputs plots inspired by (or other info related to) figure A1 from Carrera et. al. 2015

abod_mass.py: calculates predicted mass of planetesimals from overdensities as in Abod et. al. 2019

spacetime.py: calculates space vs. time vs. particle density colorplots
