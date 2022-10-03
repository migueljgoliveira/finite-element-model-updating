ECHO OFF
cls

del displacements_to_strains.pyd >nul

python -m numpy.f2py -c -m displacements_to_strains displacements_to_strains.f90 --fcompiler=intelvem --opt=/heap-arrays:0 1> displacements_to_strains.log 2>&1

rename displacements_to_strains.*.pyd displacements_to_strains.pyd >nul
