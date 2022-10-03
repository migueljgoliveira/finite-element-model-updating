ECHO OFF
cls

del femu.pyd 2>nul

python -m numpy.f2py -c -m femu femu.f90 --fcompiler=intelvem --opt=/heap-arrays:0 1> femu.log 2>&1

rename femu.*.pyd femu.pyd >nul
