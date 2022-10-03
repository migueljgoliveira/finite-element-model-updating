ECHO OFF
cls

del vfm.pyd 2>nul

copy vfm.f90 vfm_compile.f90 >nul

type ludwik_isotropic_hardening.f90 >> vfm_compile.f90

python -m numpy.f2py -c -m vfm vfm_compile.f90 --fcompiler=intelvem --opt=/heap-arrays:0 1> vfm.log 2>&1

rename vfm.*.pyd vfm.pyd 2>nul

del vfm_compile.f90 >nul
