
fc := gfortran
fflags := -pipe -O0 -fbacktrace -finit-integer=2147483648 -finit-real=nan -finit-character=35 -fmax-stack-var-size=512 -fstack-protector-all

ldflags := -llapack -lblas -lsymspg

# How to set destination for the .mod files
modflag := -J


# 'Out of tree' build path
bldpath := bld/g/d
