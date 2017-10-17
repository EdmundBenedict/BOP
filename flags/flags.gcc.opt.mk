
fc := gfortran
fflags := -pipe -O3 -march=corei7-avx -flto -funroll-loops -funsafe-loop-optimizations -funsafe-math-optimizations

ldflags := -llapack -lblas -lsymspg

# How to set destination for the .mod files
modflag := -J


# 'Out of tree' build path
bldpath := bld/g/o
