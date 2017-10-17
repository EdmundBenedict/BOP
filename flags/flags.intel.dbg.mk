fc := ifort
fflags := -O0 -g -traceback

ldflags := -mkl -lsymspg

# How to set destination for the .mod files
modflag := -module


# 'Out of tree' build path
bldpath := bld/i/d
