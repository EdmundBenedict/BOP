fc := mpif90
fflags := -DMPI -Difort_read_bug -O0 -g -traceback

ldflags := -mkl -lsymspg

# How to set destination for the .mod files
modflag := -module

# 'Out of tree' build path
bldpath := bld/mi/d
