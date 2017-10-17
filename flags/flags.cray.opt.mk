fc := ftn
fflags := -DNOSYM -Difort_read_bug -O3 -xCORE-AVX-I

ldflags := 

# How to set destination for the .mod files
modflag := -module

# 'Out of tree' build path
bldpath := bld/ci/o
