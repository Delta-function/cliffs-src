###########################################################
# Makefile for Cliffs, Cartesian and Geo, 2D and 1D       #
# Author:      Elena Tolkova 03/2013-09/2014              #
###########################################################
F95	=  gfortran
EXE	=  Cliffs
MOD	=  global_modules.f
SOL	=  cliffs_main.f cliffs.f
GRD	=  MasterGrid.f NestedGrid.f
INI	=  allabout.f get_initial_state.f apply_initial_conditions.f
FED	= get_boundary_feed.f feed_children_grid.f errhandle.f
OUT	= output_sea_state.f output_maxwave.f output_gages.f output_feed.f
OPN	= opennc4sea.f opennc4max.f opennc4gages.f opennc4feed.f
OPT         =  -m64 -fopenmp -I/usr/local/gfortran/include/ -I/usr/local/include/ -L/usr/local/gfortran/lib/ -L/usr/local/lib/ -lnetcdf -lnetcdff -lcurl

$(EXE):  $(MOD) $(SOL) $(GRD) $(INI) $(FED) $(OUT) $(OPN)
	$(F95) $(MOD) $(SOL) $(GRD) $(INI) $(FED) $(OUT) $(OPN) $(OPT) -o $(EXE)
