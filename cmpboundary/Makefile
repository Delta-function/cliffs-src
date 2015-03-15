###########################################################
# Makefile for ComputeBoundaryInput                       #
# Author:      Elena Tolkova 03/2015                      #
###########################################################
F95	=  gfortran
EXE	=  ComputeBoundaryInput
SOL	=  cmpboundary.f 
GRD	=  ReadGridLattice.f ReadSeaLattice.f
INI	=  get_sea_state.f scan4wave.f InterpolateOntoBndr.f
OPN	=  errhandle.f opennc4feed.f output_bndr.f
OPT         =  -m64 -I/usr/local/gfortran/include/ -I/usr/local/include/ -L/usr/local/gfortran/lib/ -L/usr/local/lib/ -lnetcdf -lnetcdff

$(EXE):  $(SOL) $(GRD) $(INI) $(OPN)
	$(F95) $(SOL) $(GRD) $(INI) $(OPN) $(OPT) -o $(EXE)
