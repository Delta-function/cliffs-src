###########################################################
# Makefile for depth_ssl: node-wise changes to a depth map
#       to limit depth variations. 
# Improves stability when using VTCS-2 scheme
# Author:  Elena Tolkova, e.tolkova@gmail.com   30/09/2014   
###########################################################
F95         =  gfortran
OBJ         =  depth_ssl
SRC         =  depth_ssl.f errhandle.f 
OPT         =  -O -I/usr/local/gfortran/include/ -I/usr/local/include/ -L/usr/local/gfortran/lib/ -L/usr/local/lib/ -lnetcdf -lnetcdff

$(OBJ):  $(SRC)
	$(F95) $(SRC) $(OPT) -o $(OBJ)
