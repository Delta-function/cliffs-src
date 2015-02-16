c                           CLIFFS GLOBAL VARIABLES                                        
c Cliffs is an open-source model for tsunami propagation and inundation computations       
c in the framework of the shallow-water equations. Cliffs implements:                      
c (1) VTCS-2 finite-difference scheme in an open ocean and on the open boundary            
c	 as in (Titov and Synolakis, 1995, 1998);                                            
c (2) dimensional splitting as in (Titov and Gonzalez, 1997; Titov and Synolakis, 1998);   
c (3) reflection and inundation computations as in (Tolkova, 2014);                        
c (4) data flows as in curvilinear MOST and MOST-4 (Burwell and Tolkova, 2008).            
c                               REFERENCE:                                                             
c E. Tolkova. Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting 
c Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314                     
c                                                                                          
c     Copyright (C) 2014, Elena Tolkova                                                    
c     For conditions of distribution and use, see copyright notice in cliffs_main.f        

	module MASTER
      real*8, parameter :: PI=3.1415926536, EarthR=6371009, grav=9.8
      real*4, parameter :: NAN=-1.e32
! Master grid arrays
      integer nXn,nYn,nXYn                 
      real*8, dimension(:,:), allocatable :: cel,xvel,yvel,dep !(nXn,nYn)
      real*8, dimension(:,:), allocatable :: hmax, velmax !(nXn,nYn)
      real*8, dimension(:), allocatable :: Xcrd, Ycrd !(nXn),(nYn)
      real*8, dimension(:,:), allocatable :: s1
      real*8, dimension(:), allocatable :: s2
      real*8, dimension(:), allocatable :: zeta
      logical xrun,yrun
! Master Grid input
      integer nTn,nXsrc,nYsrc
      real*8, dimension(:), allocatable :: tinput,Xsrc,Ysrc
      real*8, dimension(:,:), allocatable :: usrc,vsrc,qsrc
      real*8, dimension(:,:,:), allocatable :: west,east ! (nYn,3,nTn)
      real*8, dimension(:,:,:), allocatable :: south,north ! (nXn,3,nTn)	
      real*8, dimension(:,:), allocatable :: edge1,edge2 !(max(nXn,nYn),3)
	end module MASTER
	
	module NESTED
      integer nests  ! up to 99 nested grids
      integer nXnst(99), nYnst(99), maxlen                
      integer, dimension(:,:), allocatable :: WE_ny,SN_nx 		!(nYn,nests),(nXn,nests)
      integer, dimension(:), allocatable :: W_nx,E_nx,S_ny,N_ny 	!(nests)
      real*8, dimension(:,:), allocatable :: WE_cy, SN_cx 		!(nYn,nests),(nXn,nests)
      real*8, dimension(:), allocatable :: W_cx,E_cx,S_cy,N_cy 	!(nests)
	end module NESTED

      module PARAMETERS
      integer seaout,gout,maxout,lonsub,latsub
      real*8 cuke,dwall,ground,crough,dt,celmin
      integer steps_total,bndout,itopo,quake
      logical freeze ! true to stop computations when bndr feed ends
      logical cartesian  ! true if cartesian coordinate system
      character*200 MBathyFile
      integer Ngages
      integer, dimension(:), allocatable :: Igages, Jgages
      end module PARAMETERS
