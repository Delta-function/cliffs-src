C  Compute_Boundary_Input is an open-source code for constructing boundary input 
C  into a target area given lon x lat x time dataset of sea states in a larger area.
C  The resulting boundary input can be used with Cliffs, to drive tsunami simulation 
C  in the target area.
C
C           Copyright (c) 2015, Elena Tolkova                                              
C           under the terms of FreeBSD License                                             
C                                                                                          
C  Redistribution and use in source and binary forms, with or without                       
C  modification, are permitted provided that the following conditions are met:              
C                                                                                          
C   1. Redistributions of source code must retain the above copyright notice, this          
C     list of conditions and the following disclaimer.                                     
C   2. Redistributions in binary form must reproduce the above copyright notice,            
C     this list of conditions and the following disclaimer in the documentation            
C     and/or other materials provided with the distribution.                               
C                                                                                          
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND          
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED            
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                   
C  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR          
C  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES           
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;             
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND              
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT               
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS            
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                             

	module VARS
      real*8, dimension(:), allocatable :: xnst, ynst !(mm),(nn)
      real*8, dimension(:), allocatable :: Xsrc, Ysrc, tsrc !(nXn),(nYn)
      integer mm,nn,nXn,nYn,nt
      end module VARS
      
      module Interp
      integer, dimension(:), allocatable :: WE_ny,SN_nx !(nYn),(nXn)
      real*8, dimension(:), allocatable :: WE_cy, SN_cx !(nYn),(nXn)
      integer W_nx,E_nx,S_ny,N_ny
      real*8 W_cx,E_cx,S_cy,N_cy
      real*8, dimension(:,:,:), allocatable :: west,east ! (nn,3,nt)
      real*8, dimension(:,:,:), allocatable :: south,north ! (mm,3,nt)	
      end module Interp	
	
      program Compute_Boundary_Input
	use VARS
	use Interp
      implicit none
      include 'netcdf.inc'

      character*200 pathin,casetitle,grid,fnm,units
      integer error,i,status, feedid(4), met
      real*8 cuke	
      integer irec,ix1,ix2,iy1,iy2,nx,ny,nnt,it,ic1,ic2,ic3
      real*8 x,y,xx,yy
      real*4, dimension(:,:,:), allocatable :: src
      real*8, dimension(:), allocatable :: tim !(nnt)
                       
      irec=iargc()
      if (irec.ne.4) then
      write(*,*) 'Help: */ComputeBoundaryInput <full path to h-file> 
     & <full path to target grid> <output dir/casetitle> <units(m/cm)>'
         goto 99
      endif  

      pathin=' '
      call getarg(1,pathin)      
      grid=' '
      call getarg(2,grid)
      casetitle=' '
      call getarg(3,casetitle)
      units=' '
      call getarg(4,units)
      met=-1
      if(trim(units).eq.'m') met=1
      
      write(*,*)'Still sea threshold (m) ? '
      read(*,*) cuke
      write(*,*) 'Still sea threshold (m): ', cuke      
      
! read Master grid lattice      
	write(*,*) 'Target area bathy/topo file:  ',trim(grid)
	call ReadGrid(grid,error)
	if (error.ne.0) goto 99
! done read Master grid lattice      
! read source grid lattice      
	write(*,*) 'Sea state file:  ',trim(pathin)
	call ReadSea(pathin,nt,nXn,nYn,tsrc,Xsrc,Ysrc,error)
	if (error.ne.0) goto 99
! done read source grid lattice      

! check nested configuration
	x=xnst(1)
	xx=xnst(mm)
	y=ynst(1)
	yy=ynst(nn)
      if((x.lt.Xsrc(1)).or.(xx.gt.Xsrc(nXn)).or.
     & 	(y.lt.Ysrc(1)).or.(yy.gt.Ysrc(nYn))) then
      	error=-1
		write(*,*) 'Target area lon/lat limits: ',x,xx,y,yy		
      	write(*,*) 'not within source grid limits'
      	goto 99
      endif
	
! keep left-next parent cell and interpolation coefficient for each child edge node
      allocate(WE_ny(nn),SN_nx(mm),WE_cy(nn),SN_cx(mm))
! North and South
	call pinpoint(Ysrc,nYn,ynst(1),S_ny,S_cy)
	call pinpoint(Ysrc,nYn,ynst(nn),N_ny,N_cy)		
	do i=1,mm
     	   call pinpoint(Xsrc,nXn,xnst(i),SN_nx(i),SN_cx(i))
     	enddo
! East and West
	call pinpoint(Xsrc,nXn,xnst(1),W_nx,W_cx)
	call pinpoint(Xsrc,nXn,xnst(mm),E_nx,E_cx)		
	do i=1,nn
     	   call pinpoint(Ysrc,nYn,ynst(i),WE_ny(i),WE_cy(i))
     	enddo
     	
! scan target area within source grid for the wave
      ix1=W_nx
      ix2=E_nx+1
	iy1=S_ny
	iy2=N_ny+1	
	nx=ix2-ix1+1
	ny=iy2-iy1+1
	call scan4wave(pathin,nx,ny,nt,ix1,iy1,it,cuke,error)     	
     	nnt=nt-it+1 !tsrc length from the start of the signal
     	allocate(tim(nnt))
     	tim=tsrc(it:nt)
     	deallocate(tsrc)
     	
! allocate boundary input arrays
	allocate(west(nn,3,nnt),east(nn,3,nnt),stat=status)
	if(status.gt.0) then
		write(*,*) 'memory allocation unsuccessful'
		goto 99
	endif
 	allocate(south(mm,3,nnt),north(mm,3,nnt),stat=status)
	if(status.gt.0) then
		write(*,*) 'memory allocation unsuccessful'
		goto 99
	endif
	
	
! get sea state
      ic1=index(pathin,'/',.true.)+1
      ic2=index(pathin,'h',.true.)
      ic3=len_trim(pathin)
      if((ic2.ge.ic3).or.(ic1.ge.ic2)) then
      	write(*,*) 'Wrong naming for input dataset'
      	goto 99
      endif
	allocate(src(nx,ny,nnt),stat=status)
	if(status.gt.0) then
		write(*,*) 'memory allocation unsuccessful'
		goto 99
	endif

!  Displacement	  
	call get_sea_state(pathin,src,nx,ny,nnt,ix1,iy1,it,error,met)
	if (error.ne.0) goto 99
	call interpolate(src,nx,ny,nnt,mm,nn,3)	
! zonal/x current
	pathin(ic2:ic2)='u'
	call get_sea_state(pathin,src,nx,ny,nnt,ix1,iy1,it,error,met)
	if (error.ne.0) goto 99
	call interpolate(src,nx,ny,nnt,mm,nn,1)	
! meridional/y current	  
	pathin(ic2:ic2)='v'
	call get_sea_state(pathin,src,nx,ny,nnt,ix1,iy1,it,error,met)
	if (error.ne.0) goto 99
	call interpolate(src,nx,ny,nnt,mm,nn,2)	

	deallocate(src)
	
! open netcdf files for boundary input
      ic1=index(grid,'/',.true.)+1
      ic2=len_trim(grid)-3
	fnm=' '
	fnm=trim(casetitle)//'_'//grid(ic1:ic2)
	call opennc4feed(fnm,mm,nn,feedid)
	call output_bndr(feedid(1),nn,west,nnt,tim)	
	call output_bndr(feedid(2),nn,east,nnt,tim)	
	call output_bndr(feedid(3),mm,south,nnt,tim)	
	call output_bndr(feedid(4),mm,north,nnt,tim)	
	
C	! save initial sea state at t(it)

 99  	stop
      end program Compute_Boundary_Input


! .....................................

      subroutine pinpoint(xx,nx,xp,ip,cp)
      integer nx,ip,ipp,med
      real*8 xx(nx),xp,cp

      ip=1
      ipp=nx
      do while (ip+2.le.ipp)
         med=(ip+ipp)/2
         if (xp.gt.xx(med)) then
            ip=med
         else
            ipp=med
         endif
      enddo
      if (nx.gt.1) then
		cp=(xx(ip+1)-xp)/(xx(ip+1)-xx(ip))
	else
		cp=0.5
	endif
	end subroutine pinpoint


            

