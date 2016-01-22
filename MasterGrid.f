C  Read Master grid bathy/topo and compute cell spacing and depression     
C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine MasterGrid(indir,fnmgrid,cartesian,error)

      use MASTER, only: nXn,nYn,dep,Xcrd,Ycrd,s1,s2,
     &	 zeta,PI,EarthR,nXYn
      implicit none
      include 'netcdf.inc'

      character*200 fname,indir,fnmgrid
	logical cartesian
      integer ix1,ix2,i,j,error
      real*8 alpha,beta   
      integer ncid,xcrdtyp,ycrdtyp,topotyp
      real*4, dimension(:,:), allocatable :: r4d2temp
	real*8, dimension(:), allocatable :: r8d1temp
	real*4, dimension(:), allocatable :: r4d1temp
      integer, parameter :: idxcr=1,idycr=2,idbth=3
      integer ndims,nvar,ngatts,nmore   ! added for rivers

      fname=' '
      read(5,2) fname
      ix1=index(fname,'/',.true.)+1
      ix2=index(fname,'.nc')-1
      fnmgrid=' '
      fnmgrid=trim(fname(ix1:ix2))
      fname = trim(indir)//trim(fname)

	write(9,*) 'Master bathy/topo file:  ',trim(fname)
      call errhandle(nf_open(trim(fname),0,ncid))
      call errhandle(nf_inq_dimlen(ncid,1,nXn))
      call errhandle(nf_inq_dimlen(ncid,2,nYn))
	write(9,*) 'grid size ',nXn,nYn

      allocate(Xcrd(nXn),Ycrd(nYn))
      allocate(dep(nXn,nYn))
      nXYn=maxval( (/nXn, nYn/) )
      allocate(r4d1temp(nXYn))

      call errhandle(nf_inq_vartype(ncid,idycr,ycrdtyp))
      call errhandle(nf_inq_vartype(ncid,idxcr,xcrdtyp))
      call errhandle(nf_inq_vartype(ncid,idbth,topotyp))

      if (xcrdtyp .eq. NF_FLOAT) then
         call errhandle(nf_get_var_real(ncid,idxcr,r4d1temp))
         Xcrd(1:nXn)=r4d1temp(1:nXn)
      else if (xcrdtyp .eq. NF_DOUBLE) then
         call errhandle(nf_get_var_double(ncid,idxcr,Xcrd))
      else 
      	error=-1
      	write(9,*) 'x-vector unexpected type',xcrdtyp
		goto 99
      endif

      if (ycrdtyp .eq. NF_FLOAT) then
         call errhandle(nf_get_var_real(ncid,idycr,r4d1temp))
         Ycrd(1:nYn)=r4d1temp(1:nYn)
      else if (ycrdtyp .eq. NF_DOUBLE) then
         call errhandle(nf_get_var_double(ncid,idycr,Ycrd))
      else 
      	error=-1
      	write(9,*) 'y-vector unexpected type',ycrdtyp
      	goto 99
      endif
      
	deallocate(r4d1temp)
      
      if (topotyp .eq. NF_FLOAT) then
         allocate(r4d2temp(nXn,nYn))
         call errhandle(nf_get_var_real(ncid,idbth,r4d2temp))
         dep(1:nXn,1:nYn)=r4d2temp(1:nXn,1:nYn)
         deallocate(r4d2temp)      
      else if (topotyp .eq. NF_DOUBLE) then
         call errhandle(nf_get_var_double(ncid,idbth,dep))
      else 
      	error=-1
      	write(9,*) 'bathy/topo unexpected type',topotyp
          	goto 99
      endif
		 
	beta=PI*EarthR/180.
	if(nXn.gt.2) then ! distances between nodes along a parallel
	      allocate(s1(nXn-1,nYn),r8d1temp(nXn-1))
		if(cartesian) then
			r8d1temp(1:(nXn-1))=Xcrd(2:nXn)-Xcrd(1:(nXn-1))		
			do j=1,nYn
				s1(1:(nXn-1),j)=r8d1temp(1:(nXn-1))
			enddo
		else
			r8d1temp(1:(nXn-1))=beta*(Xcrd(2:nXn)-Xcrd(1:(nXn-1)))
			do j=1,nYn
            		alpha=COS(PI/180.*Ycrd(j))     
				s1(1:(nXn-1),j)=alpha*r8d1temp(1:(nXn-1))
			enddo
		endif
		deallocate(r8d1temp)
		write(9,*) 'range of x/lon spacing, m:',minval(s1), maxval(s1)
	endif
	if(nYn.gt.2) then ! distances between nodes along a meridian, and scaling factors (dg/ds)/g/8 
      	allocate(s2(nYn-1),zeta(nYn))
		if(cartesian) then
			s2(1:(nYn-1))=Ycrd(2:nYn)-Ycrd(1:(nYn-1))
			zeta(1:nYn)=0
		else
			s2(1:(nYn-1))=beta*(Ycrd(2:nYn)-Ycrd(1:(nYn-1)))
			alpha=0.125/EarthR
			do i=1,nYn
				zeta(i)=alpha*TAN(PI/180.*Ycrd(i))
			enddo
		endif
		write(9,*) 'range of y/lat spacing, m:',minval(s2), maxval(s2)
	endif

C for rivers, y-dim only: extra bathy line - ggs
	call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
	if(nvar.gt.idbth) 
     &	call errhandle(nf_get_var_double(ncid,nvar,zeta))
		
      call errhandle(nf_close(ncid))

 2	format(a)
 99  	continue
      end subroutine MasterGrid
