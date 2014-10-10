C  read whatever initial conditions are provided                                           
C  and syncronize with the boundary input if any                                            
C                                                                                          
C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine get_initial_state(srcname,indir,irec,
     & hinput,uinput,vinput,error)

      use MASTER, only: nTn,nXsrc,nYsrc,tinput,Xsrc,Ysrc,usrc,vsrc,qsrc
      implicit none
      include 'netcdf.inc'
      
      character*200 fname,fname2,dname,srcname,indir
      integer i,j,n,nc,error,istart(3),icount(3)
      integer*4 ncids(3),tmid,idvar(3)
      integer xxid,yyid
      integer ncerr,ndims,ngatts,nvar,nmore,dlen,ncid
      integer ivardims,ivartype
      logical hinput,uinput,vinput
      real*4, dimension(:,:,:), allocatable :: sea
      real*8, dimension(:), allocatable :: tsrc
      integer ntsrc,irec,isrc
      real*8 maxu,minu,maxv,minv,maxq,minq,nanum

      error=0
      fname=' '
      fname2=' '
      fname=trim(indir)//trim(srcname)

!  Displacement	  
      fname2=trim(fname)//'_h.nc'
      write(9,*) 'Initial displacement: ', trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncids(3))
      if (ncerr .ne. NF_NOERR) then
      	write(9,*) 'not found'
      	else
		hinput=.true.
      	endif
! zonal/x current
      fname2=trim(fname)//'_u.nc'
      write(9,*) 'Initial x-velocity: ',trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncids(1))
      if (ncerr .ne. NF_NOERR) then
      	write(9,*) 'not found'
      	else
		uinput=.true.
      	endif
! meridional/y current	  
      fname2=trim(fname)//'_v.nc'
      write(9,*) 'Initial y-velocity: ',trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncids(2))
      if (ncerr .ne. NF_NOERR) then
      	write(9,*) 'not found'
      	else
		vinput=.true.
      	endif
     	
	if(hinput) then
     		ncid=ncids(3)
     		else if(uinput) then
     			ncid=ncids(1)
     			else if(vinput) then
     				ncid=ncids(2)
     				else 
     					if(nTn<1) error=-1
     			      	goto 99
     	end if	
      	
      call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
      do i=1,ndims
		call errhandle(nf_inq_dim(ncid,i,dname,dlen))
        	if (dname(1:3).eq.'tim') then
           		ntsrc = dlen
        else if((dname(1:3).eq.'xxx').or.(dname(1:3).eq.'lon')) then
           		nXsrc = dlen
        else if((dname(1:3).eq.'yyy').or.(dname(1:3).eq.'lat')) then
           		nYsrc = dlen
        	endif
      end do
      do i=1,nvar
      	call errhandle(nf_inq_varname(ncid,i,dname))
      	if (dname(1:3).eq.'tim') then
      		tmid = i
        else if((dname(1:3).eq.'xxx').or.(dname(1:3).eq.'lon')) then
      		xxid = i
        else if ((dname(1:3).eq.'yyy').or.(dname(1:3).eq.'lat')) then
      		yyid = i
        	endif
      end do

      allocate(Xsrc(nXsrc),Ysrc(nYsrc),tsrc(ntsrc))
! get lon,lat,time
      call errhandle(nf_get_var_double(ncid, xxid,Xsrc))
      call errhandle(nf_get_var_double(ncid, yyid,Ysrc))
      call errhandle(nf_get_var_double(ncid, tmid,tsrc))
     
! find when to start
	if(nTn.lt.2) then
      	allocate(tinput(1))
		irec=1
		isrc=1
		tinput(1)=tsrc(1)	
	else
		call pinpoint(tsrc,ntsrc,tinput(irec),isrc,maxq)
	endif
		
	write(9,*) 'Initial conditions read at ',
     & 	tsrc(isrc),' applied at ',tinput(irec)
! get initial conditions
	deallocate(tsrc)
	allocate(sea(nXsrc,nYsrc,1))
	istart=(/1,1,isrc/)
	icount=(/nXsrc,nYsrc,1/)
      nanum=-1.0e+10
      maxu=-99999	
      maxv=-99999	
      maxq=-99999	
      minu=-99999	
      minv=-99999	
      minq=-99999	
	allocate(usrc(nXsrc,nYsrc),vsrc(nXsrc,nYsrc),qsrc(nXsrc,nYsrc))	
	
	if(uinput) then
		call errhandle(nf_inq(ncids(1),ndims,nvar,ngatts,nmore))
		call errhandle(nf_get_vara_real(ncids(1), nvar,
     & istart,icount,sea))
         	call errhandle(nf_close(ncids(1)))
		usrc(1:nXsrc,1:nYsrc)=sea(1:nXsrc,1:nYsrc,1)
		where(usrc < nanum) usrc=0
		maxu=maxval(usrc)
		minu=-maxval(-usrc)
	endif
	
	if(vinput) then
		call errhandle(nf_inq(ncids(2),ndims,nvar,ngatts,nmore))
		call errhandle(nf_get_vara_real(ncids(2), nvar,
     & istart,icount,sea))
         	call errhandle(nf_close(ncids(2)))
		vsrc(1:nXsrc,1:nYsrc)=sea(1:nXsrc,1:nYsrc,1)
		where(vsrc < nanum) vsrc=0
		maxv=maxval(vsrc)
		minv=-maxval(-vsrc)
	endif

	if(hinput) then
		call errhandle(nf_inq(ncids(3),ndims,nvar,ngatts,nmore))
		call errhandle(nf_get_vara_real(ncids(3), nvar,
     & istart,icount,sea))
         	call errhandle(nf_close(ncids(3)))
		qsrc(1:nXsrc,1:nYsrc)=sea(1:nXsrc,1:nYsrc,1)
		where(qsrc < nanum) qsrc=0
		maxq=maxval(qsrc)
		minq=-maxval(-qsrc)
	endif

	deallocate(sea)
     
      write(9,*)' Max initial u/v/q',maxu,maxv,maxq
      write(9,*)' Min initial u/v/q',minu,minv,minq

 99   call flush(9)      
      end subroutine get_initial_state
