C  This file is part of ComputeBoundaryInput                                                             
C                                                                                          
C     Copyright (C) 2015, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cmpboundary.f        

      subroutine scan4wave(fname,nx,ny,nt,ix1,iy1,it,cuke,error)

      implicit none
      include 'netcdf.inc'
      
      character*200 fname
      integer error,ncid,ndims,nvar,ngatts,nmore,ncerr
      integer it,ix1,iy1,nx,ny,nt, icount(3), istart(3)
      real*8 nanum,cuke,maxu
      real*4, dimension(:,:,:), allocatable :: sea	
      
      error=0
      nanum=-1.0e+10	
      ncerr = nf_open(trim(fname),0,ncid)
      if (ncerr .ne. NF_NOERR) then
      	write(9,*) 'not found'
      	error=-1
      	goto 99
      endif
      call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
	
	allocate(sea(nx,ny,1))
	icount=(/nx,ny,1/)
	do it=1,nt	
		istart=(/ix1,iy1,it/)
		call errhandle(nf_get_vara_real(ncid, nvar,
     & istart,icount,sea))
		where(sea < nanum) sea=0
		maxu=maxval(sea)
		if(maxu.gt.cuke) exit
		maxu=maxval(-sea)
		if(maxu.gt.cuke) exit
	enddo
      deallocate(sea)
	
	if(it.ge.nt) then
		write(*,*) 'no signal above treshold'
		error=-1
		goto 99
	else
		write(*,*) 'signal detected at step ', it
	endif
      call errhandle(nf_close(ncid))

 99   continue      
      end subroutine scan4wave
