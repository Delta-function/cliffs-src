C  This file is part of ComputeBoundaryInput                                                             
C                                                                                          
C     Copyright (C) 2015, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cmpboundary.f        

      subroutine get_sea_state(fname,sea,nx,ny,nnt,ix1,iy1,it,error,m)

      implicit none
      include 'netcdf.inc'
      
      character*200 fname
      integer error,ncid,ndims,nvar,ngatts,nmore,ncerr,m
      real*8 nanum
      real*4 sea(nx,ny,nnt)
      integer nx,ny,nnt,ix1,iy1,it, icount(3), istart(3)
      
      nanum=-1.0e+10	
      ncerr = nf_open(trim(fname),0,ncid)
      if (ncerr .ne. NF_NOERR) then
      	write(9,*) 'not found'
      	error=-1
      	goto 99
      endif
      call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
      
	icount=(/nx, ny, nnt/)
	istart=(/ix1,iy1,it/)
	call errhandle(nf_get_vara_real(ncid, nvar,
     & istart,icount,sea))
	where(sea < nanum) sea=0
	if(m.LT.0) sea=0.01*sea
      call errhandle(nf_close(ncid))
      write(*,*) trim(fname),' is read'

 99   continue      
      end subroutine get_sea_state
