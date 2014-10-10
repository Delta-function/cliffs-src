C This file is part of Cliffs.                                                  !

      subroutine output_maxwave(ncid,idval)
      use MASTER, only: nXn,nYn,hmax,dep,NAN,grav
      use PARAMETERS, only: celmin
      implicit none
      include 'netcdf.inc'

      integer i,j
      real*4 hh(nXn,nYn)
      real*8 aa
      integer ncid,idval

	hh=NAN
	do j=1,nYn 
		do i=1,nXn
			aa=hmax(i,j)
			if(aa.gt.celmin) hh(i,j)=(aa**2)/grav-dep(i,j)
		enddo
	enddo		
      call errhandle(nf_put_var_real(ncid,idval,hh))
      call errhandle(nf_sync(ncid))

      end subroutine output_maxwave
