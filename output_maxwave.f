C This file is part of Cliffs.                                                  !

      subroutine output_maxwave(ncid,idvals)
      use MASTER, only: nXn,nYn,hmax,dep,NAN,grav,velmax
      use PARAMETERS, only: celmin
      implicit none
      include 'netcdf.inc'

      integer i,j
      real*4 hh(nXn,nYn),vv(nXn,nYn)
      real*8 aa
      integer ncid,idvals(3)

	hh=NAN
	vv=NAN
	do j=1,nYn 
		do i=1,nXn
			aa=hmax(i,j)
			if(aa.gt.celmin) then
				hh(i,j)=(aa**2)/grav-dep(i,j)
				vv(i,j)=sqrt(velmax(i,j))
			endif
		enddo
	enddo		
      call errhandle(nf_put_var_real(ncid,idvals(1),hh))
      call errhandle(nf_put_var_real(ncid,idvals(2),vv))
      call errhandle(nf_sync(ncid))

      end subroutine output_maxwave
