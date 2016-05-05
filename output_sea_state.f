C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine output_sea_state(ncid,idtim,idvar,nx,ny,ishot,t)
      use MASTER, only: cel,xvel,yvel,dep,NAN,xrun,yrun,grav
      use PARAMETERS, only: celmin,lonsub,latsub
      implicit none
      include 'netcdf.inc'

      integer nx,ny,ishot,i,j,ii,jj
      real*8 aa,t,tm(1)
      real*4 usub(nx,ny),vsub(nx,ny),hsub(nx,ny)
      integer ncid(3),idtim(3),idvar(3)
      integer tstart(1),tcount(1),vstart(3),vcount(3)

	do j=1,ny
		jj=latsub*(j-1)+1
		do i=1,nx
			ii=lonsub*(i-1)+1
			aa=cel(ii,jj)
			if(aa.gt.celmin) then
				hsub(i,j)=(aa**2)/grav-dep(ii,jj)
				usub(i,j)=xvel(ii,jj)
				vsub(i,j)=yvel(ii,jj)
			else
				hsub(i,j)=nan
				usub(i,j)=nan
				vsub(i,j)=nan
			endif
		enddo
	enddo

      tm(1)=t
      tstart(1)=ishot
      tcount(1)=1
      vstart=(/ 1,1,ishot /)
      vcount=(/ nx,ny,1 /)
      
      call errhandle(nf_put_vara_double(ncid(1),idtim(1),
     & tstart,tcount,tm))
      call errhandle(nf_put_vara_real(ncid(1),idvar(1),
     &     vstart,vcount,hsub))
      call errhandle(nf_sync(ncid(1)))    
     		if(xrun) then
      call errhandle(nf_put_vara_double(ncid(2),idtim(2),
     & tstart,tcount,tm))
      call errhandle(nf_put_vara_real(ncid(2),idvar(2),
     &     vstart,vcount,usub))
      call errhandle(nf_sync(ncid(2)))
     		endif
     		if(yrun) then
      call errhandle(nf_put_vara_double(ncid(3),idtim(3),
     & tstart,tcount,tm))
      call errhandle(nf_put_vara_real(ncid(3),idvar(3),
     &     vstart,vcount,vsub))
      call errhandle(nf_sync(ncid(3)))
		endif
      write(9,*) 'screenshot',ishot,'at',t,'sec'          
      write(9,*)' Max elevation ',maxval(hsub)
      call flush(9)
	
	ishot=ishot+1
	
      end subroutine output_sea_state
