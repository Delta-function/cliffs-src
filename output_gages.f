C  This file is part of Cliffs.                                                             
                                                                                          
      subroutine output_gages(gncid,gtim_id,gvar_id,
     &                       irec,t,Ngages)

      use MASTER, only: dep,cel,xvel,yvel,grav,nan
      use PARAMETERS, only: Igages,Jgages,celmin
      implicit none
      include 'netcdf.inc'

      real*8 gages(Ngages,1),ugg(Ngages,1),vgg(Ngages,1),t
      integer irec
      integer Ngages,j
      integer gncid,gtim_id,gvar_id(3)
      real*8 tau(1),aa
      integer tstart(1),tcount(1),vstart(2),vcount(2)
      
		do j=1,Ngages
			aa=cel(Igages(j),Jgages(j))
			if(aa.gt.celmin) then
				gages(j,1)=(aa**2)/grav - dep(Igages(j),Jgages(j))
			else
				gages(j,1)=nan
			endif
		end do
		forall(j=1:Ngages)
			ugg(j,1)=xvel(Igages(j),Jgages(j))
			vgg(j,1)=yvel(Igages(j),Jgages(j))
		end forall

      irec=irec+1
      tau(1)=t
      tstart(1)=irec
      tcount(1)=1
      call errhandle(nf_put_vara_double(gncid,gtim_id,
     &      tstart,tcount,tau))
      vstart=(/ 1,irec /)
      vcount=(/ Ngages,1 /)
      call errhandle(nf_put_vara_double(gncid,gvar_id(1),
     &      vstart,vcount,gages))
      if(gvar_id(2).gt.0) then
      call errhandle(nf_put_vara_double(gncid,gvar_id(2),
     &      vstart,vcount,ugg))
      end if
      if(gvar_id(3).gt.0) then
      call errhandle(nf_put_vara_double(gncid,gvar_id(3),
     &      vstart,vcount,vgg))
      end if
      call errhandle(nf_sync(gncid))

      end subroutine output_gages
