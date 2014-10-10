C  This file is part of Cliffs.                                                             
                                                                                          
      subroutine output_gages(gncid,gtim_id,gvar_id,
     &                       irec,t,Ngages,gages)
      implicit none
      include 'netcdf.inc'

      real*8 gages(Ngages,1),t
      integer irec
      integer Ngages,i
      integer gncid,gtim_id,gvar_id
      double precision tau(1)
      integer tstart(1),tcount(1),vstart(2),vcount(2)
      
    	irec=irec+1
      tau(1)=t
      tstart(1)=irec
      tcount(1)=1
      call errhandle(nf_put_vara_double(gncid,gtim_id,
     &      tstart,tcount,tau))
      vstart=(/ 1,irec /)
      vcount=(/ Ngages,1 /)
      call errhandle(nf_put_vara_double(gncid,gvar_id,
     &      vstart,vcount,gages))
      call errhandle(nf_sync(gncid))

      end subroutine output_gages
