C  This file is part of ComputeBoundaryInput                                                             

	subroutine output_bndr(ncid,nn,uvq,nt,time)

      implicit none
      include 'netcdf.inc'
      
      real*8 uvq(nn,3,nt),time(nt)
      integer nn,ncid,nt,timeid,varid
      integer tstart(1),tcount(1),start(3),count(3)
	
	timeid=2
	varid=1
      tstart(1)=1
      tcount(1)=nt
      call errhandle(nf_put_vara_double(ncid,timeid,tstart,
     & tcount,time))
      count =(/ nn,3,nt /)
      start =(/ 1,1,1 /)
      call errhandle(nf_put_vara_double(ncid,varid,start,count,uvq))
      call errhandle(nf_sync(ncid))

      end subroutine output_bndr
