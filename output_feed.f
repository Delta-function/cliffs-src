C  This file is part of Cliffs.                                                             

	subroutine output_feed(ncid,nn,mm,uvq,it,time)

      implicit none
      include 'netcdf.inc'
      
      real*8 uvq(mm,3),time,var(nn,3,1)
      integer nn,mm,ncid,it,timeid,varid
      integer tstart(1),tcount(1),start(3),count(3)
	
	timeid=2
	varid=1
      tstart(1)=it
      tcount(1)=1
      call errhandle(nf_put_vara_double(ncid,timeid,tstart,
     & tcount,time))
      count =(/ nn,3,1 /)
      start =(/ 1,1,it /)
      var(1:nn,1:3,1)=uvq(1:nn,1:3)
      call errhandle(nf_put_vara_double(ncid,varid,start,count,var))
      call errhandle(nf_sync(ncid))

      end subroutine output_feed
