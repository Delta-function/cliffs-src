C  This file is part of ComputeBoundaryInput                                                             
C                                                                                          
C     Copyright (C) 2015, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cmpboundary.f        

      subroutine ReadSea(fname,error)
      use VARS, only: nt,nXn,nYn,tsrc,Xsrc,Ysrc

      implicit none
      include 'netcdf.inc'
      
      character*200 fname,dname
      integer i,error
      integer ndims,ngatts,nvar,nmore,dlen,ncid
      integer xxid,yyid,ncerr,tmid
    
      error=0
      
      ncerr = nf_open(trim(fname),0,ncid)
      if (ncerr .ne. NF_NOERR) then
      	write(*,*) 'not found'
      	error=1
      	goto 99
      endif
      	
      call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
      do i=1,ndims
		call errhandle(nf_inq_dim(ncid,i,dname,dlen))
        	if ((dname(1:3).eq.'tim').or.(dname(1:3).eq.'TIM')) then
           		nt = dlen
        else if((dname(1:3).eq.'xxx').or.(dname(1:3).eq.'LON')
     &   .or.(dname(1:3).eq.'lon')) then
           		nXn = dlen
        else if((dname(1:3).eq.'yyy').or.(dname(1:3).eq.'LAT')
     &   .or.(dname(1:3).eq.'lat')) then
           		nYn = dlen
        	endif
      end do
      do i=1,nvar
      	call errhandle(nf_inq_varname(ncid,i,dname))
      	if ((dname(1:3).eq.'tim').or.(dname(1:3).eq.'TIM')) then
      		tmid = i
        else if((dname(1:3).eq.'xxx').or.(dname(1:3).eq.'LON')
     &   .or.(dname(1:3).eq.'lon')) then
      		xxid = i
        else if ((dname(1:3).eq.'yyy').or.(dname(1:3).eq.'LAT')
     &   .or.(dname(1:3).eq.'lat')) then
      		yyid = i
        	endif
      end do

      allocate(Xsrc(nXn),Ysrc(nYn),tsrc(nt))
! get lon,lat,time
      call errhandle(nf_get_var_double(ncid, xxid,Xsrc))
      call errhandle(nf_get_var_double(ncid, yyid,Ysrc))
      call errhandle(nf_get_var_double(ncid, tmid,tsrc))
      call errhandle(nf_close(ncid))

      write (*,*) 'Nodes in Sea State Dataset (lon/lat/time):', 
     &	 nXn,nYn,nt
      write (*,*) ' Longitude:',Xsrc(1),' to',Xsrc(nXn)
      write (*,*) '  Latitude:',Ysrc(1),' to',Ysrc(nYn)
      write (*,*) '      Time:',tsrc(1),' to',tsrc(nt)

 99   continue      
      end subroutine ReadSea
