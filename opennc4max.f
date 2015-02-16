C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine opennc4max(ncfn,ncid,idvals,cartesian)
      use MASTER, only: Xcrd,Ycrd,nXn,nYn
      implicit none
      include 'netcdf.inc'

      character*200 ncfn
      character*20 xname,yname
      integer ncid,xxid,yyid
      integer idx,idy,idvals(2)
      integer xdim(1),ydim(1),vdims(2)
      logical cartesian

	if(cartesian) then
		xname='xxx'
		yname='yyy'
	else
		xname='lon'
		yname='lat'
	endif

      call errhandle(nf_create(trim(ncfn), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid,trim(xname), nXn, idx))
      call errhandle(nf_def_dim(ncid,trim(yname), nYn, idy))
! variables
      xdim(1) = idx
      ydim(1) = idy
      call errhandle(nf_def_var(ncid,trim(xname),NF_DOUBLE,1,
     &       xdim, xxid))
      call errhandle(nf_def_var(ncid,trim(yname),NF_DOUBLE,1,
     &	 ydim, yyid))
	vdims(1) = idx
      vdims(2) = idy
      call errhandle(nf_def_var(ncid,'MaxE',NF_REAL,2,
     &        vdims,idvals(1)))
      call errhandle(nf_def_var(ncid,'MaxV',NF_REAL,2,
     &        vdims,idvals(2)))
! attributes
      call errhandle(nf_put_att_text(ncid, idvals(1), 'name', 23,
     &       'Maximal Water Elevation'))
      call errhandle(nf_put_att_text(ncid,idvals(1),'units',1,'m'))
      call errhandle(nf_put_att_text(ncid, idvals(2), 'name', 15,
     &       'Maximal Current'))
      call errhandle(nf_put_att_text(ncid,idvals(2),'units',3,'m/s'))
	if(cartesian) then
      call errhandle(nf_put_att_text(ncid,xxid,'name',6,'x-axis'))
      call errhandle(nf_put_att_text(ncid,xxid,'units',6,'meters'))
      call errhandle(nf_put_att_text(ncid,yyid,'name',6,'y-axis'))
      call errhandle(nf_put_att_text(ncid,yyid,'units',6,'meters'))
     		else
      call errhandle(nf_put_att_text(ncid,xxid,'name',9,'longitude'))
      call errhandle(nf_put_att_text(ncid,xxid, 
     & 'units', 12, 'degrees East'))
      call errhandle(nf_put_att_text(ncid,yyid,'name',8,'latitude'))
      call errhandle(nf_put_att_text(ncid,yyid, 
     & 'units', 13,'degrees North'))
	endif
      
      call errhandle(nf_enddef(ncid))
	call errhandle(nf_put_var_double(ncid, xxid, Xcrd))
      call errhandle(nf_put_var_double(ncid, yyid, Ycrd))
      call errhandle(nf_sync(ncid))
      
      write(9,*) 'Maxwave output to: ',trim(ncfn)
      end subroutine opennc4max
       
