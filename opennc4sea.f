C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine opennc4sea(lon,lat,fname,who,ncid,tmid,
     &                  varid,nx,ny,cartesian)
      implicit none
      include 'netcdf.inc'

      character*200 fname
      character*20 xname,yname
      integer ncid,nx,ny,x1_id,x2_id,tim_id,who
      integer xxid,yyid,tmid,varid,xdim(1),var_dims(3)
      real*8  lon(nx),lat(ny)
      logical cartesian

	if(cartesian) then
		xname='xxx'
		yname='yyy'
	else
		xname='lon'
		yname='lat'
	endif

      call errhandle(nf_create(trim(fname), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid,trim(xname), nx, x1_id))
      call errhandle(nf_def_dim(ncid,trim(yname), ny, x2_id))
      call errhandle(nf_def_dim(ncid, 'time', NF_UNLIMITED, tim_id))
! variables
      xdim(1) = x1_id
      call errhandle(nf_def_var(ncid,trim(xname),NF_DOUBLE,1,
     & xdim,xxid))
      xdim(1) = x2_id
      call errhandle(nf_def_var(ncid,trim(yname),NF_DOUBLE,1,
     & xdim,yyid))
      xdim(1) = tim_id
      call errhandle(nf_def_var(ncid,'time',NF_DOUBLE,1,
     & xdim,tmid))

      var_dims = (/ x1_id, x2_id, tim_id /)

      select case(who)
      	case(1)
        call errhandle(nf_def_var(ncid,'ha',NF_REAL,3,var_dims,varid))
        call errhandle(nf_put_att_text(ncid,varid,'name',20,
     &         'Surface Displacement'))
        call errhandle(nf_put_att_text(ncid,varid,'units',6,'meters'))
      	case(2)
        call errhandle(nf_def_var(ncid,'ua',NF_REAL,3,var_dims,varid))
        call errhandle(nf_put_att_text(ncid,varid,'name',18,
     &         'U (zonal) Velocity'))
        call errhandle(nf_put_att_text(ncid,varid,'units',13,
     &        'meters/second'))
      	case(3)
        call errhandle(nf_def_var(ncid,'va',NF_REAL,3,var_dims,varid))
        call errhandle(nf_put_att_text(ncid, varid, 'name', 23,
     &         'V (meridional) Velocity'))
        call errhandle(nf_put_att_text(ncid, varid, 'units', 13,
     &        'meters/second'))
      end select

! attributes
	if(cartesian) then
      call errhandle(nf_put_att_text(ncid,xxid,'name',8,'subset_x'))
      call errhandle(nf_put_att_text(ncid,xxid,'units',6,'meters'))
      call errhandle(nf_put_att_text(ncid,yyid,'name',8,'subset_y'))
      call errhandle(nf_put_att_text(ncid,yyid,'units',6,'meters'))
 	else
 	      call errhandle(nf_put_att_text(ncid, xxid, 'name', 
     & 16, 'subset_longitude'))
      	call errhandle(nf_put_att_text(ncid, xxid, 'units', 
     & 12, 'degrees_east'))
      	call errhandle(nf_put_att_text(ncid, yyid, 'name', 
     & 15, 'subset_latitude'))
      	call errhandle(nf_put_att_text(ncid, yyid, 'units', 
     & 13, 'degrees_north'))
	endif    
      call errhandle(nf_put_att_text(ncid, tmid, 'name', 4, 'time'))
      call errhandle(nf_put_att_text(ncid, tmid, 'units', 
     & 7, 'seconds'))
      
      call errhandle(nf_enddef(ncid))
      call errhandle(nf_put_var_double(ncid, xxid, lon))
      call errhandle(nf_put_var_double(ncid, yyid, lat))
      
      write(9,*) 'Screenshot output to: ',trim(fname)
      end subroutine opennc4sea
       
