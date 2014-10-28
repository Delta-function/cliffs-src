C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine opennc4gages(ncfn,gncid,timid,varid,
     &      Ngages,lon,lat,cartesian,xrun,yrun)
      implicit none
      include 'netcdf.inc'

      character*200 ncfn
      integer ncerr
      integer gncid,timid,varid(3),lonid,latid
      integer idpoint,idtime
      integer Ngages
      integer xdim(1),var_dims(2),Mdims(2)
      real*8  lon(1), lat(1)
      logical cartesian,xrun,yrun
      
      varid=(/ -1, -1, -1 /)
      ncerr = nf_create(trim(ncfn)//'_gages.nc', NF_CLOBBER, gncid)
      call errhandle(ncerr)
! dimensions
      call errhandle(nf_def_dim(gncid, 'point', Ngages, idpoint))
      call errhandle(nf_def_dim(gncid, 'time', NF_UNLIMITED, idtime))
! variables
      xdim(1) = idpoint
      if(cartesian) then
      	ncerr = nf_def_var(gncid,'xxx',NF_DOUBLE,1,xdim,lonid)
      	ncerr = nf_def_var(gncid,'yyy',NF_DOUBLE,1,xdim,latid)
      else
      	ncerr = nf_def_var(gncid,'lon',NF_DOUBLE,1,xdim,lonid)
      	ncerr = nf_def_var(gncid,'lat',NF_DOUBLE,1,xdim,latid)
	endif
      xdim(1)=idtime
      ncerr = nf_def_var(gncid, 'time', NF_DOUBLE,1,xdim,timid)
      Mdims(1)=idpoint
      Mdims(2)=idtime
      call errhandle(nf_def_var(gncid,'gage', NF_DOUBLE,2,
     &            Mdims,varid(1)))
      if(xrun) then      
      	call errhandle(nf_def_var(gncid,'u', NF_DOUBLE,2,
     & 		Mdims,varid(2)))
      end if
      if(yrun) then      
      	call errhandle(nf_def_var(gncid,'v', NF_DOUBLE,2,
     & 		Mdims,varid(3)))
      end if
! attributes
      ncerr = nf_put_att_text(gncid,varid(1),'name', 26,
     &       'wave elevation at the gage')
      ncerr = nf_put_att_text(gncid,varid(1),'units',6,'meters')
      if(xrun) then      
      ncerr = nf_put_att_text(gncid,varid(2),'name', 18,
     &       'x (zonal) velocity')
      ncerr = nf_put_att_text(gncid,varid(2),'units',3,'m/s')
      end if
      if(yrun) then      
      ncerr = nf_put_att_text(gncid,varid(3),'name', 23,
     &       'y (meridional) velocity')
      ncerr = nf_put_att_text(gncid,varid(3),'units',3,'m/s')
      end if
      if(cartesian) then
      	ncerr = nf_put_att_text(gncid,lonid,'name',6,'x-axis')
      	ncerr = nf_put_att_text(gncid,lonid,'units',6,'meters')
      	ncerr = nf_put_att_text(gncid,latid,'name',6,'y-axis')
      	ncerr = nf_put_att_text(gncid,latid,'units',6,'meters')
      else
            ncerr = nf_put_att_text(gncid,lonid,'name',9,'longitude')
      	ncerr = nf_put_att_text(gncid,lonid,'units',12,'degrees_east')
      	ncerr = nf_put_att_text(gncid,latid,'name',8,'latitude')
      	ncerr = nf_put_att_text(gncid,latid,'units',13,'degrees_north')
	endif
      ncerr = nf_put_att_text(gncid,timid,'name',4,'time')
      ncerr = nf_put_att_text(gncid,timid,'units',35,'seconds')
      ncerr = nf_put_att_text(gncid,NF_GLOBAL,'model',6,'Cliffs')
      call errhandle(ncerr)      

      call errhandle(nf_enddef(gncid))
! gage coordinates
      call errhandle(nf_put_var_double(gncid,lonid,lon))
      call errhandle(nf_put_var_double(gncid,latid,lat))

	write(9,*) 'Time-series output to  '//trim(ncfn)//'_gages.nc'
      end subroutine opennc4gages
       
