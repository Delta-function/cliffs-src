C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine opennc4feed(bname,nX,nY,bndrids)
      implicit none
      include 'netcdf.inc'

      character*200 bname,fname
      integer  ncid,bndrids(4)
      integer idpnt,iduvq,idtim,valsid,timeid
      integer nX,nY,npnt,nuvq,xdims(3),tdims(1)

 	nuvq=3     
      fname=' '
            
      fname=trim(bname)//'_west.nc'
	npnt=nY
      call errhandle(nf_create(trim(fname), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid, 'pnt', npnt, idpnt))
      call errhandle(nf_def_dim(ncid, 'uvq', nuvq, iduvq))
      call errhandle(nf_def_dim(ncid, 'tim', NF_UNLIMITED, idtim))
! variables
      xdims = (/ idpnt, iduvq, idtim /)
      tdims(1) = idtim      
      call errhandle(nf_def_var(ncid, 'vals', NF_DOUBLE, 3, xdims, 
     &                  valsid))
      call errhandle(nf_def_var(ncid, 'time', NF_DOUBLE, 1, tdims, 
     &                  timeid))

      call errhandle(nf_enddef(ncid))
      bndrids(1)=ncid

      fname=trim(bname)//'_east.nc'
	npnt=nY
      call errhandle(nf_create(trim(fname), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid, 'pnt', npnt, idpnt))
      call errhandle(nf_def_dim(ncid, 'uvq', nuvq, iduvq))
      call errhandle(nf_def_dim(ncid, 'tim', NF_UNLIMITED, idtim))
! variables
      xdims = (/ idpnt, iduvq, idtim /)
      tdims(1) = idtim      
      call errhandle(nf_def_var(ncid, 'vals', NF_DOUBLE, 3, xdims, 
     &                  valsid))
      call errhandle(nf_def_var(ncid, 'time', NF_DOUBLE, 1, tdims, 
     &                  timeid))
      call errhandle(nf_enddef(ncid))
      bndrids(2)=ncid

      fname=trim(bname)//'_south.nc'
	npnt=nX
      call errhandle(nf_create(trim(fname), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid, 'pnt', npnt, idpnt))
      call errhandle(nf_def_dim(ncid, 'uvq', nuvq, iduvq))
      call errhandle(nf_def_dim(ncid, 'tim', NF_UNLIMITED, idtim))
! variables
      xdims = (/ idpnt, iduvq, idtim /)
      tdims(1) = idtim      
      call errhandle(nf_def_var(ncid, 'vals', NF_DOUBLE, 3, xdims, 
     &                  valsid))
      call errhandle(nf_def_var(ncid, 'time', NF_DOUBLE, 1, tdims, 
     &                  timeid))
      call errhandle(nf_enddef(ncid))
      bndrids(3)=ncid

      fname=trim(bname)//'_north.nc'
	npnt=nX
      call errhandle(nf_create(trim(fname), NF_CLOBBER, ncid))
! dimensions
      call errhandle(nf_def_dim(ncid, 'pnt', npnt, idpnt))
      call errhandle(nf_def_dim(ncid, 'uvq', nuvq, iduvq))
      call errhandle(nf_def_dim(ncid, 'tim', NF_UNLIMITED, idtim))
! variables
      xdims = (/ idpnt, iduvq, idtim /)
      tdims(1) = idtim      
      call errhandle(nf_def_var(ncid, 'vals', NF_DOUBLE, 3, xdims, 
     &                  valsid))
      call errhandle(nf_def_var(ncid, 'time', NF_DOUBLE, 1, tdims, 
     &                  timeid))
      call errhandle(nf_enddef(ncid))
      bndrids(4)=ncid
      
      end subroutine opennc4feed
