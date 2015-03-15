C  This file is part of ComputeBoundaryInput                                                             
C                                                                                          
C     Copyright (C) 2015, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine ReadGrid(fname,error)
      use VARS, only: mm,nn,ynst,xnst
      implicit none
      include 'netcdf.inc'

      character*200 fname
      integer error,nXYn
      integer ncid,xnsttyp,ynsttyp
	real*4, dimension(:), allocatable :: r4d1temp
      integer, parameter :: idxcr=1,idycr=2

	error=0
      call errhandle(nf_open(trim(fname),0,ncid))
      call errhandle(nf_inq_dimlen(ncid,1,mm))
      call errhandle(nf_inq_dimlen(ncid,2,nn))
	write(*,*) 'grid size ',mm,nn

      allocate(xnst(mm),ynst(nn))
      nXYn=maxval( (/mm, nn/) )
      allocate(r4d1temp(nXYn))

      call errhandle(nf_inq_vartype(ncid,idycr,ynsttyp))
      call errhandle(nf_inq_vartype(ncid,idxcr,xnsttyp))

      if (xnsttyp .eq. NF_FLOAT) then
         call errhandle(nf_get_var_real(ncid,idxcr,r4d1temp))
         xnst(1:mm)=r4d1temp(1:mm)
      else if (xnsttyp .eq. NF_DOUBLE) then
         call errhandle(nf_get_var_double(ncid,idxcr,xnst))
      else 
      	error=-1
      	write(*,*) 'x-vector unexpected type',xnsttyp
		goto 99
      endif

      if (ynsttyp .eq. NF_FLOAT) then
         call errhandle(nf_get_var_real(ncid,idycr,r4d1temp))
         ynst(1:nn)=r4d1temp(1:nn)
      else if (ynsttyp .eq. NF_DOUBLE) then
         call errhandle(nf_get_var_double(ncid,idycr,ynst))
      else 
      	error=-1
      	write(*,*) 'y-vector unexpected type',ynsttyp
      	goto 99
      endif
	deallocate(r4d1temp)
      call errhandle(nf_close(ncid))
      
 99  	continue
      end subroutine ReadGrid
