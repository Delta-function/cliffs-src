C  read boundary input into Master grid, and find start time (wave above threshold)         
C                                                                                         
C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine get_boundary_feed(srcname,indir,cuke,it,error)

      use MASTER, only: 
     & nXn,nYn,nTn,tinput,west,east,south,north,xrun,yrun
      implicit none
      include 'netcdf.inc'

      character*200 fname,fname2,srcname,indir
	integer nT,nY,nX
	real*8 cuke,t
	real*8 w1,w2,w3,w4,wave
      integer error,it,status
      integer ncid, ncerr

      error=0
      fname=' '
      fname=trim(indir)//trim(srcname)
	write(9,*) ' Boundary Input:'
	if (xrun) then

      fname2=' '
      fname2=trim(fname)//'_west.nc'
      write(9,*) trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncid)
      if (ncerr.ne.NF_NOERR) then
      	error=-1
      	write(9,*) 'not found'
      	goto 99
      	endif

      call errhandle(nf_inq_dimlen(ncid,1,nY))
      call errhandle(nf_inq_dimlen(ncid,3,nTn))
      if (nY.ne.nYn) then
      	error=-1
      	write(9,*) 'size mismatch: expected ',nYn,' provided ',nY
      	goto 99
      endif
      if (nTn.lt.2) then
      	error=-1
      	write(9,*) 'unsufficient time-series length: ',nTn,' samples'
      	goto 99
      endif
      
	allocate(west(nYn,3,nTn),east(nYn,3,nTn),stat=status) 
	if(status.gt.0) then
		error=-1
		write(9,*) 'memory allocation unsuccessful'
		goto 99
		endif
	allocate(tinput(nTn)) 
      	
      call errhandle(nf_get_var_double(ncid,1,west))
      call errhandle(nf_get_var_double(ncid,2,tinput))
      call errhandle(nf_close(ncid))
      
      fname2=trim(fname)//'_east.nc'
      write(9,*) trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncid)
      if (ncerr .ne. NF_NOERR) then
      	error=-1
      	write(9,*) 'not found'
      	goto 99
      	endif
      
      call errhandle(nf_inq_dimlen(ncid,1,nY))
      call errhandle(nf_inq_dimlen(ncid,3,nT))
      if ((nY .ne. nYn).OR.(nT.NE.nTn)) then
      	error=-1
      write(9,*) 'size mismatch: expected ',nYn,nTn,' provided ',nY,nT
      	goto 99
      endif
      
      call errhandle(nf_get_var_double(ncid, 1, east))
      call errhandle(nf_close(ncid))
      write(9,*) 'x-bndr input ', nY,'x',nT,' read successfully'  
      endif  ! nXn>2
      
      if(yrun) then
      
      fname2=trim(fname)//'_south.nc'
      write(9,*) trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncid)
      if (ncerr .ne. NF_NOERR) then
      	error=-1
      	write(9,*) 'not found'
      	goto 99
      	endif
      
      ncerr=nf_inq_dimlen(ncid,1,nX)
      call errhandle(ncerr)
      ncerr=nf_inq_dimlen(ncid,3,nT)
      call errhandle(ncerr)
      
      if(.not.xrun) then  ! get time
      	nTn=nT
      	allocate(tinput(nTn))
      	call errhandle(nf_get_var_double(ncid,2,tinput))
      endif
      if ((nX.ne.nXn).OR.(nT.ne.nTn)) then
		error=-1
      write(9,*) 'size mismatch: expected ',nXn,nTn,' provided ',nX,nT
      	goto 99
      endif
      
	allocate(south(nXn,3,nTn),north(nXn,3,nTn),stat=status)   
	if(status.gt.0) then
		error=-1
		write(9,*) 'memory allocation unsuccessful'
		goto 99
		endif
      
      ncerr = nf_get_var_double(ncid,1,south)
      call errhandle(nf_close(ncid))
      
      fname2=trim(fname)//'_north.nc'
      write(9,*) trim(fname2)
      ncerr = nf_open(trim(fname2),0,ncid)  
      if (ncerr.ne.NF_NOERR) then
      	error=-1
      	write(9,*) 'not found'
      	goto 99
      	endif
         
      call errhandle(nf_inq_dimlen(ncid,1,nX))
      call errhandle(nf_inq_dimlen(ncid,3,nT))
      if ((nX.NE.nXn).OR.(nT.NE.nTn)) then
		error=-1
      write(9,*) 'size mismatch: expected ',nXn,nTn,' provided ',nX,nT
      	goto 99
      endif      
      call errhandle(nf_get_var_double(ncid,1,north))
      call errhandle(nf_close(ncid))
      
      write(9,*) 'y-bndr input ', nX,'x',nT,' read successfully'  
      endif  ! nYn>2
          
	w1=0
	w2=0
	w3=0
	w4=0
	do it=1,nTn
		if(xrun) then
		w1=max(maxval(west(1:nYn,3,it)),-minval(west(1:nYn,3,it)))
		w2=max(maxval(east(1:nYn,3,it)),-minval(east(1:nYn,3,it)))
		endif
		if(yrun) then
		w3=max(maxval(south(1:nXn,3,it)),-minval(south(1:nXn,3,it)))
		w4=max(maxval(north(1:nXn,3,it)),-minval(north(1:nXn,3,it)))
		endif
		wave=max(w1,w2,w3,w4)
		if (wave.ge.cuke) exit
	end do
	if(it.ge.nTn) then
		error=-1
		write(9,*) 'bndr signal below threshold'
	else				
		write(9,*) 'feed duration (s):',tinput(1),' to',tinput(nTn)
		write(9,*) 'exceeds threshold at step',it,
     &   	',',tinput(it),'s, with ',wave,' m'
     	endif
      
 99   call flush(9)
      end subroutine get_boundary_feed
