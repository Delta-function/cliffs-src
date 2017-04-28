C  Cliffs is an open-source model for tsunami propagation and inundation computations       
C  in the framework of the shallow-water equations. Cliffs implements:                      
C  (1) VTCS-2 finite-difference scheme in an open ocean and on the open boundary            
C 	 as in (Titov and Synolakis, 1995, 1998);                                            
C  (2) dimensional splitting as in (Titov and Gonzalez, 1997; Titov and Synolakis, 1998);   
C  (3) reflection and inundation computations as in (Tolkova, 2014);                        
C  (4) data flows similar to curvilinear MOST and MOST-4 (Burwell and Tolkova, 2008).       
C                               REFERENCE:                                                              
C  E. Tolkova. Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting 
C  Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314                     
C                                                                                          
C  FEATURES: Cartezian or Geophysical (lon/lat) coordinates; 2D or 1D domains;              
C           grid nesting with one-way coupling; initial conditions or/and boundary forcing;
C           Open MP; NetCDF format of I/O                                                  
C
C           Copyright (c) 2014, Elena Tolkova                                              
C           under the terms of FreeBSD License                                             
C                                                                                          
C  Redistribution and use in source and binary forms, with or without                       
C  modification, are permitted provided that the following conditions are met:              
C                                                                                          
C   1. Redistributions of source code must retain the above copyright notice, this          
C     list of conditions and the following disclaimer.                                     
C   2. Redistributions in binary form must reproduce the above copyright notice,            
C     this list of conditions and the following disclaimer in the documentation            
C     and/or other materials provided with the distribution.                               
C                                                                                          
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND          
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED            
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                   
C  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR          
C  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES           
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;             
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND              
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT               
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS            
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                             

      program CLIFFS_MAIN
	use MASTER
	use NESTED
	use PARAMETERS
	
      implicit none
      include 'netcdf.inc'

      character*200 casetitle,ingrids(99),fnm,Bname,bprfx,aprfx
      character*200 indir,casename,outdir,srcdir,pathout,params
      character*200 notes, anotes
      integer ishot, bcnt(99), gcnt
      integer nc,mc,ii1,ii2,jj1,jj2
      integer j,istep,error,ns,mm,nn,i,irec,mx
      logical hinput,uinput,vinput,bndr_input,nestedgrids,arv(99),iarv  
      integer ncid(3),idtime(3),idvar(3),idmax(2),ncmaxid
      integer gncid,gidvar(3),gtim_id
      integer nxsub,nysub,nsid(4),feedid(4,99)
      integer date1time(8),date2time(8),dd(8)
      character*10 b(3)
	integer omp_get_num_threads, omp_get_thread_num,iargc
      real*8 t,cc,aa
      real*8, dimension(:), allocatable :: gglon,gglat
      
      irec=iargc()
      if (irec.LT.5) then
      write(*,*) 'Help: */Cliffs <output dir/prefix> <source dir/>
     &<boundary prefix, 0 for no input> <area prefix, 0 for no input>
     & <param dir/param file> [<short notes>]'
         goto 999
      endif  

      pathout=' '
      casetitle=' '
      outdir=' '
      call getarg(1,pathout)
      nc=index(pathout,'/',.true.)+1
      mc=index(pathout,' ')-1
      casetitle=pathout(nc:mc)
      outdir=pathout(1:(nc-1))
      
      srcdir=' '
      call getarg(2,srcdir)
      bprfx=' '
      call getarg(3,bprfx)
      aprfx=' '
      call getarg(4,aprfx)
      
      pathout=' '
      params=' '
      indir=' '
      call getarg(5,pathout)
      nc=index(pathout,'/',.true.)+1
      mc=index(pathout,' ')-1
      params=pathout(nc:mc)
      indir=pathout(1:(nc-1))
      
      notes=' '
      mc=0
      do i=6,irec
	      anotes=' '
      	call getarg(i,anotes)
      	nc=index(anotes,' ')
      	notes((mc+1):(mc+nc))=anotes(1:nc)
      	mc=mc+nc
      enddo
            
      casename=' '
      casename = trim(outdir)//trim(casetitle)//'_log.txt'
      write(*,*)'*',trim(casename),'*'
      open(unit=9,file=trim(casename),status='unknown',form='formatted')
      call date_and_time(b(1), b(2), b(3), date1time)
      write(9,15) date1time(1),date1time(2),date1time(3),
     &             date1time(5),date1time(6),date1time(7)
 15  	format (i4,'/',i2,'/',i2,i6,':',i2,':',i2)
      write(9,*) trim(notes)
!$OMP PARALLEL
	if (omp_get_thread_num().eq.0)
     &	write(9,*) 'Parallel threads available : ', 
     &		OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
      write(9,*) 'Case Title: ',trim(casetitle)
      write(9,*) 'Input Data: ',trim(bprfx),', ',trim(aprfx)
      
! Read parameter file 
      call allabout(indir,params,ingrids,error)
      if (error.NE.0) goto 99
	nestedgrids=.false.	
	if(nests.gt.0) nestedgrids=.true.
	arv=.false.   ! wave has not yet reached inner grids, if any  
	  
! Impose vertical wall at depth dwall, if no topography
	if (itopo.eq.0) then
		do j=1,nYn 
			do i=1,nXn 
				if(dep(i,j).le.dwall) dep(i,j)=-999
			enddo
		enddo
	endif

! Take care of 1D and 2D cases
	if (nXn.gt.2) then
		ii1=2
		ii2=nXn-1
		xrun=.true.
		nc=nYn
	else
		ii1=1
		ii2=1
		lonsub=1
		xrun=.false.
		nc=1
	end if
	if (nYn.gt.2) then
		jj1=2
		jj2=nYn-1
		yrun=.true.
		mc=nXn
	else
		jj1=1
		jj2=1
		latsub=1
		yrun=.false.
		mc=1
	end if
	mx=max(nc,mc) ! max length (nodes) of feeding boundary

! read input
	if((trim(bprfx).eq.'0').and.(trim(aprfx).eq.'0')) then
		write(9,*) 'no input of any kind'
		goto 99
		endif
      nTn=0   
      hinput=.false.
      uinput=.false.
      vinput=.false.
      bndr_input=.false.  
	if(trim(bprfx).ne.'0') then
		Bname=trim(bprfx)//'_'//trim(MBathyFile)
		call get_boundary_feed(Bname,srcdir,cuke,irec,error)
		bndr_input=.true.
 		endif
 	if (error.ne.0) goto 99
 	if(trim(aprfx).ne.'0') then
		call get_initial_state(aprfx,srcdir,irec,
     &	hinput,uinput,vinput,error)
      	endif
 	if (error.ne.0) goto 99
      if (nTn.LT.2) freeze=.false.
 	t=tinput(irec)   ! start time
     	
!  allocate ocean state vars
      allocate(cel(nXn,nYn),xvel(nXn,nYn),yvel(nXn,nYn))
      allocate(velmax(nXn,nYn),hmax(nXn,nYn))
      cel=dep
      xvel=0
      yvel=0
      hmax=0
      velmax=0
	
	if (uinput) call apply_initial_conditions(nXsrc,nYsrc,
     &	Xsrc,Ysrc,usrc,nXn,nYn,Xcrd,Ycrd,xvel)
	if (vinput) call apply_initial_conditions(nXsrc,nYsrc,
     &	Xsrc,Ysrc,vsrc,nXn,nYn,Xcrd,Ycrd,yvel)
     	if(hinput) then
		if(quake.eq.1) then ! deform bottom
			qsrc=-qsrc
			call apply_initial_conditions(nXsrc,nYsrc,Xsrc,Ysrc,qsrc,
     &		nXn,nYn,Xcrd,Ycrd,dep)
     		else ! deform surface
 			call apply_initial_conditions(nXsrc,nYsrc,Xsrc,Ysrc,qsrc,
     &		nXn,nYn,Xcrd,Ycrd,cel)
     			end if     			
      endif      
	if (hinput.or.uinput.or.vinput) then
      	deallocate(usrc,vsrc,qsrc)  
      	deallocate(Xsrc,Ysrc)  
      endif

	do j=1,nYn   ! switch from height h to celerity sqrt(gh)
		do i=1,nXn 
			if(cel(i,j).le.ground) then
				cel(i,j)=0
			else
				cel(i,j)=sqrt(grav*cel(i,j))
			endif
		enddo
	enddo
	celmin=sqrt(grav*ground)

! open netcdf files for boundary feed into children grids
	if (nests.GE.1) then
	fnm=' '
	do ns=1,nests
		mm=nXnst(ns)
		nn=nYnst(ns)
		fnm=trim(outdir)//trim(casetitle)//'_'//trim(ingrids(ns))
		call opennc4feed(fnm,mm,nn,nsid)
		feedid(1:4,ns)=nsid(1:4)
	enddo
	endif

      casename=' '
      casename = trim(outdir)//trim(casetitle)
      fnm = ' '
      fnm=trim(casename)//'_maxwave.nc'
      call opennc4max(fnm,ncmaxid,idmax,cartesian)
      if (seaout.le.steps_total) then		! open netcdf sea screenshot files
      	nxsub = (nXn-1)/lonsub+1
      	nysub = (nYn-1)/latsub+1
      	allocate(gglon(nxsub),gglat(nysub))
		gglon=Xcrd(1:nXn:lonsub)
		gglat=Ycrd(1:nYn:latsub)
      	fnm = ' '
      	fnm=trim(casename)//'_sea_h.nc'
      	call opennc4sea(gglon,gglat,fnm,1,ncid(1),
     &      idtime(1),idvar(1),nxsub,nysub,cartesian)
     		if(xrun) then
      		fnm= trim(casename)//'_sea_u.nc'
      		call opennc4sea(gglon,gglat,fnm,2,ncid(2),
     &      	idtime(2),idvar(2),nxsub,nysub,cartesian)
     		endif
     		if(yrun) then
      		fnm=trim(casename)//'_sea_v.nc'
      		call opennc4sea(gglon,gglat,fnm,3,ncid(3),
     &      	idtime(3),idvar(3),nxsub,nysub,cartesian)
      	endif
      	deallocate(gglon,gglat)
 	      ishot = 1
      	call output_sea_state(ncid,idtime,idvar,
     &     		nxsub,nysub,ishot,t)
     	end if
     
! open netcdf gages output file
      bcnt=0
      gcnt=0
	if(Ngages.ge.1) then
		allocate(gglon(Ngages),gglat(Ngages))
		forall(i=1:Ngages)
			gglon(i)=Xcrd(Igages(i))
			gglat(i)=Ycrd(Jgages(i))
		end forall	
	      call opennc4gages(casename,gncid,gtim_id,
     &	gidvar,Ngages,gglon,gglat,cartesian,xrun,yrun)
		deallocate(gglon,gglat)
	      call output_gages(gncid,gtim_id,gidvar,gcnt,t,Ngages)
      end if     
      deallocate(Xcrd,Ycrd)
      
	allocate(edge1(mx,3),edge2(mx,3))  ! to hold each-step bndr feed of state vars
	edge1=0
	edge2=0

! Compute boundary feed into nested grids at the start of computations - added on 2016/07/23
      if(nestedgrids)then
       	do ns = 1,nests
           	  nsid(1:4)=feedid(1:4,ns)
           	  iarv=arv(ns)
              call feed_children_grid(ns,nsid,bcnt(ns),t,ground,iarv)
              arv(ns)=iarv
           	end do
      endif
            
! Time Loop
      do istep = 1,steps_total
      t=t+dt 
      do while (bndr_input.AND.t.ge.tinput(irec)) ! set tinput(irec-1) <= t < tinput(irec)
          	if(irec.eq.nTn) then
          		bndr_input=.false.
			edge1=0
			edge2=0
			write(9,*) 'bndr feed ended at step', istep
			if(freeze) goto 95
		else
          		irec=irec+1
          	endif
      end do 

	if (xrun) then	! computations in x/lon direction
      	if(bndr_input) then
          		cc=(tinput(irec)-t)/(tinput(irec)-tinput(irec-1))		
      		edge1(1:nYn,1:3)=west(1:nYn,1:3,irec-1)*cc
     & 			+west(1:nYn,1:3,irec)*(1-cc)
      		edge2(1:nYn,1:3)=east(1:nYn,1:3,irec-1)*cc
     & 			+east(1:nYn,1:3,irec)*(1-cc)
      	endif          
!$OMP PARALLEL DO
		do j=jj1,jj2
           		call cliffs(1,j,nXn)
		end do
!$OMP END PARALLEL DO
	endif ! if nXn>2, xrun

	if (yrun) then	! computations in y/lat direction
      	if(bndr_input) then
          		cc=(tinput(irec)-t)/(tinput(irec)-tinput(irec-1))		
      		edge1(1:nXn,1:3)=south(1:nXn,1:3,irec-1)*cc
     & 			+south(1:nXn,1:3,irec)*(1-cc)
      		edge2(1:nXn,1:3)=north(1:nXn,1:3,irec-1)*cc
     & 			+north(1:nXn,1:3,irec)*(1-cc)
		endif            
!$OMP PARALLEL DO
	      do j=ii1,ii2
           		call cliffs(2,j,nYn)
      	end do
!$OMP END PARALLEL DO
	endif ! if nYn>2, yrun
	
! update max wave height
	do j=1,nYn 
		do i=1,nXn
			aa=cel(i,j)
			if(aa.gt.hmax(i,j)) hmax(i,j)=aa
			aa=xvel(i,j)**2+yvel(i,j)**2
			if(aa.gt.velmax(i,j)) velmax(i,j)=aa
		enddo
	enddo		

! Compute boundary feed into nested grids
      if(nestedgrids.and.mod(istep,bndout).EQ.0)then
       	do ns = 1,nests
           	  nsid(1:4)=feedid(1:4,ns)
           	  iarv=arv(ns)
              call feed_children_grid(ns,nsid,bcnt(ns),t,ground,iarv)
              arv(ns)=iarv
           	end do
      endif

! Output screenshot and gage time histories
      if(mod(istep,seaout).EQ.0)
     &  call output_sea_state(ncid,idtime,idvar,nxsub,nysub,ishot,t)
	
	if((Ngages.ge.1).and.(mod(istep,gout).EQ.0))then        
	      call output_gages(gncid,gtim_id,gidvar,gcnt,t,Ngages)
	end if ! Ngages>0
! Update max wave	
	if(mod(istep,maxout).EQ.0)
     &	call output_maxwave(ncmaxid,idmax)

      end do			! end of Time Loop
      
 95	call output_maxwave(ncmaxid,idmax)
! close netcdf files 
	if (Ngages.gt.0) call errhandle(nf_close(gncid))
	if (seaout.le.steps_total) then
      	call errhandle(nf_close(ncid(1)))
      	if(xrun) call errhandle(nf_close(ncid(2)))
      	if(yrun) call errhandle(nf_close(ncid(3)))
      end if
      call errhandle(nf_close(ncmaxid))
      call date_and_time(b(1), b(2), b(3), date2time)
      write(9,15) date2time(1),date2time(2),date2time(3),
     &             date2time(5),date2time(6),date2time(7)
     	dd=date2time-date1time
     	nc=3600*(24*dd(3)+dd(5))+60*dd(6)+dd(7)
     	mc=mod(nc,60)
     	nc=(nc-mc)/60
      write(9,*) 'Run finished in ',nc,'min ',mc,'sec'
     
 99	close(unit=9,status='keep')

999   stop
      end program Cliffs_main
