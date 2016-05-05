!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! depth_ssl: node-wise changes to a depth map to limit depth variations.                   ! 
! Improves stability when using VTCS-2 scheme                                              !
!                                                                                          !
!           Copyright (c) 2014, Elena Tolkova                                              !
!           under the terms of FreeBSD License                                             !
!                                                                                          !
! Redistribution and use in source and binary forms, with or without                       !
! modification, are permitted provided that the following conditions are met:              !
!                                                                                          !
!  1. Redistributions of source code must retain the above copyright notice, this          !
!     list of conditions and the following disclaimer.                                     !
!  2. Redistributions in binary form must reproduce the above copyright notice,            !
!     this list of conditions and the following disclaimer in the documentation            !
!     and/or other materials provided with the distribution.                               !
!                                                                                          !
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND          !
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED            !
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                   !
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR          !
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES           !
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;             !
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND              !
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT               !
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS            !
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                             !
!                                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program depth_ssl  

      implicit none
      include 'netcdf.inc'

      character*200 fname,inname,outdir,outname,xname,yname,att
      character*10 query
      integer xid_dim,yid_dim,mx_dim,my_dim
      integer ncid,ncid_out,var1,var2,var3,var4
      integer xcrdtyp,ycrdtyp,topotyp,mtxtyp
      real*4, dimension(:,:), allocatable :: bb
      real*8, dimension(:,:), allocatable :: bbb
      real*4, dimension(:), allocatable :: xxx
      real*8, dimension(:), allocatable :: yyy
      integer inargc,nXn,nYn,nc1,nc2,n
      integer ndims,nvar,ngatts,nmore
      real*8 alpha, depmin, GM(3,3)
      integer, parameter :: idxcr=1,idycr=2,idbth=3,idmtx=4

      inargc = iargc()
      if (inargc.lt.3) then
         write(*,*) 'Help: ./depth_ssl <full path to input bathy file>
     &  <OutputDir/> <n (0-default params)> [OutputName]' 
         stop
      endif  

      inname=' '
      outname=' '
      outdir=' '
      fname=' '
      query=' '
      
      call getarg(1,inname)
      call getarg(2,outdir)
      call getarg(3,query)
      
 	if(query(1:1).eq.'0') then
 		alpha=2.0
 		depmin=0.1
 	else
		write(*,*) 'ssl (1.5 - 2.5)? '
		read(*,*) alpha
		if (alpha.lt.1.5) alpha=1.5
		if (alpha.gt.2.5) alpha=2.5
		write(*,*) 'minimal flow depth (m)? '
		read(*,*) depmin
		if (depmin.le.0) then
			write(*,*) 'I do not think so'
			stop
		endif
	endif

      if(inargc.gt.3) then
	      call getarg(4,fname)
	      outname=trim(outdir)//trim(fname)
	else
	      nc1=index(inname,'/',.true.)+1
	      nc2=index(inname,'.nc')-1
	      outname=trim(outdir)//inname(nc1:nc2)//'_ssl.nc'
	endif
	write(*,*) 'processing ',trim(inname),' with'
	write(*,*) 'ssl=',alpha,' min depth=',depmin
	write(*,*) 'output to ', trim(outname)
      
      xname=' '
      yname=' '
      call errhandle(nf_open(trim(inname),0,ncid))
      call errhandle(nf_inq_dim(ncid,1,xname,nXn))
      call errhandle(nf_inq_dim(ncid,2,yname,nYn)) 
! define output      	
	call errhandle(nf_create(trim(outname), NF_CLOBBER, ncid_out))
! define dimensions
      call errhandle(nf_def_dim(ncid_out,trim(xname),nXn,xid_dim))
      call errhandle(nf_def_dim(ncid_out,trim(yname),nYn,yid_dim))
	call errhandle(nf_inq(ncid,ndims,nvar,ngatts,nmore))
	if(ndims.gt.2) then 
      	call errhandle(nf_def_dim(ncid_out,'Mrow',3,mx_dim))
      	call errhandle(nf_def_dim(ncid_out,'Mclm',3,my_dim))
	endif
! copy global attributes
	do n=1,ngatts
		call errhandle(nf_inq_attname(ncid,NF_GLOBAL,n,att))	
		call errhandle(nf_copy_att(ncid,NF_GLOBAL,
     &		trim(att),ncid_out,NF_GLOBAL))
	enddo
! define variables and write atts	
	call copy_vardefine(ncid,ncid_out,idxcr,var1)
	call copy_vardefine(ncid,ncid_out,idycr,var2)
	call copy_vardefine(ncid,ncid_out,idbth,var3)
	if(ndims.gt.2) then 
      	call copy_vardefine(ncid,ncid_out,idmtx,var4)
	endif
! leave define mode
      call errhandle(nf_enddef(ncid_out))
	
! transfer data	
	write(*,*) 'grid size',nXn,nYn
	
      call errhandle(nf_inq_vartype(ncid, idycr, ycrdtyp))
      call errhandle(nf_inq_vartype(ncid, idxcr, xcrdtyp))
      call errhandle(nf_inq_vartype(ncid, idbth, topotyp))

      if (xcrdtyp .eq. NF_FLOAT) then
      	allocate(xxx(nXn))
         	call errhandle(nf_get_var_real(ncid,idxcr,xxx))
		call errhandle(nf_put_var_real(ncid_out,var1,xxx))
		deallocate(xxx)
      else if (xcrdtyp .eq. NF_DOUBLE) then
      	allocate(yyy(nXn))
         	call errhandle(nf_get_var_double(ncid,idxcr,yyy))
		call errhandle(nf_put_var_double(ncid_out,var1,yyy))
		deallocate(yyy)
      else 
      	write(*,*) 'wrong lon/xxx type:', xcrdtyp
		goto 100
      endif

      if (ycrdtyp .eq. NF_FLOAT) then
      	allocate(xxx(nYn))
         	call errhandle(nf_get_var_real(ncid,idycr,xxx))
		call errhandle(nf_put_var_real(ncid_out,var2,xxx))
		deallocate(xxx)
      else if (ycrdtyp .eq. NF_DOUBLE) then
      	allocate(yyy(nYn))
         	call errhandle(nf_get_var_double(ncid,idycr,yyy))
		call errhandle(nf_put_var_double(ncid_out,var2,yyy))
		deallocate(yyy)
      else 
      	write(*,*) 'wrong lat/yyy type:', ycrdtyp
		goto 100
      endif

      if (topotyp .eq. NF_FLOAT) then
      	allocate(bb(nXn,nYn),bbb(nXn,nYn))
         	call errhandle(nf_get_var_real(ncid,idbth,bb))
         	bbb=bb
         	call setSSLim(nXn,nYn,bbb,alpha,depmin)   ! adjust depth map
         	bb=bbb
		call errhandle(nf_put_var_real(ncid_out,var3,bb))
		deallocate(bb,bbb)
      else if (topotyp .eq. NF_DOUBLE) then
      	allocate(bbb(nXn,nYn))
         	call errhandle(nf_get_var_double(ncid,idbth,bbb))
         	call setSSLim(nXn,nYn,bbb,alpha,depmin)   ! adjust depth map         	
		call errhandle(nf_put_var_double(ncid_out,var3,bbb))
		deallocate(bbb)
      else 
      	write(*,*) 'wrong bathy data type:', topotyp
		goto 100
      endif

	if(ndims.gt.2) then
	      call errhandle(nf_inq_vartype(ncid, idmtx, mtxtyp))	
		if (mtxtyp .eq. NF_DOUBLE) then
         		call errhandle(nf_get_var_double(ncid,idmtx,GM))
			call errhandle(nf_put_var_double(ncid_out,var4,GM))
      	else 
      		write(*,*) 'wrong O2G matrix data type:', mtxtyp
			goto 100
      	endif
      endif

 100  call errhandle(nf_close(ncid))
      call errhandle(nf_close(ncid_out))
      
      stop
      end
      
      
!********************************************************
	subroutine copy_vardefine(nc1,nc2,var1,var2)
      implicit none
      include 'netcdf.inc'
	
	integer nc1,nc2,var1,var2 !nc ids of input and output datasets and vars
	integer n,vtype,ndims,natts,dimids(10)
	character*50 vname
	      
      call errhandle(
     & nf_inq_var(nc1,var1,vname,vtype,ndims,dimids,natts))
      
      call errhandle(
     & nf_def_var(nc2,trim(vname),vtype,ndims,dimids(1:ndims),var2))
      
      do n=1,natts
      	vname=' '
      	call errhandle(nf_inq_attname(nc1,var1,n,vname))
      	call errhandle(nf_copy_att(nc1,var1,trim(vname),nc2,var2))
      enddo
      end subroutine copy_vardefine
!********************************************************

	subroutine setSSLim(nXn,nYn,bbb,alpha,depmin)
	implicit none
	
	integer nXn,nYn
	real*8 alpha,depmin,bbb(nXn,nYn),c1,c2,alp2,dep
      real*8, dimension(:), allocatable :: q,p
      logical, dimension(:), allocatable :: land
      integer nodes,n,i,j,itry,nm
	
	nm=max(nXn,nYn)
	allocate(q(nm),p(nm))
	allocate(land(nm))

	write(*,*) 'Node-wise depth adjustments in open water:'
	do itry=1,12
	
	nodes=0
	do i=1,nXn
		q(1:nYn)=bbb(i,1:nYn)
		p(1:nYn)=q(1:nYn)
		land(1:nm)=.false.
		do j=1,nYn
			if (q(j).lt.depmin) then
				land(j)=.true.
			else
				q(j)=sqrt(q(j))
			endif
		enddo
		n=nYn-1
		do j=2,n
			if(.not.(land(j-1).or.land(j).or.land(j+1))) then
				c1=q(j+1)/q(j)
				c2=q(j-1)/q(j)
				if ((c1.gt.alpha).or.(c2.gt.alpha)) then
					q(j)=(q(j-1)+q(j+1))/2
					p(j)=q(j)*q(j)
					nodes=nodes+1
				endif
			endif
		enddo
		bbb(i,1:nYn)=p(1:nYn)
	enddo
	
	do i=1,nYn
		q(1:nXn)=bbb(1:nXn,i)
		p(1:nXn)=q(1:nXn)
		land(1:nm)=.false.
		do j=1,nXn
			if (q(j).lt.depmin) then
				land(j)=.true.
			else
				q(j)=sqrt(q(j))
			endif
		enddo
		n=nXn-1
		do j=2,n
			if(.not.(land(j-1).or.land(j).or.land(j+1))) then
				c1=q(j+1)/q(j)
				c2=q(j-1)/q(j)
				if ((c1.gt.alpha).or.(c2.gt.alpha)) then
					q(j)=(q(j-1)+q(j+1))/2
					p(j)=q(j)*q(j)
					nodes=nodes+1
				endif
			endif
		enddo
		bbb(1:nXn,i)=p(1:nXn)
	enddo
	write(*,*) 'passage',itry,':',nodes,' nodes adjusted'
	if(nodes.eq.0) exit
	enddo  ! end itry loop
	
	if(nodes.ne.0) 
     &	write(*,*) 'ATTN: processing is not complete', 
     & ' - too rough bathy or too small alpha'
		
	write(*,*) 'Adjustments at coastlines:'
	alp2=alpha**2
	do itry=1,12
	
	nodes=0
	do i=1,nXn
		q(1:nYn)=bbb(i,1:nYn)
		land(1:nm)=.false.
		do j=1,nYn
			if (q(j).lt.depmin) land(j)=.true.
		enddo
		n=nYn-1
		do j=2,n
			if(.not.land(j)) then
				if(land(j+1).and.(.not.land(j-1))
     &				.or.(.not.land(j+1)).and.land(j-1)) then
					if(land(j+1)) then
						dep=q(j-1)
					else
						dep=q(j+1)
					endif
					if( dep/q(j).gt.alp2 ) then
						q(j)=dep/alp2	
						nodes=nodes+1
					endif
				endif
			endif
		enddo
		bbb(i,1:nYn)=q(1:nYn)			
	enddo		
	do i=1,nYn
		q(1:nXn)=bbb(1:nXn,i)
		land(1:nm)=.false.
		do j=1,nXn
			if (q(j).lt.depmin) land(j)=.true.
		enddo
		n=nXn-1
		do j=2,n
			if(.not.land(j)) then
				if(land(j+1).and.(.not.land(j-1))
     &				.or.(.not.land(j+1)).and.land(j-1)) then
					if(land(j+1)) then
						dep=q(j-1)
					else
						dep=q(j+1)
					endif
					if( dep/q(j).gt.alp2 ) then
						q(j)=dep/alp2	
						nodes=nodes+1
					endif
				endif
			endif
		enddo
		bbb(1:nXn,i)=q(1:nXn)			
	enddo
	write(*,*) 'passage',itry,':',nodes,' nodes adjusted'
			
	if(nodes.eq.0) exit
	enddo  ! end itry loop
	
	if(nodes.ne.0) 
     &	write(*,*) 'ATTN: processing is not complete', 
     & ' - too rough bathy or too small alpha'
			
	end subroutine setSSLim
	