C  read lattices for enclosed (children) grids,                                             
C  compute factors needed for interpolating onto child boundaries                          
C
C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine NestedGrid(indir,ingrids,error)
      use MASTER, only: nXn,nYn,Xcrd,Ycrd
      use NESTED
      implicit none
      include 'netcdf.inc'

      character*200 fname,indir,ingrids(99)
      integer mm,nn,ns,i 
      integer indx,error,xcrdtyp,ycrdtyp	    
      integer ncid(99),idxcr,idycr
	real*8 xx,yy,x,y
      real*4, dimension(:), allocatable :: r4d1temp	
	real*8, dimension(:), allocatable :: r8d1temp
	real*8, dimension(:,:), allocatable :: xnst, ynst ! (maxXns,nests),(maxYns,nests)
	integer maxXns,maxYns
	
	nXnst=0
	nYnst=0
	ingrids=' '
      do ns = 1,nests  ! read for grid sizes
        	fname=' '
        	read(5,2) fname
        	indx=index(fname,'.nc')-1
        	ingrids(ns)=fname(1:indx)
        	fname = trim(indir)//trim(fname)

        	call errhandle(nf_open(trim(fname),0,ncid(ns)))
        	call errhandle(nf_inq_dimlen(ncid(ns),1,mm))
        	call errhandle(nf_inq_dimlen(ncid(ns),2,nn))
		nXnst(ns)=mm
		nYnst(ns)=nn
        	write(9,*) 'nested Grid',ns,trim(fname), mm,nn
      enddo
	maxXns=maxval(nXnst)
      maxYns=maxval(nYnst)
      maxlen=maxval( (/maxXns, maxYns /) )
	
	allocate(xnst(maxXns,nests),ynst(maxYns,nests))
      allocate(r4d1temp(maxlen),r8d1temp(maxlen))
      
	idxcr=1
      idycr=2	
      do ns = 1,nests
         	mm=nXnst(ns)
         	nn=nYnst(ns)

      	call errhandle(nf_inq_vartype(ncid(ns), idycr, ycrdtyp))
      	call errhandle(nf_inq_vartype(ncid(ns), idxcr, xcrdtyp))
      
      	if (xcrdtyp .eq. NF_FLOAT) then
         	call errhandle(nf_get_var_real(ncid(ns),idxcr,r4d1temp))
		      xnst(1:mm,ns)=r4d1temp(1:mm)
      	else if (xcrdtyp .eq. NF_DOUBLE) then
         	call errhandle(nf_get_var_double(ncid(ns),idxcr,r8d1temp))
         		xnst(1:mm,ns)=r8d1temp(1:mm)
      	else 
      		error=1
      		goto 99
      	endif
         
      	if (ycrdtyp .eq. NF_FLOAT) then
         	call errhandle(nf_get_var_real(ncid(ns),idycr,r4d1temp))
         		ynst(1:nn,ns)=r4d1temp(1:nn)
      	else if (ycrdtyp .eq. NF_DOUBLE) then
         	call errhandle(nf_get_var_double(ncid(ns),idycr,r8d1temp))
         		ynst(1:nn,ns)=r8d1temp(1:nn)
      	else 
      		error=1
      		goto 99
      	endif
     
         	call errhandle(nf_close(ncid(ns)))
      end do
      deallocate(r4d1temp,r8d1temp)
      
	do ns=1,nests ! check nested configuration
		x=xnst(1,ns)
		xx=xnst(nXnst(ns),ns)
		y=ynst(1,ns)
		yy=ynst(nYnst(ns),ns)
		write(9,*) 'Grid ',ns,' lon/lat limits: ',x,xx,y,yy		
      	if((x.lt.Xcrd(1)).or.(xx.gt.Xcrd(nXn)).or.
     & 		(y.lt.Ycrd(1)).or.(yy.gt.Ycrd(nYn))) then
      		error=-1
      		write(9,*) 'not within Master grid limits'
      		goto 99
      	endif
	enddo ! ns
	
! keep left-next parent cell and interpolation coefficient for each child edge node
      allocate(WE_ny(maxYns,nests),W_nx(nests),E_nx(nests))
      allocate(SN_nx(maxXns,nests),S_ny(nests),N_ny(nests))
      allocate(WE_cy(maxYns,nests),W_cx(nests),E_cx(nests))
      allocate(SN_cx(maxXns,nests),S_cy(nests),N_cy(nests))
      
	do ns=1,nests
		mm=nXnst(ns)
		nn=nYnst(ns)
	! North and South
		call pinpoint(Ycrd,nYn,ynst(1,ns),S_ny(ns),S_cy(ns))
		call pinpoint(Ycrd,nYn,ynst(nn,ns),N_ny(ns),N_cy(ns))		
		do i=1,mm
     	   call pinpoint(Xcrd,nXn,xnst(i,ns),SN_nx(i,ns),SN_cx(i,ns))
     		enddo
	! East and West
		call pinpoint(Xcrd,nXn,xnst(1,ns),W_nx(ns),W_cx(ns))
		call pinpoint(Xcrd,nXn,xnst(mm,ns),E_nx(ns),E_cx(ns))		
		do i=1,nn
     	   call pinpoint(Ycrd,nYn,ynst(i,ns),WE_ny(i,ns),WE_cy(i,ns))
     		enddo
	enddo

	deallocate(xnst,ynst)
 2    format(a)

 99  	continue
      end subroutine NestedGrid
      
! .....................................

      subroutine pinpoint(xx,nx,xp,ip,cp)
      integer nx,ip,ipp,med
      real*8 xx(nx),xp,cp

      ip=1
      ipp=nx
      do while (ip+2.le.ipp)
         med=(ip+ipp)/2
         if (xp.gt.xx(med)) then
            ip=med
         else
            ipp=med
         endif
      enddo
      if (nx.gt.1) then
		cp=(xx(ip+1)-xp)/(xx(ip+1)-xx(ip))
	else
		cp=0.5
	endif
	end subroutine pinpoint
