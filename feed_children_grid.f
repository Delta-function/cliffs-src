C  boundary feed into ns-th nested grid                                                     
C  for each edge node, uses pre-computed coordinates of a left-next parent cell *_n[x/y]    
C  and interpolation weights for the left-next parent row and column *_c[x/y]               
C                                                                                          
C  This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine feed_children_grid(ns,bn_id,irec,time,hmin,iarv) 
      use MASTER, only: dep,xvel,yvel,cel,grav,nXYn
      use NESTED
      use PARAMETERS, only: cuke
	implicit none

      real*8 uvq(maxlen,3),clm(nXYn,3),clm1(nXYn,3)
      real*8 cc,time,hmin
	integer ns,irec,bn_id(4)
      integer i,k,ii,kk,ncid,len,ii1,ii2,kk1,kk2
      logical iarv  ! wave arrival flag
      
      kk1=WE_ny(1,ns)
      kk2=WE_ny(nYnst(ns),ns)+1
      ii1=SN_nx(1,ns)
      ii2=SN_nx(nXnst(ns),ns)+1
	if (.not.iarv) then
	do i=ii1,ii2
		do k=kk1,kk2
			if(dep(i,k).gt.hmin) then
		if(abs(cel(i,k)**2/grav-dep(i,k)).gt.(0.1*cuke)) then
					iarv=.true.
					exit
				endif
		if((abs(xvel(i,k))+abs(yvel(i,k)))*cel(i,k).gt.cuke) then
					iarv=.true.
					exit
				endif
			endif
		enddo
		if(iarv) then
			write(9,*) 'reached child grid ',ns,' at ',time,'sec'
			call flush(9)
			exit
		endif
	enddo
	endif
			    
      if(iarv) then
      
      irec=irec+1      
	len=nYnst(ns)    ! West and East boundaries
      
      ii=W_nx(ns)		! West
      cc=W_cx(ns) 
      do kk=kk1,kk2
		if(dep(ii,kk).gt.hmin) then    
       		clm(kk,1)=cc*xvel(ii,kk)
       		clm(kk,2)=cc*yvel(ii,kk)
       		clm(kk,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm(kk,1:3)=0
       	endif
      enddo
      ii=ii+1
      cc=1-cc 
      do kk=kk1,kk2
		if(dep(ii,kk).gt.hmin) then    
       		clm1(kk,1)=cc*xvel(ii,kk)
       		clm1(kk,2)=cc*yvel(ii,kk)
       		clm1(kk,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm1(kk,1:3)=0
       	endif
      enddo
      clm(kk1:kk2,1:3)=clm(kk1:kk2,1:3)+clm1(kk1:kk2,1:3)
      do i=1,len
       	kk=WE_ny(i,ns)
       	cc=WE_cy(i,ns)
       	uvq(i,1:3)=clm(kk,1:3)*cc+clm(kk+1,1:3)*(1-cc)
      enddo
 	ncid=bn_id(1)	
	call output_feed(ncid,len,maxlen,uvq,irec,time)	

      ii=E_nx(ns)  	! East
      cc=E_cx(ns)
      do kk=kk1,kk2
		if(dep(ii,kk).gt.hmin) then    
       		clm(kk,1)=cc*xvel(ii,kk)
       		clm(kk,2)=cc*yvel(ii,kk)
       		clm(kk,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm(kk,1:3)=0
       	endif
      enddo
      ii=ii+1
      cc=1-cc 
      do kk=kk1,kk2
		if(dep(ii,kk).gt.hmin) then    
       		clm1(kk,1)=cc*xvel(ii,kk)
       		clm1(kk,2)=cc*yvel(ii,kk)
       		clm1(kk,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm1(kk,1:3)=0
       	endif
      enddo
      clm(kk1:kk2,1:3)=clm(kk1:kk2,1:3)+clm1(kk1:kk2,1:3)
      do i=1,len
       	kk=WE_ny(i,ns)
       	cc=WE_cy(i,ns)
       	uvq(i,1:3)=clm(kk,1:3)*cc+clm(kk+1,1:3)*(1-cc)
      enddo
 	ncid=bn_id(2)	
	call output_feed(ncid,len,maxlen,uvq,irec,time)	
 
	len=nXnst(ns)		! South and North boundaries
      
      kk=S_ny(ns)		! South
      cc=S_cy(ns)
      do ii=ii1,ii2
		if(dep(ii,kk).gt.hmin) then    
       		clm(ii,1)=cc*xvel(ii,kk)
       		clm(ii,2)=cc*yvel(ii,kk)
       		clm(ii,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm(ii,1:3)=0
       	endif
      enddo
      kk=kk+1
      cc=1-cc
      do ii=ii1,ii2
		if(dep(ii,kk).gt.hmin) then    
       		clm1(ii,1)=cc*xvel(ii,kk)
       		clm1(ii,2)=cc*yvel(ii,kk)
       		clm1(ii,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm1(ii,1:3)=0
       	endif
      enddo
      clm(ii1:ii2,1:3)=clm(ii1:ii2,1:3)+clm1(ii1:ii2,1:3)
      do i=1,len
       	ii=SN_nx(i,ns)
       	cc=SN_cx(i,ns)
       	uvq(i,1:3)=clm(ii,1:3)*cc+clm(ii+1,1:3)*(1-cc)
      enddo
 	ncid=bn_id(3)	
	call output_feed(ncid,len,maxlen,uvq,irec,time)	

      kk=N_ny(ns)		! North
      cc=N_cy(ns)
      do ii=ii1,ii2
		if(dep(ii,kk).gt.hmin) then    
       		clm(ii,1)=cc*xvel(ii,kk)
       		clm(ii,2)=cc*yvel(ii,kk)
       		clm(ii,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm(ii,1:3)=0
       	endif
      enddo
      kk=kk+1
      cc=1-cc
      do ii=ii1,ii2
		if(dep(ii,kk).gt.hmin) then    
       		clm1(ii,1)=cc*xvel(ii,kk)
       		clm1(ii,2)=cc*yvel(ii,kk)
       		clm1(ii,3)=cc*(cel(ii,kk)**2/grav-dep(ii,kk))
       	else
       		clm1(ii,1:3)=0
       	endif
      enddo
      clm(ii1:ii2,1:3)=clm(ii1:ii2,1:3)+clm1(ii1:ii2,1:3)
      do i=1,len
       	ii=SN_nx(i,ns)
       	cc=SN_cx(i,ns)
       	uvq(i,1:3)=clm(ii,1:3)*cc+clm(ii+1,1:3)*(1-cc)
      enddo
 	ncid=bn_id(4)	
	call output_feed(ncid,len,maxlen,uvq,irec,time)	
	endif ! if iarv
      end subroutine feed_children_grid

    
