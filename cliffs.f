!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! Cliffs SOLVER to propagate tsunami wave through one time step in 1D with:                !
! (1) VTCS-2 finite-difference scheme in an open ocean and on the open boundary            !
!	 as in (Titov and Synolakis, 1995, 1998), modified on Dec 25, 2014                                             !
! (2) reflection and inundation computations as in (Tolkova, 2014)                         !
!                               REFERENCES:                                                !   
! V. Titov and C. Synolakis. Modeling of Breaking and Nonbreaking Long-Wave Evolution and  !
! Runup Using VTCS-2. J. of Waterway, Port, Coastal, and Ocean Eng. Vol. 121, No 6 (1995), ! 
! pp. 308-316.                                                                             !
! E. Tolkova. Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting.!
! Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314                     !
!                                                                                          !
!     Copyright (C) 2014, Elena Tolkova                                                    !
!     For conditions of distribution and use, see copyright notice in cliffs_main.f        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine cliffs(iu,irw,n)

      use MASTER, only: dep,cel,xvel,yvel,s1,s2,zeta,grav,edge1,edge2
      use PARAMETERS, only: ground,crough,dt,celmin
      implicit none

      integer iu,irw,n
      real*8 h(n),pp(n),qq(n),vv(n),dx(n-1),scl(n)
      real*8 depth(n),Pinv(n),Qinv(n),v(n),u(n)
      real*8 pqvdL(4,99), pqvdR(4,99)
      real*8 uj,vj,pj,qj,dpj,dpm,dtx,hj
      real*8 ej,e2,e1,Qpj,Qjm,Qsum
      real*8 v1,v2,d1,d2,p1,p2,q1,q2,dj      
      real*8 cc,sphr,fm,fp,flood     
      integer i,j,i1,i2,kseg,k,k1,k2
      integer lghost(n),rghost(n)
      logical land(n),water,newland(n)
      real*8, parameter :: sqr2=1.41421356
      
	if(iu.eq.1) then
		forall(i=1:n)
	      	depth(i)=dep(i,irw)
	      	h(i)=cel(i,irw)
      		u(i)=xvel(i,irw)
      		v(i)=yvel(i,irw)
			scl(i)=0
		end forall
		dx=s1(1:(n-1),irw)
     	else
		forall(i=1:n)
            	depth(i)=dep(irw,i)
            	h(i)=cel(irw,i)
            	u(i)=yvel(irw,i)
            	v(i)=xvel(irw,i)
			scl(i)=zeta(i)
		end forall
		dx=s2
	endif  
	
      lghost(1:n)=0
      rghost(1:n)=0
! Find dry areas
      land=.false.
      do j=1,n 
      	if(h(j).le.celmin) land(j)=.true.
      end do
      newland=land
     
      water=.false.
      kseg=0
      if(.not.land(1)) then
            water=.true.
            kseg=1
            lghost(1)=1		
      endif
      k1=kseg+1 ! first wet segment with left shore
       
      do i=2,n
            if(water) then
                  if(land(i)) then
                  	flood=(h(i-1)**2)/grav-depth(i-1)+depth(i)
      if((flood.gt.ground).and.(.not.land(i-1)).and.(u(i-1).ge.0)) then ! expand
                              u(i)=u(i-1)
                              v(i)=v(i-1)
                              h(i)=celmin
                              newland(i)=.false.
                        else
                              water=.false. ! close wet segment
                              rghost(kseg)=i
                        endif
                  endif
            else
                  if(.not.land(i)) then ! start wet segment
                      	water=.true.
                      	kseg=kseg+1
                   	flood=(h(i)**2)/grav-depth(i)+depth(i-1) 
                   	if((flood.gt.ground).and.(u(i).le.0)) then ! expand
                   		lghost(kseg)=max(i-2,1)
                   		newland(i-1)=.false.
                   		u(i-1)=u(i)
                   		v(i-1)=v(i)
                   		h(i-1)=celmin
                   	else
                        	lghost(kseg)=i-1
                        endif
               		if(kseg.gt.1) then
               			if(lghost(kseg).lt.rghost(kseg-1)) then ! glue to previous wet segment 
                        	lghost(kseg)=0
                        	kseg=kseg-1
                        	rghost(kseg)=0
                        	end if
                  	endif
            	endif
            endif
      enddo
      
      if(kseg.eq.0) goto 99
      if (rghost(kseg).eq.0) rghost(kseg)=n
      k2=kseg-1
      if(newland(rghost(kseg))) k2=kseg
      
      where(newland)
      	h=0
      	u=0
      	v=0
      end where
      vv=v 

!   compute Riemann invarients
      Pinv=0
      Qinv=0
!    - in open sea	
      do k=1,kseg
            i1=lghost(k)+1
            i2=rghost(k)-1
            forall (i=i1:i2)
               	Pinv(i)=u(i)+2*h(i)
               	Qinv(i)=u(i)-2*h(i)
            end forall
      enddo	
!    - on left shoreline in ghost node
      if (k1.le.kseg) then 
      do k=k1,kseg
            i=lghost(k)
            j=i+1
            pqvdL(1:4,k)=(/-Qinv(j), -Pinv(j), v(j), depth(j)/)         
      enddo
      endif !if (k1.le.kseg)
!    - on right shoreline in ghost node
      if (k2.ge.1) then
      do k=1,k2
            i=rghost(k)
            j=i-1
            pqvdR(1:4,k)=(/-Qinv(j), -Pinv(j), v(j), depth(j)/)         
      enddo
      endif !if (k2.ge.1)	     
!    - on left wet edge
	if(.not.newland(1)) then
            Pinv(1)=u(1)+2*h(1)
            Qinv(1)=u(1)-2*h(1)
        pqvdL(1:4,1)=(/Pinv(1), Qinv(1), v(1), depth(1)/)
	end if
!    - on right wet edge
	if(.not.newland(n)) then
            Pinv(n)=u(n)+2*h(n)
            Qinv(n)=u(n)-2*h(n)
        pqvdR(1:4,kseg)=(/Pinv(n), Qinv(n), v(n), depth(n)/)
	end if	      
      qq=Qinv
      pp=Pinv
                  
! time-step integration      
	do 9 k=1,kseg
 
		i1=lghost(k)+1
      	i2=rghost(k)-1
        	if(i1.gt.i2) goto 9
        
		do i=i1,i2
        		j=i-1
        		if(i==i1) then
           			p1=pqvdL(1,k)
           			q1=pqvdL(2,k)
           			v1=pqvdL(3,k)
           			d1=pqvdL(4,k)
        		else
           			p1=Pinv(j)
           			q1=Qinv(j)
           			v1=v(j)
           			d1=depth(j)
        		end if
        		j=i+1;
        		if(i==i2) then
           			p2=pqvdR(1,k)
           			q2=pqvdR(2,k)
           			v2=pqvdR(3,k)
           			d2=pqvdR(4,k)
        		else
           			p2=Pinv(j)
           			q2=Qinv(j)
           			v2=v(j)
           			d2=depth(j)
        		end if
        		uj=u(i)
        		vj=v(i)
        		hj=h(i)**2/grav
!  Almost Manning - Preferred friction model     
       		cc=grav*crough*uj*sqrt((uj**2+vj**2)/hj)/hj  
!   Drag force 
      !  		cc=grav*crough*uj*sqrt(uj**2+vj**2)/hj
!  Manning friction - SLOW
      !  		cc=grav*crough*uj*sqrt(uj**2+vj**2)/hj**1.3333
!  Linear friction 
      !  		cc=grav*crough*uj/hj

        		pj=Pinv(i)
        		qj=Qinv(i)
	  		sphr=(pj-qj)*(pj+qj)*scl(i)
        		fm=dx(i-1)
        		fp=dx(i)
        		dtx=dt/(fm+fp)
        		dj=depth(i)
        		dpj=grav*(d2-dj)
        		dpm=grav*(dj-d1)        
	  	water=(dj.gt.ground).and.(d1.gt.ground).and.(d2.gt.ground)
	  		
        		ej=3*pj+qj
        		e2=3*p2+q2
        		e1=3*p1+q1
        		Qpj=(0.125*(e2+ej)*(p2-pj)-dpj)/fp 
        		Qjm=(0.125*(ej+e1)*(pj-p1)-dpm)/fm
        		if((e1.lt.0).and.(e2.gt.0).and.water) then
       Qsum=(0.25*ej*(p2-p1)-2*grav*(sqrt(dj*d2)-sqrt(dj*d1)))/(fm+fp)
        		else
        			Qsum=(Qpj+Qjm)/2
        		endif
        		Qsum=Qsum+0.25*ej*dtx*(Qjm-Qpj)
        		pp(i)=pj-dt*(Qsum+cc-sphr)       		
        
        		ej=3*qj+pj
        		e2=3*q2+p2
        		e1=3*q1+p1
        		Qpj=(0.125*(e2+ej)*(q2-qj)-dpj)/fp
        		Qjm=(0.125*(ej+e1)*(qj-q1)-dpm)/fm
        		if((e1.lt.0).and.(e2.gt.0).and.water) then
       Qsum=(0.25*ej*(q2-q1)-2*grav*(sqrt(dj*d2)-sqrt(dj*d1)))/(fm+fp)
        		else
         			Qsum=(Qpj+Qjm)/2
       		endif
        		Qsum=Qsum+0.25*ej*dtx*(Qjm-Qpj)
        		qq(i)=qj-dt*(Qsum+cc+sphr)       		
        
        		Qsum=(v2-vj)/fp-(vj-v1)/fm
        		vv(i)=vj-uj*dtx*(v2-v1-Qsum*(uj*dt+fp-fm))
        
		end do
9  	continue

      do k=1,kseg
            i1=lghost(k)+1
            i2=rghost(k)-1
            forall (i=i1:i2) 
            	h(i)=(pp(i)-qq(i))/4
           		u(i)=(pp(i)+qq(i))/2
           	end forall
      end do
      
!   check water surface angle on left shoreline
      if (k1.le.kseg) then 
      do k=k1,kseg
            i=lghost(k)+1
            j=i+1
            if(land(i)) then ! if newly included node
      		if(sqr2*h(i).gt.h(j)) h(i)=h(j)/sqr2
      		if(u(i).lt.u(j)) u(i)=u(j)
		 endif
      enddo
      endif !if (k1.le.kseg)
!   check water surface angle on right shoreline
      if (k2.ge.1) then
      do k=1,k2
            i=rghost(k)-1
            j=i-1
            if(land(i)) then ! if newly included node
      		if(sqr2*h(i).gt.h(j)) h(i)=h(j)/sqr2  
      		if(u(i).gt.u(j)) u(i)=u(j)
            endif
      enddo
      endif !if (k2.ge.1)

! Integration in edges
	flood=edge1(irw,3)+depth(1)
	if((.not.newland(1)).and.(.not.newland(2)) 
     &	.and.(flood.gt.ground)) then
        	pj=Pinv(1)
        	qj=Qinv(1)
        	p2=Pinv(2)
        	q2=Qinv(2)
	  	sphr=(pj-qj)*(pj+qj)*scl(1)
        	Qpj=(3*(q2+qj)+p2+pj)*(q2-qj)/8-grav*(depth(2)-depth(1))
        	qq(1)=qj-dt*(Qpj/dx(1)+sphr)
        	
        	pp(1)=edge1(irw,iu)+2*sqrt(grav*flood)

           	if(u(1).lt.0) then
               vv(1)=v(1)-dt*(v(2)-v(1))*u(1)/dx(1) 
           	else
               vv(1)=edge1(irw,3-iu)
           	end if
            h(1)=(pp(1)-qq(1))/4
           	u(1)=(pp(1)+qq(1))/2
      endif
      flood=depth(n)+edge2(irw,3)
	if((.not.newland(n-1)).and.(.not.newland(n))
     &	.and.(flood.gt.ground)) then
        	pj=Pinv(n)
        	qj=Qinv(n)
        	p1=Pinv(n-1)
        	q1=Qinv(n-1)
	  	sphr=(pj-qj)*(pj+qj)*scl(n)
        	Qjm=(3*(p1+pj)+q1+qj)*(pj-p1)/8-grav*(depth(n)-depth(n-1)) 
        	pp(n)=pj-dt*(Qjm/dx(n-1)-sphr)   
        	     	
        	qq(n)=edge2(irw,iu)-2*sqrt(grav*flood)

           	if(u(n).gt.0) then
               vv(n)=v(n)-dt*(v(n)-v(n-1))*u(n)/dx(n-1) 
           	else
               vv(n)=edge2(irw,3-iu)
           	end if
            h(n)=(pp(n)-qq(n))/4
           	u(n)=(pp(n)+qq(n))/2
      endif
      
      if(iu.eq.1) then
      	forall(i=1:n)
            	cel(i,irw)=h(i)
            	xvel(i,irw)=u(i)
            	yvel(i,irw)=vv(i)
            end forall
	else
      	forall(i=1:n)
	            cel(irw,i)=h(i)
      	      xvel(irw,i)=vv(i)
            	yvel(irw,i)=u(i)
            end forall
      endif        
      
 99  	continue
      end subroutine cliffs
