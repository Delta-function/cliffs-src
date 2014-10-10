C interpolates qq(nx,ny) on lattice xq&yq onto pp(kx,ky) on lattice xp&yp                  !
C This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine apply_initial_conditions(nx,ny,xq,yq,qq,
     &	kx,ky,xp,yp,pp)

	implicit none
	real*8 xq(nx),yq(ny),qq(nx,ny)
	real*8 xp(kx),yp(ky),pp(kx,ky)
	real*8 y1,y2,x1,x2,zero,xcc,ycc
      integer nx,ny,kx,ky,i,j,ii,jj	
	real*8 yy,xx,q

	zero=0.1e-20
	x1=xq(1)-zero   	! source grid limits
	x2=xq(nx)+zero
	y1=yq(1)-zero
	y2=yq(ny)+zero
	if((nx.gt.2).and.(ny.gt.2)) then	! 2D domain
	do 20 i=1,kx
		xx=xp(i)
		if((xx.ge.x1).and.(xx.le.x2)) then
			call pinpoint(xq,nx,xx,ii,xcc)
		else
			goto 20
		endif
		do 10 j=1,ky	
			yy=yp(j)
			if((yy.ge.y1).and.(yy.le.y2)) then
				call pinpoint(yq,ny,yy,jj,ycc)
			else
				goto 10
			endif
			q=(qq(ii,jj)  *xcc+qq(ii+1,jj)  *(1-xcc))*ycc+
     &		  (qq(ii,jj+1)*xcc+qq(ii+1,jj+1)*(1-xcc))*(1-ycc)			
			pp(i,j)=pp(i,j)+q
 10		continue ! j=1,ky			
 20	continue ! i=1,kx
	endif

	if((ny.eq.1).and.(nx.gt.2)) then	! 1D, x-axis
		do 30 i=1,kx
			xx=xp(i)
			if((xx.ge.x1).and.(xx.le.x2)) then
				call pinpoint(xq,nx,xx,ii,xcc)
			else
				goto 30
			endif
			q=qq(ii,1)*xcc+qq(ii+1,1)*(1-xcc)	
			pp(i,1:ky)=pp(i,1:ky)+q
 30		continue ! i=1,kx
	endif
	
	if((nx.eq.1).and.(ny.gt.2)) then	! 1D, y-axis
		do 40 j=1,ky	
			yy=yp(j)
			if((yy.ge.y1).and.(yy.le.y2)) then
				call pinpoint(yq,ny,yy,jj,ycc)
			else
				goto 40
			endif
			q=qq(1,jj)*ycc+qq(1,jj+1)*(1-ycc)
			pp(1:kx,j)=pp(1:kx,j)+q
 40		continue ! j=1,ky
 	endif			

	end subroutine apply_initial_conditions
