C  This file is part of ComputeBoundaryInput                                                             
C                                                                                          
C     Copyright (C) 2015, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cmpboundary.f        

      subroutine interpolate(sea,nx,ny,nt,mm,nn,ival) 
      use Interp
	implicit none

      real*4 sea(nx,ny,nt)
      real*8 cc,Wclm(ny,nt),Eclm(ny,nt),Srow(nx,nt),Nrow(nx,nt)
      integer i,kk,mm,nn,nx,ny,nt,ival
            
      cc=W_cx 	! West
      Wclm(1:ny,1:nt)=cc*sea(1,1:ny,1:nt)+(1-cc)*sea(2,1:ny,1:nt)  
      cc=E_cx  	! East
      Eclm(1:ny,1:nt)=cc*sea(nx-1,1:ny,1:nt)+(1-cc)*sea(nx,1:ny,1:nt)
          
      do i=1,nn
       	kk=WE_ny(i)-WE_ny(1)+1
       	cc=WE_cy(i)
       	west(i,ival,1:nt)=cc*Wclm(kk,1:nt)+(1-cc)*Wclm(kk+1,1:nt)
       	east(i,ival,1:nt)=cc*Eclm(kk,1:nt)+(1-cc)*Eclm(kk+1,1:nt)
      enddo
 
      cc=S_cy 	! South
      Srow(1:nx,1:nt)=cc*sea(1:nx,1,1:nt)+(1-cc)*sea(1:nx,2,1:nt)  
      cc=N_cy   	! North
      Nrow(1:nx,1:nt)=cc*sea(1:nx,ny-1,1:nt)+(1-cc)*sea(1:nx,ny,1:nt)
          
      do i=1,mm
       	kk=SN_nx(i)-SN_nx(1)+1
       	cc=SN_cx(i)
       	south(i,ival,1:nt)=cc*Srow(kk,1:nt)+(1-cc)*Srow(kk+1,1:nt)
       	north(i,ival,1:nt)=cc*Nrow(kk,1:nt)+(1-cc)*Nrow(kk+1,1:nt)
      enddo

      end subroutine interpolate

    
