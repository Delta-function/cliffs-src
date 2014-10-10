C This file is part of Cliffs.                                                             
C                                                                                          
C     Copyright (C) 2014, Elena Tolkova                                                    
C     For conditions of distribution and use, see copyright notice in cliffs_main.f        

      subroutine allabout(indir,params,ingrids,error)

      use PARAMETERS
      use MASTER, only: nXn,nYn
      use NESTED, only: nests
      implicit none

      character*200 indir,ingrids(99),params,fname
      integer error,num,i

	error=0
	fname=' '
	
      fname = trim(indir)//trim(params)
      write(9,*) 'Read Computational parameters: ',trim(fname)
      open(unit=5,file=trim(fname),status='old',form='formatted')
      
      read(5,*) num
      if(num.ne.1) then
      	cartesian=.false.
      	write(9,*) 'SPHERICAL (LON,LAT) COORDINATE SYSTEM'
      else
      	cartesian=.true. 
       	write(9,*) 'CARTESIAN COORDINATE SYSTEM'
	endif
      call MasterGrid(indir,MBathyFile,cartesian,error)        
      if (error.NE.0) goto 99
      read(5,*) nests
      write(9,*) 'Number of grids enclosed in Master Grid: ',nests
      if((nests.gt.0).and.((nXn.lt.3).or.(nYn.lt.3))) then
      	error=-1
      	write(9,*) 'nesting in 1D is not permitted'
      	goto 99
      endif    
      if (nests.GT.0) call NestedGrid(indir,ingrids,error)
      if (error.NE.0) goto 99   
      read(5,*) cuke
      write(9,*) 'Still sea threshold (m): ', cuke
      read(5,*) ground
      write(9,*) 'Minimal flow depth (m): ',ground
      read(5,*) crough
      write(9,*) 'Friction coefficient (drag or Manning n**2): ',crough
      read(5,*) itopo
      read(5,*) dwall   ! not used if itopo.ne.0
     	if (itopo.eq.0) then
		write(9,*) 'Vertical wall at depth (m): ',dwall      
     	else
		write(9,*) 'Land inundation enabled'
	endif      
      read(5,*) dt
      write(9,*) 'Time step (sec): ',dt
      read(5,*) steps_total
      write(9,*) 'Number of steps in time loop: ',steps_total    
      read(5,*) quake
      if (quake.eq.1) then
      	write(9,*) 'Apply initial deformation to the bottom'
      else
      	write(9,*) 'Apply initial deformation to the sea surface'
      endif
      read(5,*) num
      if(num.ne.0) then
      	freeze=.false.
      	else
      	freeze=.true. 
      endif     	
      write(9,*) 'Stop computations when input stops (true/false) ',
     & 			freeze
      read(5,*)  seaout
      write(9,*)'Save screenshots every',seaout*dt,'s'
      read(5,*) lonsub
      write(9,*) 'Subsample screenshots in lon/x, nodes', lonsub
      read(5,*) latsub
      write(9,*) 'Subsample screenshots in lat/y, nodes', latsub
      read(5,*)  bndout
      write(9,*) 'Save feed into nested grids every', bndout*dt,'s'
      read(5,*)  maxout
      write(9,*) 'Update maxwave every', maxout*dt,'s'
    	read(5,*) Ngages
      if(Ngages.ge.1) then
      	read(5,*) gout
      	allocate(Igages(Ngages),Jgages(Ngages))
      	do i=1,Ngages
      		read(5,*) Igages(i),Jgages(i)
      	enddo
      	write(9,*) 'Save wave height at',Ngages,'gages every',
     & 		gout*dt,'s'
!      	write(9,*) (Igages(i), Jgages(i), i=1,Ngages)
      endif           	
      close(unit=5,status='keep')

 99  	call flush(9) 
      end subroutine allabout
