!-------------------------------------------------------------------------------------------------!
! Copyright 2025 André Gonçalves, Emerson Parazzi Lyra, Sayali Ramdas Chavan,                     ! 
! Philip L. Llewellyn, Luis Fernando Mercier Franco, and Yann Magnin.                             !
!                                                                                                 !
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software   !
! and associated documentation files (the “Software”), to deal in the Software without            !
! restriction, including without limitation the rights to use, copy, modify, merge, publish,      !
! distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the   !
! Software is furnished to do so, subject to the following conditions:                            !
!                                                                                                 !
! The above copyright notice and this permission notice shall be included in all copies or        !
! substantial portions of the Software.                                                           !
!                                                                                                 !
! THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,             !
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR        !
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR   !
! ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,     !
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN         !
! THE SOFTWARE.                                                                                   !
!-------------------------------------------------------------------------------------------------!

module parameters
implicit none
      real*8, parameter    :: boltzmann  = 1.38064852d-23
      real*8, parameter    :: avogadro   = 6.0221409d23
      real*8, parameter    :: tol        = 3.592d-4
      integer*8, parameter :: maxiter    = 2.0d5
      integer*8            :: printevery = 10**5
      real*8, parameter    :: r          = avogadro*boltzmann 
      real*8               :: rt        
      real*8               :: temperature
      real*8               :: error
      real*8               :: deltaHinf
      real*8               :: kconst 
      real*8               :: emax,emin,deltaenergy
      real*8               :: totalloading  
      real*8, allocatable  :: pressure(:),exploading(:),calcloading(:) 
      real*8, allocatable  :: energylevels(:)
      real*8, allocatable  :: fun(:),correction(:),totaltheta(:) 
      integer*8            :: mpoints               
      
      contains 

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

      subroutine read_parameters()
      implicit none
            integer*8     :: i 
             
            open(1,file="input.dat")
            mpoints=0
            do
              read(1,*,end=100,err=200)
              mpoints=mpoints+1
            enddo
200 write(*,*) 'READING ERROR'  ;stop
100 continue
            close(1)
            mpoints=mpoints-5
            open(1,file="input.dat")
            read(1,*)             
            read(1,*) temperature
            read(1,*)             
            read(1,*) deltaHinf
            read(1,*)  

            allocate(pressure(mpoints))
            allocate(exploading(mpoints)) 

            totalloading = 0.d0
                  do i = 1,mpoints
                     read(1,*) pressure(i),exploading(i)
                     totalloading = totalloading+exploading(i) 
                  enddo
            close(1)

            allocate(energylevels(mpoints))  
            allocate(calcloading(mpoints)) 
            allocate(fun(mpoints))  
            allocate(correction(mpoints)) 
            allocate(totaltheta(mpoints))            

            energylevels(:)=0.d0                                 
            rt = r*temperature
            emax = 40000.d0
            emin = 5000.d0  
            deltaenergy = (emax-emin)/dble(mpoints) 
                  do i = 1,mpoints
                     energylevels(i) = emin+deltaenergy*(i-1)
                  enddo       
      end subroutine read_parameters
      
end module parameters

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine langmuir(p_i,eps_i,theta_i)
use parameters
implicit none
      real*8              :: theta_i
      real*8              :: p_i
      real*8              :: eps_i 

      theta_i = p_i+kconst*dexp(-eps_i/rt)
      theta_i = p_i/theta_i
      
end subroutine langmuir

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine energydistribution(prefactor)
use parameters
implicit none
      real*8              :: theta_ij,prefactor(1),adsHcalc
      integer*8           :: i,ii,j,jj,jk,iter,adsHloc
      
      kconst = prefactor(1)*1.d9
      
      call read_parameters()
      
      !TO READ THE INITIAL GUESS FROM A FILE, UNCOMMENT FOLLOWING LINES:
      !open(2,file='initialguess.dat')
      !do i = 1,mpoints
      !   read(2,*) fun(i)
      !enddo
      !close(2) 
      
      fun(:)=exploading(mpoints)/mpoints/deltaenergy

      
      error = 1.d1
      iter = 0
      
      do while (iter < maxiter)
         iter = iter+1
         calcloading(:) = 0.d0
         do i = 1,mpoints
            do j = 1,mpoints   
               call langmuir(pressure(i),energylevels(j),theta_ij)          
               calcloading(i) = calcloading(i)+fun(j)*theta_ij*deltaenergy
            enddo
         enddo

         correction(1:mpoints) = 0.d0
         totaltheta(1:mpoints) = 0.d0
         error = 0.d0

         do jk = 1,mpoints

            error = 0.d0

            do jj = 1,mpoints
               call langmuir(pressure(jj),energylevels(jk),theta_ij)              
               correction(jk) = correction(jk)+theta_ij*deltaenergy*exploading(jj)/calcloading(jj)
               totaltheta(jk) = totaltheta(jk)+theta_ij*deltaenergy
               error = error+(calcloading(jj)-exploading(jj))*(calcloading(jj)-exploading(jj)) 
            enddo
            fun(jk) = fun(jk)*correction(jk)/totaltheta(jk)   
         enddo
         error = error**0.5d0
         write(*,*)'iter,error',iter,error
         if(mod(iter,printevery)==0)then
            write(*,'(1A15,5x,1I12,5x,1ES20.6)')'iteration,error',iter,error 
            open(1,file = 'energydistribution.dat')
            write(1,*) '#Energy (kJ/mol),fun(energy),calculated isotherm (mmol/g)'      
            do i = 1,mpoints
               write(1,'(3ES20.6)') energylevels(i)*1.d-3,fun(i),calcloading(i)*1.d3
            enddo
         

            
            close(1) 
            call int_riemman()
            call int_trapez()
            call int_gauss()
         endif                     
      enddo

      ii = mpoints
      do while (fun(ii).le.fun(ii-1))
         ii=ii-1
         adsHcalc = energylevels(ii)
      enddo
      adsHcalc = energylevels(ii)  
      error = error + abs(-adsHcalc*1.d-3+deltaHinf)

      call outputs()
      deallocate(fun,energylevels,calcloading,exploading,pressure,correction,totaltheta)
end subroutine energydistribution

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine outputs()
use parameters
implicit none

      integer*8           :: i

      open(1,file = 'energydistribution.dat')
      write(1,*) '#Energy (kJ/mol),fun(energy),calculated isotherm (mmol/g)'      
      do i = 1,mpoints
         write(1,'(3ES20.6)') energylevels(i)*1.d-3,fun(i),calcloading(i)*1.d3
      enddo
      close(1)

      open(10,file = 'restart.dat')
      do i = 1,mpoints
         write(10,'(1ES20.6)') fun(i)
      enddo
      close(10)
      
end subroutine outputs

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine int_riemman()
use parameters
implicit none

      real*8              :: sat
      integer*8           :: jk

      sat = 0.d0
      do jk = 1,mpoints 
        sat = sat+fun(jk)*deltaenergy 
      enddo
      open(2,file='satdensrieamman.dat')
      write(2,*)'Saturation density --riemman--          (mmol/g):',sat*1.d3
      close(2)

      write(*,'(A49,5x,1ES20.6)')'Saturation density --riemman--          (mmol/g):',sat*1.d3
end subroutine int_riemman

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine int_trapez()
use parameters
implicit none

      real*8              :: sat,smin,smax
      integer*8           :: jk

      sat = 0.d0
      do jk = 0,mpoints-1

        if(fun(jk)<fun(jk+1))then
          smin=fun(jk)
          smax=fun(jk+1)
        else
          smin=fun(jk+1)
          smax=fun(jk)
        endif

        sat = sat + deltaenergy*smin + 0.5d0*deltaenergy*dabs(smax-smin) 

      enddo
      open(2,file='satdenstrapez.dat')
      write(2,*)'Saturation density --trapez--           (mmol/g):',sat*1.d3
      close(2)      

      write(*,'(A49,5x,1ES20.6)')'Saturation density --trapez--           (mmol/g):',sat*1.d3
      
end subroutine int_trapez

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine int_gauss()
use parameters
implicit none

      real*8                 :: sat,sat1,sat2
      integer*8              :: jk,i

      real*8, dimension(7)   :: pnt,wht
      real*8, dimension(2,7) :: shp

      pnt(1) =0.d0
      wht(1) =0.417959183673469d0

      pnt(2) =0.405845151377397d0
      wht(2)=0.381830050505119d0

      pnt(3) =-0.405845151377397d0
      wht(3)=0.381830050505119d0

      pnt(4) =0.741531185599394d0
      wht(4)=0.279705391489277d0

      pnt(5) =-0.741531185599394d0
      wht(5)=0.279705391489277d0

      pnt(6) =0.949107912342759d0
      wht(6)=0.129484966168870d0

      pnt(7) =-0.949107912342759d0
      wht(7)=0.129484966168870d0

      do i=1,7
        shp(1,i)=0.5d0*(1.d0-pnt(i))
        shp(2,i)=0.5d0*(1.d0+pnt(i))
      enddo

      sat = 0.d0
      sat1=0
      sat2=0
      do jk = 0,mpoints-1
        do i = 1,7
          sat=sat + shp(1,i)*wht(i)*0.5d0*deltaenergy*fun(jk)
          sat=sat + shp(2,i)*wht(i)*0.5d0*deltaenergy*fun(jk+1)
        enddo
      enddo
      open(2,file='satdensgauss.dat')
      write(2,*)'Saturation density --gauss--           (mmol/g):',sat*1.d3
      close(2)  
      write(*,'(A49,5x,1ES20.6)')'Saturation density --gauss--            (mmol/g):',sat*1.d3  
      
end subroutine int_gauss

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
program main
use parameters
implicit double precision(A-Z)
integer K,N,LS
real*8, allocatable ::  X1(:),X2(:),X3(:),XK(:)
real*8, allocatable ::  XMIN(:),XMAX(:),EMAXNM(:)
character DUM*40,NFILE*40
    
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Reading the initial estimates                                        C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(1,file='nelder_mead.dat')
      read(1,*) DUM    
      read(1,*) N 
      read(1,*) DUM
      allocate(X1(N),X2(N),X3(N),XK(N))
      allocate(XMIN(N),XMAX(N),EMAXNM(N))      
      do  K=1,N
         read(1,*) X1(K)
      enddo    
      read(1,*) DUM
      do  K=1,N
         read(1,*) XMIN(K),XMAX(K)
      enddo
      close(1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Minimization of the object function                                  C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(*,*) ' Minimizing the objective function ... '
      EFAC = 1E-3
      call energydistribution(X1)
      Y    = error
      write(*,*)'error',error
      do K=1,N
          XK(K)   = 1.0
          EMAXNM(K) = 0.01*(XMAX(K)-XMIN(K))
          X2(K)   = X1(K)
          X3(K)   = X1(K)
      enddo
      YSTAR = Y
      LS    = 0
  100 YTEST = YSTAR
      do K=1,N
         if ((X1(K) .ge. XMAX(K)) .and. (XK(K) .gt. 0.0)) then 
             go to 1010 
         endif
         if ((X1(K) .le. XMIN(K)) .and. (XK(K) .lt. 0.0)) then 
             go to 1010 
         end if
         X2(K) = X1(K)+XK(K)*EMAXNM(K)
         if (X2(K) .gt. XMAX(K)) then 
             X2(K) = XMAX(K) 
         end if
         if (X2(K) .lt. XMIN(K)) then 
             X2(K) = XMIN(K) 
         end if
         call energydistribution(X2)
         Y1 = error
         if (Y1 .ge. YTEST) then 
             go to 1010  
         end if
         XK(K) = 2.0*XK(K)
         if (XK(K) .gt. 1.0) then
             XK(K) = 1.0  
         end if
         if (XK(K) .lt. -1.0) then 
             XK(K) = -1.0 
         end if
         go to 103
 1010    if ((X1(K) .le. XMIN(K)) .and. (XK(K) .gt. 0.0)) then 
            go to 1011 
         end if
         if ((X1(K) .ge. XMAX(K)) .and. (XK(K) .LT. 0.0)) then  
            go to 1011 
         end if
         X2(K) = X1(K)-XK(K)*EMAXNM(K)
         if (X2(K) .lt. XMIN(K)) then 
             X2(K) = XMIN(K) 
         end if
         if (X2(K) .gt. XMAX(K)) then 
             X2(K) = XMAX(K) 
         end if
         call energydistribution(X2)
         Y1 = error;
         if (Y1 .ge. YTEST) then 
             go to 1011 
         end if
         XK(K) = -1.0*XK(K)
         go to 103
 1011    X2(K) = X1(K)
         XK(K) = 0.5*XK(K)
         if (XK(K) .lt. 0.0) then 
             go to 1015 
         end if
         if (XK(K) .lt. 0.5*EFAC) then 
             XK(K) = 0.5*EFAC 
         end if
         go to 1018
 1015    if (ABS(XK(K)) .lt. 0.5*EFAC) then
             XK(K) = -0.5*EFAC 
         end if
 1018    Y1    = YTEST
  103    YTEST = Y1
      enddo  
      if (LS .eq. 1) then 
          go to 350 
      end if
      YBASE = YTEST
      if (YBASE .lt. YSTAR) then 
          go to 300 
      end if
      DO  K=1,N
         if (ABS(XK(K)) .gt. 0.5*EFAC) then 
             go to 100
         end if
      enddo
      go to 600
 300  DO K=1,N
         X1(K) = X3(K)+2.0*(X2(K)-X3(K))
         if (X1(K) .lt. XMIN(K)) then 
             X1(K) = XMIN(K)
         end if
         if (X1(K) .gt. XMAX(K)) then 
             X1(K) = XMAX(K) 
         end if
         X3(K) = X2(K)
         X2(K) = X1(K)
      enddo
      call energydistribution(X1)
      YSTAR = error
      LS    = 1
      go to 100
 350  LS    = 0
      if (YTEST .ge. YBASE) then 
          go to 500 
      end if
      YBASE = YTEST
      go to 300
 500  do K=1,N
         X1(K) = X3(K)
         X2(K) = X1(K)
      enddo
      YSTAR = YBASE
      go to 100
      call energydistribution(X1)
 600  Y = error
 
      write(*,'(A9,E15.7)') 'FUN  = ',Y
      do K=1,N
         write(*,'(A4,I1,A3,E15.7)') ' X[',K,'] = ',X1(K)
      enddo
     write(*,*) ' The calculation is ended ...'   
end program     
