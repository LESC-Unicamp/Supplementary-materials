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
      real*8, parameter   :: boltzmann = 1.38064852d-23
      real*8, parameter   :: avogadro  = 6.0221409d23
      real*8, parameter   :: tol       = 3.592d-4
      integer*8           :: printevery= 10**5
      real*8, parameter   :: r         = avogadro*boltzmann 
      real*8              :: rt        
      real*8              :: temperature
      real*8              :: vaporpressure
      real*8              :: molarheatvap
      real*8              :: kconst 
      real*8              :: emax,emin,deltaenergy
      real*8              :: totalloading  
      real*8, allocatable :: pressure(:),exploading(:),calcloading(:) 
      real*8, allocatable :: energylevels(:)
      real*8, allocatable :: fun(:),correction(:),totaltheta(:) 
      integer*8           :: mpoints                 
      
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
            mpoints=mpoints-7


            open(1,file="input.dat")
            read(1,*)             
            read(1,*) temperature
            read(1,*)             
            read(1,*) vaporpressure
            read(1,*)             
            read(1,*) molarheatvap
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
            kconst = vaporpressure*dexp(molarheatvap/rt) 
            emax = 50000.d0
            emin = 20000.d0  
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

subroutine energydistribution()
use parameters
implicit none
      real*8              :: theta_ij
      real*8              :: error
      integer*8           :: i,j,jj,jk,iter
      
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
      
      do while (error > tol)
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

      call outputs()
      write(*,*) ''
      write(*,*) ''
      call int_riemman()
      call int_trapez()
      call int_gauss()

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

      real*8                 :: sat,sat1
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
      sat1=0.d00
      do jk = 0,mpoints-1
        do i = 1,7
          sat=sat + shp(1,i)*wht(i)*0.5d0*deltaenergy*fun(jk)
          sat=sat + shp(2,i)*wht(i)*0.5d0*deltaenergy*fun(jk+1)
          if(jk.lt.mpoints*0.33)then
             sat1=sat1 + shp(1,i)*wht(i)*0.5d0*deltaenergy*fun(jk)
             sat1=sat1 + shp(2,i)*wht(i)*0.5d0*deltaenergy*fun(jk+1)                      
          endif          
        enddo
      enddo
      open(2,file='satdensgauss.dat')
      write(2,*)'Saturation density --gauss--           (mmol/g):',sat*1.d3,sat1*1.d3 
      close(2)  
      write(*,'(A49,5x,1ES20.6)')'Saturation density --gauss--            (mmol/g):',sat*1.d3,sat1*1.d3  
      
end subroutine int_gauss

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
program main
use parameters
implicit none
      call energydistribution()
      deallocate(fun,energylevels,calcloading,exploading,pressure)
end program         
