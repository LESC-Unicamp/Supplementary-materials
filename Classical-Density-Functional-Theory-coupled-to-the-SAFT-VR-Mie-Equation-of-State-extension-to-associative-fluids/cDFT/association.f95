!-------------------------------------------------------------------------------------------------!
! Copyright 2026 André de Freitas Gonçalves, Nathan Barros de Souza, and Luís Fernando Mercier    !
! Franco                                                                                          !
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
subroutine assoc(j,n0,n2,n2vec,n3,xai,ddn0,ddn2,ddn2vec,ddn3)
use globalvar
implicit none
 
integer*4 :: i,j,k,a,b,ii,jj,bb,s,iter
real*8, allocatable :: xtrial(:),dx(:),xsiold(:),xsi(:),dxi(:),dxiold(:)
real*8, allocatable :: deltaab(:,:),ddeltaab(:,:,:),k_ij(:,:)
real*8, allocatable :: sum0x(:),dif(:)
real*8, allocatable :: dxdn0(:),dxdn2(:),dxdn2vec(:),dxdn3(:)
real*8              :: rhostar,tstar,zeta,rho,xai(sites),dxsold(sites)
real*8              :: delta,ddelta(4),diff,temp,temp2,rhostarn3
real*8              :: n0,n2,n2vec,n3,ddn0,ddn2,ddn2vec,ddn3
real*8              :: theta,thetan0,thetan2,thetan3,thetan2vec,dzetadn2,dzetadn2vec


allocate(xsi(sites),dxi(sites),dxiold(sites))
allocate(deltaab(sites,sites))
allocate(ddeltaab(sites,sites,4))
allocate(sum0x(sites),dif(sites))
allocate(xsiold(sites))
allocate(xtrial(sites))
allocate(dx(sites))
allocate(dxdn0(sites),dxdn2(sites),dxdn2vec(sites),dxdn3(sites))

zeta=1.d0-n2vec*n2vec/n2/n2
rho=n0*zeta
rhostar = rho*sigma_fluid_fluid(1)**3.d00
rhostarn3 = n3*6.d00/pi
tstar = 1.d0/epsilon_fluid_fluid(1)
dzetadn2 = 2.d0*n2vec*n2vec/(n2*n2*n2)
dzetadn2vec = -2.d0*n2vec/(n2*n2)  

!Calculating delta_abij
deltaab(:,:)=0.d0
ddeltaab(:,:,:)=0.d0
do a=1,sites
   do b=1,sites
      if(a.eq.b)then
         delta=0.d00
         ddelta(:)=0.d00 
       else
         call calcdelta(n0,n2,n2vec,n3,tstar,rhostar,rhostarn3,delta,ddelta) 
       endif
       deltaab(a,b)=delta
       do k=1,4
          ddeltaab(a,b,k)=ddelta(k) 
       enddo
   enddo
enddo

!-----------------------------------------------------------------------
!Analytical solution for xai,dxdn0,dxdn2,dxdn2vec (Xa = Xb, pure water)
xsi(:)=1.d00
   temp = rho*deltaab(1,2)*4.d0
   xsi(1) = -1.d0+(1.d0+2.d0*temp)**0.5d0
   xsi(1) = xsi(1)/temp
   xsi(2) = xsi(1) 
   if(xsi(1).gt.1.d00)then  
      xsi(1)=1.d00 
      xsi(2)=xsi(1)
   endif   
   if(xsi(1).lt.1.d-20.or.abs(rho).lt.1.d-15)then
      xsi(:) = 0.99999
      dxdn0(:)=1.d-100
      dxdn2(:)=1.d-100   
      dxdn2vec(:)=1.d-100
      dxdn3(:)=1.d-100
      ddn0=1.d-100
      ddn2=1.d-100
      ddn2vec=1.d-100
      ddn3=1.d-100 
      go to 132     
   else
   temp2 = 1.d0+2.d0*temp
   temp2 = temp2**0.5d0
   temp2 = 1.d0/temp2
!dxdn0 
   dxdn0(1) = -xsi(1)*(1.d0/n0+ddeltaab(1,2,1)/deltaab(1,2))
   dxdn0(1) = dxdn0(1) + temp2*(1.d0/n0+ddeltaab(1,2,1)/deltaab(1,2))
   dxdn0(2) = dxdn0(1)   
!dxdn2   
   dxdn2(1) = -xsi(1)*(dzetadn2/zeta+ddeltaab(1,2,2)/deltaab(1,2))
   dxdn2(1) = dxdn2(1) + temp2*(dzetadn2/zeta+ddeltaab(1,2,2)/deltaab(1,2))
   dxdn2(2) = dxdn2(1) 
!dxdn2vec 
   dxdn2vec(1) =-xsi(1)*(dzetadn2vec/zeta+ddeltaab(1,2,3)/deltaab(1,2))
   dxdn2vec(1) = dxdn2vec(1) + temp2*(dzetadn2vec/zeta+ddeltaab(1,2,3)/deltaab(1,2))
   dxdn2vec(2) = dxdn2vec(1)   
!dxdn3   
   dxdn3(1) = ddeltaab(1,2,4)/deltaab(1,2)*(-xsi(1)+temp2)
   dxdn3(2) = dxdn3(1)           
   endif

!-----------------------------------------------------------------------
theta = 0.d0
thetan0 = 0.d0
thetan2 = 0.d0
thetan2vec = 0.d0
thetan3 = 0.d0
do i = 1,sites
   theta = theta + ns(i)*(dlog(xsi(i))-xsi(i)*0.5d0+0.5d0)
   thetan0 = thetan0 + ns(i)*(dxdn0(i)*1.d0/xsi(i)-dxdn0(i)*0.5d0) 
   thetan2 = thetan2 + ns(i)*(dxdn2(i)*1.d0/xsi(i)-dxdn2(i)*0.5d0)
   thetan2vec = thetan2vec + ns(i)*(dxdn2vec(i)*1.d0/xsi(i)-dxdn2vec(i)*0.5d0) 
   thetan3 = thetan3 + ns(i)*(dxdn3(i)*1.d0/xsi(i)-dxdn3(i)*0.5d0)      
enddo
ddn0 = theta*zeta+rho*thetan0
ddn2 = theta*n0*dzetadn2+rho*thetan2
ddn2vec = theta*n0*dzetadn2vec+rho*thetan2vec
ddn3 = thetan3*rho
132 continue
xai(:) = xsi(:)


deallocate(xsi,dxi,dxiold)
deallocate(deltaab)
deallocate(ddeltaab)
deallocate(sum0x,dif)
deallocate(xsiold)
deallocate(xtrial)
deallocate(dx)
deallocate(dxdn0,dxdn2,dxdn2vec,dxdn3)
end subroutine

subroutine calcdelta(n0,n2,n2vec,n3,tstar,rhostar,rhostarn3,delta,ddelta)   
use globalvar
implicit none
real*8 :: FABIJ,REPS,delta,ddelta(4),dfun(4),fun
real*8 :: tstar,rhostar,rhostarn3,n0,n2,n2vec,n3,rc,rd,t1,t2,t3,t4,t5

      rc = rcbysigma(1)*sigma_fluid_fluid(1)
      rd = rdbysigma(1)*sigma_fluid_fluid(1)
      call kernel (n0,n2,n2vec,n3,tstar,rhostar,rhostarn3,fun,dfun)
      REPS     = EPSHBIJ(1)/kB/T
      FABIJ    = EXP(REPS)-1.D0
      t1 = (rc+2.d0*rd)/dij(1,1)
      t2 = rc**3.d0+3.d0*rc*rc*rd-4.d0*rd**3.d0
      t3 = rc+2.d0*rd-dij(1,1)
      t4 = 22.d0*rd*rd-5.d0*rc*rd-7.d0*rd*dij(1,1)-8.d0*rc*rc+rc*dij(1,1)+dij(1,1)*dij(1,1) 
      t5 = 72.d0*rd*rd*sigma_fluid_fluid(1)**3.d0 
      KABIJ(1) = 4.d0*pi*dij(1,1)*dij(1,1)*(dlog(t1)*6.d0*t2+t3*t4)/t5
      delta    = sigma_fluid_fluid(1)**3.d00*FABIJ*KABIJ(1)*fun
      ddelta(:) = sigma_fluid_fluid(1)**3.d00*FABIJ*KABIJ(1)*dfun(:)  
end subroutine

subroutine kernel (n0,n2,n2vec,n3,tstar,rhostar,rhostarn3,fun,dfun)
use globalvar
implicit none

real*8 :: dfun(4),fun
real*8 :: dfundn0,dfundn2,dfundn2vec,dfundn3
real*8 :: zeta,dzetadn2,dzetadn2vec,rhostarn3
real*8 :: tstar,rhostar,n0,n2,n2vec,n3,temp,tempn2
   
      zeta = 1.d0-n2vec*n2vec/(n2*n2)
      dzetadn2 = 2.d0*n2vec*n2vec/(n2*n2*n2)
      dzetadn2vec = -2.d0*n2vec/(n2*n2)  
      fun = 0.D0
      dfundn0=0.d0
      dfundn2=0.d0
      dfundn2vec=0.d0 
      dfundn3 = 0.d0   
      temp = 1.d0/(1.d0-n3)
      tempn2 = 2.d0*n2*zeta+n2*n2*dzetadn2
      fun = temp+0.25d0*n2*zeta*dij(1,1)*temp*temp+dij(1,1)*dij(1,1)*0.25d0/18.d0*temp**3.d0*n2*n2*zeta    
      dfundn0 = 0.d0
      dfundn2 = 0.25d0*dij(1,1)*temp*temp*(zeta+n2*dzetadn2)+0.25d0/18.d0*dij(1,1)*dij(1,1)*temp**3.d0*tempn2 
      dfundn2vec = 0.25d0*dij(1,1)*temp*temp*n2*dzetadn2vec+                 &
                    0.25d0/18.d0*dij(1,1)*dij(1,1)*temp**3.d0*n2*n2*dzetadn2vec
      dfundn3 = temp*temp+n2*zeta*dij(1,1)*0.5d0*temp**3.d0+1.d0/24.d0*dij(1,1)*dij(1,1)*n2*n2*zeta*temp**4.d0
      
      dfun(1) = dfundn0
      dfun(2) = dfundn2
      dfun(3) = dfundn2vec
      dfun(4) = dfundn3
end subroutine
