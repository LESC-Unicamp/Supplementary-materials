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
module globalvar
        implicit none
        integer*8            :: nc
        integer*8            :: n_bins             
        integer*8, parameter :: ext_bins            = 8192
        real*8, parameter    :: kB                  = 1.38064852D-23
        real*8, parameter    :: avogadro            = 6.0221409d23        
        real*8, parameter    :: pi                  = 4.d0*datan(1.d0)
        real*8, parameter    :: solid_density       = 0.114d0 
        real*8, parameter    :: interlayer_spacing  = 3.35d0
        real*8, parameter    :: epsilon_solid_solid = 28.0d00
        real*8, parameter    :: sigma_solid_solid   = 3.4d0
        real*8, parameter    :: width               = 7.0d0         !Width of slit-pore          
        real*8               :: z0
        real*8               :: beta        
        integer*8            :: contrl
        real*8               :: T,inv_T
        real*8               :: chem_pot
        real*8               :: delta_z
        real*8               :: mole_fraction
        real*8               :: total_bulk_density         
        real*8, allocatable  :: bulk_density_i(:)  
        real*8, allocatable  :: density_calc(:,:) 
        real*8, allocatable  :: epsilon_fluid_fluid(:)
        real*8, allocatable  :: sigma_fluid_fluid(:) 
        real*8, allocatable  :: epsilon_solid_fluid(:) 
        real*8, allocatable  :: sigma_solid_fluid(:) 
        real*8, allocatable  :: di(:),rcbysigma(:),rdbysigma(:)  
        real*8, allocatable  :: m(:)  
        real*8, allocatable  :: dxstoren0(:),dxstoren2(:),dxstoren2vec(:)
        real*8, allocatable  :: n0global(:),n2global(:),n2vecglobal(:),n3global(:)  
        integer, allocatable :: ns(:)                      
        real*8               :: chemical_potential 
        real*8, allocatable  :: exc_chempot(:) 
        real*8               :: ideal_gas_chemical_potential
        real*8               :: length
        real*8               :: inv_length  
        real*8               :: lambdar,lambdaa     
        real*8, allocatable  :: c1(:,:,:),c2(:,:,:),c3(:,:,:),c4(:,:,:),psi(:) 
        real*8, allocatable  :: c1a2(:,:,:),c2a2(:,:,:),c3a2(:,:,:),c4a2(:,:,:)          
        real*8, allocatable  :: lambdaij(:,:,:),c_ij(:,:),F_alpha(:,:,:)
        real*8, allocatable  :: dij(:,:),dij3(:,:),x0ij(:,:),eps(:,:),epshbij(:),kabij(:),sigmaij3(:,:)
        integer*8            :: tests,sites
end module globalvar

program cdft
        use globalvar
        implicit none
        integer*8            :: k,i,j,ii,half_n_bins,iteration
        real*8, parameter    :: tolerance = 1.5d-7
        real*8, parameter    :: alpha     = 0.008d0
        real*8               :: avg_error
        real*8               :: z
        real*8               :: CIJ,LAIJ,LRIJ,F(6)        
        real*8, allocatable  :: solid_fluid_potential(:)
        real*8, allocatable  :: density_z(:) 
        real*8, allocatable  :: density_star(:,:)
        real*8, allocatable  :: dF_hard_sphere(:,:),dF_At(:,:),dF_chain(:,:)
        real*8, allocatable  :: dF_assoc(:,:)
        real*8, allocatable  :: density_old(:,:)
        real*8, allocatable  :: density_new(:,:)
        real*8, allocatable  :: density(:,:)
        real*8, allocatable  :: error(:)        
     
        open (1,file="input.dat")
        read (1,*) T
        read (1,*) nc
        allocate(psi(nc))
        allocate(rcbysigma(nc),rdbysigma(nc))
        allocate(epsilon_fluid_fluid(nc))
        allocate(epsilon_solid_fluid(nc))
        allocate(sigma_fluid_fluid(nc))          
        allocate(sigma_solid_fluid(nc)) 
        allocate(bulk_density_i(nc)) 
        allocate(lambdaij(2,nc,nc))
        allocate(m(nc)) 
        allocate(epshbij(nc),kabij(nc))
        allocate(n0global(ext_bins),n2global(ext_bins),n3global(ext_bins))
        allocate(n2vecglobal(ext_bins))             
        do i=1,nc
           read (1,*) bulk_density_i(i), sigma_fluid_fluid(i),         &
                      epsilon_fluid_fluid(i),lambdaij(1,i,i),          &
                      lambdaij(2,i,i),m(i),psi(i),epshbij(i),          &
                      rcbysigma(i),rdbysigma(i),sites
          allocate(ns(sites))
          read (1,*) (ns(j),j=1,sites)
           sigma_solid_fluid(i) = 0.5d0*(sigma_solid_solid             &
                                +sigma_fluid_fluid(i))
           epsilon_solid_fluid(i) = dsqrt(epsilon_solid_solid          &
                               *epsilon_fluid_fluid(i))
           !For changing the value of solid-fluid parameters directly, uncomment the next two lines:    
           !epsilon_solid_fluid(i) = 91.717d00
           !sigma_solid_fluid(i)=3.2165d00                    
        enddo
        close(1)
        epshbij(:) = epshbij(:)*kB
        lambdaa = lambdaij(1,1,1)        
        lambdar = lambdaij(2,1,1)
        
        allocate(di(nc))
        allocate(density_z(nc))        
        allocate(solid_fluid_potential(nc))
        allocate(exc_chempot(nc))
        allocate(error(nc)) 
        allocate(dxstoren0(sites),dxstoren2(sites),dxstoren2vec(sites))
        
        allocate(c1(2,nc,nc))
        allocate(c2(2,nc,nc))            
        allocate(c3(2,nc,nc))
        allocate(c4(2,nc,nc))
        
        allocate(c1a2(3,nc,nc))
        allocate(c2a2(3,nc,nc))            
        allocate(c3a2(3,nc,nc))
        allocate(c4a2(3,nc,nc))          

        allocate(dij(nc,nc))  
        allocate(dij3(nc,nc))
        allocate(sigmaij3(nc,nc)) 
        allocate(F_alpha(nc,nc,6))       
        allocate(eps(nc,nc))
        allocate(x0ij(nc,nc))
        allocate(c_ij(nc,nc))
        
        inv_T=1.d00/T
        beta=inv_T/kB
        epsilon_solid_fluid=epsilon_solid_fluid*inv_T
        epsilon_fluid_fluid=epsilon_fluid_fluid*inv_T

        do i=1,nc
           do j=1,nc
              lambdaij(1,i,j) = 3.d00+dsqrt((lambdaij(1,i,i)-3.d00)    &
                                *(lambdaij(1,j,j)-3.d00))          
              lambdaij(2,i,j) = 3.d00+dsqrt((lambdaij(2,i,i)-3.d00)    &
                                *(lambdaij(2,j,j)-3.d00))   
           enddo
        enddo 
        di(:)                = sigma_fluid_fluid(:) 
    
        call cmix
        call combrule
        
        do i=1,nc
           do j=1,nc
              CIJ=c_ij(i,j)
              LAIJ=lambdaij(1,i,j)
              LRIJ=lambdaij(2,i,j)
              call fun_alpha(CIJ,LRIJ,LAIJ,F)
              F_alpha(i,j,:) = F
           enddo
        enddo                             
      
        n_bins               = dint((width-dij(1,1))/(width+1.0d00*dij(nc,nc))*      &
                               dble(ext_bins-1))+1

        allocate(density_old(nc,n_bins))
        allocate(density_star(nc,n_bins))
        allocate(density_new(nc,n_bins))        
        allocate(density_calc(nc,n_bins))
        allocate(dF_hard_sphere(nc,n_bins))
        allocate(dF_At(nc,n_bins))
        allocate(dF_chain(nc,n_bins))
        allocate(dF_assoc(nc,n_bins))        
        allocate(density(nc,n_bins))

        dxstoren0(:)=1.d0
        dxstoren2(:)=1.d0
        dxstoren2vec(:)=1.d0                     
                               
        if(mod(dble(ext_bins-n_bins),2.d0).ne.0) then
           n_bins = dint(n_bins-1.d0)
        endif    
        
        half_n_bins          = dint(n_bins*0.5d00)
        delta_z              = (width-dij(1,1))/dble(n_bins-1)
        z0                   = 0.5*dij(1,1)
        z = z0

        call hard_sphere (0,dF_hard_sphere)
        call monomers (0,dF_At)
        !call chain (0,dF_chain)
        call association (0,dF_assoc)        

        do i =1,nc
            exc_chempot(i)=dF_hard_sphere(i,half_n_bins)+dF_At(i,half_n_bins)+dF_assoc(i,half_n_bins)
            write(*,*) 'exc_chempot',exc_chempot(i)
        enddo

        z=z0   
        do k=1,n_bins
           call calc_initial_distribution(z,density_z)
           do ii=1,nc
              density_old(ii,k) = density_z(ii)         
           enddo               
           z              = z+delta_z
        enddo
 
        !To start from a known denity profile specified in the file "profile.dat", uncomment line 221:      
        open(1,file="profile.dat")
        z = z0
        do k=1,n_bins
           !read(1,*) z,density_old(1,k)
           z = z+delta_z
        end do 
        close(1)
        
        iteration = 0
        error(:) = 1000
        avg_error = 1000
        do while (avg_error.gt.tolerance)
           write(*,*) 'error(i) = ', (error(i),i=1,nc)
           z                 = z0
           density_calc(:,:)      = density_old(:,:)
           error(:)             = 0.d0
           
        call hard_sphere (1,dF_hard_sphere)
        call association (1,dF_assoc)           
        call monomers (1,dF_At)
        !call chain (1,dF_chain)     
           do k=1,n_bins
              do j=1,nc
                 call calc_steele(z,solid_fluid_potential)
                 density_star(j,k)   = bulk_density_i(j)               &
                                   *dexp(exc_chempot(j)                &
                                   -dF_hard_sphere(j,k) -              &
                                   dF_At(j,k)-                        &
                                   !dF_chain(j,k)-                     &
                                   dF_assoc(j,k)-                     &                                    
                                   solid_fluid_potential(j))
                 density_new(j,k)  = (1.d0-alpha)*density_old(j,k)     &
                                   +alpha*density_star(j,k)                                  
                 error(j)          = error(j)+dabs(density_new(j,k)-   &
                                   density_old(j,k))
                 avg_error         = avg_error+error(j)

              end do
              z              = z+delta_z
           end do

           z=z0
           avg_error         = avg_error*0.5d00/dble(n_bins)
           density(:,:)        = density_new(:,:)
           density_old(:,:)    = density_new(:,:)
           close(1)
        end do

       close(3)

        open(1,file="profile.dat")
        z   = z0
        do k=1,n_bins
            write(1,*) z,(density(j,k),j=1,nc)
            z = z+delta_z
        end do 
        close(1)

        deallocate(epsilon_fluid_fluid)
        deallocate(epsilon_solid_fluid)
        deallocate(sigma_fluid_fluid)          
        deallocate(sigma_solid_fluid) 
        deallocate(bulk_density_i)  
        deallocate(di)
        deallocate(psi)
        deallocate(density_z)
        deallocate(density_calc)
        deallocate(density_star)
        deallocate(dF_hard_sphere)
        deallocate(density)
        deallocate(density_old)
        deallocate(density_new)
        deallocate(solid_fluid_potential)
        deallocate(exc_chempot)
        deallocate(error) 
        
        deallocate(c1)
        deallocate(c2)            
        deallocate(c3)
        deallocate(c4) 
        deallocate(lambdaij)          
end program cdft

subroutine calc_steele(z,potential)
        use globalvar
        implicit none
        integer*8 :: i
        real*8 :: z
        real*8 :: potential(nc)           

        do i = 1,nc
                       potential(i) = 2.d0*pi*solid_density            &
                          *epsilon_solid_fluid(i)                      &
                          *sigma_solid_fluid(i)**2.d00                 &
                          *interlayer_spacing                          &
                          *(0.4d0*(sigma_solid_fluid(i)/z)**10.d0      &
                          -(sigma_solid_fluid(i)/z)**4.d0              &
                          -sigma_solid_fluid(i)**4.d0/3.d0             &
                          /interlayer_spacing                          &
                          /(z+0.61d0*interlayer_spacing)**3.d0         &
                          +0.4d0*(sigma_solid_fluid(i)                 &
                          /(width-z))**10.d0                           &
                          -(sigma_solid_fluid(i)/(width-z))**4.d0      &
                          -sigma_solid_fluid(i)**4.d0/3.d0             &
                          /interlayer_spacing                          &
                          /(width-z                                    &
                          +0.61d0*interlayer_spacing)**3.d0)
        enddo            
end subroutine calc_steele

subroutine calc_initial_distribution(z,density_z)
        use globalvar
        implicit none
        real*8 :: z
        real*8 :: potential(nc) 
        real*8 :: density_z(nc) 
        integer*8 :: i

        call calc_steele(z,potential)
        do i=1,nc
           density_z(i)  = bulk_density_i(i)*dexp(-                    &
                           potential(i))
        enddo
end subroutine calc_initial_distribution

function density_map(z) result (dens)
        use globalvar
        implicit none
        integer*8 :: k,i
        real*8    :: z,diam
        real*8    :: dens(nc)

        do i=1,nc
           diam = dij(i,i)
           if (z < diam*0.5d00) then
              dens(i) = density_calc(i,k)           
           else if (z > width-0.5d00*diam) then
              dens(i) = density_calc(i,k)              
           else
              k    = dint((z-z0)/delta_z)
              dens(i) = density_calc(i,k)
           end if
        end do
end function density_map

subroutine hard_sphere (ct,dF_hard_sphere)
        use globalvar
	implicit none
        integer*8                      :: i,j,k,ii,kk
        integer                        :: ct       
        integer*8, parameter           :: cplxp = 1
        integer*8, parameter           :: cplxn = 0       
        integer*8                      :: dif_bins
        integer*8                         :: z_index
        real*8, dimension (nc,n_bins)     :: dF_hard_sphere
        real*8, dimension (nc,ext_bins)   :: rho_ft
        real*8, dimension (ext_bins)   :: Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain 
        real*8, dimension (ext_bins)   :: n0,n1,n2,n3,n1_vec,n2_vec
        real*8, dimension (ext_bins)   :: d3,recn3,dphidn0,dphidn1,dphidn2,dphidn3,dphidn1_vec,dphidn2_vec
        real*8, dimension (ext_bins)   :: c0ft,c1ft,c2ft,c3ft,c1_vecft,c2_vecft
        real*8, dimension (ext_bins)   :: n0_sum,n1_sum,n2_sum,n3_sum,n1_vec_sum,n2_vec_sum 
   
	rho_ft(:,:) = 0.d00
        dif_bins  = dint((ext_bins-n_bins)/2.d00)
        if(ct.eq.1)then
          do i=1,nc
	     z_index   = 1           
             do j=dif_bins+1,n_bins+dif_bins
                rho_ft(i,j) = density_calc(i,z_index)
                z_index   = z_index+1
             enddo 
          enddo
          z_index   = 1
        else
          do i=1,nc
             rho_ft(i,1:ext_bins)=bulk_density_i(i)
          enddo
        endif          

        n0_sum(:)=0.d00
	n1_sum(:)=0.d00
	n2_sum(:)=0.d00
	n3_sum(:)=0.d00
	n1_vec_sum(:)=0.d00
	n2_vec_sum(:)=0.d00
        do i=1,nc
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain)
           call fft(rho_ft(i,:),Fw0,n0,ext_bins,cplxn)
           call fft(rho_ft(i,:),Fw1,n1,ext_bins,cplxn)
           call fft(rho_ft(i,:),Fw2,n2,ext_bins,cplxn)
           call fft(rho_ft(i,:),Fw3,n3,ext_bins,cplxn)
           call fft(rho_ft(i,:),Fw1_vec,n1_vec,ext_bins,cplxp) 
           call fft(rho_ft(i,:),Fw2_vec,n2_vec,ext_bins,cplxp)
           n0_sum(:)=n0_sum(:)+n0(:)*m(i)
           n1_sum(:)=n1_sum(:)+n1(:)*m(i)
           n2_sum(:)=n2_sum(:)+n2(:)*m(i)
           n3_sum(:)=n3_sum(:)+n3(:)*m(i)
           n1_vec_sum(:)=n1_vec_sum(:)+n1_vec(:)*m(i)
           n2_vec_sum(:)=n2_vec_sum(:)+n2_vec(:)*m(i)
        enddo
        do j=1,ext_bins
          if(abs(n3_sum(j)).lt.1d-20)then
             n3_sum(j) = 1d-20
          endif          
          recn3(j)       = 1d00-n3_sum(j)
          dphidn0(j)     = -dlog(recn3(j))
          dphidn1(j)     = n2_sum(j)/(recn3(j))
          dphidn2(j)     = n1_sum(j)/recn3(j)+(n2_sum(j)*n2_sum(j)     &
                           -n2_vec_sum(j)*n2_vec_sum(j))*(n3_sum(j)+   &
                           dlog(recn3(j))*recn3(j)*recn3(j))/          &
                           (12d00*pi*n3_sum(j)*n3_sum(j)*recn3(j)*     &
                           recn3(j))        
          d3(j)          = n3_sum(j)*(n3_sum(j)*n3_sum(j)-5.d00*       &
                           n3_sum(j)+2d00)+dlog(recn3(j))*2.d00*       &
                           recn3(j)*recn3(j)*recn3(j)                     
          dphidn3(j)     = n0_sum(j)/recn3(j)+(n1_sum(j)*n2_sum(j)-    &
                           n1_vec_sum(j)*n2_vec_sum(j))/(recn3(j)*     &
                           recn3(j))-(n2_sum(j)*n2_sum(j)*n2_sum(j)    &
                           -3d00*n2_sum(j)*n2_vec_sum(j)*              & 
                           n2_vec_sum(j))*d3(j)/(36d00*pi*n3_sum(j)*   &
                           n3_sum(j)*n3_sum(j)*recn3(j)*recn3(j)*      &
                           recn3(j))            
          dphidn1_vec(j) = -n2_vec_sum(j)/recn3(j)
          dphidn2_vec(j) = -n1_vec_sum(j)/recn3(j)-(1d00/6d00)*        &
                           n2_sum(j)*n2_vec_sum(j)*(n3_sum(j)+         &
                           dlog(recn3(j))*recn3(j)*recn3(j))/(pi*      &
                           recn3(j)*recn3(j)*n3_sum(j)*n3_sum(j))       
        enddo

        do i=1,nc
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain) 
           call fft(dphidn0,Fw0,c0ft,ext_bins,cplxn)
           call fft(dphidn1,Fw1,c1ft,ext_bins,cplxn)
           call fft(dphidn2,Fw2,c2ft,ext_bins,cplxn)
           call fft(dphidn3,Fw3,c3ft,ext_bins,cplxn)
           call fft(dphidn1_vec,-Fw1_vec,c1_vecft,ext_bins,cplxp) 
           call fft(dphidn2_vec,-Fw2_vec,c2_vecft,ext_bins,cplxp)
           k=1
           do ii=dif_bins+1,n_bins+dif_bins
              dF_hard_sphere(i,k) = m(i)*(c0ft(ii)+c1ft(ii)+c2ft(ii)+  &
                                    c3ft(ii)+c1_vecft(ii)+c2_vecft(ii))
              k=k+1
           enddo
        enddo
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine CLIJ(LIJ,C)
      implicit none
      integer*8 :: j,k
      real*8    :: OM(4,4),C(4), LIJ

      OM(1,1) = 0.81096D0
      OM(1,2) = 1.7888D0
      OM(1,3) = -37.578D0
      OM(1,4) = 92.284D0
      OM(2,1) = 1.0205D0
      OM(2,2) = -19.341D0
      OM(2,3) = 151.26D0
      OM(2,4) = -463.50D0
      OM(3,1) = -1.9057D0
      OM(3,2) = 22.845D0
      OM(3,3) = -228.14D0
      OM(3,4) = 973.92D0
      OM(4,1) = 1.0885D0
      OM(4,2) = -6.1962D0
      OM(4,3) = 106.98D0
      OM(4,4) = -677.64D0
      do j=1,4
         C(j) = OM(j,1) 
         do k=2,4
            C(j) = C(j)+OM(j,k)/LIJ**(k-1) 
         enddo
      enddo
     
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine cmix
      use globalvar  
      implicit none   
      integer*8 :: i,j,k
      real*8    :: OM(4,4),C(4),LIJa,LIJr
      
      do i=1,nc
         do j=1,nc
            LIJa = lambdaij(1,i,j)
            LIJr = lambdaij(2,i,j)
            c_ij(i,j)=LIJr*((LIJr/LIJa)**(LIJa/(LIJr-LIJa)))/(LIJr-LIJa)
         enddo
      enddo 
          
end subroutine       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE FUNMIE(EPSIJ,X,LRIJ,LAIJ,temp,FMIE)
      use globalvar
      implicit none
      real*8 ::  FMIE,U,EPSIJ,X,LRIJ,LAIJ,temp,UMIE

      temp=T
      U       = UMIE(EPSIJ,X,LRIJ,LAIJ)
      FMIE    = 1.D0-EXP(-U)

      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION alphaij(CIJ,LRIJ,LAIJ)
      implicit none
      real*8 :: alphaij,CIJ,LRIJ,LAIJ

      alphaij = CIJ*(1.d00/(LAIJ-3.D00)-1.d00/(LRIJ-3.D00))
      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE fun_alpha(CIJ,LRIJ,LAIJ,F)
      use globalvar
      real*8 :: CIJ,LRIJ,LAIJ,alpha,alphaij
      real*8 :: F(6),Faux(6),PHI(6,7)

      PHI(1,1) = 7.5365557D00
      PHI(1,2) = -37.60463D00
      PHI(1,3) = 71.745953D00
      PHI(1,4) = -46.83552D00
      PHI(1,5) = -2.467982D00
      PHI(1,6) = -0.50272D00
      PHI(1,7) = 8.0956883D00
      PHI(2,1) = -359.44D00
      PHI(2,2) = 1825.6D00
      PHI(2,3) = -3168.D00
      PHI(2,4) = 1884.2D00
      PHI(2,5) = -0.82376D00
      PHI(2,6) = -3.1935D00
      PHI(2,7) = 3.7090D00
      PHI(3,1) = 1550.9D00
      PHI(3,2) = -5070.1D00
      PHI(3,3) = 6534.6D00
      PHI(3,4) = -3288.7D00
      PHI(3,5) = -2.7171D00
      PHI(3,6) = 2.0883D00
      PHI(3,7) = 0.D00
      PHI(4,1) = -1.19932D00
      PHI(4,2) = 9.063632D00
      PHI(4,3) = -17.9482D00
      PHI(4,4) = 11.34027D00
      PHI(4,5) = 20.52142D00
      PHI(4,6) = -56.6377D00
      PHI(4,7) = 40.53683D00
      PHI(5,1) = -1911.28D00
      PHI(5,2) = 21390.175D00
      PHI(5,3) = -51320.7D00
      PHI(5,4) = 37064.54D00
      PHI(5,5) = 1103.742D00
      PHI(5,6) = -3264.61D00
      PHI(5,7) = 2556.181D00
      PHI(6,1) = 9236.90d00
      PHI(6,2) = -129430D00
      PHI(6,3) = 357230D00
      PHI(6,4) = -315530D00
      PHI(6,5) = 1390.2D00
      PHI(6,6) = -4518.2D00
      PHI(6,7) = 4241.6D00  
      F(:)=0.d00 
      Faux(:)=0.d00     
      alpha = alphaij(CIJ,LRIJ,LAIJ)
      do k=1,6
         do n=1,4
            Faux(k)=Faux(k)+phi(k,n)*alpha**(dble(n-1))           
         enddo
         do n=5,7         
            F(k)=F(k)+phi(k,n)*alpha**(dble(n-4))
         enddo 
         F(k)=Faux(k)/(1.d00+F(k))
      enddo 

end subroutine              


subroutine diam(sigij,epsij,lrij,laij,dkj)
      use globalvar
      implicit DOUBLE PRECISION (A-Z)
      integer*8   ::  k
      real*8      ::  W(5),X(5),FMIE,DKJ
      
      temp=T
      A    = 0.D0 
      B    = 1.D0
      W(1) = 0.2955242247D0
      W(2) = 0.2692667193D0
      W(3) = 0.2190863625D0
      W(4) = 0.1494513491D0
      W(5) = 0.0666713443D0
      X(1) = 0.1488743389D0
      X(2) = 0.4333953941D0
      X(3) = 0.6794095682D0
      X(4) = 0.8650633666D0
      X(5) = 0.9739065285D0

      XM = 0.5D0*(B+A)
      XR = 0.5D0*(B-A)
      DD = 0.D0
      do k=1,5
         DX = XR*X(k)
         CALL FUNMIE(EPSIJ,XM+DX,LRIJ,LAIJ,temp,FMIE)
         DP = FMIE
         CALL FUNMIE(EPSIJ,XM-DX,LRIJ,LAIJ,temp,FMIE)
         DM = FMIE
         DD = DD+W(K)*(DP+DM)
      enddo
      DKJ = XR*DD*SIGIJ
     
      end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION UMIE(EPSIJ,X,LRIJ,LAIJ)
      IMPLICIT DOUBLE PRECISION (A-Z)

      CIJ  = LRIJ*((LRIJ/LAIJ)**(LAIJ/(LRIJ-LAIJ)))/(LRIJ-LAIJ)
      UMIE = CIJ*EPSIJ*(X**(-LRIJ)-X**(-LAIJ))

      END
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine combrule
      use globalvar
      implicit none
      integer    :: i,j
      real*8     :: sigij,epsij,laij,lrij,dmie
      real*8     :: dkj(nc)
      real*8     :: sigmaij

      do i=1,nc
         sigij    = di(i)
         epsij    = epsilon_fluid_fluid(i)
         laij     = lambdaij(1,i,i)
         lrij     = lambdaij(2,i,i)
         call diam(sigij,epsij,lrij,laij,dkj(i))
      enddo
      do i=1,nc
         do j=1,nc
            dij(i,j)     = (dkj(i)+dkj(j))*0.5d00
            dij3(i,j)    = dij(i,j)*dij(i,j)*dij(i,j) 
            sigmaij      = (di(i)+di(j))*0.5d00
            sigmaij3(i,j)= sigmaij**3.d00  
            x0ij(i,j)    = sigmaij/dij(i,j)             
            eps(i,j)     = dsqrt(di(i)**3.d00*di(j)**3.d00*            & 
                           epsilon_fluid_fluid(i)*                     & 
                           epsilon_fluid_fluid(j))/sigmaij**3.d00
         enddo
      enddo    
      end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                      
	subroutine a1scalc(rhos,Zetaeff,dZetaeff,mon,ep,dia3,lij,a1s,da1s) 
	use globalvar
	implicit none
	real*8 :: Zetaeff,rhos,ep,lij,dia3,mon,t3,t4,dZetaeff,a1s,da1s 
        		
        
	t3         =   (1.d00-Zetaeff*0.5d00)/(1.d00-         &
                               Zetaeff)**3.d00
	t4         =   0.5d00*(5.d00-2.d00*Zetaeff)/          &
                                (1.d00-Zetaeff)**4.d00	
	a1s        =  -2.d00*rhos*pi*ep*    &
                      dia3*(1.d00-Zetaeff*0.5d00)/    &
                      ((1.d00-Zetaeff)**3.d00*               &
                      (lij-3.d00))
	da1s       =  -2.d00*mon*t3*pi*ep*       &
                       dia3/(lij-3.d00)- &
                       2.d00*rhos*pi*ep*    &
                       dia3/(lij-3.d00)  &
                       *t4*dZetaeff
               
	end subroutine 	

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine Bijcalc(x0,lij,rhos,Zeta,dZeta,mon,ep,dia3,Bij,dBij)
	use globalvar
	implicit none
	real*8 :: t1,t2,dt1,dt2,Zeta,dZeta,Iij,Jij,Bij,dBij,x0,lij,ep,dia3,rhos,mon

                 t1  = (1.d00-Zeta*0.5d00)/(1.d00-Zeta)**3.d00
                 t2  = 4.5d00*Zeta*(1.d00+Zeta)/(1.d00-Zeta)**3.d00
                 dt1 = 0.5d00*(5.d00-2.d00*Zeta)/(1.d00-Zeta)**4.d00
                 dt2 = 0.5d00*(9.d00+9.d00*Zeta*Zeta+36.d00*Zeta)/     &
                       (1.d00-Zeta)**4.d00		
                 Iij        =  -(x0**(3.d00-lij)-1.d00)/(lij-3.d00)
                 Jij        =  -(x0**(4.d00-lij)*(lij-3.d00)-         &
                                x0**(3.d00-lij)*(lij-4.d00)-1.d00)/    &
                                ((lij-3.d00)*(lij-4.d00)) 
                 Bij        =   2.d00*pi*ep*dia3*       &
                                rhos*(Iij*t1-Jij*t2)
                 dBij       =   2.d00*pi*dia3*ep*       &
                                ((Iij*t1-Jij*t2)*mon+          &
                                rhos*(Iij*dt1*dZeta-      &
                                Jij*dt2*dZeta)) 
                                 
	end subroutine 	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            
	subroutine a1scalc_chain(rhos,Zetaeff,dZetaeff,dZetaeffrhos,ddZetaeffrhos,mon,ep,dia3,lij,a1s,da1s,da1srhos,dda1srhos) 
	use globalvar
	implicit none
	real*8 :: Zetaeff,rhos,ep,lij,dia3,mon,t3,t4,dZetaeff,a1s,da1s,dt3,dt4 
        real*8 :: dt3rhos,dt4rhos,dZetarhos,ddZetarhos,da1srhos,dda1srhos
        real*8 :: ddZetaeffrhos,dZetaeffrhos        		
        
	t3         =   (1.d00-Zetaeff*0.5d00)/(1.d00-         &
                               Zetaeff)**3.d00
	t4         =   0.5d00*(5.d00-2.d00*Zetaeff)/          &
                                (1.d00-Zetaeff)**4.d00	
	a1s        =  -2.d00*rhos*pi*ep*    &
                      dia3*(1.d00-Zetaeff*0.5d00)/    &
                      ((1.d00-Zetaeff)**3.d00*               &
                      (lij-3.d00))
	da1s       =  -2.d00*mon*t3*pi*ep*       &
                       dia3/(lij-3.d00)- &
                       2.d00*rhos*pi*ep*    &
                       dia3/(lij-3.d00)  &
                       *t4*dZetaeff
        dt3        =  t4*dZetaeff
        dt4        =  0.5d00*(-6.d00*Zetaeff+18.d00)/((1.d00-          &
                      Zetaeff)**5.d00)
        dt4        =  dt4*dZetaeff
      
        dt3rhos    =  t4*dZetaeffrhos
        dt4rhos    =  0.5d00*(-6.d00*Zetaeff+18.d00)/((1.d00-          &
                      Zetaeff)**5.d00)
        dt4rhos    =  dt4rhos*dZetaeffrhos
        da1srhos   =  -2.d00*pi*ep*dia3/(lij-3.d00)*(t3+rhos*t4        &
                      *dZetaeffrhos)
	dda1srhos  =  -2.d00*pi*ep*dia3/(lij-3.d00)*(mon*t4*dZetaeffrhos+      &
                      t4*dZetaeff+rhos*(dZetaeffrhos*dt4+t4*           &
                      ddZetaeffrhos))               
	end subroutine 	

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine Bijcalc_chain(x0,lij,rhos,Zeta,dZeta,dZetarhos,ddZetarhos,mon,ep,dia3,Bij,dBij,dBijrhos,ddBijrhos)
	use globalvar
	implicit none
	real*8 :: t1,t2,dt1,dt2,Zeta,dZeta,Iij,Jij,Bij,dBij,x0,lij,ep,dia3,rhos,mon
        real*8 :: dt1rhos,dt2rhos,ddt1rhos,ddt2rhos
        real*8 :: dBijrhos,ddBijrhos,dZetarhos,ddZetarhos

                 t1  = (1.d00-Zeta*0.5d00)/(1.d00-Zeta)**3.d00
                 t2  = 4.5d00*Zeta*(1.d00+Zeta)/(1.d00-Zeta)**3.d00
                 dt1 = 0.5d00*(5.d00-2.d00*Zeta)/(1.d00-Zeta)**4.d00
                 dt2 = 0.5d00*(9.d00+9.d00*Zeta*Zeta+36.d00*Zeta)/     &
                       (1.d00-Zeta)**4.d00		
                 Iij        =  -(x0**(3.d00-lij)-1.d00)/(lij-3.d00)
                 Jij        =  -(x0**(4.d00-lij)*(lij-3.d00)-         &
                                x0**(3.d00-lij)*(lij-4.d00)-1.d00)/    &
                                ((lij-3.d00)*(lij-4.d00)) 
                 Bij        =   2.d00*pi*ep*dia3*       &
                                rhos*(Iij*t1-Jij*t2)
                 dBij       =   2.d00*pi*dia3*ep*       &
                                ((Iij*t1-Jij*t2)*mon+          &
                                rhos*(Iij*dt1*dZeta-      &
                                Jij*dt2*dZeta)) 
                 dt1rhos = dZetarhos*dt1 
                 dt2rhos = dZetarhos*dt2
                 ddt1rhos = ddZetarhos*dt1+dZeta*dZetarhos*0.5d00*     &
                            (-6.d00*Zeta+18.d00)/((1.d0-Zeta)**5.d00)                    
                 ddt2rhos = ddZetarhos*dt2+dZeta*dZetarhos*(63.d00*     &
                            Zeta+9.d00*Zeta*Zeta+36.d00)/(1.d00-Zeta)**5.d00                            
                 dBijrhos  = 2.d00*pi*ep*dia3*(Iij*t1-Jij*t2+rhos*(Iij &
                             *dt1rhos-Jij*dt2rhos))
                 ddBijrhos = 2.d00*pi*ep*dia3*(mon*(Iij*dt1rhos-Jij*   &
                             dt2rhos)+(Iij*dt1*dZeta-Jij*dt2*dZeta)+   &
                             rhos*(Iij*ddt1rhos-Jij*ddt2rhos))  
end subroutine 	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine monomers (ct,dF_At)
	use globalvar	
	implicit none
        integer*8, parameter  :: cplxp = 1
        integer*8, parameter  :: cplxn = 0 
        integer   ::   ct 
        integer*8 ::   i,j,k,ii,jj,kk,jjj
        integer*8 ::   z_index,dif_bins
        real*8    ::   Abig(nc,ext_bins)
        real*8    ::   d_A1(nc,ext_bins),d_A2(nc,ext_bins),d_A3(nc,ext_bins)
        real*8    ::   d_A(nc,ext_bins),dF_A(nc,ext_bins)       
        real*8    ::   dF_At(nc,n_bins),da1(nc,ext_bins),da2(nc,ext_bins),da3(nc,ext_bins)
        real*8    ::   rhoseg(nc,ext_bins)
        real*8    ::   rhoseg_orig(nc,ext_bins)
        real*8    ::   rho_ft(nc,ext_bins),temp(nc,ext_bins),rho_orig(nc,ext_bins)
        real*8    ::   chempottest(nc,ext_bins)        
        real*8    ::   rhosegtotal(ext_bins),rhosegtotal_orig(ext_bins)
        real*8    ::   rhototal(ext_bins),rhototal_orig(ext_bins)
        real*8    ::   x(nc,ext_bins)
        real*8    ::   xseg(nc,ext_bins)
        real*8    ::   sumxseg,sumxseg_a2
        real*8    ::   d1(nc,ext_bins),a1(ext_bins)
        real*8    ::   d1_a2(nc,ext_bins),a2(ext_bins)  
        real*8    ::   d1_a3(nc,ext_bins),a3(ext_bins)      
        real*8    ::   Zeta,Zeta2,Zeta3,Zeta4,dZeta,t1,t2,dt1,dt2
        real*8    ::   a1ij,da1ij,a2ij,da2ij,a2ijtemp,a3ij,da3ij,a3ijtemp
        real*8    ::   p1(nc),p1_a2(nc),d2(nc),d2_a2(nc),d2_a3(nc)
        real*8    ::   Zetaeff,dZetaeff,Zetaavg,dZetaavg
        real*8    ::   a1s,da1s,Bij,dBij
        real*8    ::   K_hs,dK_hs,Xiij,dXiij,ep,dia3,lij,C(4)
        real*8, dimension (ext_bins) :: Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec
        real*8, dimension (ext_bins) :: Fdisp,Fchain 


        rho_ft(:,:) = 0.d00
        dif_bins  = dint((ext_bins-n_bins)/2.d00)
        if(ct.eq.1)then
          do i=1,nc
	     z_index   = 1           
             do j=dif_bins+1,n_bins+dif_bins
                rho_ft(i,j) = density_calc(i,z_index)
                z_index   = z_index+1
             enddo 
          enddo
          z_index   = 1
        else
          do i=1,nc
             rho_ft(i,1:ext_bins)=bulk_density_i(i)
          enddo
        endif
             
        do i=1,nc
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain)
           call fft(rho_ft(i,:),Fdisp,temp(i,:),ext_bins,cplxn)
           rho_orig(i,:)=rho_ft(i,:)        
           rho_ft(i,:)=temp(i,:)
        enddo

        rhosegtotal(:) = 0.d0
        rhoseg(:,:) = 0.d0
        rhototal(:) = 0.d0
        rhosegtotal_orig(:) = 0.d0
        rhoseg_orig(:,:) = 0.d0
        rhototal_orig(:) = 0.d0

        do i=1,nc
	   do j=1,ext_bins
	      rhoseg(i,j)    = rho_ft(i,j)*m(i)
              rhosegtotal(j) = rhosegtotal(j)+rho_ft(i,j)*m(i)
              rhototal(j)    = rhototal(j)+rho_ft(i,j)
	      rhoseg_orig(i,j)    = rho_orig(i,j)*m(i)
              rhosegtotal_orig(j) = rhosegtotal_orig(j)+rho_orig(i,j)*m(i)
              rhototal_orig(j)    = rhototal_orig(j)+rho_orig(i,j)
           enddo
        enddo
        
        d1(:,:)=0.d00   
        d2(:)=0.d00  
        a1(:)=0.d00
        d1_a2(:,:)=0.d00   
        d2_a2(:)=0.d00  
        a2(:)=0.d00  
        d1_a3(:,:)=0.d00   
        d2_a3(:)=0.d00  
        a3(:)=0.d00         
        do j = 1,ext_bins
           p1(:)=0.d0
           p1_a2(:)=0.d00
           sumxseg=0.d00
           sumxseg_a2=0.d00
           if(rhosegtotal(j).lt.1d-20)then
              rhosegtotal(j) = 1d-20 
              rhototal(j) = 1d-20
           endif    
           do i = 1,nc
              xseg(i,j) = rhoseg(i,j)/rhosegtotal(j)
              x(i,j)    = rho_ft(i,j)/rhototal(j)
           enddo
           do i=1,nc             
              do ii =1,nc
                 sumxseg = sumxseg + xseg(i,j)*xseg(ii,j)*dij3(i,ii)
                 p1(i)      =   p1(i)+xseg(ii,j)*dij3(i,ii)
                 sumxseg_a2 = sumxseg_a2 + xseg(i,j)*xseg(ii,j)*sigmaij3(i,ii)
                 p1_a2(i)      =   p1_a2(i)+xseg(ii,j)*sigmaij3(i,ii)             
              enddo
           enddo

           Zeta = pi*(1.d00/6.d00)*sumxseg*rhosegtotal(j)
           Zetaavg = pi*(1.d00/6.d00)*sumxseg_a2*rhosegtotal(j) 
           Zeta2 = Zeta*Zeta
           Zeta3 = Zeta2*Zeta
           Zeta4 = Zeta3*Zeta
           K_hs = 1.d00+4.d00*Zeta+4.d00*Zeta**2.d00-            &
                  4.d00*Zeta**3.d00+Zeta**4.d00
           K_hs = (1.d00-Zeta)**4.d00/K_hs
           d2(:)=0.d00
           d2_a2(:)=0.d00
           d2_a3(:)=0.d00           
           do jjj = 1,nc
                 dZeta      =   m(jjj)*pi*1.d00/6.d00*(2.d00*p1(jjj)-    &
                                sumxseg) 
                 dZetaavg   =   m(jjj)*pi*1.d00/6.d00*(2.d00*p1_a2(jjj)- &
                                sumxseg_a2)
                 dK_hs      =   -4.d00*K_hs/(1.d00-Zeta)*dZeta-K_hs**  &
                                2.d00/(1.d00-Zeta)**4.d00*dZeta*       &
                                (4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*  &
                                Zeta3)
           a1(j)=0.d00
           a2(j)=0.d00
           a3(j)=0.d00                      
           d2(:)=0.d00
           d2_a2(:)=0.d00 
           d2_a3(:)=0.d00                     
           do jj=1,nc
                                                                                                                                
              do kk =1,nc
                 da1ij=0.d0
                 ep = eps(jj,kk)
                 dia3 = dij3(jj,kk)
                 Xiij       =   F_alpha(jj,kk,1)*Zetaavg+F_alpha(jj,kk,2)* &
                                Zetaavg**5.d00+F_alpha(jj,kk,3)*         &
                                Zetaavg**8.d00
                 dXiij      =   F_alpha(jj,kk,1)*dZetaavg+5.d00*       &
                                F_alpha(jj,kk,2)*Zetaavg**4.d00*       &
                                dZetaavg+8.d00*F_alpha(jj,kk,3)*       &
                                Zetaavg**7.d00*dZetaavg                                
                 lij = lambdaij(1,jj,kk)
                 call CLIJ(lij,C)                 
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 call a1scalc(rhosegtotal(j),Zetaeff,dZetaeff,m(jjj),   &
                      ep,dia3,lambdaij(1,jj,kk),a1s,da1s)                                
                 call Bijcalc(x0ij(jj,kk),lambdaij(1,jj,kk),           &
                      rhosegtotal(j),Zeta,dZeta,m(jjj),ep,dia3,Bij,dBij) 
                 a1ij       =   x0ij(jj,kk)**lij*(a1s+Bij)             
                 da1ij      =   x0ij(jj,kk)**lij*(da1s+dBij) 
                 
                 lij = lambdaij(2,jj,kk)
                 call CLIJ(lij,C)                 
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 call a1scalc(rhosegtotal(j),Zetaeff,dZetaeff,m(jjj),ep,dia3,lambdaij(2,jj,kk),a1s,da1s)                                
                 call Bijcalc(x0ij(jj,kk),lambdaij(2,jj,kk),rhosegtotal(j),Zeta,dZeta,m(jjj),ep,dia3,Bij,dBij) 
                 a1ij       =   c_ij(jj,kk)*(a1ij-x0ij(jj,kk)**lij*(a1s+Bij))             
                 da1ij      =   c_ij(jj,kk)*(da1ij-x0ij(jj,kk)**lij*(da1s+dBij))                  
                 a1(j)      =   a1(j)+xseg(jj,j)*xseg(kk,j)*a1ij
                 d1(jjj,j)      =   d1(jjj,j)+xseg(jj,j)*xseg(kk,j)*da1ij
                 d2(jj)     =   d2(jj)+xseg(kk,j)*a1ij 
                   

                 lij = 2.d00*lambdaij(1,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 call a1scalc(rhosegtotal(j),Zetaeff,dZetaeff,m(jjj),ep,dia3,lij,a1s,da1s)                                
                 call Bijcalc(x0ij(jj,kk),lij,rhosegtotal(j),Zeta,dZeta,m(jjj),ep,dia3,Bij,dBij) 
                 a2ij       =   x0ij(jj,kk)**lij*(a1s+Bij)             
                 da2ij      =   x0ij(jj,kk)**lij*(da1s+dBij)
                 
                 lij = 2.d00*lambdaij(2,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 call a1scalc(rhosegtotal(j),Zetaeff,dZetaeff,m(jjj),ep,dia3,lij,a1s,da1s)                                
                 call Bijcalc(x0ij(jj,kk),lij,rhosegtotal(j),Zeta,dZeta,m(jjj),ep,dia3,Bij,dBij) 
                 a2ij       =   a2ij+x0ij(jj,kk)**lij*(a1s+Bij)             
                 da2ij      =   da2ij+x0ij(jj,kk)**lij*(da1s+dBij)                  

                 lij = lambdaij(1,jj,kk)+lambdaij(2,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 call a1scalc(rhosegtotal(j),Zetaeff,dZetaeff,m(jjj),ep,dia3,lij,a1s,da1s)                                
                 call Bijcalc(x0ij(jj,kk),lij,rhosegtotal(j),Zeta,dZeta,m(jjj),ep,dia3,Bij,dBij) 
                 a2ijtemp   =   a2ij-2.d00*x0ij(jj,kk)**lij*(a1s+Bij)             
                 a2ij       =   0.5d00*K_hs*(1.d00+Xiij)*ep*c_ij(jj,kk)**2.d00*a2ijtemp
                 da2ij      =   da2ij-2.d00*x0ij(jj,kk)**lij*(da1s+dBij)
                 da2ij      =   da2ij*K_hs*(1.d00+Xiij)+dK_hs*(1.d00+Xiij)*a2ijtemp+K_hs*dXiij*a2ijtemp
                 da2ij      =   (0.5d00*ep*c_ij(jj,kk)**2.d00)*da2ij
                                 
                 a2(j)      =   a2(j)+xseg(jj,j)*xseg(kk,j)*a2ij
                 d1_a2(jjj,j)   =   d1_a2(jjj,j)+xseg(jj,j)*xseg(kk,j)*da2ij
                 d2_a2(jj)  =   d2_a2(jj)+xseg(kk,j)*a2ij

                 a3ij       =   -ep**3.d00*F_alpha(jj,kk,4)*Zetaavg*   &
                                dexp(F_alpha(jj,kk,5)*Zetaavg+         &
                                F_alpha(jj,kk,6)*Zetaavg**2.d00)
                 a3ijtemp   =   -ep**3.d00*F_alpha(jj,kk,4)*           &
                                dexp(F_alpha(jj,kk,5)*Zetaavg+         &
                                F_alpha(jj,kk,6)*Zetaavg**2.d00)
                 da3ij      =   dZetaavg*(a3ij*(F_alpha(jj,kk,5)+2.d00*F_alpha(jj,kk,6)*Zetaavg)+a3ijtemp)
                 a3(j)      =   a3(j)+xseg(jj,j)*xseg(kk,j)*a3ij
                 d1_a3(jjj,j)   =   d1_a3(jjj,j)+xseg(jj,j)*xseg(kk,j)*da3ij
                 d2_a3(jj)  =   d2_a3(jj)+xseg(kk,j)*a3ij   
                 
              enddo
           enddo
           enddo

           do ii=1,nc
                 da1(ii,j)    =   m(ii)*(2.d00*d2(ii)-2.d00*a1(j)+     &
                                  rhosegtotal(j)*d1(ii,j)/m(ii))/         &
                                  rhosegtotal(j)
                 da2(ii,j)    =   m(ii)*(2.d00*d2_a2(ii)-2.d00*a2(j)+     &
                                  rhosegtotal(j)*d1_a2(ii,j)/m(ii))/         &
                                  rhosegtotal(j) 
                 da3(ii,j)    =   m(ii)*(2.d00*d2_a3(ii)-2.d00*a3(j)+     &
                                  rhosegtotal(j)*d1_a3(ii,j)/m(ii))/         &
                                  rhosegtotal(j)                                                     
                 da1(ii,j)    =   (rhosegtotal(j)*da1(ii,j)+a1(j)*     &
                                  (m(ii)-rhosegtotal(j)/rhototal(j)))  &
                                  *1.d00/rhototal(j)   
                 da2(ii,j)    =   (rhosegtotal(j)*da2(ii,j)+a2(j)*     &
                                  (m(ii)-rhosegtotal(j)/rhototal(j)))  &
                                  *1.d00/rhototal(j) 
                 da3(ii,j)    =   (rhosegtotal(j)*da3(ii,j)+a3(j)*     &
                                  (m(ii)-rhosegtotal(j)/rhototal(j)))  &
                                  *1.d00/rhototal(j)                               
                 Abig(ii,j) =    (a1(j)+a2(j)+a3(j))*rhosegtotal(j)/rhototal(j)

                 d_A1(ii,j)  =    a1(j)*rhosegtotal(j)/rhototal(j)+    &
                                  rhototal(j)*da1(ii,j) 
                 d_A2(ii,j)  =    a2(j)*rhosegtotal(j)/rhototal(j)+    &
                                  rhototal(j)*da2(ii,j) 
                 d_A3(ii,j)  =    a3(j)*rhosegtotal(j)/rhototal(j)+    &
                                  rhototal(j)*da3(ii,j) 

                 d_A(ii,j)   =    d_A1(ii,j)+d_A2(ii,j)+d_A3(ii,j) 
           enddo
        enddo

        do i=1,nc
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain) 
           call fft(d_A(i,1:ext_bins),Fdisp,dF_A(i,1:ext_bins),ext_bins,cplxn)
           k=1             
           do ii=dif_bins+1,n_bins+dif_bins
              dF_At(i,k) = dF_A(i,ii)               
              k=k+1   
           enddo
        enddo
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine chain (ct,dF_chain)
	use globalvar	
	implicit none
        integer*8, parameter  :: cplxp = 1
        integer*8, parameter  :: cplxn = 0 
        integer   ::   ct 
        integer*8 ::   i,j,k,ii,jj,kk,jjj
        integer*8 ::   z_index,dif_bins
        real*8    ::   Csi(nc,ext_bins),dCsin2(nc,ext_bins),dCsin2vec(nc,ext_bins)
        real*8    ::   n2(nc,ext_bins),n2_vec(nc,ext_bins),n0(nc,ext_bins)
        real*8    ::   achain1(nc,ext_bins),achain2(nc,ext_bins)
        real*8    ::   achain3(nc,ext_bins),achain4(nc,ext_bins,nc)
        real*8    ::   dF_chain1(nc,ext_bins),dF_chain2(nc,ext_bins)
        real*8    ::   dF_chain3(nc,ext_bins),dF_chain4(nc,ext_bins)
        real*8    ::   dF_chain(nc,n_bins),chain4(ext_bins)
        real*8    ::   rhoseg(nc,ext_bins)
        real*8    ::   rhoseg_orig(nc,ext_bins)
        real*8    ::   rho_ft(nc,ext_bins),temp(nc,ext_bins),rho_orig(nc,ext_bins)
        real*8    ::   chempottest(nc,ext_bins)        
        real*8    ::   rhosegtotal(ext_bins),rhosegtotal_orig(ext_bins)
        real*8    ::   rhototal(ext_bins),rhototal_orig(ext_bins)
        real*8    ::   x(nc,ext_bins)
        real*8    ::   xseg(nc,ext_bins)
        real*8    ::   sumxseg,sumxseg_a2
        real*8    ::   k0,k1,k2,k3,dk0,dk1,dk2,dk3,dBK_hs
        real*8    ::   dZetarhos,dK_hsrhos,dXiijrhos,dZetaeffrhos,dZetaavgrhos
        real*8    ::   ddZetarhos,ddK_hsrhos,ddK1,ddK2,ddXiijrhos,ddZetaeffrhos,ddZetaavgrhos
        real*8    ::   da1ijrhos,da2ijrhos,a2ijchain,a2ijtemp
        real*8    ::   dda1ijrhos,dda2ijrhos,da2ijchain,da2ijtemp
        real*8    ::   dTa2ijrhos,ddTa2ijrhos,da2byXi,dda2byXi
        real*8    ::   g1ijMie,dg1ijMie,g2ijMie,dg2ijMie,g2ijMCA,dg2ijMCA
        real*8    ::   thetaij,alpha,gamacij,dgamacij
        real*8    ::   gijHS,dgijHS,dg1,dg2
        real*8    ::   gijMie(nc,nc),dgijMie(nc,nc)
        real*8    ::   Zeta,Zetainv,Zeta2,Zeta3,Zeta4,dZeta,t1,t2,dt1,dt2
        real*8    ::   a1ij,da1ij,a2ij,da2ij,total_bulk_dens,total_bulk_dens_seg
        real*8    ::   p1(nc),p1_a2(nc),d2(nc),d2_a2(nc),d2_a3(nc)
        real*8    ::   Zetaeff,dZetaeff,Zetaavg,dZetaavg,invXi
        real*8    ::   a1s,da1s,Bij,dBij,da1srhos,dBijrhos,dda1srhos,ddBijrhos
        real*8    ::   K_hs,dK_hs,Xiij,dXiij,ep,dia3,lij,alphaij,C(4)
        real*8, dimension (ext_bins) :: Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec
        real*8, dimension (ext_bins) :: Fdisp,Fchain 


        rho_ft(:,:) = 0.d00
        dif_bins  = dint((ext_bins-n_bins)/2.d00)
        total_bulk_dens=0.d00
        total_bulk_dens_seg=0.d00
        if(ct.eq.1)then
          do i=1,nc
	     z_index   = 1           
             do j=dif_bins+1,n_bins+dif_bins
                rho_ft(i,j) = density_calc(i,z_index)
                z_index   = z_index+1
             enddo 
          enddo
          z_index   = 1
        else
          do i=1,nc
             rho_ft(i,1:ext_bins)=bulk_density_i(i)
             total_bulk_dens=total_bulk_dens+bulk_density_i(i)
             total_bulk_dens_seg=total_bulk_dens_seg+m(i)*bulk_density_i(i)
          enddo
        endif
          
        do i=1,nc
           rho_orig(i,:)=rho_ft(i,:)  
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain)
           call fft(rho_ft(i,:),Fchain,temp(i,:),ext_bins,cplxn)
           rho_ft(i,:)=temp(i,:)
         
           do j=1,ext_bins
	      if(rho_ft(i,j).lt.1d-17)then
	         rho_ft(i,j)=1d-20
	      endif           
           enddo                             
        enddo
        
        rhosegtotal(:) = 0.d0
        rhoseg(:,:) = 0.d0
        rhototal(:) = 0.d0
        rhosegtotal_orig(:) = 0.d0
        rhoseg_orig(:,:) = 0.d0
        rhototal_orig(:) = 0.d0

        do i=1,nc
	   do j=1,ext_bins
	      rhoseg(i,j)    = rho_ft(i,j)*m(i)
              rhosegtotal(j) = rhosegtotal(j)+rho_ft(i,j)*m(i)
              rhototal(j)    = rhototal(j)+rho_ft(i,j)
	      rhoseg_orig(i,j)    = rho_orig(i,j)*m(i)
              rhosegtotal_orig(j) = rhosegtotal_orig(j)+rho_orig(i,j)*m(i)
              rhototal_orig(j)    = rhototal_orig(j)+rho_orig(i,j)
           enddo
        enddo
   
        do j = 1,ext_bins
           p1(:)=0.d0
           p1_a2(:)=0.d00
           sumxseg=0.d00
           sumxseg_a2=0.d00
           if(rhosegtotal(j).lt.1d-17)then
              rhosegtotal(j) = 1d-20 
              rhototal(j) = 1d-20
           endif    
           do i = 1,nc
              xseg(i,j) = rhoseg(i,j)/rhosegtotal(j)
              x(i,j)    = rho_ft(i,j)/rhototal(j)
              if(ct.eq.0)then
                xseg(i,j) = m(i)*bulk_density_i(i)/total_bulk_dens_seg
                x(i,j)    = bulk_density_i(i)/total_bulk_dens 
              endif             
           enddo             

           do i=1,nc             
              do ii =1,nc
                 sumxseg = sumxseg + xseg(i,j)*xseg(ii,j)*dij3(i,ii)
                 p1(i)      =   p1(i)+xseg(ii,j)*dij3(i,ii)
                 sumxseg_a2 = sumxseg_a2 + xseg(i,j)*xseg(ii,j)*sigmaij3(i,ii)
                 p1_a2(i)      =   p1_a2(i)+xseg(ii,j)*sigmaij3(i,ii)                 
              enddo
           enddo
           
           if(rhosegtotal(j).lt.1d-20)then
              sumxseg= 1d-20 
              p1 = 1d-20
              sumxseg_a2= 1d-20 
              p1_a2 = 1d-20

           endif             
                    
           Zeta = pi*(1.d00/6.d00)*sumxseg*rhosegtotal(j)
           Zetaavg = pi*(1.d00/6.d00)*sumxseg_a2*rhosegtotal(j) 
           Zeta2 = Zeta*Zeta
           Zeta3 = Zeta2*Zeta
           Zeta4 = Zeta3*Zeta
           Zetainv = 1.d00/(1.d00-Zeta)
           dZetarhos  =   pi*1.d00/6.d00*sumxseg  !ATTENTION: SINGLE COMPONENT           
           dZetaavgrhos = pi*1.d00/6.d00*sumxseg_a2 !ATTENTION: SINGLE COMPONENT           
           K_hs = 1.d00+4.d00*Zeta+4.d00*Zeta**2.d00-                  &
                  4.d00*Zeta**3.d00+Zeta**4.d00
           K_hs = (1.d00-Zeta)**4.d00/K_hs
           dK_hsrhos  =   -4.d00*K_hs/(1.d00-Zeta)*dZetarhos-K_hs**  &
                          2.d00/(1.d00-Zeta)**4.d00*dZetarhos*       &
                          (4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*      &
                          Zeta3)
           ddK1       =   1.d00+4.d00*Zeta+4.d00*Zeta2-4.d00*Zeta3+  &
                          Zeta4  
           ddK2       =   4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*Zeta3

           k0  =  -dlog(1.d00-Zeta)+(42.d00*Zeta-39.d00*Zeta2+         &
                  9.d00*Zeta3-2.d00*Zeta4)*Zetainv**3.d00/6.d00 
           k1  =  (Zeta4+6.d00*Zeta2-12.d00*Zeta)*Zetainv**3.d00*0.5d00 
           k2  =  -3.d00*Zeta2*0.125d00*Zetainv**2.d00 
           k3  =  (-Zeta4+3.d00*Zeta2+3.d00*Zeta)*Zetainv**3.d00/6.d00
        
           do jjj = 1,nc
                 dZeta      =   m(jjj)*pi*1.d00/6.d00*(2.d00*p1(jjj)-    &
                                sumxseg) 
                 dZetaavg   =   m(jjj)*pi*1.d00/6.d00*(2.d00*p1_a2(jjj)- &
                                sumxseg_a2)
                 dK_hs      =   -4.d00*K_hs/(1.d00-Zeta)*dZeta-K_hs**  &
                                2.d00/(1.d00-Zeta)**4.d00*dZeta*       &
                                (4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*  &
                                Zeta3)
                                                         
           dk0 =  (Zetainv+(21.d00+3.d00*Zeta-6.d00*Zeta2-             &
                  4.d00*Zeta3+Zeta4)*Zetainv**4.d00/3.d00)*dZeta   
           dk1 =  (-Zeta4+4.d00*Zeta3+6.d00*Zeta2-12.d00*Zeta-         &
                  12.d00)*0.5d00*Zetainv**4.d00*dZeta
           dk2 =  -6.d00*Zeta*0.125d00*Zetainv**3.d00*dZeta 
           dk3 =  (Zeta4-4.d00*Zeta3+3.d00*Zeta2+12.d00*Zeta+3.d00)    &
                  *Zetainv**4.d00*(1.d00/6.d00)*dZeta                                 
! Chain terms
                 ddZetarhos =   -m(jjj)*pi*1.d00/3.d00*(sumxseg-p1(jjj))/ &
                                rhosegtotal(j)
                 ddZetaavgrhos = -m(jjj)*pi*1.d00/3.d00*(sumxseg_a2-p1_a2(jjj))/ &
                                rhosegtotal(j)

                 ddK_hsrhos =   ddZetarhos*(-4.d00*K_hs/(1.d00-Zeta)-K_hs**  &
                                2.d00/(1.d00-Zeta)**4.d00*                   &
                                (4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*        &
                                Zeta3)) 
                 dBK_hs     =   (4.d00+8.d00*Zeta-12.d00*Zeta2+4.d00*        &
                                Zeta3)*(2.d00*K_hs/(1.d00-Zeta)**4.d00* &
                                dK_hs+4.d00*K_hs**2.d00/(1.d00-         &
                                Zeta)**5.d00*dZeta)+K_hs**2.d00/        &
                                (1.d00-Zeta)**4.d00*dZeta*(8.d00-24.d0* &
                                Zeta+12.d0*Zeta2)  
                 ddK_hsrhos =   ddK_hsrhos+dZetarhos*(-4.d00*dK_hs      &
                                /(1.d00-Zeta)-4.d00*K_hs/(1.d00-Zeta)   &
                                **2.d00*dZeta-dBK_hs) 

              da2ijrhos   = 0.d00 
              dda2byXi    = 0.d00
              da2ijchain  = 0.d00
              dda2ijrhos  = 0.d00
              dTa2ijrhos  = 0.d00
              ddTa2ijrhos = 0.d00
              dda1srhos   = 0.d00
              ddBijrhos   = 0.d00
              dg1ijMie    = 0.d00
              dg2ijMie    = 0.d00
              dgijMie(:,:)= 0.d00
              
            do jj=1,nc     
              do kk =1,nc
                 ep = eps(jj,kk)
                 dia3 = dij3(jj,kk)
                 Xiij       =   F_alpha(jj,kk,1)*Zetaavg+F_alpha(jj,kk,2)* &
                                Zetaavg**5.d00+F_alpha(jj,kk,3)*         &
                                Zetaavg**8.d00
                 dXiij      =   F_alpha(jj,kk,1)*dZetaavg+5.d00*       &
                                F_alpha(jj,kk,2)*Zetaavg**4.d00*       &
                                dZetaavg+8.d00*F_alpha(jj,kk,3)*       &
                                Zetaavg**7.d00*dZetaavg                               
                 dXiijrhos  =   (F_alpha(jj,kk,1)+5.d00*               &
                                F_alpha(jj,kk,2)*Zetaavg**4.d00+       &
                                8.d00*F_alpha(jj,kk,3)*                &
                                Zetaavg**7.d00)*dZetaavgrhos 
                 ddXiijrhos =   F_alpha(jj,kk,1)*ddZetaavgrhos+5.d00*  &
                                F_alpha(jj,kk,2)*(4.d00*Zetaavg**3.d00 &
                                *dZetaavg*dZetaavgrhos+Zetaavg**4.d00  &
                                *ddZetaavgrhos)+8.d00*                 &
                                F_alpha(jj,kk,3)*(7.d00*Zetaavg**6.d00 &
                                *dZetaavg*dZetaavgrhos+Zetaavg**7.d00* &
                                ddZetaavgrhos)                                  
                 lij = lambdaij(1,jj,kk)
                 call CLIJ(lij,C)                 
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 ddZetaeffrhos  = (C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*ddZetarhos+dZeta*     &
                                (2.d00*C(2)+6.d00*C(3)*                &
                                Zeta+12.d00*C(4)*                      & 
                                Zeta2)*dZetarhos
                 dZetaeffrhos   = (C(1)+2.d00*C(2)                     &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*dZetarhos          


                 call a1scalc_chain(rhosegtotal(j),Zetaeff,dZetaeff,   &
                      dZetaeffrhos,ddZetaeffrhos,m(jjj),ep,dia3,        &
                      lambdaij(1,jj,kk),a1s,da1s,da1srhos,dda1srhos)                                
                 call Bijcalc_chain(x0ij(jj,kk),lambdaij(1,jj,kk),     &
                      rhosegtotal(j),Zeta,dZeta,dZetarhos,ddZetarhos,  &
                      m(jjj),ep,dia3,Bij,dBij,dBijrhos,ddBijrhos)                     
                 a1ij       =   x0ij(jj,kk)**lij*(a1s+Bij)             
                 da1ij      =   x0ij(jj,kk)**lij*(da1s+dBij) 
                 da1ijrhos  =   x0ij(jj,kk)**lij*(da1srhos+dBijrhos)   
                 dda1ijrhos =   x0ij(jj,kk)**lij*(dda1srhos+ddBijrhos)
                 g1ijMie    =   -c_ij(jj,kk)*lij*a1ij
                 dg1ijMie   =   -c_ij(jj,kk)*lij*(da1ij/               &
                                rhosegtotal(j)-m(jjj)*a1ij/             &
                                rhosegtotal(j)**2.d00)                 


                 lij = lambdaij(2,jj,kk)
                 call CLIJ(lij,C)                 
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 ddZetaeffrhos  = (C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*ddZetarhos+dZeta*     &
                                (2.d00*C(2)+6.d00*C(3)*                &
                                Zeta+12.d00*C(4)*                      & 
                                Zeta2)*dZetarhos
                 dZetaeffrhos   = (C(1)+2.d00*C(2)                     &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*dZetarhos   
                 call a1scalc_chain(rhosegtotal(j),Zetaeff,dZetaeff,   &
                      dZetaeffrhos,ddZetaeffrhos,m(jjj),ep,dia3,        &
                      lambdaij(2,jj,kk),a1s,da1s,da1srhos,dda1srhos)                              
                 call Bijcalc_chain(x0ij(jj,kk),lambdaij(2,jj,kk),     &
                      rhosegtotal(j),Zeta,dZeta,dZetarhos,ddZetarhos,  &
                      m(jjj),ep,dia3,Bij,dBij,dBijrhos,ddBijrhos)
                 a1ij       =   c_ij(jj,kk)*(a1ij-x0ij(jj,kk)**lij*    &
                                (a1s+Bij))             
                 da1ij      =   c_ij(jj,kk)*(da1ij-x0ij(jj,kk)**lij*   &
                                (da1s+dBij))
                 da1ijrhos  =   c_ij(jj,kk)*(da1ijrhos-x0ij(jj,kk)**   &
                                lij*(da1srhos+dBijrhos))
                 dda1ijrhos =   c_ij(jj,kk)*(dda1ijrhos-x0ij(jj,kk)**  &
                                lij*(dda1srhos+ddBijrhos))
                 g1ijMie    =   (g1ijMie+c_ij(jj,kk)*lij*x0ij(jj,kk)** &
                                lij*(a1s+Bij))*1.d00/rhosegtotal(j)+   &
                                3.d00*da1ijrhos
                 g1ijMie    =   0.5d00*g1ijMie/(pi*ep*dia3) 
                 dg1ijMie   =   dg1ijMie+c_ij(jj,kk)*lij*((x0ij(jj,kk) &
                                **lij*(da1s+dBij))/rhosegtotal(j)-     &
                                m(jjj)*x0ij(jj,kk)**lij*(a1s+Bij)/     &
                                rhosegtotal(j)**2.d00)+3.d00*dda1ijrhos
                 dg1ijMie   =   0.5d00*dg1ijMie/(pi*ep*dia3)
                 
        
                 lij = 2.d00*lambdaij(1,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 ddZetaeffrhos  = (C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*ddZetarhos+dZeta*     &
                                (2.d00*C(2)+6.d00*C(3)*                &
                                Zeta+12.d00*C(4)*                      & 
                                Zeta2)*dZetarhos
                 dZetaeffrhos   = (C(1)+2.d00*C(2)                     &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*dZetarhos 
                 call a1scalc_chain(rhosegtotal(j),Zetaeff,dZetaeff,   &
                      dZetaeffrhos,ddZetaeffrhos,m(jjj),ep,dia3,lij,    &
                      a1s,da1s,da1srhos,dda1srhos)                                
                 call Bijcalc_chain(x0ij(jj,kk),lij,rhosegtotal(j),    &
                      Zeta,dZeta,dZetarhos,ddZetarhos,m(jjj),ep,dia3,   &
                      Bij,dBij,dBijrhos,ddBijrhos) 
                 a2ij       =   x0ij(jj,kk)**lij*(a1s+Bij) 
                 a2ijchain  =   -lij*0.5d00*x0ij(jj,kk)**lij*(a1s+Bij)
                 da2ijchain =   -lij*0.5d00*x0ij(jj,kk)**lij*(da1s+dBij)                                                  
                 da2ij      =   x0ij(jj,kk)**lij*(da1s+dBij)
                 dTa2ijrhos  =   x0ij(jj,kk)**lij*(da1srhos+dBijrhos)   
                 ddTa2ijrhos =   x0ij(jj,kk)**lij*(dda1srhos+ddBijrhos)
                 
                 lij = 2.d00*lambdaij(2,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 ddZetaeffrhos  = (C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*ddZetarhos+dZeta*     &
                                (2.d00*C(2)+6.d00*C(3)*                &
                                Zeta+12.d00*C(4)*                      & 
                                Zeta2)*dZetarhos
                 dZetaeffrhos   = (C(1)+2.d00*C(2)                     &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*dZetarhos 
                 call a1scalc_chain(rhosegtotal(j),Zetaeff,dZetaeff,   &
                      dZetaeffrhos,ddZetaeffrhos,m(jjj),ep,dia3,lij,    &
                      a1s,da1s,da1srhos,dda1srhos)                                
                 call Bijcalc_chain(x0ij(jj,kk),lij,rhosegtotal(j),    &
                      Zeta,dZeta,dZetarhos,ddZetarhos,m(jjj),ep,dia3,   &
                      Bij,dBij,dBijrhos,ddBijrhos) 
                 a2ij       =   a2ij+x0ij(jj,kk)**lij*(a1s+Bij)  
                 a2ijchain  =   a2ijchain-0.5d00*lij*x0ij(jj,kk)**lij*(a1s+Bij)  
                 da2ijchain  =  da2ijchain-0.5d00*lij*x0ij(jj,kk)**lij*(da1s+dBij)                                              
                 da2ij      =   da2ij+x0ij(jj,kk)**lij*(da1s+dBij)                  
                 dTa2ijrhos  =   dTa2ijrhos+x0ij(jj,kk)**lij*(da1srhos+dBijrhos)   
                 ddTa2ijrhos =   ddTa2ijrhos+x0ij(jj,kk)**lij*(dda1srhos+ddBijrhos)
                 lij = lambdaij(1,jj,kk)+lambdaij(2,jj,kk)
                 call CLIJ(lij,C)  
                 Zetaeff    =   C(1)*Zeta+C(2)*Zeta2+                  &
                                C(3)*Zeta3+C(4)*Zeta4
                 dZetaeff   =   dZeta*(C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3) 
                 ddZetaeffrhos  = (C(1)+2.d00*C(2)                 &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*ddZetarhos+dZeta*     &
                                (2.d00*C(2)+6.d00*C(3)*                &
                                Zeta+12.d00*C(4)*                      & 
                                Zeta2)*dZetarhos
                 dZetaeffrhos   = (C(1)+2.d00*C(2)                     &
                                *Zeta+3.d00*C(3)*Zeta2+4.d00*          &
                                C(4)*Zeta3)*dZetarhos 
                 call a1scalc_chain(rhosegtotal(j),Zetaeff,dZetaeff,   &
                      dZetaeffrhos,ddZetaeffrhos,m(jjj),ep,dia3,lij,    &
                      a1s,da1s,da1srhos,dda1srhos)                                
                 call Bijcalc_chain(x0ij(jj,kk),lij,rhosegtotal(j),    &
                      Zeta,dZeta,dZetarhos,ddZetarhos,m(jjj),ep,dia3,   &
                      Bij,dBij,dBijrhos,ddBijrhos) 
                 a2ijtemp   =   a2ij-2.d00*x0ij(jj,kk)**lij*(a1s+Bij)             
                 a2ij       =   0.5d00*K_hs*(1.d00+Xiij)*ep*c_ij(jj,kk)**2.d00*a2ijtemp
                 a2ijchain  =   a2ijchain+lij*x0ij(jj,kk)**lij*(a1s+Bij)  
                 da2ijchain =   da2ijchain+lij*x0ij(jj,kk)**lij*(da1s+dBij)                                    
                 da2ij      =   da2ij-2.d00*x0ij(jj,kk)**lij*(da1s+dBij)
                 da2ijtemp  =   da2ij 
                 da2ij      =   da2ij*K_hs*(1.d00+Xiij)+dK_hs*(1.d00+Xiij)*a2ijtemp+K_hs*dXiij*a2ijtemp
                 da2ij      =   (0.5d00*ep*c_ij(jj,kk)**2.d00)*da2ij

                 dTa2ijrhos  =  dTa2ijrhos-2.d00*x0ij(jj,kk)**lij*(da1srhos+dBijrhos)   
                 ddTa2ijrhos =  ddTa2ijrhos-2.d00*x0ij(jj,kk)**lij*(dda1srhos+ddBijrhos)
                 
                 da2ijrhos   =  0.5d00*ep*c_ij(jj,kk)**2.d00*         &
                                (dK_hsrhos*(1.d00+Xiij)*a2ijtemp+     &
                                K_hs*a2ijtemp*dXiijrhos+K_hs*         &
                                (1.d00+Xiij)*dTa2ijrhos)
                 dda2ijrhos  =  0.5d00*ep*c_ij(jj,kk)**2.d00*         &
                                ((1.d00+Xiij)*dK_hs*dTa2ijrhos+       &
                                a2ijtemp*(1.d00+Xiij)*ddK_hsrhos+     &
                                dK_hs*a2ijtemp*dXiijrhos)
                 dda2ijrhos  =  dda2ijrhos+0.5d00*ep*c_ij(jj,kk)**    &
                                2.d00*(dK_hsrhos*dXiij*a2ijtemp+      &
                                ddXiijrhos*K_hs*a2ijtemp+K_hs*        &
                                dXiij*dTa2ijrhos)
                 dda2ijrhos =   dda2ijrhos+0.5d00*ep*c_ij(jj,kk)**    &
                                2.d00*(dK_hsrhos*(1.d00+Xiij)*        &
                                da2ijtemp+K_hs*dXiijrhos*da2ijtemp+   &
                                K_hs*(1.d00+Xiij)*ddTa2ijrhos)

                 invXi      =   1.d00/(1.d00+Xiij) 
                 da2byXi    =   -a2ij*(invXi**2.d00)*dXiijrhos+invXi*da2ijrhos                    
                 dda2byXi   =   dXiij*(-invXi**2.d00*da2ijrhos+       &
                                2.d00*a2ij*invXi**3.d00*dXiijrhos)
                 dda2byXi   =   dda2byXi-a2ij*invXi**2.d00*ddXiijrhos &
                                -da2ij*invXi**2.d00*dXiijrhos
                 dda2byXi   =   dda2byXi+invXi*dda2ijrhos  

                 g2ijMCA    =   3.d00*da2byXi+ep/rhosegtotal(j)*K_hs*  &
                                c_ij(jj,kk)**2.d00*a2ijchain 
                 dg2ijMCA   =   3.d00*dda2byXi+ep*K_hs*c_ij(jj,        &
                                kk)**2.d00*(1.d00/rhosegtotal(j)*      &
                                da2ijchain-1.d00/rhosegtotal(j)**      &
                                2.d00*m(jjj)*a2ijchain)+ep*dK_hs*       &
                                c_ij(jj,kk)**2.d00*1.d00/              &
                                rhosegtotal(j)*a2ijchain
                    
                 g2ijMCA    =   g2ijMCA/(2.d00*pi*ep*ep*dia3)   
                 dg2ijMCA   =   dg2ijMCA/(2.d00*pi*ep*ep*dia3)  

                 thetaij    =   dexp(ep)-1.d00
                 alpha      =   alphaij(c_ij(jj,kk),lambdaij(2,jj,kk),lambdaij(1,jj,kk))                 
                 gamacij    =   10.d00*(-dtanh(10.d00*                  &
                                (0.57d00-alpha))+1.d00)*thetaij*         &
                                Zetaavg*dexp(Zetaavg*(-6.7d00)+Zetaavg*    &
                                Zetaavg*(-8d00))
                 dgamacij   =   10.d00*(-dtanh(10.d00*                  &
                                (0.57d00-alpha))+1.d00)*thetaij*         &
                                dZetaavg*dexp(Zetaavg*(-6.7d00)+Zetaavg*   &
                                Zetaavg*(-8d00))+gamacij*dZetaavg*       &
                                ((-6.7d00)+2.d00*Zetaavg*(-8d00)) 
                 g2ijMie    =   (1.d00+gamacij)*g2ijMCA
                 dg2ijMie   =   g2ijMCA*dgamacij+(1.d00+gamacij)*      &
                                dg2ijMCA
                 gijHS      =   dexp(k0+k1*x0ij(jj,kk)+k2*x0ij(jj,kk)  &
                                **2.d00+k3*x0ij(jj,kk)**3.d00)       
                 dgijHS     =   gijHS*(dk0+x0ij(jj,kk)*dk1+x0ij(jj,kk) &
                                **2.d00*dk2+x0ij(jj,kk)**3.d00*dk3)  
                 gijMie(jj,kk)  = gijHS*dexp(ep*g1ijMie/gijHS+ep*ep*   &
                                g2ijMie/gijHS) 
                 dg1            = ep/(gijHS**2.d00)*(dg1ijMie*gijHS-     &
                                g1ijMie*dgijHS) 
                 dg2            = ep*ep/(gijHS**2.d00)*(dg2ijMie*gijHS-  &
                                g2ijMie*dgijHS)                  
                 dgijMie(jj,kk) = dgijHS*dexp(ep/gijHS*(g1ijMie+ep*  &
                                g2ijMie))+gijMie(jj,kk)*(dg1+dg2)

              enddo
           enddo
           

                 achain1(jjj,j) = -(m(jjj)-1.d00)*dlog(gijMie(jjj,jjj))

           do ii=1,nc                 
                 achain4(jjj,j,ii) = -(m(ii)-1.d00)*(dgijMie(ii,ii))*rho_orig(ii,j)/gijMie(ii,ii)
           enddo
         enddo
        enddo

        dF_chain(:,:)  = 0.d00
        dF_chain1(:,:) = 0.d00
        dF_chain2(:,:) = 0.d00
        dF_chain3(:,:) = 0.d00
        dF_chain4(:,:) = 0.d00
        do i=1,nc
           call fft_weightfun (psi(i),dij(i,i),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain) 
           dF_chain4(:,:)=0.d00
           chain4(:)=0.d00 
           do jj=1,nc
             call fft(achain4(i,1:ext_bins,jj),Fchain,chain4(1:ext_bins),ext_bins,cplxn)
             dF_chain4(i,1:ext_bins)=dF_chain4(i,1:ext_bins)+chain4(1:ext_bins)
           enddo

           k=1          
           do ii=dif_bins+1,n_bins+dif_bins
              df_chain(i,k) = achain1(i,ii)+dF_chain4(i,ii)            
              k=k+1   
           enddo
        enddo
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine association (ct,dF_assoc)
	use globalvar	
	implicit none
        integer*8, parameter  :: cplxp = 1
        integer*8, parameter  :: cplxn = 0 
        integer   ::   ct 
        integer*8 ::   i,j,k,ii,jj,kk,jjj
        integer*8 ::   z_index,dif_bins   
        real*8    ::   rho_ft(nc,ext_bins),temp(nc,ext_bins),rho_orig(nc,ext_bins)    
        real*8, dimension (ext_bins) :: Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec
        real*8, dimension (ext_bins) :: c0ft,c2ft,c2_vecft,c3ft       
        real*8, dimension (ext_bins) :: Fdisp,Fchain
        real*8, dimension (ext_bins) :: dphidn0,dphidn2,dphidn2_vec,dphidn3
        real*8, dimension (nc,n_bins)     :: dF_assoc 
        real*8                            :: xaitotal(ext_bins,2),xa(2)   
        real*8                            :: ddn0,ddn2,ddn2vec,ddn3   

        rho_ft(:,:) = 0.d00
        dif_bins  = dint((ext_bins-n_bins)/2.d00)
        if(ct.eq.1)then
          do i=1,nc
	     z_index   = 1           
             do j=dif_bins+1,n_bins+dif_bins
                rho_ft(i,j) = density_calc(i,z_index)
                z_index   = z_index+1
             enddo 
          enddo
          z_index   = 1
        else
          do i=1,nc
             rho_ft(i,1:ext_bins)=bulk_density_i(i)
          enddo
        endif

           call fft_weightfun (psi(1),dij(1,1),Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain) 
           call fft(rho_ft(1,:),Fw0,n0global,ext_bins,cplxn)
           call fft(rho_ft(1,:),Fw2,n2global,ext_bins,cplxn)
           call fft(rho_ft(1,:),Fw3,n3global,ext_bins,cplxn)
           call fft(rho_ft(1,:),Fw2_vec,n2vecglobal,ext_bins,cplxp) 
 
        do j = 1,ext_bins
           if(n0global(j).lt.1.d-100.or.n2global(j).lt.1.d-100)then
              xa(:)=0.9999
              n0global(j) = 1.d-100
              n2global(j) = 1.d-100 
              n3global(j) = 1.d-100             
              dphidn0(j) = 1.d-100
              dphidn2(j) = 1.d-100    
              dphidn2_vec(j) = 1.d-100   
              dphidn3(j) = 1.d-100                                     
           else
           call assoc(j,n0global(j),n2global(j),n2vecglobal(j),n3global(j),xa,ddn0,ddn2,ddn2vec,ddn3)
           dphidn0(j) = ddn0                 
           dphidn2(j) = ddn2          
           dphidn2_vec(j) = ddn2vec
           dphidn3(j) = ddn3
           xaitotal(j,1) = xa(1)
           endif

        enddo 


        call fft(dphidn0,Fw0,c0ft,ext_bins,cplxn)
        call fft(dphidn2,Fw2,c2ft,ext_bins,cplxn)
        call fft(dphidn3,Fw3,c3ft,ext_bins,cplxn)        
        call fft(dphidn2_vec,-Fw2_vec,c2_vecft,ext_bins,cplxp)
        k=1
        do ii = dif_bins+1,n_bins+dif_bins         
           dF_assoc(1,k) = (c0ft(ii)+c2ft(ii)+c2_vecft(ii)+c3ft(ii))
           k=k+1
        enddo

end subroutine      

subroutine fft_weightfun (psi_i,dii,Fw0,Fw1,Fw2,Fw3,Fw1_vec,Fw2_vec,Fdisp,Fchain)
   use globalvar
   implicit none
   integer*8           :: i,k,half_ext_bins
   real*8              :: csi,csi0,dii,psi_i
   real*8              :: surf,vol,inv_R
   real*8              :: j0,j2,j0disp,j2disp,R,j0_ch,j2_ch 
   real*8              :: Fw0(ext_bins),Fw1(ext_bins),Fw2(ext_bins),Fw3(ext_bins)  
   real*8              :: Fw1_vec(ext_bins),Fw2_vec(ext_bins),Fdisp(ext_bins),Fchain(ext_bins)


   R=dii*0.5d00
   length=(ext_bins-1d00)*delta_z
   inv_length=1.d00/length
   csi0=2.d00*pi*R*inv_length
   half_ext_bins=int(ext_bins*0.5d00)
   surf=4.d00*pi*R**2.d00
   vol=4.d00*pi*R*R*R/3.d00
   inv_R=1.d00/R

   Fw0(1)     = 1d00 
   Fw1(1)     = R
   Fw2(1)     = surf      
   Fw3(1)     = vol
   Fw1_vec(1) = 0d00 
   Fw2_vec(1) = 0d00
   Fdisp(1)   = 1d00
   Fchain(1)  = 1d00
   do k=1,half_ext_bins-1
      csi=dabs(csi0*k)
      j0           = dble(sin(csi)/csi)
      j2           = j0*(3.d00/csi**2.d00-1.d00)-3.d00*                &
                     dble(dcos(csi)/csi**2.d00)
      j0_ch        = dble(sin(2.d0*csi)/(2.d0*csi))
      j2_ch        = j0_ch*(3.d00/(2.d00*csi)**2.d00-1.d00)-3.d00*        &
                     dble(dcos(2.d00*csi)/(2.d00*csi)**2.d00)                        
      j0disp       = dble(sin(2.d00*csi*psi_i)/(2.d00*csi*psi_i))
      j2disp       = j0disp*(3.d00/(2.d00*csi*psi_i)**2.d00-1.d00)-              &
                     3.d00*dble(dcos((2.d00*csi*psi_i))/(2.d00*csi*psi_i)**2.d00)
      Fw0(k+1)     = j0
      Fw1(k+1)     = R*j0
      Fw2(k+1)     = surf*j0  
      Fw3(k+1)     = vol*(j0+j2) 
      Fw1_vec(k+1) = -Fw3(k+1)*k*0.5d00*inv_R*inv_length 
      Fw2_vec(k+1) = -2.d00*pi*Fw3(k+1)*k*inv_length
      Fdisp(k+1)   = j0disp+j2disp
      Fchain(k+1)  = j0_ch+j2_ch   
   enddo
   i=half_ext_bins+1
   do k=-half_ext_bins,-1,1
      csi=dabs(csi0*k)
      j0         = dble(sin(csi)/csi)
      j2         = j0*(3.d00/csi**2.d00-1.d00)-3.d00*                  &
                   dble(dcos(csi)/csi**2.d00)
      j0_ch        = dble(sin(2.d0*csi)/(2.d0*csi))
      j2_ch        = j0_ch*(3.d00/(2.d00*csi)**2.d00-1.d00)-3.d00*        &
                     dble(dcos(2.d00*csi)/(2.d00*csi)**2.d00)                     
      j0disp       = dble(sin(2.d00*csi*psi_i)/(2.d00*csi*psi_i))
      j2disp       = j0disp*(3.d00/(2.d00*csi*psi_i)**2.d00-1.d00)-              &
                     3.d00*dble(dcos((2.d00*csi*psi_i))/(2.d00*csi*psi_i)**2.d00)
      Fw0(i)     = j0
      Fw1(i)     = R*j0
      Fw2(i)     = surf*j0      
      Fw3(i)     = vol*(j0+j2) 
      Fw1_vec(i) = -Fw3(i)*k*0.5d00*inv_R*inv_length  
      Fw2_vec(i) = -2.d00*pi*Fw3(i)*k*inv_length
      Fdisp(i)   = j0disp+j2disp
      Fchain(i)  = j0_ch+j2_ch  
      i=i+1
   enddo
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine fft (data_vec1,data_vec2,conv,nn,cplx)
use globalvar
implicit none
integer*8           :: nn,i,j,ii,cplx,isig
real*8              :: nn_inv
real*8              :: data_vec1 (ext_bins),data_vec2 (ext_bins)
real*8              :: conv (ext_bins)
real*8, allocatable :: data_vec1ft (:),data_vec2ft (:)
real*8, allocatable :: conv_ft (:)

allocate(data_vec1ft(2*ext_bins))
allocate(data_vec2ft(2*ext_bins))
allocate(conv_ft(2*ext_bins))
data_vec1ft(:)=0.d00
data_vec2ft(:)=0.d00
conv_ft(:)=0.d00

if(cplx.eq.1)then
   do i=1,nn
      data_vec1ft(2*i-1)=data_vec1(i)
      data_vec2ft(2*i)=data_vec2(i)
   enddo
else
   do i=1,nn
      data_vec1ft(2*i-1)=data_vec1(i)
      data_vec2ft(2*i-1)=data_vec2(i)  
   enddo
endif

nn_inv=dble(1.d00/nn)



isig=1
call four1 (data_vec1ft,nn,isig)


do ii=1,nn
    conv_ft(2*ii-1)=data_vec1ft(2*ii-1)*data_vec2ft(2*ii-1)- &
                  data_vec2ft(2*ii)*data_vec1ft(2*ii)
    conv_ft(2*ii)=data_vec1ft(2*ii-1)*data_vec2ft(2*ii)+     &
                  data_vec1ft(2*ii)*data_vec2ft(2*ii-1)
enddo
isig=-1 
call four1 (conv_ft,nn,isig)

do j=1,nn
  conv(j)=conv_ft(2*j-1)*nn_inv 
enddo

deallocate(data_vec1ft)
deallocate(data_vec2ft)
deallocate(conv_ft)

end subroutine

SUBROUTINE four1(dat,nn,isig)
INTEGER*8 isig,nn
REAL*8 dat(2*nn)
!Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
!data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as −1.
!data is a complex array of length nn or, equivalently, a real array of length 2*nn . nn
!MUST be an integer power of 2 (this is not checked for!).
INTEGER*8 i,istep,j,m,mmax,n
REAL*8 tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
!Double precision for the trigonometric recurrences. 
n=2*nn
j=1
!This is the bit-reversal section of the routine.
do i=1,n,2
if(j.gt.i)then
!Exchange the two complex numbers.
tempr=dat(j)
tempi=dat(j+1)
dat(j)=dat(i)
dat(j+1)=dat(i+1)
dat(i)=tempr
dat(i+1)=tempi
endif
m=n/2
1 if ((m.ge.2).and.(j.gt.m)) then
j=j-m
m=m/2
goto 1
endif
j=j+m
enddo 
mmax=2
!Here begins the Danielson-Lanczos section of the routine.
!Outer loop executed log 2 nn times.
2 if (n.gt.mmax) then
istep=2*mmax
!Initialize for the trigonometric recurrence.
theta=6.28318530717959d00/(isig*mmax)
wpr=-2.d00*dsin(0.5d00*theta)**2.d00
wpi=dsin(theta)
wr=1.d00
wi=0.d00
!Here are the two nested inner loops.
do  m=1,mmax,2
do  i=m,n,istep
j=i+mmax
!This is the Danielson-Lanczos formula:
tempr=(wr)*dat(j)-(wi)*dat(j+1)
tempi=(wr)*dat(j+1)+(wi)*dat(j)
dat(j)=dat(i)-tempr
dat(j+1)=dat(i+1)-tempi
dat(i)=dat(i)+tempr
dat(i+1)=dat(i+1)+tempi
enddo 
!Trigonometric recurrence.
wtemp=wr
wr=wr*wpr-wi*wpi+wr
wi=wi*wpr+wtemp*wpi+wi
enddo 
mmax=istep
goto 2
endif
return
END
