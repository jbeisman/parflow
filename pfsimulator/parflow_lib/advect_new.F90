!BHEADER**********************************************************************
!
!  Copyright (c) 1995-2009, Lawrence Livermore National Security,
!  LLC. Produced at the Lawrence Livermore National Laboratory. Written
!  by the Parflow Team (see the CONTRIBUTORS file)
!  <parflow@lists.llnl.gov> CODE-OCE!08-103. All rights reserved.
!
!  This file is part of Parflow. For details, see
!  http://www.llnl.gov/casc/parflow
!
!  Please read the COPYRIGHT file or Our Notice and the LICENSE file
!  for the GNU Lesser General Public License.
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Publi!License (as published
!  by the Free Software Foundation) version 2.1 dated February 1999.
!
!  This program is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
!  and conditions of the GNU General Publi!License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
!  USA
!**********************************************************************EHEADER
! advect: this file contains 5 primary subroutines to simulate advective 
! transport, advect_upwind, advect_highorder, advect_transverse, advect_limit,
! and advect_computeconcen
!
!
!-------------------------------------------------------------------------------
!    advect_upwind: computes first-order upwind fluxes, updates concentrations
!    can handle transient saturations
!    
!    the new low-order concentrations should always be positive and monotone 
!    in time, although small errors (~ 0.3%) may arise in the case of transient 
!    saturations 
!-------------------------------------------------------------------------------
      subroutine advect_upwind(s,sn,uedge,vedge,wedge,phi, &
          dlo,dhi,hx,dt,old_sat,sat, &
          iteration,num_iterations,fx,fy,fz,smin,smax)

      implicit none
      integer,  parameter   :: dp = selected_real_kind(15)
      integer,  intent(in)  :: dlo(3), dhi(3)
      real(dp), intent(in)  :: hx(3), dt
      integer,  intent(in)  :: iteration,num_iterations
      real(dp), intent(in)  :: s(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(dp), intent(out) :: sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(dp), intent(in)  :: uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(dp), intent(in)  :: vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(dp), intent(in)  :: wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(in)  :: phi(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2)
      real(dp), intent(out) :: fx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(out) :: fy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(out) :: fz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(in)  :: old_sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(dp), intent(in)  :: sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1)
      real(dp), intent(out) :: smin(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out) :: smax(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2)
 
      integer i,j,k
      integer is,ie,js,je,ks,ke
      integer ii,jj,kk
      real(dp) dx,dy,dz,dx_inv,dy_inv,dz_inv
      real(dp) sat_diff,iter,num_iter,num_iter_inv

      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)
      dx_inv = 1.0_dp/dx
      dy_inv = 1.0_dp/dy
      dz_inv = 1.0_dp/dz

      iter=DBLE(iteration)
      num_iter=DBLE(num_iterations)
      num_iter_inv = 1.0_dp/num_iter

      do k=ks,ke+1
        do j=js,je+1
          do i=is,ie+1
     
          if (uedge(i,j,k) .ge. 0.0)then
            ii =  i-1
          else
            ii =  i
          endif
    
          if (vedge(i,j,k) .ge. 0.0)then
            jj =  j-1
          else
            jj =  j
          endif

          if (wedge(i,j,k) .ge. 0.0)then
            kk =  k-1
          else
            kk =  k
          endif

          !these values give a first order scheme
          fx(i,j,k)=uedge(i,j,k)*s(ii,j,k) 
          fy(i,j,k)=vedge(i,j,k)*s(i,jj,k) 
          fz(i,j,k)=wedge(i,j,k)*s(i,j,kk)

          smin(i,j,k) = min(s(i,j,k),s(ii,j,k),s(i,jj,k),s(i,j,kk))
          smax(i,j,k) = max(s(i,j,k),s(ii,j,k),s(i,jj,k),s(i,j,kk))

          enddo
        enddo
      enddo
    
      !! first-order transport 
      do k=ks,ke
        do j=js,je
          do i=is,ie
            sat_diff=(sat(i,j,k) - old_sat(i,j,k)) * num_iter_inv

            sn(i,j,k) = ((iter*sat_diff + old_sat(i,j,k))*phi(i,j,k)*s(i,j,k) + &
            ((dt*dx_inv)*(fx(i,j,k) - fx(i+1,j,k)) + (dt*dy_inv)*(fy(i,j,k)-fy(i,j+1,k)) + &
            (dt*dz_inv)*(fz(i,j,k)-fz(i,j,k+1)))) / (((iter+1.0_dp)*sat_diff + old_sat(i,j,k))*phi(i,j,k))

          enddo
        enddo
      enddo

      ! clear memory for resuse in other advect subroutines
      fx = 0.0_dp
      fy = 0.0_dp
      fz = 0.0_dp

      return
      end subroutine advect_upwind


!-------------------------------------------------------------------------------
!    advect_highorder: computes high-order anti-diffusive flux corrections to 
!    be added to the low-order concentrations computed with advect_upwind
!    
!    employs the monotone-centered limiter to ensure montonicity in each 
!    primary dir
!
!    Lax-Wendrof style centered approximation 
!-------------------------------------------------------------------------------
    subroutine advect_highorder(s,uedge,vedge,wedge, &
          dlo,dhi,hx,dt,sx,sy,sz)

      implicit none
      integer,  parameter    :: dp = selected_real_kind(15)
      integer,  intent (in)  :: dlo(3), dhi(3)
      real(dp), intent (in)  :: hx(3), dt
      real(dp), intent (in)  :: s(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(dp), intent (in)  :: uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(dp), intent (in)  :: vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(dp), intent (in)  :: wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(dp), intent (inout) :: sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)  
      real(dp), intent (inout) :: sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent (inout) :: sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      integer i,j,k
      integer is,ie,js,je,ks,ke
      integer ks1,ke1
      integer js1,je1
      integer is1,ie1
      integer ii2,ii3,jj2,jj3,kk2,kk3
      real(dp) dx,dy,dz,dx_inv,dy_inv,dz_inv
      real(dp) half
      real(dp) rx,ry,rz,mclimit,limx,thetax,thetay,thetaz,limy,limz

      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)
      dx_inv = 1.0_dp/dx
      dy_inv = 1.0_dp/dy
      dz_inv = 1.0_dp/dz
      half = 0.5_dp

      !! make smaller and slightly faster for lower D problems
      if((ie-is) == 0) then
        is1 = is+1
        ie1 = ie-2
      else
        is1 = is
        ie1 = ie
      endif

      if((je-js) == 0) then
        js1 = js+1
        je1 = je-2
      else
        js1 = js
        je1 = je
      endif

      if((ke-ks) == 0) then
        ks1 = ks+1
        ke1 = ke-2
      else
        ks1 = ks
        ke1 = ke
      endif

      !! compute second order corrections
      do k=ks1-1,ke1+2
        do j=js1-1,je1+2
          do i=is1-1,ie1+2
     
          if (uedge(i,j,k) .ge. 0.0)then
            ii2 = i
            ii3 = i-1
          else
            ii2 = i-1
            ii3 = i+1
          endif
    
          if (vedge(i,j,k) .ge. 0.0)then
            jj2 = j
            jj3 = j-1
          else
            jj2 = j-1
            jj3 = j+1
          endif

          if (wedge(i,j,k) .ge. 0.0)then
            kk2 = k
            kk3 = k-1
          else
            kk2 = k-1
            kk3 = k+1
          endif
    
          !gradient
          rx=s(i,j,k)-s(i-1,j,k)
          ry=s(i,j,k)-s(i,j-1,k)
          rz=s(i,j,k)-s(i,j,k-1)
    
          !second order - LW flux - monotone centered limiter
          thetax=(s(ii3,j,k)-s(ii3-1,j,k))/(s(i,j,k)-s(i-1,j,k))
          limx=mclimit(thetax)
          sx(i,j,k)=half*abs(uedge(i,j,k))*(1.0_dp-(dt*dx_inv)*abs(uedge(i,j,k)))*rx*limx
    
          thetay=(s(i,jj3,k)-s(i,jj3-1,k))/(s(i,j,k)-s(i,j-1,k))
          limy=mclimit(thetay)
          sy(i,j,k) = half*abs(vedge(i,j,k))*(1.0_dp-(dt*dy_inv)*abs(vedge(i,j,k)))*ry*limy
    
          thetaz=(s(i,j,kk3)-s(i,j,kk3-1))/(s(i,j,k)-s(i,j,k-1))
          limz=mclimit(thetaz)
          sz(i,j,k)= half*abs(wedge(i,j,k))*(1.0_dp-(dt*dz_inv)*abs(wedge(i,j,k)))*rz*limz
    
          enddo
        enddo
      enddo

      return
      end subroutine advect_highorder


!-------------------------------------------------------------------------------
!    advect_transverse: computes approximate solution to Riemann problems 
!    emanating from neighboring interfaces 
!    
!    computes transverse corrections to both the first-order increment waves
!    and second-order correction waves
!-------------------------------------------------------------------------------
      subroutine advect_transverse(s,uedge,vedge,wedge, &
          dlo,dhi,hx,dt,vx,wx,uy,wy,uz,vz, &
          sx,sy,sz)

      implicit none
      integer,  parameter     :: dp = selected_real_kind(15)
      integer,  intent(in)    :: dlo(3), dhi(3)
      real(dp), intent(in)    :: hx(3), dt
      real(dp), intent(in)    :: s(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(dp), intent(in)    :: uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(dp), intent(in)    :: vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(dp), intent(in)    :: wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(out)   :: vx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: wx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp), intent(out)   :: uy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: wy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: uz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: vz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)  
      real(dp), intent(inout) :: sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)  
      real(dp), intent(inout) :: sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(inout) :: sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
 
      integer i,j,k
      integer is,ie,js,je,ks,ke
      integer ks1,ke1
      integer js1,je1
      integer is1,ie1
      integer ii2,ii3,jj2,jj3,kk2,kk3
      integer iiz,iiz2,jjx,jjx2,kky,kky2,iiy,jjz,kkx
      real(dp) dx,dy,dz,dx_inv,dy_inv,dz_inv
      real(dp) half,third
      real(dp) rx,ry,rz
      real(dp) transvel
      !real(dp) minmod4,minmod2,median

      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)
      dx_inv = 1.0_dp/dx
      dy_inv = 1.0_dp/dy
      dz_inv = 1.0_dp/dz
      half = 0.5_dp
      third = 1.0_dp/3.0_dp


      !! make smaller and faster for lower D problems
      if((ie-is) == 0) then
        is1 = is+1
        ie1 = ie-2
      else
        is1 = is
        ie1 = ie
      endif

      if((je-js) == 0) then
        js1 = js+1
        je1 = je-2
      else
        js1 = js
        je1 = je
      endif

      if((ke-ks) == 0) then
        ks1 = ks+1
        ke1 = ke-2
      else
        ks1 = ks
        ke1 = ke
      endif

      !! main loop
      !! compute increment fluxes and second order corrections
      !! compute transverse riemann problem, move fluxes appropriately
      do k=ks1-1,ke1+2
        do j=js1-1,je1+2
          do i=is1-1,ie1+2
     
          if (uedge(i,j,k) .ge. 0.0)then
            ii2 = i
            ii3 = i-1
          else
            ii2 = i-1
            ii3 = i+1
          endif
    
          if (vedge(i,j,k) .ge. 0.0)then
            jj2 = j
            jj3 = j-1
          else
            jj2 = j-1
            jj3 = j+1
          endif

          if (wedge(i,j,k) .ge. 0.0)then
            kk2 = k
            kk3 = k-1
          else
            kk2 = k-1
            kk3 = k+1
          endif

    
          !! transverse velocities - 0 if both don't have same sign
          vx(i,j,k) = (dt*dx_inv)*transvel(vedge(ii2,j,k),vedge(ii2,j+1,k)) 
          wx(i,j,k) = (dt*dx_inv)*transvel(wedge(ii2,j,k),wedge(ii2,j,k+1)) 
          uy(i,j,k) = (dt*dy_inv)*transvel(uedge(i,jj2,k),uedge(i+1,jj2,k)) 
          wy(i,j,k) = (dt*dy_inv)*transvel(wedge(i,jj2,k),wedge(i,jj2,k+1))
          uz(i,j,k) = (dt*dz_inv)*transvel(uedge(i,j,kk2),uedge(i+1,j,kk2))
          vz(i,j,k) = (dt*dz_inv)*transvel(vedge(i,j,kk2),vedge(i,j+1,kk2))
    
          !gradient
          rx=s(i,j,k)-s(i-1,j,k)
          ry=s(i,j,k)-s(i,j-1,k)
          rz=s(i,j,k)-s(i,j,k-1)
    
          
          if (vx(i,j,k) .gt. 0.0)then
            jjx = j+1
            jjx2 = j+1
          else
            jjx = j
            jjx2 = j-1
          endif
        
          if (wx(i,j,k) .gt. 0.0)then
            kkx = k+1
          else
            kkx = k
          endif
        
          if (uy(i,j,k) .gt. 0.0)then
            iiy = i+1
          else
            iiy = i
          endif
        
          if (wy(i,j,k) .gt. 0.0)then
            kky = k+1
            kky2 = k+1
          else
            kky = k
            kky2 = k-1
          endif
        
          if (uz(i,j,k) .gt. 0.0)then
            iiz = i+1
            iiz2 = i+1
          else
            iiz = i
            iiz2 = i-1
          endif
        
          if (vz(i,j,k) .gt. 0.0)then
            jjz = j+1
          else
            jjz = j
          endif
    
          ! these sweeps add transverse propogation information to the first order upwind method 
          !x sweep
          sy(ii2,jjx,k) = sy(ii2,jjx,k)-half*rx*uedge(i,j,k)*vx(i,j,k)
          sz(ii2,j,kkx)=sz(ii2,j,kkx) - half*rx*uedge(i,j,k)*wx(i,j,k) + &
          (third)*uedge(i,j,k)*abs(vx(i,j,k))*wx(i,j,k)*rx 
          sz(ii2,jjx2,kkx)=sz(ii2,jjx2,kkx) - (third)*uedge(i,j,k) * &
          abs(vx(i,j,k))*wx(i,j,k)*rx 
    
          !y sweep
          sz(i,jj2,kky) = sz(i,jj2,kky)-half*ry*vedge(i,j,k)*wy(i,j,k)
          sx(iiy,jj2,k)=sx(iiy,jj2,k) - half*ry*vedge(i,j,k)*uy(i,j,k) + &
          (third)*vedge(i,j,k)*abs(wy(i,j,k))*uy(i,j,k)*ry 
          sx(iiy,jj2,kky2)=sx(iiy,jj2,kky2) - (third)*vedge(i,j,k) * &
          abs(wy(i,j,k))*uy(i,j,k)*ry 
    
          !z sweep
          sx(iiz,j,kk2) = sx(iiz,j,kk2)-half*rz*wedge(i,j,k)*uz(i,j,k)
          sy(i,jjz,kk2)=sy(i,jjz,kk2) - half*rz*wedge(i,j,k)*vz(i,j,k) + &
          (third)*wedge(i,j,k)*abs(uz(i,j,k))*vz(i,j,k)*rz 
          sy(iiz2,jjz,kk2)=sy(iiz2,jjz,kk2) - (third)*wedge(i,j,k) * &
          abs(uz(i,j,k))*vz(i,j,k)*rz
    
          ! Add second order transverse propogation information
          !x sweep
          sy(i,jjx,k)      = sy(i,jjx,k) + vx(i,j,k)*sx(i,j,k)
          sy(i-1,jjx,k)    = sy(i-1,jjx,k) - vx(i,j,k)*sx(i,j,k)
          sz(i,j,kkx)      = sz(i,j,kkx) + (1.0_dp - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
          sz(i,jjx2,kkx)   = sz(i,jjx2,kkx) + abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)
          sz(i-1,j,kkx)    = sz(i-1,j,kkx) - (1.0_dp - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
          sz(i-1,jjx2,kkx) = sz(i-1,jjx2,kkx) - abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)
    
          !y sweep
          sz(i,j,kky)      = sz(i,j,kky) + wy(i,j,k)*sy(i,j,k)
          sz(i,j-1,kky)    = sz(i,j-1,kky) - wy(i,j,k)*sy(i,j,k)
          sx(iiy,j,k)      = sx(iiy,j,k) + (1.0_dp - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
          sx(iiy,j,kky2)   = sx(iiy,j,kky2) + abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)
          sx(iiy,j-1,k)    = sx(iiy,j-1,k) - (1.0_dp - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
          sx(iiy,j-1,kky2) = sx(iiy,j-1,kky2) - abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)
    
          !z sweep
          sx(iiz,j,k)      = sx(iiz,j,k) + uz(i,j,k)*sz(i,j,k)
          sx(iiz,j,k-1)    = sx(iiz,j,k-1) - uz(i,j,k)*sz(i,j,k)
          sy(i,jjz,k)      = sy(i,jjz,k) + (1.0_dp - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
          sy(iiz2,jjz,k)   = sy(iiz2,jjz,k) + abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k)
          sy(i,jjz,k-1)    = sy(i,jjz,k-1) - (1.0_dp - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
          sy(iiz2,jjz,k-1) = sy(iiz2,jjz,k-1) - abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k) 
          
          enddo
        enddo
      enddo

      return
      end subroutine advect_transverse


!-------------------------------------------------------------------------------
!    advect_computeconcen: adds high-order corrections from advect_high-order
!    and transverse corrections from advect_transverse to low-order solution 
!    from advect_upwind 
!
!    Also accounts for transient saturations    
!-------------------------------------------------------------------------------
      subroutine advect_computeconcen(sn,phi, &
          dlo,dhi,hx,dt,old_sat,sat, &
          iteration,num_iterations, &
          smin,smax,sx,sy,sz)

      implicit none
      integer,  parameter     :: dp = selected_real_kind(15)
      integer,  intent(in)    :: dlo(3), dhi(3)
      real(dp), intent(in)    :: hx(3), dt
      integer,  intent(in)    :: iteration,num_iterations
      real(dp), intent(inout) :: sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(dp), intent(in)    :: phi(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(dp), intent(in)    :: smin(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(in)    :: smax(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(in)    :: sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)  
      real(dp), intent(in)    :: sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(in)    :: sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(in)    :: old_sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(dp), intent(in)    :: sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1)
 
      integer i,j,k
      integer is,ie,js,je,ks,ke
      real(dp) dx,dy,dz,dx_inv,dy_inv,dz_inv
      real(dp) sat_diff,iter,num_iter,num_iter_inv

      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)
      dx_inv = 1.0_dp/dx
      dy_inv = 1.0_dp/dy
      dz_inv = 1.0_dp/dz

      iter=DBLE(iteration)
      num_iter=DBLE(num_iterations)
      num_iter_inv = 1.0_dp/num_iter


      do k=ks,ke
        do j=js,je
          do i=is,ie

            sat_diff=(sat(i,j,k)-old_sat(i,j,k)) * num_iter_inv

            sn(i,j,k)=sn(i,j,k) + ((dt*dx_inv)*(sx(i,j,k) - sx(i+1,j,k)) &
            + (dt*dy_inv)*(sy(i,j,k)-sy(i,j+1,k)) + &
            (dt*dz_inv)*(sz(i,j,k)-sz(i,j,k+1))) / &
            (((iter+1.0_dp)*sat_diff + old_sat(i,j,k))*phi(i,j,k))

          enddo
        enddo
      enddo


      do k=ks,ke
        do j=js,je
          do i=is,ie  

           !! cutoff any values that violate min/max
           !! this typically doesn't happen if the limiter was called
           !! but transient saturations can sometimes cause small errors, about 0.3%
           !! this keeps that from happening 
           if (sn(i,j,k) .lt. smin(i,j,k)) sn(i,j,k) = smin(i,j,k)
           if (sn(i,j,k) .gt. smax(i,j,k)) sn(i,j,k) = smax(i,j,k)  

          enddo
        enddo
      enddo 

      return
      end subroutine advect_computeconcen


!-------------------------------------------------------------------------------
!    advect_limit: multi-dimensional limiter to ensure siolution is TVD and that
!    min-max is not violated
!
!    modeled after Zalesak FCT limiter    
!-------------------------------------------------------------------------------
      subroutine advect_limit(sn,sx,sy,sz,dlo,dhi,hx,dt,&
                       p_plus,p_minus,q_plus,q_minus,r_plus,r_minus)
      implicit none
      integer,  parameter     :: dp = selected_real_kind(15)
      integer,  intent(in)    :: dlo(3), dhi(3)
      real(dp), intent(in)    :: dt,hx(3)
      real(dp), intent(in)    :: sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(dp), intent(inout) :: sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp), intent(inout) :: sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)
      real(dp), intent(inout) :: sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)
      real(dp), intent(out)   :: p_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: p_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp), intent(out)   :: q_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: q_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp), intent(out)   :: r_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp), intent(out)   :: r_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)

      integer  is,ie,js,je,ks,ke,i,j,k 
      real(dp) dx_inv,dy_inv,dz_inv

      !! initialize to 0
      p_plus = 0.0_dp
      p_minus = 0.0_dp
      q_plus = 0.0_dp
      q_minus = 0.0_dp
      r_plus = 0.0_dp
      r_minus = 0.0_dp

      dx_inv = 1.0_dp/hx(1)
      dy_inv = 1.0_dp/hx(2)
      dz_inv = 1.0_dp/hx(3)

      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
       
      
      do k=ks-1,ke+1
        do j=js-1,je+1
          do i=is-1,ie+1 
            
            p_plus(i,j,k) = (dt*dx_inv)*(max(0.0_dp,sx(i,j,k)) - min(0.0_dp,sx(i+1,j,k))) + &
            (dt*dy_inv)*(max(0.0_dp,sy(i,j,k)) - min(0.0_dp,sy(i,j+1,k))) + &
            (dt*dz_inv)*(max(0.0_dp,sz(i,j,k)) - min(0.0_dp,sz(i,j,k+1)))

            p_minus(i,j,k) = (dt*dx_inv)*(max(0.0_dp,sx(i+1,j,k)) - min(0.0_dp,sx(i,j,k))) + &
            (dt*dy_inv)*(max(0.0_dp,sy(i,j+1,k)) - min(0.0_dp,sy(i,j,k))) + &
            (dt*dz_inv)*(max(0.0_dp,sz(i,j,k+1)) - min(0.0_dp,sz(i,j,k)))

            q_plus(i,j,k) =  max(sn(i-1,j,k),sn(i,j,k),sn(i+1,j,k), sn(i,j-1,k),sn(i,j+1,k), sn(i,j,k-1),sn(i,j,k+1)) - sn(i,j,k)

            q_minus(i,j,k) = sn(i,j,k) - min(sn(i-1,j,k),sn(i,j,k),sn(i+1,j,k), sn(i,j-1,k),sn(i,j+1,k), sn(i,j,k-1),sn(i,j,k+1))

            if (p_plus(i,j,k) .gt. 0.0_dp) then
              r_plus(i,j,k) = min(q_plus(i,j,k)/p_plus(i,j,k),1.0_dp)
            else
              r_plus(i,j,k) = 0.0_dp
            endif 
   
            if (p_minus(i,j,k) .gt. 0.0_dp) then
              r_minus(i,j,k) = min(q_minus(i,j,k)/p_minus(i,j,k),1.0_dp)
            else
              r_minus(i,j,k) = 0.0_dp
            endif 
          
          enddo
        enddo
      enddo


      do k=ks,ke+1
        do j=js,je+1
          do i=is,ie+1
       
            if (sx(i,j,k) .gt. 0.0_dp) then
              sx(i,j,k) = sx(i,j,k) * min(r_plus(i,j,k),r_minus(i-1,j,k))
            else
              sx(i,j,k) = sx(i,j,k) * min(r_plus(i-1,j,k),r_minus(i,j,k))
            endif
            
            if (sy(i,j,k) .gt. 0.0_dp) then
               sy(i,j,k) = sy(i,j,k) * min(r_plus(i,j,k),r_minus(i,j-1,k))
            else
               sy(i,j,k) = sy(i,j,k) * min(r_plus(i,j-1,k),r_minus(i,j,k))
            endif
      
            if (sz(i,j,k) .gt. 0.0_dp) then
               sz(i,j,k) = sz(i,j,k) * min(r_plus(i,j,k),r_minus(i,j,k-1))
            else
               sz(i,j,k) = sz(i,j,k) * min(r_plus(i,j,k-1),r_minus(i,j,k))
            endif

          enddo
        enddo
      enddo
   
      return
      end subroutine advect_limit




!      !! dispersion calculations -- likely need to be updated for boundary conditions
!      subroutine disperse(sn,uedge,vedge,wedge,dlo,dhi,hx,dt) 
!      implicit none
!
!      integer, parameter :: dp = selected_real_kind(15)            
!      integer i,j,k,ie,is,je,js,ke,ks
!      integer dlo(3), dhi(3)
!      real(dp)  hx(3), dt, al, at,dx,dy,dz
!      
!      real(dp) uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
!      real(dp) vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
!      real(dp) wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
!      real(dp) v_xx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_xy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_xz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) v_yy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_yx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_yz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) v_zz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_zx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) v_zy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
!      real(dp) abs_vx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) abs_vy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) abs_vz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) d_xx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_xy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_xz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) d_yy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_yx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_yz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) d_zz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_zx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) d_zy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      real(dp) dc_x(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) dc_y(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
!      real(dp) dc_z(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
!      
!      al = .00625
!      at = .000625
!      
!      is = dlo(1)
!      ie = dhi(1)
!      js = dlo(2)
!      je = dhi(2)
!      ks = dlo(3)
!      ke = dhi(3)
!      dx = hx(1)
!      dy = hx(2)
!      dz = hx(3)
!  
!      !velocities
!      do k=ks,ke+1
!      do j=js,je+1
!      do i=is,ie+1
!      
!      d_xx(i,j,k) = 0.0
!      d_xy(i,j,k) = 0.0
!      d_xz(i,j,k) = 0.0
!      
!      d_yy(i,j,k) = 0.0
!      d_yx(i,j,k) = 0.0
!      d_yz(i,j,k) = 0.0
!      
!      d_zz(i,j,k) = 0.0
!      d_zx(i,j,k) = 0.0
!      d_zy(i,j,k) = 0.0
!      
!     
!      
!     
!     v_xx(i,j,k)=uedge(i,j,k)
!     v_xy(i,j,k) = sum(vedge(i-1:i,j:j+1,k))/4.0
!     v_xz(i,j,k) = sum(wedge(i-1:i,j,k:k+1))/4.0
!    
!     v_yy(i,j,k)=vedge(i,j,k) 
!     v_yx(i,j,k) = sum(uedge(i:i+1,j-1:j,k))/4.0
!     v_yz(i,j,k) = sum(wedge(i,j-1:j,k:k+1))/4.0
!     
!     v_zz(i,j,k)=wedge(i,j,k)
!     v_zx(i,j,k) = sum(uedge(i:i+1,j,k-1:k))/4.0
!     v_zy(i,j,k) = sum(vedge(i,j:j+1,k-1:k))/4.0
!     
!     abs_vx(i,j,k) = (v_xx(i,j,k)**2 + v_xy(i,j,k)**2 + v_xz(i,j,k)**2)**0.5
!     abs_vy(i,j,k) = (v_yy(i,j,k)**2 + v_yx(i,j,k)**2 + v_yz(i,j,k)**2)**0.5
!     abs_vz(i,j,k) = (v_zz(i,j,k)**2 + v_zx(i,j,k)**2 + v_zy(i,j,k)**2)**0.5 
!     
!     enddo
!     enddo
!     enddo
!      
!     
!      do i=is,ie+1
!      do j=js,je
!      do k=ks,ke
!      
!      if (abs_vx(i,j,k) /= 0.0) THEN
!      
!      d_xx(i,j,k) = at*abs_vx(i,j,k) + (al-at)*v_xx(i,j,k)*v_xx(i,j,k)/abs_vx(i,j,k)
!      else 
!      d_xx(i,j,k) = 0.0 
!      endif
!      enddo
!      enddo
!      enddo
!      
!      do i=is,ie
!      do j=js,je+1
!      do k=ks,ke
!      
!      if (abs_vy(i,j,k) /= 0.0) THEN
!      
!      d_yy(i,j,k) = at*abs_vy(i,j,k) + (al-at)*v_yy(i,j,k)*v_yy(i,j,k)/abs_vy(i,j,k)
!      else 
!      d_yy(i,j,k) = 0.0 
!      endif
!      enddo
!      enddo
!      enddo
!      
!      
!      do i=is,ie
!      do j=js,je
!      do k=ks,ke+1
!      
!      if (abs_vz(i,j,k) /= 0.0) THEN
!      
!      d_zz(i,j,k) = at*abs_vz(i,j,k) + (al-at)*v_zz(i,j,k)*v_zz(i,j,k)/abs_vz(i,j,k)
!      else 
!      d_zz(i,j,k) = 0.0 
!      endif
!      enddo
!      enddo
!      enddo
!    
!    
!      do k=ks,ke
!      do j=js,je
!      do i=is,ie
!       sn(i,j,k) = sn(i,j,k) + dt*(d_xx(i,j,k)*(sn(i-1,j,k)-sn(i,j,k))/(dy*dz) - d_xx(i+1,j,k)*(sn(i,j,k)-sn(i+1,j,k))/(dy*dz) + &
!       d_yy(i,j,k)*(sn(i,j-1,k)-sn(i,j,k))/(dx*dz) - d_yy(i,j+1,k)*(sn(i,j,k)-sn(i,j+1,k))/(dx*dz) + &
!       d_zz(i,j,k)*(sn(i,j,k-1)-sn(i,j,k))/(dx*dy) - d_zz(i,j,k+1)*(sn(i,j,k)-sn(i,j,k+1))/(dx*dy))
!      enddo
!      enddo
!      enddo
!      end subroutine disperse

      real(selected_real_kind(15)) function mclimit(theta)
      implicit none
      real(selected_real_kind(15)) theta   
      mclimit = max(0.0d0,min((1.0d0 + theta)/2.0d0,2.0d0,2.0d0*theta))
      end function mclimit 
  
      real(selected_real_kind(15)) function transvel(a,b)
      implicit none
      real(selected_real_kind(15)) a,b,prod,avg 
      prod = a*b
      if (prod .gt. 0.0d0) then
      avg = (abs(a) + abs(b))/2.0d0
      transvel = sign(avg,a)
      else
      transvel = 0.0d0
      endif
      end function transvel 
  
  
!      real(selected_real_kind(15)) function minmod4(a,b,c,d)
!      implicit none
!      real(selected_real_kind(15)) a,b,c,d,one
!      one = 1.0d0
!      minmod4 = (1.0d0/2.0d0)*(sign(one,a)+sign(one,b))*(1.0d0/2.0d0)*&
!      (sign(one,a)+sign(one,c))*(1.0d0/2.0d0)*(sign(one,a)+sign(one,d)) &
!      * min(abs(a),abs(b),abs(c),abs(d))   
!      end function minmod4
!  
!      real(selected_real_kind(15)) function minmod2(a,b)
!      implicit none
!      real(selected_real_kind(15)) a,b,one
!      one = 1.0d0
!      minmod2 = (1.0d0/2.0d0)*(sign(one,a)+sign(one,b))*min(abs(a),abs(b))   
!      end function minmod2
!  
!      real(selected_real_kind(15)) function median(a,b,c)
!      implicit none
!      real(selected_real_kind(15)) a,b,c,minmod2
!      median = a + minmod2(b-a,c-a)    
!      end function median
  





