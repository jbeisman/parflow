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

!---------------------------------------------------------------------
!    advect:
!    Godunov advection routine, Lax-Wendroff fluxes with
!    monotonized-centered flux limiter + multi-dimensional
!    limiter to ensure min-max adherence
!---------------------------------------------------------------------

      subroutine advect(s,sn,uedge,vedge,wedge,phi, &
          dlo,dhi,hx,dt,order,old_sat,sat, &
          iteration,num_iterations,gnx,gny,gnz, &
          gx,gy,gz,fx,fy,fz,vx,wx,uy,wy,uz,vz, &
          stemp,smin,smax,sx,sy,sz,sxtemp,sytemp,sztemp)

      implicit none
      integer, parameter :: dp = selected_real_kind(15)
      integer dlo(3), dhi(3)
      real(dp) hx(3), dt
      integer order
      integer iteration,num_iterations
      integer gnx,gny,gnz,gx,gy,gz

      real(dp) s(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(dp) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(dp) stemp(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(dp) uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(dp) vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(dp) wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(dp) phi(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2)
      real(dp) fx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) fy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) fz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) vx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) wx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp) uy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) wy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) uz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) vz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)  
      real(dp) smin(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2) 
      real(dp) smax(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+2) 
      real(dp) sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)  
      real(dp) sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) sxtemp(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) sytemp(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) sztemp(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) old_sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(dp) sat(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 

 
      integer i,j,k
      integer is,ie,js,je,ks,ke
      integer ks1,ks2,ke1,ke2,js1,js2,je1,je2
      integer is1,is2,ie1,ie2
      integer ii,ii2,ii3,jj,jj2,jj3,kk,kk2,kk3
      integer iiz,iiz2,jjx,jjx2,kky,kky2,iiy,jjz,kkx
      real(dp) dx,dy,dz,dx_inv,dy_inv,dz_inv
      real(dp) half,third
      real(dp) rx,ry,rz,mclimit,limx,thetax,thetay,thetaz,limy,limz
      real(dp) transvel,sat_diff,iter,num_iter,num_iter_inv
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

      iter=DBLE(iteration)
      num_iter=DBLE(num_iterations)
      num_iter_inv = 1.0_dp/num_iter

     !! initialization
      stemp = s

      fx = 0.0_dp
      fy = 0.0_dp
      fz = 0.0_dp
      
      sx = 0.0_dp
      sy = 0.0_dp
      sz = 0.0_dp
  
      vx = 0.0_dp
      wx = 0.0_dp
      uy = 0.0_dp
      wy = 0.0_dp
      uz = 0.0_dp
      vz = 0.0_dp


      sxtemp = 0.0_dp
      sytemp = 0.0_dp
      sztemp = 0.0_dp


      !! make smaller and faster for lower D problems
      if((ie-is) == 0) then
        is1 = is+1
        ie1 = ie-2
        is2 = is+1
        ie2 = ie-1
      else
        is1 = is
        ie1 = ie
        is2 = is
        ie2 = ie
      endif

      if((je-js) == 0) then
        js1 = js+1
        je1 = je-2
        js2 = js+1
        je2 = je-1
      else
        js1 = js
        je1 = je
        js2 = js
        je2 = je
      endif

      if((ke-ks) == 0) then
        ks1 = ks+1
        ke1 = ke-2
        ks2 = ks+1
        ke2 = ke-1
      else
        ks1 = ks
        ke1 = ke
        ks2 = ks
        ke2 = ke
      endif

      !! main loop
      !! compute increment fluxes and second order corrections
      !! compute transverse riemann problem, move fluxes appropriately
      do k=ks1-1,ke1+2
        do j=js1-1,je1+2
          do i=is1-1,ie1+2
     
          if (uedge(i,j,k) .ge. 0.0)then
            ii =  i-1
            ii2 = i
            ii3 = i-1
          else
            ii =  i
            ii2 = i-1
            ii3 = i+1
          endif
    
          if (vedge(i,j,k) .ge. 0.0)then
            jj =  j-1
            jj2 = j
            jj3 = j-1
          else
            jj =  j
            jj2 = j-1
            jj3 = j+1
          endif

          if (wedge(i,j,k) .ge. 0.0)then
            kk =  k-1
            kk2 = k
            kk3 = k-1
          else
            kk =  k
            kk2 = k-1
            kk3 = k+1
          endif

          smin(i,j,k) = min(s(i,j,k),s(ii,j,k),s(i,jj,k),s(i,j,kk))
          smax(i,j,k) = max(s(i,j,k),s(ii,j,k),s(i,jj,k),s(i,j,kk))
    
          !these values give a first order scheme
          fx(i,j,k)=fx(i,j,k)+uedge(i,j,k)*s(ii,j,k) 
          fy(i,j,k)=fy(i,j,k)+vedge(i,j,k)*s(i,jj,k) 
          fz(i,j,k)=fz(i,j,k)+wedge(i,j,k)*s(i,j,kk)
    
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
          sytemp(ii2,jjx,k) = sytemp(ii2,jjx,k)-half*rx*uedge(i,j,k)*vx(i,j,k)
          sztemp(ii2,j,kkx)=sztemp(ii2,j,kkx) - half*rx*uedge(i,j,k)*wx(i,j,k) + &
          (third)*uedge(i,j,k)*abs(vx(i,j,k))*wx(i,j,k)*rx 
          sztemp(ii2,jjx2,kkx)=sztemp(ii2,jjx2,kkx) - (third)*uedge(i,j,k) * &
          abs(vx(i,j,k))*wx(i,j,k)*rx 
    
          !y sweep
          sztemp(i,jj2,kky) = sztemp(i,jj2,kky)-half*ry*vedge(i,j,k)*wy(i,j,k)
          sxtemp(iiy,jj2,k)=sxtemp(iiy,jj2,k) - half*ry*vedge(i,j,k)*uy(i,j,k) + &
          (third)*vedge(i,j,k)*abs(wy(i,j,k))*uy(i,j,k)*ry 
          sxtemp(iiy,jj2,kky2)=sxtemp(iiy,jj2,kky2) - (third)*vedge(i,j,k) * &
          abs(wy(i,j,k))*uy(i,j,k)*ry 
    
          !z sweep
          sxtemp(iiz,j,kk2) = sxtemp(iiz,j,kk2)-half*rz*wedge(i,j,k)*uz(i,j,k)
          sytemp(i,jjz,kk2)=sytemp(i,jjz,kk2) - half*rz*wedge(i,j,k)*vz(i,j,k) + &
          (third)*wedge(i,j,k)*abs(uz(i,j,k))*vz(i,j,k)*rz 
          sytemp(iiz2,jjz,kk2)=sytemp(iiz2,jjz,kk2) - (third)*wedge(i,j,k) * &
          abs(uz(i,j,k))*vz(i,j,k)*rz
    
          ! Add second order transverse propogation information
          !x sweep
          sytemp(i,jjx,k)      = sytemp(i,jjx,k) + vx(i,j,k)*sx(i,j,k)
          sytemp(i-1,jjx,k)    = sytemp(i-1,jjx,k) - vx(i,j,k)*sx(i,j,k)
          sztemp(i,j,kkx)      = sztemp(i,j,kkx) + (1.0_dp - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
          sztemp(i,jjx2,kkx)   = sztemp(i,jjx2,kkx) + abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)
          sztemp(i-1,j,kkx)    = sztemp(i-1,j,kkx) - (1.0_dp - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
          sztemp(i-1,jjx2,kkx) = sztemp(i-1,jjx2,kkx) - abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)
    
          !y sweep
          sztemp(i,j,kky)      = sztemp(i,j,kky) + wy(i,j,k)*sy(i,j,k)
          sztemp(i,j-1,kky)    = sztemp(i,j-1,kky) - wy(i,j,k)*sy(i,j,k)
          sxtemp(iiy,j,k)      = sxtemp(iiy,j,k) + (1.0_dp - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
          sxtemp(iiy,j,kky2)   = sxtemp(iiy,j,kky2) + abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)
          sxtemp(iiy,j-1,k)    = sxtemp(iiy,j-1,k) - (1.0_dp - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
          sxtemp(iiy,j-1,kky2) = sxtemp(iiy,j-1,kky2) - abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)
    
          !z sweep
          sxtemp(iiz,j,k)      = sxtemp(iiz,j,k) + uz(i,j,k)*sz(i,j,k)
          sxtemp(iiz,j,k-1)    = sxtemp(iiz,j,k-1) - uz(i,j,k)*sz(i,j,k)
          sytemp(i,jjz,k)      = sytemp(i,jjz,k) + (1.0_dp - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
          sytemp(iiz2,jjz,k)   = sytemp(iiz2,jjz,k) + abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k)
          sytemp(i,jjz,k-1)    = sytemp(i,jjz,k-1) - (1.0_dp - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
          sytemp(iiz2,jjz,k-1) = sytemp(iiz2,jjz,k-1) - abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k) 
          enddo
        enddo
      enddo
    
    
      !! first-order transport 
      do k=ks2-1,ke2+1
        do j=js2-1,je2+1
          do i=is2-1,ie2+1

            sat_diff=(sat(i,j,k) - old_sat(i,j,k)) * num_iter_inv

            stemp(i,j,k)= ((iter*sat_diff + old_sat(i,j,k))*phi(i,j,k)*s(i,j,k) + &
            ((dt*dx_inv)*(fx(i,j,k) - fx(i+1,j,k)) + (dt*dy_inv)*(fy(i,j,k)-fy(i,j+1,k)) + &
            (dt*dz_inv)*(fz(i,j,k)-fz(i,j,k+1)))) / (((iter+1.0_dp)*sat_diff + old_sat(i,j,k))*phi(i,j,k))

          enddo
        enddo
      enddo

     
      sn(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)) = stemp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

     
      if (order == 2) then

        !! remove high order and transverse waves at
        !! global boundary
        if (is == gx) then
          sxtemp(is,js:je+1,ks:ke+1) = 0.0_dp
          sytemp(is,js:je+1,ks:ke+1) = 0.0_dp
          sztemp(is,js:je+1,ks:ke+1) = 0.0_dp
        endif

        if (js == gy) then
          sxtemp(is:ie+1,js,ks:ke+1) = 0.0_dp
          sytemp(is:ie+1,js,ks:ke+1) = 0.0_dp
          sztemp(is:ie+1,js,ks:ke+1) = 0.0_dp
        endif

        if (ks == gz) then
          sxtemp(is:ie+1,js:je+1,ks) = 0.0_dp
          sytemp(is:ie+1,js:je+1,ks) = 0.0_dp
          sztemp(is:ie+1,js:je+1,ks) = 0.0_dp
        endif

        if (ie == gnx) then
          sxtemp(ie,js:je+1,ks:ke+1) = 0.0_dp
          sytemp(ie,js:je+1,ks:ke+1) = 0.0_dp
          sztemp(ie,js:je+1,ks:ke+1) = 0.0_dp
        endif

        if (je == gny) then
          sxtemp(is:ie+1,je,ks:ke+1) = 0.0_dp
          sytemp(is:ie+1,je,ks:ke+1) = 0.0_dp
          sztemp(is:ie+1,je,ks:ke+1) = 0.0_dp
        endif

        if (ke == gnz) then
          sxtemp(is:ie+1,js:je+1,ke) = 0.0_dp
          sytemp(is:ie+1,js:je+1,ke) = 0.0_dp
          sztemp(is:ie+1,js:je+1,ke) = 0.0_dp
        endif

        !!gather 2nd order anti-diffusive correction and transverse terms
        sx=sx+sxtemp
        sy=sy+sytemp
        sz=sz+sztemp

        !!call limiter to enforce min/max
        call limit(smax,smin,stemp,sx,sy,sz,fx,fy,fz,&
          dlo,dhi,dt,dx_inv,dy_inv,dz_inv,vx,wx,uy,wy,uz,vz)

        do k=ks,ke
          do j=js,je
            do i=is,ie

              sat_diff=(sat(i,j,k)-old_sat(i,j,k)) * num_iter_inv

              sn(i,j,k)=sn(i,j,k) + ((dt*dx_inv)*(fx(i,j,k)*sx(i,j,k) - fx(i+1,j,k)*sx(i+1,j,k)) &
              + (dt*dy_inv)*(fy(i,j,k)*sy(i,j,k)-fy(i,j+1,k)*sy(i,j+1,k)) + &
              (dt*dz_inv)*(fz(i,j,k)*sz(i,j,k)-fz(i,j,k+1)*sz(i,j,k+1))) / &
              (((iter+1.0_dp)*sat_diff + old_sat(i,j,k))*phi(i,j,k))

            enddo
          enddo
        enddo

      endif


      do k=ks,ke
        do j=js,je
          do i=is,ie  

            !!cutoff any values that violate min/max
           if (sn(i,j,k) .lt. smin(i,j,k)) sn(i,j,k) = smin(i,j,k)
           if (sn(i,j,k) .gt. smax(i,j,k)) sn(i,j,k) = smax(i,j,k)  

          enddo
        enddo
      enddo 
       
      !! needs work
      !!call disperse(sn,uedge,vedge,wedge,lo,hi,dlo,dhi,hx,dt)

      return
      end subroutine advect

    

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
  


      subroutine limit(smax,smin,sn,sx,sy,sz,cx,cy,cz,dlo,dhi,dt,dx_inv,&
                       dy_inv,dz_inv,p_plus,p_minus,q_plus,q_minus,r_plus,r_minus)
      implicit none
      integer, parameter :: dp = selected_real_kind(15)
      integer  dlo(3), dhi(3)
      integer  is,ie,js,je,ks,ke,i,j,k 
      real(dp) dt,dx_inv,dy_inv,dz_inv
      real(dp) smin(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) smax(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) sx(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3) 
      real(dp) sy(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)
      real(dp) sz(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+3)
      real(dp) cx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) cy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) cz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) p_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) p_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp) q_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) q_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp) r_plus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(dp) r_minus(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(dp) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)

      !! clear memory
      cx = 0.0_dp
      cy = 0.0_dp
      cz = 0.0_dp
      p_plus = 0.0_dp
      p_minus = 0.0_dp
      q_plus = 0.0_dp
      q_minus = 0.0_dp
      r_plus = 0.0_dp
      r_minus = 0.0_dp

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

            q_plus(i,j,k) = max(smax(i,j,k),maxval(sn(i-1:i+1,j-1:j+1,k-1:k+1))) - sn(i,j,k)

            q_minus(i,j,k) = sn(i,j,k) - min(smin(i,j,k),minval(sn(i-1:i+1,j-1:j+1,k-1:k+1)))

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
       
            if (sx(i,j,k) .ge. 0.0_dp) then
              cx(i,j,k) = min(r_plus(i,j,k),r_minus(i-1,j,k))
            else
              cx(i,j,k) = min(r_plus(i-1,j,k),r_minus(i,j,k))
            endif
            
            if (sy(i,j,k) .ge. 0.0_dp) then
              cy(i,j,k) = min(r_plus(i,j,k),r_minus(i,j-1,k))
            else
              cy(i,j,k) = min(r_plus(i,j-1,k),r_minus(i,j,k))
            endif
      
            if (sz(i,j,k) .ge. 0.0_dp) then
              cz(i,j,k) = min(r_plus(i,j,k),r_minus(i,j,k-1))
            else
              cz(i,j,k) = min(r_plus(i,j,k-1),r_minus(i,j,k))
            endif

          enddo
        enddo
      enddo
   
      return
      end subroutine limit


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




