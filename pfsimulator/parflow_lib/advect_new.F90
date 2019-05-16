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
!  for the GNU Lesser General Publi!License.
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

!*********************************************************************
!    advect, slopexy, slopez
!    
!    Godunov advection routine
!    
!*********************************************************************




!---------------------------------------------------------------------
!    advect:
!    Godunov advection routine
!---------------------------------------------------------------------

      subroutine advect(s,sn,uedge,vedge,wedge,phi, &
          lo,hi,dlo,dhi,hx,dt,order, &
          old_sat,sat,iteration,num_iterations) 
      implicit none
    
!    ::: argument declarations

      integer lo(3), hi(3)
      integer dlo(3), dhi(3)
      real(selected_real_kind(8))  hx(3), dt
      integer order
      integer iteration,num_iterations

      real(selected_real_kind(8)) s(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(selected_real_kind(8)) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(selected_real_kind(8)) uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(selected_real_kind(8)) phi(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2)
      real(selected_real_kind(8)) cx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) cy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) cz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) fx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) fy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) fz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) vx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) wx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(selected_real_kind(8)) uy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) wy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) uz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) vz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)  
      real(selected_real_kind(8)) smin(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) smax(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) sx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)  
      real(selected_real_kind(8)) sy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) sz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) sxtemp(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) sytemp(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) sztemp(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) old_sat(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 
      real(selected_real_kind(8)) sat(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3) 

 
      integer is,ie,js,je,ks,ke,n,ii,jj,kk,ii2,jj2,kk2,ii3,jj3,kk3,ii4,jj4,kk4,ii5,jj5,kk5
      integer i,j,k,z,ii6,jj6,kk6,jjx,jjx2,jjx3,jjx4,jjx5,jjx6,kkx,kkx2,kkx3,kkx4,kkx5,kkx6
      integer iiy,iiy2,iiy3,iiy4,iiy5,iiy6,kky,kky2,kky3,kky4,kky5,kky6
      integer iiz,iiz2,iiz3,iiz4,iiz5,iiz6,jjz,jjz2,jjz3,jjz4,jjz5,jjz6,t
      real(selected_real_kind(8)) dx,dy,dz
      real(selected_real_kind(8)) half
      real(selected_real_kind(8)) rx,ry,rz,mclimit,limx,thetax,thetay,thetaz,limy,limz
      real(selected_real_kind(8)) transvel,sat_diff,iter,num_iter
      
      real(selected_real_kind(8)) minmod4,minmod2,median
      integer flag

      print*,"inside advection"
     
      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)
      half = 0.5d0
      iter=real(iteration)
      num_iter=real(num_iterations)
   
!!initialize some stuff
print*,"initiailizing"

     fx = 0.0d0
     fy = 0.0d0
     fz = 0.0d0

     cx = 0.0d0 
     cy = 0.0d0
     cz = 0.0d0
     
     sx = 0.0d0
     sy = 0.0d0
     sz = 0.0d0

     vx = 0.0d0
     wx = 0.0d0
     uy = 0.0d0
     wy = 0.0d0
     uz = 0.0d0
     vz = 0.0d0

     sxtemp = 0.0d0
     sytemp = 0.0d0
     sztemp = 0.0d0

    do i=is-2,ie+2
    do j=js-2,je+2
    do k=ks-2,ke+2
     smin(i,j,k) = minval(s(i-1:i+1,j-1:j+1,k-1:k+1))
     smax(i,j,k) = maxval(s(i-1:i+1,j-1:j+1,k-1:k+1))
    enddo
    enddo
    enddo    
!  do i=is-2,ie+3
!    do j=js,je
!    do k=ks,ke
!    print*,i,j,k,uedge(i,j,k),s(i,j,k)
!    print*,vedge(i,j,k)
!    print*,wedge(i,j,k)
!  enddo
!enddo
!enddo
   print*,"done init"


  
  
      flag = 0
     do k=ks-1,ke+2
     do j=js-1,je+2
     do i=is-1,ie+2

    !! uedge(i,j,k) = 0.25
    !! vedge(i,j,k) = 0.0
    !! wedge(i,j,k) = 0.0
       
     if (uedge(i,j,k) .ge. 0.0)then
     ii =  i-1
     ii2 = i
     ii3 = i+1
     ii4 = i+1
     ii5 = i+1
     ii6 = i-1
     else
     ii =  i
     ii2 = i-1
     ii3 = i
     ii4 = i
     ii5 = i-1
     ii6 = i+1
     endif

    if (vedge(i,j,k) .ge. 0.0)then
     jj =  j-1
     jj2 = j
     jj3 = j+1
     jj4 = j+1
     jj5 = j+1
     jj6 = j-1
     else
     jj =  j
     jj2 = j-1
     jj3 = j
     jj4 = j
     jj5 = j-1
     jj6 = j+1
     endif


     if (wedge(i,j,k) .ge. 0.0)then
     kk =  k-1
     kk2 = k
     kk3 = k+1
     kk4 = k+1
     kk5 = k+1
     kk6 = k-1
     else
     kk =  k
     kk2 = k-1
     kk3 = k
     kk4 = k
     kk5 = k-1
     kk6 = k+1
     endif


    if (flag == 0) print*,"before transvel"
     !! Transverse velocities
     vx(i,j,k) = (dt/dx)*transvel(vedge(ii2,j,k),vedge(ii2,j+1,k)) 
     wx(i,j,k) = (dt/dx)*transvel(wedge(ii2,j,k),wedge(ii2,j,k+1)) 
     uy(i,j,k) = (dt/dy)*transvel(uedge(i,jj2,k),uedge(i+1,jj2,k)) 
     wy(i,j,k) = (dt/dy)*transvel(wedge(i,jj2,k),wedge(i,jj2,k+1))
     uz(i,j,k) = (dt/dz)*transvel(uedge(i,j,kk2),uedge(i+1,j,kk2))
     vz(i,j,k) = (dt/dz)*transvel(vedge(i,j,kk2),vedge(i,j+1,kk2))
if (flag == 0) print*,"after transvel"

!! not all of these logic flags are necessary anymore
   if (vx(i,j,k) .gt. 0.0)then
    jjx =  j-1
    jjx2 = j
    jjx3 = j+1
    jjx4 = j+1
    jjx5 = j+1
    jjx6 = j-1
    else
    jjx =  j
    jjx2 = j-1
    jjx3 = j
    jjx4 = j
    jjx5 = j-1
    jjx6 = j+1
    endif

   if (wx(i,j,k) .gt. 0.0)then
    kkx =  k-1
    kkx2 = k
    kkx3 = k+1
    kkx4 = k+1
    kkx5 = k+1
    kkx6 = k-1
    else
    kkx =  k
    kkx2 = k-1
    kkx3 = k
    kkx4 = k
    kkx5 = k-1
    kkx6 = k+1
    endif


   if (uy(i,j,k) .gt. 0.0)then
    iiy =  i-1
    iiy2 = i
    iiy3 = i+1
    iiy4 = i+1
    iiy5 = i+1
    iiy6 = i-1
    else
    iiy =  i
    iiy2 = i-1
    iiy3 = i
    iiy4 = i
    iiy5 = i-1
    iiy6 = i+1
    endif

    if (wy(i,j,k) .gt. 0.0)then
    kky =  k-1
    kky2 = k
    kky3 = k+1
    kky4 = k+1
    kky5 = k+1
    kky6 = k-1
    else
    kky =  k
    kky2 = k-1
    kky3 = k
    kky4 = k
    kky5 = k-1
    kky6 = k+1
    endif


    if (uz(i,j,k) .gt. 0.0)then
    iiz =  i-1
    iiz2 = i
    iiz3 = i+1
    iiz4 = i+1
    iiz5 = i+1
    iiz6 = i-1
    else
    iiz =  i
    iiz2 = i-1
    iiz3 = i
    iiz4 = i
    iiz5 = i-1
    iiz6 = i+1
    endif

    if (vz(i,j,k) .gt. 0.0)then
    jjz =  j-1
    jjz2 = j
    jjz3 = j+1
    jjz4 = j+1
    jjz5 = j+1
    jjz6 = j-1
    else
    jjz =  j
    jjz2 = j-1
    jjz3 = j
    jjz4 = j
    jjz5 = j-1
    jjz6 = j+1
    endif

if (flag == 0) print*,"before fluxes"
     !these values give a first order scheme
     fx(i,j,k)=fx(i,j,k)+uedge(i,j,k)*s(ii,j,k) 
     fy(i,j,k)=fy(i,j,k)+vedge(i,j,k)*s(i,jj,k) 
     fz(i,j,k)=fz(i,j,k)+wedge(i,j,k)*s(i,j,kk)

   if (order == 2) then  
     !concentration gradient
     rx=s(i,j,k)-s(i-1,j,k)
     ry=s(i,j,k)-s(i,j-1,k)
     rz=s(i,j,k)-s(i,j,k-1)


    ! these sweeps add transverse propogation information to the first order upwind method 
    !x sweep
!    sytemp(ii2,jjx3,k) = sytemp(ii2,jjx3,k)-half*rx*uedge(i,j,k)*vx(i,j,k)
!    sztemp(ii2,j,kkx4)=sztemp(ii2,j,kkx4) - half*rx*uedge(i,j,k)*wx(i,j,k) + (1.0/3.0)*uedge(i,j,k)*abs(vx(i,j,k))*wx(i,j,k)*rx 
!    sztemp(ii2,jjx5,kkx4)=sztemp(ii2,jjx5,kkx4) - (1.0/3.0)*uedge(i,j,k)*abs(vx(i,j,k))*wx(i,j,k)*rx 


!   !y sweep
!    sztemp(i,jj2,kky3) = sztemp(i,jj2,kky3)-half*ry*vedge(i,j,k)*wy(i,j,k)
!    sxtemp(iiy4,jj2,k)=sxtemp(iiy4,jj2,k) - half*ry*vedge(i,j,k)*uy(i,j,k) + (1.0/3.0)*vedge(i,j,k)*abs(wy(i,j,k))*uy(i,j,k)*ry 
!    sxtemp(iiy4,jj2,kky5)=sxtemp(iiy4,jj2,kky5) - (1.0/3.0)*vedge(i,j,k)*abs(wy(i,j,k))*uy(i,j,k)*ry 

!   !z sweep
!    sxtemp(iiz3,j,kk2) = sxtemp(iiz3,j,kk2)-half*rz*wedge(i,j,k)*uz(i,j,k)
!    sytemp(i,jjz4,kk2)=sytemp(i,jjz4,kk2) - half*rz*wedge(i,j,k)*vz(i,j,k) + (1.0/3.0)*wedge(i,j,k)*abs(uz(i,j,k))*vz(i,j,k)*rz 
!    sytemp(iiz5,jjz4,kk2)=sytemp(iiz5,jjz4,kk2) - (1.0/3.0)*wedge(i,j,k)*abs(uz(i,j,k))*vz(i,j,k)*rz 
    
     if (flag == 0) print*,"after s_temp assignement"

     !second order - LW flux - 
     
     thetax=(s(ii6,j,k)-s(ii6-1,j,k))/(s(i,j,k)-s(i-1,j,k))
     limx=mclimit(thetax)
     sx(i,j,k)=half*abs(uedge(i,j,k))*(1.0-(dt/dx)*abs(uedge(i,j,k)))*rx*limx
    

     thetay=(s(i,jj6,k)-s(i,jj6-1,k))/(s(i,j,k)-s(i,j-1,k))
     limy=mclimit(thetay)
     sy(i,j,k) = half*abs(vedge(i,j,k))*(1.0-(dt/dy)*abs(vedge(i,j,k)))*ry*limy
    
     
   
     thetaz=(s(i,j,kk6)-s(i,j,kk6-1))/(s(i,j,k)-s(i,j,k-1))
     limz=mclimit(thetaz)
     sz(i,j,k)= half*abs(wedge(i,j,k))*(1.0-(dt/dz)*abs(wedge(i,j,k)))*rz*limz

     if (flag == 0) print*,"after 2nd order fluxes"
    
     




     ! Add second order transverse propogation information
    !x sweep
!    sytemp(i,jjx3,k)      = sytemp(i,jjx3,k) + vx(i,j,k)*sx(i,j,k)
!    sytemp(i-1,jjx3,k)    = sytemp(i-1,jjx3,k) - vx(i,j,k)*sx(i,j,k)
!    sztemp(i,j,kkx4)      = sztemp(i,j,kkx4) + (1.0 - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
!    sztemp(i,jjx5,kkx4)   = sztemp(i,jjx5,kkx4) + abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)
!    sztemp(i-1,j,kkx4)    = sztemp(i-1,j,kkx4) - (1.0 - abs(vx(i,j,k)))*wx(i,j,k)*sx(i,j,k)
!    sztemp(i-1,jjx5,kkx4) = sztemp(i-1,jjx5,kkx4) - abs(vx(i,j,k))*wx(i,j,k)*sx(i,j,k)


!   !y sweep
!    sztemp(i,j,kky3)      = sztemp(i,j,kky3) + wy(i,j,k)*sy(i,j,k)
!    sztemp(i,j-1,kky3)    = sztemp(i,j-1,kky3) - wy(i,j,k)*sy(i,j,k)
!    sxtemp(iiy4,j,k)      = sxtemp(iiy4,j,k) + (1.0 - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
!    sxtemp(iiy4,j,kky5)   = sxtemp(iiy4,j,kky5) + abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)
!    sxtemp(iiy4,j-1,k)    = sxtemp(iiy4,j-1,k) - (1.0 - abs(wy(i,j,k)))*uy(i,j,k)*sy(i,j,k)
!    sxtemp(iiy4,j-1,kky5) = sxtemp(iiy4,j-1,kky5) - abs(wy(i,j,k))*uy(i,j,k)*sy(i,j,k)


!   !z sweep
!    sxtemp(iiz3,j,k)      = sxtemp(iiz3,j,k) + uz(i,j,k)*sz(i,j,k)
!    sxtemp(iiz3,j,k-1)    = sxtemp(iiz3,j,k-1) - uz(i,j,k)*sz(i,j,k)
!    sytemp(i,jjz4,k)      = sytemp(i,jjz4,k) + (1.0 - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
!    sytemp(iiz5,jjz4,k)   = sytemp(iiz5,jjz4,k) + abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k)
!    sytemp(i,jjz4,k-1)    = sytemp(i,jjz4,k-1) - (1.0 - abs(uz(i,j,k)))*vz(i,j,k)*sz(i,j,k)
!    sytemp(iiz5,jjz4,k-1) = sytemp(iiz5,jjz4,k-1) - abs(uz(i,j,k))*vz(i,j,k)*sz(i,j,k)
     if (flag == 0) print*,"after s_temp 2nd order"
  
  endif
flag = 1
    enddo
    enddo
    enddo

  print*,"done with main loop"  
 

    
         !! first-order transport 
     do k=ks,ke
     do j=js,je
     do i=is,ie  
     sat_diff=(sat(i,j,k)-old_sat(i,j,k))/num_iter
     sn(i,j,k)=((iter*sat_diff + old_sat(i,j,k))*phi(i,j,k)*s(i,j,k) + ((dt/dx)*(fx(i,j,k) - fx(i+1,j,k)) + &
     (dt/dy)*(fy(i,j,k)-fy(i,j+1,k))+ (dt/dz)*(fz(i,j,k)-fz(i,j,k+1)))) & 
     /(((iter+1.0)*sat_diff + old_sat(i,j,k))*phi(i,j,k))

     if (sn(i,j,k) .lt. smin(i,j,k)) sn(i,j,k) = smin(i,j,k)
     if (sn(i,j,k) .gt. smax(i,j,k)) sn(i,j,k) = smax(i,j,k)
    
    
    
    
     enddo
     enddo
     enddo
     
     
  if (order == 2) then
  
  
  
  !!Accumulate 2nd order anti-diffusive correction and transverse terms
    sx=sx+sxtemp
    sy=sy+sytemp
    sz=sz+sztemp
    print*,"before limit"
   !!call limiter to enforce min/max
    !!call limit(smax,smin,sn,sx,sy,sz,cx,cy,cz,lo,hi,dlo,dhi,dt,dx,dy,dz)
        print*,"after limit"
    cx = 1.0d0
    cy = 1.0d0
    cz = 1.0d0 
      
         
    do k=ks,ke
    do j=js,je
    do i=is,ie
   
   sat_diff=(sat(i,j,k)-old_sat(i,j,k))/num_iter
  
   sn(i,j,k)=sn(i,j,k) + ((dt/dx)*(cx(i,j,k)*sx(i,j,k) - cx(i+1,j,k)*sx(i+1,j,k)) & 
   + (dt/dy)*(cy(i,j,k)*sy(i,j,k)-cy(i,j+1,k)*sy(i,j+1,k)) + & 
   (dt/dz)*(cz(i,j,k)*sz(i,j,k)-cz(i,j,k+1)*sz(i,j,k+1)))/(((iter+1.0)*sat_diff + old_sat(i,j,k))*phi(i,j,k))
  


     if (sn(i,j,k) .lt. smin(i,j,k)) sn(i,j,k) = smin(i,j,k)
     if (sn(i,j,k) .gt. smax(i,j,k)) sn(i,j,k) = smax(i,j,k)

     enddo
     enddo
     enddo
     

endif 
 
!!call disperse(sn,uedge,vedge,wedge,lo,hi,dlo,dhi,hx,dt) 


     return
     end subroutine advect 

    

!! dispersion calculations

      subroutine disperse(sn,uedge,vedge,wedge,lo,hi,dlo,dhi,hx,dt) 
      implicit none
      
            
      integer lo(3), hi(3),i,j,k,ie,is,je,js,ke,ks
      integer dlo(3), dhi(3)
      real(selected_real_kind(8))  hx(3), dt, al, at,dx,dy,dz
      
      real(selected_real_kind(8)) uedge(dlo(1)-2:dhi(1)+3,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) vedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+3,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) wedge(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+3) 
      real(selected_real_kind(8)) v_xx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_xy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_xz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) v_yy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_yx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_yz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) v_zz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_zx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) v_zy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)
      real(selected_real_kind(8)) abs_vx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) abs_vy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) abs_vz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) d_xx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_xy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_xz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) d_yy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_yx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_yz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) d_zz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_zx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) d_zy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      real(selected_real_kind(8)) dc_x(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) dc_y(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) dc_z(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1)
      
      
      
      al = .00625
      at = .000625
      
      
      
      
      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
      dx = hx(1)
      dy = hx(2)
      dz = hx(3)

!velocities
      do k=ks,ke+1
      do j=js,je+1
      do i=is,ie+1
      
      d_xx(i,j,k) = 0.0
      d_xy(i,j,k) = 0.0
      d_xz(i,j,k) = 0.0
      
      d_yy(i,j,k) = 0.0
      d_yx(i,j,k) = 0.0
      d_yz(i,j,k) = 0.0
      
      d_zz(i,j,k) = 0.0
      d_zx(i,j,k) = 0.0
      d_zy(i,j,k) = 0.0
      
     
      
     
     v_xx(i,j,k)=uedge(i,j,k)
     v_xy(i,j,k) = sum(vedge(i-1:i,j:j+1,k))/4.0
     v_xz(i,j,k) = sum(wedge(i-1:i,j,k:k+1))/4.0
    
     v_yy(i,j,k)=vedge(i,j,k) 
     v_yx(i,j,k) = sum(uedge(i:i+1,j-1:j,k))/4.0
     v_yz(i,j,k) = sum(wedge(i,j-1:j,k:k+1))/4.0
     
     v_zz(i,j,k)=wedge(i,j,k)
     v_zx(i,j,k) = sum(uedge(i:i+1,j,k-1:k))/4.0
     v_zy(i,j,k) = sum(vedge(i,j:j+1,k-1:k))/4.0
     
     abs_vx(i,j,k) = (v_xx(i,j,k)**2 + v_xy(i,j,k)**2 + v_xz(i,j,k)**2)**0.5
     abs_vy(i,j,k) = (v_yy(i,j,k)**2 + v_yx(i,j,k)**2 + v_yz(i,j,k)**2)**0.5
     abs_vz(i,j,k) = (v_zz(i,j,k)**2 + v_zx(i,j,k)**2 + v_zy(i,j,k)**2)**0.5 
     
     enddo
     enddo
     enddo
      
     
      do i=is,ie+1
      do j=js,je
      do k=ks,ke
      
      if (abs_vx(i,j,k) /= 0.0) THEN
      
      d_xx(i,j,k) = at*abs_vx(i,j,k) + (al-at)*v_xx(i,j,k)*v_xx(i,j,k)/abs_vx(i,j,k)
      else 
      d_xx(i,j,k) = 0.0 
      endif
      enddo
      enddo
      enddo
      
      do i=is,ie
      do j=js,je+1
      do k=ks,ke
      
      if (abs_vy(i,j,k) /= 0.0) THEN
      
      d_yy(i,j,k) = at*abs_vy(i,j,k) + (al-at)*v_yy(i,j,k)*v_yy(i,j,k)/abs_vy(i,j,k)
      else 
      d_yy(i,j,k) = 0.0 
      endif
      enddo
      enddo
      enddo
      
      
      do i=is,ie
      do j=js,je
      do k=ks,ke+1
      
      if (abs_vz(i,j,k) /= 0.0) THEN
      
      d_zz(i,j,k) = at*abs_vz(i,j,k) + (al-at)*v_zz(i,j,k)*v_zz(i,j,k)/abs_vz(i,j,k)
      else 
      d_zz(i,j,k) = 0.0 
      endif
      enddo
      enddo
      enddo
   
  
      do k=ks,ke
      do j=js,je
      do i=is,ie
       sn(i,j,k) = sn(i,j,k) + dt*(d_xx(i,j,k)*(sn(i-1,j,k)-sn(i,j,k))/(dy*dz) - d_xx(i+1,j,k)*(sn(i,j,k)-sn(i+1,j,k))/(dy*dz) + &
       d_yy(i,j,k)*(sn(i,j-1,k)-sn(i,j,k))/(dx*dz) - d_yy(i,j+1,k)*(sn(i,j,k)-sn(i,j+1,k))/(dx*dz) + &
       d_zz(i,j,k)*(sn(i,j,k-1)-sn(i,j,k))/(dx*dy) - d_zz(i,j,k+1)*(sn(i,j,k)-sn(i,j,k+1))/(dx*dy))
    enddo
    enddo
    enddo

  end subroutine disperse



      subroutine limit(smax,smin,sn,sx,sy,sz,cx,cy,cz,lo,hi,dlo,dhi,dt,dx,dy,dz)
      implicit none
      integer lo(3), hi(3)
      integer dlo(3), dhi(3)
      integer is,ie,js,je,ks,ke,i,j,k 
      real(selected_real_kind(8)) dt,dx,dy,dz
      real(selected_real_kind(8)) smin(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) smax(dlo(1)-2:dhi(1)+2,dlo(2)-2:dhi(2)+2,dlo(3)-2:dhi(3)+2) 
      real(selected_real_kind(8)) sx(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2) 
      real(selected_real_kind(8)) sy(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(selected_real_kind(8)) sz(dlo(1)-1:dhi(1)+2,dlo(2)-1:dhi(2)+2,dlo(3)-1:dhi(3)+2)
      real(selected_real_kind(8)) cx(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) cy(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) cz(dlo(1):dhi(1)+1,dlo(2):dhi(2)+1,dlo(3):dhi(3)+1) 
      real(selected_real_kind(8)) p_plus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(selected_real_kind(8)) p_minus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1)
      real(selected_real_kind(8)) q_plus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(selected_real_kind(8)) q_minus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1)
      real(selected_real_kind(8)) r_plus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1) 
      real(selected_real_kind(8)) r_minus(dlo(1)-1:dhi(1)+1,dlo(2)-1:dhi(2)+1,dlo(3)-1:dhi(3)+1)
      real(selected_real_kind(8)) sn(dlo(1)-3:dhi(1)+3,dlo(2)-3:dhi(2)+3,dlo(3)-3:dhi(3)+3)

       print*, "inside limit"
      is = dlo(1)
      ie = dhi(1)
      js = dlo(2)
      je = dhi(2)
      ks = dlo(3)
      ke = dhi(3)
       
      
      do k=ks-1,ke+1
      do j=js-1,je+1
      do i=is-1,ie+1 
            
      p_plus(i,j,k) = (dt/dx)*(max(0.0,sx(i,j,k)) - min(0.0,sx(i+1,j,k))) + (dt/dy)*(max(0.0,sy(i,j,k)) - min(0.0,sy(i,j+1,k))) &
      + (dt/dz)*(max(0.0,sz(i,j,k)) - min(0.0,sz(i,j,k+1)))

      p_minus(i,j,k) = (dt/dx)*(max(0.0,sx(i+1,j,k)) - min(0.0,sx(i,j,k))) + (dt/dy)*(max(0.0,sy(i,j+1,k)) - min(0.0,sy(i,j,k))) &
      + (dt/dz)*(max(0.0,sz(i,j,k+1)) - min(0.0,sz(i,j,k)))

      q_plus(i,j,k) = max(smax(i,j,k),maxval(sn(i-1:i+1,j-1:j+1,k-1:k+1))) - sn(i,j,k)

      q_minus(i,j,k) = sn(i,j,k) - min(smin(i,j,k),minval(sn(i-1:i+1,j-1:j+1,k-1:k+1)))

      if (p_plus(i,j,k) .gt. 0.0) then
      r_plus(i,j,k) = min(q_plus(i,j,k)/p_plus(i,j,k),1.0)
      else
      r_plus(i,j,k) = 0.0
      endif 
   
      if (p_minus(i,j,k) .gt. 0.0) then
      r_minus(i,j,k) = min(q_minus(i,j,k)/p_minus(i,j,k),1.0)
      else
      r_minus(i,j,k) = 0.0
      endif 
      
      enddo
      enddo
      enddo

       do k=ks,ke+1
       do j=js,je+1
       do i=is,ie+1
       
       if (sx(i,j,k) .ge. 0.0) then
       cx(i,j,k) = min(r_plus(i,j,k),r_minus(i-1,j,k))
       else
       cx(i,j,k) = min(r_plus(i-1,j,k),r_minus(i,j,k))
       endif
       
       if (sy(i,j,k) .ge. 0.0) then
       cy(i,j,k) = min(r_plus(i,j,k),r_minus(i,j-1,k))
       else
       cy(i,j,k) = min(r_plus(i,j-1,k),r_minus(i,j,k))
       endif

       if (sz(i,j,k) .ge. 0.0) then
       cz(i,j,k) = min(r_plus(i,j,k),r_minus(i,j,k-1))
       else
       cz(i,j,k) = min(r_plus(i,j,k-1),r_minus(i,j,k))
       endif

      enddo
      enddo
      enddo
   return
   end subroutine limit


    real(selected_real_kind(8)) function mclimit(theta)
    implicit none
    real(selected_real_kind(8)) theta   
    mclimit = max(0.0,min((1.0 + theta)/2.0,2.0,2.0*theta))
    end function mclimit 

    real(selected_real_kind(8)) function transvel(a,b)
    implicit none
    real(selected_real_kind(8)) a,b,prod 
    prod = a*b
    if (prod .gt. 0.0) then
    transvel = sign(min(abs(a),abs(b)),a)
    else
    transvel = 0.0
    endif
    end function transvel 


    real(selected_real_kind(8)) function minmod4(a,b,c,d)
    implicit none
    real(selected_real_kind(8)) a,b,c,d,one
    one = 1.0
    minmod4 = (1.0/2.0)*(sign(one,a)+sign(one,b))*(1.0/2.0)*(sign(one,a)+sign(one,c))*(1.0/2.0)*(sign(one,a)+sign(one,d)) &
    * min(abs(a),abs(b),abs(c),abs(d))   
    end function minmod4

    real(selected_real_kind(8)) function minmod2(a,b)
    implicit none
    real(selected_real_kind(8)) a,b,one
    one = 1.0
    minmod2 = (1.0/2.0)*(sign(one,a)+sign(one,b))*min(abs(a),abs(b))   
    end function minmod2

    real(selected_real_kind(8)) function median(a,b,c)
    implicit none
    real(selected_real_kind(8)) a,b,c,minmod2
    
    median = a + minmod2(b-a,c-a)    
    end function median




