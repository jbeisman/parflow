/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/

/*****************************************************************************
* C to Fortran interfacing macros
*
*****************************************************************************/

/*advect_new.f90 */
#if defined(_CRAYMPP)
#define ADVECT_UPWIND ADVECT_UPWIND
#elif defined(__bg__)
#define ADVECT_UPWIND advect_upwind
#else
#define ADVECT_UPWIND advect_upwind_
#endif

#define CALL_ADVECT_UPWIND(s, sn, uedge, vedge, wedge, phi,\
                    dlo, dhi, hx, dt,\
                    old_sat, sat, iteration, num_iterations,fx,fy,fz,smin,smax)\
             ADVECT_UPWIND(s, sn, uedge, vedge, wedge, phi,\
                    dlo, dhi, hx, &dt,\
                    old_sat, sat, &iteration, &num_iterations,fx,fy,fz,smin,smax)

void ADVECT_UPWIND(double *s, double *sn,
            double *uedge, double *vedge, double *wedge, double *phi,
            int *dlo, int *dhi, double *hx, double *dt,
            double *old_sat, double *sat, int *iteration, int *num_iterations,
            double *fx, double *fy,double *fz, double *smin, double *smax);


#if defined(_CRAYMPP)
#define ADVECT_HIGHORDER ADVECT_HIGHORDER
#elif defined(__bg__)
#define ADVECT_HIGHORDER advect_highorder
#else
#define ADVECT_HIGHORDER advect_highorder_
#endif

#define CALL_ADVECT_HIGHORDER(s, uedge, vedge, wedge,\
                    dlo, dhi, hx, dt,sx,sy,sz)\
             ADVECT_HIGHORDER(s, uedge, vedge, wedge,\
                    dlo, dhi, hx, &dt,sx,sy,sz)

void ADVECT_HIGHORDER(double *s, double *uedge, double *vedge, double *wedge, 
            int *dlo, int *dhi, double *hx, double *dt,
            double *sx, double *sy,double *sz);



#if defined(_CRAYMPP)
#define ADVECT_TRANSVERSE ADVECT_TRANSVERSE
#elif defined(__bg__)
#define ADVECT_TRANSVERSE advect_transverse
#else
#define ADVECT_TRANSVERSE advect_transverse_
#endif

#define CALL_ADVECT_TRANSVERSE(s, uedge, vedge, wedge,\
                    dlo, dhi, hx, dt,vx,wx,uy,wy,uz,vz,sx,sy,sz)\
             ADVECT_TRANSVERSE(s, uedge, vedge, wedge,\
                    dlo, dhi, hx, &dt,vx,wx,uy,wy,uz,vz,sx,sy,sz)

void ADVECT_TRANSVERSE(double *s, double *uedge, double *vedge, double *wedge, 
            int *dlo, int *dhi, double *hx, double *dt,double *vx,double *wx,
            double *uy,double *wy, 
            double *uz,double *vz,
            double *sx, double *sy,double *sz);


#if defined(_CRAYMPP)
#define ADVECT_COMPUTECONCEN ADVECT_COMPUTECONCEN
#elif defined(__bg__)
#define ADVECT_COMPUTECONCEN advect_computeconcen
#else
#define ADVECT_COMPUTECONCEN advect_computeconcen_
#endif

#define CALL_ADVECT_COMPUTECONCEN(sn, phi,\
                    dlo, dhi, hx, dt,old_sat,sat,\
                    iteration,num_iterations,smin,smax,sx,sy,sz)\
             ADVECT_COMPUTECONCEN(sn, phi,\
                    dlo, dhi, hx, &dt,old_sat,sat,\
                    &iteration,&num_iterations,smin,smax,sx,sy,sz)

void ADVECT_COMPUTECONCEN(double *sn, double *phi, 
            int *dlo, int *dhi, double *hx, double *dt,
            double *old_sat, double *sat,int *iteration, int *num_iterations,
             double *smin, double *smax,
            double *sx, double *sy,double *sz);



#if defined(_CRAYMPP)
#define ADVECT_LIMIT ADVECT_LIMIT
#elif defined(__bg__)
#define ADVECT_LIMIT advect_limit
#else
#define ADVECT_LIMIT advect_limit_
#endif

#define CALL_ADVECT_LIMIT(sn,sx,sy,sz,dlo,dhi,hx,dt,\
                       p_plus,p_minus,q_plus,q_minus,r_plus,r_minus)\
             ADVECT_LIMIT(sn,sx,sy,sz,dlo,dhi,hx,&dt,\
                       p_plus,p_minus,q_plus,q_minus,r_plus,r_minus)

void ADVECT_LIMIT(double *sn, double *sx, double *sy, double *sz, 
            int *dlo, int *dhi, double *hx, double *dt,
             double *p_plus,double *p_minus,
            double *q_plus,double *q_minus, 
            double *r_plus,double *r_minus);

/* sadvect.f */
#if defined(_CRAYMPP)
#define SADVECT SADVECT
#else define (__bg__)
#define SADVECT sadvect
#else
#define SADVECT sadvect_
#endif

#define CALL_SADVECT(s, sn, uedge, vedge, wedge, betaedge, phi, \
                     viscosity, density, gravity, \
                     slx, sly, slz, \
                     lohi, dlohi, hx, dt, \
                     sbot, stop, sbotp, sfrt, sbck, sleft, sright, sfluxz, \
                     dxscr, dyscr, dzscr, dzfrm) \
  SADVECT(s, sn, uedge, vedge, wedge, betaedge, phi, \
          viscosity, density, &gravity, \
          slx, sly, slz, \
          lohi, dlohi, hx, &dt, \
          sbot, stop, sbotp, sfrt, sbck, sleft, sright, sfluxz, \
          dxscr, dyscr, dzscr, dzfrm)

void SADVECT(double *s, double *sn,
             double *uedge, double *vedge, double *wedge, double *betaedge, double *phi,
             double *viscosity, double *density, double *gravity,
             double *slx, double *sly, double *slz,
             int *lohi, int *dlohi, double *hx, double *dt,
             double *sbot, double *stop, double *sbotp,
             double *sfrt, double *sbck,
             double *sleft, double *sright, double *sfluxz,
             double *dxscr, double *dyscr, double *dzscr, double *dzfrm);

/* sk: ftest.f90*/
#define FTEST ftest_

#define CALL_FTEST(outflow_log) \
  FTEST(outflow_log);

void FTEST(double *outflow_log);
