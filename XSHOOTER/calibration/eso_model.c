/* Mathematical kernel from ESO X-shooter pipeline 3.8.12.
 * Copyright (C) 2006 European Southern Observatory.
 * GPL-2.0-or-later; see ESO_COPYING.
 * Source: xsh/xsh_model_kernel.c, xsh_model_kernel.h, xsh_model_ref_ind.h.
 * https://ftp.eso.org/pub/dfs/pipelines/instruments/xshooter/xshoo-kit-3.8.12-3.tar.gz
 * Mathematical function bodies are unchanged. CPL I/O and unused functions
 * are omitted; the last four functions expose the kernel for offline table
 * preparation only. This file is not part of the ScopeSim runtime.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NIR_FIX 36
#define mm2um 1000.0
#define mm2nm 1000000.0
#define XSH_TRACE_MSG(x)
static const int vectordim=4;
typedef double DOUBLE;


typedef DOUBLE vec[4];
typedef DOUBLE cvec[4];
typedef DOUBLE mat[4][4];

typedef int detloc[3];
typedef int XY[3];

typedef struct {
  DOUBLE best,min,max;
  int flag;
  char  name[20];
} ann_all_par;

typedef struct
  {
    int counter;
    double x;
    double y;
    int arm;
    double wave;
    double flux;
    int slit_pos;
    int order;
  } coord ;

struct xs_3
{

  int arm;
  double chipxpix,chipypix;
  double xsize_corr, ysize_corr;
  int ASIZE, BSIZE, SIZE;
  int morder, morder_min, morder_max;
  DOUBLE blaze_pad;
  DOUBLE temper, t_ir_p2, t_ir_p3;
  DOUBLE mues, nues, taues, slit_scale, es_s, es_w, es_x, es_xf, es_y, es_y_tot;
  DOUBLE fcol;
  DOUBLE mup1, cmup1, nup1, taup1;
  DOUBLE rind, rind2, rind3, rind_a1, rind_a2, rind_a3, rind_b1, rind_b2, rind_lm1, rind_lm2, rind_lm3, rind_wta0, rind_wta1, rind_wta2, rind_wtb0, rind_wtb1, rind_wtb2;
  DOUBLE rind_delta;
  DOUBLE mup2, nup2, taup2;
  DOUBLE mup3, nup3, taup3;
  DOUBLE mup4, nup4, taup4;
  DOUBLE mup5, nup5, taup5;
  DOUBLE mup6, nup6, taup6;

  DOUBLE mug, nug, taug, sg, grat_alpha, grat_beta;
  DOUBLE mud, nud, taud, pix, pix_X, pix_Y, offx, offy, flipx, flipy, xpospix, ypospix, fdet;

  detloc chippix;
  double xdet,ydet;

  double chipx, chipy, chipxup, chipyup, chipxdown, chipydown, chiprot;


  double pc_x_xx, pc_x_x1, pc_x_yy, pc_x_y1, pc_x_xy, pc_x_x3, pc_x_x2y, pc_x_y2x, pc_x_y3;
  double pc_y_xx, pc_y_x1, pc_y_yy, pc_y_y1, pc_y_xy, pc_y_x3, pc_y_x2y, pc_y_y2x, pc_y_y3;

  double pc4_x_xy3,pc4_x_x3y,pc4_x_x2y2,pc4_x_x4,pc4_x_y4;
  double pc4_y_xy3,pc4_y_x3y,pc4_y_x2y2,pc4_y_x4,pc4_y_y4;


  double d2_x3, d2_x2, d2_x1;
  double d2_y3x3, d2_y3x2, d2_y3x1,d2_y3x0,d2_y2x3, d2_y2x2, d2_y2x1,d2_y2x0,d2_y1x3, d2_y1x2, d2_y1x1,d2_y1x0;

  double ca_x0,ca_x1,ca_y0,ca_y1;


  mat e_slit;
  mat mup_ir_cor, mup_ir_cor_out, toprism1, toprism2, toprism3, toprism4, toprism5, toprism6, prism_out, ret_prism2, ret_prism1, ret_prism_out1,ret_prism4, ret_prism3, ret_prism_out2,ret_prism6, ret_prism5, ret_prism_out3,todetector;
  mat tograt, grat_out;
  DOUBLE slit[10];


  double config_mjd;
};
typedef struct xs_3 xsh_xs_3;
/* $Id: xsh_model_ref_ind.h,v 1.4 2010-07-28 11:17:30 amodigli Exp $
 *
 *   This file is part of the ESO X-shooter Pipeline
 *   Copyright (C) 2006 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */



#ifndef XSH_MODEL_REF_IND_H
#define XSH_MODEL_REF_IND_H

static double schott_sf6_manu [6][7] = {
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,260.0},
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,270.0},
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,280.0},
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,290.0},
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,300.0},
  {1.72618475,0.392383376,1.05068036,0.0134871947,0.0569318095,118.557185,310.0}
  };

static double silica_zemax [6][7] = {
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,260.0},
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,270.0},
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,280.0},
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,290.0},
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,300.0},
  {0.6961663,0.4079426,0.8974794,0.004679148,0.01351206,97.934,310.0}};

static double znse [6][7] = {
  {4.43472613286,0.437534439541,2.94296707239,0.20047046355,0.383301991611,47.9937969433,95.0},
  {4.43542230714,0.438059580351,2.9265319319,0.200506374102,0.383453912766,47.8502966353,100.0},
  {4.43612895652,0.438586230317,2.91141326949,0.200541768113,0.383612197906,47.7172014285,105.0},
  {4.43604895636,0.439902006841,2.88537551027,0.200541339195,0.38368776333,47.4933353795,110.0},
  {4.43675967885,0.440440921835,2.87838606989,0.200575497514,0.383853032135,47.4256517629,115.0},
  {4.4374705012,0.440979834519,2.87139508026,0.200609780065,0.384017392306,47.3579544185,120.0}};
#endif



void xsh_3_detpix(struct xs_3 *p_xs_3)
;
void xsh_3_eval(DOUBLE  lambda,
	       int morder,
	       DOUBLE** ref_ind,
	       struct xs_3* p_xs_3)
;
void xsh_3_init(struct xs_3 * p_xs_3)
;
static void xsh_copyvector(vec a, vec b)
;
void xsh_matrixforvector(vec a,mat B,vec c)
;
double xsh_model_sellmeier_ext(int arm, double temper, double lam_sqr)
;
void xsh_multiplymatrix(mat A, mat B, mat C)
;
void xsh_normall(vec a)
;
void xsh_normz(vec a)
;
static void xsh_nullmatrix(mat A)
;
static void xsh_nullvector(vec a)
;
void xsh_ref_ind_read(int arm, DOUBLE** ref_ind_T,DOUBLE temper)
;
void xsh_refract(vec b, DOUBLE rind, vec a)
;
void
xsh_rotationmatrix(mat A, const char axis, const DOUBLE angle)
;
void
xsh_rotin(mat A,const DOUBLE x_angle, const DOUBLE y_angle,const DOUBLE z_angle)
;
void xsh_subtractvectors(vec a, vec b)
;
void xsh_transpose(mat A ,mat B)
;
void xsh_3_detpix(struct xs_3 *p_xs_3)
{

  double xdet2, ydet2;

  double xdet_c, xdet_c2, xdet_c3, ydet_c, ydet_c2, ydet_c3;
  double d2_y1=0.0;
  double d2_y2=0.0;
  double d2_y3=0.0;
  double pc_x, pc_y, dstn_x, dstn_y;
  double xdet_c4, ydet_c4;

  int array_size_fix=0;




  p_xs_3->xpospix=0.0;
  p_xs_3->ypospix=0.0;
  p_xs_3->chippix[0]=-1;
  p_xs_3->chippix[1]=0;
  p_xs_3->chippix[2]=0;

#ifndef FCCDD_FLAG



  ydet2=p_xs_3->ydet*p_xs_3->ydet;
  double ydet3=ydet2*p_xs_3->ydet;
  xdet2=p_xs_3->xdet*p_xs_3->xdet;
  double xdet3=xdet2*p_xs_3->xdet;
  xdet_c=p_xs_3->xdet-p_xs_3->chipx;

  xdet_c2=xdet_c*xdet_c;
  xdet_c3=xdet_c2*xdet_c;
  xdet_c4=xdet_c3*xdet_c;
  ydet_c=p_xs_3->ydet-p_xs_3->chipy;

  ydet_c2=ydet_c*ydet_c;
  ydet_c3=ydet_c2*ydet_c;
  ydet_c4=ydet_c3*ydet_c;


  pc_x=p_xs_3->chipx+
    (p_xs_3->pc_x_xx*xdet_c2)+
    (p_xs_3->pc_x_x1*xdet_c)+
    (p_xs_3->pc_x_yy*ydet_c2)+
    (p_xs_3->pc_x_y1*ydet_c)+
    (p_xs_3->pc_x_xy*xdet_c*ydet_c)+
    (p_xs_3->pc_x_x3*xdet_c3)+
    (p_xs_3->pc_x_x2y*xdet_c2*ydet_c)+
    (p_xs_3->pc_x_y2x*xdet_c*ydet_c2)+
    (p_xs_3->pc_x_y3*ydet_c3)+
    (p_xs_3->pc4_x_xy3*xdet_c*ydet_c3)+
    (p_xs_3->pc4_x_x3y*xdet_c3*ydet_c)+
    (p_xs_3->pc4_x_x2y2*xdet_c2*ydet_c2)+
    (p_xs_3->pc4_x_x4*xdet_c4)+
    (p_xs_3->pc4_x_y4*ydet_c4);
  pc_y=p_xs_3->chipy+
    (p_xs_3->pc_y_xx*xdet_c2)+
    (p_xs_3->pc_y_x1*xdet_c)+
    (p_xs_3->pc_y_yy*ydet_c2)+
    (p_xs_3->pc_y_y1*ydet_c)+
    (p_xs_3->pc_y_xy*xdet_c*ydet_c)+
    (p_xs_3->pc_y_x3*xdet_c3)+
    (p_xs_3->pc_y_x2y*xdet_c2*ydet_c)+
    (p_xs_3->pc_y_y2x*xdet_c*ydet_c2)+
    (p_xs_3->pc_y_y3*ydet_c3)+
    (p_xs_3->pc4_y_xy3*xdet_c*ydet_c3)+
    (p_xs_3->pc4_y_x3y*xdet_c3*ydet_c)+
    (p_xs_3->pc4_y_x2y2*xdet_c2*ydet_c2)+
    (p_xs_3->pc4_y_x4*xdet_c4)+
    (p_xs_3->pc4_y_y4*ydet_c4);






















  p_xs_3->xdet=pc_x;
  p_xs_3->ydet=pc_y;


  if (p_xs_3->arm<2) {
    dstn_y=(p_xs_3->d2_x3*ydet_c3)+(p_xs_3->d2_x2*ydet_c2)+(p_xs_3->d2_x1*ydet_c);

    p_xs_3->ydet=p_xs_3->ydet+dstn_y;







  }
  else if (p_xs_3->arm==2) {

    if (p_xs_3->d2_x3!=0.0) {
      dstn_x=(p_xs_3->d2_x3*xdet_c3)+(p_xs_3->d2_x2*xdet_c2)+(p_xs_3->d2_x1*xdet_c);

      p_xs_3->xdet=p_xs_3->xdet+dstn_x;
      d2_y3=(p_xs_3->d2_y3x3*xdet_c3)+(p_xs_3->d2_y3x2*xdet_c2)+(p_xs_3->d2_y3x1*xdet_c)+p_xs_3->d2_y3x0;
      d2_y2=(p_xs_3->d2_y2x3*xdet_c3)+(p_xs_3->d2_y2x2*xdet_c2)+(p_xs_3->d2_y2x1*xdet_c)+p_xs_3->d2_y2x0;
      d2_y1=(p_xs_3->d2_y1x3*xdet_c3)+(p_xs_3->d2_y1x2*xdet_c2)+(p_xs_3->d2_y1x1*xdet_c)+p_xs_3->d2_y1x0;




      dstn_y=(d2_y3*ydet_c3)+(d2_y2*ydet_c2)+(d2_y1*ydet_c);
      p_xs_3->ydet=p_xs_3->ydet+dstn_y;
    }

    else {
      dstn_y=(p_xs_3->ca_x1*(xdet_c-p_xs_3->ca_x0)+p_xs_3->ca_y1*(ydet_c-p_xs_3->ca_y0))*p_xs_3->es_y_tot;
      p_xs_3->ydet=p_xs_3->ydet+dstn_y;
    }
  }
#endif
#ifdef CCCDD_FLAG

  h_det=sqrt(p_xs_3->xdet*p_xs_3->xdet+p_xs_3->xdet*p_xs_3->xdet);

  p_xs_3->xdet=(p_xs_3->fdet*p_xs_3->xdet/h_det)*sin(h_det/p_xs_3->fdet);
  p_xs_3->ydet=(p_xs_3->fdet*p_xs_3->ydet/h_det)*sin(h_det/p_xs_3->fdet);
#endif




















    p_xs_3->chippix[0]=1;


    if(p_xs_3->arm!=2) {
      p_xs_3->xpospix=0.5*(p_xs_3->chipxpix+p_xs_3->xsize_corr)-(p_xs_3->ydet-p_xs_3->chipy)/p_xs_3->pix_X;
      p_xs_3->ypospix=0.5*(p_xs_3->chipypix+p_xs_3->ysize_corr)-(p_xs_3->xdet-p_xs_3->chipx)/p_xs_3->pix_Y;
    }
    else {
      p_xs_3->xpospix=0.5*(p_xs_3->chipxpix+p_xs_3->xsize_corr)+(p_xs_3->xdet-p_xs_3->chipx)/p_xs_3->pix_X;
      p_xs_3->ypospix=0.5*(p_xs_3->chipypix+p_xs_3->ysize_corr)+(p_xs_3->ydet-p_xs_3->chipy)/p_xs_3->pix_Y;
      array_size_fix=NIR_FIX;
    }


    if (p_xs_3->pix_X>0.029) {
      p_xs_3->xpospix=p_xs_3->xpospix+0.25;
    }
    if (p_xs_3->pix_Y>0.029) {
      p_xs_3->ypospix=p_xs_3->ypospix+0.25;
    }

    if (p_xs_3->xpospix>p_xs_3->chipxpix-1+array_size_fix) {
      p_xs_3->xpospix=-1.0;
      p_xs_3->ypospix=-1.0;
      p_xs_3->chippix[0]=6;
    }
    if (p_xs_3->xpospix<0) {
      p_xs_3->xpospix=-1.0;
      p_xs_3->ypospix=-1.0;
      p_xs_3->chippix[0]=5;
    }
    if (p_xs_3->ypospix<0) {
      p_xs_3->xpospix=-1.0;
      p_xs_3->ypospix=-1.0;
      p_xs_3->chippix[0]=8;
    }
    if (p_xs_3->ypospix>p_xs_3->chipypix-1) {
      p_xs_3->xpospix=-1.0;
      p_xs_3->ypospix=-1.0;
      p_xs_3->chippix[0]=7;
    }

      p_xs_3->chippix[1]=(int)(p_xs_3->xpospix-0.5)+1;
      p_xs_3->chippix[2]=(int)(p_xs_3->ypospix-0.5)+1;























  return;
}
void xsh_3_eval(DOUBLE  lambda,
	       int morder,
	       DOUBLE** ref_ind,
	       struct xs_3* p_xs_3)
{
  DOUBLE lambda_um,lamf,rind_up,rind_down,xx,yy,dndT;
  vec vres,slitpos2,slitpos,detproj;

  mat mech={ {1,0,0,0},
	     {morder*(p_xs_3->sg),1,0,0},
	     {0,0,1,0},
	     {0,0,0,-1}};
  mat chrom_corr={ {1,0,0,0},
		   {0,1,0,0},
		   {0,0,1,0},
		   {0,0,0,1}};



  vec vres1={lambda,
	      0,
	      0,
	      -p_xs_3->fcol};


  XSH_TRACE_MSG(vres1)


  slitpos[0]=0.0;
  slitpos[1]=1.0*p_xs_3->es_x;

  if (p_xs_3->config_mjd>2455409.0) {

    if (p_xs_3->arm==0) {
      if (lambda>0.000425) {
	slitpos[2]=-1.0*p_xs_3->es_y_tot*(1.0+(lambda-0.000425)*p_xs_3->offx);
      }
      else{
	slitpos[2]=-1.0*p_xs_3->es_y_tot*(1.0+(lambda-0.000425)*p_xs_3->offy);
      }
    }
    else if (p_xs_3->arm==1) {
      if (lambda>0.000650) {
	slitpos[2]=1.0*(p_xs_3->es_y_tot)*(1.0+(lambda-0.000900)*p_xs_3->offx);
      }
      else{
	slitpos[2]=1.0*(p_xs_3->es_y_tot)*(1.0+(lambda-0.000018)*p_xs_3->offy);
      }
    }
    else if (p_xs_3->arm==2) {
      if (lambda>0.001400) {
	slitpos[2]=1.0*(p_xs_3->es_y_tot-p_xs_3->es_y)*(1.0+(lambda-0.001600)*p_xs_3->offx)+p_xs_3->es_y;
      }
      else{
	slitpos[2]=1.0*(p_xs_3->es_y_tot)*(1.0+(lambda-0.000990)*p_xs_3->offy);
      }
    }
  }
  else {
    if (p_xs_3->arm==0) {
      slitpos[2]=-1.0*p_xs_3->es_y_tot;
    }
    else {
      slitpos[2]=1.0*p_xs_3->es_y_tot;
    }
  }
  slitpos[3]=0.0;
  xsh_matrixforvector(slitpos2,p_xs_3->e_slit,slitpos);
  xsh_subtractvectors(vres1,slitpos2);
  xsh_copyvector(vres,vres1);
  xsh_normall(vres);
  if (p_xs_3->arm!=2 && p_xs_3->config_mjd<2455409.0) {

    xsh_rotin(chrom_corr,vres[2]*(lambda-p_xs_3->offy)*p_xs_3->offx,0.0,0.0);
    xsh_matrixforvector(vres,chrom_corr,vres);
  }

  XSH_TRACE_MSG(vres)



  if (p_xs_3->arm==2) xsh_matrixforvector(vres,p_xs_3->mup_ir_cor,vres);

  xsh_matrixforvector(vres,p_xs_3->toprism1,vres);


  lambda_um=lambda* mm2um;
  lamf=lambda_um*lambda_um;

  if (p_xs_3->arm==1) {


    p_xs_3->rind=sqrt(ref_ind[0][0]*lamf/(lamf-ref_ind[0][3])+ref_ind[0][1]*lamf/(lamf-ref_ind[0][4])+ref_ind[0][2]*lamf/(lamf-ref_ind[0][5])+1.0);


    double lam_ig2=0.04963984;
    double H=71.0867E-06;
    double G=-50.2451E-06;
    double R;

    R=lamf/(lamf-lam_ig2);
    dndT=(G*R+H*R*R)/(2.0*p_xs_3->rind);




    p_xs_3->rind+=dndT*(p_xs_3->temper-293.15);

  }
  else {

    p_xs_3->rind=xsh_model_sellmeier_ext(p_xs_3->arm,p_xs_3->temper,lamf);
  }







  xsh_refract(vres,p_xs_3->rind,vres);



  xsh_matrixforvector(vres,p_xs_3->toprism2,vres);



  xsh_refract(vres,1.0/(p_xs_3->rind),vres);




  xsh_matrixforvector(vres,p_xs_3->toprism1,vres);




 XSH_TRACE_MSG(vres)
  if (p_xs_3->arm==2) {
 XSH_TRACE_MSG(vres)


    xsh_matrixforvector(vres,p_xs_3->toprism3,vres);




    rind_up=sqrt(ref_ind[1][0]*lamf/(lamf-ref_ind[1][3]*ref_ind[1][3])+ref_ind[1][1]*lamf/(lamf-ref_ind[1][4]*ref_ind[1][4])+ref_ind[1][2]*lamf/(lamf-ref_ind[1][5]*ref_ind[1][5])+1.0);


    rind_down=sqrt(ref_ind[0][0]*lamf/(lamf-ref_ind[0][3]*ref_ind[0][3])+ref_ind[0][1]*lamf/(lamf-ref_ind[0][4]*ref_ind[0][4])+ref_ind[0][2]*lamf/(lamf-ref_ind[0][5]*ref_ind[0][5])+1.0);


    p_xs_3->rind2=rind_down+((p_xs_3->t_ir_p2-ref_ind[0][6])/(ref_ind[1][6]-ref_ind[0][6]))*(rind_up-rind_down);


    xsh_refract(vres,p_xs_3->rind2,vres);



    xsh_matrixforvector(vres,p_xs_3->toprism4,vres);



    xsh_refract(vres,1.0/(p_xs_3->rind2),vres);




    xsh_matrixforvector(vres,p_xs_3->toprism3,vres);


 XSH_TRACE_MSG(vres)


    xsh_matrixforvector(vres,p_xs_3->toprism5,vres);
 XSH_TRACE_MSG(vres)



    rind_up=sqrt(ref_ind[1][0]*lamf/(lamf-ref_ind[1][3]*ref_ind[1][3])+ref_ind[1][1]*lamf/(lamf-ref_ind[1][4]*ref_ind[1][4])+ref_ind[1][2]*lamf/(lamf-ref_ind[1][5]*ref_ind[1][5])+1.0);


    rind_down=sqrt(ref_ind[0][0]*lamf/(lamf-ref_ind[0][3]*ref_ind[0][3])+ref_ind[0][1]*lamf/(lamf-ref_ind[0][4]*ref_ind[0][4])+ref_ind[0][2]*lamf/(lamf-ref_ind[0][5]*ref_ind[0][5])+1.0);


    p_xs_3->rind3=rind_down+((p_xs_3->t_ir_p3-ref_ind[0][6])/(ref_ind[1][6]-ref_ind[0][6]))*(rind_up-rind_down);


    xsh_refract(vres,p_xs_3->rind3,vres);
 XSH_TRACE_MSG(vres)




    xsh_matrixforvector(vres,p_xs_3->toprism6,vres);
 XSH_TRACE_MSG(vres)




    xsh_refract(vres,1.0/(p_xs_3->rind3),vres);
 XSH_TRACE_MSG(vres)






    xsh_matrixforvector(vres,p_xs_3->toprism5,vres);
 XSH_TRACE_MSG(vres)



  }









  xsh_matrixforvector(vres,p_xs_3->tograt,vres);
 XSH_TRACE_MSG(vres)


  p_xs_3->grat_alpha=atan(vres[1]/vres[3]);



  xsh_matrixforvector(vres,mech,vres);
 XSH_TRACE_MSG(vres)




  xsh_normz(vres);
 XSH_TRACE_MSG(vres)


  p_xs_3->grat_beta=atan(vres[1]/vres[3]);



  xsh_matrixforvector(vres,p_xs_3->grat_out,vres);
 XSH_TRACE_MSG(vres)




  if (p_xs_3->arm==2) {

    xsh_matrixforvector(vres,p_xs_3->ret_prism6,vres);
 XSH_TRACE_MSG(vres)





    xsh_refract(vres,p_xs_3->rind3,vres);
 XSH_TRACE_MSG(vres)





    xsh_matrixforvector(vres,p_xs_3->ret_prism5,vres);
 XSH_TRACE_MSG(vres)





    xsh_refract(vres,1.0/(p_xs_3->rind3),vres);
 XSH_TRACE_MSG(vres)




    xsh_matrixforvector(vres,p_xs_3->ret_prism6,vres);
 XSH_TRACE_MSG(vres)





    xsh_matrixforvector(vres,p_xs_3->ret_prism4,vres);
 XSH_TRACE_MSG(vres)




    xsh_refract(vres,p_xs_3->rind2,vres);
 XSH_TRACE_MSG(vres)




    xsh_matrixforvector(vres,p_xs_3->ret_prism3,vres);
 XSH_TRACE_MSG(vres)




    xsh_refract(vres,1.0/(p_xs_3->rind2),vres);
 XSH_TRACE_MSG(vres)




    xsh_matrixforvector(vres,p_xs_3->ret_prism4,vres);
 XSH_TRACE_MSG(vres)


  }


  xsh_matrixforvector(vres,p_xs_3->ret_prism2,vres);
 XSH_TRACE_MSG(vres)




  xsh_refract(vres,p_xs_3->rind,vres);
 XSH_TRACE_MSG(vres)




  xsh_matrixforvector(vres,p_xs_3->ret_prism1,vres);
 XSH_TRACE_MSG(vres)




  xsh_refract(vres,1.0/(p_xs_3->rind),vres);
 XSH_TRACE_MSG(vres)




  xsh_matrixforvector(vres,p_xs_3->ret_prism2,vres);
 XSH_TRACE_MSG(vres)


  if (p_xs_3->arm==2) xsh_matrixforvector(vres,p_xs_3->mup_ir_cor_out,vres);
 XSH_TRACE_MSG(vres)


#ifdef FCCDD_FLAG
  xsh_matrixforvector(vres,p_xs_3->todetector,vres);
 XSH_TRACE_MSG(vres)

#endif


  xx=(p_xs_3->flipx)*(p_xs_3->fdet)/((vres[3]/vres[1])*cos(p_xs_3->nud)-sin(p_xs_3->nud));
  yy=(p_xs_3->flipy)*(p_xs_3->fdet)/((vres[3]/vres[2])*cos(p_xs_3->mud)-sin(p_xs_3->mud));




#ifdef CCCDD_FLAG

    xx=(p_xs_3->flipx)*(p_xs_3->fdet)*vres[1];
    yy=(p_xs_3->flipy)*(p_xs_3->fdet)*vres[2];
#endif

#ifdef XCCCDD_FLAG

    xx=(p_xs_3->flipx)*(p_xs_3->fdet)*2.0*vres[1];
    yy=(p_xs_3->flipy)*(p_xs_3->fdet)*2.0*vres[2];
#endif


  detproj[0]=0.0;
  detproj[1]=xx*cos(-p_xs_3->taud)+yy*sin(-p_xs_3->taud);
  detproj[2]=-xx*sin(-p_xs_3->taud)+yy*cos(-p_xs_3->taud);
  detproj[3]=0.0;
  xsh_copyvector(vres,detproj);
 XSH_TRACE_MSG(vres)


  p_xs_3->morder=(float)(morder);

  p_xs_3->xdet=vres[1];
  p_xs_3->ydet=vres[2];
}
void xsh_3_init(struct xs_3 * p_xs_3)
{

#ifdef FCCDD_FLAG
  xsh_rotin(p_xs_3->todetector,-0.062,0.01,0.0);
#endif



 xsh_rotin(p_xs_3->e_slit,p_xs_3->mues,p_xs_3->nues,p_xs_3->taues);


 xsh_rotin(p_xs_3->mup_ir_cor,p_xs_3->cmup1,0.0,0.0);
 xsh_rotin(p_xs_3->mup_ir_cor_out,-p_xs_3->cmup1,0.0,0.0);
 xsh_rotin(p_xs_3->toprism1,p_xs_3->mup1,p_xs_3->nup1,p_xs_3->taup1);
 xsh_rotin(p_xs_3->toprism2,p_xs_3->mup2,p_xs_3->nup2,p_xs_3->taup2);







 xsh_rotin(p_xs_3->toprism3,p_xs_3->mup3,p_xs_3->nup3,p_xs_3->taup3);
 xsh_rotin(p_xs_3->toprism4,p_xs_3->mup4,p_xs_3->nup4,p_xs_3->taup4);






 xsh_rotin(p_xs_3->toprism5,p_xs_3->mup5,p_xs_3->nup5,p_xs_3->taup5);
 xsh_rotin(p_xs_3->toprism6,p_xs_3->mup6,p_xs_3->nup6,p_xs_3->taup6);






 xsh_rotin(p_xs_3->tograt,p_xs_3->mug,p_xs_3->nug,p_xs_3->taug);
 xsh_transpose(p_xs_3->grat_out,p_xs_3->tograt);





 xsh_transpose(p_xs_3->ret_prism2,p_xs_3->toprism1);
 xsh_transpose(p_xs_3->ret_prism1,p_xs_3->toprism2);
 xsh_multiplymatrix(p_xs_3->ret_prism_out1,p_xs_3->toprism2,p_xs_3->toprism1);
 xsh_transpose(p_xs_3->ret_prism4,p_xs_3->toprism3);
 xsh_transpose(p_xs_3->ret_prism3,p_xs_3->toprism4);
 xsh_multiplymatrix(p_xs_3->ret_prism_out2,p_xs_3->toprism4,p_xs_3->toprism3);
 xsh_transpose(p_xs_3->ret_prism6,p_xs_3->toprism5);
 xsh_transpose(p_xs_3->ret_prism5,p_xs_3->toprism6);
 xsh_multiplymatrix(p_xs_3->ret_prism_out3,p_xs_3->toprism6,p_xs_3->toprism5);


 if (p_xs_3->arm!=2) {
   p_xs_3->chipxup=p_xs_3->chipx+p_xs_3->pix_Y*(0.5*p_xs_3->chipypix);
   p_xs_3->chipxdown=p_xs_3->chipx-p_xs_3->pix_Y*(0.5*p_xs_3->chipypix);
   p_xs_3->chipyup=p_xs_3->chipy+p_xs_3->pix_X*(0.5*p_xs_3->chipxpix);
   p_xs_3->chipydown=p_xs_3->chipy-p_xs_3->pix_X*(0.5*p_xs_3->chipxpix);
 }
 else {
   p_xs_3->chipxup=p_xs_3->chipx+p_xs_3->pix_X*(0.5*p_xs_3->chipxpix+NIR_FIX);
   p_xs_3->chipxdown=p_xs_3->chipx-p_xs_3->pix_X*(0.5*p_xs_3->chipxpix);
   p_xs_3->chipyup=p_xs_3->chipy+p_xs_3->pix_Y*(+0.5*p_xs_3->chipypix);
   p_xs_3->chipydown=p_xs_3->chipy-p_xs_3->pix_Y*(0.5*p_xs_3->chipypix);
 }
}
static void xsh_copyvector(vec a, vec b)
{
  int i;
  for(i=0; i<vectordim; i++) a[i]=b[i];

}
void xsh_matrixforvector(vec a,mat B,vec c)
{
  register int i,j;
  vec temp;
  xsh_nullvector(temp);
  for(j=0; j<vectordim; j++) {
    double c_j = c[j];
    for(i=0; i<vectordim; i++) {
        temp[i] +=B[i][j]*c_j;
    }
  }
  xsh_copyvector(a,temp);
}
double xsh_model_sellmeier_ext(int arm, double temper, double lam_sqr)
{
  typedef DOUBLE sellmeier_mat[3][5];
  int ii, jj;
  double ref_ind,top,bottom;
  double T[5];
  sellmeier_mat S_UVB={{1.10127,     -4.94251E-05, 5.27414E-07,  -1.59700E-09, 1.75949E-12},
		       {1.78752E-05, 4.76391E-05,  -4.49019E-07, 1.44546E-09,  -1.57223E-12},
		       {0.793552,    -1.27815E-03, 1.84595E-05,  -9.20275E-08, 1.48829E-10}};
  sellmeier_mat L_UVB={{-0.08906,    9.08730E-06,  -6.53638E-08, 7.77072E-11,  6.84605E-14},
		       {0.297562,    -859578E-04,  6.59069E-06,  -1.09482E-08, 7.85145E-13},
		       {9.34454,     -7.09788E-03, 1.01968E-04,  -5.07660E-07, 8.21348E-10}};

  sellmeier_mat S_NIR={{0.104940144,5.40E-06,3.23E-08,1.83E-13,-3.60E-14},
		   {0.996335571,-3.29E-06,1.48E-08,-3.01E-11,4.99E-14},
		   {0.832484961,6.38E-04,-2.40E-06,6.10E-10,4.77E-12}};
  sellmeier_mat L_NIR={{-3.07E-03,-2.56E-05,6.21E-07,-2.54E-09,2.88E-12},
		   {9.40E-02,-1.59E-06,1.28E-08,1.82E-12,-3.01E-14},
		   {9.598633568,3.15E-03,-1.22E-05,5.48E-09,1.96E-11}};
  T[0]=1.0;
  T[1]=temper;
  T[2]=temper*temper;
  T[3]=T[2]*temper;
  T[4]=T[3]*temper;
  ref_ind=1.0;
  if (arm==0) {
    for (ii=0;ii<3;ii++) {
      top=bottom=0.0;
      for (jj=0;jj<5;jj++) {
	top+=S_UVB[ii][jj]*T[jj];
	bottom+=L_UVB[ii][jj]*T[jj];
      }
      ref_ind+=(top*lam_sqr)/(lam_sqr-(bottom*bottom));
    }
  }
  else if (arm==2) {
    for (ii=0;ii<3;ii++) {
      top=bottom=0.0;
      for (jj=0;jj<5;jj++) {
	top+=S_NIR[ii][jj]*T[jj];
	bottom+=L_NIR[ii][jj]*T[jj];
      }
      ref_ind+=(top*lam_sqr)/(lam_sqr-(bottom*bottom));
    }
  }
  return sqrt(ref_ind);
}
void xsh_multiplymatrix(mat A, mat B, mat C)
{
  register int i,j,k;
  xsh_nullmatrix(A);
  for(i=0; i<vectordim; i++) {
    for(j=0; j<vectordim; j++) {
      for(k=0; k<vectordim; k++){
	  A[i][j] +=B[i][k]*C[k][j];
	}
    }
  }
}
void xsh_normall(vec a)
{
   DOUBLE mag;
   mag=1./sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3]);
   a[1]=a[1]*mag;
   a[2]=a[2]*mag;
   a[3]=a[3]*mag;
}
void xsh_normz(vec a)
{
  double x2_sqr = a[2]*a[2];
  if (x2_sqr > 1.0) {
    a[1]=0.0;
    a[2]=1.0;
    a[3]=0.0;
    return;
  }
  double x1_sqr = a[1]*a[1];
  if (x1_sqr > 1.0) {
    a[1]=1.0;
    a[2]=0.0;
    a[3]=0.0;
    return;
  }
  double sum_sqrs = x1_sqr + x2_sqr;
  if (sum_sqrs > 1.0) {
    double norm = 1./sqrt (sum_sqrs);
    a[1] *=norm;
    a[2] *=norm;
    a[3]=0.0;
  }
  else{
    if (a[3]>=0) a[3]=sqrt(1.0- sum_sqrs);
    else a[3]=-sqrt(1.0-sum_sqrs);
  }

  return;
}
static void xsh_nullmatrix(mat A)
{

  register int i,j;
   for(i=0; i<vectordim; i++) {
    for(j=0; j<vectordim; j++) {
      A[i][j]=0.0;
    }
   }



}
static void xsh_nullvector(vec a)
{

  int i;
  for(i=0; i<vectordim; i++) a[i]=0.0;


}
void xsh_ref_ind_read(int arm, DOUBLE** ref_ind_T,DOUBLE temper)
{
  int ii,jj,flag;



  ii=2;
  flag=0;
  for (ii=2; ii<8; ii+=1) {
    for (jj=0;jj<7;jj+=1){
      if (arm==0) {
	ref_ind_T[ii][jj]=silica_zemax[ii-2][jj];
      }
      else if (arm==1) {
	ref_ind_T[ii][jj]=schott_sf6_manu[ii-2][jj];
      }
      else {
	ref_ind_T[ii][jj]=znse[ii-2][jj];
      }
    }

    if (temper>=ref_ind_T[ii][6])
    {

      flag=1;
      for (jj=0; jj<7; jj++) {
	ref_ind_T[0][jj]=ref_ind_T[ii][jj];
      }
    }
    else if (temper<ref_ind_T[ii][6] && flag==1)
    {

      flag=2;
      for (jj=0; jj<7; jj++) {
	ref_ind_T[1][jj]=ref_ind_T[ii][jj];
      }
    }
  }

  if (flag!=2) {
    printf("******* Temperature out of range! ******* %lf \n", temper);
    for (jj=0; jj<7; jj++) {
      ref_ind_T[1][jj]=ref_ind_T[7][jj];
    }
  }
  return;
}
void xsh_refract(vec b, DOUBLE rind, vec a)
{
  DOUBLE phi, pheta, y_signflag, z_signflag, sin_refrangle;
  y_signflag=1.0;
  z_signflag=1.0;
  phi=atan(sqrt(a[1]*a[1]+a[2]*a[2])/a[3]);
  pheta=atan(a[1]/a[2]);
  if (a[2]<0.0) y_signflag=-1.0;
  if (a[3]<0.0) z_signflag=-1.0;

  sin_refrangle = sin(phi)/rind;
  b[0]=a[0];
  b[1]=z_signflag*y_signflag*sin(pheta)*sin_refrangle;
  b[2]=z_signflag*y_signflag*cos(pheta)*sin_refrangle;
  b[3]=z_signflag*cos(asin(sin_refrangle));
  return;
}
void
xsh_rotationmatrix(mat A, const char axis, const DOUBLE angle)
{
  xsh_nullmatrix(A);
  double sin_a=sin(angle);
  double cos_a=cos(angle);


  switch(axis){

  case 'x':{

  A[0][0]=A[1][1]=1;
  A[2][2]=A[3][3]=cos_a;
  A[2][3]=sin_a;
  A[3][2]= -A[2][3];
  break;
  }


  case 'y':{

  A[0][0]=A[2][2]=1;
  A[1][1]=A[3][3]=cos_a;
  A[1][3]=-sin_a;
  A[3][1]=-A[1][3];
  break;
  }


  case 'z':{
  A[0][0]=A[3][3]=1;
  A[1][1]=A[2][2]=cos_a;
  A[1][2]=sin_a;
  A[2][1]=-A[1][2];
  break;
  }
  printf("Error on creating rotation matrix\n");
 }

}
void
xsh_rotin(mat A,const DOUBLE x_angle, const DOUBLE y_angle,const DOUBLE z_angle)
{
  mat temp;
  mat X,Y,Z;
  xsh_rotationmatrix(X,'x',x_angle);
  xsh_rotationmatrix(Y,'y',y_angle);
  xsh_rotationmatrix(Z,'z',z_angle);

  xsh_multiplymatrix(temp,Y,X);
  xsh_multiplymatrix(A,Z,temp);
}
void xsh_subtractvectors(vec a, vec b)
{
  int i;
  for(i=0; i<vectordim; i++) a[i]-=b[i];
}
void xsh_transpose(mat A ,mat B)
{
   register int i,j;
   for(i=0; i<vectordim; i++) {
    for(j=0; j<vectordim; j++) {
      A[i][j]=B[j][i];
    }
   }
}
static struct xs_3 model;
static double indices[8][7];
static double *index_rows[8];
void reset_model(void) {memset(&model, 0, sizeof(model));}
int set_parameter(const char *name, double value) {struct xs_3 *p=&model;
if (!strcmp(name, "arm")) {p->arm=value; return 1;}
if (!strcmp(name, "ca_x0")) {p->ca_x0=value; return 1;}
if (!strcmp(name, "ca_x1")) {p->ca_x1=value; return 1;}
if (!strcmp(name, "ca_y0")) {p->ca_y0=value; return 1;}
if (!strcmp(name, "ca_y1")) {p->ca_y1=value; return 1;}
if (!strcmp(name, "chiprot")) {p->chiprot=value; return 1;}
if (!strcmp(name, "chipx")) {p->chipx=value; return 1;}
if (!strcmp(name, "chipy")) {p->chipy=value; return 1;}
if (!strcmp(name, "cmup1")) {p->cmup1=value; return 1;}
if (!strcmp(name, "config_mjd")) {p->config_mjd=value; return 1;}
if (!strcmp(name, "d2_x1")) {p->d2_x1=value; return 1;}
if (!strcmp(name, "d2_x2")) {p->d2_x2=value; return 1;}
if (!strcmp(name, "d2_x3")) {p->d2_x3=value; return 1;}
if (!strcmp(name, "d2_y1x0")) {p->d2_y1x0=value; return 1;}
if (!strcmp(name, "d2_y1x1")) {p->d2_y1x1=value; return 1;}
if (!strcmp(name, "d2_y1x2")) {p->d2_y1x2=value; return 1;}
if (!strcmp(name, "d2_y1x3")) {p->d2_y1x3=value; return 1;}
if (!strcmp(name, "d2_y2x0")) {p->d2_y2x0=value; return 1;}
if (!strcmp(name, "d2_y2x1")) {p->d2_y2x1=value; return 1;}
if (!strcmp(name, "d2_y2x2")) {p->d2_y2x2=value; return 1;}
if (!strcmp(name, "d2_y2x3")) {p->d2_y2x3=value; return 1;}
if (!strcmp(name, "d2_y3x0")) {p->d2_y3x0=value; return 1;}
if (!strcmp(name, "d2_y3x1")) {p->d2_y3x1=value; return 1;}
if (!strcmp(name, "d2_y3x2")) {p->d2_y3x2=value; return 1;}
if (!strcmp(name, "d2_y3x3")) {p->d2_y3x3=value; return 1;}
if (!strcmp(name, "es_s")) {p->es_s=value; return 1;}
if (!strcmp(name, "es_w")) {p->es_w=value; return 1;}
if (!strcmp(name, "es_x")) {p->es_x=value; return 1;}
if (!strcmp(name, "es_y")) {p->es_y=value; return 1;}
if (!strcmp(name, "fcol")) {p->fcol=value; return 1;}
if (!strcmp(name, "fdet")) {p->fdet=value; return 1;}
if (!strcmp(name, "flipx")) {p->flipx=value; return 1;}
if (!strcmp(name, "flipy")) {p->flipy=value; return 1;}
if (!strcmp(name, "mud")) {p->mud=value; return 1;}
if (!strcmp(name, "mues")) {p->mues=value; return 1;}
if (!strcmp(name, "mug")) {p->mug=value; return 1;}
if (!strcmp(name, "mup1")) {p->mup1=value; return 1;}
if (!strcmp(name, "mup2")) {p->mup2=value; return 1;}
if (!strcmp(name, "mup3")) {p->mup3=value; return 1;}
if (!strcmp(name, "mup4")) {p->mup4=value; return 1;}
if (!strcmp(name, "mup5")) {p->mup5=value; return 1;}
if (!strcmp(name, "mup6")) {p->mup6=value; return 1;}
if (!strcmp(name, "nud")) {p->nud=value; return 1;}
if (!strcmp(name, "nues")) {p->nues=value; return 1;}
if (!strcmp(name, "nug")) {p->nug=value; return 1;}
if (!strcmp(name, "nup1")) {p->nup1=value; return 1;}
if (!strcmp(name, "nup2")) {p->nup2=value; return 1;}
if (!strcmp(name, "nup3")) {p->nup3=value; return 1;}
if (!strcmp(name, "nup4")) {p->nup4=value; return 1;}
if (!strcmp(name, "nup5")) {p->nup5=value; return 1;}
if (!strcmp(name, "nup6")) {p->nup6=value; return 1;}
if (!strcmp(name, "offx")) {p->offx=value; return 1;}
if (!strcmp(name, "offy")) {p->offy=value; return 1;}
if (!strcmp(name, "pc4_x_x2y2")) {p->pc4_x_x2y2=value; return 1;}
if (!strcmp(name, "pc4_x_x3y")) {p->pc4_x_x3y=value; return 1;}
if (!strcmp(name, "pc4_x_x4")) {p->pc4_x_x4=value; return 1;}
if (!strcmp(name, "pc4_x_xy3")) {p->pc4_x_xy3=value; return 1;}
if (!strcmp(name, "pc4_x_y4")) {p->pc4_x_y4=value; return 1;}
if (!strcmp(name, "pc4_y_x2y2")) {p->pc4_y_x2y2=value; return 1;}
if (!strcmp(name, "pc4_y_x3y")) {p->pc4_y_x3y=value; return 1;}
if (!strcmp(name, "pc4_y_x4")) {p->pc4_y_x4=value; return 1;}
if (!strcmp(name, "pc4_y_xy3")) {p->pc4_y_xy3=value; return 1;}
if (!strcmp(name, "pc4_y_y4")) {p->pc4_y_y4=value; return 1;}
if (!strcmp(name, "pc_x_x1")) {p->pc_x_x1=value; return 1;}
if (!strcmp(name, "pc_x_x2y")) {p->pc_x_x2y=value; return 1;}
if (!strcmp(name, "pc_x_x3")) {p->pc_x_x3=value; return 1;}
if (!strcmp(name, "pc_x_xx")) {p->pc_x_xx=value; return 1;}
if (!strcmp(name, "pc_x_xy")) {p->pc_x_xy=value; return 1;}
if (!strcmp(name, "pc_x_y1")) {p->pc_x_y1=value; return 1;}
if (!strcmp(name, "pc_x_y2x")) {p->pc_x_y2x=value; return 1;}
if (!strcmp(name, "pc_x_y3")) {p->pc_x_y3=value; return 1;}
if (!strcmp(name, "pc_x_yy")) {p->pc_x_yy=value; return 1;}
if (!strcmp(name, "pc_y_x1")) {p->pc_y_x1=value; return 1;}
if (!strcmp(name, "pc_y_x2y")) {p->pc_y_x2y=value; return 1;}
if (!strcmp(name, "pc_y_x3")) {p->pc_y_x3=value; return 1;}
if (!strcmp(name, "pc_y_xx")) {p->pc_y_xx=value; return 1;}
if (!strcmp(name, "pc_y_xy")) {p->pc_y_xy=value; return 1;}
if (!strcmp(name, "pc_y_y1")) {p->pc_y_y1=value; return 1;}
if (!strcmp(name, "pc_y_y2x")) {p->pc_y_y2x=value; return 1;}
if (!strcmp(name, "pc_y_y3")) {p->pc_y_y3=value; return 1;}
if (!strcmp(name, "pc_y_yy")) {p->pc_y_yy=value; return 1;}
if (!strcmp(name, "pix")) {p->pix=value; return 1;}
if (!strcmp(name, "sg")) {p->sg=value; return 1;}
if (!strcmp(name, "slit[0]")) {p->slit[0]=value; return 1;}
if (!strcmp(name, "slit[1]")) {p->slit[1]=value; return 1;}
if (!strcmp(name, "slit[2]")) {p->slit[2]=value; return 1;}
if (!strcmp(name, "slit[3]")) {p->slit[3]=value; return 1;}
if (!strcmp(name, "slit[4]")) {p->slit[4]=value; return 1;}
if (!strcmp(name, "slit[5]")) {p->slit[5]=value; return 1;}
if (!strcmp(name, "slit[6]")) {p->slit[6]=value; return 1;}
if (!strcmp(name, "slit[7]")) {p->slit[7]=value; return 1;}
if (!strcmp(name, "slit[8]")) {p->slit[8]=value; return 1;}
if (!strcmp(name, "slit_scale")) {p->slit_scale=value; return 1;}
if (!strcmp(name, "t_ir_p2")) {p->t_ir_p2=value; return 1;}
if (!strcmp(name, "t_ir_p3")) {p->t_ir_p3=value; return 1;}
if (!strcmp(name, "taud")) {p->taud=value; return 1;}
if (!strcmp(name, "taues")) {p->taues=value; return 1;}
if (!strcmp(name, "taug")) {p->taug=value; return 1;}
if (!strcmp(name, "taup1")) {p->taup1=value; return 1;}
if (!strcmp(name, "taup2")) {p->taup2=value; return 1;}
if (!strcmp(name, "taup3")) {p->taup3=value; return 1;}
if (!strcmp(name, "taup4")) {p->taup4=value; return 1;}
if (!strcmp(name, "taup5")) {p->taup5=value; return 1;}
if (!strcmp(name, "taup6")) {p->taup6=value; return 1;}
if (!strcmp(name, "temper")) {p->temper=value; return 1;}
return 0;}

void initialize_model(int nx, int ny, double xcorr, double ycorr) {
    model.ASIZE=nx; model.BSIZE=ny; model.SIZE=ny;
    model.chipxpix=nx; model.chipypix=ny;
    model.xsize_corr=xcorr; model.ysize_corr=ycorr;
    model.pix_X=model.pix; model.pix_Y=model.pix;
    for (int i=0; i<8; i++) index_rows[i]=indices[i];
    double temperature = model.arm==2 ? model.t_ir_p2 : model.temper;
    xsh_ref_ind_read(model.arm, index_rows, temperature);
    xsh_3_init(&model);
}
void evaluate(int n, double *wave_nm, double *slit_as, double eta_as, int order, double *x, double *y) {
    double original_es_x=model.es_x;
    model.morder=order;
    model.es_x=original_es_x + eta_as * model.slit_scale;
    for (int i=0; i<n; i++) {
        model.es_y_tot=model.es_y+slit_as[i]*model.slit_scale;
        xsh_3_eval(wave_nm[i]/1000000.0, order, index_rows, &model);
        xsh_3_detpix(&model);
        // Unclipped version of ESO's pixel transformation, pixel centers at 0.5.
        if (model.arm != 2) {
            x[i]=0.5*(model.chipxpix+model.xsize_corr)-(model.ydet-model.chipy)/model.pix_X;
            y[i]=0.5*(model.chipypix+model.ysize_corr)-(model.xdet-model.chipx)/model.pix_Y;
        } else {
            x[i]=0.5*(model.chipxpix+model.xsize_corr)+(model.xdet-model.chipx)/model.pix_X;
            y[i]=0.5*(model.chipypix+model.ysize_corr)+(model.ydet-model.chipy)/model.pix_Y;
        }
    }
    model.es_x=original_es_x;
}
