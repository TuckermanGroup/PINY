/*=======================================================================================*/
#include "standard_include.h"
#include "../proto_defs/proto_friend_lib_entry.h"


#define SQR(a) ((a)*(a))
#define SQR6(a) ((a)*(a)*(a)*(a)*(a)*(a))
#define SQR7(a) ((a)*(a)*(a)*(a)*(a)*(a)*(a))
#define SQR8(a) ((a)*(a)*(a)*(a)*(a)*(a)*(a)*(a))
#define SQR9(a) ((a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a))
#define SQR10(a) ((a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a))
#define SQR11(a) ((a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a)*(a))
#define SQR(a) ((a)*(a))
#define fact 45.0/sqrt(1120.0)
#define PI 3.14159265358979323846264338327
double length(double *);
void vprod(double *, double *, double *);
double sprod(double *,double * );
double g202(double );
double g2022(double );
double g224(double ,double ,double);
double g2242(double ,double ,double,double,double);
double v000(double);
double v022(double);
double v224(double);
double H2H2pot(double , double , double , double );
double H2H2pot2(double , double , double , double,double ,double);
double Dv000DR(double );
double Dv022DR(double );
double Dv224DR(double );
double Dg202Dth(double );
double Dg202Dcosth(double);
double Dg224Dcosth(double ,double,double,int);
double Dg224Dsinth(double ,double,double,int);
double Dg224Dphi(double ,double,double);
double Dg224Dphit(double ,double,double);
double Dg224D2phi(double theta1,double theta2,double phi);
double dRdX1(double ,double ,double ,double ,double );

// ----------------------------------------------------SUBROUTINES---------------

// -----------------------------to get internal coordinates----------------------

void H2_pot_for(double *r1, double *r2, double *r3, double *r4, double *ff1,double *ff2, double *ff3, double *ff4,
                double *E_pot)
{
 int i;
 double *vFH,*vHH,*vCC,lvCC,*FHC,*FHH,*HHC;
 double rr,R,theta1,theta2,phi,E_pot_now;
 double *n1,*n2,*n3;
 double thetFHtemp,fact2,sthetHH,sthetFH;
 double temp_r3_r4;
 double rFH,rHH,ln1,ln2,thetFH,thetHH,thetHHtemp;
 double testX;
 double sign,signFH,signHH,test,tau;
 double dVdX1R, dVdX1theta1,dVdX1theta2,dVdx1phi;
 double dVdX1costheta1,dVdX1costheta2,dVdX1sintheta1,dVdX1sintheta2,dVdX1phi,dVdX12phi,V_tot_dx1;
 double dVdX2theta1,dVdX2theta2,dVdx2phi;
 double dtheta1_dx,dtheta2_dx,dtheta1_dx2,dtheta2_dx2;
 double dcostheta1_dx1,dcostheta2_dx1,dsintheta1_dx1,dsintheta2_dx1;
 double dcosphi_dx1,dRdX;
 double dn1n2, dn1n2av,dphidx,totX1,totX2 ;
 double dn1n2x2,dn1n2avx2,dphidx2;
 double dn1n2a,tag,dn1n2av2;
 double c1,c2,c3,h1,h2,h3,f1,f2,f3;
 double delta1;
 double dVdR,dVdcostheta1,dVdcostheta2,dVdsintheta1,dVdsintheta2,dVdphi;
 double dVdX2R,dVdX2costheta1,dVdX2costheta2,dVdX2sintheta1,dVdX2sintheta2,dVdX2phi,V_tot_dx2;
 double dVdX3R,dVdX3costheta1,dVdX3costheta2,dVdX3sintheta1,dVdX3sintheta2,dVdX3phi,V_tot_dx3;
 double dVdX4R,dVdX4costheta1,dVdX4costheta2,dVdX4sintheta1,dVdX4sintheta2,dVdX4phi,V_tot_dx4;
 double dRdY,dVdY1R,dVdY1costheta1,dVdY1costheta2,dVdY1sintheta1,dVdY1sintheta2,dVdY1phi,V_tot_dy1;
 double dVdY2R,dVdY2costheta1,dVdY2costheta2,dVdY2sintheta1,dVdY2sintheta2,dVdY2phi,V_tot_dy2;
 double dVdY3R,dVdY3costheta1,dVdY3costheta2,dVdY3sintheta1,dVdY3sintheta2,dVdY3phi,V_tot_dy3;
 double dVdY4R,dVdY4costheta1,dVdY4costheta2,dVdY4sintheta1,dVdY4sintheta2,dVdY4phi,V_tot_dy4;
 double dRdZ,dVdZ1R,dVdZ1costheta1,dVdZ1costheta2,dVdZ1sintheta1,dVdZ1sintheta2,dVdZ1phi,V_tot_dz1;
 double dVdZ2R,dVdZ2costheta1,dVdZ2costheta2,dVdZ2sintheta1,dVdZ2sintheta2,dVdZ2phi,V_tot_dz2;
 double dVdZ3R,dVdZ3costheta1,dVdZ3costheta2,dVdZ3sintheta1,dVdZ3sintheta2,dVdZ3phi,V_tot_dz3;
 double dVdZ4R,dVdZ4costheta1,dVdZ4costheta2,dVdZ4sintheta1,dVdZ4sintheta2,dVdZ4phi,V_tot_dz4;

// -----------------------------to get internal coordinates----------------------

 vFH = (double *) cmalloc(3*sizeof(double));
 vHH = (double *) cmalloc(3*sizeof(double));
 vCC = (double *) cmalloc(3*sizeof(double));
 FHC = (double *) cmalloc(3*sizeof(double));
 FHH = (double *) cmalloc(3*sizeof(double));
 HHC = (double *) cmalloc(3*sizeof(double));

 n1 = (double *) cmalloc(3*sizeof(double));
 n2 = (double *) cmalloc(3*sizeof(double));
 n3 = (double *) cmalloc(3*sizeof(double));

// Center of mass coordinates:
for(i=0;i<=2;i++)
{
  FHC[i]=0.5*(r1[i]+r2[i]);
  HHC[i]=0.5*(r3[i]+r4[i]);
//printf("\n r1[i] %15.14f, r2[i] %15.14f, r3[i] %15.14f, r4[i] %15.14f\n",r1[i],r2[i],r3[i], r4[i]);
}

for(i=0;i<=2;i++)
{
  vFH[i]=r2[i]-r1[i];
  vHH[i]=r4[i]-r3[i];
  vCC[i]=HHC[i]-FHC[i];
}
  vprod(vCC,vHH,n1);
  
  vprod(vCC,vFH,n2);
  vprod(n2,n1,n3);
  sign=sprod(n3,vCC);
// GET vector lengths

for(i=0;i<=2;i++)
{
}
 // GET normal vectors

  rFH=length(vFH);
  rHH=length(vHH);
//rHH is eq. to Fortran code
  rr= length(vCC);
  ln1=length(n1);
  ln2=length(n2);
// GET angles
  thetFHtemp=-sprod(vFH,vCC)/(rFH*length(vCC));
  thetFH=acos(sprod(vFH,vCC)/(rFH*length(vCC)));
  thetHHtemp=sprod(vHH,vCC)/(rHH*length(vCC));
  sthetHH=ln1/(rHH*length(vCC));
  sthetFH=ln2/(rFH*length(vCC));
  thetHH=acos(thetHHtemp);
// IF for acos(1)
//  thetHHtemp=1;
  test=sprod(n1,n2)/(ln1*ln2);
  if(test<=-1.0) test=-1.0;
  if(test>=1.0)  test=1.0;
  tau=acos(test);
  if(sign<=0.0) tau=2.0*PI-tau;
//  if(sign<=0.0) tau=2.0*M_PI-tau;
// assigning final values to pointers:
  R=rr;
  theta1=thetHH;
//  theta2= M_PI-thetFH;
  theta2= PI-thetFH;

  if(sign==0) 
  {
    test=0;
    tau=0;
    phi=tau;
  }
  else 
  {
    phi=tau;
  }
// make simpler the vectorproducts
c1 = (r3[0]+r4[0]-r2[0]-r1[0])/2;
c2 = (r3[1]+r4[1]-r2[1]-r1[1])/2;
c3 = (r3[2]+r4[2]-r2[2]-r1[2])/2;
h1 = (r4[0]-r3[0]);
h2 = (r4[1]-r3[1]);
h3 = (r4[2]-r3[2]);

f1 = (r2[0]-r1[0]);
f2 = (r2[1]-r1[1]);
f3 = (r2[2]-r1[2]);
 
//    COMPUTING DERIVATIVES

// FIRST things we always need: 
/* old with transformation to angles first and then use ================================
   of arcus functions
dVdR = Dv000DR(rr) + Dv022DR(rr) * g202(theta1) +  Dv022DR(rr) * g202(theta2) + \
            Dv224DR(rr) * g224(theta1,theta2,phi);
dVdcostheta1 =  v022(rr)*Dg202Dcosth(theta1) +v224(rr)* Dg224Dcosth(theta1,theta2,phi,1);
dVdcostheta2 = -( v022(rr)*Dg202Dcosth(theta2) +v224(rr)* Dg224Dcosth(theta2,theta1,phi,2));

dVdsintheta1 = v224(rr)* Dg224Dsinth(theta1,theta2,phi,1);
dVdsintheta2 = v224(rr)* Dg224Dsinth(theta2,theta1,phi,2);
dVdphi       = v224(rr)* Dg224Dphit(theta1,theta2,phi);
*/
dVdR = Dv000DR(rr) + Dv022DR(rr) * g2022(thetHHtemp) +  Dv022DR(rr) * g2022(thetFHtemp) + \
            Dv224DR(rr) * g2242(thetHHtemp,thetFHtemp,sthetHH,sthetFH,test);
dVdcostheta1 =  v022(rr)*Dg202Dcosth(theta1) +v224(rr)* Dg224Dcosth(theta1,theta2,phi,1);
dVdcostheta2 = -( v022(rr)*Dg202Dcosth(theta2) +v224(rr)* Dg224Dcosth(theta2,theta1,phi,2));

dVdsintheta1 = v224(rr)* Dg224Dsinth(theta1,theta2,phi,1);
dVdsintheta2 = v224(rr)* Dg224Dsinth(theta2,theta1,phi,2);
dVdphi       = v224(rr)* Dg224Dphit(theta1,theta2,phi);
// derivatives wrt. x1
dRdX = dRdX1(r1[0],r2[0],r3[0],r4[0],rr);

//temporary variables:

temp_r3_r4 = (r3[0]-r4[0])/2;

//f(ln1==0)
//{
//  ln1=0.001;
//}
dVdX1R = dVdR * dRdX;

dVdX1costheta1 = ((temp_r3_r4*rHH*rr - rHH*sprod(vHH,vCC)*dRdX)/(rHH*rr*rHH*rr) )* dVdcostheta1;

dVdX1costheta2 = ( ( (-r3[0]/2-r4[0]/2+r1[0])*rFH*rr -((sprod(vFH,vCC))*( ((-r2[0]+r1[0])*rr/rFH) + rFH*dRdX) ))/(rFH*rFH*rr*rr )) * dVdcostheta2;

dVdX1sintheta1 =(( (((h1*c3-h3*c1)*h3+(h1*c2-h2*c1)*h2)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdX) )/(rHH*rr*rHH*rr)) * dVdsintheta1;


dVdX1sintheta2 = (((( ((-f3*c1+f1*c3)*(0.5*f3-c3)) + ((-f1*c2+f2*c1)*(-0.5*f2+c2))  )*rFH*rr/ln2) - ln2*(((-r2[0]+r1[0])*rr/rFH ) + rFH*dRdX))/(rFH*rFH*rr*rr) )* dVdsintheta2;


dVdX1phi = (((-c3*c3*h1 +0.5*c3*h1*f3+h3*c3*(0.5*f1+c1) -c1*h3*f3 -h2*f2*c1 +h2*c2*(0.5*f1+c1) +\
               0.5*c2*h1*f2 -c2*c2*h1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c3*h1-c1*h3)*h3-(c1*h2-c2*h1)*h2)*ln2/(2*ln1) +( ((c3*f1-c1*f3)*(-c3+0.5*f3)+(c1*f2-c2*f1)*(-0.5*f2+c2))*ln1/(ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdX1phi=0;
    dVdX12phi=0;
  }

if((ln1==0) )
{
  dVdX1sintheta1=0;
  dVdX1phi =0;
}
if(ln2==0)
{
  dVdX1sintheta2=0;
  dVdX1phi=0;
}
V_tot_dx1 = dVdX1R + dVdX1costheta1 + dVdX1costheta2 + dVdX1sintheta1 + dVdX1sintheta2 + dVdX1phi;
// correct result
 

// ----------------------------------------------------- DVdx2
dRdX = dRdX1(r1[0],r2[0],r3[0],r4[0],rr);

//temporary variables:

temp_r3_r4 = (r3[0]-r4[0])/2;

dVdX2R = dVdR * dRdX;

dVdX2costheta1 = ((temp_r3_r4*rHH*rr - rHH*sprod(vHH,vCC)*dRdX)/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdX2costheta2 = (( ( (r3[0]/2+r4[0]/2-r2[0])*rFH*rr -((sprod(vFH,vCC))*( ((r2[0]-r1[0])*rr/rFH) + rFH*dRdX) ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdX2sintheta1 =(( (((h1*c3-h3*c1)*h3+(h1*c2-h2*c1)*h2)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdX) )/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdX2sintheta2 = (((( ((-f3*c1+f1*c3)*(0.5*f3+c3)) + ((-f1*c2+f2*c1)*(-0.5*f2-c2))  )*rFH*rr/ln2) - ln2*(((r2[0]-r1[0])*rr/rFH ) + rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdX2phi = (((c3*c3*h1 +0.5*c3*h1*f3-h3*c3*(-0.5*f1+c1) -c1*h3*f3 -h2*f2*c1 -h2*c2*(-0.5*f1+c1) +\
               0.5*c2*h1*f2 +c2*c2*h1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c3*h1-c1*h3)*h3-(c1*h2-c2*h1)*h2)*ln2/(2*ln1) +( ((c3*f1-c1*f3)*(c3+0.5*f3)+(c1*f2-c2*f1)*(-0.5*f2-c2))*ln1/(ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
//dVdX12phi= v224(rr)*Dg224D2phi(theta1,theta2,phi) * (4*(sprod(n1,n2))/(ln1*ln2)*dcosphi_dx1);

  if(sign==0)
  {
    dVdX2phi=0;
  }

if((ln1==0) )
{
  dVdX2sintheta1=0;
  dVdX2phi =0;
}
if(ln2==0)
{
  dVdX2sintheta2=0;
  dVdX2phi=0;
}

V_tot_dx2 = dVdX2R + dVdX2costheta1 + dVdX2costheta2 + dVdX2sintheta1 + dVdX2sintheta2 + dVdX2phi;
//correct result


// ----------------------------------------------------- DVdx3
dRdX = dRdX1(r3[0],r4[0],r2[0],r1[0],rr);

dVdX3R = dVdR * dRdX;

dVdX3costheta1 = (((r1[0]/2+r2[0]/2-r3[0])*rHH*rr - sprod(vHH,vCC)*((r3[0]-r4[0])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdX3costheta2 = (( ( (r2[0]/2-r1[0]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdX3sintheta1 =(( (((c3*h1-c1*h3)*(-c3-0.5*h3)+(c1*h2-c2*h1)*(0.5*h2+c2))*rHH*rr/(ln1)) - (ln1*( (r3[0]-r4[0])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdX3sintheta2 = (((( ((c3*f1-c1*f3)*(-f3)) + ((c1*f2-c2*f1)*(f2))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdX3phi = (((-c3*c3*f1 -c3*f3*(-c1+0.5*h1) -0.5*h3*c3*f1 +c1*h3*f3 +c1*h2*f2 -0.5*h2*c2*f1 -\
               c2*f2*(0.5*h1-c1) -c2*c2*f1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c3*h1-c1*h3)*(-c3-0.5*h3)+(c1*h2-c2*h1)*(0.5*h2+c2))*ln2/(ln1) +( (((c3*f1-c1*f3)*(-f3)) + ((c1*f2-c2*f1)*(f2))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdX3phi=0;
  }

if((ln1==0) )
{
  dVdX3sintheta1=0;
  dVdX3phi =0;
}
if(ln2==0)
{
  dVdX3sintheta2=0;
  dVdX3phi=0;
}


V_tot_dx3 = dVdX3R + dVdX3costheta1 + dVdX3costheta2 + dVdX3sintheta1 + dVdX3sintheta2 + dVdX3phi;
// correct -checked


// ----------------------------------------------------- DVdx4
dRdX = dRdX1(r3[0],r4[0],r2[0],r1[0],rr);

dVdX4R = dVdR * dRdX;

dVdX4costheta1 = (((-r1[0]/2-r2[0]/2+r4[0])*rHH*rr - sprod(vHH,vCC)*((r4[0]-r3[0])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdX4costheta2 = (( ( (r2[0]/2-r1[0]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdX4sintheta1 =(( (((c3*h1-c1*h3)*(+c3-0.5*h3)+(c1*h2-c2*h1)*(0.5*h2-c2))*rHH*rr/(ln1)) - (ln1*( (r4[0]-r3[0])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdX4sintheta2 = (((( ((c3*f1-c1*f3)*(-f3)) + ((c1*f2-c2*f1)*(f2))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdX4phi = (((c3*c3*f1 -c3*f3*(c1+0.5*h1) -0.5*h3*c3*f1 +c1*h3*f3 +c1*h2*f2 -0.5*h2*c2*f1 -\
               c2*f2*(0.5*h1+c1) +c2*c2*f1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c3*h1-c1*h3)*(+c3-0.5*h3)+(c1*h2-c2*h1)*(0.5*h2-c2))*ln2/(ln1) +( (((c3*f1-c1*f3)*(-f3)) + ((c1*f2-c2*f1)*(f2))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdX4phi=0;
  }

if((ln1==0) )
{
  dVdX4sintheta1=0;
  dVdX4phi =0;
}
if(ln2==0)
{
  dVdX4sintheta2=0;
  dVdX4phi=0;
}


V_tot_dx4 = dVdX4R + dVdX4costheta1 + dVdX4costheta2 + dVdX4sintheta1 + dVdX4sintheta2 + dVdX4phi;
// checked and working

// --------------------------------------------    dVdY1:

dRdY = dRdX1(r1[1],r2[1],r3[1],r4[1],rr);

dVdY1R = dVdR * dRdY;

dVdY1costheta1 =  dVdcostheta1 * (((r3[1]/2-r4[1]/2)*rHH*rr - rHH*sprod(vHH,vCC)*dRdY)/(rHH*rr*rHH*rr) ) ;

dVdY1costheta2 = dVdcostheta2 * ( ( (-r3[1]/2-r4[1]/2+r1[1])*rFH*rr -((sprod(vFH,vCC))*( ((-r2[1]+r1[1])*rr/rFH) + rFH*dRdY) ))/(rFH*rFH*rr*rr )) ;


dVdY1sintheta1 =dVdsintheta1 * (( (((c3*h2-h3*c2)*h3+(h2*c1-c2*h1)*h1)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdY) )/(rHH*rr*rHH*rr)) ;


dVdY1sintheta2 = dVdsintheta2 * (((( ((f3*c2-c3*f2)*(-0.5*f3+c3)) + ((-f1*c2+f2*c1)*(0.5*f1-c1))  )*rFH*rr/ln2) - ln2*(((-r2[1]+r1[1])*rr/rFH ) + rFH*dRdY))/(rFH*rFH*rr*rr) );


dVdY1phi=dVdphi * ((-c2*h3*f3 +h3*c3*(0.5*f2+c2) + 0.5*c3*h2*f3 -c3*c3*h2 -c1*c1*h2 +0.5*c1*h2*f1 +\
                h1*c1*(0.5*f2+c2) -c2*h1*f1 )*ln1*ln2 - \
               sprod(n1,n2)*((((c3*h2-h3*c2)*h3+(h2*c1-c2*h1)*h1) )*ln2/(2*ln1) +( ( ((f3*c2-c3*f2)*(-0.5*f3+c3)) + ((-f1*c2+f2*c1)*(0.5*f1-c1)))*ln1/(ln2)   )) )        \
              /  (ln1*ln1*ln2*ln2);

  if(sign==0)
  {
    dVdY1phi=0;
  }

if((ln1==0) )
{
  dVdY1sintheta1=0;
  dVdY1phi =0;
}
if(ln2==0)
{
  dVdY1sintheta2=0;
  dVdY1phi=0;
}

V_tot_dy1 = dVdY1R + dVdY1costheta1 + dVdY1costheta2 + dVdY1sintheta1 + dVdY1sintheta2 + dVdY1phi;

// correct result


// ----------------------------------------------------dVdy2

dRdY = dRdX1(r1[1],r2[1],r3[1],r4[1],rr);

dVdY2R = dVdR * dRdY;

dVdY2costheta1 =  dVdcostheta1 * (((r3[1]/2-r4[1]/2)*rHH*rr - rHH*sprod(vHH,vCC)*dRdY)/(rHH*rr*rHH*rr) ) ;

dVdY2costheta2 = dVdcostheta2 * ( ( (r3[1]/2+r4[1]/2-r2[1])*rFH*rr -((sprod(vFH,vCC))*( ((r2[1]-r1[1])*rr/rFH) + rFH*dRdY) ))/(rFH*rFH*rr*rr )) ;


dVdY2sintheta1 =dVdsintheta1 * (( (((c3*h2-h3*c2)*h3+(h2*c1-c2*h1)*h1)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdY) )/(rHH*rr*rHH*rr)) ;


dVdY2sintheta2 = dVdsintheta2 * (((( ((f3*c2-c3*f2)*(-0.5*f3-c3)) + ((-f1*c2+f2*c1)*(0.5*f1+c1))  )*rFH*rr/ln2) - ln2*(((r2[1]-r1[1])*rr/rFH ) + rFH*dRdY))/(rFH*rFH*rr*rr) );


dVdY2phi=dVdphi * ((-c2*h3*f3 -h3*c3*(-0.5*f2+c2) + 0.5*c3*h2*f3 +c3*c3*h2 +c1*c1*h2 +0.5*c1*h2*f1-\
                h1*c1*(-0.5*f2+c2) -c2*h1*f1 )*ln1*ln2 - \
               sprod(n1,n2)*((((c3*h2-h3*c2)*h3+(h2*c1-c2*h1)*h1) )*ln2/(2*ln1) +( ( ((f3*c2-c3*f2)*(-0.5*f3-c3)) + ((-f1*c2+f2*c1)*(0.5*f1+c1)))*ln1/(ln2)   )) )        \
              /  (ln1*ln1*ln2*ln2);

  if(sign==0)
  {
    dVdY2phi=0;
  }

if((ln1==0) )
{
  dVdY2sintheta1=0;
  dVdY2phi =0;
}
if(ln2==0)
{
  dVdY2sintheta2=0;
  dVdY2phi=0;
}


V_tot_dy2 = dVdY2R + dVdY2costheta1 + dVdY2costheta2 + dVdY2sintheta1 + dVdY2sintheta2 + dVdY2phi;
// checked -- working

// ----------------------------------------------------- DVdy3
dRdX = dRdX1(r3[1],r4[1],r2[1],r1[1],rr);

dVdY3R = dVdR * dRdX;

dVdY3costheta1 = (((r1[1]/2+r2[1]/2-r3[1])*rHH*rr - sprod(vHH,vCC)*((r3[1]-r4[1])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdY3costheta2 = (( ( (r2[1]/2-r1[1]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdY3sintheta1 =(( (((c2*h3-c3*h2)*(c3+0.5*h3)+(c1*h2-c2*h1)*(-0.5*h1-c1))*rHH*rr/(ln1)) - (ln1*( (r3[1]-r4[1])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdY3sintheta2 = (((( ((c2*f3-c3*f2)*(f3)) + ((c1*f2-c2*f1)*(-f1))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdY3phi = (((c2*h3*f3 -0.5*h3*c3*f2 -c3*f3*(0.5*h2-c2) -c3*c3*f2 -c1*c1*f2 -c1*f1*(0.5*h2-c2) -\
               0.5*h1*c1*f2 +c2*h1*f1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c2*h3-c3*h2)*(c3+0.5*h3)+(c1*h2-c2*h1)*(-0.5*h1-c1))*ln2/(ln1) +( (((c2*f3-c3*f2)*(f3)) + ((c1*f2-c2*f1)*(-f1))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdY3phi=0;
  }

if((ln1==0) )
{
  dVdY3sintheta1=0;
  dVdY3phi =0;
}
if(ln2==0)
{
  dVdY3sintheta2=0;
  dVdY3phi=0;
}


V_tot_dy3 = dVdY3R + dVdY3costheta1 + dVdY3costheta2 + dVdY3sintheta1 + dVdY3sintheta2 + dVdY3phi;
// correct and checked


// ----------------------------------------------------- DVdy4
dRdX = dRdX1(r3[1],r4[1],r2[1],r1[1],rr);

dVdY4R = dVdR * dRdX;

dVdY4costheta1 = (((r4[1]-r2[1]/2-r1[1]/2)*rHH*rr - sprod(vHH,vCC)*((r4[1]-r3[1])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdY4costheta2 = (( ( (r2[1]/2-r1[1]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdY4sintheta1 =(( (((c2*h3-c3*h2)*(-c3+0.5*h3)+(c1*h2-c2*h1)*(-0.5*h1+c1))*rHH*rr/(ln1)) - (ln1*( (r4[1]-r3[1])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdY4sintheta2 = (((( ((c2*f3-c3*f2)*(f3)) + ((c1*f2-c2*f1)*(-f1))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdY4phi = (((c2*h3*f3 -0.5*h3*c3*f2 -c3*f3*(0.5*h2+c2) +c3*c3*f2 +c1*c1*f2 -c1*f1*(0.5*h2+c2) -\
               0.5*h1*c1*f2 +c2*h1*f1  )*ln1*ln2 - \
               sprod(n1,n2)*(((c2*h3-c3*h2)*(-c3+0.5*h3)+(c1*h2-c2*h1)*(-0.5*h1+c1))*ln2/(ln1) +( (((c2*f3-c3*f2)*(f3)) + ((c1*f2-c2*f1)*(-f1))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdY4phi=0;
  }

if((ln1==0) )
{
  dVdY4sintheta1=0;
  dVdY4phi =0;
}
if(ln2==0)
{
  dVdY4sintheta2=0;
  dVdY4phi=0;
}


V_tot_dy4 = dVdY4R + dVdY4costheta1 + dVdY4costheta2 + dVdY4sintheta1 + dVdY4sintheta2 + dVdY4phi;
// checked --correct



// ----------------------------------------------------dVdz1

dRdZ = dRdX1(r1[2],r2[2],r3[2],r4[2],rr);

dVdZ1R = dVdR * dRdZ;

dVdZ1costheta1 =  dVdcostheta1 * (((r3[2]/2-r4[2]/2)*rHH*rr - rHH*sprod(vHH,vCC)*dRdZ)/(rHH*rr*rHH*rr) ) ;

dVdZ1costheta2 = dVdcostheta2 * ( ( (-r3[2]/2-r4[2]/2+r1[2])*rFH*rr -((sprod(vFH,vCC))*( ((-r2[2]+r1[2])*rr/rFH) + rFH*dRdZ) ))/(rFH*rFH*rr*rr )) ;

dVdZ1sintheta1 =dVdsintheta1 * (( (((h2*c3-h3*c2)*(-h2)+(h3*c1-h1*c3)*h1)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdZ) )/(rHH*rr*rHH*rr)) ;

dVdZ1sintheta2 = dVdsintheta2 * (((( ((f3*c2-c3*f2)*(0.5*f2-c2)) + ((f1*c3-f3*c1)*(-0.5*f1+c1))  )*rFH*rr/ln2) - ln2*(((-r2[2]+r1[2])*rr/rFH ) + rFH*dRdZ))/(rFH*rFH*rr*rr) );

dVdZ1phi = dVdphi * ((-c2*h3*c2 +0.5*h3*c2*f2 + h2*c2*(c3+0.5*f3) -c3*h2*f2 -\
                       c3*h1*f1 +h1*c1*(c3+0.5*f3) + 0.5*c1*h3*f1 -c1*c1*h3 )*ln1*ln2 - \
               sprod(n1,n2)*((( ((h2*c3-h3*c2)*(-h2)+(h3*c1-h1*c3)*h1) ) )*ln2/(2*ln1) + \
           (((f3*c2-c3*f2)*(+0.5*f2-c2)) + ((f1*c3-f3*c1)*(-0.5*f1+c1))  )*ln1/(ln2)   ) )   \
              /  (ln1*ln1*ln2*ln2);
  if(sign==0)
  {
    dVdZ1phi=0;
  }

if((ln1==0) )
{
  dVdZ1sintheta1=0;
  dVdZ1phi =0;
}
if(ln2==0)
{
  dVdZ1sintheta2=0;
  dVdZ1phi=0;
}


V_tot_dz1 = dVdZ1R + dVdZ1costheta1 + dVdZ1costheta2 + dVdZ1sintheta1 + dVdZ1sintheta2 + dVdZ1phi;
// correct result


// ----------------------------------------------------dVdz2

dRdZ = dRdX1(r1[2],r2[2],r3[2],r4[2],rr);

dVdZ2R = dVdR * dRdZ;

dVdZ2costheta1 =  dVdcostheta1 * (((r3[2]/2-r4[2]/2)*rHH*rr - rHH*sprod(vHH,vCC)*dRdZ)/(rHH*rr*rHH*rr) ) ;

dVdZ2costheta2 = dVdcostheta2 * ( ( (r3[2]/2+r4[2]/2-r2[2])*rFH*rr -((sprod(vFH,vCC))*( ((r2[2]-r1[2])*rr/rFH) + rFH*dRdZ) ))/(rFH*rFH*rr*rr )) ;

dVdZ2sintheta1 =dVdsintheta1 * (( (((h2*c3-h3*c2)*(-h2)+(h3*c1-h1*c3)*h1)*rHH*rr/(2*ln1)) - (ln1*rHH*dRdZ) )/(rHH*rr*rHH*rr)) ;

dVdZ2sintheta2 = dVdsintheta2 * (((( ((f3*c2-c3*f2)*(0.5*f2+c2)) + ((f1*c3-f3*c1)*(-0.5*f1-c1))  )*rFH*rr/ln2) - ln2*(((r2[2]-r1[2])*rr/rFH ) + rFH*dRdZ))/(rFH*rFH*rr*rr) );

dVdZ2phi = dVdphi * ((c2*h3*c2 +0.5*h3*c2*f2 - h2*c2*(c3-0.5*f3) -c3*h2*f2 -\
                       c3*h1*f1 -h1*c1*(c3-0.5*f3) + 0.5*c1*h3*f1 +c1*c1*h3 )*ln1*ln2 - \
               sprod(n1,n2)*((( ((h2*c3-h3*c2)*(-h2)+(h3*c1-h1*c3)*h1) ) )*ln2/(2*ln1) + \
           (((f3*c2-c3*f2)*(+0.5*f2+c2)) + ((f1*c3-f3*c1)*(-0.5*f1-c1))  )*ln1/(ln2)   ) )   \
              /  (ln1*ln1*ln2*ln2);
  if(sign==0)
  {
    dVdZ2phi=0;
  }

if((ln1==0) )
{
  dVdZ2sintheta1=0;
  dVdZ2phi =0;
}
if(ln2==0)
{
  dVdZ2sintheta2=0;
  dVdZ2phi=0;
}


V_tot_dz2 = dVdZ2R + dVdZ2costheta1 + dVdZ2costheta2 + dVdZ2sintheta1 + dVdZ2sintheta2 + dVdZ2phi;
// correct result



// ----------------------------------------------------- DVdz3
dRdX = dRdX1(r3[2],r4[2],r2[2],r1[2],rr);

dVdZ3R = dVdR * dRdX;

dVdZ3costheta1 = (((r1[2]/2+r2[2]/2-r3[2])*rHH*rr - sprod(vHH,vCC)*((r3[2]-r4[2])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdZ3costheta2 = (( ( (r2[2]/2-r1[2]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdZ3sintheta1 =(( (((c2*h3-c3*h2)*(-c2-0.5*h2)+(c3*h1-c1*h3)*(0.5*h1+c1))*rHH*rr/(ln1)) - (ln1*( (r3[2]-r4[2])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdZ3sintheta2 = (((( ((c2*f3-c3*f2)*(-f2)) + ((c3*f1-c1*f3)*(f1))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdZ3phi = (((-c2*c2*f3 -c2*f2*(-c3+0.5*h3) -0.5*h2*c2*f3 +c3*h2*f2 +c3*h1*f1 -0.5*h1*c1*f3 -\
               c1*f1*(0.5*h3-c3) -c1*c1*f3  )*ln1*ln2 - \
               sprod(n1,n2)*(((c2*h3-c3*h2)*(-c2-0.5*h2)+(c3*h1-c1*h3)*(0.5*h1+c1))*ln2/(ln1) +( (((c2*f3-c3*f2)*(-f2)) + ((c3*f1-c1*f3)*(f1))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdZ3phi=0;
  }

if((ln1==0) )
{
  dVdZ3sintheta1=0;
  dVdZ3phi =0;
}
if(ln2==0)
{
  dVdZ3sintheta2=0;
  dVdZ3phi=0;
}


V_tot_dz3 = dVdZ3R + dVdZ3costheta1 + dVdZ3costheta2 + dVdZ3sintheta1 + dVdZ3sintheta2 + dVdZ3phi;
//correct --checked

// ----------------------------------------------------- DVdz4
dRdX = dRdX1(r3[2],r4[2],r2[2],r1[2],rr);

dVdZ4R = dVdR * dRdX;

dVdZ4costheta1 = (((-r1[2]/2-r2[2]/2+r4[2])*rHH*rr - sprod(vHH,vCC)*((r4[2]-r3[2])/rHH*rr +rHH*dRdX))/(rHH*rr*rHH*rr)) * dVdcostheta1;

dVdZ4costheta2 = (( ( (r2[2]/2-r1[2]/2)*rFH*rr -((sprod(vFH,vCC))*rFH*dRdX ))/(rFH*rFH*rr*rr ))) * dVdcostheta2 ;

dVdZ4sintheta1 =(( (((c2*h3-c3*h2)*(+c2-0.5*h2)+(c3*h1-c1*h3)*(0.5*h1-c1))*rHH*rr/(ln1)) - (ln1*( (r4[2]-r3[2])*rr/rHH+rHH*dRdX) ))/(rHH*rr*rHH*rr)) *dVdsintheta1 ;

dVdZ4sintheta2 = (((( ((c2*f3-c3*f2)*(-f2)) + ((c3*f1-c1*f3)*(f1))  )*rFH*rr/(2*ln2)) - ln2*(rFH*dRdX))/(rFH*rFH*rr*rr)) * dVdsintheta2;

dVdZ4phi = (((c2*c2*f3 -c2*f2*(c3+0.5*h3) -0.5*h2*c2*f3 +c3*h2*f2 +c3*h1*f1 -0.5*h1*c1*f3 -\
               c1*f1*(0.5*h3+c3) +c1*c1*f3  )*ln1*ln2 - \
               sprod(n1,n2)*(((c2*h3-c3*h2)*(+c2-0.5*h2)+(c3*h1-c1*h3)*(0.5*h1-c1))*ln2/(ln1) +( (((c2*f3-c3*f2)*(-f2)) + ((c3*f1-c1*f3)*(f1))   )*ln1/(2*ln2)   ) ))        \
              /  (ln1*ln1*ln2*ln2)) * dVdphi;
  if(sign==0)
  {
    dVdZ4phi=0;
  }

if((ln1==0) )
{
  dVdZ4sintheta1=0;
  dVdZ4phi =0;
}
if(ln2==0)
{
  dVdZ4sintheta2=0;
  dVdZ4phi=0;
}


V_tot_dz4 = dVdZ4R + dVdZ4costheta1 + dVdZ4costheta2 + dVdZ4sintheta1 + dVdZ4sintheta2 + dVdZ4phi;
// correct --checked

ff1[0]=-V_tot_dx1;
ff1[1]=-V_tot_dy1;
ff1[2]=-V_tot_dz1;
ff2[0]=-V_tot_dx2;
ff2[1]=-V_tot_dy2;
ff2[2]=-V_tot_dz2;
ff3[0]=-V_tot_dx3;
ff3[1]=-V_tot_dy3;
ff3[2]=-V_tot_dz3;
ff4[0]=-V_tot_dx4;
ff4[1]=-V_tot_dy4;
ff4[2]=-V_tot_dz4;


   E_pot_now = H2H2pot(R, theta1, theta2, phi);
   *E_pot = E_pot_now;

/* =================================================================*/
/* Free memory */
 
// printf("\n rFH %15.14f, rHH %15.14f, rr %15.14f ",rFH,rHH,rr);
// printf("\n thetHHtemp %15.14f, thetHH %15.14f, thetFHtemp %15.14f, thetFH %15.14f ",thetHHtemp,thetHH,thetFHtemp,thetFH);
// printf("\n sthetHH %15.14f, sthetFH %15.14f, phi %15.14f, test %15.14f ",sthetHH,sthetFH,phi,test);

   cfree(vFH);
   cfree(vHH);
   cfree(vCC);
   cfree(FHC);
   cfree(FHH);
   cfree(HHC);

   cfree(n1);
   cfree(n2);
   cfree(n3);

/* =================================================================*/
}  /* end function */
/* =================================================================*/




// -----------------------------------------FUNCTIONS NEEDED---------------------
// computes length of a vector
double length(double *vec)
{
 int i,j,k;
 double  length2;
 length2=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
 return length2;
}
// computes normal vector of two given vectors
void vprod(double *vec1, double *vec2, double *vec3)
{
 int i,j,k;
 double vect1[3],vect2[3],vect3[3];
 vec3[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
 vec3[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
 vec3[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
 return;
}
// computes the scalar product
double sprod(double *vec1, double *vec2)
{
 int i;
 double sprod2=0;
 for(i=0;i<=2;i++)
 {
   sprod2 += vec1[i]*vec2[i];
 }
 return sprod2;
}
// 

double g202(double theta)
{
 double g202;
 
   g202 = 2.5*( 3*cos(theta)*cos(theta) -1   );

 return g202;
}

double g2022(double theta)
{
 double g202;
 
   g202 = 2.5*( 3*theta*theta -1   );

 return g202;
}

double g2242(double theta1,double theta2, double sin1, double sin2, double phi)
{
  double g224;

  g224 = fact * ( 2.0*  (3.0*theta1*theta1 -1.0 )*( 3.0*theta2*theta2 -1.0 )    - 16.0* sin1*theta1*sin2*theta2*phi + sin1*sin1*sin2*sin2*(2*phi*phi - 1 ) );

  return g224;
}
double g224(double theta1,double theta2, double phi)
{
  double g224;

  g224 = fact * ( 2.0*  (3.0*cos(theta1)*cos(theta1) -1.0 )*( 3.0*cos(theta2)*cos(theta2) -1.0 )    - 16.0* sin(theta1)*cos(theta1)*sin(theta2)*cos(theta2)*cos(phi) + sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*cos(2.0*phi) );


  return g224;
}


double v000(double rr)
{
//  Computes the radial component v000 of the H2-H2 potential.
//  This is the isotropic contribution.  Distance and energies are in 
//  atomic units (bohr and hartree).
//**********************************************************************
   double v000;
   double alpha=0.655914, beta=1.018447, gamma=8.070618E-02;
   double rc=9.034308, c6=13.076837, c8=80.700360, c10=3687.082967, fcut;
   double rr6,rr8,rr10;
   if(rr<rc) 
   {
     fcut = exp(-((rc/rr-1.0)*(rc/rr-1.0)));
   } 
   else
   {
         fcut = 1.0;
   }
   v000 = exp(alpha - beta*rr - gamma*SQR(rr)) - \
          fcut*(c6/SQR6(rr) + c8/SQR8(rr) + c10/(SQR10(rr)));
//   v000 = exp(alpha - beta*rr - gamma*pow(rr,2)) - \
//          fcut*(c6/(pow(rr,6)) + c8/(pow(rr,8)) + c10/(pow(rr,10)));
   return v000;
}

// ***************************************DERIVATIVE   d(v000) / dR ************DERIVATIVE
double Dv000DR(double rr)
{
  double Dv000DR;
  double alpha=0.655914, beta=1.018447, gamma=8.070618E-02;
  double rc=9.034308, c6=13.076837, c8=80.700360, c10=3687.082967, fcut,dfcut;
   if(rr<rc)
   {
     fcut = exp(-((rc/rr-1.0)*(rc/rr-1.0)));
     dfcut = exp(-( (rc/rr-1.0) * (rc/rr-1.0) ))*2 *(rc/rr-1) * (rc/(rr*rr));
   Dv000DR = exp(alpha - beta*rr-gamma*rr*rr ) * (-beta -2*gamma*rr) - \
             dfcut * (c6/(pow(rr,6)) + c8/(pow(rr,8)) + c10/(pow(rr,10))) +\
             fcut  * (6*c6/(pow(rr,7)) + 8*c8/(pow(rr,9)) + 10*c10/(pow(rr,11))); 
   }
   else
   {
         fcut = 1.0;
      Dv000DR = exp(alpha - beta*rr-gamma*rr*rr )*(-beta -2*gamma*rr) +\
                fcut*(6*c6/(pow(rr,7)) + 8*c8/(pow(rr,9)) + 10*c10/(pow(rr,11)));
   }
   return Dv000DR;
}

double v022(double rr)
//  Computes the radial components, v022 and v202, of the H2-H2 
//  potential.  v022 = v202 due to symmetry.
//  Atomic units (bohr and hartree).
//***********************************************************************
{
  double rr6,rr8,rr10, v022;
  double alpha=-3.428227, beta=0.715011, gamma=0.100120, rc=8.422755;
  double c6=0.288396, c8=8.242595, c10=210.984075, fcut;
  if (rr<rc) 
  {
    fcut = exp(-(rc/rr-1.0)*(rc/rr-1.0));
  }
  else
  {
    fcut = 1.0;
  }
  v022 = exp(alpha - beta*rr - gamma*rr*rr)-fcut*(c6/pow(rr,6) + c8/pow(rr,8) + c10/pow(rr,10));
  return v022;
}


// ***************************************DERIVATIVE   d(v022) / dR ************DERIVATIVE
double Dv022DR(double rr)
{
  double Dv022DR;
  double alpha=-3.428227, beta=0.715011, gamma=0.100120, rc=8.422755;
  double c6=0.288396, c8=8.242595, c10=210.984075, fcut,dfcut;

   if(rr<rc)
   {
     fcut = exp(-((rc/rr-1.0)*(rc/rr-1.0)));
     dfcut = 2*exp(-((rc/rr-1.0)*(rc/rr-1.0))) *(rc/rr-1)*(rc/(rr*rr));
   Dv022DR = exp(alpha - beta*rr-gamma*rr*rr )*(-beta -2*gamma*rr) - \
             dfcut*(c6/(pow(rr,6)) + c8/(pow(rr,8)) + c10/(pow(rr,10))) +\
             fcut*(6*c6/(pow(rr,7)) + 8*c8/(pow(rr,9)) + 10*c10/(pow(rr,11)));
   }
   else
   {
         fcut = 1.0;
      Dv022DR = exp(alpha - beta*rr-gamma*rr*rr )*(-beta -2*gamma*rr) +\
                fcut*(6*c6/(pow(rr,7)) + 8*c8/(pow(rr,9)) + 10*c10/(pow(rr,11)));
   }
   
   return Dv022DR;
}

//      end
//***********************************************************************
double v224(double rr)
//     This is the so-called quadrupole-quadrupole interaction term.
//***********************************************************************
{
  double rr5, v224, epsi=0.135269;
  v224 = epsi/pow(rr,5);
  return v224;
}

// ***************************************DERIVATIVE   d(v224) / dR ************DERIVATIVE
double Dv224DR(double rr)
{
  double Dv224DR, epsi=0.135269;
  Dv224DR= -( 5*epsi)/pow(rr,6);
  return Dv224DR;
}

double Dg202Dth(double theta)
{
  double Dg202th;
  Dg202th = 15 * cos(theta);
  return Dg202th;
}
double Dg202Dcosth(double theta)
{
  double Dg202th;
  Dg202th = 15 * cos(theta);
  return Dg202th;
}

double Dg224Dcosth(double theta1,double theta2,double phi,int tag)
{
  double Dg224Dcosth;

  Dg224Dcosth = fact * ( 36*cos(theta1)*cos(theta2)*cos(theta2) -12*cos(theta1) \
                         -16* sin(theta1) *sin(theta2)*cos(theta2)*cos(phi)       );
  return Dg224Dcosth;
}

double Dg224Dsinth(double theta1,double theta2,double phi,int tag)
{
  double Dg224Dsinth, temp1,temp2;
  Dg224Dsinth = fact * (-16* cos(theta1) *sin(theta2)*cos(theta2)*cos(phi) +\
                        2*sin(theta1)*sin(theta2)*sin(theta2)*cos(2*phi)   );
  return Dg224Dsinth;
}

double Dg224Dphi(double theta1,double theta2,double phi)
{
  double Dg224Dphi;
  Dg224Dphi = fact * ( -16 * sin(theta1)*cos(theta1)*sin(theta2)*cos(theta2) + \
                      sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*cos(2*phi)   );
  return Dg224Dphi;
}
double Dg224Dphit(double theta1,double theta2,double phi)
{
  double Dg224Dphit;
  Dg224Dphit = fact * ( -16 * sin(theta1)*cos(theta1)*sin(theta2)*cos(theta2) + \
                     4*sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*cos(phi)   );
  return Dg224Dphit;
}

double Dg224D2phi(double theta1,double theta2,double phi)
{
  double Dg2242Dphi;
  Dg2242Dphi = fact * ( -16 * sin(theta1)*cos(theta1)*sin(theta2)*cos(theta2)*cos(phi) + \
                      sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)   );
  return Dg2242Dphi;
}


double dRdX1(double r1,double r2,double r3,double r4,double rr)
{
  double ddrdx;
  ddrdx = (r1 + r2 - r3 -r4)/(4*rr);
  return ddrdx;
}

// ###############################################    V_H2H2   ################################

double H2H2pot2(double rr, double theta1, double theta2, double sin1, double sin2, double phi12)
{
  double pot;
  pot = v000(rr)+ v022(rr) * g2022(theta2)+ v022(rr) * g2022(theta1) + \
        v224(rr) * g2242(theta1, theta2, sin1,sin2,phi12);
  return pot;

}


double H2H2pot(double rr, double theta1, double theta2, double phi12)
{
  double pot;
  pot = v000(rr) + v022(rr) *2.5*(3.0*(cos(theta2))*(cos(theta2))-1.0) + \
                   v022(rr) *2.5*(3.0*(cos(theta1))*(cos(theta1))-1.0) + \
        v224(rr) * fact * ( 2.0*(3.0*(cos(theta1))*(cos(theta1))-1.0) * (3.0*(cos(theta2))*(cos(theta2))-1.0)-  \
       16.0 * sin(theta1)*cos(theta1)* sin(theta2)*cos(theta2)*cos(phi12)  + \
           (sin(theta1))*(sin(theta1))* (sin(theta2))*(sin(theta2))*cos(2.0*phi12));

  return pot;
}

