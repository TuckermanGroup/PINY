#include <math.h>

/*============*/
/* prototypes */

/* helpers */
double length(const double *);
inline double sprod(const double *, const double *);
void vprod(const double *, const double *, double *);
double dRdX(const double, const double, const double, const double, const double);

/* radial potential compotents */
void V000(const double, double *, double *);
void V022(const double, double *, double *);
void V224(const double, double *, double *);

/* angular potential components */
void G202_trig(const double, double *, double *);
void G224_trig(const double, const double, const double, const double, const double,
               double *, double *, double *, double *, double *, double *);


/*======================================================================================*/
void H2_pot_for(const double *r1, const double *r2, const double *r3, const double *r4,
                double *ff1, double *ff2, double *ff3, double *ff4,
                double *E_pot) {
/*======================================================================================*/

    /*=======================*/
    /* variable declarations */

    unsigned int i;

    /* internal coordinates */
    double R;
    double sintheta1, sintheta2;
    double costheta1, costheta2;
    double cosphi;
    double sign;
    double rFH, rHH, ln1, ln2;

    /* transformation to internal coordinates */
    double FHC[3];   /* COM of 1st H2 */
    double HHC[3];   /* COM of 2nd H2 */
    double vFH[3];   /* bond vector of 1st H2 */
    double vHH[3];   /* bond vector of 2nd H2 */
    double vCC[3];   /* COMs relative position */
    double n1[3];
    double n2[3];
    double n3[3];

    /* some helpers */
    double c1, c2, c3, h1, h2, h3, f1, f2, f3;

    /* potential energy and its derivatives in internal coordinates */
    double V000R, dV000dRR, V022R, dV022dRR, V224R, dV224dRR;
    double dVdR, dVdcostheta1, dVdcostheta2, dVdsintheta1, dVdsintheta2, dVdphi;
    double G202theta1, dG202dcostheta1, G202theta2, dG202dcostheta2;
    double G224, dG224dcostheta1, dG224dcostheta2, dG224dsintheta1, dG224dsintheta2, dG224dcosphi;

    /* derivative transformation */
    double dRdX1, dcostheta1dX1, dcostheta2dX1, dsintheta1dX1, dsintheta2dX1, dphidX1;
    double dRdX2, dcostheta1dX2, dcostheta2dX2, dsintheta1dX2, dsintheta2dX2, dphidX2;
    double dRdX3, dcostheta1dX3, dcostheta2dX3, dsintheta1dX3, dsintheta2dX3, dphidX3;
    double dRdX4, dcostheta1dX4, dcostheta2dX4, dsintheta1dX4, dsintheta2dX4, dphidX4;
    double dRdY1, dcostheta1dY1, dcostheta2dY1, dsintheta1dY1, dsintheta2dY1, dphidY1;
    double dRdY2, dcostheta1dY2, dcostheta2dY2, dsintheta1dY2, dsintheta2dY2, dphidY2;
    double dRdY3, dcostheta1dY3, dcostheta2dY3, dsintheta1dY3, dsintheta2dY3, dphidY3;
    double dRdY4, dcostheta1dY4, dcostheta2dY4, dsintheta1dY4, dsintheta2dY4, dphidY4;
    double dRdZ1, dcostheta1dZ1, dcostheta2dZ1, dsintheta1dZ1, dsintheta2dZ1, dphidZ1;
    double dRdZ2, dcostheta1dZ2, dcostheta2dZ2, dsintheta1dZ2, dsintheta2dZ2, dphidZ2;
    double dRdZ3, dcostheta1dZ3, dcostheta2dZ3, dsintheta1dZ3, dsintheta2dZ3, dphidZ3;
    double dRdZ4, dcostheta1dZ4, dcostheta2dZ4, dsintheta1dZ4, dsintheta2dZ4, dphidZ4;


    /*======================*/
    /* internal coordinates */

    /* [r1, r2, r3, r4] -> [R, theta1, theta2, phi] */

    /* vectors */
    for(i=0; i<3; ++i) {
        FHC[i] = 0.5 * (r1[i] + r2[i]);
        HHC[i] = 0.5 * (r3[i] + r4[i]);
        vFH[i] = r2[i] - r1[i];
        vHH[i] = r4[i] - r3[i];
        vCC[i] = HHC[i] - FHC[i];
    }

    /* normal vectors */
    vprod(vCC, vHH, n1);
    vprod(vCC, vFH, n2);
    vprod(n2, n1, n3);
    sign = sprod(n3, vCC);

    /* vector lengths */
    R = length(vCC);
    rFH = length(vFH);
    rHH = length(vHH);
    ln1 = length(n1);
    ln2 = length(n2);

    /* angles */
    sintheta1 = ln1 / (rHH * R);
    sintheta2 = ln2 / (rFH * R);
    costheta1 = sprod(vHH, vCC) / (rHH * R);
    costheta2 = -sprod(vFH, vCC) / (rFH * R);
    cosphi = sprod(n1, n2) / (ln1 * ln2);
    if (cosphi < -1.0)
        cosphi = -1.0;
    if (cosphi > 1.0)
        cosphi = 1.0;
    if (sign == 0) {
        cosphi = 0;
    }


    /*==============================================================*/
    /* potential energy and its derivatives in internal coordinates */

    V000(R, &V000R, &dV000dRR);
    V022(R, &V022R, &dV022dRR);
    V224(R, &V224R, &dV224dRR);
    G202_trig(costheta1, &G202theta1, &dG202dcostheta1);
    G202_trig(costheta2, &G202theta2, &dG202dcostheta2);
    G224_trig(costheta1, costheta2, sintheta1, sintheta2, cosphi,
              &G224, &dG224dcostheta1, &dG224dcostheta2,
              &dG224dsintheta1, &dG224dsintheta2, &dG224dcosphi);

    *E_pot = V000R + V022R * (G202theta1 + G202theta2) + V224R * G224;

    dVdR = dV000dRR + dV022dRR * (G202theta1 + G202theta2) + dV224dRR * G224;

    dVdcostheta1 = V022R * dG202dcostheta1 + V224R * dG224dcostheta1;
    dVdcostheta2 = -(V022R * dG202dcostheta2 + V224R * dG224dcostheta2);

    dVdsintheta1 = V224R * dG224dsintheta1;
    dVdsintheta2 = V224R * dG224dsintheta2;

    dVdphi = V224R * dG224dcosphi;


    /*=====================================================================*/
    /* coeeficients for derivative transformation to Cartesian coordinates */

    /* some helpers */
    c1 = (r3[0] + r4[0] - r2[0] - r1[0]) / 2;
    c2 = (r3[1] + r4[1] - r2[1] - r1[1]) / 2;
    c3 = (r3[2] + r4[2] - r2[2] - r1[2]) / 2;
    h1 = r4[0] - r3[0];
    h2 = r4[1] - r3[1];
    h3 = r4[2] - r3[2];
    f1 = r2[0] - r1[0];
    f2 = r2[1] - r1[1];
    f3 = r2[2] - r1[2];

    /* X1 */
    dRdX1 = dRdX(r1[0], r2[0], r3[0], r4[0], R);
    dcostheta1dX1 = ((0.5 * (r3[0] - r4[0]) * rHH * R - rHH * sprod(vHH, vCC) * dRdX1) /
                     (rHH * R * rHH * R));
    dcostheta2dX1 = (((-0.5*r3[0] - 0.5*r4[0] + r1[0]) * rFH * R -
                      ((sprod(vFH, vCC)) * (((-r2[0] + r1[0]) * R / rFH) + rFH * dRdX1))) /
                     (rFH * rFH * R * R));
    dsintheta1dX1 = (((((h1*c3 - h3*c1) * h3 + (h1*c2 - h2*c1) * h2) * rHH * R / (2*ln1)) - ln1*rHH*dRdX1) /
                     (rHH*R*rHH*R));
    dsintheta2dX1 = ((((((-f3*c1 + f1*c3) * (0.5*f3 - c3)) + ((-f1*c2 + f2*c1) * (-0.5*f2 + c2))) *
                       rFH * R / ln2) - ln2 * (((-r2[0] + r1[0]) * R / rFH) + rFH * dRdX1)) /
                      (rFH*rFH*R*R));
    dphidX1 = (((-c3*c3*h1 + 0.5*c3*h1*f3 + h3*c3*(0.5*f1 + c1) - c1*h3*f3 - h2*f2*c1 +
                 h2*c2*(0.5*f1 + c1) + 0.5*c2*h1*f2 - c2*c2*h1) * ln1 * ln2 -
                sprod(n1,n2) * (((c3*h1 - c1*h3)*h3 - (c1*h2-c2*h1)*h2) * ln2 / (2*ln1) +
                (((c3*f1 - c1*f3)*(-c3 + 0.5*f3) + (c1*f2 - c2*f1)*(-0.5*f2 + c2)) * ln1 / ln2))) /
               (ln1*ln1*ln2*ln2));

    /* X2 */
    dRdX2 = dRdX(r1[0], r2[0], r3[0], r4[0], R);
    dcostheta1dX2 = ((0.5 * (r3[0] - r4[0]) * rHH * R - rHH * sprod(vHH, vCC) * dRdX2) /
                     (rHH * R * rHH * R));
    dcostheta2dX2 = ((((r3[0] / 2 + r4[0] / 2 - r2[0]) * rFH * R -
                      ((sprod(vFH, vCC)) * (((r2[0] - r1[0]) * R / rFH) + rFH * dRdX2))) /
                     (rFH * rFH * R * R)));
    dsintheta1dX2 = (((((h1*c3 - h3*c1) * h3 + (h1*c2 - h2*c1) * h2) * rHH * R / (2*ln1)) - ln1*rHH*dRdX2) /
                     (rHH*R*rHH*R));
    dsintheta2dX2 = ((((((-f3*c1 + f1*c3) * (0.5*f3 + c3)) + ((-f1*c2 + f2*c1) * (-0.5*f2 - c2))) *
                       rFH * R / ln2) - ln2 * (((r2[0] - r1[0]) * R / rFH) + rFH * dRdX2)) /
                      (rFH*rFH*R*R));
    dphidX2 = (((c3*c3*h1 + 0.5*c3*h1*f3 - h3*c3*(-0.5*f1 + c1) - c1*h3*f3 - h2*f2*c1 -
                 h2*c2*(-0.5*f1 + c1) + 0.5*c2*h1*f2 + c2*c2*h1) * ln1 * ln2 -
                sprod(n1, n2) * (((c3*h1 - c1*h3)*h3 - (c1*h2 - c2*h1)*h2) * ln2 / (2*ln1) +
                (((c3*f1 - c1*f3) * (c3 + 0.5*f3) + (c1*f2 - c2*f1) * (-0.5*f2 - c2)) * ln1 / ln2))) /
               (ln1*ln1*ln2*ln2));

    /* X3 */
    dRdX3 = dRdX(r3[0], r4[0], r2[0], r1[0], R);
    dcostheta1dX3 = (((0.5*r1[0] + 0.5*r2[0] - r3[0]) * rHH * R - sprod(vHH, vCC) * 
                      ((r3[0] - r4[0]) / rHH * R + rHH*dRdX3)) / (rHH*R*rHH*R));
    dcostheta2dX3 = ((((r2[0] / 2 - r1[0] / 2) * rFH * R - ((sprod(vFH, vCC)) * rFH * dRdX3)) /
                      (rFH*rFH*R*R)));
    dsintheta1dX3 = (((((c3*h1 - c1*h3) * (-c3 - 0.5*h3) + (c1*h2 - c2*h1) * (0.5*h2 + c2)) * rHH * R / ln1) -
                      (ln1 * ((r3[0] - r4[0]) * R / rHH + rHH*dRdX3))) / (rHH*R*rHH*R));
    dsintheta2dX3 = ((((((c3*f1 - c1*f3) * (-f3)) + ((c1*f2 - c2*f1) * f2)) * rFH * R / (2*ln2)) -
                      ln2*rFH*dRdX3) / (rFH*rFH*R*R));
    dphidX3 = (((-c3*c3*f1 - c3*f3*(-c1 + 0.5*h1) - 0.5*h3*c3*f1 + c1*h3*f3 + c1*h2*f2 - 0.5*h2*c2*f1 -
                 c2*f2*(0.5*h1 - c1) - c2*c2*f1) * ln1 * ln2 -
                sprod(n1, n2) * (((c3*h1 - c1*h3) * (-c3 - 0.5*h3) +
                (c1*h2 - c2*h1) * (0.5*h2 + c2)) * ln2 / ln1 +
                ((((c3*f1 - c1*f3) * (-f3)) + ((c1*f2 - c2*f1) * f2)) * ln1 / (2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    /* X4 */
    dRdX4 = dRdX(r3[0], r4[0], r2[0], r1[0], R);
    dcostheta1dX4 = (((-0.5*r1[0] - 0.5*r2[0] + r4[0]) * rHH * R - sprod(vHH, vCC) *
                      ((r4[0] - r3[0]) / rHH * R + rHH * dRdX4)) / (rHH*R*rHH*R));
    dcostheta2dX4 = ((((0.5*r2[0] - 0.5*r1[0]) * rFH * R - ((sprod(vFH, vCC)) * rFH * dRdX4)) /
                      (rFH*rFH*R*R)));
    dsintheta1dX4 = (((((c3*h1 - c1*h3) * (c3 - 0.5*h3) + (c1*h2 - c2*h1) * (0.5*h2 - c2)) * rHH * R / ln1) -
                      (ln1 * ((r4[0] - r3[0]) * R / rHH + rHH * dRdX4))) / (rHH*R*rHH*R));
    dsintheta2dX4 = ((((((c3*f1 - c1*f3) * (-f3)) + ((c1*f2 - c2*f1) * f2)) * rFH * R / (2*ln2)) -
                      ln2*rFH*dRdX4) / (rFH*rFH*R*R));
    dphidX4 = (((c3*c3*f1 - c3*f3*(c1 + 0.5*h1) - 0.5*h3*c3*f1 + c1*h3*f3 + c1*h2*f2 - 0.5*h2*c2*f1 -
                 c2*f2*(0.5*h1 + c1) + c2*c2*f1) * ln1 * ln2 -
                sprod(n1, n2) * (((c3*h1 - c1*h3) * (c3 - 0.5*h3) +
                (c1*h2 - c2*h1) * (0.5*h2 - c2)) * ln2 / ln1 +
                ((((c3*f1 - c1*f3) * (-f3)) + ((c1*f2 - c2*f1) * f2)) * ln1 / (2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    /* Y1 */
    dRdY1 = dRdX(r1[1], r2[1], r3[1], r4[1], R);
    dcostheta1dY1 = ((0.5*(r3[1] - r4[1])*rHH*R - rHH*sprod(vHH, vCC)*dRdY1) / (rHH*R*rHH*R));
    dcostheta2dY1 = (((-0.5*r3[1] - 0.5*r4[1] + r1[1]) * rFH * R -
                      ((sprod(vFH, vCC)) * (((-r2[1] + r1[1]) * R / rFH) + rFH*dRdY1))) /
                     (rFH*rFH*R*R));
    dsintheta1dY1 = (((((c3*h2 - h3*c2)*h3 + (h2*c1 - c2*h1)*h1) * rHH * R / (2*ln1)) - (ln1*rHH*dRdY1)) /
                     (rHH*R*rHH*R));
    dsintheta2dY1 = (((((((f3*c2 - c3*f2) * (-0.5*f3 + c3)) +
                         ((-f1*c2 + f2*c1)*(0.5*f1 - c1))) * rFH * R / ln2) -
                        ln2 * (((-r2[1] + r1[1]) * R / rFH) + rFH*dRdY1)) / (rFH*rFH*R*R)));
    dphidY1 = (((-c2*h3*f3 + h3*c3*(0.5*f2 + c2) + 0.5*c3*h2*f3 - c3*c3*h2 - c1*c1*h2 +
                 0.5*c1*h2*f1 + h1*c1*(0.5*f2 + c2) - c2*h1*f1) * ln1 * ln2 -
                sprod(n1, n2) * ((((c3*h2 - h3*c2)*h3 + (h2*c1 - c2*h1)*h1))*ln2 / (2*ln1) +
                ((((f3*c2 - c3*f2)*(-0.5*f3 + c3)) + ((-f1*c2 + f2*c1)*(0.5*f1 - c1))) * ln1 / ln2))) /
               (ln1*ln1*ln2*ln2));

    /* Y2 */
    dRdY2 = dRdX(r1[1], r2[1], r3[1], r4[1], R);
    dcostheta1dY2 = ((0.5*(r3[1] - r4[1])*rHH*R - rHH*sprod(vHH, vCC)*dRdY2) / (rHH*R*rHH*R));
    dcostheta2dY2 = (((0.5*r3[1] + 0.5*r4[1] - r2[1]) * rFH * R - ((sprod(vFH, vCC)) *
                      (((r2[1] - r1[1]) * R / rFH) + rFH*dRdY2))) / (rFH*rFH*R*R));
    dsintheta1dY2 = (((((c3*h2 - h3*c2)*h3 + (h2*c1 - c2*h1)*h1) * rHH * R / (2*ln1)) -
                      (ln1*rHH*dRdY2)) / (rHH*R*rHH*R));
    dsintheta2dY2 = ((((((f3*c2 - c3*f2) * (-0.5*f3 - c3)) + ((-f1*c2 + f2*c1)*(0.5*f1 + c1))) *
                       rFH * R / ln2) - ln2*(((r2[1]-r1[1])*R/rFH ) + rFH*dRdY2)) /
                     (rFH*rFH*R*R));
    dphidY2 = (((-c2*h3*f3 - h3*c3*(-0.5*f2 + c2) + 0.5*c3*h2*f3 + c3*c3*h2 + c1*c1*h2 +
                 0.5*c1*h2*f1 - h1*c1*(-0.5*f2 + c2) - c2*h1*f1) * ln1 * ln2 -
                sprod(n1, n2) * ((((c3*h2 - h3*c2)*h3 + (h2*c1 - c2*h1)*h1)) * ln2 / (2*ln1) +
                ((((f3*c2 - c3*f2)*(-0.5*f3 - c3)) + ((-f1*c2 + f2*c1)*(0.5*f1 + c1))) * ln1 / ln2))) /
               (ln1*ln1*ln2*ln2));

    /* Y3 */
    dRdY3 = dRdX(r3[1], r4[1], r2[1], r1[1], R);
    dcostheta1dY3 = (((0.5*r1[1] + 0.5*r2[1] - r3[1])*rHH*R - sprod(vHH, vCC) * ((r3[1] - r4[1]) /
                      rHH*R +rHH*dRdY3)) / (rHH*R*rHH*R));
    dcostheta2dY3 = (((0.5*(r2[1] - r1[1])*rFH*R - sprod(vFH, vCC)*rFH*dRdY3) / (rFH*rFH*R*R)));
    dsintheta1dY3 = (((((c2*h3 - c3*h2) * (c3 + 0.5*h3) + (c1*h2 - c2*h1) * (-0.5*h1 - c1)) * rHH*R / ln1) -
                      (ln1*((r3[1] - r4[1]) * R / rHH + rHH*dRdY3))) / (rHH*R*rHH*R));
    dsintheta2dY3 = ((((((c2*f3 - c3*f2)*f3) + ((c1*f2 - c2*f1)*(-f1))) * rFH * R / (2*ln2)) -
                      ln2*rFH*dRdY3) / (rFH*rFH*R*R));
    dphidY3 = (((c2*h3*f3 - 0.5*h3*c3*f2 - c3*f3*(0.5*h2 - c2) - c3*c3*f2 - c1*c1*f2 - c1*f1*(0.5*h2 - c2) -
                 0.5*h1*c1*f2 + c2*h1*f1) * ln1 * ln2 -
                sprod(n1, n2) * (((c2*h3 - c3*h2)*(c3 + 0.5*h3) + (c1*h2 - c2*h1)*(-0.5*h1 - c1)) * ln2/ln1 + 
                ((((c2*f3 - c3*f2)*f3) + ((c1*f2 - c2*f1)*(-f1))) * ln1 / (2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    /* Y4 */
    dRdY4 = dRdX(r3[1], r4[1], r2[1], r1[1], R);
    dcostheta1dY4 = (((r4[1] - 0.5*r2[1] - 0.5*r1[1])*rHH*R - sprod(vHH, vCC)*((r4[1] - r3[1]) /
                      rHH*R + rHH*dRdY4)) / (rHH*R*rHH*R));
    dcostheta2dY4 = ((((0.5*r2[1] - 0.5*r1[1])*rFH*R - sprod(vFH, vCC)*rFH*dRdY4) /
                      (rFH*rFH*R*R)));
    dsintheta1dY4 = (((((c2*h3 - c3*h2)*(-c3 + 0.5*h3) + (c1*h2 - c2*h1)*(-0.5*h1 + c1)) * rHH * R/ ln1) -
                      (ln1 * ((r4[1] - r3[1]) * R / rHH + rHH*dRdY4))) / (rHH*R*rHH*R));
    dsintheta2dY4 = ((((((c2*f3-c3*f2)*(f3)) + ((c1*f2-c2*f1)*(-f1)))*rFH*R/(2*ln2)) -
                      ln2*(rFH*dRdY4))/(rFH*rFH*R*R));
    dphidY4 = (((c2*h3*f3 - 0.5*h3*c3*f2 - c3*f3*(0.5*h2 + c2) + c3*c3*f2 + c1*c1*f2 - c1*f1*(0.5*h2 + c2) -
                 0.5*h1*c1*f2 + c2*h1*f1) * ln1 * ln2 -
                sprod(n1, n2) * (((c2*h3 - c3*h2) * (-c3 + 0.5*h3) +
                (c1*h2 - c2*h1) * (-0.5*h1 + c1)) * ln2 / ln1 +
                ((((c2*f3 - c3*f2)*f3) + (c1*f2 - c2*f1)*(-f1)) * ln1 / (2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    /* Z1 */
    dRdZ1 = dRdX(r1[2], r2[2], r3[2], r4[2], R);
    dcostheta1dZ1 = ((0.5*(r3[2] - r4[2])*rHH*R - rHH*sprod(vHH, vCC)*dRdZ1) / (rHH*R*rHH*R));
    dcostheta2dZ1 = (((-0.5*r3[2] - 0.5*r4[2] + r1[2])*rFH*R - ((sprod(vFH,vCC)) *
                      ((-r2[2] + r1[2]) * R / rFH + rFH*dRdZ1))) / (rFH*rFH*R*R));
    dsintheta1dZ1 = (((((h2*c3 - h3*c2)*(-h2) + (h3*c1 - h1*c3)*h1) * rHH * R / (2*ln1)) -
                      (ln1*rHH*dRdZ1)) / (rHH*R*rHH*R));
    dsintheta2dZ1 = ((((((f3*c2 - c3*f2)*(0.5*f2 - c2)) + ((f1*c3 - f3*c1)*(-0.5*f1 + c1)))*rFH*R/ln2) -
                      ln2 * (((-r2[2] + r1[2])*R/rFH) + rFH*dRdZ1)) / (rFH*rFH*R*R));
    dphidZ1 = (((-c2*h3*c2 + 0.5*h3*c2*f2 + h2*c2*(c3 + 0.5*f3) - c3*h2*f2 -
                 c3*h1*f1 + h1*c1*(c3 + 0.5*f3) + 0.5*c1*h3*f1 -c1*c1*h3) * ln1 * ln2 -
                sprod(n1, n2) * (((((h2*c3 - h3*c2)*(-h2) + (h3*c1 - h1*c3)*h1))) * ln2/(2*ln1) +
                (((f3*c2 - c3*f2)*(0.5*f2 - c2)) + ((f1*c3 - f3*c1)*(-0.5*f1 + c1))) * ln1 / ln2)) /
               (ln1*ln1*ln2*ln2));

    /* Z2 */
    dRdZ2 = dRdX(r1[2], r2[2], r3[2], r4[2], R);
    dcostheta1dZ2 = (((0.5*r3[2] - 0.5*r4[2])*rHH*R - rHH*sprod(vHH, vCC)*dRdZ2) / (rHH*R*rHH*R));
    dcostheta2dZ2 = (((0.5*r3[2] + 0.5*r4[2] - r2[2])*rFH*R - ((sprod(vFH,vCC)) *
                      ((r2[2] - r1[2]) * R / rFH + rFH*dRdZ2))) / (rFH*rFH*R*R));
    dsintheta1dZ2 = (((((h2*c3 - h3*c2)*(-h2) + (h3*c1 - h1*c3)*h1) * rHH * R / (2*ln1)) - (ln1*rHH*dRdZ2)) /
                     (rHH*R*rHH*R));
    dsintheta2dZ2 = ((((((f3*c2 - c3*f2)*(0.5*f2 + c2)) + ((f1*c3 - f3*c1)*(-0.5*f1 - c1))) * rFH*R/ln2) -
                      ln2 * ((r2[2] - r1[2]) * R / rFH + rFH*dRdZ2)) / (rFH*rFH*R*R));
    dphidZ2 = (((c2*h3*c2 + 0.5*h3*c2*f2 - h2*c2*(c3 - 0.5*f3) - c3*h2*f2 -
                 c3*h1*f1 - h1*c1*(c3 - 0.5*f3) + 0.5*c1*h3*f1 + c1*c1*h3) * ln1 * ln2 -
                sprod(n1, n2) * (((((h2*c3 - h3*c2)*(-h2) + (h3*c1 - h1*c3)*h1))) * ln2 / (2*ln1) +
                (((f3*c2 - c3*f2)*(0.5*f2 + c2)) + ((f1*c3 - f3*c1)*(-0.5*f1 - c1))) * ln1 / ln2)) /
               (ln1*ln1*ln2*ln2));

    /* Z3 */
    dRdZ3 = dRdX(r3[2], r4[2], r2[2], r1[2], R);
    dcostheta1dZ3 = (((0.5*r1[2] + 0.5*r2[2] - r3[2])*rHH*R -
                      sprod(vHH,vCC) * ((r3[2] - r4[2]) / rHH * R + rHH*dRdZ3)) /
                     (rHH*R*rHH*R));
    dcostheta2dZ3 = ((0.5*(r2[2] - r1[2])*rFH*R - ((sprod(vFH,vCC))*rFH*dRdZ3)) /
                     (rFH*rFH*R*R));
    dsintheta1dZ3 = (((((c2*h3 - c3*h2)*(-c2 - 0.5*h2) + (c3*h1 - c1*h3)*(0.5*h1 + c1))*rHH*R/ln1) -
                    (ln1 * ((r3[2]-r4[2])*R/rHH + rHH*dRdZ3))) / (rHH*R*rHH*R));
    dsintheta2dZ3 = ((((((c2*f3 - c3*f2)*(-f2)) + ((c3*f1 - c1*f3)*f1)) * rFH*R / (2*ln2)) -
                    ln2 * (rFH*dRdZ3)) / (rFH*rFH*R*R));
    dphidZ3 = (((-c2*c2*f3 - c2*f2*(-c3 + 0.5*h3) - 0.5*h2*c2*f3 + c3*h2*f2 + c3*h1*f1 - 0.5*h1*c1*f3 -
                 c1*f1*(0.5*h3 - c3) - c1*c1*f3) * ln1 * ln2 -
                sprod(n1, n2) * (((c2*h3 - c3*h2)*(-c2 - 0.5*h2) + (c3*h1 - c1*h3)*(0.5*h1 + c1))*ln2/ln1 +
                ((((c2*f3 - c3*f2)*(-f2)) + ((c3*f1 - c1*f3)*(f1)))*ln1/(2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    /* Z4 */
    dRdZ4 = dRdX(r3[2], r4[2], r2[2], r1[2], R);
    dcostheta1dZ4 = (((-0.5*r1[2] - 0.5*r2[2] + r4[2])*rHH*R -
                      sprod(vHH, vCC)*((r4[2] - r3[2])/rHH*R + rHH*dRdZ4)) / (rHH*R*rHH*R));
    dcostheta2dZ4 = (((0.5*(r2[2] - r1[2])*rFH*R -
                      ((sprod(vFH, vCC))*rFH*dRdZ4)) / (rFH*rFH*R*R)));
    dsintheta1dZ4 = (((((c2*h3 - c3*h2)*(+c2 - 0.5*h2)+(c3*h1 - c1*h3)*(0.5*h1 - c1))*rHH*R/ln1) -
                    (ln1 * ((r4[2] - r3[2])*R/rHH + rHH*dRdZ4))) / (rHH*R*rHH*R));
    dsintheta2dZ4 = ((((((c2*f3 - c3*f2)*(-f2)) + ((c3*f1 - c1*f3)*f1))*rFH*R/(2*ln2)) -
                    ln2*rFH*dRdZ4) / (rFH*rFH*R*R));
    dphidZ4 = (((c2*c2*f3 - c2*f2*(c3 + 0.5*h3) - 0.5*h2*c2*f3 + c3*h2*f2 + c3*h1*f1 - 0.5*h1*c1*f3 -
                 c1*f1*(0.5*h3+c3) + c1*c1*f3) * ln1 * ln2 -
                sprod(n1, n2) * (((c2*h3 - c3*h2)*(c2 - 0.5*h2) + (c3*h1 - c1*h3)*(0.5*h1 - c1)) * ln2 / ln1 +
                ((((c2*f3 - c3*f2)*(-f2)) + ((c3*f1 - c1*f3)*f1)) * ln1 / (2*ln2)))) /
               (ln1*ln1*ln2*ln2));

    if ((sign == 0) || (ln1 == 0) || (ln2 == 0)) {
        dphidX1 = dphidX2 = dphidX3 = dphidX4 = 0;
        dphidY1 = dphidY2 = dphidY3 = dphidY4 = 0;
        dphidZ1 = dphidZ2 = dphidZ3 = dphidZ4 = 0;
    }

    if (ln1 == 0) {
        dsintheta1dX1 = dsintheta1dX2 = dsintheta1dX3 = dsintheta1dX4 = 0;
        dsintheta1dY1 = dsintheta1dY2 = dsintheta1dY3 = dsintheta1dY4 = 0;
        dsintheta1dZ1 = dsintheta1dZ2 = dsintheta1dZ3 = dsintheta1dZ4 = 0;
    }

    if(ln2 == 0) {
        dsintheta2dX1 = dsintheta2dX2 = dsintheta2dX3 = dsintheta2dX4 = 0;
        dsintheta2dY1 = dsintheta2dY2 = dsintheta2dY3 = dsintheta2dY4 = 0;
        dsintheta2dZ1 = dsintheta2dZ2 = dsintheta2dZ3 = dsintheta2dZ4 = 0;
    }


    /*==========================================================*/
    /* compose Cartesian forces and store them in output arrays */

    ff1[0] = -(dRdX1*dVdR +
               dcostheta1dX1*dVdcostheta1 + dcostheta2dX1*dVdcostheta2 +
               dsintheta1dX1*dVdsintheta1 + dsintheta2dX1*dVdsintheta2 + dphidX1*dVdphi);

    ff2[0] = -(dRdX2*dVdR + 
               dcostheta1dX2*dVdcostheta1 + dcostheta2dX2*dVdcostheta2 +
               dsintheta1dX2*dVdsintheta1 + dsintheta2dX2*dVdsintheta2 + dphidX2*dVdphi);

    ff3[0] = -(dRdX3*dVdR +
               dcostheta1dX3*dVdcostheta1 + dcostheta2dX3*dVdcostheta2 +
               dsintheta1dX3*dVdsintheta1 + dsintheta2dX3*dVdsintheta2 + dphidX3*dVdphi);

    ff4[0] = -(dRdX4*dVdR +
               dcostheta1dX4*dVdcostheta1 + dcostheta2dX4*dVdcostheta2 +
               dsintheta1dX4*dVdsintheta1 + dsintheta2dX4*dVdsintheta2 + dphidX4*dVdphi);

    ff1[1] = -(dRdY1*dVdR +
               dcostheta1dY1*dVdcostheta1 + dcostheta2dY1*dVdcostheta2 +
               dsintheta1dY1*dVdsintheta1 + dsintheta2dY1*dVdsintheta2 + dphidY1*dVdphi);

    ff2[1] = -(dRdY2*dVdR +
               dcostheta1dY2*dVdcostheta1 + dcostheta2dY2*dVdcostheta2 +
               dsintheta1dY2*dVdsintheta1 + dsintheta2dY2*dVdsintheta2 + dphidY2*dVdphi);

    ff3[1] = -(dRdY3*dVdR +
               dcostheta1dY3*dVdcostheta1 + dcostheta2dY3*dVdcostheta2 +
               dsintheta1dY3*dVdsintheta1 + dsintheta2dY3*dVdsintheta2 + dphidY3*dVdphi);

    ff4[1] = -(dRdY4*dVdR +
               dcostheta1dY4*dVdcostheta1 + dcostheta2dY4*dVdcostheta2 +
               dsintheta1dY4*dVdsintheta1 + dsintheta2dY4*dVdsintheta2 + dphidY4*dVdphi);

    ff1[2] = -(dRdZ1*dVdR +
               dcostheta1dZ1*dVdcostheta1 + dcostheta2dZ1*dVdcostheta2 +
               dsintheta1dZ1*dVdsintheta1 + dsintheta2dZ1*dVdsintheta2 + dphidZ1*dVdphi);

    ff2[2] = -(dRdZ2*dVdR +
               dcostheta1dZ2*dVdcostheta1 + dcostheta2dZ2*dVdcostheta2 +
               dsintheta1dZ2*dVdsintheta1 + dsintheta2dZ2*dVdsintheta2 + dphidZ2*dVdphi);

    ff3[2] = -(dRdZ3*dVdR +
               dcostheta1dZ3*dVdcostheta1 + dcostheta2dZ3*dVdcostheta2 +
               dsintheta1dZ3*dVdsintheta1 + dsintheta2dZ3*dVdsintheta2 + dphidZ3*dVdphi);

    ff4[2] = -(dRdZ4*dVdR +
               dcostheta1dZ4*dVdcostheta1 + dcostheta2dZ4*dVdcostheta2 +
               dsintheta1dZ4*dVdsintheta1 + dsintheta2dZ4*dVdsintheta2 + dphidZ4*dVdphi);
}


double length(const double *vec) {
    /* returns the length of a 3D vector */

    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}


inline double sprod(const double *vec1, const double *vec2) {
    /* returns the scalar product of two 3D vectors */

    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}


void vprod(const double *vec1, const double *vec2, double *vec3) {
    /* computes the vector product of two 3D vectors */

    vec3[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec3[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec3[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


double dRdX(const double r1, const double r2, const double r3, const double r4, const double R) {

    return (r1 + r2 - r3 - r4) / (4 * R);
}


void V000(const double R, double *V000, double *dV000dR) {
    /* Computes the radial component V000 and its derivative.
     * This is the isotropic contribution.
     */

    /* potential parameters */
    const double alpha = 0.655914;
    const double beta = 1.018447;
    const double gamma = 8.070618E-02;
    const double rc = 9.034308;
    const double c6 = 13.076837;
    const double c8 = 80.700360;
    const double c10 = 3687.082967;

    const double R2 = R * R;
    const double R6 = R2 * R2 * R2;
    const double R8 = R6 * R2;
    const double Rinv = 1 / R;
    const double cRinv6 = c6 / R6;
    const double cRinv8 = c8 / R8;
    const double cRinv10 = c10 / (R8 * R2);
    const double b = exp(alpha - beta*R - gamma*R2);

    double a, fcut, dfcut;

    if (R < rc) {
        a = (rc/R - 1.0);
        fcut = exp(-a*a);
        dfcut = 2 * fcut * a * (rc/R2);
        *dV000dR = b * (-beta - 2*gamma*R) -
                   dfcut * (cRinv6 + cRinv8 + cRinv10) +
                   fcut * (6*cRinv6*Rinv + 8*cRinv8*Rinv + 10*cRinv10*Rinv);
    } else {
        fcut = 1.0;
        *dV000dR = b * (-beta -2*gamma*R) +
                   6*cRinv6*Rinv + 8*cRinv8*Rinv + 10*cRinv10*Rinv;
    }

    *V000 = b - fcut * (cRinv6 + cRinv8 + cRinv10);
}


void V022(const double R, double* V022, double *dV022dR) {
    /* Computes the radial components, V022 and V202, of the H2-H2 potential.
     * V022 = V202 due to symmetry.
     */

    /* potential parameters */
    const double alpha = -3.428227;
    const double beta = 0.715011;
    const double gamma = 0.100120;
    const double rc = 8.422755;
    const double c6 = 0.288396;
    const double c8 = 8.242595;
    const double c10 = 210.984075;

    const double R2 = R * R;
    const double R6 = R2 * R2 * R2;
    const double R8 = R6 * R2;
    const double Rinv = 1 / R;
    const double cRinv6 = c6 / R6;
    const double cRinv8 = c8 / R8;
    const double cRinv10 = c10 / (R8 * R2);
    const double b = exp(alpha - beta*R - gamma*R2);

    double a, fcut, dfcut;

    if (R < rc) {
        a = (rc/R - 1);
        fcut = exp(-a*a);
        dfcut = 2 * fcut * a * (rc/R2);
        *dV022dR = b * (-beta - 2*gamma*R) -
                   dfcut * (cRinv6 + cRinv8 + cRinv10) +
                   fcut * (6*cRinv6*Rinv + 8*cRinv8*Rinv + 10*cRinv10*Rinv);
    } else {
        fcut = 1;
        *dV022dR = b * (-beta - 2*gamma*R) +
                   6*cRinv6*Rinv + 8*cRinv8*Rinv + 10*cRinv10*Rinv;
    }

    *V022 = b - fcut * (cRinv6 + cRinv8 + cRinv10);
}


void V224(const double R, double *V224, double *dV224dR) {
    /* This is the so-called quadrupole-quadrupole interaction term. */

    /* potential parameter */
    const double epsi = 0.135269;

    const double R2 = R * R;
    const double R5 = R2 * R2 * R;
    const double R6 = R2 * R2 * R2;

    *V224 = epsi/R5;
    *dV224dR = -5.0 * epsi / R6;
}


void G202_trig(const double costheta, double *G202, double *dG202dcostheta) {

    *G202 = 2.5 * (3 * costheta * costheta - 1);
    *dG202dcostheta =  15.0 * costheta;
}


void G224_trig(const double ct1, const double ct2,
               const double st1, const double st2, const double cp,
               double *G224,
               double *dG224dct1, double *dG224dct2,
               double *dG224dst1, double *dG224dst2,
               double *dG224dcp) {

    /* potential parameter */
    const double fact = 45.0 / (4.0 * sqrt(70.0));

    const double ct1sq = ct1 * ct1;
    const double ct2sq = ct2 * ct2;
    const double st1sq = st1 * st1;
    const double st2sq = st2 * st2;
    const double cos2phi = (2*cp*cp - 1);

    *G224 = fact * (2.0 * (3.0*ct1sq - 1.0) * (3.0*ct2sq - 1.0) -
            16.0*st1*ct1*st2*ct2*cp + st1sq*st2sq*cos2phi);

    *dG224dct1 = fact * (36*ct1*ct2sq - 12*ct1 - 16*st1*st2*ct2*cp);

    *dG224dct2 = fact * (36*ct2*ct1sq - 12*ct2 - 16*st2*st1*ct1*cp);

    *dG224dst1 = fact * (-16*ct1*st2*ct2*cp + 2*st1*st2sq*cos2phi);

    *dG224dst2 = fact * (-16*ct2*st1*ct1*cp + 2*st2*st1*st1*cos2phi);

    *dG224dcp = fact * (-16*st1*ct1*st2*ct2 + 4*st1sq*st2sq*cp);
}
