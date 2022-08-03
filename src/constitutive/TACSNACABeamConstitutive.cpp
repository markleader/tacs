#include "TACSNACABeamConstitutive.h"

TACSNACABeamConstitutive::TACSNACABeamConstitutive( TACSMaterialProperties *properties,
                                                    int naca,
                                                    TacsScalar chord_init,
                                                    TacsScalar twist_init,
                                                    TacsScalar wall_init,
                                                    int chord_dv,
                                                    int twist_dv,
                                                    int wall_dv,
                                                    TacsScalar chord_lb,
                                                    TacsScalar chord_ub,
                                                    TacsScalar twist_lb,
                                                    TacsScalar twist_ub,
                                                    TacsScalar wall_lb,
                                                    TacsScalar wall_ub,
                                                    int _npts,
                                                    int _use_cm,
                                                    TacsScalar _xrot,
                                                    TacsScalar _yrot ){
  props = properties;
  props->incref();

  TacsScalar m = (naca/1000)/100.0;
  TacsScalar p = ((naca%1000)/100)/10.0;
  TacsScalar tt = (naca%100)/100.0;

  chord = chord_init;
  twist = twist_init;
  wall = wall_init;
  chordDV = chord_dv;
  twistDV = twist_dv;
  wallDV = wall_dv;
  chordLb = chord_lb;
  chordUb = chord_ub;
  twistLb = twist_lb;
  twistUb = twist_ub;
  wallLb = wall_lb;
  wallUb = wall_ub;
  npts = _npts;
  use_cm = _use_cm;
  xrot = _xrot;
  yrot = _yrot;

  double a1 = 0.2969;
  double a2 = 0.1260;
  double a3 = 0.3516;
  double a4 = 0.2843;
  double a5 = 0.1015;

  TacsScalar yc_pts[npts];
  TacsScalar xl_pts[npts];
  TacsScalar xu_pts[npts];
  TacsScalar yl_pts[npts];
  TacsScalar yu_pts[npts];
}

TACSNACABeamConstitutive::~TACSNACABeamConstitutive(){
  props->decref();
}

int TACSNACABeamConstitutive::getDesignVarNums( int elemIndex,
                                                int dvLen,
                                                int dvNums[] ){
  int index = 0;
  if (chordDV >= 0){
    if (dvNums && dvLen > index){
      dvNums[index] = chordDV;
    }
    index++;
  }
  if (twistDV >= 0){
    if (dvNums && dvLen > index){
      dvNums[index] = twistDV;
    }
    index++;
  }
  if (wallDV >= 0){
    if (dvNums && dvLen > index){
      dvNums[index] = wallDV;
    }
    index++;
  }
  return index;
}


int TACSNACABeamConstitutive::setDesignVars( int elemIndex,
                                             int dvLen,
                                             const TacsScalar dvs[] ){
  int index = 0;
  if (chordDV >= 0){
    chord = dvs[index];
    index++;
  }
  if (twistDV >= 0){
    twist = dvs[index];
    index++;
  }
  if (wallDV >= 0){
    wall = dvs[index];
    index++;
  }
  return index;
}

int TACSNACABeamConstitutive::getDesignVars( int elemIndex,
                                             int dvLen,
                                             TacsScalar dvs[] ){
  int index = 0;
  if (chordDV >= 0){
    dvs[index] = chord;
    index++;
  }
  if (twistDV >= 0){
    dvs[index] = twist;
    index++;
  }
  if (wallDV >= 0){
    dvs[index] = wall;
    index++;
  }
  return index;
}

int TACSNACABeamConstitutive::getDesignVarRange( int elemIndex, int dvLen,
                                                 TacsScalar lb[],
                                                 TacsScalar ub[] ){
  int index = 0;
  if (chordDV >= 0){
    lb[index] = chordLb;
    ub[index] = chordUb;
    index++;
  }
  if (twistDV >= 0){
    lb[index] = twistLb;
    ub[index] = twistUb;
    index++;
  }
  if (wallDV >= 0){
    lb[index] = wallLb;
    ub[index] = wallUb;
    index++;
  }
  return index;
}

void TACSNACABeamConstitutive::computeNACAPoints( TacsScalar yc_pts[],
                                                  TacsScalar xl_pts[],
                                                  TacsScalar xu_pts[],
                                                  TacsScalar yl_pts[],
                                                  TacsScalar yu_pts[] ){

  TacsScalar x = 0.0;
  TacsScalar dx = 1.0/(npts-1);
  TacsScalar yt;
  TacsScalar dyc_dx;
  TacsScalar theta;
  for ( int i = 0; i < npts; i++ ){
    yt = 5.0*tt*(a1*x**0.5 - a2*x - a3*x**2 + a4*x**3 - a5*x**4);
    if ( (m > 0.0) && (p > 0.0) ){
      if ( x <= p ){
        yc_pts[i] = (m/p**2)*(2*p*x - x**2);
        dyc_dx = 2*m*(p - x)/(p**2);
      }
      else {
        yc_pts[i] = m*((1 - 2.0*p) + 2*p*x - x**2)/((1 - p)**2);
        dyc_dx = 2*m*(p - x)/((1 - p)**2);
      }
      theta = atan(dyc_dx);
      xl_pts[i] = x + yt*sin(theta);
      xu_pts[i] = x - yt*sin(theta);
      yl_pts[i] = yc_pts[i] - yt*cos(theta);
      yu_pts[i] = yc_pts[i] + yt*cos(theta);
    }
    else {
      yc_pts[i] = 0.0;
      xl_pts[i] = x;
      xu_pts[i] = x;
      yl_pts[i] = -yt;
      yu_pts[i] = yt;
    }
  x += dx;
  }
}

void TACSNACABeamConstitutive::computeNACAPoint( TacsScalar x,
                                                 TacsScalar yc,
                                                 TacsScalar xl,
                                                 TacsScalar xu,
                                                 TacsScalar yl,
                                                 TacsScalar yu ){

  TacsScalar yt;
  TacsScalar dyc_dx;
  TacsScalar theta;

  yt = 5.0*tt*(a1*x**0.5 - a2*x - a3*x**2 + a4*x**3 - a5*x**4);
  if ( (m > 0.0) && (p > 0.0) ){
    if ( x <= p ){
      yc = (m/p**2)*(2*p*x - x**2);
      dyc_dx = 2*m*(p - x)/(p**2);
    }
    else {
      yc = m*((1 - 2.0*p) + 2*p*x - x**2)/((1 - p)**2);
      dyc_dx = 2*m*(p - x)/((1 - p)**2);
    }
    theta = atan(dyc_dx);
    xl = x + yt*sin(theta);
    xu = x - yt*sin(theta);
    yl = yc_pts[i] - yt*cos(theta);
    yu = yc_pts[i] + yt*cos(theta);
  }
  else {
    yc = 0.0;
    xl = x;
    xu = x;
    yl = -yt;
    yu = yt;
  }
}

void TACSNACABeamConstitutive::computeDs( TacsScalar x,
                                          TacsScalar dsl,
                                          TacsScalar dsu ){

  TacsScalar yt = 5*tt*(a1*x**0.5 - a2*x - a3*x**2 + a4*x**3 - a5*x**4);
  TacsScalar dyt_dx = 5*tt*(0.5*a1*x**(-0.5) - a2 - 2.0*a3*x + 3.0*a4*x**2 - 4.0*a5*x**3);

  TacsScalar dyc_dx;
  TacsScalar d2yc_dx2;
  if ( (m > 0.0) && (p > 0.0) ){
    if ( x <= p ){
      dyc_dx = 2*m*(p - x)/(p**2);
      d2yc_dx2 = -2.0*m/p**2;
    }
    else {
      dyc_dx = 2*m*(p - x)/((1 - p)**2);
      d2yc_dx2 = -2.0*m/((1.0 - p)**2);
    }
  }
  else {
    dyc_dx = 0.0;
    d2yc_dx2 = 0.0;
  }

  // Lower curve
  TacsScalar dxl_dx = 1.0 + dyt_dx*sin(theta) + yt*cos(theta)*dtheta_dx;
  TacsScalar dyl_dx = dyc_dx - dyt_dx*cos(theta) + yt*sin(theta)*dtheta_dx;
  dsl = sqrt((dxl_dx)**2 + (dyl_dx)**2);

  // Upper curve
  TacsScalar dxu_dx = 1.0 - dyt_dx*sin(theta) - yt*cos(theta)*dtheta_dx;
  TacsScalar dyu_dx = dyc_dx + dyt_dx*cos(theta) - yt*sin(theta)*dtheta_dx;
  dsu = sqrt((dxu_dx)**2 + (dyu_dx)**2);
}

TacsScalar TACSNACABeamConstitutive::computeArea(){
  TacsScalar dx = 1.0/(npts-1);
  TacsScalar x = 0.5*dx;
  TacsScalar A = 0.0;
  TacsScalar dsl;
  TacsScalar dsu;
  for ( int i = 0; i < npts-1; i++ ){
    computeDs(x, dsl, dsu)
    A += wall*dx*(dlu+dsu);
  }
  A *= chord**2
  return A
}

TacsScalar TACSNACABeamConstitutive::computeTorsionConstant(){
  x = 0.0
  dx = 1.0/(npts-1);

  // Compute the enclosed area
  TacsScalar A_enc = 0.0;
  TacsScalar yc;
  TacsScalar xl; // dummy
  TacsScalar xu; // dummy
  TacsScalar xl_i;
  TacsScalar xu_i;
  TacsScalar yl; // dummy
  TacsScalar yu; // dummy
  TacsScalar xl_ip1;
  TacsScalar xu_ip1;
  TacsScalar dxi;  // variable dx width for integration
  for ( int i = 0; i < npts-1; i++ ){
    computeNACAPoint(x, yc, xl_i, xu_i, yl, yu);  // evaluate xl, xu at x[i]
    computeNACAPoint(x+dx, yc, xl_ip1, xu_ip1, yl, yu); // evaluate xl, xu at x[i+1]
    computeNACAPoint(x+0.5*dx, yc, xl, xu, yl, yu);  // evaluate yc at the midpoint
    dxi = 0.5*(xl_i+xu_i) - 0.5*(xl_ip1+xu_ip1);
    // ***** left off here
  }


  // # Compute the enclosed area
  //   x = np.linspace(0.0, 1.0, npts)
  //   A_enc = 0.0
  //   for i in range(npts-1):

  //       # Get the points on the curve
  //       zli,zui, _, _, _ = naca(x[i], m, p, tt)
  //       _, _, yl, yu, _ = naca(0.5*(x[i] + x[i+1]), m, p, tt)
  //       zlip1, zuip1, _, _, _ = naca(x[i+1], m, p, tt)
  //       dz = 0.5*(zli+zui) - 0.5*(zlip1+zuip1)

  //       A_enc += dz*(yu-yl)

  //   A = area(m, p, tt, tw, npts=npts)
  //   S = A/tw  // A = S*tw -> S = A/tw
    J = 4.0*(A_enc**2)*tw/S
    J *= chord**2
    return J
}

void TACSNACABeamConstitutive::computeCentroid(){

}

void TACSNACABeamConstitutive::computeMoments(){

}

void TACSNACABeamConstitutive::evalMassMoments( int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                TacsScalar moments[] ){
  TacsScalar rho = props->getDensity();
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1))/4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1))/64.0;

  moments[0] = rho * A;
  moments[1] = 0.0;
  moments[2] = 0.0;
  moments[3] = rho * Ia;
  moments[4] = rho * Ia;
  moments[5] = 0.0;
}

void TACSNACABeamConstitutive::addMassMomentsDVSens( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar scale[],
                                                     int dvLen,
                                                     TacsScalar dfdx[] ){
  TacsScalar rho = props->getDensity();
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;

  int index = 0;
  if (innerDV >= 0){
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1))/16.0;

    dfdx[index] += rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (wallDV >= 0){
    TacsScalar dA = M_PI * d0/2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0)/16.0;

    dfdx[index] += rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
}

TacsScalar TACSNACABeamConstitutive::evalSpecificHeat( int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[] ){
  if (props){
    return props->getSpecificHeat();
  }
  return 0.0;
}


TacsScalar TACSNACABeamConstitutive::evalDensity( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[] ){
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar rho = props->getDensity();
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1))/4.0;

  return rho * A;
}

void TACSNACABeamConstitutive::addDensityDVSens( int elemIndex,
                                                 TacsScalar scale,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 int dvLen,
                                                 TacsScalar dfdx[] ){
  int index = 0;
  TacsScalar d0 = inner + wall;
  TacsScalar rho = props->getDensity();

  if (innerDV >= 0){
    TacsScalar dA = M_PI * wall / 2.0;
    dfdx[index] += scale * rho * dA;
    index++;
  }
  if (wallDV >= 0){
    TacsScalar dA = M_PI * d0/2.0;
    dfdx[index] += scale * rho * dA;
    index++;
  }
}

void TACSNACABeamConstitutive::evalStress( int elemIndex,
                                           const double pt[],
                                           const TacsScalar X[],
                                           const TacsScalar e[],
                                           TacsScalar s[] ){
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5*E/(1.0 + nu);
  TacsScalar kcorr = 2.0*(1.0 + nu)/(4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1))/4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1))/64.0;

  s[0] = E * A * e[0];
  s[1] = 2.0 * G * Ia * e[1];
  s[2] = E * Ia * e[2];
  s[3] = E * Ia * e[3];
  s[4] = kcorr * G * A * e[4];
  s[5] = kcorr * G * A * e[5];
}

void TACSNACABeamConstitutive::evalTangentStiffness( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     TacsScalar C[] ){
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5*E/(1.0 + nu);
  TacsScalar kcorr = 2.0*(1.0 + nu)/(4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1))/4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1))/64.0;

  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES*sizeof(TacsScalar));

  C[0] = E * A;
  C[6] = 2.0 * G * Ia;
  C[11] = E * Ia;
  C[15] = E * Ia;
  C[18] = kcorr * G * A;
  C[20] = kcorr * G * A;
}

void TACSNACABeamConstitutive::addStressDVSens( int elemIndex,
                                                TacsScalar scale,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar e[],
                                                const TacsScalar psi[],
                                                int dvLen,
                                                TacsScalar dfdx[] ){
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5*E/(1.0 + nu);
  TacsScalar kcorr = 2.0*(1.0 + nu)/(4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;

  int index = 0;
  if (innerDV >= 0){
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1))/16.0;

    dfdx[index] += scale * (E * dA * e[0] * psi[0] +
                            2.0 * G * dIa * e[1] * psi[1] +
                            E * dIa * e[2] * psi[2] +
                            E * dIa * e[3] * psi[3] +
                            kcorr * G * dA * e[4] * psi[4] +
                            kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
  if (wallDV >= 0){
    TacsScalar dA = M_PI * d0/2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0)/16.0;

    dfdx[index] += scale * (E * dA * e[0] * psi[0] +
                            2.0 * G * dIa * e[1] * psi[1] +
                            E * dIa * e[2] * psi[2] +
                            E * dIa * e[3] * psi[3] +
                            kcorr * G * dA * e[4] * psi[4] +
                            kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
}

TacsScalar TACSNACABeamConstitutive::evalFailure( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar e[] ){
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3]; // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  return props->vonMisesFailure3D(s0);
}

TacsScalar TACSNACABeamConstitutive::evalFailureStrainSens( int elemIndex,
                                                            const double pt[],
                                                            const TacsScalar X[],
                                                            const TacsScalar e[],
                                                            TacsScalar sens[] ){
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3]; // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  TacsScalar s0d[6];
  TacsScalar fail = props->vonMisesFailure3DStressSens(s0, s0d);

  TacsScalar e0d[6];
  props->evalStress3D(s0d, e0d);

  sens[0] = e0d[0];
  sens[2] = r0 * e0d[0];
  sens[3] = r0 * e0d[0];
  sens[1] = r0 * e0d[3];
  sens[5] = e0d[4];
  sens[4] = e0d[5];

  return fail;
}

void TACSNACABeamConstitutive::addFailureDVSens( int elemIndex,
                                                 TacsScalar scale,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar e[],
                                                 int dvLen,
                                                 TacsScalar dfdx[] ){
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3]; // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  TacsScalar s0d[6];
  props->vonMisesFailure3DStressSens(s0, s0d);

  TacsScalar e0d[6];
  props->evalStress3D(s0d, e0d);

  int index = 0;
  if (innerDV >= 0){
    TacsScalar dr0 = 0.5;
    dfdx[index] += scale * (e0d[0] * dr0 * (e[2] + e[3]) +
                            e0d[3] * dr0 * e[1]);
    index++;
  }
  if (wallDV >= 0){
    TacsScalar dr0 = 1.0;
    dfdx[index] += scale * (e0d[0] * dr0 * (e[2] + e[3]) +
                            e0d[3] * dr0 * e[1]);

    index++;
  }
}

TacsScalar TACSNACABeamConstitutive::evalDesignFieldValue( int elemIndex,
                                                           const double pt[],
                                                           const TacsScalar X[],
                                                           int index ){
  if (index == 0){
    return inner;
  }
  else if (index == 1){
    return wall;
  }
  return 0.0;
}

const char* TACSNACABeamConstitutive::constName = "TACSNACABeamConstitutive";

/*
  Return the constitutive name
*/
const char* TACSNACABeamConstitutive::getObjectName(){
  return constName;
}
