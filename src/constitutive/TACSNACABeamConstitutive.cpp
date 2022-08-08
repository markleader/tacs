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
                                                    TacsScalar _yrot,
                                                    TacsScalar _zrot ){
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
  yrot = _yrot;
  zrot = _zrot;

  double a1 = 0.2969;
  double a2 = 0.1260;
  double a3 = 0.3516;
  double a4 = 0.2843;
  double a5 = 0.1015;
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
                                                  TacsScalar yl_pts[],
                                                  TacsScalar yu_pts[],
                                                  TacsScalar zl_pts[],
                                                  TacsScalar zu_pts[] ){

  TacsScalar x = 0.0;
  TacsScalar dx = 1.0/(npts-1);
  TacsScalar yt, dyc_dx, theta;
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
      yl_pts[i] = yc_pts[i] - yt*cos(theta);
      yu_pts[i] = yc_pts[i] + yt*cos(theta);
      zl_pts[i] = x + yt*sin(theta);
      zu_pts[i] = x - yt*sin(theta);
    }
    else {
      yc_pts[i] = 0.0;
      yl_pts[i] = -yt;
      yu_pts[i] = yt;
      zl_pts[i] = x;
      zu_pts[i] = x;
    }
  x += dx;
  }
}

void TACSNACABeamConstitutive::computeNACAPoint( TacsScalar x,
                                                 TacsScalar yc,
                                                 TacsScalar yl,
                                                 TacsScalar yu,
                                                 TacsScalar zl,
                                                 TacsScalar zu ){

  TacsScalar yt, dyc_dx, theta;

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
    yl = yc_pts[i] - yt*cos(theta);
    yu = yc_pts[i] + yt*cos(theta);
    zl = x + yt*sin(theta);
    zu = x - yt*sin(theta);
  }
  else {
    yc = 0.0;
    yl = -yt;
    yu = yt;
    zl = x;
    zu = x;
  }
}

void TACSNACABeamConstitutive::computeDs( TacsScalar x,
                                          TacsScalar dsl,
                                          TacsScalar dsu ){

  TacsScalar yt = 5*tt*(a1*x**0.5 - a2*x - a3*x**2 + a4*x**3 - a5*x**4);
  TacsScalar dyt_dx = 5*tt*(0.5*a1*x**(-0.5) - a2 - 2.0*a3*x + 3.0*a4*x**2 - 4.0*a5*x**3);

  TacsScalar dyc_dx, d2yc_dx2;
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
  TacsScalar dzl_dx = 1.0 + dyt_dx*sin(theta) + yt*cos(theta)*dtheta_dx;
  TacsScalar dyl_dx = dyc_dx - dyt_dx*cos(theta) + yt*sin(theta)*dtheta_dx;
  dsl = sqrt((dzl_dx)**2 + (dyl_dx)**2);

  // Upper curve
  TacsScalar dzu_dx = 1.0 - dyt_dx*sin(theta) - yt*cos(theta)*dtheta_dx;
  TacsScalar dyu_dx = dyc_dx + dyt_dx*cos(theta) - yt*sin(theta)*dtheta_dx;
  dsu = sqrt((dzu_dx)**2 + (dyu_dx)**2);
}

TacsScalar TACSNACABeamConstitutive::computeArea(){
  TacsScalar dx = 1.0/(npts-1);
  TacsScalar x = 0.5*dx;
  TacsScalar A = 0.0;
  TacsScalar dsl, dsu;
  for ( int i = 0; i < npts-1; i++ ){
    computeDs(x, dsl, dsu)
    A += wall*dx*(dlu+dsu);
  }
  A *= chord**2
  return A
}

TacsScalar TACSNACABeamConstitutive::computeTorsionConstant(){
  TacsScalar x = 0.0;
  TacsScalar dx = 1.0/(npts-1);

  // Compute the enclosed area
  TacsScalar A_enc = 0.0;
  TacsScalar yc;
  TacsScalar zl, zu; // dummy
  TacsScalar zl_i, zu_i;
  TacsScalar yl, yu;
  TacsScalar zl_ip1, zu_ip1;
  TacsScalar dzi;  // variable dx width for integration
  for ( int i = 0; i < npts-1; i++ ){
    computeNACAPoint(x, yc, yl, yu, zl_i, zu_i);  // evaluate xl, xu at x[i]
    computeNACAPoint(x+dx, yc, yl, yu, zl_ip1, zu_ip1,); // evaluate xl, xu at x[i+1]
    computeNACAPoint(x+0.5*dx, yc, yl, yu, zl, zu);  // evaluate yc at the midpoint
    dzi = 0.5*(zl_i+zu_i) - 0.5*(zl_ip1+zu_ip1);
    A_enc += dzi*(yu-yl);
  }

  TacsScalar A = computeArea();
  TacsScalar S = A/wall;
  TacsScalar J = 4.0*(A_enc**2)*tw/S
  J *= chord**2
  return J
}

void TACSNACABeamConstitutive::computeRelCentroid( TacsScalar ystar,
                                                   TacsScalar zstar ){
  TacsScalar x = 0.0;
  TacsScalar dx = 1.0/(npts-1);
  ystar = 0.0;
  zstar = 0.0;

  TacsScalar A = computeArea();
  TacsScalar dsl, dsu;
  TacsScalar ystar, xl, xu, yl, yu;

  for ( int i = 0; i < npts-1; i++ ){
    computeDs(x+0.5*dx, dsl, dsu);
    computeNACAPoint(x+0.5*dx, ystar, xl, xu, yl, yu);

    // Add the contribution from the centroid of the upper dA at this x
    ystar += wall*dsu*dx*yu;  // dA = wall*dsu*dx
    zstar += wall*dsu*dx*zu;

    // Add the contribution from the centroid of the lower dA at this x
    ystar += wall*dsl*dx*yl;  // dA = wall*dsl*dx
    zstar += wall*dsl*dx*zl;

    x += dx;
  }

  ystar /= A;
  zstar /= A;
}

void TACSNACABeamConstitutive::computeRefMoments( TacsScalar Iyy,
                                                  TacsScalar Izz,
                                                  TacsScalar Iyz ){
  TacsScalar x = 0.0;
  TacsScalar dx = 1.0/(npts-1);
  TacsScalar A = computeArea();
  TacsScalar Iyy = 0.0;
  TacsScalar Izz = 0.0;
  TacsScalar Iyz = 0.0;
  TacsScalar dsl, dsu;
  TacsScalar yc, xl, xu, yl, yu;
  TacsScalar dA;

  // Compute the moments of inertia about the leading edge
  for ( int i = 0; i < npts-1; i++ ){
    computeDs(x+0.5*dx, dsl, dsu);
    computeNACAPoint(x+0.5*dx, yc, xl, xu, yl, yu);

    // Add the inertia contributions from the upper curve dA
    dA = wall*dsu*dx;
    Iyy += dA*zu**2;
    Izz += dA*yu**2;
    Iyz += dA*yu*zu;

    // Add the inertia contributions from the lower curve dA
    dA = tw*dsl*dx;
    Iyy += dA*zl**2;
    Izz += dA*yl**2;
    Iyz += dA*yl*zl;
  }

  // Get the centroidal moments of inertia
  TacsScalar ystar, zstar;
  computeRelCentroid(ystar, zstar);
  TacsScalar Iyy_c = Iyy - A*zstar**2;
  TacsScalar Izz_c = Izz - A*ystar**2;
  TacsScalar Iyz_c = Iyz - A*ystar*zstar;

  // Apply parallel axis theorem if the moments are not about the centroid
  TacsScalar dy = yrot - ystar;
  TacsScalar dz = zrot - zstar;
  if (!use_cm){
    Iyy = Iyy_c + A*dz**2;
    Izz = Izz_c + A*dy**2;
    Iyz = Iyz_c + A*dy*dz;
  }

  Iyy *= chord**4;
  Izz *= chord**4;
  Iyz *= chord**4;
}

void TACSNACABeamConstitutive::computeMoments( TacsScalar Iyy,
                                               TacsScalar Izz,
                                               TacsScalar Iyz ){
  TacsScalar Iyy0, Izz0, Iyz0;
  computeRefMoments(Iyy0, Izz0, Iyz0);

  // Rotate the moments of inertia by the twist angle of the airfoil
  TacsScalar ystar, zstar;
  computeRelCentroid(ystar, zstar);
  TacsScalar dz = 0.0;
  if (!use_cm){
    dz = zrot - zstar;
  }
  // Additional term for dy???

  Iyy = c2*Iyy0 + s2*Izz0 + A*(-dz*chord*c)**2;
  Izz = s2*Iyy0 + c2*Izz0 + A*(dz*chord*s)**2;
  Iyz = (Izz0 - Iyy0)*s*c + A*(dz*chord*c)*(-dz*chord*c);  // *** Update for nonzero Iyz0
}

void TACSNACABeamConstitutive::evalMassMoments( int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                TacsScalar moments[] ){
  TacsScalar rho = props->getDensity();
  TacsScalar A = computeArea();
  TacsScalar ystar, zstar;
  TacsScalar Iyy, Izz, Iyz;
  computeRelCentroid(ystar, zstar);
  computeMoments(Iyy, Izz, Iyz);

  TacsScalar dy = 0.0;
  TacsScalar dz = 0.0;
  if (!use_cm){
    dy = -(zstar - zrot)*sin(twist) + (ystar - yrot)*cos(twist);
    dz = (zstar - zrot)*cos(twist) + (ystar - yrot)*sin(twist);
  }

  moments[0] = rho*A;
  moments[1] = rho*A*chord*dy;
  moments[2] = rho*A*chord*dz;
  moments[3] = rho*Iyy;
  moments[4] = rho*Izz;
  moments[5] = rho*Iyz;
}

void TACSNACABeamConstitutive::addMassMomentsDVSens( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar scale[],
                                                     int dvLen,
                                                     TacsScalar dfdx[] ){
  TacsScalar rho = props->getDensity();
  TacsScalar A = computeArea();
  TacsScalar ystar, zstar;
  TacsScalar Iyy0, Izz0, Iyz0;
  TacsScalar Iyy, Izz, Iyz;
  computeRelCentroid(ystar, zstar);
  computeRefMoments(Iyy0, Izz0, Iyz0);
  computeMoments(Iyy, Izz, Iyz);
  TacsScalar dA, dystar, dzstar, dAy, dAz, dIyy, dIzz, dIyz;

  TacsScalar dz = 0.0;
  if (!use_cm){
    dz = zrot - zstar;
  }

  TacsScalar s = sin(theta);
  TacsScalar c = cos(theta);
  TacsScalar s2 = s**2;
  TacsScalar c2 = c**2;

  int index = 0;
  if (chordDV >= 0){
    dA = (2.0/chord)*A;
    dIyy = (4.0/chord)*(Iyy0*c2 + Izz0*s2) + 4.0*A*chord*(-dz*chord*c)**2;
    dIzz = (4.0/chord)*(Iyy0*s2 + Izz0*c2) + 4.0*A*chord*(-dz*chord*s)**2;
    dIyz = (4.0/chord)*(Izz0 - Iyy0)*c*s - 4.0*A*chord*s*c*(-dz*chord)**2;
    dystar = ystar/chord;
    dzstar = zstar/chord;
    dAy = dA*ystar + A*dystar;
    dAz = dA*zstar + A*dzstar;

    dfdx[index] += rho * ( scale[0]*dA + scale[1]*dAy + scale[2]*dAz );
    dfdx[index] += rho * ( scale[3]*dIyy + scale[4]*dIzz + scale[5]*dIyz );
    index++;
  }
  if (twistDV >= 0){
    TacsScalar ddy, ddz;
    TacsScalar ds2 = 2.0*s*c;
    TacsScalar dc2 = -2.0*s*c;
    TacsScalar dcs = c2 - s2;
    ddy = -(zstar - zrot)*cos(twist) - (ystar - yrot)*sin(twist);
    ddz = -(zstar - zrot)*sin(twist) + (ystar - yrot)*cos(twist);

    dA = 0.0;
    dAy = rho*A*chord*ddy;
    dAz = rho*A*chord*ddz;

    dIyy = dc2*Iyy0 + ds2*Izz0 + A*dc2*(dz*chord)**2;
    dIzz = ds2*Iyy0 + dc2*Izz0 + A*ds2*(dz*chord)**2;
    dIyz = (Izz0 - Iyy0)*dcs - A*dcs*(dz*chord)**2;

    dfdx[index] += rho * ( scale[0]*dA + scale[1]*dAy + scale[2]*dAz );
    dfdx[index] += rho * ( scale[3]*dIyy + scale[4]*dIzz + scale[5]*dIyz );
    index++;
  }
  // *** Left off here
  if (wallDV >= 0){
    dA = A/wall;
    dAy = 0.0; // *** TODO: dA*dy + A*ddy
    dAz = 0.0; // *** TODO: dA*dz + A*ddz
    dIyy = Iyy/wall;
    dIzz = Izz/wall;
    dIyz = Iyz/wall;

    dfdx[index] += rho * ( scale[0]*dA + scale[1]*dAy + scale[2]*dAz );
    dfdx[index] += rho * ( scale[3]*dIyy + scale[4]*dIzz + scale[5]*dIyz );
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
