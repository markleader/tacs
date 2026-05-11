#ifndef TACS_BEAM2_PINNED_H
#define TACS_BEAM2_PINNED_H

#include "TACSBeamCentrifugalForce.h"
#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamElementTransform.h"
#include "TACSBeamInertialForce.h"
#include "TACSBeamTraction.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "a2d.h"

class TACSBeam2Pinned : public TACSElement {
 public:
  static const int offset = 3;
  static const int vars_per_node = 6;
  static const int num_nodes = 2;
  static const int num_variables = vars_per_node * num_nodes;

  TACSBeam2Pinned(TACSBeamTransform *_transform, TACSBeamConstitutive *_con) {
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  ~TACSBeam2Pinned() {
    if (transform) {
      transform->decref();
    }

    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSBeam2Pinned"; }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return num_nodes; }

  ElementLayout getLayoutType() { return TACS_LINE_ELEMENT; }

  ElementType getElementType() { return TACS_BEAM_OR_SHELL_ELEMENT; }

  int getNumQuadraturePoints() {
    return TACSBeamLinearQuadrature::getNumQuadraturePoints();
  }

  double getQuadratureWeight(int n) {
    return TACSBeamLinearQuadrature::getQuadratureWeight(n);
  }

  double getQuadraturePoint(int n, double pt[]) {
    return TACSBeamLinearQuadrature::getQuadraturePoint(n, pt);
  }

  int getNumElementFaces() {
    return TACSBeamLinearQuadrature::getNumElementFaces();
  }

  int getNumFaceQuadraturePoints(int face) {
    return TACSBeamLinearQuadrature::getNumFaceQuadraturePoints(face);
  }

  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return TACSBeamLinearQuadrature::getFaceQuadraturePoint(face, n, pt,
                                                            tangent);
  }

  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
    return con->getDesignVarNums(elemIndex, dvLen, dvNums);
  }

  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    return con->setDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    return con->getDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]) {
    return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
  }

  TACSElement *createElementTraction(int faceIndex, const TacsScalar t[]) {
    return new TACSBeamTraction<vars_per_node, TACSBeamLinearQuadrature,
                                TACSBeamBasis<2> >(t);
  }

  TACSElement *createElementInertialForce(const TacsScalar inertiaVec[]) {
    return new TACSBeamInertialForce<vars_per_node, TACSBeamLinearQuadrature,
                                     TACSBeamBasis<2> >(transform, con,
                                                        inertiaVec);
  }

  TACSElement *createElementCentrifugalForce(const TacsScalar omegaVec[],
                                             const TacsScalar rotCenter[],
                                             const bool first_order = false) {
    return new TACSBeamCentrifugalForce<vars_per_node,
                                        TACSBeamLinearQuadrature,
                                        TACSBeamBasis<2> >(transform, con,
                                                           omegaVec,
                                                           rotCenter);
  }

  void computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                       const TacsScalar vars[], const TacsScalar dvars[],
                       TacsScalar *Te, TacsScalar *Pe) {
    *Te = 0.0;
    *Pe = 0.0;

    TacsScalar uaxial[3];
    TacsScalar X0[3], t1[3], T[9];
    TacsScalar detXd = 0.0;
    computeAxis(Xpts, X0, t1, T, &detXd);

    const TacsScalar L = 2.0 * detXd;
    if (TacsRealPart(L) <= 0.0) {
      return;
    }

    TacsScalar e[6] = {0.0};
    computeAxialStrain(vars, t1, L, &e[0]);

    for (int n = 0; n < getNumQuadraturePoints(); n++) {
      double pt[3];
      const double weight = getQuadraturePoint(n, pt);

      TacsScalar X[3];
      interpPosition(pt[0], Xpts, X);

      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X, e, s);
      *Pe += 0.5 * weight * detXd * e[0] * s[0];

      const TacsScalar N0 = 0.5 * (1.0 - pt[0]);
      const TacsScalar N1 = 0.5 * (1.0 + pt[0]);
      TacsScalar density = con->evalDensity(elemIndex, pt, X);
      TacsScalar v[3];
      for (int i = 0; i < 3; i++) {
        v[i] = N0 * dvars[i] + N1 * dvars[vars_per_node + i];
      }
      *Te += 0.5 * weight * detXd * density *
             (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
  }

  void addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[]) {
    TacsScalar X0[3], t1[3], T[9];
    TacsScalar detXd = 0.0;
    computeAxis(Xpts, X0, t1, T, &detXd);

    const TacsScalar L = 2.0 * detXd;
    if (TacsRealPart(L) <= 0.0) {
      return;
    }

    TacsScalar e[6] = {0.0};
    computeAxialStrain(vars, t1, L, &e[0]);

    for (int n = 0; n < getNumQuadraturePoints(); n++) {
      double pt[3];
      const double weight = getQuadraturePoint(n, pt);

      TacsScalar X[3];
      interpPosition(pt[0], Xpts, X);

      const TacsScalar N0 = 0.5 * (1.0 - pt[0]);
      const TacsScalar N1 = 0.5 * (1.0 + pt[0]);
      const TacsScalar bscale = weight * detXd / L;

      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X, e, s);
      for (int i = 0; i < 3; i++) {
        res[i] -= bscale * s[0] * t1[i];
        res[vars_per_node + i] += bscale * s[0] * t1[i];
      }

      const TacsScalar density = con->evalDensity(elemIndex, pt, X);
      for (int i = 0; i < 3; i++) {
        const TacsScalar accel =
            N0 * ddvars[i] + N1 * ddvars[vars_per_node + i];
        res[i] += weight * detXd * density * N0 * accel;
        res[vars_per_node + i] += weight * detXd * density * N1 * accel;
      }
    }
  }

  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]) {
    TacsScalar X0[3], t1[3], T[9];
    TacsScalar detXd = 0.0;
    computeAxis(Xpts, X0, t1, T, &detXd);

    const TacsScalar L = 2.0 * detXd;
    if (TacsRealPart(L) <= 0.0) {
      return;
    }

    TacsScalar e[6] = {0.0};
    computeAxialStrain(vars, t1, L, &e[0]);

    for (int n = 0; n < getNumQuadraturePoints(); n++) {
      double pt[3];
      const double weight = getQuadraturePoint(n, pt);

      TacsScalar X[3];
      interpPosition(pt[0], Xpts, X);

      const TacsScalar N[2] = {0.5 * (1.0 - pt[0]), 0.5 * (1.0 + pt[0])};
      const TacsScalar bscale = weight * detXd / L;

      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X, e, s);

      if (res) {
        for (int i = 0; i < 3; i++) {
          res[i] -= bscale * s[0] * t1[i];
          res[vars_per_node + i] += bscale * s[0] * t1[i];
        }
      }

      TacsScalar C[21];
      con->evalTangentStiffness(elemIndex, pt, X, C);
      const TacsScalar kaxial = alpha * bscale * C[0] / L;
      for (int node_i = 0; node_i < 2; node_i++) {
        const TacsScalar sign_i = (node_i == 0) ? -1.0 : 1.0;
        const int row0 = node_i * vars_per_node;
        for (int node_j = 0; node_j < 2; node_j++) {
          const TacsScalar sign_j = (node_j == 0) ? -1.0 : 1.0;
          const int col0 = node_j * vars_per_node;
          const TacsScalar scale = sign_i * sign_j * kaxial;
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              mat[(col0 + j) + (row0 + i) * num_variables] +=
                  scale * t1[i] * t1[j];
            }
          }
        }
      }

      const TacsScalar density = con->evalDensity(elemIndex, pt, X);
      if (res) {
        for (int i = 0; i < 3; i++) {
          const TacsScalar accel =
              N[0] * ddvars[i] + N[1] * ddvars[vars_per_node + i];
          res[i] += weight * detXd * density * N[0] * accel;
          res[vars_per_node + i] += weight * detXd * density * N[1] * accel;
        }
      }
      const TacsScalar mscale = gamma * weight * detXd * density;
      for (int node_i = 0; node_i < 2; node_i++) {
        const int row0 = node_i * vars_per_node;
        for (int node_j = 0; node_j < 2; node_j++) {
          const int col0 = node_j * vars_per_node;
          const TacsScalar mass = mscale * N[node_i] * N[node_j];
          for (int i = 0; i < 3; i++) {
            mat[(col0 + i) + (row0 + i) * num_variables] += mass;
          }
        }
      }
    }
  }

  int evalPointQuantity(int elemIndex, int quantityType, double time, int n,
                        double pt[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar *detXdval,
                        TacsScalar *quantity) {
    TacsScalar X0[3], t1[3], T[9];
    TacsScalar detXd = 0.0;
    computeAxis(Xpts, X0, t1, T, &detXd);

    if (detXdval) {
      *detXdval = detXd;
    }

    TacsScalar X[3];
    interpPosition(pt[0], Xpts, X);

    if (quantityType == TACS_ELEMENT_DENSITY) {
      if (quantity) {
        *quantity = con->evalDensity(elemIndex, pt, X);
      }
      return 1;
    } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
      if (quantity) {
        const TacsScalar N0 = 0.5 * (1.0 - pt[0]);
        const TacsScalar N1 = 0.5 * (1.0 + pt[0]);
        for (int i = 0; i < 3; i++) {
          quantity[i] = N0 * vars[i] + N1 * vars[vars_per_node + i];
        }
      }
      return 3;
    } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
      if (quantity) {
        TacsScalar mass_moment[6];
        con->evalMassMoments(elemIndex, pt, X, mass_moment);
        const TacsScalar density = mass_moment[0];
        const TacsScalar n1[3] = {T[1], T[4], T[7]};
        const TacsScalar n2[3] = {T[2], T[5], T[8]};
        for (int i = 0; i < 3; i++) {
          quantity[i] = density * X[i] - mass_moment[1] * n1[i] -
                        mass_moment[2] * n2[i];
        }
      }
      return 3;
    } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
      if (quantity) {
        TacsScalar moments[6];
        con->evalMassMoments(elemIndex, pt, X, moments);
        const TacsScalar density = moments[0];
        const TacsScalar I0[6] = {
            0.0,
            0.0,
            0.0,
            moments[4] - moments[2] * moments[2] / density,
            -moments[5] + moments[1] * moments[2] / density,
            moments[3] - moments[1] * moments[1] / density};
        mat3x3SymmTransform(T, I0, quantity);

        const TacsScalar n1[3] = {T[1], T[4], T[7]};
        const TacsScalar n2[3] = {T[2], T[5], T[8]};
        TacsScalar dXcg[3];
        for (int i = 0; i < 3; i++) {
          dXcg[i] = X[i] - (moments[1] * n1[i] + moments[2] * n2[i]) / density;
        }

        quantity[0] += density * (dXcg[1] * dXcg[1] + dXcg[2] * dXcg[2]);
        quantity[1] += -density * dXcg[0] * dXcg[1];
        quantity[2] += -density * dXcg[0] * dXcg[2];
        quantity[3] += density * (dXcg[0] * dXcg[0] + dXcg[2] * dXcg[2]);
        quantity[4] += -density * dXcg[1] * dXcg[2];
        quantity[5] += density * (dXcg[0] * dXcg[0] + dXcg[1] * dXcg[1]);
      }
      return 6;
    }

    const TacsScalar L = 2.0 * detXd;
    if (TacsRealPart(L) <= 0.0) {
      return 0;
    }

    TacsScalar e[6] = {0.0};
    computeAxialStrain(vars, t1, L, &e[0]);
    if (quantityType == TACS_ELEMENT_STRAIN) {
      if (quantity) {
        for (int i = 0; i < 6; i++) {
          quantity[i] = e[i];
        }
      }
      return 6;
    }

    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X, e, s);

    if (quantityType == TACS_ELEMENT_STRESS) {
      if (quantity) {
        for (int i = 0; i < 6; i++) {
          quantity[i] = s[i];
        }
      }
      return 6;
    } else if (quantityType == TACS_FAILURE_INDEX) {
      if (quantity) {
        *quantity = con->evalFailure(elemIndex, pt, X, e);
      }
      return 1;
    } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
      if (quantity) {
        *quantity = e[0] * s[0];
      }
      return 1;
    }

    return 0;
  }

 private:
  static void interpPosition(TacsScalar xi, const TacsScalar Xpts[],
                             TacsScalar X[]) {
    const TacsScalar N0 = 0.5 * (1.0 - xi);
    const TacsScalar N1 = 0.5 * (1.0 + xi);
    for (int i = 0; i < 3; i++) {
      X[i] = N0 * Xpts[i] + N1 * Xpts[3 + i];
    }
  }

  void computeAxis(const TacsScalar Xpts[], TacsScalar X0[], TacsScalar t1[],
                   TacsScalar T[], TacsScalar *detXd) const {
    for (int i = 0; i < 3; i++) {
      X0[i] = 0.5 * (Xpts[i] + Xpts[3 + i]);
    }

    TacsScalar Xxi[3];
    for (int i = 0; i < 3; i++) {
      Xxi[i] = 0.5 * (Xpts[3 + i] - Xpts[i]);
    }

    TacsScalar norm = sqrt(Xxi[0] * Xxi[0] + Xxi[1] * Xxi[1] +
                           Xxi[2] * Xxi[2]);
    *detXd = norm;
    if (TacsRealPart(norm) > 0.0) {
      t1[0] = Xxi[0] / norm;
      t1[1] = Xxi[1] / norm;
      t1[2] = Xxi[2] / norm;
    } else {
      t1[0] = 1.0;
      t1[1] = 0.0;
      t1[2] = 0.0;
    }

    transform->computeTransform(Xxi, T);
  }

  static void computeAxialStrain(const TacsScalar vars[], const TacsScalar t1[],
                                 TacsScalar L, TacsScalar *eps) {
    *eps = 0.0;
    for (int i = 0; i < 3; i++) {
      *eps += t1[i] * (vars[vars_per_node + i] - vars[i]);
    }
    *eps /= L;
  }

  TACSBeamTransform *transform;
  TACSBeamConstitutive *con;
};

#endif  // TACS_BEAM2_PINNED_H
