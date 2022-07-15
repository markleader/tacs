/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H
#define TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H

#include "TACSConstitutive.h"
#include "TACSMaterialProperties.h"

/*
  This is the base class for the phase change material constitutive objects.
*/
class TACSPhaseChangeMaterialConstitutive : public TACSConstitutive {
 public:

  TACSPhaseChangeMaterialConstitutive( TACSMaterialProperties *solid_props,
                                       TACSMaterialProperties *liquid_props,
                                       TacsScalar _lh, TacsScalar _mt,
                                       TacsScalar _t=1.0, int _tNum=-1,
                                       TacsScalar _tlb=0.0, TacsScalar _tub=1.0 );
  ~TACSPhaseChangeMaterialConstitutive();

  // Retrieve the global design variable numbers
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  // Set the element design variable from the design vector
  int setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  // Get the element design variables values
  int getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange( int elemIndex, int dvLen,
                         TacsScalar lb[], TacsScalar ub[] );

  // Evaluate the temperature at a given element
  TacsScalar evalTemperature( int elemIndex, const double pt[],
                              const TacsScalar X[], const TacsScalar U );

  // Evaluate the material's phase
  int evalPhase( const TacsScalar U);

  // Evaluate the material density
  TacsScalar evalDensity( int elemIndex, const double pt[],
                          const TacsScalar X[] );

  // Add the derivative of the density
  void addDensityDVSens( int elemIndex, TacsScalar scale,
                         const double pt[], const TacsScalar X[],
                         int dvLen, TacsScalar dfdx[] );

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat( int elemIndex, const double pt[],
                               const TacsScalar X[], const TacsScalar U );

  // Evaluate the thermal strain
  void evalThermalStrain( int elemIndex, const double pt[],
                          const TacsScalar X[], TacsScalar theta,
                          TacsScalar strain[], const TacsScalar U );

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux( int elemIndex, const double pt[],
                     const TacsScalar X[], const TacsScalar grad[],
                     TacsScalar flux[], const TacsScalar U );

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux( int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar Kc[],
                            const TacsScalar U );

  // Add the derivative of the heat flux
  void addHeatFluxDVSens( int elemIndex, TacsScalar scale,
                          const double pt[], const TacsScalar X[],
                          const TacsScalar grad[], const TacsScalar psi[],
                          int dvLen, TacsScalar dfdx[], const TacsScalar U );

  // Extra info about the constitutive class
  const char *getObjectName();

 protected:
  // Materiial properties class
  TACSMaterialProperties *solid_props;
  TACSMaterialProperties *liquid_props;

 private:
  // Store information about the design variable
  TacsScalar lh, mt, t, tlb, tub;
  int tNum;

  static const char *psName;
};

#endif // TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H
