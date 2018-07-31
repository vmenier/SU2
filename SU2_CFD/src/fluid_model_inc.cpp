/*!
 * fluid_model_inc.cpp
 * \brief Source of the incompressible fluid models.
 * \author T. Economon
 * \version 6.0.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/fluid_model.hpp"

CConstantDensity::CConstantDensity() : CFluidModel() {
  Density = 0.0;
  Cp      = 0.0;
  Cv      = 0.0;
}

CConstantDensity::CConstantDensity(su2double val_Density,
                                   su2double val_Cp) : CFluidModel() {
  Density = val_Density;
  Cp      = val_Cp;
  Cv      = val_Cp;
}

CConstantDensity::~CConstantDensity(void) { }

void CConstantDensity::SetTDState_T (su2double val_Temperature) {
  
  /*--- Density is constant and thermodynamic pressure is
   not required for incompressible, constant density flows,
   but the energy equation can still be computed as a
   decoupled equation. Hence, we update the value.
   Note Cp = Cv (gamma = 1). ---*/
  
  Temperature = val_Temperature;
  
}

CIncIdealGas::CIncIdealGas() : CFluidModel() {
  Pressure        = 0.0;
  Gamma           = 0.0;
  Gas_Constant    = 0.0;
  Cp              = 0.0;
  Cv              = 0.0;
}

CIncIdealGas::CIncIdealGas(su2double val_Cp,
                           su2double val_gas_constant,
                           su2double val_operating_pressure) : CFluidModel() {

  /*--- In the incompressible ideal gas model, the thermodynamic pressure
    is decoupled from the governing equations and held constant. The 
    density is therefore only a function of temperature variations. ---*/

  Gas_Constant = val_gas_constant;
  Pressure     = val_operating_pressure;
  Gamma        = 1.0;
  Cp           = val_Cp;
  Cv           = Cp;

}

CIncIdealGas::~CIncIdealGas(void) { }

void CIncIdealGas::SetTDState_T(su2double val_temperature) {

 /*--- The EoS only depends upon temperature. ---*/

  Temperature  = val_temperature;
  Density      = Pressure/(Temperature*Gas_Constant);

}
