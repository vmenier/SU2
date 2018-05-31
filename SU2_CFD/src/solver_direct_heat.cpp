/*!
 * \file solution_direct_heat.cpp
 * \brief Main subrotuines for solving the heat equation
 * \author F. Palacios, T. Economon
 * \version 6.0.1 "Falcon"
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

#include "../include/solver_structure.hpp"

CHeatSolver::CHeatSolver(void) : CSolver() { }

CHeatSolver::CHeatSolver(CGeometry *geometry, CConfig *config) : CSolver() {
  
  unsigned short iDim, iVar, nLineLets;
  unsigned long iPoint;
  
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
  nDim    =       geometry->GetnDim();
  node    =       new CVariable*[nPoint];
  nVar    =       1;
  
  Residual     = new su2double[nVar]; Residual_RMS = new su2double[nVar];
  Solution     = new su2double[nVar];
  Res_Sour     = new su2double[nVar];
  Residual_Max = new su2double[nVar]; Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Point to point stiffness matrix (only for triangles)---*/
  
  StiffMatrix_Elem = new su2double*[nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++) {
    StiffMatrix_Elem[iVar] = new su2double [nDim+1];
  }
  
  StiffMatrix_Node = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    StiffMatrix_Node[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of matrix structures ---*/
  
  StiffMatrixSpace.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  StiffMatrixTime.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Linear Elasticity)." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  /*--- Initialization of linear solver structures ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Heat coefficient for all of the markers ---*/
  
  CHeat = new su2double[config->GetnMarker_All()];
  Total_CHeat = 0.0;
  
  /*--- Check for a restart (not really used), initialize from zero otherwise ---*/
  
  bool restart = (config->GetRestart());
  if (!restart) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Zero initial condition for testing source terms & forcing BCs ---*/
      
      Solution[0] = 0.0;
      Solution[1] = 0.0;
      
      node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
      
      /*--- Copy solution to old containers if using dual time ---*/
      
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) {
        node[iPoint]->Set_Solution_time_n();
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {
        node[iPoint]->Set_Solution_time_n();
        node[iPoint]->Set_Solution_time_n1();
      }
      
    }
  } else {
    
    cout << "Heat restart file not currently configured!!" << endl;
    
    string mesh_filename = config->GetSolution_FlowFileName();
    ifstream restart_file;
    
    char *cstr; cstr = new char [mesh_filename.size()+1];
    strcpy (cstr, mesh_filename.c_str());
    restart_file.open(cstr, ios::in);
    
    if (restart_file.fail()) {
      SU2_MPI::Error("There is no Heat restart file", CURRENT_FUNCTION);
    }
    unsigned long index;
    string text_line;
    
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      getline(restart_file, text_line);
      istringstream point_line(text_line);
      point_line >> index >> Solution[0] >> Solution[1];
      node[iPoint] = new CHeatVariable(Solution, nDim, nVar, config);
    }
    restart_file.close();
  }
  
}

CHeatSolver::~CHeatSolver(void) {
  
  unsigned short iVar;
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] StiffMatrix_Elem[iVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] StiffMatrix_Node[iVar];
  
  delete [] StiffMatrix_Elem;
  delete [] StiffMatrix_Node;
  
}

void CHeatSolver::Preprocessing(CGeometry *geometry,
                                CSolver **solver_container,
                                CConfig   *config,
                                unsigned short iMesh,
                                unsigned short iRKStep,
                                unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint;
  
  /*--- Set residuals and matrix entries to zero ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysSol.SetBlock_Zero(iPoint);
    LinSysAux.SetBlock_Zero(iPoint);
    LinSysRes.SetBlock_Zero(iPoint);
  }
  
  /*--- Zero out the entries in the various matrices ---*/
  StiffMatrixSpace.SetValZero();
  StiffMatrixTime.SetValZero();
  Jacobian.SetValZero();
  
}

void CHeatSolver::Source_Residual(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics, CNumerics *second_numerics,
                                  CConfig   *config,
                                  unsigned short iMesh) {

  if (config->GetUnsteady_Simulation() != STEADY) {

    unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
    su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, d[3] = {0.0,0.0,0.0}, Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
    su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
    unsigned short iDim;

    /*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
    if (config->GetUnsteady_Simulation() == TIME_STEPPING) Time_Num = config->GetDelta_UnstTimeND();
    else Time_Num = 1E30;
    
    /*--- Loop through elements to compute contributions from the matrix
     blocks involving time. These contributions are also added to the
     Jacobian w/ the time step. Spatial source terms are also computed. ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      /*--- Get node numbers and their coordinate vectors ---*/
      
      Point_0 = geometry->elem[iElem]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      
      /*--- Compute area and volume ---*/
      
      if (nDim == 2) {
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
        }
        Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
      }
      else {
        Point_3 = geometry->elem[iElem]->GetNode(3);
        Coord_3 = geometry->node[Point_3]->GetCoord();
        for (iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = Coord_0[iDim]-Coord_2[iDim];
          b[iDim] = Coord_1[iDim]-Coord_2[iDim];
          c[iDim] = Coord_3[iDim]-Coord_2[iDim];
        }
        d[0] = a[1]*b[2]-a[2]*b[1];
        d[1] = -(a[0]*b[2]-a[2]*b[0]);
        d[2] = a[0]*b[1]-a[1]*b[0];
        Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
      }
      
      /*--- Block contributions to the Jacobian (includes time step) ---*/
      
      if (nDim == 2) { StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local/Time_Num); }
      else { StiffMatrix_Node[0][0] = (2.0/20.0)*(Volume_Local/Time_Num); }
      Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      if (nDim == 3) Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
      
      if (nDim == 2) { StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local/Time_Num); }
      else { StiffMatrix_Node[0][0] = (1.0/20.0)*(Volume_Local/Time_Num); }
      
      Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      if (nDim == 3) {
        Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
        Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
      }
      
    }
    
  }
  
}

void CHeatSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                  CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0, total_index, iPoint;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3 = NULL;
  
  if (nDim == 2 ) {
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      Point_0 = geometry->elem[iElem]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
      
      numerics->SetCoord(Coord_0, Coord_1, Coord_2);
      numerics->ComputeResidual(StiffMatrix_Elem, config);
      
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      
    }
  }
  
  if (nDim == 3 ) {
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      Point_0 = geometry->elem[iElem]->GetNode(0);   Coord_0 = geometry->node[Point_0]->GetCoord();
      Point_1 = geometry->elem[iElem]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
      Point_2 = geometry->elem[iElem]->GetNode(2);   Coord_2 = geometry->node[Point_2]->GetCoord();
      Point_3 = geometry->elem[iElem]->GetNode(3);  Coord_3 = geometry->node[Point_3]->GetCoord();
      
      numerics->SetCoord(Coord_0, Coord_1, Coord_2, Coord_3);
      numerics->ComputeResidual(StiffMatrix_Elem, config);
      
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0]; StiffMatrixSpace.AddBlock(Point_0, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][1]; StiffMatrixSpace.AddBlock(Point_0, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2]; StiffMatrixSpace.AddBlock(Point_0, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][3]; StiffMatrixSpace.AddBlock(Point_0, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][0]; StiffMatrixSpace.AddBlock(Point_1, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][1]; StiffMatrixSpace.AddBlock(Point_1, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][2]; StiffMatrixSpace.AddBlock(Point_1, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[1][3]; StiffMatrixSpace.AddBlock(Point_1, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0]; StiffMatrixSpace.AddBlock(Point_2, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][1]; StiffMatrixSpace.AddBlock(Point_2, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2]; StiffMatrixSpace.AddBlock(Point_2, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][3]; StiffMatrixSpace.AddBlock(Point_2, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][0]; StiffMatrixSpace.AddBlock(Point_3, Point_0, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][1]; StiffMatrixSpace.AddBlock(Point_3, Point_1, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][2]; StiffMatrixSpace.AddBlock(Point_3, Point_2, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node);
      StiffMatrix_Node[0][0] = StiffMatrix_Elem[3][3]; StiffMatrixSpace.AddBlock(Point_3, Point_3, StiffMatrix_Node); Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node);
      
    }
  }
  
  if (config->GetUnsteady_Simulation() != STEADY) {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      total_index = iPoint*nVar;
      LinSysSol[total_index] = node[iPoint]->GetSolution(0);
      LinSysAux[total_index] = 0.0;
    }
    
    StiffMatrixSpace.MatrixVectorProduct(LinSysSol, LinSysAux);
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      total_index = iPoint*nVar;
      Residual[0] = LinSysAux[total_index];
      LinSysRes.SubtractBlock(iPoint, Residual);
    }
  }
  
}

void CHeatSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

void CHeatSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex, total_index;
  su2double Twall;
  
  /*--- Identify the boundary ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature ---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag);
  
  /*--- Set the solution at the boundary nodes and zero the residual ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    Solution[0] = Twall;
    
    node[iPoint]->SetSolution(Solution);
    node[iPoint]->SetSolution_Old(Solution);

    /*--- Unsteady solution, the equation is solved in terms of increments ---*/
    if (config->GetUnsteady_Simulation() != STEADY) Residual[0] = 0.0;
    
    LinSysRes.SetBlock(iPoint, Residual);
    LinSysSol.SetBlock(iPoint, Residual);

    total_index = iPoint*nVar;
    Jacobian.DeleteValsRowi(total_index);
    
  }
  
}


void CHeatSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                       unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3 = 0;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, d[3] = {0.0,0.0,0.0}, Area_Local = 0.0, Volume_Local = 0.0, Time_Num;
  su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
  unsigned short iDim, iVar;
  su2double TimeJac = 0.0;
  
  /*--- Numerical time step (this system is uncoditional stable... a very big number can be used) ---*/
  Time_Num = config->GetDelta_UnstTimeND();
  
  /*--- Loop through elements to compute contributions from the matrix
   blocks involving time. These contributions are also added to the
   Jacobian w/ the time step. Spatial source terms are also computed. ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    /*--- Get node numbers and their coordinate vectors ---*/
    
    Point_0 = geometry->elem[iElem]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
    Point_1 = geometry->elem[iElem]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
    Point_2 = geometry->elem[iElem]->GetNode(2);  Coord_2 = geometry->node[Point_2]->GetCoord();
    
    /*--- Compute area and volume ---*/

    if (nDim == 2) {
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = Coord_0[iDim]-Coord_2[iDim];
        b[iDim] = Coord_1[iDim]-Coord_2[iDim];
      }
      Area_Local = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
    }
    else {
      Point_3 = geometry->elem[iElem]->GetNode(3);
      Coord_3 = geometry->node[Point_3]->GetCoord();
      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = Coord_0[iDim]-Coord_2[iDim];
        b[iDim] = Coord_1[iDim]-Coord_2[iDim];
        c[iDim] = Coord_3[iDim]-Coord_2[iDim];
      }
      d[0] = a[1]*b[2]-a[2]*b[1]; d[1] = -(a[0]*b[2]-a[2]*b[0]); d[2] = a[0]*b[1]-a[1]*b[0];
      Volume_Local = fabs(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])/6.0;
    }
    
    /*--- Block contributions to the Jacobian (includes time step) ---*/
        
    if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) TimeJac = 1.0/Time_Num;
    if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) TimeJac = 3.0/(2.0*Time_Num);
      
    if (nDim == 2) { StiffMatrix_Node[0][0] = (2.0/12.0)*(Area_Local*TimeJac); }
    else { StiffMatrix_Node[0][0] = (2.0/20.0)*(Volume_Local*TimeJac); }
    
    Jacobian.AddBlock(Point_0, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_0, StiffMatrix_Node);
    Jacobian.AddBlock(Point_1, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_1, StiffMatrix_Node);
    Jacobian.AddBlock(Point_2, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node);
    if (nDim == 3) { Jacobian.AddBlock(Point_3, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_2, StiffMatrix_Node); }
      
    if (nDim == 2) { StiffMatrix_Node[0][0] = (1.0/12.0)*(Area_Local*TimeJac); }
    else { StiffMatrix_Node[0][0] = (1.0/20.0)*(Volume_Local*TimeJac); }
    
    Jacobian.AddBlock(Point_0, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_1, StiffMatrix_Node);
    Jacobian.AddBlock(Point_0, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_2, StiffMatrix_Node);
    Jacobian.AddBlock(Point_1, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_0, StiffMatrix_Node);
    Jacobian.AddBlock(Point_1, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_2, StiffMatrix_Node);
    Jacobian.AddBlock(Point_2, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_0, StiffMatrix_Node);
    Jacobian.AddBlock(Point_2, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_1, StiffMatrix_Node);
    if (nDim == 3) {
      Jacobian.AddBlock(Point_0, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_0, Point_3, StiffMatrix_Node);
      Jacobian.AddBlock(Point_1, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_1, Point_3, StiffMatrix_Node);
      Jacobian.AddBlock(Point_2, Point_3, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_2, Point_3, StiffMatrix_Node);
      Jacobian.AddBlock(Point_3, Point_0, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_0, StiffMatrix_Node);
      Jacobian.AddBlock(Point_3, Point_1, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_1, StiffMatrix_Node);
      Jacobian.AddBlock(Point_3, Point_2, StiffMatrix_Node); StiffMatrixTime.AddBlock(Point_3, Point_2, StiffMatrix_Node);
    }
    
  }
  
  unsigned long iPoint, total_index;
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  
  /*--- loop over points ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Solution at time n-1, n and n+1 ---*/
    
    U_time_nM1 = node[iPoint]->GetSolution_time_n1();
    U_time_n   = node[iPoint]->GetSolution_time_n();
    U_time_nP1 = node[iPoint]->GetSolution();
    
    /*--- Compute residual ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        LinSysSol[total_index] = ( U_time_nP1[iVar] - U_time_n[iVar] );
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        LinSysSol[total_index] = ( U_time_nP1[iVar] - (4.0/3.0)*U_time_n[iVar] + (1.0/3.0)*U_time_nM1[iVar] );
    }
  }
  
  /*--- Contribution to the residual ---*/
  
  StiffMatrixTime.MatrixVectorProduct(LinSysSol, LinSysAux);
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    total_index = iPoint*nVar;
    Residual[0] = LinSysAux[total_index];
    LinSysRes.SubtractBlock(iPoint, Residual);
  }
  
}

void CHeatSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysSol[total_index] = 0.0;
    }
    
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      if (config->GetUnsteady_Simulation() == STEADY) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
      else node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
    }
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*---  Compute the residual Ax-f ---*/
  
  Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Compute the residual ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
      AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

CHeatSolverFVM::CHeatSolverFVM(void) : CSolver() {

  ConjugateVar = NULL;
}

CHeatSolverFVM::CHeatSolverFVM(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  unsigned short iVar, iDim, nLineLets, iMarker;
  unsigned long iPoint, iVertex;

  int rank = MASTER_NODE;

  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool heat_equation      = config->GetKind_Solver() == HEAT_EQUATION_FVM;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Dimension of the problem --> temperature is the only conservative variable ---*/

  nVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  nMarker = config->GetnMarker_All();

  CurrentMesh = iMesh;
  /*--- Define some auxiliar vector related with the residual ---*/

  Residual      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Residual_Max  = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Res_Conv      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliar vector related with the solution ---*/

  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the primitive flow solution ---*/

  Primitive_Flow_i = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) Primitive_Flow_i[iVar] = 0.0;
  Primitive_Flow_j = new su2double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) Primitive_Flow_j[iVar] = 0.0;

  /*--- Jacobians and vector structures for implicit computations ---*/

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

  /*--- Initialization of the structure of the whole Jacobian ---*/

  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (heat equation) MG level: " << iMesh << "." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  if (config->GetExtraOutput()) {
    if (nDim == 2) { nOutputVariables = 13; }
    else if (nDim == 3) { nOutputVariables = 19; }
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
    OutputHeadingNames = new string[nOutputVariables];
  }

  /*--- Computation of gradients by least squares ---*/

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
  }

  Heat_Flux = new su2double[nMarker];
  AvgTemperature = new su2double[nMarker];
  Surface_Areas = new su2double[config->GetnMarker_HeatFlux()];

  for(iMarker = 0; iMarker < nMarker; iMarker++) {
    Heat_Flux[iMarker] = 0.0;
    AvgTemperature[iMarker] = 0.0;
  }
  for(iMarker = 0; iMarker < config->GetnMarker_HeatFlux(); iMarker++) {
    Surface_Areas[iMarker] = 0.0;
  }

  Set_Heatflux_Areas(geometry, config);

  /*--- Non-dimensionalization of heat equation */

  su2double Temperature_FreeStream = config->GetInc_Temperature_Init();
  config->SetTemperature_FreeStream(Temperature_FreeStream);
  su2double Temperature_Ref = 0.0;

  if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
    Temperature_Ref = 1.0;
  }
  else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
    Temperature_Ref = Temperature_FreeStream;
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    Temperature_Ref = config->GetInc_Temperature_Ref();
  }
  config->SetTemperature_Ref(Temperature_Ref);

  config->SetTemperature_FreeStreamND(config->GetTemperature_FreeStream()/config->GetTemperature_Ref());
  if (rank == MASTER_NODE) {
    cout << "Weakly coupled heat solver's freestream temperature: " << config->GetTemperature_FreeStreamND() << endl;
  }

  su2double Temperature_Solid_Freestream_ND = config->GetTemperature_Freestream_Solid()/config->GetTemperature_Ref();
  if (heat_equation && (rank == MASTER_NODE)) {
    cout << "Heat solver freestream temperature in case for solids: " << Temperature_Solid_Freestream_ND << endl;
  }

  /*--- Store the value of the temperature and the heat flux density at the boundaries,
   used for IO with a donor cell ---*/
  unsigned short nConjVariables = 4;

  ConjugateVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    ConjugateVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      ConjugateVar[iMarker][iVertex] = new su2double [nConjVariables];
      for (iVar = 0; iVar < nConjVariables ; iVar++) {
        ConjugateVar[iMarker][iVertex][iVar] = 0.0;
      }
      ConjugateVar[iMarker][iVertex][0] = config->GetTemperature_FreeStreamND();
    }
  }

  /*--- If the heat solver runs stand-alone, we have to set the reference values ---*/
  if(heat_equation) {
    su2double rho_cp = config->GetDensity_Solid()*config->GetSpecific_Heat_Cp_Solid();
    su2double thermal_diffusivity_solid = config->GetThermalConductivity_Solid() / rho_cp;
    config->SetThermalDiffusivity_Solid(thermal_diffusivity_solid);
  }

    for (iPoint = 0; iPoint < nPoint; iPoint++)
      if (flow)
        node[iPoint] = new CHeatFVMVariable(config->GetTemperature_FreeStreamND(), nDim, nVar, config);
      else
        node[iPoint] = new CHeatFVMVariable(Temperature_Solid_Freestream_ND, nDim, nVar, config);

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
}

CHeatSolverFVM::~CHeatSolverFVM(void) { }


void CHeatSolverFVM::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;
  bool center = (config->GetKind_ConvNumScheme_Heat() == SPACE_CENTERED);

  if (center) {
    SetUndivided_Laplacian(geometry, config);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CHeatSolverFVM::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) { }

void CHeatSolverFVM::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool heat_equation = config->GetKind_Solver() == HEAT_EQUATION_FVM;

  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = 0;

  if (flow) {

    if (config->GetKind_Turb_Model() == SA || config->GetKind_Turb_Model() == SA_NEG) {
      if (nDim == 2) skipVars += 6;
      if (nDim == 3) skipVars += 8;
    }
    else if (config->GetKind_Turb_Model() == SST ) {
      if (nDim == 2) skipVars += 7;
      if (nDim == 3) skipVars += 9;
    }
    else {
      if (nDim == 2) skipVars += 5;
      if (nDim == 3) skipVars += 7;
    }
  }
  else if (heat_equation) {

    if (nDim == 2) skipVars += 2;
    if (nDim == 3) skipVars += 3;
  }
  else {
    cout << "WARNING: Finite volume heat solver's restart routine could not load data." << endl;
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    if (rank == MASTER_NODE) {
      cout << endl << "The solution file " << restart_filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][HEAT_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][HEAT_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][HEAT_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][HEAT_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][HEAT_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][HEAT_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_HEAT_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}


void CHeatSolverFVM::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, jPoint, iEdge;
  su2double *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;

  Diff = new su2double[nVar];

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Solution differences ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);

    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

    /*--- Both points inside the domain, or both in the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

  }

  /*--- MPI parallelization ---*/

  Set_MPI_Undivided_Laplacian(geometry, config);

  delete [] Diff;

}

void CHeatSolverFVM::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Only copy conserved variables - no transformation necessary. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;

    }

  }

}

void CHeatSolverFVM::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  su2double *V_i, *V_j, Temp_i, Temp_j;
  unsigned long iEdge, iPoint, jPoint;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  if(flow) {

    nVarFlow = solver_container[FLOW_SOL]->GetnVar();

      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

        /*--- Points in edge ---*/
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

        /*--- Primitive variables w/o reconstruction ---*/
        V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
        V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();

        Temp_i = node[iPoint]->GetSolution(0);
        Temp_j = node[jPoint]->GetSolution(0);

        numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
        numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());

        numerics->SetPrimitive(V_i, V_j);
        numerics->SetTemperature(Temp_i, Temp_j);

        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        LinSysRes.AddBlock(iPoint, Residual);
        LinSysRes.SubtractBlock(jPoint, Residual);

        /*--- Implicit part ---*/

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
      }
  }
}

void CHeatSolverFVM::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {

  su2double *V_i, *V_j, Temp_i, Temp_i_Corrected, Temp_j, Temp_j_Corrected, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j,
          **Temp_i_Grad, **Temp_j_Grad, Project_Temp_i_Grad, Project_Temp_j_Grad, Non_Physical = 1.0;
  unsigned short iDim, iVar;
  unsigned long iEdge, iPoint, jPoint;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool muscl = (config->GetMUSCL_Heat());

  if(flow) {

    nVarFlow = solver_container[FLOW_SOL]->GetnVar();

      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

        /*--- Points in edge ---*/
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

        /*--- Primitive variables w/o reconstruction ---*/
        V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
        V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();

        Temp_i = node[iPoint]->GetSolution(0);
        Temp_j = node[jPoint]->GetSolution(0);

        /* Second order reconstruction */
        if (muscl) {

            for (iDim = 0; iDim < nDim; iDim++) {
              Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
              Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
            }

            Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
            Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
            Temp_i_Grad = node[iPoint]->GetGradient();
            Temp_j_Grad = node[jPoint]->GetGradient();

            /*Loop to correct the flow variables*/
            for (iVar = 0; iVar < nVarFlow; iVar++) {

              /*Apply the Gradient to get the right temperature value on the edge */
              Project_Grad_i = 0.0; Project_Grad_j = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                  Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
                  Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
              }

              Primitive_Flow_i[iVar] = V_i[iVar] + Project_Grad_i;
              Primitive_Flow_j[iVar] = V_j[iVar] + Project_Grad_j;
            }

            /* Correct the temperature variables */
            Project_Temp_i_Grad = 0.0; Project_Temp_j_Grad = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
                Project_Temp_i_Grad += Vector_i[iDim]*Temp_i_Grad[0][iDim]*Non_Physical;
                Project_Temp_j_Grad += Vector_j[iDim]*Temp_j_Grad[0][iDim]*Non_Physical;
            }

            Temp_i_Corrected = Temp_i + Project_Temp_i_Grad;
            Temp_j_Corrected = Temp_j + Project_Temp_j_Grad;

            numerics->SetPrimitive(Primitive_Flow_i, Primitive_Flow_j);
            numerics->SetTemperature(Temp_i_Corrected, Temp_j_Corrected);
        }

        else {

          numerics->SetPrimitive(V_i, V_j);
          numerics->SetTemperature(Temp_i, Temp_j);
        }

        numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        LinSysRes.AddBlock(iPoint, Residual);
        LinSysRes.SubtractBlock(jPoint, Residual);

        /*--- Implicit part ---*/

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
        }
  }

}

void CHeatSolverFVM::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  su2double laminar_viscosity, Prandtl_Lam, Prandtl_Turb, eddy_viscosity_i, eddy_viscosity_j,
      thermal_diffusivity_i, thermal_diffusivity_j, Temp_i, Temp_j, **Temp_i_Grad, **Temp_j_Grad;
  unsigned long iEdge, iPoint, jPoint;

  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  laminar_viscosity = config->GetMu_ConstantND();
  Prandtl_Lam = config->GetPrandtl_Lam();
  Prandtl_Turb = config->GetPrandtl_Turb();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    Temp_i_Grad = node[iPoint]->GetGradient();
    Temp_j_Grad = node[jPoint]->GetGradient();
    numerics->SetConsVarGradient(Temp_i_Grad, Temp_j_Grad);

    /*--- Primitive variables w/o reconstruction ---*/
    Temp_i = node[iPoint]->GetSolution(0);
    Temp_j = node[jPoint]->GetSolution(0);
    numerics->SetTemperature(Temp_i, Temp_j);

    /*--- Eddy viscosity to compute thermal conductivity ---*/
    if (flow) {
      eddy_viscosity_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
      eddy_viscosity_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity();

      thermal_diffusivity_i = (laminar_viscosity/Prandtl_Lam) + (eddy_viscosity_i/Prandtl_Turb);
      thermal_diffusivity_j = (laminar_viscosity/Prandtl_Lam) + (eddy_viscosity_j/Prandtl_Turb);
    }
    else {
      thermal_diffusivity_i = config->GetThermalDiffusivity_Solid();
      thermal_diffusivity_j = config->GetThermalDiffusivity_Solid();
    }

    numerics->SetThermalDiffusivity(thermal_diffusivity_i,thermal_diffusivity_j);

    /*--- Compute residual, and Jacobians ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobians ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
  }
}

void CHeatSolverFVM::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }

void CHeatSolverFVM::Set_Heatflux_Areas(CGeometry *geometry, CConfig *config) {

  unsigned short iMarker, iMarker_HeatFlux, Monitoring, iDim;
  unsigned long iPoint, iVertex;
  string HeatFlux_Tag, Marker_Tag;

  su2double *Local_Surface_Areas, Local_HeatFlux_Areas_Monitor, Area, *Normal;
  Local_Surface_Areas = new su2double[config->GetnMarker_HeatFlux()];

  for ( iMarker_HeatFlux = 0; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {
    Local_Surface_Areas[iMarker_HeatFlux] = 0.0;
  }
  Local_HeatFlux_Areas_Monitor = 0.0;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    for ( iMarker_HeatFlux = 0; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {

      HeatFlux_Tag = config->GetMarker_HeatFlux_TagBound(iMarker_HeatFlux);
      Marker_Tag = config->GetMarker_All_TagBound(iMarker);

      if (Marker_Tag == HeatFlux_Tag) {

        Local_Surface_Areas[iMarker_HeatFlux] = 0.0;

        for( iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if(geometry->node[iPoint]->GetDomain()) {

            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

            Local_Surface_Areas[iMarker_HeatFlux] += Area;

            if(Monitoring == YES) {
              Local_HeatFlux_Areas_Monitor += Area;
            }
          }
        }
      }
    }
  }
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(Local_Surface_Areas, Surface_Areas, config->GetnMarker_HeatFlux(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_HeatFlux_Areas_Monitor, &Total_HeatFlux_Areas_Monitor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for( iMarker_HeatFlux = 0; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {
      Surface_Areas[iMarker_HeatFlux] = Local_Surface_Areas[iMarker_HeatFlux];
    }
    Total_HeatFlux_Areas_Monitor = Local_HeatFlux_Areas_Monitor;
#endif

  Total_HeatFlux_Areas = 0.0;
  for( iMarker_HeatFlux = 0; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {
    Total_HeatFlux_Areas += Surface_Areas[iMarker_HeatFlux];
  }

  delete[] Local_Surface_Areas;
}

void CHeatSolverFVM::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {

  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iDim;
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, laminar_viscosity, thermal_diffusivity, Twall, dTdn, Prandtl_Lam;
  //su2double Prandtl_Turb;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  Prandtl_Lam = config->GetPrandtl_Lam();
//  Prandtl_Turb = config->GetPrandtl_Turb();
  laminar_viscosity = config->GetMu_ConstantND();
  //Prandtl_Turb = config->GetPrandtl_Turb();
  //laminar_viscosity = config->GetViscosity_FreeStreamND(); // TDE check for consistency for CHT

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);

        Coord_i = geometry->node[iPoint]->GetCoord();
        Coord_j = geometry->node[Point_Normal]->GetCoord();
        dist_ij = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
        dist_ij = sqrt(dist_ij);

        dTdn = -(node[Point_Normal]->GetSolution(0) - Twall)/dist_ij;

        if(flow) {
          thermal_diffusivity = laminar_viscosity/Prandtl_Lam;
        }
        else
          thermal_diffusivity = config->GetThermalDiffusivity_Solid();

        Res_Visc[0] = thermal_diffusivity*dTdn*Area;

        if(implicit) {

          Jacobian_i[0][0] = -thermal_diffusivity/dist_ij * Area;
        }

        LinSysRes.SubtractBlock(iPoint, Res_Visc);
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    }
  }
}

void CHeatSolverFVM::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                                     unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double Wall_HeatFlux, Area, *Normal;

  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  if(config->GetIntegrated_HeatFlux()) {

    unsigned short iMarker_HeatFlux;
    string HeatFlux_Tag, Marker_Tag;

    // Find out which heat flux wall to get the right surface area

    for ( iMarker_HeatFlux = 0; iMarker_HeatFlux < config->GetnMarker_HeatFlux(); iMarker_HeatFlux++ ) {

      HeatFlux_Tag = config->GetMarker_HeatFlux_TagBound(iMarker_HeatFlux);
      Marker_Tag = config->GetMarker_All_TagBound(val_marker);

      if (Marker_Tag == HeatFlux_Tag) {
        Wall_HeatFlux = Wall_HeatFlux / Surface_Areas[iMarker_HeatFlux];
      }
    }
  }

  if(flow) {
    Wall_HeatFlux = Wall_HeatFlux/(config->GetViscosity_Ref()*config->GetSpecific_Heat_Cp()*config->GetTemperature_Ref());
  }
  else {
    Wall_HeatFlux = Wall_HeatFlux/(config->GetDensity_Solid()*config->GetSpecific_Heat_Cp_Solid()*config->GetTemperature_Ref());
  }

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      Res_Visc[0] = 0.0;

      Res_Visc[0] = Wall_HeatFlux * Area;

      /*--- Viscous contribution to the residual at the wall ---*/

      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }

  }
}

void CHeatSolverFVM::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *Flow_Dir,  Vel_Mag;
  su2double *V_inlet, *V_domain;

  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  bool viscous              = config->GetViscous();
  bool grid_movement        = config->GetGrid_Movement();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];

  su2double *Coord_i, *Coord_j, Area, dist_ij, laminar_viscosity, thermal_diffusivity, Twall, dTdn, Prandtl_Lam;
  //su2double Prandtl_Turb;
  Prandtl_Lam = config->GetPrandtl_Lam();
//  Prandtl_Turb = config->GetPrandtl_Turb();
  laminar_viscosity = config->GetMu_ConstantND();
  //laminar_viscosity = config->GetViscosity_FreeStreamND(); //TDE check for consistency with CHT

  Twall = config->GetTemperature_FreeStreamND();

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      if(flow) {

        /*--- Normal vector for this vertex (negate for outward convention) ---*/

        conv_numerics->SetNormal(Normal);

        /*--- Retrieve solution at this boundary node ---*/

        V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

        /*--- Retrieve the specified velocity for the inlet. ---*/

        Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag)/config->GetVelocity_Ref();
        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

        V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];

        conv_numerics->SetPrimitive(V_domain, V_inlet);

        if (grid_movement)
          conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

        conv_numerics->SetTemperature(node[iPoint]->GetSolution(0), config->GetInlet_Ttotal(Marker_Tag)/config->GetTemperature_Ref());

        /*--- Compute the residual using an upwind scheme ---*/

        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        /*--- Update residual value ---*/

        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/

        if (implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }

      /*--- Viscous contribution ---*/

      if (viscous) {

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);

        Coord_i = geometry->node[iPoint]->GetCoord();
        Coord_j = geometry->node[Point_Normal]->GetCoord();
        dist_ij = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
        dist_ij = sqrt(dist_ij);

        dTdn = -(node[Point_Normal]->GetSolution(0) - Twall)/dist_ij;

        thermal_diffusivity = laminar_viscosity/Prandtl_Lam;

        Res_Visc[0] = thermal_diffusivity*dTdn*Area;

        if(implicit) {

          Jacobian_i[0][0] = -thermal_diffusivity/dist_ij * Area;
        }
        /*--- Viscous contribution to the residual at the wall ---*/

        LinSysRes.SubtractBlock(iPoint, Res_Visc);
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CHeatSolverFVM::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain;

  bool flow                 = (config->GetKind_Solver() != HEAT_EQUATION);
  bool grid_movement        = config->GetGrid_Movement();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  su2double *Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Normal vector for this vertex (negate for outward convention) ---*/

        geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
        for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

        if(flow) {
            conv_numerics->SetNormal(Normal);

            /*--- Retrieve solution at this boundary node ---*/

            V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

            /*--- Retrieve the specified velocity for the inlet. ---*/

            V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
            for (iDim = 0; iDim < nDim; iDim++)
              V_outlet[iDim+1] = solver_container[FLOW_SOL]->node[Point_Normal]->GetPrimitive(iDim+1);

            conv_numerics->SetPrimitive(V_domain, V_outlet);

            if (grid_movement)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

            conv_numerics->SetTemperature(node[iPoint]->GetSolution(0), node[Point_Normal]->GetSolution(0));

            /*--- Compute the residual using an upwind scheme ---*/

            conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

            /*--- Update residual value ---*/

            LinSysRes.AddBlock(iPoint, Residual);

            /*--- Jacobian contribution for implicit integration ---*/

            if (implicit)
              Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        }

        // viscous contribution is still missing...
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CHeatSolverFVM::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint, total_index;
  unsigned short iDim, iVar, iMarker;

  su2double Area, rho_cp_solid,
      Temperature_Ref, Tinterface, T_Conjugate, Tnormal_Conjugate, Conductance, HeatFluxDensity, HeatFluxValue;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

  su2double *Normal = new su2double[nDim];

  Temperature_Ref       = config->GetTemperature_Ref();
  rho_cp_solid          = config->GetDensity_Solid()*config->GetSpecific_Heat_Cp_Solid();

  if (flow) {

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);

            T_Conjugate = GetConjugateHeatVariable(iMarker, iVertex, 0)/Temperature_Ref;

            node[iPoint]->SetSolution_Old(&T_Conjugate);
            LinSysRes.SetBlock_Zero(iPoint, 0);
            node[iPoint]->SetRes_TruncErrorZero();

            if (implicit) {
              for (iVar = 0; iVar < nVar; iVar++) {
                total_index = iPoint*nVar+iVar;
                Jacobian.DeleteValsRowi(total_index);
              }
            }
          }
        }
      }
    }
  }
  else {

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);

            Tinterface          = node[iPoint]->GetSolution(0);
            Tnormal_Conjugate   = GetConjugateHeatVariable(iMarker, iVertex, 3)/Temperature_Ref;
            Conductance         = GetConjugateHeatVariable(iMarker, iVertex, 2)/rho_cp_solid;

            HeatFluxDensity     = Conductance*(Tinterface - Tnormal_Conjugate);

            HeatFluxValue       = HeatFluxDensity * Area;

            Res_Visc[0] = -HeatFluxValue;
            LinSysRes.SubtractBlock(iPoint, Res_Visc);

            if (implicit) {

              Jacobian_i[0][0] = Conductance*Area;
              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            }
          }
        }
      }
    }
  }
}

void CHeatSolverFVM::Heat_Fluxes(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iDim;
  su2double *Coord, *Coord_Normal, *Normal, Area, dist, Twall, dTdn, cp_fluid, rho_cp_solid,
      thermal_conductivity, thermal_diffusivity;
  string Marker_Tag, HeatFlux_Tag;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));

#ifdef HAVE_MPI
  su2double MyAllBound_HeatFlux, MyAllBound_AvgTemperature;
#endif

  cp_fluid = config->GetSpecific_Heat_Cp();
  rho_cp_solid = config->GetSpecific_Heat_Cp_Solid()*config->GetDensity_Solid();

  AllBound_HeatFlux = 0.0;
  AllBound_AvgTemperature = 0.0;

  for ( iMarker = 0; iMarker < nMarker; iMarker++ ) {

    AvgTemperature[iMarker] = 0.0;

    Boundary = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    Heat_Flux[iMarker] = 0.0;

    if ( Boundary == ISOTHERMAL ) {

      Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

      for( iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if(geometry->node[iPoint]->GetDomain()) {

          iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          Coord = geometry->node[iPoint]->GetCoord();
          Coord_Normal = geometry->node[iPointNormal]->GetCoord();

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

          dist = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) dist += (Coord_Normal[iDim]-Coord[iDim])*(Coord_Normal[iDim]-Coord[iDim]);
          dist = sqrt(dist);

          dTdn = (Twall - node[iPointNormal]->GetSolution(0))/dist;

          if(flow) {
            thermal_diffusivity = config->GetViscosity_FreeStreamND()/config->GetPrandtl_Lam();
            thermal_conductivity = thermal_diffusivity*config->GetViscosity_Ref()*cp_fluid;
          }
          else {
            thermal_conductivity = config->GetThermalDiffusivity_Solid()*rho_cp_solid;
          }

          Heat_Flux[iMarker] += thermal_conductivity*dTdn*config->GetTemperature_Ref()*Area;
        }
      }
    }
    else if ( Boundary == CHT_WALL_INTERFACE || Boundary == HEAT_FLUX ) {

      for( iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++ ) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if(geometry->node[iPoint]->GetDomain()) {

          iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          Twall = node[iPoint]->GetSolution(0);

          Coord = geometry->node[iPoint]->GetCoord();
          Coord_Normal = geometry->node[iPointNormal]->GetCoord();

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

          dist = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) dist += (Coord_Normal[iDim]-Coord[iDim])*(Coord_Normal[iDim]-Coord[iDim]);
          dist = sqrt(dist);

          dTdn = (Twall - node[iPointNormal]->GetSolution(0))/dist;

          if(flow) {
            thermal_diffusivity = config->GetViscosity_FreeStreamND()/config->GetPrandtl_Lam();
            thermal_conductivity = thermal_diffusivity*config->GetViscosity_Ref()*cp_fluid;
          }
          else {
            thermal_conductivity = config->GetThermalDiffusivity_Solid()*rho_cp_solid;
          }

          Heat_Flux[iMarker] += thermal_conductivity*dTdn*config->GetTemperature_Ref()*Area;

          /*--- We do only aim to compute averaged temperatures on the (interesting) heat flux walls ---*/
          if ( Boundary == HEAT_FLUX ) {

            AvgTemperature[iMarker] += Twall*config->GetTemperature_Ref()*Area;
          }

        }
      }
    }

    if (Monitoring == YES) {

    AllBound_HeatFlux += Heat_Flux[iMarker];
    AllBound_AvgTemperature += AvgTemperature[iMarker];
    }
  }

#ifdef HAVE_MPI
  MyAllBound_HeatFlux = AllBound_HeatFlux;
  MyAllBound_AvgTemperature = AllBound_AvgTemperature;
  SU2_MPI::Allreduce(&MyAllBound_HeatFlux, &AllBound_HeatFlux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_AvgTemperature, &AllBound_AvgTemperature, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (Total_HeatFlux_Areas_Monitor != 0.0) {
    Total_AvgTemperature = AllBound_AvgTemperature/Total_HeatFlux_Areas_Monitor;
  }
  else {
    Total_AvgTemperature = 0.0;
  }


  Total_HeatFlux = AllBound_HeatFlux;
}

void CHeatSolverFVM::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                               unsigned short iMesh, unsigned long Iteration) {

  unsigned short iDim, iMarker;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  su2double *Normal, Area, Vol, laminar_viscosity, eddy_viscosity, thermal_diffusivity, Prandtl_Lam, Prandtl_Turb, Mean_ProjVel, Mean_BetaInc2, Mean_DensityInc, Mean_SoundSpeed, Lambda;
  su2double Global_Delta_Time, Global_Delta_UnstTimeND = 0.0, Local_Delta_Time = 0.0, Local_Delta_Time_Inv, Local_Delta_Time_Visc, CFL_Reduction, K_v = 0.25;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  laminar_viscosity = config->GetMu_ConstantND();
  Prandtl_Lam = config->GetPrandtl_Lam();
  Prandtl_Turb = config->GetPrandtl_Turb();

  thermal_diffusivity = config->GetThermalDiffusivity_Solid();

  /*--- Compute spectral radius based on thermal conductivity ---*/

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  CFL_Reduction = config->GetCFLRedCoeff_Turb();

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- get the edge's normal vector to compute the edge's area ---*/
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Inviscid contribution ---*/

    if (flow) {
      Mean_ProjVel = 0.5 * (solver_container[FLOW_SOL]->node[iPoint]->GetProjVel(Normal) + solver_container[FLOW_SOL]->node[jPoint]->GetProjVel(Normal));
      Mean_BetaInc2 = 0.5 * (solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2() + solver_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
      Mean_DensityInc = 0.5 * (solver_container[FLOW_SOL]->node[iPoint]->GetDensity() + solver_container[FLOW_SOL]->node[jPoint]->GetDensity());
      Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);

      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    }

    /*--- Viscous contribution ---*/

    thermal_diffusivity = config->GetThermalDiffusivity_Solid();
    if(flow) {
      eddy_viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
      thermal_diffusivity = laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb;
    }

    Lambda = thermal_diffusivity*Area*Area;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

      /*--- Inviscid contribution ---*/

      if (flow) {
        Mean_ProjVel = solver_container[FLOW_SOL]->node[iPoint]->GetProjVel(Normal);
        Mean_BetaInc2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
        Mean_DensityInc = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
        Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);

        Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
        if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }

      /*--- Viscous contribution ---*/

      thermal_diffusivity = config->GetThermalDiffusivity_Solid();
      if(flow) {
        eddy_viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        thermal_diffusivity = laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb;
      }

      Lambda = thermal_diffusivity*Area*Area;
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);

    }
  }

  /*--- Each element uses their own speed, steady state simulation ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();

    if (Vol != 0.0) {

      if(flow) {
        Local_Delta_Time_Inv = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
        Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      }
      else {
        Local_Delta_Time_Inv = config->GetMax_DeltaTime();
        Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
        //Local_Delta_Time_Visc = 100.0*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      }

      /*--- Time step setting method ---*/

      if (config->GetKind_TimeStep_Heat() == MINIMUM)
        Local_Delta_Time = min(Local_Delta_Time_Inv, Local_Delta_Time_Visc);
      else if (config->GetKind_TimeStep_Heat() == CONVECTIVE)
        Local_Delta_Time = Local_Delta_Time_Inv;
      else if (config->GetKind_TimeStep_Heat() == VISCOUS)
        Local_Delta_Time = Local_Delta_Time_Visc;
      else if (config->GetKind_TimeStep_Heat() == BYFLOW)
        Local_Delta_Time = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();

      /*--- Min-Max-Logic ---*/

      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();

      node[iPoint]->SetDelta_Time(CFL_Reduction*Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }

  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        cout << "Using unsteady time: " << config->GetDelta_UnstTimeND() << endl;
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
  }
}

void CHeatSolverFVM::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}


void CHeatSolverFVM::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, *local_Res_TruncError;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));


  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = node[iPoint]->GetResTruncError();

    /*--- Read the volume ---*/

    Vol = geometry->node[iPoint]->GetVolume();

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    if (node[iPoint]->GetDelta_Time() != 0.0) {

      if(flow) {
        Delta = Vol / node[iPoint]->GetDelta_Time();
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
      else {
        Delta = Vol / node[iPoint]->GetDelta_Time();
        Jacobian.AddVal2Diag(iPoint, Delta);
      }

    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/

  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
    }
  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CHeatSolverFVM::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++) {
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        }
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }

}

void CHeatSolverFVM::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }
}

void CHeatSolverFVM::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;

  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];

#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];

        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }

        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;

    }

  }

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;

}

void CHeatSolverFVM::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution_Fine, *Solution;

  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/

  if (restart && (ExtIter == 0)) {

    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][HEAT_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][HEAT_SOL]->node[iPoint]->SetSolution(Solution);
      }
      solver_container[iMesh][HEAT_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
  }

  /*--- The value of the solution for the first iteration of the dual time ---*/

  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][HEAT_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][HEAT_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }

    if ((restart && (long)ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/

      solver_container[MESH_0][HEAT_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), true);

      /*--- Push back this new solution to time level N. ---*/

      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][HEAT_SOL]->node[iPoint]->Set_Solution_time_n();
        }
      }
    }
  }
}

void CHeatSolverFVM::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Local variables ---*/

  unsigned short iVar, jVar;
  unsigned long iPoint;

  su2double *U_time_n, *U_time_nP1, *U_time_nM1;
  su2double Volume_nP1, TimeStep;

  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!grid_movement) {

    /*--- Loop over all nodes (excluding halos) ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->node[iPoint]->GetVolume();

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/

      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
}
