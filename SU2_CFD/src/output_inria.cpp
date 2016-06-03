#include "../include/output_structure.hpp"

enum BCVAR  { bcMach, bcTemp, bcPres, bcDens };

void COutput::SetInriaRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  //ofstream restart_file;
  string filename;
  
	unsigned long OutSol,i, npoin = geometry->GetGlobal_nPointDomain();
	int VarTyp[GmfMaxTyp];
	double bufDbl[GmfMaxTyp];
	char OutNam[1024], BasNam[1024];
	char *ptr=NULL;
	
	
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
	
	/*--- Get output name *.solb ---*/
	
	strcpy(BasNam, filename.c_str());
	ptr = strstr(BasNam,".dat");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	ptr = strstr(BasNam,".solb");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	sprintf(OutNam, "%s.solb", BasNam);
	
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iZone));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if ((fem) && (config->GetWrt_Dynamic())) {
	filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  /*--- Open the restart file and write the solution. ---*/
	
	OutSol = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutSol ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  /*--- Write the restart file ---*/

	for (iVar = 0; iVar < nVar_Total; iVar++) {
		VarTyp[iVar] = GmfSca;
	}
	
	npoin = geometry->GetGlobal_nPointDomain();
	
	//printf("SET KWD : %d %d %d %d\n", OutSol, GmfSolAtVertices, npoin, nVar_Total);
	//printf("VarTyp = ");
	//for (iVar = 0; iVar < nVar_Total; iVar++) 
	//	printf("%d ",VarTyp[iVar]);
	//printf("\n");
	
	if ( !GmfSetKwd(OutSol, GmfSolAtVertices, npoin, nVar_Total, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}

  for (iPoint = 0; iPoint < npoin; iPoint++) {
	
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
			bufDbl[iVar] = SU2_TYPE::GetValue(Data[iVar][iPoint]);
    }

		GmfSetLin(OutSol, GmfSolAtVertices, bufDbl);
	}
  
	if ( !GmfCloseMesh(OutSol) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file %s.", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}  
	
}




/*

	Write solutions of interest : mach, dens, pres etc.

*/
void COutput::WriteInriaOutputs(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  //ofstream restart_file;
  string filename;
  
	unsigned long OutMach, OutPres, i, npoin = geometry->GetGlobal_nPointDomain();
	int VarTyp[GmfMaxTyp];
	double bufDbl[GmfMaxTyp];
	char OutNam[1024], BasNam[1024];
	char *ptr=NULL;
	
	int NbrVar, idxVar;
	
	/* Get indices of mach, pres, etc. in the solution array */
	unsigned short *TagBc;
	TagBc = new unsigned short [100];
	
	idxVar=0;
  idxVar += nVar_Consv;
  
  if (!config->GetLow_MemoryOutput()) {
    
    if (config->GetWrt_Limiters()) {
      idxVar += nVar_Consv;
    }
    if (config->GetWrt_Residuals()) {
      idxVar += nVar_Consv;
    }
		
    if (config->GetKind_Regime() == FREESURFACE) {
      TagBc[bcDens] = idxVar;
      idxVar++;
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      TagBc[bcPres] = idxVar;
      TagBc[bcTemp] = idxVar+1;
      TagBc[bcMach] = idxVar+3;
      idxVar += 4;
    }
		
	}
	
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
	
	/*--- Get output name *.solb ---*/
	

	
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);
	
  /*--- Open the restart file and write the solution. ---*/
	sprintf(OutNam, "mach.solb");
	OutMach = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutMach ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
	sprintf(OutNam, "pres.solb");
	OutPres = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutPres ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
	
  /*--- Write MACH ---*/
	
	npoin = geometry->GetGlobal_nPointDomain();
		
	NbrVar = 1;
	VarTyp[0]  = GmfSca;
	
	if ( !GmfSetKwd(OutPres, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write Mach");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  for (iPoint = 0; iPoint < npoin; iPoint++) {
		
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
			iVar = TagBc[bcMach];
			bufDbl[0] = SU2_TYPE::GetValue(Data[iVar][iPoint]);
			GmfSetLin(OutMach, GmfSolAtVertices, bufDbl);
    }
		
	}
	
	
  /*--- Write PRES ---*/
	
	npoin = geometry->GetGlobal_nPointDomain();
		
	NbrVar = 1;
	VarTyp[0]  = GmfSca;
	
	if ( !GmfSetKwd(OutPres, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write pressure");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  for (iPoint = 0; iPoint < npoin; iPoint++) {
		iVar = TagBc[bcPres];
		bufDbl[0] = SU2_TYPE::GetValue(Data[iVar][iPoint]);
		GmfSetLin(OutPres, GmfSolAtVertices, bufDbl);
	}
		
	/*--- Close files ---*/
  
	if ( !GmfCloseMesh(OutMach) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}  
	
	if ( !GmfCloseMesh(OutPres) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
	delete [] TagBc;
	
}




void COutput::SetInriaMesh(CConfig *config, CGeometry *geometry) {
  
  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
  unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], iNode, nElem;
  unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, iPeriodic, nPeriodic = 0, VTK_Type, nMarker_;
  su2double *center, *angles, *transl;
  ofstream output_file;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag, str;
  string::size_type position;
	
	unsigned short nMarker = config->GetnMarker_All();
	unsigned long cptElem = 0, nTri=0, nLin=0, nQua=0;

	int Dim;
	int OutMsh,i;
  int iVer,iTri,iEfr,iTet;
	double bufDbl[8];
	char OutNam[2014];
	int bufInt[8];
	
	unsigned long *PointSurface=NULL;
	unsigned long nPointSurface=0;
	
	CPrimalGrid* bnd = NULL;

  /*--- Read the name of the output and input file ---*/
  str = config->GetMesh_Out_FileName();

  //strcpy (out_file, str.c_str());
  //strcpy (cstr, out_file);
  //output_file.precision(15);
  //output_file.open(cstr, ios::out);

	sprintf(OutNam, "%s.meshb", str.c_str());
	
	Dim = nDim;
	if ( !(OutMsh = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,Dim)) ) {
    printf("  ## ERROR: Cannot open mesh file %s ! \n",OutNam);
		return;
  }
  
	/*--- Write vertices ---*/
	
	GmfSetKwd(OutMsh, GmfVertices, nGlobal_Poin);
	
  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++)
			bufDbl[iDim] = SU2_TYPE::GetValue(Coords[iDim][iPoint]);
		
		if ( nDim == 2 ) {
			GmfSetLin(OutMsh, GmfVertices,bufDbl[0],bufDbl[1],0); 
		}
		else {
			GmfSetLin(OutMsh, GmfVertices,bufDbl[0],bufDbl[1],bufDbl[2],0); 
		}
  }

	/*--- Write 2D elements ---
	Note: in 3D, triangles/quads are boundary markers
	*/
	
	if ( nDim == 2 ){
		
		/*--- Write triangles ---*/
		
		GmfSetKwd(OutMsh, GmfTriangles, nGlobal_Tria);
  	for (iElem = 0; iElem < nGlobal_Tria; iElem++) {
  	  iNode = iElem*N_POINTS_TRIANGLE;
			GmfSetLin(OutMsh, GmfTriangles,Conn_Tria[iNode+0],Conn_Tria[iNode+1],Conn_Tria[iNode+2], 0);  
  	}	

		/*--- Write quadrilaterals ---*/
		
		if ( nGlobal_Quad > 0  ) {
			GmfSetKwd(OutMsh, GmfQuadrilaterals, nGlobal_Quad);
			for (iElem = 0; iElem < nGlobal_Quad; iElem++) {
  		  iNode = iElem*N_POINTS_QUADRILATERAL;
				GmfSetLin(OutMsh, GmfQuadrilaterals,Conn_Quad[iNode+0],Conn_Quad[iNode+1],Conn_Quad[iNode+2], Conn_Quad[iNode+3], 0);  
  		}
		}
	
	}
	
	/*--- Write tetrahedra ---*/
	
	
	
	if ( nGlobal_Tetr > 0  ) {
		GmfSetKwd(OutMsh, GmfTetrahedra, nGlobal_Tetr);
		for (iElem = 0; iElem < nGlobal_Tetr; iElem++) {
	    iNode = iElem*N_POINTS_TETRAHEDRON;
			GmfSetLin(OutMsh, GmfTetrahedra,Conn_Tetr[iNode+0],Conn_Tetr[iNode+1],Conn_Tetr[iNode+2], Conn_Tetr[iNode+3], 0); 
	  }
	}
	
	/*--- Write hexahedra ---*/
	
	if ( nGlobal_Hexa > 0 ) {
		GmfSetKwd(OutMsh, GmfHexahedra, nGlobal_Hexa);
		for (iElem = 0; iElem < nGlobal_Hexa; iElem++) {
	    iNode = iElem*N_POINTS_HEXAHEDRON;
			GmfSetLin(OutMsh, GmfHexahedra,Conn_Hexa[iNode+0],Conn_Hexa[iNode+1], Conn_Hexa[iNode+2], Conn_Hexa[iNode+3], Conn_Hexa[iNode+4],Conn_Hexa[iNode+5],Conn_Hexa[iNode+6], Conn_Hexa[iNode+7],  0); 
	  }
	}
	
	/*--- Write prisms ---*/
	
	if ( nGlobal_Pris > 0 ) {
		GmfSetKwd(OutMsh, GmfPrisms, nGlobal_Pris);
		for (iElem = 0; iElem < nGlobal_Pris; iElem++) {
	    iNode = iElem*N_POINTS_PRISM;
			GmfSetLin(OutMsh, GmfPrisms,Conn_Pris[iNode+0],Conn_Pris[iNode+1], Conn_Pris[iNode+2], Conn_Pris[iNode+3], Conn_Pris[iNode+4],Conn_Pris[iNode+5],  0); 
	  }
	}
	
	/*--- Write pyramids ---*/
	
	if ( nGlobal_Pyra > 0 ) {
		GmfSetKwd(OutMsh, GmfPyramids, nGlobal_Pyra);
		for (iElem = 0; iElem < nGlobal_Pyra; iElem++) {
	    iNode = iElem*N_POINTS_PYRAMID;
	  	GmfSetLin(OutMsh, GmfPyramids,Conn_Pyra[iNode+0],Conn_Pyra[iNode+1], Conn_Pyra[iNode+2], Conn_Pyra[iNode+3], Conn_Pyra[iNode+4],0); 
		}
	}
	
	
	/* --- Boundary elements ---*/
	
	/*--- Get surface points ---*/
	
	nPointSurface = 0;
	PointSurface = new unsigned long[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (geometry->node[iPoint]->GetBoundary()) {
      PointSurface[nPointSurface] = iPoint;
      nPointSurface++;
    }
	}
	
	/*--- Count elements ---*/
	
	nLin = nTri = nQua = 0;
	
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
			bnd = geometry->bound[iMarker][iElem];
			switch ( bnd->GetVTK_Type() ) {
				case LINE:          nLin++; break;
				case TRIANGLE:      nTri++; break;
				case QUADRILATERAL: nQua++; break;
			}
		}
	}
	
	
	/*--- Write edges ---*/
	
	if ( nLin > 0 ) {
		
		GmfSetKwd(OutMsh, GmfEdges, nLin);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != LINE ) 
					continue;
				
				cptElem++;

				
				GmfSetLin(OutMsh, GmfEdges,bnd->GetNode(0)+1,bnd->GetNode(1)+1,iMarker); 	
			}
		}
	
		if ( cptElem != nLin ) {
			cout << "  !! Error Inria output:  Inconsistent number of edges\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	/*--- Write triangles ---*/
	
	if ( nTri > 0 ) {
		
		GmfSetKwd(OutMsh, GmfTriangles, nTri);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != TRIANGLE ) 
					continue;
				
				
				if ( cptElem < 10 ) 
					cout << bnd->GetNode(0)+1 << " " << bnd->GetNode(1)+1 << " " << bnd->GetNode(2)+1 << endl;
				
				cptElem++;
				
				GmfSetLin(OutMsh, GmfTriangles,PointSurface[bnd->GetNode(0)]+1,PointSurface[bnd->GetNode(1)]+1, PointSurface[bnd->GetNode(2)]+1,iMarker); 	
			}
		}
	
		if ( cptElem != nTri ) {
			cout << "  !! Error Inria output:  Inconsistent number of triangles\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	/*--- Write quadrilaterals ---*/
	
	if ( nQua > 0 ) {
		
		GmfSetKwd(OutMsh, GmfQuadrilaterals, nQua);
		
		cptElem = 0;
		
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
				
				bnd = geometry->bound[iMarker][iElem];
				
				if ( bnd->GetVTK_Type() != QUADRILATERAL ) 
					continue;
				
				cptElem++;
				
				GmfSetLin(OutMsh, GmfTriangles,bnd->GetNode(0)+1,bnd->GetNode(1)+1, bnd->GetNode(2)+1, bnd->GetNode(3)+1,iMarker); 	
			}
		}
	
		if ( cptElem != nQua ) {
			cout << "  !! Error Inria output:  Inconsistent number of quadrilaterals\n" << endl;
			exit(EXIT_FAILURE);
		}
		
	}
	
	if ( PointSurface )
		delete [] PointSurface;
	
	GmfCloseMesh(OutMsh);
  
}

