PROGRAM DiffusionEquationWithLinearSource

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Program parmaeters
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP/2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP/2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: coordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: regionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: basisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: generatedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: meshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: decompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: geometricFieldUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: equationsSetFieldUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: dependentFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: materialsFieldUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: sourceFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: analyticFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: equationsSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: problemUserNumber=1
  
  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  
  !CMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem,worldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,sourceField,analyticField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh  
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver, linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: equationsSetIndex
  INTEGER(CMISSIntg) :: firstNodeNumber,lastNodeNumber
  INTEGER(CMISSIntg) :: err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  !Intialise OpenCMISS
  CALL cmfe_Initialise(worldCoordinateSystem,worldRegion,err)
  
  !Set the random seeds so we can test multi process
  CALL cmfe_RandomSeedsSet(9999,err)
  
  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,numberOfComputationalNodes,err)
  CALL cmfe_ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,computationalNodeNumber,err)

  !Set output on
  CALL cmfe_OutputSetOn("DiffusionWithLinearSource",err)

  numberOfGlobalXElements=2
  numberOfGlobalYElements=4
  numberOfGlobalZElements=4

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(coordinateSystemUserNumber,coordinateSystem,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(regionUserNumber,worldRegion,region,err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_LabelSet(region,"diffusion_equation_linear_source",err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------   
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(basisUserNumber,basis,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a biquadratic Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(basis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],err)    
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[3,3],err)
  ELSE
    !Set the basis to be a triquadratic Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(basis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(basis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
      & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],err)    
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[3,3,3],err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(generatedMeshUserNumber,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  IF(numberOfGlobalZElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MeshUserNumber,Mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(decompositionUserNumber,Mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(decomposition,numberOfComputationalNodes,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(geometricFieldUserNumber,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)
  
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  !Set the equations set to be a linear source diffusion problem
  CALL cmfe_EquationsSet_CreateStart(equationsSetUserNumber,region,geometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE],equationsSetFieldUserNumber, &
    & equationsSetField,equationsSet,err)  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,dependentFieldUserNumber,dependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MaterialsFieldUserNumber,materialsField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 4,-1.0_CMISSRP,err)
  
  !Create the equations set source field variables
  CALL cmfe_Field_Initialise(sourceField,err)
  CALL cmfe_EquationsSet_SourceCreateStart(equationsSet,SourceFieldUserNumber,sourceField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(equationsSet,err)
  CALL cmfe_Field_ComponentValuesInitialise(sourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,0.0_CMISSRP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTICAL FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(analyticField,err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSet,CMFE_EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_THREE_DIM_1, & 
    & analyticFieldUserNumber,analyticField,err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSet,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !----------------------------------------------------------------------------------------------------------- 
  
  !Create the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(problemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE, &
    & CMFE_PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOP
  !----------------------------------------------------------------------------------------------------------- 
  
  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(controlLoop,0.0_CMISSRP,1.0000_CMISSRP,0.005_CMISSRP,err)
  !Set the output
  CALL cmfe_ControlLoop_OutputTypeSet(controlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVERs
  !-----------------------------------------------------------------------------------------------------------
  
  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  CALL cmfe_Solver_DynamicLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolver,1000,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem solver equations
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------
  
  !Solve the problem
  CALL cmfe_Problem_Solve(problem,err)

  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(dependentField,"diffusion_equation_linear_source",err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------
  
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"diffusion_equation_linear_source","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"diffusion_equation_linear_source","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  CALL cmfe_Finalise(err)
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM DiffusionEquationWithLinearSource
