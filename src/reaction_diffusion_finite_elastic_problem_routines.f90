!> \file
!> \author Vijay Rajagopal
!> \brief This module handles reaction diffusion mechanics problem routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all reaction diffusion equation routines.
MODULE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_ROUTINES
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  
  PUBLIC REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE
  
  PUBLIC REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP

  PUBLIC REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET

  PUBLIC REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE
  
  PUBLIC REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP


  CONTAINS


  !
  !================================================================================================================================
  !
  !>Sets/changes the equation subtype for a reaction-diffusion-mechanics problem type of a multiphysics problem class.
  SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER, INTENT(IN) :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
        PROBLEM%CLASS=PROBLEM_MULTI_PHYSICS_CLASS
        PROBLEM%TYPE=PROBLEM_REACTION_DIFFUSION_FINITE_ELASTIC_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="The specified problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid or implemented for a reaction-diffusion-mechanics problem type of the multi-physics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT

    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !
  !>Sets up the reaction-diffusion-mechanics problem.
  SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a bioelectric domain equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT, &
      & MECH_CONTROL_LOOP,REACDIFF_CONTROL_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: REACDIFF_SOLVERS,MECH_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(CELLML_EQUATIONS)
    NULLIFY(CONTROL_LOOP)
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(MECH_CONTROL_LOOP)
    NULLIFY(REACDIFF_CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(REACDIFF_SOLVERS)
    NULLIFY(MECH_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    
    CALL ENTERS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing????
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a reaction diffusion and mechanics problem type."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)

          !Set up a time control loop and a simple and load increment subloops
          !create start sets up a root control loop - passed back as CONTROL_LOOP
          CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
          !setting up the root control as a time loop.
          CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          !create two sub loops
          CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
          CALL CONTROL_LOOP_LABEL_SET(CONTROL_LOOP,"Time Loop",ERR,ERROR,*999)

          !get the first subloop of the root loop.
          CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_TYPE_SET(REACDIFF_CONTROL_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
          CALL CONTROL_LOOP_LABEL_SET(ReacDiff_CONTROL_LOOP,"Simple Reaction Diffusion",ERR,ERROR,*999)

          !get the second subloop of the root loop and set up mechanics load increment type
          CALL CONTROL_LOOP_GET(CONTROL_LOOP,[2,CONTROL_LOOP_NODE],MECH_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_TYPE_SET(MECH_CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
          CALL CONTROL_LOOP_LABEL_SET(MECH_CONTROL_LOOP,"Mechanics Load Loop",ERR,ERROR,*999)


        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          !CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          !CALL CONTROL_LOOP_CREATE_FINISH(REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          !CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[2,CONTROL_LOOP_NODE],MECH_CONTROL_LOOP,ERR,ERROR,*999)
          !CALL CONTROL_LOOP_CREATE_FINISH(MECH_CONTROL_LOOP,ERR,ERROR,*999)

          CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP_ROOT,ERR,ERROR,*999)
            
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a reaction-diffusion-mechanics problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        !Get the control loop
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        !get the two sub loops
        CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
        CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[2,CONTROL_LOOP_NODE],MECH_CONTROL_LOOP,ERR,ERROR,*999)
        SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Start the solvers creation
          CALL SOLVERS_CREATE_START(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_CREATE_START(MECH_CONTROL_LOOP,MECH_SOLVERS,ERR,ERROR,*999)
          SELECT CASE(PROBLEM%SUBTYPE)
          CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
            CALL SOLVERS_NUMBER_SET(REACDIFF_SOLVERS,3,ERR,ERROR,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"First ODE solver",ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the third solver to be a differential-algebraic equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Second ODE solver",ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !Set the solver for the mech_solvers to be a nonlinear equations solver
            NULLIFY(SOLVER)
            CALL SOLVERS_NUMBER_SET(MECH_SOLVERS,1,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(MECH_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear Solver",ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a reaction-diffusion-mechanics problem type of the multi-physics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          CALL CONTROL_LOOP_SOLVERS_GET(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(MECH_CONTROL_LOOP,MECH_SOLVERS,ERR,ERROR,*999)
          !Finish the solvers creation
          CALL SOLVERS_CREATE_FINISH(REACDIFF_SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_CREATE_FINISH(MECH_SOLVERS,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for this multi-physics reaction-diffusion-mechanics problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[2,CONTROL_LOOP_NODE],MECH_CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM%SUBTYPE)
          CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MECH_CONTROL_LOOP,MECH_SOLVERS,ERR,ERROR,*999)
            !Create the solver equations for the second (parabolic) in the reaction diffusion control loop
            !and the first (nonlinear) solver in the mechanics control loop
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(MECH_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)

          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a reaction-diffusion-mechanics problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loops and solvers
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[2,CONTROL_LOOP_NODE],MECH_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(MECH_CONTROL_LOOP,MECH_SOLVERS,ERR,ERROR,*999)

          SELECT CASE(PROBLEM%SUBTYPE)
          CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
            !Get the solver equations for the second (parabolic) solver in first subloop
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)       
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)  
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)       
            !Get the solver equations for the first (nonlinear) solver in second subloop
            CALL SOLVERS_SOLVER_GET(MECH_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
    
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a reaction-diffusion-mechanics problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a reaction-diffusion-mechanics problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          !Get the solver
          CALL CONTROL_LOOP_SOLVERS_GET(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)

          IF(PROBLEM%SUBTYPE==PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE) THEN
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Create the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(REACDIFF_CONTROL_LOOP,REACDIFF_SOLVERS,ERR,ERROR,*999)
          SELECT CASE(PROBLEM%SUBTYPE)
          CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
            !Get the CellML equations for the second DAE solver
            NULLIFY(SOLVER)
            NULLIFY(CELLML_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,3,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for reaction-diffusion-mechanics problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for reaction-diffusion-mechanics problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
          & " is invalid for a reaction-diffusion-mechanics problem setup."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
  !>Performs pre-solve actions for reaction-diffusion-mechanics problems.
  SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to perform the pre-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          !CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            SELECT CASE(PROBLEM%SUBTYPE)
            CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
              SELECT CASE(SOLVER%GLOBAL_NUMBER)
              CASE(1)
                !CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT/2.0_DP,ERR,ERROR,*999)
              CASE(2)
                !Do nothing
              CASE(3)
                !CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME+TIME_INCREMENT/2.0_DP,CURRENT_TIME+TIME_INCREMENT, &
                 ! & ERR,ERROR,*999)
              CASE(4)
                !Do nothing
              CASE DEFAULT
                LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for a Strang split reaction-diffusion problem."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is invalid for a reaction-diffusion problem type."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Control loop problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE")
    RETURN
999 CALL ERRORS("REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE")
    RETURN 1
    
  END SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(SOLVERS_TYPE), POINTER :: REACDIFF_SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(CONTROL_LOOP_TYPE), POINTER :: REACDIFF_CONTROL_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    NULLIFY(REACDIFF_CONTROL_LOOP)
    NULLIFY(REACDIFF_SOLVERS)
    NULLIFY(SOLVER)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        ! Eventually we may want to do different things depending on problem type/subtype
        ! too but for now we can just check the loop type.
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
            !get the first sub loop containing reaction diffusion solvers
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,[1,CONTROL_LOOP_NODE],REACDIFF_CONTROL_LOOP,ERR,ERROR,*999)

            SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
              CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
                REACDIFF_SOLVERS=>REACDIFF_CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(REACDIFF_SOLVERS)) THEN
                  CALL SOLVERS_SOLVER_GET(REACDIFF_SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  SELECT CASE(SOLVER%GLOBAL_NUMBER)
                  CASE(1)
                    CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT/2.0_DP,ERR,ERROR,*999)
                  CASE(2)
                    !do nothing
                  CASE(3)
                    CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME+TIME_INCREMENT/2.0_DP,CURRENT_TIME+TIME_INCREMENT, &
                      & ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for a Strang split reaction-diffusion problem."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
                ENDIF

              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                  & " is not valid for REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT

          CASE DEFAULT
            !do nothing
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE REACDIFF_FINITE_ELASTIC_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Sets up the reaction diffusion problem post solve.
  SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STRANGSPLIT_REACDIFF_FINITE_ELASTIC_SUBTYPE)
              !OUTPUT SOLUTIONS AT EACH TIME STEP
              !CALL REACTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a reaction diffusion type of a classical field problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE")
    RETURN
999 CALL ERRORS("REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE",ERR,ERROR)
    CALL EXITS("REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE")
    RETURN 1
  END SUBROUTINE REACTION_DIFFUSION_FINITE_ELASTIC_POST_SOLVE


END MODULE REACTION_DIFFUSION_FINITE_ELASTIC_PROBLEM_ROUTINES
