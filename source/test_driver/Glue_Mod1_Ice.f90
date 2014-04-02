!..................................................................................................................................
! LICENSING                                                                                                                         
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    This file is the "glue code" example for the new FAST modularization.  
!
!    This code embodies the formulation and numerical methods for "Predictor-Corrector Loose Coupling."   The methods
!    employed here share similarities with those described in Gasmi et al. (2013).   However, the method used here is a 
!    "symmetric" predictor-corrector approach (order of module UpdateStates does not matter).    Further, the method reduces
!    to an explicit method when only the prediction step is taken (pc_max = 1 below).   However, for all modules, inputs
!    and outputs are stored for up to three steps, allowing up to quadratic interpolation or exptrapolation of input and 
!    output data.
!
!    The test problem is a simple two-degree-of-freedom damped oscillator, where each "mass" is treated by a module; see
!    Gasmi et al. (2013) for details.
!
!    Three fourth-order explicit numerical time integrators are included: Runge-Kutta (RK4), Adams-Bashforth (AB4), and 
!    Adams-Bashforth-Moulton (ABM4).    RK4 and ABM4 have an implcit dependence on other-module data.
!
!    Numerical experiments have shown that, if quadratic interpolation of inputs and outpus is employed, order of accuracy of 
!    the methods with pc_max predictor-corrector iterations are as follows (with Mod1 & Ice using same integrator):
!
!    RK4, PC1: third order
!    RK4, PC2: third order (but more accurate than PC1)
!
!    AB4, PC1: fourth order
!    AB4, PC2: fourth order (should be identical to PC1; no implicit dependence on other-module data)
!     
!    ABM4, PC1: third order
!    ABM4, PC2: fourth order
!
!    NOTE: These convergence results can be obtained only when the multi-step methods have their first three steps initialized
!          with the exact benchmark solution.
!	
!    References:
!
!    Gasmi, A., M. A. Sprague, J. M. Jonkman, and W. B. Jones, Numerical stability and accuracy of temporally coupled 
!    multi-physics modules in wind turbine CAE tools. In proceedings of the 32nd ASME Wind Energy Symposium, 51st AIAA 
!    Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, Grapevine, TX, January 7-10, 
!    2013.   Also published as NREL Report No. CP-2C00-57298.   Available in pdf format at:
!    http://www.nrel.gov/docs/fy13osti/57298.pdf
!
!**********************************************************************************************************************************
PROGRAM MAIN

   USE Module1
   USE Module1_Types

   USE IceDyn
   USE IceDyn_Types

   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   CHARACTER(1024)                    :: OutFileName      ! Name of the output file

   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                         :: t_initial        ! time at initialization
   REAL(DbKi)                         :: t_final          ! time at simulation end 
   REAL(DbKi)                         :: t_global         ! global-loop time marker

   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter

   INTEGER(IntKi)                     :: pc_max           ! 1:explicit loose; 2:pc loose
   INTEGER(IntKi)                     :: pc               ! counter for pc iterations

   INTEGER(IntKi)                     :: Mod1_interp_order     ! order of interpolation/extrapolation
   INTEGER(IntKi)                     :: ID_interp_order      ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InitInputType)           :: Mod1_InitInput
   TYPE(Mod1_ParameterType)           :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousStateDeriv
   TYPE(Mod1_InitOutputType)          :: Mod1_InitOutput
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType)          :: Mod1_OtherState

   TYPE(Mod1_InputType),Dimension(:),Allocatable   :: Mod1_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_InputTimes

   TYPE(Mod1_OutputType),Dimension(:),Allocatable  :: Mod1_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_OutputTimes

   TYPE(Mod1_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(Mod1_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module IceLoad deived data typed needed in pc-coupling; predicted states

   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState_pred
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState_pred
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState_pred

   ! Module IceLoad Derived-types variables; see Registry_Module2.txt

   TYPE(ID_InitInputType)           :: ID_InitInput
   TYPE(ID_ParameterType)           :: ID_Parameter
   TYPE(ID_ContinuousStateType)     :: ID_ContinuousState
   TYPE(ID_ContinuousStateType)     :: ID_ContinuousStateDeriv
   TYPE(ID_InitOutputType)          :: ID_InitOutput
   TYPE(ID_DiscreteStateType)       :: ID_DiscreteState
   TYPE(ID_ConstraintStateType)     :: ID_ConstraintState
   TYPE(ID_OtherStateType)          :: ID_OtherState

   ! Module IceLoad deived data typed needed in pc-coupling; predicted states

   TYPE(ID_ContinuousStateType)     :: ID_ContinuousState_pred
   TYPE(ID_DiscreteStateType)       :: ID_DiscreteState_pred
   TYPE(ID_ConstraintStateType)     :: ID_ConstraintState_pred

   TYPE(ID_InputType),Dimension(:),Allocatable    :: ID_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ID_InputTimes

   TYPE(ID_OutputType),Dimension(:),Allocatable  :: ID_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ID_OutputTimes

   TYPE(ID_InputType)   :: u2    ! local variable for extrapolated inputs
   TYPE(ID_OutputType)  :: y2    ! local variable for extrapolated outputs

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization
   CHARACTER(200)                     :: Frmt
   CHARACTER(1)                       :: Dlim = TAB
   Integer(IntKi)                     :: Un

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.
   t_final   = 250.

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.0125

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   Mod1_interp_order = 2 
   ID_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(Mod1_Input(Mod1_interp_order + 1)) 
   Allocate(Mod1_InputTimes(Mod1_interp_order + 1)) 

   Allocate(Mod1_Output(Mod1_interp_order + 1)) 
   Allocate(Mod1_OutputTimes(Mod1_interp_order + 1)) 

   ! Module2: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(ID_Input(ID_interp_order + 1)) 
   Allocate(ID_InputTimes(ID_interp_order + 1)) 

   Allocate(ID_Output(ID_interp_order + 1)) 
   Allocate(ID_OutputTimes(ID_interp_order + 1)) 

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL Mod1_Init( Mod1_InitInput, Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), dt_global, Mod1_InitOutput, ErrStat, ErrMsg )


   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod1_Input)
   do i = 1, Mod1_interp_order + 1  
      Mod1_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod1_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod1_interp_order
     Call Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  0, Errstat, ErrMsg)
     Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  0, Errstat, ErrMsg)
   enddo

   ID_InitInput%InputFile = 'IceDyn_Input.dat'
   ID_InitInput%RootName = 'IceDyn_Mod1Test'
   
   CALL ID_Init( ID_InitInput, ID_Input(1), ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                   ID_ConstraintState, ID_OtherState, ID_Output(1), dt_global, ID_InitOutput, ErrStat, ErrMsg )

   do i = 1, ID_interp_order + 1  
      ID_InputTimes(i) = t_initial - (i - 1) * dt_global
      ID_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, ID_interp_order
     Call ID_CopyInput (ID_Input(i),  ID_Input(i+1),  0, Errstat, ErrMsg)
     Call ID_CopyOutput (ID_Output(i),  ID_Output(i+1),  0, Errstat, ErrMsg)
   enddo


   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! Open output file
   OutFileName = Trim(ID_parameter%RootName)//'.txt'
   CALL GetNewUnit (Un)
   CALL OpenFOutFile (Un, OutFileName, ErrStat, ErrMsg)
   
   ! write headers for output columns:
   Frmt = '(A8)'
   WRITE (Un, Frmt, ADVANCE = 'no') TRIM ('Time')
   
   Frmt = '(3(:,A,A10))'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, 'StrDisp', Dlim, 'IceDisp', Dlim, 'IceForce'
   
   ! CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   ! CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   ! CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(Mod1_ContinuousState%q))   
   
   Frmt = '(/ F8.3)'
   WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
   Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, Mod1_ContinuousState%q, Dlim, ID_ContinuousState%q, Dlim, ID_Output(1)%fice
   
   DO n_t_global = 0, n_t_final

      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules

      call Mod1_ID_InputOutputSolve(t_global, &
                   Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), &
                   ID_Input(1), ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                   ID_ConstraintState, ID_OtherState, ID_Output(1),  &
                   ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod1_Input_ExtrapInterp(Mod1_Input, Mod1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod1_Output_ExtrapInterp(Mod1_Output, Mod1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod1_interp_order, 1, -1
         Call Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  0, Errstat, ErrMsg)
         Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  0, Errstat, ErrMsg)
         Mod1_InputTimes(i+1) = Mod1_InputTimes(i)
         Mod1_OutputTimes(i+1) = Mod1_OutputTimes(i)
      enddo

      Call Mod1_CopyInput (u1,  Mod1_Input(1),  0, Errstat, ErrMsg)
      Call Mod1_CopyOutput (y1,  Mod1_Output(1),  0, Errstat, ErrMsg)
      Mod1_InputTimes(1) = t_global + dt_global
      Mod1_OutputTimes(1) = t_global + dt_global

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ID_Input_ExtrapInterp(ID_Input, ID_InputTimes, u2, t_global + dt_global, ErrStat, ErrMsg)

      CALL ID_Output_ExtrapInterp(ID_Output, ID_OutputTimes, y2, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Ice_Input and Ice_Output

      do i = ID_interp_order, 1, -1
         Call ID_CopyInput (ID_Input(i),  ID_Input(i+1),  0, Errstat, ErrMsg)
         Call ID_CopyOutput (ID_Output(i),  ID_Output(i+1),  0, Errstat, ErrMsg)
         ID_InputTimes(i+1) = ID_InputTimes(i)
         ID_OutputTimes(i+1) = ID_OutputTimes(i)
      enddo

      Call ID_CopyInput (u2,  ID_Input(1),  0, Errstat, ErrMsg)
      Call ID_CopyOutput (y2,  ID_Output(1),  0, Errstat, ErrMsg)
      ID_InputTimes(1) = t_global + dt_global
      ID_OutputTimes(1) = t_global + dt_global

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Module 1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod1_CopyContState   (Mod1_ContinuousState, Mod1_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod1_CopyConstrState (Mod1_ConstraintState, Mod1_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod1_CopyDiscState   (Mod1_DiscreteState,   Mod1_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod1_UpdateStates( t_global, n_t_global, Mod1_Input, Mod1_InputTimes, Mod1_Parameter, Mod1_ContinuousState_pred, &
                                 Mod1_DiscreteState_pred, Mod1_ConstraintState_pred, &
                                 Mod1_OtherState, ErrStat, ErrMsg )

         !----------------------------------------------------------------------------------------
         ! Module 2
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call ID_CopyContState   (ID_ContinuousState, ID_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call ID_CopyConstrState (ID_ConstraintState, ID_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call ID_CopyDiscState   (ID_DiscreteState,   ID_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL ID_UpdateStates( t_global, n_t_global, ID_Input, ID_InputTimes, ID_Parameter, ID_ContinuousState_pred, &
                                 ID_DiscreteState_pred, ID_ConstraintState_pred, &
                                 ID_OtherState, ErrStat, ErrMsg )

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .lt. pc_max) then

            call Mod1_ID_InputOutputSolve( t_global + dt_global, &
                                             Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
                                             Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
                                             ID_Input(1), ID_Parameter, ID_ContinuousState_pred, ID_DiscreteState_pred, &
                                             ID_ConstraintState_pred, ID_OtherState, ID_Output(1),  &
                                             ErrStat, ErrMsg)

         endif

      enddo

      ! Save all final variables 

      Call Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)

      Call ID_CopyContState   (ID_ContinuousState_pred,  ID_ContinuousState, 0, Errstat, ErrMsg)
      Call ID_CopyConstrState (ID_ConstraintState_pred,  ID_ConstraintState, 0, Errstat, ErrMsg)
      Call ID_CopyDiscState   (ID_DiscreteState_pred,    ID_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial
      
      Frmt = '(/ F8.3)'
      WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
      Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
      WRITE (Un, Frmt, ADVANCE = 'no') Dlim, Mod1_ContinuousState%q, Dlim, ID_ContinuousState%q, Dlim, ID_Output(1)%fice


   END DO


   

   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL Mod1_End(  Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), ErrStat, ErrMsg )

   do i = 2, Mod1_interp_order+1
      CALL Mod1_DestroyInput(Mod1_Input(i), ErrStat, ErrMsg )
      CALL Mod1_DestroyOutput(Mod1_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod1_InputTimes)
   DEALLOCATE(Mod1_OutputTimes)

   CALL ID_End(  ID_Input(1), ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                   ID_ConstraintState, ID_OtherState, ID_Output(1), ErrStat, ErrMsg )
   
   do i = 2, ID_interp_order+1
      CALL ID_DestroyInput(ID_Input(i), ErrStat, ErrMsg )
      CALL ID_DestroyOutput(ID_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(ID_InputTimes)
   DEALLOCATE(ID_Output)

END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_ID_InputOutputSolve(time, &
                   Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, &
                   ID_Input, ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                   ID_ConstraintState, ID_OtherState, ID_Output,  &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in 
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................
 
   USE Module1
   USE Module1_Types

   USE IceDyn
   USE IceDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details
 
   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_ParameterType),       INTENT(IN   ) :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType), INTENT(IN   ) :: Mod1_ContinuousState
   TYPE(Mod1_DiscreteStateType),   INTENT(IN   ) :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType), INTENT(INOUT) :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType),      INTENT(INOUT) :: Mod1_OtherState
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! IceLoad Derived-types variables; see Registry_Module2.txt

   TYPE(ID_InputType),           INTENT(INOUT) :: ID_Input
   TYPE(ID_ParameterType),       INTENT(IN   ) :: ID_Parameter
   TYPE(ID_ContinuousStateType), INTENT(IN   ) :: ID_ContinuousState
   TYPE(ID_DiscreteStateType),   INTENT(IN   ) :: ID_DiscreteState
   TYPE(ID_ConstraintStateType), INTENT(INOUT) :: ID_ConstraintState
   TYPE(ID_OtherStateType),      INTENT(INOUT) :: ID_OtherState
   TYPE(ID_OutputType),          INTENT(INOUT) :: ID_Output

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds

   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules; could be placed in a separate routine.
   ! Note that Module2 has direct feedthrough, but Module1 does not. Thus, Module1 should be called first.

   CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

   call ID_InputSolve( Mod1_Input, Mod1_Output, ID_Input, ID_Output, ErrStat, ErrMsg)

   CALL ID_CalcOutput( time, ID_Input, ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                ID_ConstraintState, ID_OtherState, ID_Output, ErrStat, ErrMsg )

   call Mod1_InputSolve( Mod1_Input, Mod1_Output, ID_Input, ID_Output, ErrStat, ErrMsg)
 
END SUBROUTINE Mod1_ID_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_InputSolve( Mod1_Input, Mod1_Output,  ID_Input, ID_Output, ErrStat, ErrMsg)
 
   USE Module1
   USE Module1_Types

   USE IceDyn
   USE IceDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(  OUT) :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(IN   ) :: Mod1_Output

   ! Module2 Derived-types variables; see Registry_Module2.txt
 
   TYPE(ID_InputType),           INTENT(IN   ) :: ID_Input
   TYPE(ID_OutputType),          INTENT(IN   ) :: ID_Output

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg  = ''

   Mod1_Input%fc   = ID_Output%fice

 
END SUBROUTINE Mod1_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_InputSolve( Mod1_Input,  Mod1_Output, ID_Input,  ID_Output, ErrStat, ErrMsg)
 
   USE Module1
   USE Module1_Types

   USE IceDyn
   USE IceDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(IN   ) :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(IN   ) :: Mod1_Output

   ! Module2 Derived-types variables; see Registry_Module2.txt
 
   TYPE(ID_InputType),           INTENT(  OUT) :: ID_Input
   TYPE(ID_OutputType),          INTENT(IN   ) :: ID_Output

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   ErrStat = ErrID_None
   ErrMsg  = ''   
   
   ID_Input%q    = Mod1_Output%q
   ID_Input%dqdt = Mod1_Output%dqdt
 
 
END SUBROUTINE ID_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
