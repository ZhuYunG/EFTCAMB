!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2019 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 10p5_Horndeski.f90
!! This file contains the definition of the Horndeski full mapping model.

!> @author Cheng-Zhi Dou

module EFTCAMB_FM_horndeski

    use precision
    use constants, only : c, const_pi, kappa, Mpc, G
    use IniFile
    use AMLutils
    use equispaced_linear_interpolation_1D
    use EFTCAMB_cache
    use EFT_def
    use EFTCAMB_mixed_algorithms, only : double_NaN
    use EFTCAMB_rootfind
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_parametrizations_1D
    use EFTCAMB_power_law_parametrizations_1D
    use EFTCAMB_constant_parametrization_1D
    use EFTCAMB_taylorseries_parametrizations_1D
    use EFTCAMB_padeseries_parametrizations_1D
    use EFTCAMB_fourier_parametrizations_1D
    use EFTCAMB_exponential_parametrizations_2_1D
    use EFTCAMB_double_exponential_parametrizations_1D
    use EFTCAMB_cosine_parametrizations_1D
    use EFTCAMB_axion_parametrizations_1D

    implicit none

    private

    public EFTCAMB_Horndeski

    type, extends(EFTCAMB_full_model) :: EFTCAMB_Horndeski

        !-----------------------------
        ! 时间网格设置：x = ln a
        !-----------------------------
        integer  :: interpolation_num_points = 800
        real(dl) :: x_initial = log(1.d-8)
        real(dl) :: x_final   = 0._dl
        real(dl) :: omegaLambda = 0._dl

        !-----------------------------
        ! 模型选择：Omega(a)、Lambda(a) 的参数化类型
        ! 在 read_model_selection 里从 ini 读入：
        !   Omega_model, Lambda_model
        ! 在 allocate_model_selection 里根据它们 allocate
        !   power_law_parametrization_1D / taylorseries / padeseries 等
        !-----------------------------
        integer :: Omega_model  = 0
        integer :: Lambda_model = 0

        ! 这里的父类型和 5e 保持一致
        class( parametrized_function_1D ), allocatable :: Omega
        class( parametrized_function_1D ), allocatable :: Lambda

        !-----------------------------
        ! 解出来的 EFT 背景函数，在 x 网格上等间隔采样：
        !   EFTc      : c(a) a^2 / m0^2
        !   EFTLambda : \Lambda(a) a^2 / m0^2
        ! 这两个是在 EFTCAMBHorndeskiSolveBackgroundEquations 的 output()里被填充的，在 InitBackground 里用 initialize() 设置网格.
        !-----------------------------
        type(equispaced_linear_interpolate_function_1D) :: EFTc
        type(equispaced_linear_interpolate_function_1D) :: EFTLambda

        !（可选）将来如果要在整套代码里用到 rho_DE, w_DE, 可以把它们加进来并在 output() 里填值：
        !type(equispaced_linear_interpolate_function_1D) :: rhoDE
        !type(equispaced_linear_interpolate_function_1D) :: wDE

        !（可选）如果你想保存整条 H^2(x) 解，也可以保留这个数组；
        ! 当前的 SolveBackgroundEquations 没有写它，可以之后在 output() 里加上：
        ! self%H2(ind) = H2
        real(dl), allocatable :: H2(:)
        real(dl), allocatable :: wDE(:)
        real(dl), allocatable :: OmegaDE(:)
        real(dl), allocatable :: Hphys(:)
        real(dl), allocatable :: rhoDE(:)
        real(dl), allocatable :: pDE(:)
        real(dl), allocatable :: cEFT(:)


    contains

        !-----------------------------
        ! 配置与参数初始化相关
        !-----------------------------
        procedure :: read_model_selection            => EFTCAMBHorndeskiReadModelSelection
        procedure :: allocate_model_selection        => EFTCAMBHorndeskiAllocateModelSelection
        procedure :: init_model_parameters           => EFTCAMBHorndeskiInitModelParameters
        procedure :: init_model_parameters_from_file => EFTCAMBHorndeskiInitModelParametersFromFile

        procedure :: compute_param_number            => EFTCAMBHorndeskiComputeParametersNumber
        procedure :: feedback                        => EFTCAMBHorndeskiFeedback
        procedure :: parameter_names                 => EFTCAMBHorndeskiParameterNames
        procedure :: parameter_names_latex           => EFTCAMBHorndeskiParameterNamesLatex
        procedure :: parameter_values                => EFTCAMBHorndeskiParameterValues

        !-----------------------------
        ! 背景 / EFT 函数接口（供 EFTCAMB 主框架调用）
        !-----------------------------
        procedure :: compute_background_EFT_functions  => EFTCAMBHorndeskiBackgroundEFTFunctions
        procedure :: compute_secondorder_EFT_functions => EFTCAMBHorndeskiSecondOrderEFTFunctions
        procedure :: compute_adotoa                    => EFTCAMBHorndeskiComputeAdotoa
        procedure :: compute_H_derivs                  => EFTCAMBHorndeskiComputeHubbleDer

        ! ----------------------------
        procedure :: initialize_background           => EFTCAMBHorndeskiInitBackground               !< subroutine that initializes the background of 5e.
        procedure :: solve_background_equations      => EFTCAMBHorndeskiSolveBackgroundEquations     !< subroutine that solves the 5e background equations.
        procedure :: find_initial_conditions         => EFTCAMBHorndeskiFindInitialConditions        !< subroutine that solves the background equations several time to determine the values of the initial conditions.

    end type EFTCAMB_Horndeski


contains

    subroutine EFTCAMBHorndeskiReadModelSelection( self, Ini )

        implicit none

        class(EFTCAMB_Horndeski) :: self
        type(TIniFile)          :: Ini

        self%Omega_model = Ini_Read_Int_File(Ini, 'Omega_model', 0)
        self%Lambda_model = Ini_Read_Int_File(Ini, 'Lambda_model', 0)

        ! self%num_points = Ini_Read_Int_File(Ini, 'RandomHorndeski_num_points', self%num_points)
        ! self%x_initial  = Ini_Read_Double_File(Ini, 'RandomHorndeski_loga_initial', self%x_initial)
        ! self%x_final    = Ini_Read_Double_File(Ini, 'RandomHorndeski_loga_final', self%x_final)
        ! self%initial_condition_mode = Ini_Read_Int_File(Ini, 'RandomHorndeski_initial_mode', self%initial_condition_mode)
    end subroutine EFTCAMBHorndeskiReadModelSelection

    subroutine EFTCAMBHorndeskiAllocateModelSelection( self, Ini )

        implicit none

        class(EFTCAMB_Horndeski) :: self
        type(TIniFile)          :: Ini

        if ( allocated(self%Omega) ) deallocate(self%Omega)
        select case ( self%Omega_model)
            case(1)
                allocate( power_law_parametrization_1D::self%Omega )
            case(2)
                allocate( taylorseries_parametrization_1D::self%Omega )
            case(3)
                allocate( padeseries_parametrization_1D::self%Omega )
            case default
                write(*,'(a,I3)') 'No model corresponding to Omega_model =', self%Omega_model
                write(*,'(a)')    'Please select an appropriate model.'
                stop
        end select

        if ( allocated(self%Lambda) ) deallocate(self%Lambda)
        select case ( self%Lambda_model)
            case(1)
                allocate( power_law_parametrization_1D::self%Lambda )
            case(2)
                allocate( taylorseries_parametrization_1D::self%Lambda )
            case(3)
                allocate( padeseries_parametrization_1D::self%Lambda )
            case default
                write(*,'(a,I3)') 'No model corresponding to Lambda_model =', self%Lambda_model
                write(*,'(a)')    'Please select an appropriate model.'
                stop
        end select

        ! ------ 给函数起名（可选但很有用，用于 feedback） ------
        call self%Omega%set_name(   'Omega',   '\Omega(a)'   )
        call self%Lambda%set_name( 'Lambda', '\Lambda(a)' )

        ! ------ 从 ini 读函数结构设置：比如节点个数、loga 区间等 ------
        call self%Omega%init_func_from_file(   Ini )
        call self%Lambda%init_func_from_file(  Ini )

    end subroutine EFTCAMBHorndeskiAllocateModelSelection

    !> 从 array 中按顺序读入 Ω(a)、Λ(a) 的参数值
    subroutine EFTCAMBHorndeskiInitModelParameters( self, array )

        implicit none

        class(EFTCAMB_Horndeski) :: self
        real(dl), dimension(self%parameter_number), intent(in) :: array

        real(dl), allocatable :: temp(:)
        integer :: num_params_function
        integer :: num_params_temp
        integer :: i

        num_params_temp = 1   ! 当前 array 的读指针，从 1 开始

        ! ---------- Omega(a) ----------
        if ( self%Omega_model > 0 ) then
            num_params_function = self%Omega%parameter_number
            allocate( temp(num_params_function) )

            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp + 1
            end do

            call self%Omega%init_parameters( temp )
            deallocate( temp )
        end if

        ! ---------- Lambda(a) ----------
        if ( self%Lambda_model > 0 ) then
            num_params_function = self%Lambda%parameter_number
            allocate( temp(num_params_function) )

            do i = 1, num_params_function
                temp(i)         = array(num_params_temp)
                num_params_temp = num_params_temp + 1
            end do

            call self%Lambda%init_parameters( temp )
            deallocate( temp )
        end if

        ! ---------- 额外常数 omegaLambda ----------
        self%omegaLambda    = array(num_params_temp)
        num_params_temp = num_params_temp + 1

        ! ---------- 安全检查 ----------
        if ( num_params_temp-1 /= self%parameter_number ) then
            write(*,*) 'In EFTCAMBHorndeskiInitModelParameters:'
            write(*,*) 'Length of array and self%parameter_number do not coincide.'
            write(*,*) 'Used elements    :', num_params_temp-1
            write(*,*) 'Expected elements:', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

    end subroutine EFTCAMBHorndeskiInitModelParameters


    !> 从 ini 文件读取本次运行要用的 Omega 和 Lambda 的参数值
    subroutine EFTCAMBHorndeskiInitModelParametersFromFile( self, Ini )

        implicit none

        class(EFTCAMB_Horndeski) :: self
        type(TIniFile)           :: Ini

        ! Omega(a) 的参数（如果模型号非 0）
        if ( self%Omega_model > 0 ) then
            call self%Omega%init_from_file( Ini )
        end if

        ! Lambda(a) 的参数
        if ( self%Lambda_model > 0 ) then
            call self%Lambda%init_from_file( Ini )
        end if

       ! 读入额外的常数 omegaLambda
        self%omegaLambda = Ini_Read_Double_File( Ini, 'omegaLambda', 0.7_dl)

    end subroutine EFTCAMBHorndeskiInitModelParametersFromFile

    !> 计算 Horndeski 模型总共有多少个参数
    subroutine EFTCAMBHorndeskiComputeParametersNumber( self )

        implicit none

        class(EFTCAMB_Horndeski) :: self

        self%parameter_number = 0

        if ( self%Omega_model  > 0 ) &
            self%parameter_number = self%parameter_number + self%Omega%parameter_number

        if ( self%Lambda_model > 0 ) &
            self%parameter_number = self%parameter_number + self%Lambda%parameter_number

        self%parameter_number = self%parameter_number + 1

    end subroutine EFTCAMBHorndeskiComputeParametersNumber

    !> Subroutine that prints on the screen feedback information about the Horndeski model.
    subroutine EFTCAMBHorndeskiFeedback( self, print_params )

        implicit none

        class(EFTCAMB_Horndeski) :: self         !< the base class
        logical, optional        :: print_params !< optional flag that decides whether to print
                                                !! numerical values of the parameters.

        !----------------------------
        ! 一些总体信息
        !----------------------------
        write(*,*)
        write(*,'(a,a)')   '   Model               =  ', self%name
        write(*,'(a,I3)')  '   Number of params    =',  self%parameter_number

        !----------------------------
        ! Horndeski EFT 函数的 model 选择信息
        !----------------------------
        write(*,*)
        if ( self%Omega_model /= 0 ) then
            write(*,'(a,I3)') '   Omega_model         =', self%Omega_model
        end if
        if ( self%Lambda_model /= 0 ) then
            write(*,'(a,I3)') '   Lambda_model        =', self%Lambda_model
        end if

        !----------------------------
        ! 具体函数对象的 feedback
        !----------------------------
        write(*,*)
        if ( self%Omega_model > 0 ) then
            call self%Omega%feedback( print_params )
        end if

        if ( self%Lambda_model > 0 ) then
            call self%Lambda%feedback( print_params )
        end if

    end subroutine EFTCAMBHorndeskiFeedback


    subroutine EFTCAMBHorndeskiParameterNames( self, i, name )

        implicit none

        class(EFTCAMB_Horndeski) :: self
        integer     , intent(in) :: i
        character(*), intent(out):: name

        integer :: NOmega, NLambda
        integer :: j

        ! 累积数量：前 NOmega 个是 Omega，之后的是 Lambda
        NOmega  = 0
        NLambda = 0

        if ( self%Omega_model  > 0 ) NOmega  = self%Omega%parameter_number
        if ( self%Lambda_model > 0 ) NLambda = NOmega + self%Lambda%parameter_number

        ! 合法性检查
        if ( i <= 0 .or. i > self%parameter_number ) then
            write(*,'(a,I3)') 'No parameter corresponding to: ', i
            write(*,'(a,I3)') 'Total number of parameters is: ', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

        ! omegade 固定放在最后一个槽
        if ( i == self%parameter_number ) then
            name = 'omegaLambda'
            return
        end if

        ! 来自 Omega(a)
        if ( self%Omega_model > 0 .and. i <= NOmega ) then
            do j = 1, self%Omega%parameter_number
                if ( i == j ) call self%Omega%parameter_names( j, name )
            end do
            return
        end if

        ! 来自 Lambda(a)
        if ( self%Lambda_model > 0 .and. i <= NLambda ) then
            do j = 1, self%Lambda%parameter_number
                if ( i-NOmega == j ) call self%Lambda%parameter_names( j, name )
            end do
            return
        end if

    end subroutine EFTCAMBHorndeskiParameterNames

    !> Subroutine that returns the LaTeX name of the i-th parameter of the Horndeski model.
    subroutine EFTCAMBHorndeskiParameterNamesLatex( self, i, latexname )

        implicit none

        class(EFTCAMB_Horndeski) :: self        !< the base class
        integer     , intent(in) :: i           !< the index of the parameter
        character(*), intent(out):: latexname   !< LaTeX name of the i-th parameter

        integer :: NOmega, NLambda

        !----------------------------------------------------------
        ! 1) 累积各个函数的参数数目:
        !    约定：前 NOmega 个属于 Omega，之后 NLambda-NOmega 个属于 Lambda
        !----------------------------------------------------------
        NOmega  = 0
        NLambda = 0

        if ( self%Omega_model  > 0 ) NOmega  = self%Omega%parameter_number
        if ( self%Lambda_model > 0 ) NLambda = NOmega + self%Lambda%parameter_number

        !----------------------------------------------------------
        ! 2) 检查合法性
        !----------------------------------------------------------
        if ( i <= 0 .or. i > self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if


        if ( i == self%parameter_number ) then
            latexname = '\Omega_{\rm Lambda}'
            return
        end if


        !----------------------------------------------------------
        ! 3) 来自 Omega(a) 的参数
        !----------------------------------------------------------
        if ( self%Omega_model > 0 .and. i <= NOmega ) then
            call self%Omega%parameter_names_latex( i, latexname )
            return
        end if

        !----------------------------------------------------------
        ! 4) 来自 Lambda(a) 的参数
        !----------------------------------------------------------
        if ( self%Lambda_model > 0 .and. i <= NLambda ) then
            call self%Lambda%parameter_names_latex( i - NOmega, latexname )
            return
        end if

        ! 理论上不应该走到这里，如果走到了说明上面的逻辑有遗漏:
        write(*,*) 'EFTCAMB error in EFTCAMBHorndeskiParameterNamesLatex: index mapping failed.'
        call MpiStop('EFTCAMB error')

    end subroutine EFTCAMBHorndeskiParameterNamesLatex


    !> Subroutine that returns the value of the i-th parameter of the Horndeski model.
    subroutine EFTCAMBHorndeskiParameterValues( self, i, value )

        implicit none

        class(EFTCAMB_Horndeski) :: self        !< the base class
        integer     , intent(in) :: i           !< the index of the parameter
        real(dl)    , intent(out):: value       !< the output value of the i-th parameter

        integer :: NOmega, NLambda

        !----------------------------------
        ! 1) 计算 Omega / Lambda 的累计个数
        !----------------------------------
        NOmega  = 0
        NLambda = 0

        if ( self%Omega_model  > 0 ) NOmega  = self%Omega%parameter_number
        if ( self%Lambda_model > 0 ) NLambda = NOmega + self%Lambda%parameter_number

        !----------------------------------
        ! 2) 合法性检查
        !----------------------------------
        if ( i <= 0 .or. i > self%parameter_number ) then
            write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
            write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
            call MpiStop('EFTCAMB error')
        end if

        if ( i == self%parameter_number ) then
            value = self%omegaLambda
            return
        end if

        !----------------------------------
        ! 3) 来自 Omega(a) 的参数
        !----------------------------------
        if ( self%Omega_model > 0 .and. i <= NOmega ) then
            call self%Omega%parameter_value( i, value )
            return
        end if

        !----------------------------------
        ! 4) 来自 Lambda(a) 的参数
        !----------------------------------
        if ( self%Lambda_model > 0 .and. i <= NLambda ) then
            call self%Lambda%parameter_value( i - NOmega, value )
            return
        end if

        ! 理论上到不了这里
        write(*,*) 'EFTCAMB error in EFTCAMBHorndeskiParameterValues: index mapping failed.'
        call MpiStop('EFTCAMB error')

    end subroutine EFTCAMBHorndeskiParameterValues

    !> Subroutine that solves the Horndeski background equations driven by {Omega(a), Lambda(a)}. ini is today.
    subroutine EFTCAMBHorndeskiSolveBackgroundEquations( self, params_cache, H2_ini, H02, only_solve, success, outroot )

    use, intrinsic :: ieee_arithmetic
    implicit none

    class(EFTCAMB_Horndeski)                  :: self          !< the base class
    type(EFTCAMB_parameter_cache), intent(in) :: params_cache  !< cosmological background parameters
    real(dl),                   intent(in)    :: H2_ini        !< initial value of y = H^2 at a=1 (x=0)
    real(dl),                   intent(out)   :: H02           !< value of H^2 today (a=1)
    logical,        optional                  :: only_solve    !< if .true. only solves for H^2, no EFT tables
    logical,                   intent(out)    :: success       !< whether the integration completed successfully
    character(len=*), optional, intent(in)    :: outroot       !< root for debug output files

    ! ODE system size: here we only evolve y = H^2
    integer, parameter :: num_eq = 1
    integer, parameter :: unit_bg  = 33
    real(dl) :: y(num_eq), ydot(num_eq)

    ! odepack quantities:
    integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
    real(dl) :: rtol, atol, t1, t2
    real(dl), allocatable :: rwork(:)
    integer , allocatable :: iwork(:)

    ! background quantities shared between derivs / output:
    real(dl) :: a, a2
    real(dl) :: grhob_t, grhoc_t, grhor_t, grhog_t
    real(dl) :: grhonu_tot, gpinu_tot
    real(dl) :: grhonu, gpinu, grhormass_t
    real(dl) :: grho_matter, gpres_matter
    real(dl) :: H2, Hdot
    real(dl) :: Omega, Omegap, Omegapp
    real(dl) :: Lambda_a2, Lambda_a2_prime, LLambda
    real(dl) :: Lambda, Lambda_prime, LLambda_prime, LLambda_dot
    real(dl) :: ca2_over_m0sq, cdot_a2_over_m0sq, cdot
    real(dl) :: dOmega_dN, dOmegap_dN, d_factor_dN
    real(dl) :: drho_ma2_dN, dLambda_a2_dN, dC_dN
    real(dl) :: omega_r_t, omega_m_t, omega_nu_t, omega_DE_t, omega_tot_t
    integer  :: nu_i

    ! A(a) 相关：方程前面的系数 A = 1 + Omega + 0.5 a Omega'
    real(dl) :: Acoef, Bcoef, Cterm
    real(dl) :: Acrit, Amin_grid, Atmp
    logical  :: hit_A_singularity, ok_A

    logical :: loc_only_solve

    ! write(*,*) 'DEBUG: entered EFTCAMBHorndeskiSolveBackgroundEquations'
    ! call flush(6)


    ! ---------------------------------------------------------------
    ! 0) digest the input flags
    ! ---------------------------------------------------------------
    if ( present(only_solve) ) then
        loc_only_solve = only_solve
    else
        loc_only_solve = .False.
    end if

    success            = .True.
    H02                = H2_ini
    Acrit              = 1.d-8    ! 运行时兜底阈值：|A| < Acrit 认为接近奇点
    Amin_grid          = 1.d-3    ! 预检查阈值：网格上 |A| 必须都 > Amin_grid
    hit_A_singularity  = .False.

    if (DebugEFTCAMB) then
        if ( present(outroot) .and. len_trim(outroot) > 0 ) then
            open(unit=unit_bg, file=trim(outroot)//'Horndeski_solution.dat', status='replace', &
                 action='write', form='formatted')
        else
            open(unit=unit_bg, file='Horndeski_solution.dat', status='replace', &
                 action='write', form='formatted')
        end if
        write(unit_bg,'(A)') '# x  a  z  H2  Hdot  grho_m  gpres_m  Omega  Omegap  Omegapp  ' // &
                             'ca2_over_m0sq  Lambda_a2  omega_r  omega_m  omega_nu  omega_DE  omega_tot  ' // &
                             'wDE  OmegaDE  rho_m/(3H2)  P_m/(3H2)  Lambda_a2/(3H2)  Acoef  Bcoef  Cterm  ' // &
                             'rhoDE_real  c_real'
    end if

    ! ---------------------------------------------------------------
    ! 0.5) 预检查：在整个背景 x-grid 上扫描 A(a)，排除明显不好的 Omega(a)
    !
    ! 要求：对所有网格点 i，
    !   |A(a_i)| > Amin_grid
    !
    ! 也可以在这里加上 A(a_i) > 0 的条件（ghost-free），看后面想不想加先验。
    ! ---------------------------------------------------------------
    ok_A = .True.

    ! do i = 1, self%EFTc%num_points
    !     a      = exp( self%EFTc%x(i) )
    !     Omega  = self%Omega%value(a)
    !     Omegap = self%Omega%first_derivative(a)

    !     Atmp = 1._dl + Omega + 0.5_dl * a * Omegap

    !     if ( abs(Atmp) < Amin_grid ) then
    !         ok_A = .False.
    !         exit
    !     end if
    !     ! 如果想强制 A>0，可以改成：
    !     ! if ( Atmp <= 0._dl .or. abs(Atmp) < Amin_grid ) then ...
    ! end do

    ! if ( .not. ok_A ) then
    !     if ( DebugEFTCAMB ) then
    !         write(*,'(a,es15.6)') 'Horndeski: A(a) too small on background grid, aborting. Atmp = ', Atmp
    !         close(unit_bg)
    !     end if
    !     success = .False.
    !     return
    ! end if

    ! ---------------------------------------------------------------
    ! 1) set initial conditions
    !
    ! 约定：背景 x-grid 满足
    !   self%EFTc%x( self%EFTc%num_points ) = 0.d0   <=>   a = 1 today
    ! 所以这里把初值放在 a=1（今天）这个点上：
    !   y(1) = H^2(a=1) = H2_ini
    ! ---------------------------------------------------------------
    t1   = self%EFTc%x( self%EFTc%num_points )   ! 这个应该是 0.d0
    y(1) = H2_ini

    ! ---------------------------------------------------------------
    ! 2) Initialize DLSODA
    ! ---------------------------------------------------------------
    itol = 1
    rtol = 1.d-10
    atol = 1.d-14

    itask  = 1
    istate = 1
    iopt   = 1

    LRN = 20 + 16*num_eq
    LRS = 22 + 9*num_eq + num_eq**2
    LRW = max(LRN, LRS)

    LIS = 20 + num_eq
    LIN = 20
    LIW = max(LIS, LIN)

    allocate(rwork(LRW))
    allocate(iwork(LIW))

    ! optional lsoda input:
    rwork(5) = 0._dl  ! initial step size (0 => chosen by solver)
    rwork(6) = 0._dl  ! max step size (0 => infinite)
    rwork(7) = 0._dl  ! min step size (0 => 0)

    iwork(5) = 0      ! no extra printing
    iwork(6) = 1000   ! max internal steps per call
    iwork(7) = 0      ! max warning messages
    iwork(8) = 0      ! max order Adams
    iwork(9) = 0      ! max order BDF

    call XSETF(0)     ! suppress odepack printing

    JacobianMode = 1  ! 1 = full Jacobian provided

    ! ---------------------------------------------------------------
    ! 3) 在 a=1 处（x=0）先写入 EFT 表的最后一个点
    !    ind = num_points  对应 a=1
    ! ---------------------------------------------------------------
    if ( .not. loc_only_solve ) then
        call output( num_eq, self%EFTc%num_points, t1, y )
    end if

    ! ---------------------------------------------------------------
    ! 4) 从今天 (a=1, x=0) 向过去积分：
    !    i = num_points, num_points-1, ..., 2
    !    每一步从 x(i) 到 x(i-1)
    !
    ! 网格要求：self%EFTc%x(:) 单调递增，最后一个点 = 0
    ! 例如：x(1) = ln(a_min), x(num_points) = 0
    ! ---------------------------------------------------------------
    do i = self%EFTc%num_points, 2, -1

        t1 = self%EFTc%x(i)
        t2 = self%EFTc%x(i-1)

        call DLSODA( derivs, num_eq, y, t1, t2, itol, rtol, atol, &
                     itask, istate, iopt, rwork, LRW, iwork, LIW, &
                     jacobian, JacobianMode )

        ! 检查 LSODA 状态
        if ( istate < 0 ) then
            if ( istate == -1 ) then
                ! -1: 需要更多步，重置后继续
                if ( DebugEFTCAMB ) then
                    write(*,'(a,i4,2es15.6)') &
                         'Horndeski LSODA istate=-1, t1, t2 = ', istate, t1, t2
                end if
                istate = 1
            else
                ! 其它负值：真正的错误
                if ( DebugEFTCAMB ) then
                    write(*,'(a,i4,2es15.6)') &
                         'Horndeski LSODA ERROR, istate, t1, t2 = ', istate, t1, t2
                    close(unit_bg)
                end if
                success = .False.
                return
            end if
        end if

        ! 再检查一下是否在 derivs 中触发了 A 奇点标记
        if ( hit_A_singularity ) then
            if ( DebugEFTCAMB ) then
                write(*,'(a,3es15.6)') 'Horndeski: Acoef too small during integration, t1,t2,A = ', t1, t2, Acoef
                close(unit_bg)
            end if
            success = .False.
            return
        end if

        ! 写入更早一个网格点 i-1 的 EFT 表
        if ( .not. loc_only_solve ) then
            call output( num_eq, i-1, t2, y )
        end if

    end do

    ! ---------------------------------------------------------------
    ! 5) 全部积分成功结束
    !
    ! 此时：
    !   - self%EFTc%y(1:num_points), self%EFTLambda%y(1:num_points),
    !     self%H2(:), self%wDE(:), self%OmegaDE(:) 等表都已填满；
    !   - H02 保持等于输入的 H2_ini（今天的 H^2）。
    !   - A(a) 在整个网格上都 > Amin_grid，且运行过程中未触发 Acrit。
    ! ---------------------------------------------------------------
    if (.not. ieee_is_finite(H02)) then
        error stop 'NaN/Inf in H02 in EFTCAMBHorndeskiSolveBackgroundEquations'
    end if

    success = .True.

    if (DebugEFTCAMB) then
        close(unit_bg)
    end if

    return

contains


        ! ================================================================
        ! 下面为 derivs / jacobian / output 三个内部子程序
        !
        ! 只要它们用了上面声明的共享变量（a, H2, Omega, Lambda_a2, ...）
        ! 并且假设 x = ln a，y(1) = H^2，逻辑就是一致的。
        ! ================================================================


        ! =================================================================
        !> RHS: dy/dx for Horndeski background, x = ln a, y = H^2
        subroutine derivs( num_eq, x, y, ydot )

            implicit none

            integer , intent(in)                     :: num_eq
            real(dl), intent(in)                     :: x
            real(dl), intent(in) , dimension(num_eq) :: y
            real(dl), intent(out), dimension(num_eq) :: ydot

            real(dl) :: Pm_a2_over_m0sq

            ! 0) x -> a
            a  = exp(x)
            a2 = a*a

            ! 1) compute background densities (完全沿用 5e 的写法)
            grhob_t = params_cache%grhob/a        ! bayron
            grhoc_t = params_cache%grhoc/a        ! cold dark matter
            grhor_t = params_cache%grhornomass/a2 ! massive nu\eutrinos
            grhog_t = params_cache%grhog/a2       ! radiation

            grhonu_tot = 0._dl
            gpinu_tot  = 0._dl
            if ( params_cache%Num_Nu_Massive /= 0 ) then
                do nu_i = 1, params_cache%Nu_mass_eigenstates
                    grhonu      = 0._dl
                    gpinu       = 0._dl
                    grhormass_t = params_cache%grhormass(nu_i)/a2
                    call params_cache%Nu_background( a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                    grhonu_tot = grhonu_tot + grhormass_t*grhonu
                    gpinu_tot  = gpinu_tot  + grhormass_t*gpinu
                end do
            end if

            grho_matter  = grhob_t + grhoc_t + grhor_t + grhog_t + grhonu_tot
            gpres_matter = gpinu_tot + (grhog_t + grhor_t)/3._dl ! P_m a^2 / m0^2

            ! 2) EFT function Omega(a) and derivatives w.r.t a
            Omega   = self%Omega%value(a)
            Omegap  = self%Omega%first_derivative(a)   ! d Omega / d a
            Omegapp = self%Omega%second_derivative(a)  ! d^2 Omega / d a^2

            ! 3) Lambda(a): 用组合形式 Lambda(a)a^2/m0^2 存在 self%Lambda%value 里
            Lambda       = self%Lambda%value(a)
            Lambda_prime = self%Lambda%first_derivative(a)  ! d/da(Lambda a^2 / m0^2)，供 output 用
            LLambda = -3._dl * self%omegaLambda * H2_ini * (1 + Lambda) * a2 ! Lambda a^2 / m0^2

            Lambda_a2 = (kappa / (c**2)) * Lambda * a2 * Mpc**2
            Lambda_a2_prime = (kappa / (c**2)) * Mpc**2 * ( 2._dl*a*Lambda + a2*Lambda_prime )

            ! 4) 系数 A,B,C in eq. (2)
            Acoef = 1._dl + Omega + 0.5_dl * a * Omegap
            Bcoef = 1._dl + Omega + 2._dl*a*Omegap + a2*Omegapp
            Cterm = gpres_matter + LLambda              ! C = P_m a^2/m0^2 + Lambda a^2/m0^2

            H2 = y(1)

            ! eq.(2): A dH^2/dN + B H^2 + C = 0
            ! --------- A 保护逻辑 ----------
            if ( abs(Acoef) < Acrit ) then
                hit_A_singularity = .True.
                call MpiStop('Horndeski: |Acoef| < Acrit, aborting background integration')
            else
                ydot(1) = -( Bcoef*H2 + Cterm ) / Acoef
                ! conformal-time derivative of \mathcal{H}: 这里 Hdot = 0.5 dH^2/dN
                Hdot    = 0.5_dl*ydot(1)
            end if

        end subroutine derivs


            ! =================================================================
        !> Jacobian matrix \partial(dy/dx)/\partial y, needed by DLSODA
        subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

            implicit none

            integer                            :: num_eq, ml, mu, nrowpd
            real(dl)                           :: x
            real(dl), dimension(num_eq)        :: y
            real(dl), dimension(nrowpd,num_eq) :: pd

            real(dl) :: dummy(1)

            ! 在当前 (x,y) 调用一下 derivs，确保 Acoef,Bcoef,H2 等已经更新
            call derivs( num_eq, x, y, dummy )

            ! 如果已经标记了 A 奇点，那这个 Jacobian 基本也没用了, 置零就行, 上层会根据 hit_A_singularity 直接退出。
            if ( hit_A_singularity ) then
                pd(1,1) = 0._dl
                return
            end if

            ! 正常情况：dy/dx = -(B y + C)/A => \partial (dy/dx)/\partial y = -(B)/A
            if ( abs(Acoef) < 1.d-20 ) then
                pd(1,1) = 0._dl
            else
                pd(1,1) = -Bcoef / Acoef
            end if

        end subroutine jacobian


        ! =================================================================
        !> Take the solution y=H^2 and compute EFT functions c(a), Lambda(a)
        !! and some auxiliary background quantities.
        subroutine output( num_eq, ind, x, y )

        implicit none

        integer , intent(in)                     :: num_eq
        integer , intent(in)                     :: ind
        real(dl), intent(in)                     :: x
        real(dl), intent(in) , dimension(num_eq) :: y

        logical :: is_open

        ! DE 相关量（组合 + 真实）
        real(dl) :: rhoDE_hat, pDE_hat
        real(dl) :: wDE_here, OmegaDE_here
        real(dl) :: rhoDE_real, pDE_real
        real(dl) :: c_real
        real(dl) :: c_real_dot
        real(dl) :: pref_m0_over_a2
        real(dl) :: H_phys, H
        real(dl) :: Lambda_dot

        ! Make sure all background quantities are set at this (x,y):
        call derivs( num_eq, x, y, ydot )

        pref_m0_over_a2 = c**2 / (kappa * a2)

        a  = exp(x)
        a2 = a*a
        H2 = y(1)          ! \mathcal{H}^2

        ! ---- 1) EFT functions: c and Lambda ---------------------------------

        ! 裸的 Lambda(a) 及其对 a 的导数
        Lambda       = self%Lambda%value(a)
        Lambda_prime = self%Lambda%first_derivative(a)
        LLambda = -3._dl * self%omegaLambda * H2_ini * (1 + Lambda) * a2
        LLambda_prime = -3._dl * self%omegaLambda * H2_ini * Lambda_prime * a2

        ! Lambda a^2 / m0^2 = (kappa/c^2) Lambda a^2 * Mpc^2
        ! Lambda_a2 = (kappa / (c**2)) * Lambda * a2 * Mpc**2

        ! d/da [ (kappa/c^2) Lambda a^2 Mpc^2 ]
        ! = (kappa/c^2) Mpc^2 (2 a Lambda + a^2 Lambda')
        ! Lambda_a2_prime = (kappa / (c**2)) * Mpc**2 * ( a2*Lambda_prime )

        ! c(a) via Raveri eq. (2):
        ! ca^2/m0^2 = 3/2 (1+Omega+aOmega') H^2 - 1/2 rho_m a^2/m0^2 + 1/2 Lambda a^2/m0^2
        ca2_over_m0sq = 1.5_dl * ( 1._dl + Omega + a*Omegap ) * H2 &
                        - 0.5_dl * grho_matter &
                        + 0.5_dl * LLambda

        ! 存 EFT 表：这里仍然是组合量 c a^2/m0^2 和 Lambda a^2/m0^2
        self%EFTc%y(ind)      = ca2_over_m0sq
        self%EFTLambda%y(ind) = LLambda

        ! ---- 2) time derivatives c_dot and Lambda_dot -----------------------

        if ( H2 > 0._dl ) then

            H = sqrt(H2)

            ! 以下求cdot表达式, 实际上是cdot a^2 / m0^2
            cdot = 1.5_dl * H * (a * H2 * (2 * Omegap + a * Omegapp) + 2 * (1 + Omega + a * Omegap) * (Hdot - H2)) &
                    + 0.5_dl * H * a * LLambda_prime &
                    + 1.5_dl * H * (grho_matter + gpres_matter)

            LLambda_dot = -3._dl * self%omegaLambda * H2_ini * Lambda_prime * a2 * a * H

            ! Lambdadot a^2 / m0^2
            self%EFTc%yp(ind)      = cdot
            self%EFTLambda%yp(ind) = LLambda_dot

        else
            self%EFTc%yp(ind)      = double_NaN
            self%EFTLambda%yp(ind) = double_NaN
        end if

        ! ---- 3) standard Omega_i for debugging ----------------------------------
        omega_r_t  = (grhog_t + grhor_t)/(3._dl*H2)
        omega_m_t  = (grhob_t + grhoc_t)/(3._dl*H2)
        omega_nu_t =  grhonu_tot/(3._dl*H2)
        omega_DE_t = (3._dl*H2 - (grhog_t + grhor_t + grhob_t + grhoc_t + grhonu_tot)) / (3._dl*H2)
        omega_tot_t= omega_r_t + omega_m_t + omega_nu_t + omega_DE_t

        ! ---- 4) effective DE density, pressure, w_DE, rho_DE (真实) ----------

        ! 组合密度:  rho_DE a^2 / m0^2 = 3 H^2 - rho_matter a^2 / m0^2
        rhoDE_hat = 3._dl*H2 - grho_matter

        ! 组合压强:  P_DE a^2 / m0^2 = -2 \dot{H} - H^2 - P_matter a^2 / m0^2
        ! 注意：Hdot 在 derivs 里定义为 Hdot = dH/d eta, 所有的H都对应文中\mathcal{H}
        pDE_hat   = -2._dl*Hdot - H2 - gpres_matter

        ! m0^2 / a^2 = c^2 / (kappa * Mpc^2) / a^2
        ! 真实暗能量密度和压强
        rhoDE_real = pref_m0_over_a2 * rhoDE_hat
        pDE_real   = pref_m0_over_a2 * pDE_hat


        ! --- 先算 Omega_DE，用组合量（始终应该是 well-defined，只要 H2>0） ---
        if ( H2 > 0._dl ) then
            OmegaDE_here = rhoDE_hat / (3._dl*H2)
        else
            OmegaDE_here = double_NaN
        end if

        ! --- 再算 w_DE，只在 rho_DE 不太接近 0 时给出数值 ---
        if ( H2 > 0._dl .and. abs(OmegaDE_here) > 1.d-4 ) then
            ! 这里用 rhoDE_hat 做阈值，并且按总密度 3H^2 缩放一下
            wDE_here = pDE_hat / rhoDE_hat
        else
            wDE_here = double_NaN
        end if

        H_phys = sqrt(H2) / a

        if ( allocated(self%Hphys) ) self%Hphys(ind) = H_phys
        if ( allocated(self%rhoDE) ) self%rhoDE(ind) = rhoDE_real
        if ( allocated(self%pDE) )   self%pDE(ind)   = pDE_real
        if ( allocated(self%cEFT) )  self%cEFT(ind)  = c_real
        ! 存到 type 里的数组中（保持原有接口）
        if ( allocated(self%H2) )      self%H2(ind)      = H2
        if ( allocated(self%wDE) )     self%wDE(ind)     = wDE_here
        if ( allocated(self%OmegaDE) ) self%OmegaDE(ind) = OmegaDE_here
        ! 如果在 type 里另外加了真实的 rhoDE/c 数组，以下可选：
        ! if ( allocated(self%rhoDE) ) self%rhoDE(ind) = rhoDE_real
        ! if ( allocated(self%cEFT) )  self%cEFT(ind)  = c_real

        ! ! ---- 5) debug 输出 ---------------------------------------------------
        ! if ( DebugEFTCAMB ) then
        !     write(*,'(200ES15.4E3)') x, a, 1._dl/a-1._dl, H2, Hdot, &
        !         grho_matter, gpres_matter, Omega, Omegap, Omegapp, &
        !         ca2_over_m0sq, Lambda_a2, omega_r_t, omega_m_t, omega_nu_t, &
        !         omega_DE_t, omega_tot_t, wDE_here, OmegaDE_here, &
        !         grho_matter/(3._dl*H2), gpres_matter/(3._dl*H2), &
        !         Lambda_a2/(3._dl*H2), Acoef, Bcoef, Cterm, &
        !         rhoDE_real, c_real
        ! end if

        if ( DebugEFTCAMB ) then
            write(unit_bg,'(200ES15.4E3)') x, a, 1._dl/a-1._dl, H2, Hdot, &
                 grho_matter, gpres_matter, Omega, Omegap, Omegapp, &
                 ca2_over_m0sq, Lambda_a2, omega_r_t, omega_m_t, omega_nu_t, &
                 omega_DE_t, omega_tot_t, wDE_here, OmegaDE_here, &
                 grho_matter/(3._dl*H2), gpres_matter/(3._dl*H2), &
                 Lambda_a2/(3._dl*H2), Acoef, Bcoef, Cterm, &
                 rhoDE_real, c_real
        end if

    end subroutine output

    end subroutine EFTCAMBHorndeskiSolveBackgroundEquations



    !> Subroutine that initializes the background of the Horndeski model, shooting method
!     subroutine EFTCAMBHorndeskiInitBackground( self, params_cache, feedback_level, success, outroot )

!     implicit none

!     class(EFTCAMB_Horndeski)                     :: self           !< the base class
!     type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< EFTCAMB parameter cache with cosmological params
!     integer                      , intent(in)    :: feedback_level !< 0=none; 1=some; 2=chatty
!     logical                      , intent(out)   :: success        !< whether background init succeeded
!     character(len=*), optional   , intent(in)    :: outroot        !< root for debug output files

!     real(dl) :: H2_ini       !< initial value of y = H^2 (or 𝓗^2) at x_initial = log(a_initial)
!     real(dl) :: H02          !< H^2 today (a=1), returned by the solver

!     !---------------------------------------------------------------
!     ! 1) some feedback
!     !---------------------------------------------------------------
!     if ( feedback_level > 1 ) then
!         write(*,'(a)') "***************************************************************"
!         write(*,'(a)') ' EFTCAMB Horndeski background solver'
!         write(*,'(a)')
!     end if

!     if ( DebugEFTCAMB .or. feedback_level > 2 ) then
!         call params_cache%print()
!     end if

!     !---------------------------------------------------------------
!     ! 2) initialize interpolating functions for EFT background
!     !    all EFT functions are sampled on the same x = log(a) grid
!     !---------------------------------------------------------------
!     call self%EFTLambda%initialize( self%interpolation_num_points, self%x_initial, self%x_final )
!     call self%EFTc%initialize     ( self%interpolation_num_points, self%x_initial, self%x_final )

!     !---------------------------------------------------------------
!     ! 2b) 分配用来存整条演化的数组: H2(a), w_DE(a), Omega_DE(a)      ! <-- NEW
!     !     网格点数与 EFTc 的插值表一致 (它在 initialize 里已经设好)
!     !---------------------------------------------------------------
!     if ( allocated(self%H2) )      deallocate(self%H2)             ! <-- NEW
!     if ( allocated(self%wDE) )     deallocate(self%wDE)            ! <-- NEW
!     if ( allocated(self%OmegaDE) ) deallocate(self%OmegaDE)        ! <-- NEW

!     allocate(self%H2     ( self%EFTc%num_points ))                 ! <-- NEW
!     allocate(self%wDE    ( self%EFTc%num_points ))                 ! <-- NEW
!     allocate(self%OmegaDE( self%EFTc%num_points ))                 ! <-- NEW

!     !（注意：真正往这三个数组里填值，是在 SolveBackgroundEquations()
!     !  的 output(...) 子程序里做的——那里我们已经约好用 ind 作为索引）

!     !---------------------------------------------------------------
!     ! 3) find initial condition H2_ini at x_initial by shooting so
!     !    that H02 matches params_cache%h0_Mpc**2
!     !---------------------------------------------------------------
!     call self%find_initial_conditions( params_cache, feedback_level, H2_ini, success )
!     if ( success ) then
!         if ( feedback_level > 1 ) write(*,'(a,E13.4)') '   initial condition H2_ini  = ', H2_ini
!     else
!         ! 如果找初始条件失败，直接返回，让上层决定怎么处理
!         return
!     end if

!     !---------------------------------------------------------------
!     ! 4) optional debug: print full background solution to files
!     !---------------------------------------------------------------
!     if ( DebugEFTCAMB ) then
!         if ( present(outroot) ) then
!             call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_1.dat', 33 )
!             call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_2.dat', 34 )
!             call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_3.dat', 35 )
!         else
!             call CreateTxtFile( 'background_Horndeski_solution_1.dat', 33 )
!             call CreateTxtFile( 'background_Horndeski_solution_2.dat', 34 )
!             call CreateTxtFile( 'background_Horndeski_solution_3.dat', 35 )
!         end if

!         ! 这里你也可以顺便把 w_DE 和 Omega_DE 写进 header，方便后面画图  ! <-- 可选
!         write(33,'(a)') '# x  a  z  H2  Hdot  grho_matter  gpres_matter  Omega  Omegap  Omegapp  wDE  OmegaDE'
!         write(34,'(a)') '# x  a  z  c_a2_over_m0sq  Lambda_a2_over_m0sq'
!         write(35,'(a)') '# x  a  z  omega_r_t  omega_m_t  omega_nu_t  omega_DE_t  omega_tot_t'

!         call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success )

!         close(33)
!         close(34)
!         close(35)

!         if ( .not. success ) return
!     end if

!     !---------------------------------------------------------------
!     ! 5) solve the background equations and store the EFT functions
!     !---------------------------------------------------------------
!     call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success )

! end subroutine EFTCAMBHorndeskiInitBackground


!     !> Subroutine that finds the initial condition H2_ini for the Horndeski background,
! !!   by requiring that the integrated H0^2 matches the desired value.
!     subroutine EFTCAMBHorndeskiFindInitialConditions( self, params_cache, feedback_level, H2_ini, success )

!     implicit none

!     class(EFTCAMB_Horndeski)                     :: self           !< the base class
!     type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< cosmological parameters
!     integer                      , intent(in)    :: feedback_level !< 0=none; 1=some; 2=chatty.
!     real(dl)                     , intent(out)   :: H2_ini         !< initial value of y = H^2 at x_initial
!     logical                      , intent(out)   :: success        !< whether the IC search succeeded

!     ! 射击 / 根寻找相关变量
!     real(dl) :: H0_wanted
!     real(dl) :: H2_guess, H2_min, H2_max
!     real(dl) :: H2_trial, H0_trial
!     real(dl) :: H2_prev,  H0_prev
!     real(dl) :: logH2_min, logH2_max, logH2
!     real(dl) :: a_ini, tol
!     integer  :: ind, num_search

!     logical  :: have_prev   ! 是否已经有一个“上一点” (H2_prev, H0_prev)
!     logical  :: solver_ok   ! helperH0 最近一次调用是否成功
!     logical  :: found_bracket

!     success       = .False.
!     solver_ok     = .False.
!     have_prev     = .False.
!     found_bracket = .False.

!     ! 绝对容差（H0^2 级别）
!     tol        = 1.d-12
!     num_search = 40

!     ! 初始 a
!     a_ini     = exp( self%x_initial )

!     ! 想要复现的今天的 H0^2
!     H0_wanted = params_cache%h0_Mpc**2

!     ! 一点输出，方便看逻辑
!     if (feedback_level > 2) then
!         write(*,'(a)') '-------------------------------------------------------'
!         write(*,'(a)') ' Horndeski IC finder: scanning H2_ini'
!         write(*,'(a,1p,e12.4)') '   x_initial = ', self%x_initial
!         write(*,'(a,1p,e12.4)') '   a_initial = ', a_ini
!         write(*,'(a,1p,e12.4)') '   H0_wanted^2 = ', H0_wanted
!     end if

!     !------------------------------------------------------------------
!     ! 1) 构造一个大致猜测 H2_guess
!     !    这里先用简单的近似：~ H0^2 / max(a_ini^3, 1e-8)
!     !    后面需要的话可以再用更物理的 GR 近似给个更精确的 guess。
!     !------------------------------------------------------------------
!     H2_guess = H0_wanted / max( a_ini**3, 1.d-8 )

!     if ( H2_guess <= 0._dl ) H2_guess = H0_wanted

!     ! 在 log10(H2) 空间上给一个很宽的扫描区间 [1e-4, 1e+4] * H2_guess
!     H2_min = H2_guess * 1.d-4
!     H2_max = H2_guess * 1.d+4
!     if (H2_min <= 0._dl) H2_min = H2_guess * 1.d-6

!     logH2_min = log10(H2_min)
!     logH2_max = log10(H2_max)

!     if (feedback_level > 2) then
!         write(*,'(a,1p,3e12.4)') '   H2_guess, H2_min, H2_max = ', H2_guess, H2_min, H2_max
!         write(*,'(a)') '   Scan table: index, H2_try, H0^2(H2_try), solver_ok'
!     end if

!     !------------------------------------------------------------------
!     ! 2) 在 [H2_min, H2_max] 上扫描，寻找 “成功积分且有符号变化” 的一对点
!     !------------------------------------------------------------------
!     do ind = 0, num_search

!         ! 扫描点：在 logH2 空间等间隔
!         logH2    = logH2_min + real(ind,dl)/real(num_search,dl) * (logH2_max - logH2_min)
!         H2_trial = 10._dl**logH2

!         ! 调用 helperH0 得到对应的 H0^2，并通过 solver_ok 知道这次积分有没有崩
!         H0_trial = helperH0( H2_trial )

!         if (feedback_level > 2) then
!             write(*,'(a,i3,1p,2e14.6,l2)') '   scan: ', ind, H2_trial, H0_trial, solver_ok
!         end if

!         ! 如果这次积分失败（DLSODA 报错或者 H02 始终为 0），直接跳过，不拿它做 bracket
!         if (.not. solver_ok) cycle

!         ! 如果这是第一个成功点，就先存起来当 “前一个点”
!         if (.not. have_prev) then
!             H2_prev  = H2_trial
!             H0_prev  = H0_trial
!             have_prev = .True.
!             cycle
!         end if

!         ! 此时已经有前一个成功点 (H2_prev, H0_prev)
!         ! 检查 f(H2) = H0^2(H2) - H0_wanted 是否在 [prev, trial] 之间变号
!         if ( (H0_prev - H0_wanted) * (H0_trial - H0_wanted) <= 0._dl ) then
!             ! 找到 bracket
!             H2_min        = H2_prev
!             H2_max        = H2_trial
!             found_bracket = .True.
!             exit
!         end if

!         ! 否则，更新前一个点，继续扫
!         H2_prev = H2_trial
!         H0_prev = H0_trial

!     end do

!     ! 如果整个扫描下来都没有找到 bracket
!     if (.not. found_bracket) then
!         if (feedback_level > 1) then
!             write(*,'(a)') '   EFTCAMBHorndeskiFindInitialConditions: bracket for H2_ini not found'
!             if (have_prev) then
!                 write(*,'(a,1p,2e12.4)') '   Last successful H2, H0^2 = ', H2_prev, H0_prev
!             else
!                 write(*,'(a)') '   (No successful background integration for any H2 in scan range)'
!             end if
!         end if
!         return
!     end if

!     if (feedback_level > 2) then
!         write(*,'(a,1p,2e12.4)') '   Bracket found between H2_min, H2_max = ', H2_min, H2_max
!     end if

!     !------------------------------------------------------------------
!     ! 3) 有了 [H2_min, H2_max]，用 zbrent 细化根：
!     !    helperH0(H2_ini) ≈ H0_wanted
!     !------------------------------------------------------------------
!     H2_ini = zbrent( helperH0, H2_min, H2_max, tol, H0_wanted, success )

!     if ( success ) then
!         if ( feedback_level > 2 ) then
!             write(*,'(a,1p,e12.4)') '   H0^2 given  = ', H0_wanted
!             write(*,'(a,1p,e12.4)') '   H0^2 found  = ', helperH0(H2_ini)
!             write(*,'(a,1p,e12.4)') '   H2_ini best = ', H2_ini
!         end if
!     else
!         if ( feedback_level > 1 ) then
!             write(*,'(a)') '   EFTCAMBHorndeskiFindInitialConditions: root finder failed'
!         end if
!     end if

! contains

!     !----------------------------------------------------------------------------
!     ! helperH0(H2_ini_local)：
!     !   给定某个 H2_ini = H^2(a_ini)，调用背景方程求解器，
!     !   得到今天 a=1 时的 H0^2，并把它返回。
!     !
!     !   每次调用会更新上层局部变量 solver_ok，告诉外层这次积分是否成功。
!     !   这样既能作为 zbrent 的函数，也能在扫描阶段判断可用/不可用。
!     !----------------------------------------------------------------------------
!     function helperH0( H2_ini_local ) result( H0_sq )

!         implicit none
!         real(dl), intent(in) :: H2_ini_local
!         real(dl)             :: H0_sq
!         logical              :: success_local

!         ! 先初始化
!         H0_sq        = 0._dl
!         success_local = .True.

!         call self%solve_background_equations( params_cache, H2_ini_local, &
!                                               H02        = H0_sq,         &
!                                               only_solve = .True.,        &
!                                               success    = success_local )

!         solver_ok = success_local

!         if (feedback_level > 3) then
!             write(*,'(a,1p,e12.4,a,1p,e12.4,l2)') &
!               '      helperH0: H2_ini=', H2_ini_local, '  H0^2=', H0_sq, solver_ok
!         end if

!     end function helperH0

! end subroutine EFTCAMBHorndeskiFindInitialConditions

    subroutine EFTCAMBHorndeskiInitBackground( self, params_cache, feedback_level, success, outroot )

    implicit none

    class(EFTCAMB_Horndeski)                     :: self           !< the base class
    type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< EFTCAMB parameter cache with cosmological params
    integer                      , intent(in)    :: feedback_level !< 0=none; 1=some; 2=chatty
    logical                      , intent(out)   :: success        !< whether background init succeeded
    character(len=*), optional   , intent(in)    :: outroot        !< root for debug output files

    ! 现在的含义：H2_ini 是今天 a=1 处的 H^2（或 𝓗^2），
    ! H02 只是为了和其它模型接口一致，实际上等于 H2_ini。
    real(dl) :: H2_ini       !< initial value of y = H^2 at a=1 (x=0)
    real(dl) :: H02          !< H^2 today (a=1), returned by the solver (≈ H2_ini)

    !---------------------------------------------------------------
    ! 1) some feedback
    !---------------------------------------------------------------
    if ( feedback_level > 1 ) then
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') ' EFTCAMB Horndeski background solver'
        write(*,'(a)')
    end if

    if ( DebugEFTCAMB .or. feedback_level > 2 ) then
        call params_cache%print()
    end if

    !---------------------------------------------------------------
    ! 2) initialize interpolating functions for EFT background
    !    all EFT functions are sampled on the same x = log(a) grid
    !    要求：x(1)=x_initial=ln(a_min)<0, x(num_points)=x_final=0 => a=1 today
    !---------------------------------------------------------------
    call self%EFTLambda%initialize( self%interpolation_num_points, self%x_initial, self%x_final )
    call self%EFTc%initialize     ( self%interpolation_num_points, self%x_initial, self%x_final )

    !---------------------------------------------------------------
    ! 2b) 分配用来存整条演化的数组: H2(a), w_DE(a), Omega_DE(a)
    !---------------------------------------------------------------
    if ( allocated(self%H2) )      deallocate(self%H2)
    if ( allocated(self%wDE) )     deallocate(self%wDE)
    if ( allocated(self%OmegaDE) ) deallocate(self%OmegaDE)
    if (allocated(self%Hphys)) deallocate(self%Hphys)
    if (allocated(self%rhoDE)) deallocate(self%rhoDE)
    if (allocated(self%pDE))   deallocate(self%pDE)
    if (allocated(self%cEFT))  deallocate(self%cEFT)

    allocate(self%H2     ( self%EFTc%num_points ))
    allocate(self%wDE    ( self%EFTc%num_points ))
    allocate(self%OmegaDE( self%EFTc%num_points ))
    allocate(self%Hphys( self%EFTc%num_points ))
    allocate(self%rhoDE( self%EFTc%num_points ))
    allocate(self%pDE( self%EFTc%num_points ))
    allocate(self%cEFT( self%EFTc%num_points ))

    !---------------------------------------------------------------
    ! 3) 现在的初始条件非常简单：
    !    直接要求今天 a=1 时 H^2 = (h0_Mpc)^2
    !---------------------------------------------------------------
    call self%find_initial_conditions( params_cache, feedback_level, H2_ini, success )
    if ( success ) then
        if ( feedback_level > 1 ) write(*,'(a,E13.4)') '   initial condition H2_ini(a=1) = ', H2_ini
    else
        ! 如果找初始条件失败，直接返回，让上层决定怎么处理
        return
    end if

    !---------------------------------------------------------------
    ! 4) optional debug: print full background solution to files
    !    注意：SolveBackgroundEquations 的 output(...) 里会在 DebugEFTCAMB 时写这几个单元
    !---------------------------------------------------------------
    if ( DebugEFTCAMB ) then
        if ( present(outroot) ) then
            call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_1.dat', 33 )
            call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_2.dat', 34 )
            call CreateTxtFile( trim(outroot)//'background_Horndeski_solution_3.dat', 35 )
        else
            call CreateTxtFile( 'background_Horndeski_solution_1.dat', 33 )
            call CreateTxtFile( 'background_Horndeski_solution_2.dat', 34 )
            call CreateTxtFile( 'background_Horndeski_solution_3.dat', 35 )
        end if

        write(33,'(a)') '# x  a  z  H2  Hdot  grho_matter  gpres_matter  Omega  Omegap  Omegapp  wDE  OmegaDE'
        write(34,'(a)') '# x  a  z  c_a2_over_m0sq  Lambda_a2_over_m0sq'
        write(35,'(a)') '# x  a  z  omega_r_t  omega_m_t  omega_nu_t  omega_DE_t  omega_tot_t'

        if ( present(outroot) ) then
            call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success, outroot = outroot )
        else
            call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success )
        end if

        close(33)
        close(34)
        close(35)

        if ( .not. success ) return
    end if

    !---------------------------------------------------------------
    ! 5) solve the background equations and store the EFT functions
    !    正式跑一遍背景，填满 EFT 插值表和 H2/w_DE/Omega_DE 的数组
    !---------------------------------------------------------------
    if ( present(outroot) ) then
        call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success, outroot = outroot )
    else
        call self%solve_background_equations( params_cache, H2_ini, H02, only_solve = .False., success = success )
    end if

    end subroutine EFTCAMBHorndeskiInitBackground

    !> Subroutine that sets the initial condition H2_ini for the Horndeski background.
    !!  在新的求解方案中，我们从今天 a=1 (x=0) 往过去积分，
    !!  所以 H2_ini 就是今天的 H0^2 = (h0_Mpc)^2。
    subroutine EFTCAMBHorndeskiFindInitialConditions( self, params_cache, feedback_level, H2_ini, success )

        implicit none

        class(EFTCAMB_Horndeski)                     :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< cosmological parameters
        integer                      , intent(in)    :: feedback_level !< 0=none; 1=some; 2=chatty.
        real(dl)                     , intent(out)   :: H2_ini         !< initial value of y = H^2 at a=1 (today)
        logical                      , intent(out)   :: success        !< whether the IC assignment succeeded

        real(dl) :: H0_wanted

        success    = .False.
        H0_wanted  = params_cache%h0_Mpc**2

        ! 基本 sanity check：H0^2 必须是正数
        if ( H0_wanted <= 0._dl ) then
            if ( feedback_level > 0 ) then
                write(*,'(a,1p,e12.4)') '  Horndeski IC finder: h0_Mpc^2 <= 0, got ', H0_wanted
            end if
            return
        end if

        ! 初始条件：今天 a=1 的 H^2 就是 H0^2
        H2_ini = H0_wanted
        success = .True.

        if ( feedback_level > 2 ) then
            write(*,'(a,1p,e12.4)') '  Horndeski IC finder: set H2_ini(a=1) = ', H2_ini
        end if

    end subroutine EFTCAMBHorndeskiFindInitialConditions



    !-----------------------------------------------------------------
    !> 背景 EFT 函数：给定尺度因子 a，返回 c(a), Lambda(a) 以及它们的 dot，
    !!  并计算 EFTOmega 及其导数。
    !-----------------------------------------------------------------
    subroutine EFTCAMBHorndeskiBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horndeski)                     :: self          !< the base class
        real(dl)              , intent(in)           :: a             !< input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< UNUSED here, but kept for interface
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< timestep cache with EFT values

        real(dl) :: x, mu
        real(dl) :: a_eff, a_min
        integer  :: ind

        !---------------------------------------------------------------
        ! 1) 保护 a 的范围：
        !    我们的背景和 EFT 表是在 [a_min, 1] 上采样的，
        !    如果 CAMB 在此区间外询问 EFT 函数，就用边界值。
        !---------------------------------------------------------------
        a_min = exp( self%x_initial )   ! x_initial = ln(a_min)

        a_eff = a
        if ( a_eff < a_min ) a_eff = a_min
        if ( a_eff > 1._dl ) a_eff = 1._dl

        x = log( a_eff )

        !---------------------------------------------------------------
        ! 2) 用 EFTc 的 sampled_function 来做预计算（找到区间和样条系数），
        !    所有 EFT 函数 (c, Lambda) 共用同一条 x 网格。
        !---------------------------------------------------------------
        call self%EFTc%precompute( a_eff, ind, mu )

        !---------------------------------------------------------------
        ! 3) 从预计算的样条表中插值出 c(a), Lambda(a) 以及它们的 dot
        !    注意：现在表里存的 c(a) 和 Λ(a)以及它们的共形时间导数，已乘 a^2/m0^2。
        !---------------------------------------------------------------
        eft_cache%EFTc         = self%EFTc%value(      a_eff, index=ind, coeff=mu )
        eft_cache%EFTLambda    = self%EFTLambda%value( a_eff, index=ind, coeff=mu )
        eft_cache%EFTcdot      = self%EFTc%first_derivative(      a_eff, index=ind, coeff=mu )
        eft_cache%EFTLambdadot = self%EFTLambda%first_derivative( a_eff, index=ind, coeff=mu )

        !---------------------------------------------------------------
        ! 4) Horndeski 路径下，Omega 本身就是解析给定的 EFT 函数，
        !    这里统一在 a_eff 上评估 Ω 及其导数，保证和背景解使用同一 a 
        !---------------------------------------------------------------
        eft_cache%EFTOmegaV   = self%Omega%value(            a_eff )
        eft_cache%EFTOmegaP   = self%Omega%first_derivative( a_eff )  ! dΩ/da
        eft_cache%EFTOmegaPP  = self%Omega%second_derivative( a_eff ) ! d²Ω/da²
        eft_cache%EFTOmegaPPP = self%Omega%third_derivative(  a_eff ) ! d³Ω/da³

    end subroutine EFTCAMBHorndeskiBackgroundEFTFunctions

    !-----------------------------------------------------------------
    !> Use the precomputed H2 background solution to set adotoa.
    !! Falls back to the full-map expression if H2 is not available.
    !-----------------------------------------------------------------
    subroutine EFTCAMBHorndeskiComputeAdotoa( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horndeski)                     :: self          !< the base class
        real(dl)              , intent(in)           :: a             !< input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< EFT parameter cache
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< timestep cache

        real(dl) :: a_eff, a_min, x_eff, H2_val
        real(dl) :: x1, x2, mu, temp
        integer  :: ind

        a_min = exp( self%x_initial )
        a_eff = a
        if ( a_eff < a_min ) a_eff = a_min
        if ( a_eff > 1._dl ) a_eff = 1._dl

        if ( .not. allocated(self%H2) ) then
            call MpiStop('Horndeski: H2 array not allocated in compute_adotoa')
        end if

        x_eff = log( a_eff )
        if ( x_eff <= self%EFTc%x_initial ) then
            H2_val = self%H2(1)
        else if ( x_eff >= self%EFTc%x_final ) then
            H2_val = self%H2(self%EFTc%num_points)
        else
            ind = int( ( x_eff - self%EFTc%x_initial )/self%EFTc%grid_width ) + 1
            if ( ind < 1 ) ind = 1
            if ( ind > self%EFTc%num_points - 1 ) ind = self%EFTc%num_points - 1
            x1 = self%EFTc%x(ind)
            x2 = self%EFTc%x(ind+1)
            mu = ( x_eff - x1 )/( x2 - x1 )
            H2_val = self%H2(ind)*( 1._dl - mu ) + self%H2(ind+1)*mu
        end if

        if ( IsNaN(H2_val) .or. H2_val < 0._dl ) then
            call MpiStop('Horndeski: H2 is NaN or negative in compute_adotoa')
        end if

        if ( H2_val > 0._dl ) then
            eft_cache%adotoa = sqrt( H2_val )
        else
            eft_cache%adotoa = 0._dl
        end if

    end subroutine EFTCAMBHorndeskiComputeAdotoa

    !-----------------------------------------------------------------
    !> Override Hdot/Hdotdot with user-specified expressions.
    !-----------------------------------------------------------------
    subroutine EFTCAMBHorndeskiComputeHubbleDer( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horndeski)                     :: self          !< the base class
        real(dl)              , intent(in)           :: a             !< input scale factor
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< EFT parameter cache
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< timestep cache

        real(dl) :: H, Omega, OmegaP, OmegaPP, OmegaPPP
        real(dl) :: cdot, ca2_over_m0sq, grho_matter, gpres_matter, gpres_matter_dot
        real(dl) :: LLambda, LLambda_dot, denom1
        real(dl) :: XXX

        H             = eft_cache%adotoa
        Omega         = eft_cache%EFTOmegaV
        OmegaP        = eft_cache%EFTOmegaP
        OmegaPP       = eft_cache%EFTOmegaPP
        OmegaPPP      = eft_cache%EFTOmegaPPP
        cdot          = eft_cache%EFTcdot
        ca2_over_m0sq = eft_cache%EFTc
        grho_matter   = eft_cache%grhom_t
        gpres_matter  = eft_cache%gpresm_t
        gpres_matter_dot = eft_cache%gpresdotm_t
        LLambda       = eft_cache%EFTLambda
        LLambda_dot   = eft_cache%EFTLambdadot

        ! denom1 = 3._dl*H*(1._dl + Omega + a*OmegaP)
        ! if ( denom1 /= 0._dl ) then
        !     eft_cache%Hdot = ( cdot + 2._dl*H*ca2_over_m0sq - 3._dl*a*H**3*OmegaP - 1.5_dl*a**2*H**3*OmegaPP &
        !         & - 0.5_dl*H*grho_matter - 1.5_dl*H*gpres_matter - 0.5_dl*LLambda_dot - H*LLambda )/denom1
        ! else
        !     eft_cache%Hdot = 0._dl
        ! end if

        denom1 = 2._dl*(1._dl + Omega + 0.5_dl*a*OmegaP)
        if ( denom1 /= 0._dl ) then
            eft_cache%Hdot = ( - gpres_matter - LLambda - (1._dl + Omega + 2._dl*a*OmegaP + a**2*OmegaPP)*H**2 )/denom1
            eft_cache%Hdotdot = ( -gpres_matter_dot - 2._dl*H*gpres_matter - LLambda_dot - 2._dl*H*LLambda - (3._dl*a*OmegaP + 4._dl*a**2*OmegaPP + a**3*OmegaPPP)*H**3 - 2._dl*(1 + Omega + 2._dl*a*OmegaP + a**2*OmegaPP)*H*eft_cache%Hdot - 2._dl*(1.5_dl*a*H*OmegaP + 0.5_dl*a**2*H*OmegaPP)*eft_cache%Hdot )/denom1
        else
            eft_cache%Hdot = 0._dl
            eft_cache%Hdotdot = 0._dl
        end if
        

    end subroutine EFTCAMBHorndeskiComputeHubbleDer


    !-----------------------------------------------------------------
    !> 二阶 EFT 函数（Gammas）：目前 Horndeski 路径先全部设为 0（stub），
    !!  和 5e quintessence 的处理一致。将来若要做 full Horndeski mapping，
    !!  再在此处用 Ω, c, Λ 及其导数重建 γ_i(a)。
    !-----------------------------------------------------------------
    subroutine EFTCAMBHorndeskiSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_Horndeski)                     :: self
        real(dl)              , intent(in)           :: a
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache

        ! For now, we do not implement the full Horndeski mapping of Gamma_i.
        ! Set all second-order EFT functions to zero, as in the 5e quintessence model.

        eft_cache%EFTGamma1V  = 0._dl
        eft_cache%EFTGamma1P  = 0._dl

        eft_cache%EFTGamma2V  = 0._dl
        eft_cache%EFTGamma2P  = 0._dl

        eft_cache%EFTGamma3V  = 0._dl
        eft_cache%EFTGamma3P  = 0._dl

        eft_cache%EFTGamma4V  = 0._dl
        eft_cache%EFTGamma4P  = 0._dl
        eft_cache%EFTGamma4PP = 0._dl

        eft_cache%EFTGamma5V  = 0._dl
        eft_cache%EFTGamma5P  = 0._dl

        eft_cache%EFTGamma6V  = 0._dl
        eft_cache%EFTGamma6P  = 0._dl

    end subroutine EFTCAMBHorndeskiSecondOrderEFTFunctions


end module EFTCAMB_FM_horndeski
