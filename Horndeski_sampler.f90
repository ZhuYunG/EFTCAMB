program Horndeski_sampler
    use precision
    use IniFile
    use constants
    use EFTCAMB_cache                 ! 从 EFTCosmoMC/EFTCAMB 路径取用
    use EFTCAMB_FM_horndeski          ! 从 EFTCosmoMC/EFTCAMB 路径取用
    implicit none

    type(TIniFile)              :: Ini
    type(EFTCAMB_Horndeski)     :: Hor
    type(EFTCAMB_parameter_cache) :: cache

    integer :: n_samp, np, i, logunit, sample_id, ini_unit
    real(dl), allocatable :: params(:), pmin(:), pmax(:)
    character(len=64), allocatable :: par_names(:)
    logical :: ok, ini_error
    character(len=256) :: config_file

    config_file = 'params_EFT.ini'
    call get_command_argument(1, config_file)
    if (len_trim(config_file) == 0) config_file = 'params_EFT.ini'

    ini_unit  = 99          ! 随便挑一个没用过的 unit 号，避免 5,6 之类
    call Ini_Open_File(Ini, trim(config_file), ini_unit, ini_error)
    if (ini_error) then
        write(*,*) 'Error opening ini file: ', trim(config_file)
        stop
    end if


    n_samp = Ini_Read_Int_File(Ini, 'sampler_number', 100)

    call Hor%read_model_selection(Ini)
    call Hor%allocate_model_selection(Ini)
    call Hor%compute_param_number()
    np = Hor%parameter_number
    if (np <= 0) then
        write(*,*) 'Horndeski sampler: parameter_number <= 0, nothing to do.'
        stop
    end if

    allocate(params(np), pmin(np), pmax(np), par_names(np))

    do i = 1, np
        call Hor%parameter_names(i, par_names(i))
        pmin(i) = Ini_Read_Double_File(Ini, trim(par_names(i))//'_min')
        pmax(i) = Ini_Read_Double_File(Ini, trim(par_names(i))//'_max')
    end do

    call cache%initialize()
    call fill_cosmo_cache_from_ini(Ini, cache)

    open(newunit=logunit, file='Horndeski_sampler_log.dat', status='replace', action='write')
    write(logunit,'(a)', advance='no') '# sample_id success'
    do i = 1, np
        write(logunit,'(a)', advance='no') ' '//trim(par_names(i))
    end do
    write(logunit,*)

    do sample_id = 1, n_samp

        call draw_uniform(np, pmin, pmax, params)
        call Hor%init_model_parameters(params)
        call Hor%initialize_background(cache, 0, ok)

        if (.not. ok) then
            call log_sample(logunit, sample_id, .false., params)
            cycle
        end if

        call write_sample_file(sample_id, Hor)
        call log_sample(logunit, sample_id, .true., params)
    end do

    close(logunit)
    write(*,'(a,i0,a)') 'Horndeski sampler completed, attempted samples: ', n_samp, '.'

contains

    subroutine draw_uniform(np, pmin, pmax, params)
        integer, intent(in) :: np
        real(dl), intent(in) :: pmin(:), pmax(:)
        real(dl), intent(out) :: params(:)
        real(dl) :: u
        integer :: j

        do j = 1, np
            call random_number(u)
            params(j) = pmin(j) + u*(pmax(j) - pmin(j))
        end do
    end subroutine draw_uniform

    subroutine fill_cosmo_cache_from_ini(Ini, cache)
        type(TIniFile), intent(in) :: Ini
        type(EFTCAMB_parameter_cache), intent(inout) :: cache

        real(dl) :: hubble, h, ombh2, omch2, omnuh2, omk
        real(dl) :: temp_cmb, neff_massless
        real(dl) :: omega_gamma_h2, omega_gamma, omega_nu_massless, omega_nu_massless_h2
        real(dl) :: omega_b, omega_c, omega_nu_massive, omega_k, omega_v
        real(dl) :: omega_r_total

        hubble         = Ini_Read_Double_File(Ini, 'hubble', 70._dl)
        h              = hubble/100._dl
        ombh2          = Ini_Read_Double_File(Ini, 'ombh2', 0.0226_dl)
        omch2          = Ini_Read_Double_File(Ini, 'omch2', 0.112_dl)
        omnuh2         = Ini_Read_Double_File(Ini, 'omnuh2', 0._dl)
        omk            = Ini_Read_Double_File(Ini, 'omk', 0._dl)
        temp_cmb       = Ini_Read_Double_File(Ini, 'temp_cmb', 2.7255_dl)
        neff_massless  = Ini_Read_Double_File(Ini, 'massless_neutrinos', 3.046_dl)

        omega_b           = ombh2/(h*h)
        omega_c           = omch2/(h*h)
        omega_nu_massive  = omnuh2/(h*h)
        omega_k           = omk

        omega_gamma_h2      = 2.469e-5_dl*(temp_cmb/2.7255_dl)**4
        omega_gamma         = omega_gamma_h2/(h*h)
        omega_nu_massless_h2= omega_gamma_h2 * 0.22710731766_dl * neff_massless
        omega_nu_massless   = omega_nu_massless_h2/(h*h)

        omega_r_total = omega_gamma + omega_nu_massless
        omega_v       = max(0._dl, 1._dl - (omega_b + omega_c + omega_nu_massive + omega_r_total + omega_k))

        cache%omegab = omega_b
        cache%omegac = omega_c
        cache%omegan = omega_nu_massive
        cache%omegag = omega_gamma
        cache%omegar = omega_nu_massless
        cache%omegak = omega_k
        cache%omegav = omega_v

        cache%h0     = h
        cache%h0_Mpc = hubble*1000._dl/c

        cache%grhob       = 3._dl*cache%h0_Mpc**2*cache%omegab
        cache%grhoc       = 3._dl*cache%h0_Mpc**2*cache%omegac
        cache%grhornomass = 3._dl*cache%h0_Mpc**2*cache%omegar
        cache%grhog       = 3._dl*cache%h0_Mpc**2*cache%omegag
        cache%grhok       = 3._dl*cache%h0_Mpc**2*cache%omegak
        cache%grhov       = 3._dl*cache%h0_Mpc**2*cache%omegav

        cache%Num_Nu_Massive      = 0
        cache%Nu_mass_eigenstates = 0
        if (allocated(cache%grhormass)) deallocate(cache%grhormass)
        if (allocated(cache%nu_masses)) deallocate(cache%nu_masses)
    end subroutine fill_cosmo_cache_from_ini

    subroutine write_sample_file(idx, model)
        integer, intent(in) :: idx
        type(EFTCAMB_Horndeski), intent(in) :: model

        integer :: unit_id, npts, k
        real(dl) :: a, H_phys, H_phys2, wde, omegade, rhoDE_real
        character(len=64) :: fname

        fname = 'Horndeski_sample_'//trim(itoa(idx))//'.dat'
        open(newunit=unit_id, file=trim(fname), status='replace', action='write')
        write(unit_id,'(a)') '# a  H(a)  H2(a)  w_DE  Omega_DE  rhoDE_real'

        npts = model%EFTc%num_points
        do k = 1, npts
            a          = exp(model%EFTc%x(k))
            H_phys     = model%Hphys(k)
            H_phys2    = H_phys*H_phys
            wde        = model%wDE(k)
            omegade    = model%OmegaDE(k)
            rhoDE_real = model%rhoDE(k)
            write(unit_id,'(6ES18.8E3)') a, H_phys, H_phys2, wde, omegade, rhoDE_real
        end do

        close(unit_id)
    end subroutine write_sample_file

    subroutine log_sample(logunit, idx, success, params)
        integer, intent(in) :: logunit, idx
        logical, intent(in) :: success
        real(dl), intent(in) :: params(:)

        write(logunit,'(i10,1x,l1,*(1x,ES16.6))') idx, success, params
    end subroutine log_sample

    function itoa(k) result(str)
        integer, intent(in) :: k
        character(len=16) :: str
        write(str,'(I0)') k
    end function itoa

end program Horndeski_sampler
