!-----------------------------------------------------------------------
! scram_box: box-model driver for the CESM-SCRAM aerosol kernel.
!
! It calls the kernel through the same interface and in the same order as
! CESM (carma_model_mod, CARMA_DiagnoseBins):
!   scram_cesm_init_from_discretization  (once)
!   per step:  emissions and gas sources (operator splitting, before SCRAM,
!              as CAM applies its tendencies before SCRAM), then dilution
!              towards a background state, x <- x - (1-exp(-k dt)) (x - x_bg)
!              SOA photolysis loss exp(-0.0004 J_NO2 dt) (as CESM)
!              scram_cesm_set_met, scram_cesm_set_hygro_rh, scram_cesm_set_state
!              n_sub x scram_cesm_step(dt/n_sub)   (scram_target_dt, 600 s in CESM)
!              scram_cesm_get_state
!
! Inputs (written by the Python case compiler, scrambox.compile):
!   run.nml     scalars, process switches, species ids, file names
!   inputs.txt  arrays as keyword blocks: "<key> <n>" followed by n values
! Output:
!   output.txt  keyword blocks per output time (parsed by scrambox.results)
! The species table is read by the kernel from SCRAM_SPECIES_FILE.
!
! Usage: scram_box [run.nml]
!-----------------------------------------------------------------------
program scram_box

  use scram_shared_mod
  use scram_cesm_driver_mod, only: scram_cesm_init_from_discretization, &
       scram_cesm_set_met, scram_cesm_set_hygro_rh, scram_cesm_set_state, &
       scram_cesm_get_state, scram_cesm_step
  use time_manager, only: box_set_nstep

  implicit none

  ! ---- namelist variables
  integer :: nsize = 0, nfrac = 1, naer = 0, ngas = 0, ngroup = 1
  double precision :: duration_s = 0.d0, dt_s = 600.d0, substep_s = 600.d0
  double precision :: output_every_s = 600.d0
  integer :: with_coagulation = 1, with_condensation = 1, with_nucleation = 0
  integer :: nucleation_model = 0, dynamic_solver = 1, thermodynamic_model = 0
  integer :: use_isorropia = 1, sulfate_computation = 1, redistribution_method = 3
  integer :: rdb_core_conserv = 8, wet_diam_estimation = 0, density_mode = 0
  double precision :: fixed_density_kgm3 = 1400.d0
  double precision :: bc_core_floor_um = 0.05d0, salt_core_floor_um = 0.05d0
  double precision :: dust_core_floor_um = 0.90d0
  integer :: salt_rdb_bidirectional = 0, salt_fix_size_rdb = 0, dust_fix_size_rdb = 0
  double precision :: nucl_rate_cap_cm3_s = 1.0d3
  double precision :: nucl_h2so4_gas_fraction_max = 0.25d0
  double precision :: nucl_nh3_gas_fraction_max = 0.25d0
  integer :: sid_pmd = -1, sid_pbc = -1, sid_pna = -1, sid_pso4 = -1, sid_pnh4 = -1
  integer :: sid_pno3 = -1, sid_phcl = -1, sid_poa = -1, sid_ph2o = -1
  integer :: sid_soa(5) = -1
  integer :: gid_sulf = -1
  integer :: soa_photolysis = 1
  character(len=512) :: inputs_file = 'inputs.txt', output_file = 'output.txt'
  namelist /box_dims/ nsize, nfrac, naer, ngas, ngroup
  namelist /box_time/ duration_s, dt_s, substep_s, output_every_s
  namelist /box_proc/ with_coagulation, with_condensation, with_nucleation, &
       nucleation_model, dynamic_solver, thermodynamic_model, use_isorropia, &
       sulfate_computation, redistribution_method, rdb_core_conserv, &
       wet_diam_estimation, density_mode, fixed_density_kgm3, bc_core_floor_um, &
       salt_core_floor_um, dust_core_floor_um, salt_rdb_bidirectional, &
       salt_fix_size_rdb, dust_fix_size_rdb, nucl_rate_cap_cm3_s, &
       nucl_h2so4_gas_fraction_max, nucl_nh3_gas_fraction_max, soa_photolysis
  namelist /box_ids/ sid_pmd, sid_pbc, sid_pna, sid_pso4, sid_pnh4, sid_pno3, &
       sid_phcl, sid_poa, sid_ph2o, sid_soa, gid_sulf
  namelist /box_files/ inputs_file, output_file

  ! ---- arrays
  double precision, allocatable :: bin_bound_m(:), comp_bounds(:), density(:)
  integer, allocatable :: group_rel(:)
  double precision, allocatable :: num(:), mass(:,:), gas(:)
  double precision, allocatable :: emis_num(:), emis_mass(:,:), gas_src(:)
  double precision, allocatable :: bg_num(:), bg_mass(:,:), bg_gas(:)
  ! (NMETCOL, nmet): t, T, RH, P, J_NO2, particle emission factor, gas source factor,
  ! dilution rate (1/s)
  integer, parameter :: NMETCOL = 8
  double precision, allocatable :: met(:,:)
  integer :: nbin, nmet, rc, u, ios, istep, nstep, isub, nsub, k, iout
  double precision :: t, tk, rh, p, jno2, efac, gfac, kdil, fdil, dt_sub, next_out, survival
  double precision :: mv(NMETCOL)
  double precision :: t0_wall, t1_wall
  character(len=512) :: nml_file

  ! ---- read the namelist
  nml_file = 'run.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, nml_file)
  open(newunit=u, file=trim(nml_file), status='old', action='read', iostat=ios)
  if (ios /= 0) call fail('cannot open '//trim(nml_file))
  read(u, nml=box_dims, iostat=ios);  if (ios /= 0) call fail('bad &box_dims')
  rewind(u); read(u, nml=box_time, iostat=ios);  if (ios /= 0) call fail('bad &box_time')
  rewind(u); read(u, nml=box_proc, iostat=ios);  if (ios /= 0) call fail('bad &box_proc')
  rewind(u); read(u, nml=box_ids, iostat=ios);   if (ios /= 0) call fail('bad &box_ids')
  rewind(u); read(u, nml=box_files, iostat=ios); if (ios /= 0) call fail('bad &box_files')
  close(u)
  if (nsize <= 0 .or. nfrac <= 0 .or. naer <= 0 .or. ngas <= 0) call fail('bad dimensions')
  if (dt_s <= 0.d0 .or. duration_s < 0.d0) call fail('bad time settings')
  nbin = nsize * nfrac

  ! ---- process switches: the same module variables CAM's cfg parser sets
  with_coag = with_coagulation
  with_cond = with_condensation
  with_nucl = with_nucleation
  scram_nucl_model = nucleation_model
  scram_coag_repart_mode = 1                  ! NEAREST analytical repartition
  scram_dynamic_solver = dynamic_solver
  scram_thermodynamic_model = thermodynamic_model
  scram_use_isorropia = use_isorropia
  scram_sulfate_computation = sulfate_computation
  scram_redistribution_method = redistribution_method
  scram_rdb_core_conserv = rdb_core_conserv
  scram_wet_diam_estimation = wet_diam_estimation
  scram_tagrho = density_mode
  scram_fixed_density_kgm3 = fixed_density_kgm3
  scram_bc_core_floor_diam_um = bc_core_floor_um
  scram_salt_core_floor_diam_um = salt_core_floor_um
  scram_dust_core_floor_diam_um = dust_core_floor_um
  scram_salt_rdb_bidirectional = salt_rdb_bidirectional
  scram_salt_fix_size_rdb = salt_fix_size_rdb
  scram_dust_fix_size_rdb = dust_fix_size_rdb
  scram_nucl_rate_cap_cm3_s = nucl_rate_cap_cm3_s
  scram_nucl_h2so4_gas_fraction_max = nucl_h2so4_gas_fraction_max
  scram_nucl_nh3_gas_fraction_max = nucl_nh3_gas_fraction_max
  scram_ptop_pa = 0.d0
  scram_sid_pmd = sid_pmd;   scram_sid_pbc = sid_pbc;   scram_sid_pna = sid_pna
  scram_sid_pso4 = sid_pso4; scram_sid_pnh4 = sid_pnh4; scram_sid_pno3 = sid_pno3
  scram_sid_phcl = sid_phcl; scram_sid_poa = sid_poa;   scram_sid_ph2o = sid_ph2o
  scram_sid_soa1 = sid_soa(1); scram_sid_soa2 = sid_soa(2); scram_sid_soa3 = sid_soa(3)
  scram_sid_soa4 = sid_soa(4); scram_sid_soa5 = sid_soa(5)
  scram_gid_sulf = gid_sulf

  ! ---- read the arrays
  open(newunit=u, file=trim(inputs_file), status='old', action='read', iostat=ios)
  if (ios /= 0) call fail('cannot open '//trim(inputs_file))
  allocate(bin_bound_m(nsize+1), comp_bounds(nfrac*ngroup*2), density(naer), group_rel(naer))
  allocate(num(nbin), mass(nbin,naer), gas(ngas), emis_num(nbin), emis_mass(nbin,naer), gas_src(ngas))
  allocate(bg_num(nbin), bg_mass(nbin,naer), bg_gas(ngas))
  call read_block(u, 'bin_bound_m', bin_bound_m)
  call read_block(u, 'comp_bounds', comp_bounds)
  call read_block(u, 'density_ugum3', density)
  call read_block_int(u, 'group_rel', group_rel)
  call read_block(u, 'number_m3', num)
  call read_block2(u, 'mass_ugm3', mass)
  call read_block(u, 'gas_ugm3', gas)
  call read_block(u, 'emis_number_m3s', emis_num)
  call read_block2(u, 'emis_mass_ugm3s', emis_mass)
  call read_block(u, 'gas_source_ugm3s', gas_src)
  call read_block(u, 'bg_number_m3', bg_num)
  call read_block2(u, 'bg_mass_ugm3', bg_mass)
  call read_block(u, 'bg_gas_ugm3', bg_gas)
  call read_met(u)
  close(u)

  ! ---- initialise the kernel exactly as CESM does
  call scram_cesm_init_from_discretization(nsize, nfrac, naer, bin_bound_m, comp_bounds, &
       density, group_rel, rc, n_gas=ngas)
  if (rc /= 0) call fail('scram_cesm_init_from_discretization failed')

  open(newunit=iout, file=trim(output_file), status='replace', action='write')
  write(iout,'(A)') '# scram_box output v1'
  write(iout,'(A,5I8)') 'dims', nsize, nfrac, naer, ngas, ngroup
  call cpu_time(t0_wall)

  nstep = nint(duration_s / dt_s)
  nsub = max(1, nint(dt_s / max(substep_s, 1.d-9)))
  dt_sub = dt_s / dble(nsub)
  t = 0.d0
  call met_at(t, mv)
  call write_state(t, mv(2), mv(3), mv(4))
  next_out = output_every_s

  do istep = 1, nstep
    call box_set_nstep(istep)
    call met_at(t, mv)
    tk = mv(2); rh = mv(3); p = mv(4); jno2 = mv(5)
    efac = mv(6); gfac = mv(7); kdil = max(mv(8), 0.d0)

    ! emissions and gas sources over the step (before SCRAM, as in CESM)
    num  = num  + emis_num  * efac * dt_s
    mass = mass + emis_mass * efac * dt_s
    gas  = gas  + gas_src   * gfac * dt_s

    ! dilution (mixing with the background air mass)
    if (kdil > 0.d0) then
      fdil = 1.d0 - exp(-kdil * dt_s)
      num  = num  - fdil * (num  - bg_num)
      mass = mass - fdil * (mass - bg_mass)
      gas  = gas  - fdil * (gas  - bg_gas)
    end if

    ! SOA photolysis, the same first-order loss CESM applies before SCRAM
    if (soa_photolysis /= 0 .and. jno2 > 0.d0) then
      survival = exp(-0.0004d0 * jno2 * dt_s)
      do k = 1, 5
        if (sid_soa(k) >= 1 .and. sid_soa(k) <= naer) mass(:,sid_soa(k)) = mass(:,sid_soa(k)) * survival
      end do
    end if

    call scram_cesm_set_met(tk, rh, p)
    call scram_cesm_set_hygro_rh(min(max(rh, 0.d0), 0.98d0))
    call scram_cesm_set_state(num, mass, gas, rc)
    if (rc /= 0) call fail('scram_cesm_set_state failed')
    do isub = 1, nsub
      call scram_cesm_step(dt_sub, rc)
      if (rc /= 0) call fail('scram_cesm_step failed')
    end do
    call scram_cesm_get_state(num, mass, gas, rc)
    if (rc /= 0) call fail('scram_cesm_get_state failed')

    t = dble(istep) * dt_s
    if (t >= next_out - 1.d-6 .or. istep == nstep) then
      call write_state(t, tk, rh, p)
      next_out = next_out + output_every_s
    end if
  end do

  call cpu_time(t1_wall)
  write(iout,'(A,ES17.8E3)') 'cpu_s', t1_wall - t0_wall
  write(iout,'(A)') 'end'
  close(iout)
  write(*,'(A,I0,A,I0,A,F10.2,A)') 'scram_box: ', nstep, ' steps x ', nsub, ' substeps, cpu ', &
       t1_wall - t0_wall, ' s'

contains

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    write(0,'(A)') 'scram_box ERROR: '//trim(msg)
    write(*,'(A)') 'scram_box ERROR: '//trim(msg)
    error stop 1
  end subroutine fail

  subroutine expect_key(uu, key, n_expected)
    integer, intent(in) :: uu, n_expected
    character(len=*), intent(in) :: key
    character(len=64) :: k
    integer :: n, ios2
    read(uu, *, iostat=ios2) k, n
    if (ios2 /= 0) call fail('inputs: cannot read header of block '//key)
    if (trim(k) /= key) call fail('inputs: expected block '//key//' found '//trim(k))
    if (n /= n_expected) call fail('inputs: wrong size for block '//key)
  end subroutine expect_key

  subroutine read_block(uu, key, a)
    integer, intent(in) :: uu
    character(len=*), intent(in) :: key
    double precision, intent(out) :: a(:)
    integer :: ios2
    call expect_key(uu, key, size(a))
    read(uu, *, iostat=ios2) a
    if (ios2 /= 0) call fail('inputs: cannot read values of block '//key)
  end subroutine read_block

  subroutine read_block_int(uu, key, a)
    integer, intent(in) :: uu
    character(len=*), intent(in) :: key
    integer, intent(out) :: a(:)
    integer :: ios2
    call expect_key(uu, key, size(a))
    read(uu, *, iostat=ios2) a
    if (ios2 /= 0) call fail('inputs: cannot read values of block '//key)
  end subroutine read_block_int

  subroutine read_block2(uu, key, a)
    ! column-major: all bins of species 1, then species 2, ...
    integer, intent(in) :: uu
    character(len=*), intent(in) :: key
    double precision, intent(out) :: a(:,:)
    integer :: ios2
    call expect_key(uu, key, size(a))
    read(uu, *, iostat=ios2) a
    if (ios2 /= 0) call fail('inputs: cannot read values of block '//key)
  end subroutine read_block2

  subroutine read_met(uu)
    integer, intent(in) :: uu
    character(len=64) :: k
    integer :: ncol, ios2
    read(uu, *, iostat=ios2) k, ncol, nmet
    if (ios2 /= 0 .or. trim(k) /= 'met' .or. ncol /= NMETCOL .or. nmet < 1) &
         call fail('inputs: bad met block (expected "met 8 n")')
    allocate(met(NMETCOL, nmet))
    read(uu, *, iostat=ios2) met
    if (ios2 /= 0) call fail('inputs: cannot read met values')
  end subroutine read_met

  subroutine met_at(tt, v)
    ! piecewise-linear in time, constant outside the given range
    double precision, intent(in) :: tt
    double precision, intent(out) :: v(NMETCOL)
    integer :: i
    double precision :: w
    if (nmet == 1 .or. tt <= met(1,1)) then
      v = met(:,1)
      return
    else if (tt >= met(1,nmet)) then
      v = met(:,nmet)
      return
    end if
    i = 1
    do while (met(1,i+1) < tt)
      i = i + 1
    end do
    w = (tt - met(1,i)) / max(met(1,i+1) - met(1,i), 1.d-30)
    v = met(:,i) + w * (met(:,i+1) - met(:,i))
  end subroutine met_at

  subroutine write_state(tt, tk_o, rh_o, p_o)
    double precision, intent(in) :: tt, tk_o, rh_o, p_o
    write(iout,'(A,4ES17.8E3)') 'time', tt, tk_o, rh_o, p_o
    write(iout,'(A,I0)') 'number_m3 ', nbin
    write(iout,'(5ES17.8E3)') num
    write(iout,'(A,I0)') 'mass_ugm3 ', nbin*naer
    write(iout,'(5ES17.8E3)') mass
    write(iout,'(A,I0)') 'gas_ugm3 ', ngas
    write(iout,'(5ES17.8E3)') gas
  end subroutine write_state

end program scram_box
