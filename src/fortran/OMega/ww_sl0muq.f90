! O'Mega scattering amplitudes for the processes
!
!   allowed: (after merging flavor and color states)
!
!     e- e+ -> u/1/ dbar//1 mu- numubar
!     e- e+ -> c/1/ sbar//1 mu- numubar
!     e- e+ -> ubar//1 d/1/ mu+ numu
!     e- e+ -> cbar//1 s/1/ mu+ numu
!
! in Gauged Colorization Functor ( minimal electroweak standard model in unitarity gauge )
!
module ww_sl0muq
  use omega_kinds !NODEP!
  use omega95 !NODEP!
  use omega_parameters !NODEP!
  ! use parameters
  ! use omega_parameters_whizard
  implicit none
  private
  public :: scatter_nonzero,  scatter_diagonal_nonzero, scatter_diagonal_polarized, &
    scatter_diagonal_colored_nz, scatter_colored_nonzero
  public :: allocate_zero
  interface allocate_zero
     module procedure allocate_zero_1, allocate_zero_2
  end interface
  public :: allocate_zero_1, allocate_zero_2
  public :: number_particles, number_particles_in, number_particles_out
  public :: number_spin_states, number_spin_states_in, &
    number_spin_states_out,  spin_states, spin_states_in, spin_states_out
  public :: number_flavor_states, number_flavor_states_in, &
    number_flavor_states_out, flavor_states,  flavor_states_in, &
    flavor_states_out
  public :: number_flavor_zeros, number_flavor_zeros_in, &
    number_flavor_zeros_out, flavor_zeros, flavor_zeros_in, flavor_zeros_out
  public :: number_color_flows, color_flows, anticolor_flows
!  public :: create, reset, destroy
  ! DON'T EVEN THINK of removing the following!
  ! If the compiler complains about undeclared
  ! or undefined variables, you are compiling
  ! against an incompatible omega95 module!
  integer, dimension(6), parameter, private :: require = &
    (/ omega_spinors_2003_03_A, omega_spinor_cpls_2003_03_A, &
    omega_vectors_2003_03_A, omega_polarizations_2003_03_A, &
    omega_couplings_2003_03_A, omega_utils_2003_03_A /)

  integer, parameter, private :: n_prt     = 6
  integer, parameter, private :: n_in      = 2
  integer, parameter, private :: n_out     = 4
  integer, parameter, private :: n_cflow   = 2
  integer, parameter, private :: n_flv     = 4
  integer, parameter, private :: n_flv_in  = 1
  integer, parameter, private :: n_flv_out = 4
  integer, parameter, private :: n_hel     = 64
  integer, parameter, private :: n_hel_in  = 4
  integer, parameter, private :: n_hel_out = 16

  logical, parameter, private :: F = .false.
  logical, parameter, private :: T = .true.

  integer, dimension(n_prt), parameter, private ::  s0001 = &
    (/ -1, -1, -1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0002 = &
    (/ -1, -1, -1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0003 = &
    (/ -1, -1, -1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0004 = &
    (/ -1, -1, -1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0005 = &
    (/ -1, -1, -1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0006 = &
    (/ -1, -1, -1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0007 = &
    (/ -1, -1, -1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0008 = &
    (/ -1, -1, -1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0009 = &
    (/ -1, -1,  1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0010 = &
    (/ -1, -1,  1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0011 = &
    (/ -1, -1,  1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0012 = &
    (/ -1, -1,  1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0013 = &
    (/ -1, -1,  1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0014 = &
    (/ -1, -1,  1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0015 = &
    (/ -1, -1,  1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0016 = &
    (/ -1, -1,  1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0017 = &
    (/ -1,  1, -1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0018 = &
    (/ -1,  1, -1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0019 = &
    (/ -1,  1, -1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0020 = &
    (/ -1,  1, -1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0021 = &
    (/ -1,  1, -1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0022 = &
    (/ -1,  1, -1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0023 = &
    (/ -1,  1, -1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0024 = &
    (/ -1,  1, -1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0025 = &
    (/ -1,  1,  1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0026 = &
    (/ -1,  1,  1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0027 = &
    (/ -1,  1,  1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0028 = &
    (/ -1,  1,  1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0029 = &
    (/ -1,  1,  1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0030 = &
    (/ -1,  1,  1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0031 = &
    (/ -1,  1,  1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0032 = &
    (/ -1,  1,  1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0033 = &
    (/  1, -1, -1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0034 = &
    (/  1, -1, -1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0035 = &
    (/  1, -1, -1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0036 = &
    (/  1, -1, -1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0037 = &
    (/  1, -1, -1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0038 = &
    (/  1, -1, -1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0039 = &
    (/  1, -1, -1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0040 = &
    (/  1, -1, -1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0041 = &
    (/  1, -1,  1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0042 = &
    (/  1, -1,  1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0043 = &
    (/  1, -1,  1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0044 = &
    (/  1, -1,  1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0045 = &
    (/  1, -1,  1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0046 = &
    (/  1, -1,  1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0047 = &
    (/  1, -1,  1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0048 = &
    (/  1, -1,  1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0049 = &
    (/  1,  1, -1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0050 = &
    (/  1,  1, -1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0051 = &
    (/  1,  1, -1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0052 = &
    (/  1,  1, -1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0053 = &
    (/  1,  1, -1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0054 = &
    (/  1,  1, -1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0055 = &
    (/  1,  1, -1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0056 = &
    (/  1,  1, -1,  1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0057 = &
    (/  1,  1,  1, -1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0058 = &
    (/  1,  1,  1, -1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0059 = &
    (/  1,  1,  1, -1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0060 = &
    (/  1,  1,  1, -1,  1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0061 = &
    (/  1,  1,  1,  1, -1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0062 = &
    (/  1,  1,  1,  1, -1,  1 /)
  integer, dimension(n_prt), parameter, private ::  s0063 = &
    (/  1,  1,  1,  1,  1, -1 /)
  integer, dimension(n_prt), parameter, private ::  s0064 = &
    (/  1,  1,  1,  1,  1,  1 /)
  integer, dimension(n_prt,n_hel), parameter, private :: table_spin_states = &
    reshape ( (/ s0001, s0002, s0003, s0004, s0005, s0006, s0007, s0008, &
    s0009,  s0010, s0011, s0012, s0013, s0014, s0015, s0016, s0017, s0018, &
    s0019,  s0020, s0021, s0022, s0023, s0024, s0025, s0026, s0027, s0028, &
    s0029,  s0030, s0031, s0032, s0033, s0034, s0035, s0036, s0037, s0038, &
    s0039,  s0040, s0041, s0042, s0043, s0044, s0045, s0046, s0047, s0048, &
    s0049,  s0050, s0051, s0052, s0053, s0054, s0055, s0056, s0057, s0058, &
    s0059,  s0060, s0061, s0062, s0063, s0064 /), (/ n_prt, n_hel /) )
  integer, dimension(n_in), parameter, private :: si0001 = (/ -1, -1 /)
  integer, dimension(n_in), parameter, private :: si0002 = (/ -1,  1 /)
  integer, dimension(n_in), parameter, private :: si0003 = (/  1, -1 /)
  integer, dimension(n_in), parameter, private :: si0004 = (/  1,  1 /)
  integer, dimension(n_in,n_hel_in), parameter, private :: &
    table_spin_states_in = &
    reshape ( (/ si0001, si0002, si0003, si0004 /), (/ n_in, n_hel_in /) )
  integer, dimension(n_out), parameter, private :: so0001 = &
    (/ -1, -1, -1, -1 /)
  integer, dimension(n_out), parameter, private :: so0002 = &
    (/ -1, -1, -1,  1 /)
  integer, dimension(n_out), parameter, private :: so0003 = &
    (/ -1, -1,  1, -1 /)
  integer, dimension(n_out), parameter, private :: so0004 = &
    (/ -1, -1,  1,  1 /)
  integer, dimension(n_out), parameter, private :: so0005 = &
    (/ -1,  1, -1, -1 /)
  integer, dimension(n_out), parameter, private :: so0006 = &
    (/ -1,  1, -1,  1 /)
  integer, dimension(n_out), parameter, private :: so0007 = &
    (/ -1,  1,  1, -1 /)
  integer, dimension(n_out), parameter, private :: so0008 = &
    (/ -1,  1,  1,  1 /)
  integer, dimension(n_out), parameter, private :: so0009 = &
    (/  1, -1, -1, -1 /)
  integer, dimension(n_out), parameter, private :: so0010 = &
    (/  1, -1, -1,  1 /)
  integer, dimension(n_out), parameter, private :: so0011 = &
    (/  1, -1,  1, -1 /)
  integer, dimension(n_out), parameter, private :: so0012 = &
    (/  1, -1,  1,  1 /)
  integer, dimension(n_out), parameter, private :: so0013 = &
    (/  1,  1, -1, -1 /)
  integer, dimension(n_out), parameter, private :: so0014 = &
    (/  1,  1, -1,  1 /)
  integer, dimension(n_out), parameter, private :: so0015 = &
    (/  1,  1,  1, -1 /)
  integer, dimension(n_out), parameter, private :: so0016 = &
    (/  1,  1,  1,  1 /)
  integer, dimension(n_out,n_hel_out), parameter, private :: &
    table_spin_states_out = &
    reshape ( (/ so0001, so0002, so0003, so0004, so0005, so0006, so0007, &
    so0008, so0009, so0010, so0011, so0012, so0013, so0014, so0015, so0016 &
    /), (/ n_out, n_hel_out/) )

  integer, dimension(n_prt), parameter, private ::  f0001 = &
    (/  11, -11,   2,  -1,  13, -14 /) ! e- e+ u/1/ dbar//1 mu- numubar
  integer, parameter, private :: f0001m = 1
  integer, dimension(n_prt), parameter, private ::  f0002 = &
    (/  11, -11,   4,  -3,  13, -14 /) ! e- e+ c/1/ sbar//1 mu- numubar
  integer, parameter, private :: f0002m = 1
  integer, dimension(n_prt), parameter, private ::  f0003 = &
    (/  11, -11,  -2,   1, -13,  14 /) ! e- e+ ubar//1 d/1/ mu+ numu
  integer, parameter, private :: f0003m = 1
  integer, dimension(n_prt), parameter, private ::  f0004 = &
    (/  11, -11,  -4,   3, -13,  14 /) ! e- e+ cbar//1 s/1/ mu+ numu
  integer, parameter, private :: f0004m = 1
  integer, dimension(n_prt,n_flv), parameter, private :: table_flavor_states &
    =  reshape ( (/ f0001, f0002, f0003, f0004 /), (/ n_prt, n_flv /) )
  integer, dimension(n_in), parameter, private :: fi0001 = (/  11, -11 /) ! e- e+
  integer, parameter, private :: fi0001m = 1
  integer, dimension(n_in,n_flv_in), parameter, private :: &
    table_flavor_states_in =  reshape ( (/ fi0001 /), (/ n_in, n_flv_in /) )
  integer, dimension(n_out), parameter, private ::  fo0001 = &
    (/   2,  -1,  13, -14 /) ! u/1/ dbar//1 mu- numubar
  integer, parameter, private :: fo0001m = 1
  integer, dimension(n_out), parameter, private ::  fo0002 = &
    (/   4,  -3,  13, -14 /) ! c/1/ sbar//1 mu- numubar
  integer, parameter, private :: fo0002m = 1
  integer, dimension(n_out), parameter, private ::  fo0003 = &
    (/  -4,   3, -13,  14 /) ! cbar//1 s/1/ mu+ numu
  integer, parameter, private :: fo0003m = 1
  integer, dimension(n_out), parameter, private ::  fo0004 = &
    (/  -2,   1, -13,  14 /) ! ubar//1 d/1/ mu+ numu
  integer, parameter, private :: fo0004m = 1
  integer, dimension(n_out,n_flv_out), parameter, private :: &
    table_flavor_states_out = &
    reshape ( (/ fo0001, fo0002, fo0003, fo0004 /), (/ n_out, n_flv_out /) )

  logical, dimension(n_hel,n_flv) :: flv_hel_flag

  integer, dimension(n_prt,0), private :: table_flavor_zeros
  integer, dimension(n_in,0), private  :: table_flavor_zeros_in
  integer, dimension(n_out,0), private :: table_flavor_zeros_out

  integer, dimension(n_prt), parameter, private :: c0001 = &
    (/ 0, 0, 4, 0, 0, 0 /)
  integer, dimension(n_prt), parameter, private :: c0002 = &
    (/ 0, 0, 0, 3, 0, 0 /)
  integer, dimension(n_prt,n_cflow), parameter, private :: table_color_flows &
    = reshape ( (/ c0001, c0002 /), (/ n_prt, n_cflow /) )
  integer, dimension(n_prt), parameter, private :: a0001 = &
    (/ 0, 0, 0, 3, 0, 0 /)
  integer, dimension(n_prt), parameter, private :: a0002 = &
    (/ 0, 0, 4, 0, 0, 0 /)
  integer, dimension(n_prt,n_cflow), parameter, private :: &
    table_anticolor_flows = reshape ( (/ a0001, a0002 /), (/ n_prt, n_cflow/) )
  integer, dimension(n_cflow), parameter, private :: fn0001 = (/ 3, 0 /)
  integer, dimension(n_cflow), parameter, private :: fd0001 = (/ 1, 1 /)
  logical, dimension(n_cflow), parameter, private :: ff0001 = (/ T, F /)
  integer, dimension(n_cflow), parameter, private :: fn0002 = (/ 0, 3 /)
  integer, dimension(n_cflow), parameter, private :: fd0002 = (/ 1, 1 /)
  logical, dimension(n_cflow), parameter, private :: ff0002 = (/ F, T /)
  integer, dimension(n_cflow,n_cflow), parameter, private :: flow_num = &
    reshape ( (/ fn0001, fn0002 /), (/ n_cflow, n_cflow /) )
  integer, dimension(n_cflow,n_cflow), parameter, private :: flow_den = &
    reshape ( (/ fd0001, fd0002 /), (/ n_cflow, n_cflow /) )
  logical, dimension(n_cflow,n_cflow), parameter, private :: flow_flag = &
    reshape ( (/ ff0001, ff0002 /), (/ n_cflow, n_cflow /) )
  real(kind=omega_prec), dimension(n_cflow,n_cflow), save, private :: &
    flow_coeff
  logical, dimension(n_cflow), parameter, private :: flow_is_physical = &
    (/ T, T /)

  logical, dimension(n_cflow), parameter, private :: fc0001 = (/T, F/)
  logical, dimension(n_cflow), parameter, private :: fc0002 = (/T, F/)
  logical, dimension(n_cflow), parameter, private :: fc0003 = (/F, T/)
  logical, dimension(n_cflow), parameter, private :: fc0004 = (/F, T/)
  logical, dimension(n_cflow,n_flv), parameter, private :: flv_col_flag = &
    reshape ( (/ fc0001, fc0002, fc0003, fc0004 /), (/ n_cflow, n_flv/) )

  integer, dimension(n_flv), parameter, private :: table_symmetry = &
    (/   1, 1, 1, 1 /)

contains

  subroutine allocate_zero_1 (zero)
    integer, dimension(:,:), pointer :: zero
    allocate (zero(size(table_spin_states,dim=2), &
      size(table_flavor_states,dim=2)))
    zero = 0
  end subroutine allocate_zero_1

  subroutine allocate_zero_2 (zero)
    integer, dimension(:,:,:,:), pointer :: zero
    allocate(zero(size(table_spin_states_in,dim=2), &
      size(table_flavor_states_in,dim=2),  size(table_spin_states_out,dim=2), &
       size(table_flavor_states_out,dim=2)))
    zero = 0
  end subroutine allocate_zero_2

  subroutine create ()
    flow_coeff = &
      real (flow_num, kind=omega_prec) / real (flow_den, kind=omega_prec)
  end subroutine create
!  subroutine reset (par)
!    type(parameter_set), intent(in) :: par
!    call import_from_whizard (par)
!  end subroutine reset
!  subroutine destroy ()
!  end subroutine destroy

  pure function number_particles () result (n)
    integer :: n
    n = size (table_flavor_states, dim=1)
  end function number_particles
  pure function number_particles_in () result (n)
    integer :: n
    n = size (table_flavor_states_in, dim=1)
  end function number_particles_in
  pure function number_particles_out () result (n)
    integer :: n
    n = size (table_flavor_states_out, dim=1)
  end function number_particles_out

  pure function number_spin_states () result (n)
    integer :: n
    n = size (table_spin_states, dim=2)
  end function number_spin_states
  pure subroutine spin_states (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_spin_states
  end subroutine spin_states
  pure function number_spin_states_in () result (n)
    integer :: n
    n = size (table_spin_states_in, dim=2)
  end function number_spin_states_in
  pure subroutine spin_states_in (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_spin_states_in
  end subroutine spin_states_in
  pure function number_spin_states_out () result (n)
    integer :: n
    n = size (table_spin_states_out, dim=2)
  end function number_spin_states_out
  pure subroutine spin_states_out (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_spin_states_out
  end subroutine spin_states_out

  pure function number_flavor_states () result (n)
    integer :: n
    n = size (table_flavor_states, dim=2)
  end function number_flavor_states
  pure subroutine flavor_states (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_flavor_states
  end subroutine flavor_states
  pure function number_flavor_states_in () result (n)
    integer :: n
    n = size (table_flavor_states_in, dim=2)
  end function number_flavor_states_in
  pure subroutine flavor_states_in (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_flavor_states_in
  end subroutine flavor_states_in
  pure function number_flavor_states_out () result (n)
    integer :: n
    n = size (table_flavor_states_out, dim=2)
  end function number_flavor_states_out
  pure subroutine flavor_states_out (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_flavor_states_out
  end subroutine flavor_states_out

  function number_flavor_zeros () result (n)
    integer :: n
    n = 0
  end function number_flavor_zeros
  subroutine flavor_zeros (a)
    integer, dimension(:,:) :: a
    a = table_flavor_zeros
  end subroutine flavor_zeros
  function number_flavor_zeros_in () result (n)
    integer :: n
    n = 0
  end function number_flavor_zeros_in
  subroutine flavor_zeros_in (a)
    integer, dimension(:,:) :: a
    a = table_flavor_zeros_in
  end subroutine flavor_zeros_in
  function number_flavor_zeros_out () result (n)
    integer :: n
    n = 0
  end function number_flavor_zeros_out
  subroutine flavor_zeros_out (a)
    integer, dimension(:,:) :: a
    a = table_flavor_zeros_out
  end subroutine flavor_zeros_out

  function number_color_flows () result (n)
    integer :: n
    n = size (table_color_flows, dim=2)
  end function number_color_flows
  subroutine color_flows (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_color_flows
  end subroutine color_flows
  subroutine anticolor_flows (a)
    integer, dimension(:,:), intent(inout) :: a
    a = table_anticolor_flows
  end subroutine anticolor_flows

  subroutine calculate_amplitudes (amp, k, zero_ct, n)
    complex(kind=omega_prec), dimension(:,:,:), intent(out) :: amp
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer, parameter :: SAMPLE = 10, REPEAT = 5
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    complex(kind=omega_prec) :: bk1_1_1
    complex(kind=omega_prec) :: bk1_1_2
    complex(kind=omega_prec) :: bk1_1_3
    complex(kind=omega_prec) :: bk1_1_4
    complex(kind=omega_prec) :: bk1_1_5
    complex(kind=omega_prec) :: bk1_1_6
    complex(kind=omega_prec) :: bk1_1_7
    type(momentum) :: p1, p2, p3, p4, p5, p6
    type(spinor) :: n2_6, d1_1__4, l1_1
    type(conjspinor) :: l2b_5, u1b__1_3, l1b_2
    type(spinor) :: n1_156_111c000a000, n2_126_111c000a000, &
      d1_1__124_111c000a004
    type(conjspinor) :: l2b_125_111c000a000, u1b__1_123_111c004a000
    type(vector) :: a_12_11c00a00, wp_56_11c00a00, wm_34_11c20a01, &
      z_12_11c00a00
    type(momentum) :: p12, p123, p124, p125, p126, p156, p34, p56
    complex(kind=omega_prec) :: bk1_2_1
    complex(kind=omega_prec) :: bk1_2_2
    complex(kind=omega_prec) :: bk1_2_3
    complex(kind=omega_prec) :: bk1_2_4
    complex(kind=omega_prec) :: bk1_2_5
    complex(kind=omega_prec) :: bk1_2_6
    complex(kind=omega_prec) :: bk1_2_7
    type(spinor) :: d2_1__4
    type(conjspinor) :: u2b__1_3
    type(spinor) :: d2_1__124_112c000a004
    type(conjspinor) :: u2b__1_123_112c004a000
    type(vector) :: wm_34_22c20a01
    complex(kind=omega_prec) :: bk2_3_1
    complex(kind=omega_prec) :: bk2_3_2
    complex(kind=omega_prec) :: bk2_3_3
    complex(kind=omega_prec) :: bk2_3_4
    complex(kind=omega_prec) :: bk2_3_5
    complex(kind=omega_prec) :: bk2_3_6
    complex(kind=omega_prec) :: bk2_3_7
    type(spinor) :: l2_5, u1_1__3
    type(conjspinor) :: n2b_6, d1b__1_4
    type(spinor) :: l2_125_112c000a000, n1_134_133c002a030, &
      u1_1__123_113c000a004
    type(conjspinor) :: n2b_126_112c000a000, d1b__1_124_113c004a000
    type(vector) :: wm_56_22c00a00, wp_34_33c01a20
    type(momentum) :: p134
    complex(kind=omega_prec) :: bk2_4_1
    complex(kind=omega_prec) :: bk2_4_2
    complex(kind=omega_prec) :: bk2_4_3
    complex(kind=omega_prec) :: bk2_4_4
    complex(kind=omega_prec) :: bk2_4_5
    complex(kind=omega_prec) :: bk2_4_6
    complex(kind=omega_prec) :: bk2_4_7
    type(spinor) :: u2_1__3
    type(conjspinor) :: d2b__1_4
    type(spinor) :: n1_134_144c002a030, u2_1__123_114c000a004
    type(conjspinor) :: d2b__1_124_114c004a000
    type(vector) :: wp_34_44c01a20
    integer, dimension(n_prt) :: s
    integer :: hi, ho, h
    p1 = - k(:,1) ! incoming e-
    p2 = - k(:,2) ! incoming e+
    p3 =   k(:,3) ! outgoing u/1/
    p4 =   k(:,4) ! outgoing dbar//1
    p5 =   k(:,5) ! outgoing mu-
    p6 =   k(:,6) ! outgoing numubar
    p12 = p1 + p2
    p34 = p3 + p4
    p56 = p5 + p6
    p123 = p3 + p12
    p124 = p4 + p12
    p125 = p5 + p12
    p126 = p6 + p12
    p156 = p1 + p56
    p134 = p1 + p34
    amp = 0
    flv_hel_flag = .false.
    h = 0
    do hi = 1, n_hel_in
    do ho = 1, n_hel_out
       h = h + 1
       if (zero_ct(hi,1,ho,1) > REPEAT)  cycle
       s = table_spin_states (:,h)
       ! flavor: [-11, 11, 2, -1, 13, -14]   color: (0, 0, 4, 0, 0, 0)
       !                                            (0, 0, 0, 3, 0, 0)
       l1_1 = u (mass(11), - p1, s(1))
       l1b_2 = vbar (mass(11), - p2, s(2))
       u1b__1_3 = ubar (mass(2), p3, s(3))
       d1_1__4 = v (mass(1), p4, s(4))
       l2b_5 = ubar (mass(13), p5, s(5))
       n2_6 = v (mass(14), p6, s(6))
       a_12_11c00a00 = pr_feynman(p12, + v_ff(qlep,l1b_2,l1_1))
       z_12_11c00a00 = &
         pr_unitarity(p12,mass(23),width(23), &
         + va_ff(gnclep(1),gnclep(2),l1b_2,l1_1))
       wm_34_11c20a01 = &
         pr_unitarity(p34,mass(24),width(24), + vl_ff(gcc,u1b__1_3,d1_1__4))
       wp_56_11c00a00 = &
         pr_unitarity(p56,mass(24),width(24), + vl_ff(gcc,l2b_5,n2_6))
       u1b__1_123_111c004a000 = &
         pr_psibar(p123,mass(2),width(2), &
         - f_fv(qup,u1b__1_3,a_12_11c00a00) &
         - f_fva(gncup(1),gncup(2),u1b__1_3,z_12_11c00a00))
       d1_1__124_111c000a004 = &
         pr_psi(p124,mass(1),width(1), - f_vf(qdwn,a_12_11c00a00,d1_1__4) &
         - f_vaf(gncdwn(1),gncdwn(2),z_12_11c00a00,d1_1__4))
       l2b_125_111c000a000 = &
         pr_psibar(p125,mass(13),width(13), &
         - f_fv(qlep,l2b_5,a_12_11c00a00) &
         - f_fva(gnclep(1),gnclep(2),l2b_5,z_12_11c00a00))
       n2_126_111c000a000 = &
         pr_psi(p126,mass(14),width(14), &
         - f_vaf(gncneu(1),gncneu(2),z_12_11c00a00,n2_6))
       n1_156_111c000a000 = &
         pr_psi(p156,mass(12),width(12), + f_vlf(gcc,wp_56_11c00a00,l1_1))
       bk1_1_1 = &
         + u1b__1_123_111c004a000*( - f_vlf(gcc,wp_56_11c00a00,d1_1__4))
       bk1_1_2 = &
         + ( - f_fvl(gcc,u1b__1_3,wp_56_11c00a00))*d1_1__124_111c000a004
       bk1_1_3 = &
         + z_12_11c00a00*( &
         + kg5_kgkg((-il5z),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + kg_kgkg((-ilz),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + l5kv_vv((-rg5z),wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + t5kv_vv((-ik5z),wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + lv_kvv(ig1mkmg4z,wp_56_11c00a00,p56,wm_34_11c20a01) &
         + lv_kvv((-ig1mkpg4z),wm_34_11c20a01,p34,wp_56_11c00a00) &
         + (-1)*tv_kvv((-ig1pkpg4z),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + (-1)*tv_kvv(ig1pkmg4z,wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + tkv_vv((-ig1z),wm_34_11c20a01,p34,wp_56_11c00a00,p56))
       bk1_1_4 = &
         + a_12_11c00a00*( &
         + kg5_kgkg((-il5a),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + kg_kgkg((-ila),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + l5kv_vv((-rg5a),wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + t5kv_vv((-ik5a),wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + lv_kvv(ig1mkmg4a,wp_56_11c00a00,p56,wm_34_11c20a01) &
         + lv_kvv((-ig1mkpg4a),wm_34_11c20a01,p34,wp_56_11c00a00) &
         + (-1)*tv_kvv((-ig1pkpg4a),wp_56_11c00a00,p56,wm_34_11c20a01,p34) &
         + (-1)*tv_kvv(ig1pkmg4a,wm_34_11c20a01,p34,wp_56_11c00a00,p56) &
         + tkv_vv((-ig1a),wm_34_11c20a01,p34,wp_56_11c00a00,p56))
       bk1_1_5 = + ( - f_fvl(gcc,l2b_5,wm_34_11c20a01))*n2_126_111c000a000
       bk1_1_6 = + ( + f_fvl(gcc,l1b_2,wm_34_11c20a01))*n1_156_111c000a000
       bk1_1_7 = + l2b_125_111c000a000*( - f_vlf(gcc,wm_34_11c20a01,n2_6))
       amp(1,h,1) = &
         + bk1_1_1 + bk1_1_2 + bk1_1_3 + bk1_1_4 + bk1_1_5 + bk1_1_6 + bk1_1_7
       amp(1,h,1) = - amp(1,h,1) ! 4 vertices, 3 propagators
       flv_hel_flag(h,1) = flv_hel_flag(h,1) .or. amp(1,h,1) /= 0
       ! unit symmetry factor
       ! Number of external adjoints: 0
       ! flavor: [-11, 11, 4, -3, 13, -14]   color: (0, 0, 4, 0, 0, 0)
       !                                            (0, 0, 0, 3, 0, 0)
       u2b__1_3 = ubar (mass(4), p3, s(3))
       d2_1__4 = v (mass(3), p4, s(4))
       wm_34_22c20a01 = &
         pr_unitarity(p34,mass(24),width(24), + vl_ff(gcc,u2b__1_3,d2_1__4))
       u2b__1_123_112c004a000 = &
         pr_psibar(p123,mass(4),width(4), &
         - f_fv(qup,u2b__1_3,a_12_11c00a00) &
         - f_fva(gncup(1),gncup(2),u2b__1_3,z_12_11c00a00))
       d2_1__124_112c000a004 = &
         pr_psi(p124,mass(3),width(3), - f_vf(qdwn,a_12_11c00a00,d2_1__4) &
         - f_vaf(gncdwn(1),gncdwn(2),z_12_11c00a00,d2_1__4))
       bk1_2_1 = &
         + u2b__1_123_112c004a000*( - f_vlf(gcc,wp_56_11c00a00,d2_1__4))
       bk1_2_2 = &
         + ( - f_fvl(gcc,u2b__1_3,wp_56_11c00a00))*d2_1__124_112c000a004
       bk1_2_3 = &
         + z_12_11c00a00*( &
         + kg5_kgkg((-il5z),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + kg_kgkg((-ilz),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + l5kv_vv((-rg5z),wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + t5kv_vv((-ik5z),wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + lv_kvv(ig1mkmg4z,wp_56_11c00a00,p56,wm_34_22c20a01) &
         + lv_kvv((-ig1mkpg4z),wm_34_22c20a01,p34,wp_56_11c00a00) &
         + (-1)*tv_kvv((-ig1pkpg4z),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + (-1)*tv_kvv(ig1pkmg4z,wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + tkv_vv((-ig1z),wm_34_22c20a01,p34,wp_56_11c00a00,p56))
       bk1_2_4 = &
         + a_12_11c00a00*( &
         + kg5_kgkg((-il5a),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + kg_kgkg((-ila),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + l5kv_vv((-rg5a),wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + t5kv_vv((-ik5a),wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + lv_kvv(ig1mkmg4a,wp_56_11c00a00,p56,wm_34_22c20a01) &
         + lv_kvv((-ig1mkpg4a),wm_34_22c20a01,p34,wp_56_11c00a00) &
         + (-1)*tv_kvv((-ig1pkpg4a),wp_56_11c00a00,p56,wm_34_22c20a01,p34) &
         + (-1)*tv_kvv(ig1pkmg4a,wm_34_22c20a01,p34,wp_56_11c00a00,p56) &
         + tkv_vv((-ig1a),wm_34_22c20a01,p34,wp_56_11c00a00,p56))
       bk1_2_5 = + ( - f_fvl(gcc,l2b_5,wm_34_22c20a01))*n2_126_111c000a000
       bk1_2_6 = + ( + f_fvl(gcc,l1b_2,wm_34_22c20a01))*n1_156_111c000a000
       bk1_2_7 = + l2b_125_111c000a000*( - f_vlf(gcc,wm_34_22c20a01,n2_6))
       amp(1,h,2) = &
         + bk1_2_1 + bk1_2_2 + bk1_2_3 + bk1_2_4 + bk1_2_5 + bk1_2_6 + bk1_2_7
       amp(1,h,2) = - amp(1,h,2) ! 4 vertices, 3 propagators
       flv_hel_flag(h,2) = flv_hel_flag(h,2) .or. amp(1,h,2) /= 0
       ! unit symmetry factor
       ! Number of external adjoints: 0
       ! flavor: [-11, 11, -2, 1, -13, 14]   color: (0, 0, 0, 3, 0, 0)
       !                                            (0, 0, 4, 0, 0, 0)
       u1_1__3 = v (mass(2), p3, s(3))
       d1b__1_4 = ubar (mass(1), p4, s(4))
       l2_5 = v (mass(13), p5, s(5))
       n2b_6 = ubar (mass(14), p6, s(6))
       wp_34_33c01a20 = &
         pr_unitarity(p34,mass(24),width(24), + vl_ff(gcc,d1b__1_4,u1_1__3))
       wm_56_22c00a00 = &
         pr_unitarity(p56,mass(24),width(24), + vl_ff(gcc,n2b_6,l2_5))
       u1_1__123_113c000a004 = &
         pr_psi(p123,mass(2),width(2), - f_vf(qup,a_12_11c00a00,u1_1__3) &
         - f_vaf(gncup(1),gncup(2),z_12_11c00a00,u1_1__3))
       d1b__1_124_113c004a000 = &
         pr_psibar(p124,mass(1),width(1), &
         - f_fv(qdwn,d1b__1_4,a_12_11c00a00) &
         - f_fva(gncdwn(1),gncdwn(2),d1b__1_4,z_12_11c00a00))
       n1_134_133c002a030 = &
         pr_psi(p134,mass(12),width(12), + f_vlf(gcc,wp_34_33c01a20,l1_1))
       l2_125_112c000a000 = &
         pr_psi(p125,mass(13),width(13), - f_vf(qlep,a_12_11c00a00,l2_5) &
         - f_vaf(gnclep(1),gnclep(2),z_12_11c00a00,l2_5))
       n2b_126_112c000a000 = &
         pr_psibar(p126,mass(14),width(14), &
         - f_fva(gncneu(1),gncneu(2),n2b_6,z_12_11c00a00))
       bk2_3_1 = &
         + d1b__1_124_113c004a000*( - f_vlf(gcc,wm_56_22c00a00,u1_1__3))
       bk2_3_2 = &
         + ( - f_fvl(gcc,d1b__1_4,wm_56_22c00a00))*u1_1__123_113c000a004
       bk2_3_3 = &
         + z_12_11c00a00*( &
         + kg5_kgkg((-il5z),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + kg_kgkg((-ilz),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + l5kv_vv((-rg5z),wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + t5kv_vv((-ik5z),wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + lv_kvv(ig1mkmg4z,wp_34_33c01a20,p34,wm_56_22c00a00) &
         + lv_kvv((-ig1mkpg4z),wm_56_22c00a00,p56,wp_34_33c01a20) &
         + (-1)*tv_kvv((-ig1pkpg4z),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + (-1)*tv_kvv(ig1pkmg4z,wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + tkv_vv((-ig1z),wm_56_22c00a00,p56,wp_34_33c01a20,p34))
       bk2_3_4 = &
         + a_12_11c00a00*( &
         + kg5_kgkg((-il5a),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + kg_kgkg((-ila),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + l5kv_vv((-rg5a),wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + t5kv_vv((-ik5a),wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + lv_kvv(ig1mkmg4a,wp_34_33c01a20,p34,wm_56_22c00a00) &
         + lv_kvv((-ig1mkpg4a),wm_56_22c00a00,p56,wp_34_33c01a20) &
         + (-1)*tv_kvv((-ig1pkpg4a),wp_34_33c01a20,p34,wm_56_22c00a00,p56) &
         + (-1)*tv_kvv(ig1pkmg4a,wm_56_22c00a00,p56,wp_34_33c01a20,p34) &
         + tkv_vv((-ig1a),wm_56_22c00a00,p56,wp_34_33c01a20,p34))
       bk2_3_5 = + ( + f_fvl(gcc,l1b_2,wm_56_22c00a00))*n1_134_133c002a030
       bk2_3_6 = + n2b_126_112c000a000*( - f_vlf(gcc,wp_34_33c01a20,l2_5))
       bk2_3_7 = + ( - f_fvl(gcc,n2b_6,wp_34_33c01a20))*l2_125_112c000a000
       amp(2,h,3) = &
         + bk2_3_1 + bk2_3_2 + bk2_3_3 + bk2_3_4 + bk2_3_5 + bk2_3_6 + bk2_3_7
       amp(2,h,3) = - amp(2,h,3) ! 4 vertices, 3 propagators
       flv_hel_flag(h,3) = flv_hel_flag(h,3) .or. amp(2,h,3) /= 0
       ! unit symmetry factor
       ! Number of external adjoints: 0
       ! flavor: [-11, 11, -4, 3, -13, 14]   color: (0, 0, 0, 3, 0, 0)
       !                                            (0, 0, 4, 0, 0, 0)
       u2_1__3 = v (mass(4), p3, s(3))
       d2b__1_4 = ubar (mass(3), p4, s(4))
       wp_34_44c01a20 = &
         pr_unitarity(p34,mass(24),width(24), + vl_ff(gcc,d2b__1_4,u2_1__3))
       u2_1__123_114c000a004 = &
         pr_psi(p123,mass(4),width(4), - f_vf(qup,a_12_11c00a00,u2_1__3) &
         - f_vaf(gncup(1),gncup(2),z_12_11c00a00,u2_1__3))
       d2b__1_124_114c004a000 = &
         pr_psibar(p124,mass(3),width(3), &
         - f_fv(qdwn,d2b__1_4,a_12_11c00a00) &
         - f_fva(gncdwn(1),gncdwn(2),d2b__1_4,z_12_11c00a00))
       n1_134_144c002a030 = &
         pr_psi(p134,mass(12),width(12), + f_vlf(gcc,wp_34_44c01a20,l1_1))
       bk2_4_1 = &
         + d2b__1_124_114c004a000*( - f_vlf(gcc,wm_56_22c00a00,u2_1__3))
       bk2_4_2 = &
         + ( - f_fvl(gcc,d2b__1_4,wm_56_22c00a00))*u2_1__123_114c000a004
       bk2_4_3 = &
         + z_12_11c00a00*( &
         + kg5_kgkg((-il5z),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + kg_kgkg((-ilz),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + l5kv_vv((-rg5z),wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + t5kv_vv((-ik5z),wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + lv_kvv(ig1mkmg4z,wp_34_44c01a20,p34,wm_56_22c00a00) &
         + lv_kvv((-ig1mkpg4z),wm_56_22c00a00,p56,wp_34_44c01a20) &
         + (-1)*tv_kvv((-ig1pkpg4z),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + (-1)*tv_kvv(ig1pkmg4z,wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + tkv_vv((-ig1z),wm_56_22c00a00,p56,wp_34_44c01a20,p34))
       bk2_4_4 = &
         + a_12_11c00a00*( &
         + kg5_kgkg((-il5a),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + kg_kgkg((-ila),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + l5kv_vv((-rg5a),wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + t5kv_vv((-ik5a),wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + lv_kvv(ig1mkmg4a,wp_34_44c01a20,p34,wm_56_22c00a00) &
         + lv_kvv((-ig1mkpg4a),wm_56_22c00a00,p56,wp_34_44c01a20) &
         + (-1)*tv_kvv((-ig1pkpg4a),wp_34_44c01a20,p34,wm_56_22c00a00,p56) &
         + (-1)*tv_kvv(ig1pkmg4a,wm_56_22c00a00,p56,wp_34_44c01a20,p34) &
         + tkv_vv((-ig1a),wm_56_22c00a00,p56,wp_34_44c01a20,p34))
       bk2_4_5 = + ( + f_fvl(gcc,l1b_2,wm_56_22c00a00))*n1_134_144c002a030
       bk2_4_6 = + n2b_126_112c000a000*( - f_vlf(gcc,wp_34_44c01a20,l2_5))
       bk2_4_7 = + ( - f_fvl(gcc,n2b_6,wp_34_44c01a20))*l2_125_112c000a000
       amp(2,h,4) = &
         + bk2_4_1 + bk2_4_2 + bk2_4_3 + bk2_4_4 + bk2_4_5 + bk2_4_6 + bk2_4_7
       amp(2,h,4) = - amp(2,h,4) ! 4 vertices, 3 propagators
       flv_hel_flag(h,4) = flv_hel_flag(h,4) .or. amp(2,h,4) /= 0
       ! unit symmetry factor
       ! Number of external adjoints: 0
       if (n <= SAMPLE) then
          if (all (abs(amp(:,h,:)) <= &
            tiny(1._omega_prec)))  zero_ct(hi,1,ho,1) = zero_ct(hi,1,ho,1) + 1
       end if
    end do
    end do
  end subroutine calculate_amplitudes

  function flv (fin, fout)
    integer, intent(in) :: fin, fout
    integer :: flv
    flv = 0
    select case (fin)
    case (1)
       select case (fout)
       case (1);  flv = 1   ! (/ 11, -11, 2, -1, 13, -14 /)
       case (2);  flv = 2   ! (/ 11, -11, 4, -3, 13, -14 /)
       case (3);  flv = 4   ! (/ 11, -11, -4, 3, -13, 14 /)
       case (4);  flv = 3   ! (/ 11, -11, -2, 1, -13, 14 /)
       end select
    end select
  end function flv

  function flv_pdf (fin)
    integer, intent(in) :: fin
    integer, dimension(2) :: flv_pdf
    flv_pdf = 0
    select case (fin)
    case (1)
       flv_pdf = (/  11, -11 /)
    end select
  end function flv_pdf

  subroutine scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    real(kind=omega_prec), dimension(:,:), intent(in) :: rho_in
    real(kind=omega_prec), dimension(:,:), intent(inout) :: rho_out
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer :: fi, fo, f, hi, ho, h, c
    complex(kind=omega_prec), dimension(n_cflow,n_hel,n_flv) :: amp
    complex(kind=omega_prec), dimension(n_cflow) :: amp1, amp2
    call calculate_amplitudes (amp, p, zero_ct, n)
    rho_out = 0
    h = 0
    do hi = 1, n_hel_in
    do ho = 1, n_hel_out
       h = h + 1
     !  print *,"h = ", h
       do fi = 1, n_flv_in
       do fo = 1, n_flv_out
          f = flv (fi, fo);  if (f == 0)  cycle
          if (.not. flv_hel_flag(h,f)) cycle
          amp1 = 0
         ! print *, flow_coeff(1:n_cflow,1:n_cflow)
         ! print *, flow_flag(1:n_cflow,1:n_cflow)
         ! print *, flv_col_flag(1:n_cflow,f)
          forall (c = 1:n_cflow, flv_col_flag(c,f))
            amp1(c) = &
             ! sum (amp(:,h,f) * flow_coeff(:,c), flow_flag(:,c) .and. &
              sum (amp(:,h,f), flow_flag(:,c) .and. &
              flv_col_flag(:,f))             
          end forall
          amp2 = amp(:,h,f)
         ! print *, amp1
         ! print *, amp2
         ! print *, dot_product (amp1, amp2)
         ! print *," "
          rho_out(ho,fo) = &
            rho_out(ho,fo) + dot_product (amp1, amp2) * rho_in(hi,fi)
       end do
       end do
    end do
    end do
  end subroutine scatter_diagonal_nonzero
  
  subroutine scatter_diagonal_polarized (p, rho_in, rho_out, zero_ct, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    real(kind=omega_prec), dimension(:), intent(in) :: rho_in
    real(kind=omega_prec), dimension(:,:,:), intent(inout) :: rho_out
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer :: fi, fo, f, hi, ho, h, c
    complex(kind=omega_prec), dimension(n_cflow,n_hel,n_flv) :: amp
    complex(kind=omega_prec), dimension(n_cflow) :: amp1, amp2
    call calculate_amplitudes (amp, p, zero_ct, n)
    rho_out = 0
    h = 0
    do hi = 1, n_hel_in
    do ho = 1, n_hel_out
       h = h + 1
     !  print *,"h = ", h
       do fi = 1, n_flv_in
       do fo = 1, n_flv_out
          f = flv (fi, fo);  if (f == 0)  cycle
          if (.not. flv_hel_flag(h,f)) cycle
          amp1 = 0
          forall (c = 1:n_cflow, flv_col_flag(c,f))
            amp1(c) = &
             ! sum (amp(:,h,f) * flow_coeff(:,c), flow_flag(:,c) .and. &
              sum (amp(:,h,f), flow_flag(:,c) .and. &
              flv_col_flag(:,f))             
          end forall
          amp2 = amp(:,h,f)
          rho_out(ho,hi,fo) = &
            rho_out(ho,hi,fo) + dot_product (amp1, amp2) * rho_in(fi)
       end do
       end do
    end do
    end do
  end subroutine scatter_diagonal_polarized  

  subroutine scatter_nonzero (p, rho_in, rho_out, zero_ct, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    complex(kind=omega_prec), dimension(:,:,:,:), intent(in) :: rho_in
    complex(kind=omega_prec), dimension(:,:,:,:), intent(inout) :: rho_out
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer :: fi, fo, f, hi1, hi2, ho1, ho2, h1, h2, c
    complex(kind=omega_prec), dimension(n_cflow,n_hel,n_flv) :: amp
    complex(kind=omega_prec), dimension(n_cflow) :: amp1, amp2
    call calculate_amplitudes (amp, p, zero_ct, n)
    rho_out = 0
    h1 = 0
    do hi1 = 1, n_hel_in
    do ho1 = 1, n_hel_out
       h1 = h1 + 1
       do hi2 = 1, n_hel_in
       do ho2 = 1, n_hel_out
          h2 = h2 + 1
          do fi = 1, n_flv_in
          do fo = 1, n_flv_out
             f = flv (fi, fo);  if (f == 0)  cycle
             if (.not. (flv_hel_flag(h1,f) .and. flv_hel_flag(h2,f))) cycle
             amp1 = 0
             forall (c = 1:n_cflow, flv_col_flag(c,f))
               amp1(c) = &
                 sum (amp(:,h1,f) * flow_coeff(:,c), flow_flag(:,c) .and. &
                 flv_col_flag(:,f))
             end forall
             amp2 = amp(:,h2,f)
             rho_out(ho1,fo,ho2,fo) = &
               rho_out(ho1,fo,ho2,fo) &
               + dot_product (amp1, amp2) * rho_in(hi1,fi,hi2,fi)
          end do
          end do
       end do
       end do
    end do
    end do
  end subroutine scatter_nonzero

  subroutine scatter_diagonal_colored_nz (p, rho_in, rho_out, zero_ct, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    real(kind=omega_prec), dimension(:,:), intent(in) :: rho_in
    real(kind=omega_prec), dimension(:,:,:), intent(inout) :: rho_out
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer :: fi, fo, f, hi, ho, h, c
    complex(kind=omega_prec), dimension(n_cflow,n_hel,n_flv) :: amp
    real(kind=omega_prec), dimension(n_cflow) :: pcol
    real(kind=omega_prec) :: sqme, sum_pcol
    complex(kind=omega_prec), dimension(n_cflow) :: amp1, amp2
    call calculate_amplitudes (amp, p, zero_ct, n)
    rho_out = 0
    pcol = 0
    h = 0
    do hi = 1, n_hel_in
    do ho = 1, n_hel_out
       h = h + 1
       do fi = 1, n_flv_in
       do fo = 1, n_flv_out
          f = flv (fi, fo);  if (f == 0)  cycle
          if (.not. flv_hel_flag(h,f)) cycle
          amp1 = 0
          forall (c = 1:n_cflow, flv_col_flag(c,f))
            amp1(c) = &
              sum (amp(:,h,f) * flow_coeff(:,c), flow_flag(:,c) .and. &
              flv_col_flag(:,f))
          end forall
          amp2 = amp(:,h,f)
          sqme = dot_product (amp1, amp2) * rho_in(hi,fi)
          where (flow_is_physical)
             pcol(:) = abs (amp(:,h,f))**2
          end where
          sum_pcol = sum (pcol)
          if (sum_pcol /= 0) then
             rho_out(:,ho,fo) = rho_out(:,ho,fo) + sqme * pcol(:) / sum_pcol
          end if
       end do
       end do
    end do
    end do
  end subroutine scatter_diagonal_colored_nz

  subroutine scatter_colored_nonzero (p, rho_in, rho_out, rho_col_out, &
    zero_ct, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    complex(kind=omega_prec), dimension(:,:,:,:), intent(in) :: rho_in
    complex(kind=omega_prec), dimension(:,:,:,:), intent(inout) :: rho_out
    complex(kind=omega_prec), dimension(:,:,:,:,:), intent(inout) :: &
      rho_col_out
    integer, dimension(:,:,:,:), intent(inout) :: zero_ct
    integer, intent(in) :: n
    integer :: fi, fo, f, hi1, hi2, ho1, ho2, h1, h2, c
    complex(kind=omega_prec), dimension(n_cflow,n_hel,n_flv) :: amp
    complex(kind=omega_prec), dimension(n_cflow) :: amp1, amp2
    call calculate_amplitudes (amp, p, zero_ct, n)
    rho_out = 0
    rho_col_out = 0
    h1 = 0
    do hi1 = 1, n_hel_in
    do ho1 = 1, n_hel_out
       h1 = h1 + 1
       do hi2 = 1, n_hel_in
       do ho2 = 1, n_hel_out
          h2 = h2 + 1
          do fi = 1, n_flv_in
          do fo = 1, n_flv_out
             f = flv (fi, fo);  if (f == 0)  cycle
             if (.not. (flv_hel_flag(h1,f) .and. flv_hel_flag(h2,f))) cycle
             amp1 = 0
             forall (c = 1:n_cflow, flv_col_flag(c,f))
               amp1(c) = &
                 sum (amp(:,h1,f) * flow_coeff(:,c), flow_flag(:,c) .and. &
                 flv_col_flag(:,f))
             end forall
             amp2 = amp(:,h2,f)
             rho_out(ho1,fo,ho2,fo) = &
               rho_out(ho1,fo,ho2,fo) &
               + dot_product (amp1, amp2) * rho_in(hi1,fi,hi2,fi)
             where (flow_is_physical)
                rho_col_out(:,ho1,fo,ho2,fo) = &
                  rho_col_out(:,ho1,fo,ho2,fo) &
                  + amp(:,h1,f) * conjg (amp2(:)) * rho_in(hi1,fi,hi2,fi)
             end where
          end do
          end do
       end do
       end do
    end do
    end do
  end subroutine scatter_colored_nonzero

end module ww_sl0muq
! O'Mega revision control information:
!    Colorize.Gauge(Models.SM):
!      Gauged Colorization Functor ( minimal electroweak standard model in unitarity gauge )
!      URL: /home/sources/ohl/ml/omega/src/models.ml,v 
!      revision: 326  checked in by reuter  at 2008-08-17 06:49:19 +0200 (Sun, 17 Aug 2008) 
!    Targets.Make_Fortran():
!      NB: non-gauge vector couplings are not available yet
!      URL: /home/sources/ohl/ml/omega/src/targets.ml,v 
!      revision: 326  checked in by reuter  at 2008-08-17 06:49:19 +0200 (Sun, 17 Aug 2008) 
!    Targets.Fortran_Fermions():
!      generates Fortran95 code for Dirac fermions
!      using revision 2000_10_A of module omega95
!      URL: /home/sources/ohl/ml/omega/src/targets.ml,v 
!      revision: 326  checked in by reuter  at 2008-08-17 06:49:19 +0200 (Sun, 17 Aug 2008) 
!    DAG.Graded():
!      Graded directed Acyclical Graph 
!      representing binary or n-ary trees
!      URL: svn+ssh://jr_reuter@login.hepforge.org/hepforge/svn/whizard/branches/1.xx/omega-src/bundle/src/dAG.ml 
!      revision: 68  checked in by ohl  at 2007-11-22 12:11:19 +0100 (Thu, 22 Nov 2007) 
!    Topology.Mixed23:
!      phi**3 + phi**4 topology
!      URL: svn+ssh://jr_reuter@login.hepforge.org/hepforge/svn/whizard/branches/1.xx/omega-src/bundle/src/topology.ml 
!      revision: 68  checked in by ohl  at 2007-11-22 12:11:19 +0100 (Thu, 22 Nov 2007) 
!    Momentum.Bits():
!      Finite disjoint sums of momenta
!      using bitfields as representation.
!      URL: svn+ssh://jr_reuter@login.hepforge.org/hepforge/svn/whizard/branches/1.xx/omega-src/bundle/src/momentum.ml 
!      revision: 68  checked in by ohl  at 2007-11-22 12:11:19 +0100 (Thu, 22 Nov 2007) 
!    Fusion.Make():
!      Fusions for arbitrary topologies
!      URL: svn+ssh://jr_reuter@login.hepforge.org/hepforge/svn/whizard/branches/1.xx/omega-src/bundle/src/fusion.ml 
!      revision: 326  checked in by reuter  at 2008-08-17 06:49:19 +0200 (Sun, 17 Aug 2008) 
