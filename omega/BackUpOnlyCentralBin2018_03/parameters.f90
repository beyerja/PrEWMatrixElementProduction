! WHIZARD parameter definition file
! Automatically generated code, do not edit

module parameters

  use kinds, only: default  !NODEP!
  use file_utils, only: free_unit


  implicit none
  private

  public :: read, write, init
  public :: QED_charge, SM_sw, SM_cw

  type, public :: parameter_set
    logical :: gwidth, rwidth
    real(kind=default) :: sqrt2
include "par_def.inc"
  end type parameter_set
  
include "par_val.inc"

  interface read
     module procedure read_parameters_unit
     module procedure read_parameters_name
  end interface
  interface write
     module procedure write_parameters_unit
     module procedure write_parameters_name
  end interface
  interface init
     module procedure init_parameters
  end interface

contains

  function QED_charge (par) result (ee)
    type(parameter_set), intent(in) :: par
    real(kind=default) :: ee
include "par_qed.inc"
  end function QED_charge

  function SM_cw (par) result (cw)
    type(parameter_set), intent(in) :: par
    real(kind=default) :: sw, cw
include "par_sw.inc"
  end function SM_cw

  function SM_sw (par) result (sw)
    type(parameter_set), intent(in) :: par
    real(kind=default) :: sw, cw
include "par_sw.inc"
  end function SM_sw

  subroutine read_parameters_unit (unit, par)
    integer, intent(in) :: unit
    type(parameter_set), intent(inout) :: par

    logical :: gwidth, rwidth
    real(kind=default) :: sqrt2
include "par_def.inc"
include "par_nml.inc"
    gwidth = par%gwidth
    rwidth = par%rwidth
    sqrt2 = par%sqrt2
include "par_get.inc"
    read (unit, nml=parameter_input)
    par%gwidth = gwidth
    par%rwidth = rwidth
include "par_set.inc"

  end subroutine read_parameters_unit

  subroutine write_parameters_unit (unit, par)
    integer, intent(in) :: unit
    type(parameter_set), intent(in) :: par
    write(unit, *) '&parameter_input'
    write(unit, *) 'gwidth   = ', par%gwidth
    write(unit, *) 'rwidth   = ', par%rwidth
include "par_write.inc"
  end subroutine write_parameters_unit

  subroutine read_parameters_name (name, par)
    character(len=*), intent(in) :: name
    type(parameter_set), intent(inout) :: par
    integer :: unit
    unit = free_unit ()
    open (unit=unit, action='read', status='old', file=name)
    call read_parameters_unit (unit, par)
    close (unit=unit)
  end subroutine read_parameters_name

  subroutine write_parameters_name (name, par)
    character(len=*), intent(in) :: name
    type(parameter_set), intent(in) :: par
    integer :: unit
    unit = free_unit ()
    open (unit=unit, action='write', status='replace', file=name)
    call write_parameters_unit (unit, par)
    close (unit=unit)
  end subroutine write_parameters_name

  subroutine init_parameters (par)
    type(parameter_set), intent(inout) :: par
    par = parameters_default
  end subroutine init_parameters

end module parameters
