! $Id: omegalib.nw 1148 2009-09-02 13:58:39Z jr_reuter $
!
! Copyright (C) 2000-2002 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
! and others
!
! O'Mega is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! O'Mega is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
module omega_utils
  use omega_kinds
  use omega_vectors
  use omega_polarizations
  implicit none
  private
  public :: omega_ward_warn, omega_ward_panic
  public :: omega_slavnov_warn, omega_slavnov_panic
  public :: omega_check_arguments_warn, omega_check_arguments_panic
  public :: omega_check_helicities_warn, omega_check_helicities_panic
  private :: omega_check_helicity
  public :: omega_check_momenta_warn, omega_check_momenta_panic
  private :: check_momentum_conservation, check_mass_shell
  public :: omega_spin_sum_sqme_1, omega_sum_sqme
  public :: omega_spin_sum_sqme_1_nonzero, omega_sum_sqme_nonzero
  public :: omega_amplitude_1_nonzero, omega_amplitude_2_nonzero
  public :: omega_scatter, omega_scatter_nonzero
  public :: omega_scatter_diagonal, omega_scatter_diagonal_nonzero
  public :: omega_sum, omega_sum_nonzero, omega_nonzero
  private :: state_index
  public :: num_states
  integer, parameter, private :: MOMENTUM_TOLERANCE = 10000
  integer, parameter, private :: ON_SHELL_TOLERANCE = 1000000
  integer, parameter, private :: REPEAT = 5, SAMPLE = 10
  integer, parameter, public :: omega_utils_2003_03_A = 0
contains
  subroutine omega_ward_warn (name, m, k, e)
    character(len=*), intent(in) :: name
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    type(vector) :: ek
    real(kind=omega_prec) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e)
    abs_ek_abs_e = abs (ek) * abs (e)
    print *, name, ":", abs_eke / abs_ek_abs_e, abs (ek), abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: warning: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
    end if
  end subroutine omega_ward_warn
  subroutine omega_ward_panic (name, m, k, e)
    character(len=*), intent(in) :: name
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    type(vector) :: ek
    real(kind=omega_prec) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e)
    abs_ek_abs_e = abs (ek) * abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: panic: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
       stop
    end if
  end subroutine omega_ward_panic
  subroutine omega_slavnov_warn (name, m, k, e, phi)
    character(len=*), intent(in) :: name
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    complex(kind=omega_prec), intent(in) :: phi
    type(vector) :: ek
    real(kind=omega_prec) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e - phi)
    abs_ek_abs_e = abs (ek) * abs (e)
    print *, name, ":", abs_eke / abs_ek_abs_e, abs (ek), abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: warning: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
    end if
  end subroutine omega_slavnov_warn
  subroutine omega_slavnov_panic (name, m, k, e, phi)
    character(len=*), intent(in) :: name
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    complex(kind=omega_prec), intent(in) :: phi
    type(vector) :: ek
    real(kind=omega_prec) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e - phi)
    abs_ek_abs_e = abs (ek) * abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: panic: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
       stop
    end if
  end subroutine omega_slavnov_panic
  subroutine omega_check_arguments_warn (n, k, s)
    integer, intent(in) :: n
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, dimension(:), intent(in) :: s
    integer :: i
    i = size(k,dim=1)
    if (i /= 4) then
       print *, "O'Mega: warning: wrong # of dimensions:", i
    end if
    i = size(k,dim=2)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of momenta:", i, &
            ", expected", n
    end if
    i = size (s)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of spins:", i, &
            ", expected", n
    end if
  end subroutine omega_check_arguments_warn
  subroutine omega_check_arguments_panic (n, k, s)
    integer, intent(in) :: n
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, dimension(:), intent(in) :: s
    logical :: error
    integer :: i
    error = .false.
    i = size(k,dim=1)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of dimensions:", i
       error = .true.
    end if
    i = size(k,dim=2)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of momenta:", i, &
            ", expected", n
       error = .true.
    end if
    i = size (s)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of spins:", i, &
            ", expected", n
       error = .true.
    end if
    if (error) then
       stop
    end if
  end subroutine omega_check_arguments_panic
  function omega_check_helicity (m, smax, s) result (error)
    real(kind=omega_prec), intent(in) :: m
    integer, intent(in) :: smax, s
    logical :: error
    select case (smax)
    case (0)
       error = (s /= 0)
    case (1)
       error = (abs (s) /= 1)
    case (2)
       if (m == 0.0_omega_prec) then
          error = .not. (abs (s) == 1 .or. abs (s) == 4)
       else
          error = .not. (abs (s) <= 1 .or. abs (s) == 4)
       end if
    case (4)
       error = .true.
    case default
       error = .true.
    end select
  end function omega_check_helicity
  subroutine omega_check_helicities_warn (m, smax, s)
    real(kind=omega_prec), dimension(:), intent(in) :: m
    integer, dimension(:), intent(in) :: smax, s
    integer :: i
    do i = 1, size (m)
       if (omega_check_helicity (m(i), smax(i), s(i))) then
          print *, "O'Mega: warning: invalid helicity", s(i)
       end if
    end do
  end subroutine omega_check_helicities_warn
  subroutine omega_check_helicities_panic (m, smax, s)
    real(kind=omega_prec), dimension(:), intent(in) :: m
    integer, dimension(:), intent(in) :: smax, s
    logical :: error
    logical :: error1
    integer :: i
    error = .false.
    do i = 1, size (m)
       error1 = omega_check_helicity (m(i), smax(i), s(i))
       if (error1) then
          print *, "O'Mega: panic: invalid helicity", s(i)
          error = .true.
       end if
    end do
    if (error) then
       stop
    end if
  end subroutine omega_check_helicities_panic
  function check_momentum_conservation (k) result (error)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    logical :: error 
    error = any (abs (sum (k(:,3:), dim = 2) - k(:,1) - k(:,2)) > &
         MOMENTUM_TOLERANCE * epsilon (maxval (abs (k), dim = 2)))
    if (error) then
       print *, sum (k(:,3:), dim = 2) - k(:,1) - k(:,2)
       print *, MOMENTUM_TOLERANCE * epsilon (maxval (abs (k), dim = 2)), &
            maxval (abs (k), dim = 2)
    end if
  end function check_momentum_conservation
  function check_mass_shell (m, k) result (error)
    real(kind=omega_prec), intent(in) :: m
    real(kind=omega_prec), dimension(0:), intent(in) :: k
    real(kind=omega_prec) :: e2
    logical :: error
    e2 = k(1)**2 + k(2)**2 + k(3)**2 + m**2
    error = abs (k(0)**2 - e2) > ON_SHELL_TOLERANCE * epsilon (max (k(0)**2, e2))
    if (error) then
       print *, k(0)**2 - e2
       print *, ON_SHELL_TOLERANCE * epsilon (max (k(0)**2, e2)), max (k(0)**2, e2)
    end if
  end function check_mass_shell
  subroutine omega_check_momenta_warn (m, k)
    real(kind=omega_prec), dimension(:), intent(in) :: m
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer :: i
    if (check_momentum_conservation (k)) then
       print *, "O'Mega: warning: momentum not conserved"
    end if
    do i = 1, size(m)
       if (check_mass_shell (m(i), k(:,i))) then
          print *, "O'Mega: warning: particle #", i, "not on-shell"
       end if
    end do
  end subroutine omega_check_momenta_warn
  subroutine omega_check_momenta_panic (m, k)
    real(kind=omega_prec), dimension(:), intent(in) :: m
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    logical :: error
    logical :: error1
    integer :: i
    error = check_momentum_conservation (k)
    if (error) then
       print *, "O'Mega: panic: momentum not conserved"
    end if
    do i = 1, size(m)
       error1 = check_mass_shell (m(i), k(0:,i))
       if (error1) then
          print *, "O'Mega: panic: particle #", i, "not on-shell"
          error = .true.
       end if
    end do
    if (error) then
       stop
    end if
  end subroutine omega_check_momenta_panic
  pure function omega_spin_sum_sqme_1 &
         (amplitude_1, k, f, s_max, smask) result (amp2)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, intent(in) :: f, s_max
    logical, dimension(:), intent(in), optional :: smask
    real(kind=omega_prec) :: amp2
    interface
      pure function amplitude_1 (k, s, f) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s, f
        complex(kind=omega_prec) :: amp
      end function amplitude_1
    end interface
    complex(kind=omega_prec) :: amp
    integer :: s
    amp2 = 0
    if (present (smask)) then
      do s = 1, s_max
        if (smask(s)) then
          amp = amplitude_1 (k, s, f)
          amp2 = amp2 + amp * conjg (amp)
        end if
      end do
    else
      do s = 1, s_max
        amp = amplitude_1 (k, s, f)
        amp2 = amp2 + amp * conjg (amp)
      end do
    end if
  end function omega_spin_sum_sqme_1
  pure function omega_sum_sqme &
         (amplitude_1, k, s_max, f_max, mult, smask, fmask) result (amp2)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, intent(in) :: s_max, f_max
    integer, dimension(:), intent(in) :: mult
    logical, dimension(:), intent(in), optional :: smask, fmask
    real(kind=omega_prec) :: amp2
    interface
      pure function amplitude_1 (k, s, f) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s, f
        complex(kind=omega_prec) :: amp
      end function amplitude_1
    end interface
    complex(kind=omega_prec) :: amp
    integer :: s, f
    amp2 = 0
    if (present (smask)) then
      if (present (fmask)) then
        do s = 1, s_max
          if (smask(s)) then
            do f = 1, f_max
              if (fmask(f)) then
                amp = amplitude_1 (k, s, f)
                amp2 = amp2 + amp * conjg (amp) / mult(f)
              end if
            end do
          end if
        end do
      else
        do s = 1, s_max
          if (smask(s)) then
            do f = 1, f_max
              amp = amplitude_1 (k, s, f)
              amp2 = amp2 + amp * conjg (amp) / mult(f)
            end do
          end if
        end do
      end if
    else
      if (present (fmask)) then
        do f = 1, f_max
          if (fmask(f)) then
            do s = 1, s_max
              amp = amplitude_1 (k, s, f)
              amp2 = amp2 + amp * conjg (amp) / mult(f)
            end do
          end if
        end do
      else
        do s = 1, s_max
          do f = 1, f_max
            amp = amplitude_1 (k, s, f)
            amp2 = amp2 + amp * conjg (amp) / mult(f)
          end do
        end do
      end if
    end if
  end function omega_sum_sqme
  pure subroutine omega_spin_sum_sqme_1_nonzero &
         (amplitude_1, amp2, k, f, zero, n, smask)
    real(kind=omega_prec), intent(out) :: amp2
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, intent(in) :: f
    integer, dimension(:,:), intent(inout) :: zero
    integer, intent(in) :: n
    logical, dimension(:), intent(in), optional :: smask
    interface
      pure function amplitude_1 (k, s, f) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s, f
        complex(kind=omega_prec) :: amp
      end function amplitude_1
    end interface
    complex(kind=omega_prec) :: amp
    real(kind=omega_prec) :: dummy
    integer :: s, i
    if (n <= SAMPLE) then
      call omega_sum_sqme_nonzero &
        (amplitude_1, dummy, k, (/ (1, i = 1, size(zero,dim=2)) /), zero, n)
    end if
    amp2 = 0
    if (present (smask)) then
      do s = 1, size(zero,dim=1)
        if (smask(s)) then
          if (zero(s,f) <= REPEAT) then
            amp = amplitude_1 (k, s, f)
            amp2 = amp2 + amp * conjg (amp)
          end if
        end if
      end do
    else
      do s = 1, size(zero,dim=1)
        if (zero(s,f) <= REPEAT) then
          amp = amplitude_1 (k, s, f)
          amp2 = amp2 + amp * conjg (amp)
        end if
      end do
    end if
  end subroutine omega_spin_sum_sqme_1_nonzero
  pure subroutine omega_sum_sqme_nonzero &
         (amplitude_1, amp2, k, mult, zero, n, smask, fmask)
    real(kind=omega_prec), intent(out) :: amp2
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, dimension(:), intent(in) :: mult
    integer, dimension(:,:), intent(inout) :: zero
    integer, intent(in) :: n
    logical, dimension(:), intent(in), optional :: smask, fmask
    interface
      pure function amplitude_1 (k, s, f) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s, f
        complex(kind=omega_prec) :: amp
      end function amplitude_1
    end interface
    complex(kind=omega_prec) :: amp
    integer :: s, f
    if (n <= SAMPLE) then
      do s = 1, size(zero,dim=1)
        do f = 1, size(zero,dim=2)
          if (zero(s,f) <= REPEAT) then
            amp = amplitude_1 (k, s, f)
            if (real (amp * conjg (amp), kind=omega_prec) &
                 <= tiny (1.0_omega_prec)) then
              zero(s,f) = zero(s,f) + 1
            end if
          end if
        end do
      end do
    end if
    amp2 = 0
    if (present (smask)) then
      if (present (fmask)) then
        do s = 1, size(zero,dim=1)
          if (smask(s)) then
            do f = 1, size(zero,dim=2)
              if (fmask(f)) then
                if (zero(s,f) <= REPEAT) then
                  amp = amplitude_1 (k, s, f)
                  amp2 = amp2 + amp * conjg (amp) / mult(f)
                end if
              end if
            end do
          end if
        end do
      else
        do s = 1, size(zero,dim=1)
          if (smask(s)) then
            do f = 1, size(zero,dim=2)
              if (zero(s,f) <= REPEAT) then
                amp = amplitude_1 (k, s, f)
                amp2 = amp2 + amp * conjg (amp) / mult(f)
              end if
            end do
          end if
        end do
      end if
    else
      if (present (fmask)) then
        do f = 1, size(zero,dim=2)
          if (fmask(f)) then
            do s = 1, size(zero,dim=1)
              if (zero(s,f) <= REPEAT) then
                amp = amplitude_1 (k, s, f)
                amp2 = amp2 + amp * conjg (amp) / mult(f)
              end if
            end do
          end if
        end do
      else
        do s = 1, size(zero,dim=1)
          do f = 1, size(zero,dim=2)
            if (zero(s,f) <= REPEAT) then
              amp = amplitude_1 (k, s, f)
              amp2 = amp2 + amp * conjg (amp) / mult(f)
            end if
          end do
        end do
      end if
    end if
  end subroutine omega_sum_sqme_nonzero
  pure subroutine omega_amplitude_1_nonzero &
         (amplitude_1, amp, k, s, f, zero, n)
    complex(kind=omega_prec), intent(out) :: amp
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, intent(in) :: s, f
    integer, dimension(:,:), intent(inout) :: zero
    integer, intent(in) :: n
    interface
      pure function amplitude_1 (k, s, f) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s, f
        complex(kind=omega_prec) :: amp
      end function amplitude_1
    end interface
    integer :: i
    real(kind=omega_prec) :: dummy
    if (n <= SAMPLE) then
      call omega_sum_sqme_nonzero &
        (amplitude_1, dummy, k, (/ (1, i = 1, size(zero,dim=2)) /), zero, n)
    end if
    if (zero(s,f) < REPEAT) then
      amp = amplitude_1 (k, s, f)
    else
      amp = 0
    end if
  end subroutine omega_amplitude_1_nonzero
  pure subroutine omega_amplitude_2_nonzero &
       (amplitude_2, amp, k, s_in, f_in, s_out, f_out, zero, n)
    complex(kind=omega_prec), intent(out) :: amp
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    integer, intent(in) :: s_in, f_in, s_out, f_out
    integer, dimension(:,:,:,:), intent(inout) :: zero
    integer, intent(in) :: n
    interface
      pure function amplitude_2 (k, s_in, f_in, s_out, f_out) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s_in, f_in, s_out, f_out
        complex(kind=omega_prec) :: amp
      end function amplitude_2
    end interface
    integer :: si, fi, so, fo
    if (n <= SAMPLE) then
      do si = 1, size(zero,dim=1)
        do fi = 1, size(zero,dim=2)
          do so = 1, size(zero,dim=3)
            do fo = 1, size(zero,dim=4)
              if (zero(si,fi,so,fo) <= REPEAT) then
                amp = amplitude_2 (k, si, fi, so, fo)
                if (real (amp * conjg (amp), kind=omega_prec) &
                     <= tiny (1.0_omega_prec)) then
                  zero(si,fi,so,fo) = zero(si,fi,so,fo) + 1
                end if
              end if
            end do
          end do
        end do
      end do
    end if
    if (zero(s_in,f_in,s_out,f_out) < REPEAT) then
      amp = amplitude_2 (k, s_in, f_in, s_out, f_out)
    else
      amp = 0
    end if
  end subroutine omega_amplitude_2_nonzero
  pure subroutine omega_scatter (amplitude_2, k, rho_in, rho_out, mult)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    complex(kind=omega_prec), dimension(:,:,:,:), intent(in) :: rho_in
    complex(kind=omega_prec), dimension(:,:,:,:), intent(inout) :: rho_out
    integer, dimension(:), intent(in) :: mult
    interface
      pure function amplitude_2 (k, s_in, f_in, s_out, f_out) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s_in, f_in, s_out, f_out
        complex(kind=omega_prec) :: amp
      end function amplitude_2
    end interface
    integer :: s_in1, s_in2, f_in1, f_in2, s_out1, s_out2, f_out1, f_out2
    complex(kind=omega_prec), &
      dimension(size(rho_in,dim=1),size(rho_in,dim=2),&
                size(rho_out,dim=1),size(rho_out,dim=2)) :: a
    do s_in1 = 1, size(rho_in,dim=1)
      do f_in1 = 1, size(rho_in,dim=2)
        do s_out1 = 1, size(rho_out,dim=1)
          do f_out1 = 1, size(rho_out,dim=2)
            a(s_in1,f_in1,s_out1,f_out1) = &
              amplitude_2 (k, s_in1, f_in1, s_out1, f_out1) &
                / sqrt (real (mult(f_out1), kind=omega_prec))
          end do
        end do
      end do
    end do
    do s_out1 = 1, size(rho_out,dim=1)
      do f_out1 = 1, size(rho_out,dim=2)
        do s_out2 = 1, size(rho_out,dim=3)
          do f_out2 = 1, size(rho_out,dim=4)
            rho_out(s_out1,f_out1,s_out2,f_out2) = 0
            do s_in1 = 1, size(rho_in,dim=1)
              do f_in1 = 1, size(rho_in,dim=2)
                do s_in2 = 1, size(rho_in,dim=3)
                  do f_in2 = 1, size(rho_in,dim=4)
                    rho_out(s_out1,f_out1,s_out2,f_out2) = &
                      rho_out(s_out1,f_out1,s_out2,f_out2) &
                        + a(s_in1,f_in1,s_out1,f_out1) &
                            * rho_in(s_in1,f_in1,s_in2,f_in2) &
                            * conjg (a(s_in2,f_in2,s_out2,f_out2))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine omega_scatter
  pure subroutine omega_scatter_nonzero &
         (amplitude_2, k, rho_in, rho_out, mult, zero, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    complex(kind=omega_prec), dimension(:,:,:,:), intent(in) :: rho_in
    complex(kind=omega_prec), dimension(:,:,:,:), intent(inout) :: rho_out
    integer, dimension(:), intent(in) :: mult
    integer, dimension(:,:,:,:), intent(inout) :: zero
    integer, intent(in) :: n
    interface
      pure subroutine amplitude_2 (amp, k, s_in, f_in, s_out, f_out, zero, n)
        use omega_kinds
        implicit none
        complex(kind=omega_prec), intent(out) :: amp
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s_in, f_in, s_out, f_out
        integer, dimension(:,:,:,:), intent(inout) :: zero
        integer, intent(in) :: n
      end subroutine amplitude_2
    end interface
    integer :: s_in1, s_in2, f_in1, f_in2, s_out1, s_out2, f_out1, f_out2
    complex(kind=omega_prec), &
      dimension(size(rho_in,dim=1),size(rho_in,dim=2),&
                size(rho_out,dim=1),size(rho_out,dim=2)) :: a
    do s_in1 = 1, size(rho_in,dim=1)
      do f_in1 = 1, size(rho_in,dim=2)
        do s_out1 = 1, size(rho_out,dim=1)
          do f_out1 = 1, size(rho_out,dim=2)
            call amplitude_2 (a(s_in1,f_in1,s_out1,f_out1), &
                              k, s_in1, f_in1, s_out1, f_out1, zero, n)
            a(s_in1,f_in1,s_out1,f_out1) = &
               a(s_in1,f_in1,s_out1,f_out1) &
                 / sqrt (real (mult(f_out1), kind=omega_prec))
          end do
        end do
      end do
    end do
    do s_out1 = 1, size(rho_out,dim=1)
      do f_out1 = 1, size(rho_out,dim=2)
        do s_out2 = 1, size(rho_out,dim=3)
          do f_out2 = 1, size(rho_out,dim=4)
            rho_out(s_out1,f_out1,s_out2,f_out2) = 0
            do s_in1 = 1, size(rho_in,dim=1)
              do f_in1 = 1, size(rho_in,dim=2)
                do s_in2 = 1, size(rho_in,dim=3)
                  do f_in2 = 1, size(rho_in,dim=4)
                    rho_out(s_out1,f_out1,s_out2,f_out2) = &
                      rho_out(s_out1,f_out1,s_out2,f_out2) &
                        + a(s_in1,f_in1,s_out1,f_out1) &
                            * rho_in(s_in1,f_in1,s_in2,f_in2) &
                            * conjg (a(s_in2,f_in2,s_out2,f_out2))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine omega_scatter_nonzero
  pure subroutine omega_scatter_diagonal &
         (amplitude_2, k, rho_in, rho_out, mult)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    real(kind=omega_prec), dimension(:,:), intent(in) :: rho_in
    real(kind=omega_prec), dimension(:,:), intent(inout) :: rho_out
    integer, dimension(:), intent(in) :: mult
    interface
      pure function amplitude_2 (k, s_in, f_in, s_out, f_out) result (amp)
        use omega_kinds
        implicit none
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s_in, f_in, s_out, f_out
        complex(kind=omega_prec) :: amp
      end function amplitude_2
    end interface
    integer :: s_in, f_in, s_out, f_out
    complex(kind=omega_prec) :: a
    do s_out = 1, size(rho_out,dim=1)
      do f_out = 1, size(rho_out,dim=2)
        rho_out(s_out,f_out) = 0
        do s_in = 1, size(rho_in,dim=1)
          do f_in = 1, size(rho_in,dim=2)
            a = amplitude_2 (k, s_in, f_in, s_out, f_out)
            rho_out(s_out,f_out) = rho_out(s_out,f_out) &
              + rho_in(s_in,f_in) * real (a*conjg(a), kind=omega_prec) &
                                  / mult(f_out)
          end do
        end do
      end do
    end do
  end subroutine omega_scatter_diagonal
  pure subroutine omega_scatter_diagonal_nonzero &
       (amplitude_2, k, rho_in, rho_out, mult, zero, n)
    real(kind=omega_prec), dimension(0:,:), intent(in) :: k
    real(kind=omega_prec), dimension(:,:), intent(in) :: rho_in
    real(kind=omega_prec), dimension(:,:), intent(inout) :: rho_out
    integer, dimension(:), intent(in) :: mult
    integer, dimension(:,:,:,:), intent(inout) :: zero
    integer, intent(in) :: n
    interface
      pure subroutine amplitude_2 (amp, k, s_in, f_in, s_out, f_out, zero, n)
        use omega_kinds
        implicit none
        complex(kind=omega_prec), intent(out) :: amp
        real(kind=omega_prec), dimension(0:,:), intent(in) :: k
        integer, intent(in) :: s_in, f_in, s_out, f_out
        integer, dimension(:,:,:,:), intent(inout) :: zero
        integer, intent(in) :: n
      end subroutine amplitude_2
    end interface
    integer :: s_in, f_in, s_out, f_out
    complex(kind=omega_prec) :: a
    do s_out = 1, size(rho_out,dim=1)
      do f_out = 1, size(rho_out,dim=2)
        rho_out(s_out,f_out) = 0
        do s_in = 1, size(rho_in,dim=1)
          do f_in = 1, size(rho_in,dim=2)
            call amplitude_2 (a, k, s_in, f_in, s_out, f_out, zero, n)
            rho_out(s_out,f_out) = rho_out(s_out,f_out) &
               + rho_in(s_in,f_in) * real (a*conjg(a), kind=omega_prec) &
                                   / mult(f_out)
          end do
        end do
      end do
    end do
  end subroutine omega_scatter_diagonal_nonzero
  pure function omega_sum (omega, p, states, fixed) result (sigma)
    real(kind=omega_prec) :: sigma
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    integer, dimension(:), intent(in), optional :: states, fixed
    interface
       pure function omega (p, s) result (me)
         use omega_kinds
         implicit none
         complex(kind=omega_prec) :: me
         real(kind=omega_prec), dimension(0:,:), intent(in) :: p
         integer, dimension(:), intent(in) :: s
       end function omega
    end interface
    integer, dimension(size(p,dim=2)) :: s, nstates
    integer :: j
    complex(kind=omega_prec) :: a
    if (present (states)) then
       nstates = states
    else
       nstates = 2
    end if
    sigma = 0
    s = -1
    sum_spins: do
       if (present (fixed)) then
          !!! print *, 's = ', s, ', fixed = ', fixed, ', nstates = ', nstates, &
          !!!    ', fixed|s = ', merge (fixed, s, mask = nstates == 0)
          a = omega (p, merge (fixed, s, mask = nstates == 0))
       else
          a = omega (p, s)
       end if
       sigma = sigma + a * conjg(a)
       do j = size (p, dim = 2), 1, -1
          select case (nstates (j))
          case (3) ! massive vectors
             s(j) = modulo (s(j) + 2, 3) - 1
          case (2) ! spinors, massless vectors
             s(j) = - s(j)
          case (1) ! scalars
             s(j) = -1
          case (0) ! fized spin
             s(j) = -1
          case default ! ???
             s(j) = -1
          end select
          if (s(j) /= -1) then
             cycle sum_spins
          end if
       end do
       exit sum_spins
    end do sum_spins
    sigma = sigma / num_states (2, nstates(1:2))
  end function omega_sum
  pure function state_index (s, states) result (n)
    integer, dimension(:), intent(in) :: s
    integer, dimension(:), intent(in), optional :: states
    integer :: n
    integer :: j, p
    n = 1
    p = 1
    if (present (states)) then
       do j = size (s), 1, -1
          select case (states(j))
          case (3)
             n = n + p * (s(j) + 1) 
          case (2)
             n = n + p * (s(j) + 1) / 2
          end select
          p = p * states(j)
       end do
    else
       do j = size (s), 1, -1
          n = n + p * (s(j) + 1) / 2
          p = p * 2
       end do
    end if
  end function state_index
  pure subroutine omega_sum_nonzero (sigma, omega, p, zero, n, states, fixed)
    real(kind=omega_prec), intent(out) :: sigma
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    integer, dimension(:), intent(inout) :: zero
    integer, intent(in) :: n
    integer, dimension(:), intent(in), optional :: states, fixed
    interface
       pure function omega (p, s) result (me)
         use omega_kinds
         implicit none
         complex(kind=omega_prec) :: me
         real(kind=omega_prec), dimension(0:,:), intent(in) :: p
         integer, dimension(:), intent(in) :: s
       end function omega
    end interface
    integer, dimension(size(p,dim=2)) :: s, nstates
    integer :: j, k
    complex(kind=omega_prec) :: a
    real(kind=omega_prec) :: a2
    if (present (states)) then
       nstates = states
    else
       nstates = 2
    end if
    sigma = 0
    s = -1
    k = 1
    sum_spins: do
       if (zero (k) < REPEAT) then
          if (present (fixed)) then
             a = omega (p, merge (fixed, s, mask = nstates == 0))
          else
             a = omega (p, s)
          end if
          a2 = a * conjg(a)
          if (n <= SAMPLE .and. a2 <= tiny (1.0_omega_prec)) then
             zero (k) = zero (k) + 1
          end if
          sigma = sigma + a2
       end if
       k = k + 1
       do j = size (p, dim = 2), 1, -1
          select case (nstates (j))
          case (3) ! massive vectors
             s(j) = modulo (s(j) + 2, 3) - 1
          case (2) ! spinors, massless vectors
             s(j) = - s(j)
          case (1) ! scalars
             s(j) = -1
          case (0) ! fized spin
             s(j) = -1
          case default ! ???
             s(j) = -1
          end select
          if (s(j) /= -1) then
             cycle sum_spins
          end if
       end do
       exit sum_spins
    end do sum_spins
    sigma = sigma / num_states (2, nstates(1:2))
  end subroutine omega_sum_nonzero
  pure function num_states (n, states) result (ns)
    integer, intent(in) :: n
    integer, dimension(:), intent(in), optional :: states
    integer :: ns
    if (present (states)) then
       ns = product (states, mask = states == 2 .or. states == 3)
    else
       ns = 2**n
    end if
  end function num_states
  pure subroutine omega_nonzero (a, omega, p, s, zero, n, states)
    complex(kind=omega_prec), intent(out) :: a
    real(kind=omega_prec), dimension(0:,:), intent(in) :: p
    integer, dimension(:), intent(in) :: s
    integer, dimension(:), intent(inout) :: zero
    integer, intent(in) :: n
    integer, dimension(:), intent(in), optional :: states
    interface
       pure function omega (p, s) result (me)
         use omega_kinds
         implicit none
         complex(kind=omega_prec) :: me
         real(kind=omega_prec), dimension(0:,:), intent(in) :: p
         integer, dimension(:), intent(in) :: s
       end function omega
    end interface
    real(kind=omega_prec) :: dummy
    if (n < SAMPLE) then
       call omega_sum_nonzero (dummy, omega, p, zero, n, states)
    end if
    if (zero (state_index (s, states)) < REPEAT) then
       a = omega (p, s)
    else
       a = 0
    end if
  end subroutine omega_nonzero
end module omega_utils
