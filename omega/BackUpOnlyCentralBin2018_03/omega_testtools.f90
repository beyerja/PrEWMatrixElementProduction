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
module omega_testtools
  use omega_kinds
  implicit none
  private
  public :: print_matrix
  public :: expect
  real(kind=omega_prec), parameter, private :: TOLERANCE = 1.0e8
  interface expect
     module procedure expect_integer, expect_real, expect_complex, &
          expect_double_integer, expect_complex_integer, expect_complex_real
  end interface
  private :: expect_integer, expect_real, expect_complex, &
       expect_double_integer, expect_complex_integer, expect_complex_real
contains
  subroutine print_matrix (a)
    complex(kind=omega_prec), dimension(:,:), intent(in) :: a
    integer :: row
    do row = 1, size (a, dim=1)
       write (unit = *, fmt = "(10(tr2, f5.2, '+', f5.2, 'I'))") a(row,:)
    end do
  end subroutine print_matrix
  subroutine expect_integer (x, x0, msg)
    integer, intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    if (x == x0) then
       print *, msg, " passed"
    else
       print *, msg, " FAILED: expected ", x0, " got ", x
    end if
  end subroutine expect_integer
  subroutine expect_real (x, x0, msg)
    real(kind=omega_prec), intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    if (x == x0) then
       print *, msg, " passed exactly"
    else if (abs (x - x0) <= epsilon (x)) then
       print *, msg, " passed at machine precision"
    else if (abs (x - x0) <= TOLERANCE * epsilon (x)) then
       print *, msg, " passed at", &
            ceiling (abs (x - x0) / epsilon (x)), "* machine precision"
    else 
       print *, msg, " FAILED: expected ", x0, " got ", x, " (", &
            (x - x0) / epsilon (x), " epsilon)"
    end if
  end subroutine expect_real
  subroutine expect_complex (x, x0, msg)
    complex(kind=omega_prec), intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    if (x == x0) then
       print *, msg, " passed exactly"
    else if (abs (x - x0) <= epsilon (real(x))) then
       print *, msg, " passed at machine precision"
    else if (abs (x - x0) <= TOLERANCE * epsilon (real(x))) then
       print *, msg, " passed at", &
            ceiling (abs (x - x0) / epsilon (real(x))), "* machine precision"
    else 
       print *, msg, " FAILED: expected ", x0, " got ", x, " (", &
            (x - x0) / epsilon (real(x)), " epsilon)"
    end if
  end subroutine expect_complex
  subroutine expect_double_integer (x, x0, msg)
    real(kind=omega_prec), intent(in) :: x
    integer, intent(in) :: x0
    character(len=*), intent(in) :: msg
    call expect_real (x, real (x0, kind=omega_prec), msg)
  end subroutine expect_double_integer
  subroutine expect_complex_integer (x, x0, msg)
    complex(kind=omega_prec), intent(in) :: x
    integer, intent(in) :: x0
    character(len=*), intent(in) :: msg
    call expect_complex (x, cmplx (x0, kind=omega_prec), msg)
  end subroutine expect_complex_integer
  subroutine expect_complex_real (x, x0, msg)
    complex(kind=omega_prec), intent(in) :: x
    real(kind=omega_prec), intent(in) :: x0
    character(len=*), intent(in) :: msg
    call expect_complex (x, cmplx (x0, kind=omega_prec), msg)
  end subroutine expect_complex_real
end module omega_testtools
