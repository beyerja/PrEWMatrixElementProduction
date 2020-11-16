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
module omega_bispinors
  use omega_kinds
  use omega_constants
  implicit none
  private
  public :: operator (*), operator (+), operator (-)
  public :: abs
  type, public :: bispinor
     ! private (omegalib needs access, but DON'T TOUCH IT!)
     complex(kind=omega_prec), dimension(4) :: a
  end type bispinor 
  interface operator (*)
    module procedure spinor_product
  end interface
  private :: spinor_product
  interface operator (*)
     module procedure integer_bispinor, bispinor_integer, &
            real_bispinor, double_bispinor, &
            complex_bispinor, dcomplex_bispinor, &
            bispinor_real, bispinor_double, &
            bispinor_complex, bispinor_dcomplex 
  end interface
  private :: integer_bispinor, bispinor_integer, real_bispinor, &
       double_bispinor, complex_bispinor, dcomplex_bispinor, &
       bispinor_real, bispinor_double, bispinor_complex, bispinor_dcomplex
  interface operator (+)
     module procedure plus_bispinor
  end interface
  private :: plus_bispinor
  interface operator (-)
     module procedure neg_bispinor
  end interface
  private :: neg_bispinor
  interface operator (+)
     module procedure add_bispinor
  end interface
  private :: add_bispinor
  interface operator (-)
     module procedure sub_bispinor
  end interface
  private :: sub_bispinor
  interface abs
     module procedure abs_bispinor
  end interface
  private :: abs_bispinor
  integer, parameter, public :: omega_bispinors_2003_03_A = 0
contains
  pure function spinor_product (psil, psir) result (psilpsir)
    complex(kind=omega_prec) :: psilpsir
    type(bispinor), intent(in) :: psil, psir
    type(bispinor) :: psidum
    psidum%a(1) = psir%a(2)
    psidum%a(2) = - psir%a(1)
    psidum%a(3) = - psir%a(4)
    psidum%a(4) = psir%a(3)
    psilpsir = dot_product (conjg (psil%a), psidum%a)    
  end function spinor_product
  pure function integer_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    integer, intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function integer_bispinor
  pure function real_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    real(kind=single), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function real_bispinor
  pure function double_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    real(kind=omega_prec), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function double_bispinor      
  pure function complex_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    complex(kind=single), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function complex_bispinor
  pure function dcomplex_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    complex(kind=omega_prec), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function dcomplex_bispinor   
  pure function bispinor_integer (y, x) result (xy)
    type(bispinor) :: xy
    integer, intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function bispinor_integer
  pure function bispinor_real (y, x) result (xy)
    type(bispinor) :: xy
    real(kind=single), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function bispinor_real
  pure function bispinor_double (y, x) result (xy)
    type(bispinor) :: xy
    real(kind=omega_prec), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function bispinor_double         
  pure function bispinor_complex (y, x) result (xy)
    type(bispinor) :: xy
    complex(kind=single), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function bispinor_complex
  pure function bispinor_dcomplex (y, x) result (xy)
    type(bispinor) :: xy
    complex(kind=omega_prec), intent(in) :: x
    type(bispinor), intent(in) :: y
    xy%a = x * y%a
  end function bispinor_dcomplex         
  pure function plus_bispinor (x) result (plus_x)
    type(bispinor) :: plus_x
    type(bispinor), intent(in) :: x
    plus_x%a = x%a
  end function plus_bispinor
  pure function neg_bispinor (x) result (neg_x)
    type(bispinor) :: neg_x
    type(bispinor), intent(in) :: x
    neg_x%a = - x%a
  end function neg_bispinor
  pure function add_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    type(bispinor), intent(in) :: x, y
    xy%a = x%a + y%a
  end function add_bispinor
  pure function sub_bispinor (x, y) result (xy)
    type(bispinor) :: xy
    type(bispinor), intent(in) :: x, y
    xy%a = x%a - y%a
  end function sub_bispinor
  pure function abs_bispinor (psi) result (x)
    real(kind=omega_prec) :: x
    type(bispinor), intent(in) :: psi
    x = sqrt (dot_product (psi%a, psi%a))
  end function abs_bispinor
end module omega_bispinors                                      
