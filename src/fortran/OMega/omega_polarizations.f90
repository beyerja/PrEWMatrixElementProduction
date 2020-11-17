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
module omega_polarizations
  use omega_kinds
  use omega_constants
  use omega_vectors
  implicit none
  private
  public :: eps
  integer, parameter, public :: omega_polarizations_2003_03_A = 0
contains
  pure function eps (m, k, s) result (e)
    type(vector) :: e
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    integer, intent(in) :: s
    real(kind=omega_prec) :: kt, kabs, kabs2, sqrt2
    sqrt2 = sqrt (2.0_omega_prec)
    kabs2 = dot_product (k%x, k%x)
    e%t = 0
    e%x = 0
    if (kabs2 > 0) then
       kabs = sqrt (kabs2)
       select case (s)
       case (1)
          kt = sqrt (k%x(1)**2 + k%x(2)**2)
          if (abs(kt) <= epsilon(kt) * kabs) then 
             if (k%x(3) > 0) then           
                e%x(1) = cmplx (   1,   0, kind=omega_prec) / sqrt2
                e%x(2) = cmplx (   0,   1, kind=omega_prec) / sqrt2
             else
                e%x(1) = cmplx ( - 1,   0, kind=omega_prec) / sqrt2
                e%x(2) = cmplx (   0,   1, kind=omega_prec) / sqrt2
             end if
          else
             e%x(1) = cmplx (   k%x(3)*k%x(1)/kabs, &
                  - k%x(2), kind=omega_prec) / kt / sqrt2
             e%x(2) = cmplx (   k%x(2)*k%x(3)/kabs, &
                  k%x(1), kind=omega_prec) / kt / sqrt2
             e%x(3) = - kt / kabs / sqrt2
          end if
       case (-1)
          kt = sqrt (k%x(1)**2 + k%x(2)**2)
          if (abs(kt) <= epsilon(kt) * kabs) then
             if (k%x(3) > 0) then
                e%x(1) = cmplx (   1,   0, kind=omega_prec) / sqrt2
                e%x(2) = cmplx (   0, - 1, kind=omega_prec) / sqrt2           
             else
                e%x(1) = cmplx (  -1,   0, kind=omega_prec) / sqrt2
                e%x(2) = cmplx (   0, - 1, kind=omega_prec) / sqrt2           
             end if
          else
             e%x(1) = cmplx (   k%x(3)*k%x(1)/kabs, &
                  k%x(2), kind=omega_prec) / kt / sqrt2
             e%x(2) = cmplx (   k%x(2)*k%x(3)/kabs, &
                  - k%x(1), kind=omega_prec) / kt / sqrt2
             e%x(3) = - kt / kabs / sqrt2
          end if
       case (0)
          if (m > 0) then
             e%t = kabs / m
             e%x = k%t / (m*kabs) * k%x 
          end if
       case (3) 
          e = (0,1) * k
       case (4)
          if (m > 0) then
             e = (1 / m) * k
          else
             e = (1 / k%t) * k
          end if
       end select
    else   !!! for particles in their rest frame defined to be 
           !!! polarized along the 3-direction
       select case (s)
       case (1)
          e%x(1) = cmplx (   1,   0, kind=omega_prec) / sqrt2
          e%x(2) = cmplx (   0,   1, kind=omega_prec) / sqrt2
       case (-1)
          e%x(1) = cmplx (   1,   0, kind=omega_prec) / sqrt2
          e%x(2) = cmplx (   0, - 1, kind=omega_prec) / sqrt2
       case (0)
          if (m > 0) then
             e%x(3) = 1
          end if
       case (4)
          if (m > 0) then
             e = (1 / m) * k
          else
             e = (1 / k%t) * k
          end if
       end select
    end if
  end function eps
  !!! OLD VERSION !!!!!!
  !!! pure function eps (m, k, s) result (e)
  !!!   type(vector) :: e
  !!!   real(kind=omega_prec), intent(in) :: m
  !!!   type(momentum), intent(in) :: k
  !!!   integer, intent(in) :: s
  !!!   real(kind=omega_prec) :: kt, kabs, kabs2, sqrt2
  !!!   integer, parameter :: x = 2, y = 3, z = 1
  !!!   sqrt2 = sqrt (2.0_omega_prec)
  !!!   kabs2 = dot_product (k%x, k%x)
  !!!   e%t = 0
  !!!   e%x = 0
  !!!   if (kabs2 > 0) then
  !!!      kabs = sqrt (kabs2)
  !!!      select case (s)
  !!!      case (1)
  !!!         kt = sqrt (k%x(x)**2 + k%x(y)**2)
  !!!         e%x(x) = cmplx (   k%x(z)*k%x(x)/kabs, &
  !!!                          - k%x(y), kind=omega_prec) / kt / sqrt2
  !!!         e%x(y) = cmplx (   k%x(y)*k%x(z)/kabs, &
  !!!                            k%x(x), kind=omega_prec) / kt / sqrt2
  !!!         e%x(z) = - kt / kabs / sqrt2
  !!!      case (-1)
  !!!         kt = sqrt (k%x(x)**2 + k%x(y)**2)
  !!!         e%x(x) = cmplx (   k%x(z)*k%x(x)/kabs, &
  !!!                            k%x(y), kind=omega_prec) / kt / sqrt2
  !!!         e%x(y) = cmplx (   k%x(y)*k%x(z)/kabs, &
  !!!                          - k%x(x), kind=omega_prec) / kt / sqrt2
  !!!         e%x(z) = - kt / kabs / sqrt2
  !!!      case (0)
  !!!         if (m > 0) then
  !!!            e%t = kabs / m
  !!!            e%x = k%t / (m*kabs) * k%x 
  !!!         end if
  !!!      case (3) 
  !!!         e = (0,1) * k
  !!!      case (4)
  !!!         if (m > 0) then
  !!!            e = (1 / m) * k
  !!!         else
  !!!            e = (1 / k%t) * k
  !!!         end if
  !!!      end select
  !!!   else
  !!!      select case (s)
  !!!      case (1)
  !!!         e%x(x) = cmplx (   1,   0, kind=omega_prec) / sqrt2
  !!!         e%x(y) = cmplx (   0,   1, kind=omega_prec) / sqrt2
  !!!      case (-1)
  !!!         e%x(x) = cmplx (   1,   0, kind=omega_prec) / sqrt2
  !!!         e%x(y) = cmplx (   0, - 1, kind=omega_prec) / sqrt2
  !!!      case (0)
  !!!         if (m > 0) then
  !!!            e%x(z) = 1
  !!!         end if
  !!!      case (4)
  !!!         if (m > 0) then
  !!!            e = (1 / m) * k
  !!!         else
  !!!            e = (1 / k%t) * k
  !!!         end if
  !!!      end select
  !!!   end if
  !!! end function eps
  !!!!!!!!!!!!!!!!!!!!!!!!
end module omega_polarizations
