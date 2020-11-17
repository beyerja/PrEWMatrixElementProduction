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
module omega_tensor_polarizations
  use omega_kinds
  use omega_constants
  use omega_vectors
  use omega_tensors
  use omega_polarizations
  implicit none
  private
  public :: eps2
  integer, parameter, public :: omega_tensor_pols_2003_03_A = 0
contains
  pure function eps2 (m, k, s) result (t)
    type(tensor) :: t
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k
    integer, intent(in) :: s
    type(vector) :: ep, em, e0
    t%t = 0
    select case (s)
    case (2)
       ep = eps (m, k, 1)
       t = ep.tprod.ep
    case (1)
       ep = eps (m, k, 1)
       e0 = eps (m, k, 0)
       t = (1 / sqrt (2.0_omega_prec)) &
            * ((ep.tprod.e0) + (e0.tprod.ep))
    case (0)
       ep = eps (m, k, 1)
       e0 = eps (m, k, 0)
       em = eps (m, k, -1)
       t = (1 / sqrt (6.0_omega_prec)) &
             * ((ep.tprod.em) + (em.tprod.ep) - 2*(e0.tprod.e0))
    case (-1)
       e0 = eps (m, k, 0)
       em = eps (m, k, -1)
       t = (1 / sqrt (2.0_omega_prec)) &
             * ((em.tprod.e0) + (e0.tprod.em))
    case (-2)
       em = eps (m, k, -1)
       t = em.tprod.em
    end select
  end function eps2
end module omega_tensor_polarizations
