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
module omega_spinor_couplings
  use omega_kinds
  use omega_constants
  use omega_spinors
  use omega_vectors
  use omega_tensors
  use omega_couplings
  implicit none
  private
  public :: u, ubar, v, vbar
  private :: chi_plus, chi_minus
  public :: brs_u, brs_ubar, brs_v, brs_vbar
  public :: va_ff, v_ff, a_ff, vl_ff, vr_ff, vlr_ff, grav_ff
  public :: f_vaf, f_vf, f_af, f_vlf, f_vrf, f_vlrf
  public :: f_fva, f_fv, f_fa, f_fvl, f_fvr, f_fvlr
  public :: sp_ff, s_ff, p_ff, sl_ff, sr_ff, slr_ff
  public :: f_spf, f_sf, f_pf, f_slf, f_srf, f_slrf
  public :: f_fsp, f_fs, f_fp, f_fsl, f_fsr, f_fslr
  public :: f_gravf, f_fgrav
  public :: pr_psi, pr_psibar
  public :: pj_psi, pj_psibar
  public :: pg_psi, pg_psibar
  integer, parameter, public :: omega_spinor_cpls_2003_03_A = 0
contains
  pure function chi_plus (p) result (chi)
    complex(kind=omega_prec), dimension(2) :: chi
    type(momentum), intent(in) :: p
    real(kind=omega_prec) :: pabs
    pabs = sqrt (dot_product (p%x, p%x))
    if (pabs + p%x(3) <= 1000 * epsilon (pabs) * pabs) then
  !!! OLD VERSION !!!!!!
  !!!  if (1 + p%x(3) / pabs <= epsilon (pabs)) then
  !!!!!!!!!!!!!!!!!!!!!!
       chi = (/ cmplx ( 0.0, 0.0, kind=omega_prec), &
                cmplx ( 1.0, 0.0, kind=omega_prec) /)
    else
       chi = 1 / sqrt (2*pabs*(pabs + p%x(3))) &
            * (/ cmplx (pabs + p%x(3), kind=omega_prec), &
                 cmplx (p%x(1), p%x(2), kind=omega_prec) /)
    end if
  end function chi_plus
  pure function chi_minus (p) result (chi)
    complex(kind=omega_prec), dimension(2) :: chi
    type(momentum), intent(in) :: p
    real(kind=omega_prec) :: pabs
    pabs = sqrt (dot_product (p%x, p%x))
    if (pabs + p%x(3) <= 1000 * epsilon (pabs) * pabs) then
  !!! OLD VERSION !!!!!!!!!!!
  !!!  if (1 + p%x(3) / pabs <= epsilon (pabs)) then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
       chi = (/ cmplx (-1.0, 0.0, kind=omega_prec), &
                cmplx ( 0.0, 0.0, kind=omega_prec) /)
    else
       chi = 1 / sqrt (2*pabs*(pabs + p%x(3))) &
            * (/ cmplx (-p%x(1), p%x(2), kind=omega_prec), &
                 cmplx (pabs + p%x(3), kind=omega_prec) /)
    end if
  end function chi_minus
  pure function u (m, p, s) result (psi)
    type(spinor) :: psi
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    complex(kind=omega_prec), dimension(2) :: chi
    real(kind=omega_prec) :: pabs
    pabs = sqrt (dot_product (p%x, p%x))
    select case (s)
    case (1)
       chi = chi_plus (p)
       psi%a(1:2) = sqrt (max (p%t - pabs, 0.0_omega_prec)) * chi
       psi%a(3:4) = sqrt (p%t + pabs) * chi
    case (-1)
       chi = chi_minus (p)
       psi%a(1:2) = sqrt (p%t + pabs) * chi
       psi%a(3:4) = sqrt (max (p%t - pabs, 0.0_omega_prec)) * chi
    case default
       pabs = m ! make the compiler happy and use m
       psi%a = 0
    end select
  end function u
  pure function ubar (m, p, s) result (psibar)
    type(conjspinor) :: psibar
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    type(spinor) :: psi
    psi = u (m, p, s)
    psibar%a(1:2) = conjg (psi%a(3:4))
    psibar%a(3:4) = conjg (psi%a(1:2))
  end function ubar
  pure function v (m, p, s) result (psi)
    type(spinor) :: psi
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    complex(kind=omega_prec), dimension(2) :: chi
    real(kind=omega_prec) :: pabs
    pabs = sqrt (dot_product (p%x, p%x))
    select case (s)
    case (1)
       chi = chi_minus (p)
       psi%a(1:2) = - sqrt (p%t + pabs) * chi
       psi%a(3:4) =   sqrt (max (p%t - pabs, 0.0_omega_prec)) * chi
    case (-1)
       chi = chi_plus (p)
       psi%a(1:2) =   sqrt (max (p%t - pabs, 0.0_omega_prec)) * chi
       psi%a(3:4) = - sqrt (p%t + pabs) * chi
    case default
       pabs = m ! make the compiler happy and use m
       psi%a = 0
    end select
  end function v
  pure function vbar (m, p, s) result (psibar)
    type(conjspinor) :: psibar
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    type(spinor) :: psi
    psi = v (m, p, s)
    psibar%a(1:2) = conjg (psi%a(3:4))
    psibar%a(3:4) = conjg (psi%a(1:2))
  end function vbar
  pure function brs_u (m, p, s) result (dpsi)
      type(spinor) :: dpsi,psi
      real(kind=omega_prec), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) :: s
      type (vector)::vp
      complex(kind=omega_prec), parameter :: one = (1, 0)
      vp=p
      psi=u(m,p,s)
      dpsi=cmplx(0.0,-1.0)*(f_vf(one,vp,psi)-m*psi)
  end function brs_u
  pure function brs_v (m, p, s) result (dpsi)
      type(spinor) :: dpsi, psi
      real(kind=omega_prec), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) ::   s
      type (vector)::vp
      complex(kind=omega_prec), parameter :: one = (1, 0)
      vp=p
      psi=v(m,p,s)
      dpsi=cmplx(0.0,1.0)*(f_vf(one,vp,psi)+m*psi)
  end function brs_v
   pure function brs_ubar (m, p, s)result (dpsibar)
      type(conjspinor) :: dpsibar, psibar
      real(kind=omega_prec), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) :: s
      type (vector)::vp
      complex(kind=omega_prec), parameter :: one = (1, 0)
       vp=p
       psibar=ubar(m,p,s)
      dpsibar=cmplx(0.0,-1.0)*(f_fv(one,psibar,vp)-m*psibar)
    end function brs_ubar
   pure function brs_vbar (m, p, s) result (dpsibar)
      type(conjspinor) :: dpsibar,psibar
      real(kind=omega_prec), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) :: s
      type(vector)::vp
      complex(kind=omega_prec), parameter :: one = (1, 0)
      vp=p
      psibar=vbar(m,p,s)
     dpsibar=cmplx(0.0,1.0)*(f_fv(one,psibar,vp)+m*psibar)
  end function brs_vbar
  pure function va_ff (gv, ga, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: gv, ga
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gl, gr
    complex(kind=omega_prec) :: g13, g14, g23, g24, g31, g32, g41, g42
    gl = gv + ga
    gr = gv - ga
    g13 = psibar%a(1)*psi%a(3)
    g14 = psibar%a(1)*psi%a(4)
    g23 = psibar%a(2)*psi%a(3)
    g24 = psibar%a(2)*psi%a(4)
    g31 = psibar%a(3)*psi%a(1)
    g32 = psibar%a(3)*psi%a(2)
    g41 = psibar%a(4)*psi%a(1)
    g42 = psibar%a(4)*psi%a(2)
    j%t    =  gr * (   g13 + g24) + gl * (   g31 + g42)
    j%x(1) =  gr * (   g14 + g23) - gl * (   g32 + g41)
    j%x(2) = (gr * ( - g14 + g23) + gl * (   g32 - g41)) * (0, 1)
    j%x(3) =  gr * (   g13 - g24) + gl * ( - g31 + g42)
  end function va_ff
  pure function v_ff (gv, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: gv
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: g13, g14, g23, g24, g31, g32, g41, g42
    g13 = psibar%a(1)*psi%a(3)
    g14 = psibar%a(1)*psi%a(4)
    g23 = psibar%a(2)*psi%a(3)
    g24 = psibar%a(2)*psi%a(4)
    g31 = psibar%a(3)*psi%a(1)
    g32 = psibar%a(3)*psi%a(2)
    g41 = psibar%a(4)*psi%a(1)
    g42 = psibar%a(4)*psi%a(2)
    j%t    =   gv * (   g13 + g24 + g31 + g42)
    j%x(1) =   gv * (   g14 + g23 - g32 - g41)
    j%x(2) =   gv * ( - g14 + g23 + g32 - g41) * (0, 1)
    j%x(3) =   gv * (   g13 - g24 - g31 + g42)
  end function v_ff
  pure function a_ff (ga, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: ga
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: g13, g14, g23, g24, g31, g32, g41, g42
    g13 = psibar%a(1)*psi%a(3)
    g14 = psibar%a(1)*psi%a(4)
    g23 = psibar%a(2)*psi%a(3)
    g24 = psibar%a(2)*psi%a(4)
    g31 = psibar%a(3)*psi%a(1)
    g32 = psibar%a(3)*psi%a(2)
    g41 = psibar%a(4)*psi%a(1)
    g42 = psibar%a(4)*psi%a(2)
    j%t    =   ga * ( - g13 - g24 + g31 + g42)
    j%x(1) = - ga * (   g14 + g23 + g32 + g41)
    j%x(2) =   ga * (   g14 - g23 + g32 - g41) * (0, 1)
    j%x(3) =   ga * ( - g13 + g24 - g31 + g42)
  end function a_ff
  pure function vl_ff (gl, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: gl
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gl2
    complex(kind=omega_prec) :: g31, g32, g41, g42
    gl2 = 2 * gl
    g31 = psibar%a(3)*psi%a(1)
    g32 = psibar%a(3)*psi%a(2)
    g41 = psibar%a(4)*psi%a(1)
    g42 = psibar%a(4)*psi%a(2)
    j%t    =   gl2 * (   g31 + g42)
    j%x(1) = - gl2 * (   g32 + g41)
    j%x(2) =   gl2 * (   g32 - g41) * (0, 1)
    j%x(3) =   gl2 * ( - g31 + g42)
  end function vl_ff
  pure function vr_ff (gr, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: gr
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gr2
    complex(kind=omega_prec) :: g13, g14, g23, g24
    gr2 = 2 * gr
    g13 = psibar%a(1)*psi%a(3)
    g14 = psibar%a(1)*psi%a(4)
    g23 = psibar%a(2)*psi%a(3)
    g24 = psibar%a(2)*psi%a(4)
    j%t    = gr2 * (   g13 + g24)
    j%x(1) = gr2 * (   g14 + g23)
    j%x(2) = gr2 * ( - g14 + g23) * (0, 1)
    j%x(3) = gr2 * (   g13 - g24)
  end function vr_ff
  pure function grav_ff (g, m, kb, k, psibar, psi) result (j)
    type(tensor) :: j
    complex(kind=omega_prec), intent(in) :: g
    real(kind=omega_prec), intent(in) :: m
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    type(momentum), intent(in) :: kb, k  
    complex(kind=omega_prec) :: g2, g8, c_dum
    type(vector) :: v_dum
    type(tensor) :: t_metric
    t_metric%t = 0
    t_metric%t(0,0) = 1.0_omega_prec
    t_metric%t(1,1) = - 1.0_omega_prec
    t_metric%t(2,2) = - 1.0_omega_prec
    t_metric%t(3,3) = - 1.0_omega_prec
    g2 = g/2.0_omega_prec
    g8 = g/8.0_omega_prec
    v_dum = v_ff(g8, psibar, psi)
    c_dum = (- m) * s_ff (g2, psibar, psi) - (kb+k)*v_dum
    j = c_dum*t_metric - (((kb+k).tprod.v_dum) + &
         (v_dum.tprod.(kb+k)))
  end function grav_ff
  pure function vlr_ff (gl, gr, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=omega_prec), intent(in) :: gl, gr
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j = va_ff (gl+gr, gl-gr, psibar, psi)
  end function vlr_ff
  pure function f_vaf (gv, ga, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: gv, ga
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gl, gr
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gl = gv + ga
    gr = gv - ga
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gr * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gr * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = gl * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gl * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vaf
  pure function f_vf (gv, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: gv
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gv * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gv * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = gv * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gv * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vf
  pure function f_af (ga, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: ga
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = ga * ( - vm  * psi%a(3) + v12s * psi%a(4))
    vpsi%a(2) = ga * (   v12 * psi%a(3) - vp   * psi%a(4))
    vpsi%a(3) = ga * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = ga * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_af
  pure function f_vlf (gl, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: gl
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gl2
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gl2 = 2 * gl
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = 0
    vpsi%a(2) = 0
    vpsi%a(3) = gl2 * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gl2 * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vlf
  pure function f_vrf (gr, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: gr
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    complex(kind=omega_prec) :: gr2
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gr2 = 2 * gr
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gr2 * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gr2 * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = 0
    vpsi%a(4) = 0
  end function f_vrf
  pure function f_vlrf (gl, gr, v, psi) result (vpsi)
    type(spinor) :: vpsi
    complex(kind=omega_prec), intent(in) :: gl, gr
    type(vector), intent(in) :: v
    type(spinor), intent(in) :: psi
    vpsi = f_vaf (gl+gr, gl-gr, v, psi)
  end function f_vlrf
  pure function f_fva (gv, ga, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: gv, ga
    type(conjspinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    complex(kind=omega_prec) :: gl, gr
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gl = gv + ga
    gr = gv - ga
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    psibarv%a(1) = gl * (   psibar%a(3) * vp   + psibar%a(4) * v12)
    psibarv%a(2) = gl * (   psibar%a(3) * v12s + psibar%a(4) * vm )
    psibarv%a(3) = gr * (   psibar%a(1) * vm   - psibar%a(2) * v12)
    psibarv%a(4) = gr * ( - psibar%a(1) * v12s + psibar%a(2) * vp )
  end function f_fva
  pure function f_fv (gv, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: gv
    type(conjspinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    psibarv%a(1) = gv * (   psibar%a(3) * vp   + psibar%a(4) * v12)
    psibarv%a(2) = gv * (   psibar%a(3) * v12s + psibar%a(4) * vm )
    psibarv%a(3) = gv * (   psibar%a(1) * vm   - psibar%a(2) * v12)
    psibarv%a(4) = gv * ( - psibar%a(1) * v12s + psibar%a(2) * vp )
  end function f_fv
  pure function f_fa (ga, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: ga
    type(vector), intent(in) :: v
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    psibarv%a(1) = ga * (   psibar%a(3) * vp   + psibar%a(4) * v12)
    psibarv%a(2) = ga * (   psibar%a(3) * v12s + psibar%a(4) * vm )
    psibarv%a(3) = ga * ( - psibar%a(1) * vm   + psibar%a(2) * v12)
    psibarv%a(4) = ga * (   psibar%a(1) * v12s - psibar%a(2) * vp )
  end function f_fa
  pure function f_fvl (gl, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: gl
    type(conjspinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    complex(kind=omega_prec) :: gl2
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gl2 = 2 * gl
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    psibarv%a(1) = gl2 * (   psibar%a(3) * vp   + psibar%a(4) * v12)
    psibarv%a(2) = gl2 * (   psibar%a(3) * v12s + psibar%a(4) * vm )
    psibarv%a(3) = 0
    psibarv%a(4) = 0
  end function f_fvl
  pure function f_fvr (gr, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: gr
    type(conjspinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    complex(kind=omega_prec) :: gr2
    complex(kind=omega_prec) :: vp, vm, v12, v12s
    gr2 = 2 * gr
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    psibarv%a(1) = 0
    psibarv%a(2) = 0
    psibarv%a(3) = gr2 * (   psibar%a(1) * vm   - psibar%a(2) * v12)
    psibarv%a(4) = gr2 * ( - psibar%a(1) * v12s + psibar%a(2) * vp )
  end function f_fvr
  pure function f_fvlr (gl, gr, psibar, v) result (psibarv)
    type(conjspinor) :: psibarv
    complex(kind=omega_prec), intent(in) :: gl, gr
    type(conjspinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    psibarv = f_fva (gl+gr, gl-gr, psibar, v)
  end function f_fvlr
  pure function sp_ff (gs, gp, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gs, gp
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j =    (gs - gp) * (psibar%a(1)*psi%a(1) + psibar%a(2)*psi%a(2)) &
         + (gs + gp) * (psibar%a(3)*psi%a(3) + psibar%a(4)*psi%a(4))
  end function sp_ff
  pure function s_ff (gs, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gs
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j = gs * (psibar * psi)
  end function s_ff
  pure function p_ff (gp, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gp
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j = gp * (  psibar%a(3)*psi%a(3) + psibar%a(4)*psi%a(4) &
              - psibar%a(1)*psi%a(1) - psibar%a(2)*psi%a(2))
  end function p_ff
  pure function sl_ff (gl, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gl
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j =  2 * gl * (psibar%a(1)*psi%a(1) + psibar%a(2)*psi%a(2))
  end function sl_ff
  pure function sr_ff (gr, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gr
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j = 2 * gr * (psibar%a(3)*psi%a(3) + psibar%a(4)*psi%a(4))
  end function sr_ff
  pure function slr_ff (gl, gr, psibar, psi) result (j)
    complex(kind=omega_prec) :: j
    complex(kind=omega_prec), intent(in) :: gl, gr
    type(conjspinor), intent(in) :: psibar
    type(spinor), intent(in) :: psi
    j = sp_ff (gr+gl, gr-gl, psibar, psi)
  end function slr_ff
  pure function f_spf (gs, gp, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gs, gp
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi%a(1:2) = ((gs - gp) * phi) * psi%a(1:2)
    phipsi%a(3:4) = ((gs + gp) * phi) * psi%a(3:4)
  end function f_spf
  pure function f_sf (gs, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gs
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi%a = (gs * phi) * psi%a
  end function f_sf
  pure function f_pf (gp, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gp
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi%a(1:2) = (- gp * phi) * psi%a(1:2)
    phipsi%a(3:4) = (  gp * phi) * psi%a(3:4)
  end function f_pf
  pure function f_slf (gl, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gl
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi%a(1:2) = (2 * gl * phi) * psi%a(1:2)
    phipsi%a(3:4) = 0
  end function f_slf
  pure function f_srf (gr, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gr
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi%a(1:2) = 0
    phipsi%a(3:4) = (2 * gr * phi) * psi%a(3:4)
  end function f_srf
  pure function f_slrf (gl, gr, phi, psi) result (phipsi)
    type(spinor) :: phipsi
    complex(kind=omega_prec), intent(in) :: gl, gr
    complex(kind=omega_prec), intent(in) :: phi
    type(spinor), intent(in) :: psi
    phipsi =  f_spf (gr+gl, gr-gl, phi, psi)
  end function f_slrf
  pure function f_fsp (gs, gp, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gs, gp
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi%a(1:2) = ((gs - gp) * phi) * psibar%a(1:2)
    psibarphi%a(3:4) = ((gs + gp) * phi) * psibar%a(3:4)
  end function f_fsp
  pure function f_fs (gs, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gs
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi%a = (gs * phi) * psibar%a
  end function f_fs
  pure function f_fp (gp, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gp
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi%a(1:2) = (- gp * phi) * psibar%a(1:2)
    psibarphi%a(3:4) = (  gp * phi) * psibar%a(3:4)
  end function f_fp
  pure function f_fsl (gl, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gl
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi%a(1:2) = (2 * gl * phi) * psibar%a(1:2)
    psibarphi%a(3:4) = 0
  end function f_fsl
  pure function f_fsr (gr, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gr
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi%a(1:2) = 0
    psibarphi%a(3:4) = (2 * gr * phi) * psibar%a(3:4)
  end function f_fsr
  pure function f_fslr (gl, gr, psibar, phi) result (psibarphi)
    type(conjspinor) :: psibarphi
    complex(kind=omega_prec), intent(in) :: gl, gr
    type(conjspinor), intent(in) :: psibar
    complex(kind=omega_prec), intent(in) :: phi
    psibarphi = f_fsp (gr+gl, gr-gl, psibar, phi)
  end function f_fslr
  pure function f_gravf (g, m, kb, k, t, psi) result (tpsi)
    type(spinor) :: tpsi
    complex(kind=omega_prec), intent(in) :: g
    real(kind=omega_prec), intent(in) :: m
    type(spinor), intent(in) :: psi
    type(tensor), intent(in) :: t
    type(momentum), intent(in) :: kb, k
    complex(kind=omega_prec) :: g2, g8, t_tr
    type(vector) :: kkb
    kkb = k + kb
    g2 = g / 2.0_omega_prec
    g8 = g / 8.0_omega_prec
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)      
    tpsi = (- f_sf (g2, cmplx (m,0.0, kind=omega_prec), psi) & 
           - f_vf ((g8*m), kkb, psi)) * t_tr - &
           f_vf (g8,(t*kkb + kkb*t),psi)
  end function f_gravf
  pure function f_fgrav (g, m, kb, k, psibar, t) result (psibart)
    type(conjspinor) :: psibart
    complex(kind=omega_prec), intent(in) :: g
    real(kind=omega_prec), intent(in) :: m
    type(conjspinor), intent(in) :: psibar
    type(tensor), intent(in) :: t
    type(momentum), intent(in) :: kb, k
    type(vector) :: kkb
    complex(kind=omega_prec) :: g2, g8, t_tr
    kkb = k + kb
    g2 = g / 2.0_omega_prec
    g8 = g / 8.0_omega_prec
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)      
    psibart = (- f_fs (g2, psibar, cmplx (m, 0.0, kind=omega_prec)) &
        - f_fv ((g8 * m), psibar, kkb)) * t_tr - &
        f_fv (g8,psibar,(t*kkb + kkb*t)) 
  end function f_fgrav
  pure function pr_psi (p, m, w, psi) result (ppsi)
    type(spinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(spinor), intent(in) :: psi
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsi = (1 / cmplx (p*p - m**2, m*w, kind=omega_prec)) &
         * (- f_vf (one, vp, psi) + m * psi)
  end function pr_psi
  pure function pj_psi (p, m, w, psi) result (ppsi)
    type(spinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(spinor), intent(in) :: psi
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsi = (0, -1) * sqrt (PI / m / w) * (- f_vf (one, vp, psi) + m * psi)
  end function pj_psi
  pure function pg_psi (p, m, w, psi) result (ppsi)
    type(spinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(spinor), intent(in) :: psi
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsi = gauss(p*p, m, w) *  (- f_vf (one, vp, psi) + m * psi)
  end function pg_psi
  pure function pr_psibar (p, m, w, psibar) result (ppsibar)
    type(conjspinor) :: ppsibar
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(conjspinor), intent(in) :: psibar
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsibar = (1 / cmplx (p*p - m**2, m*w, kind=omega_prec)) &
         * (f_fv (one, psibar, vp) + m * psibar)
  end function pr_psibar
  pure function pj_psibar (p, m, w, psibar) result (ppsibar)
    type(conjspinor) :: ppsibar
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(conjspinor), intent(in) :: psibar
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsibar = (0, -1) * sqrt (PI / m / w) * (f_fv (one, psibar, vp) + m * psibar)
  end function pj_psibar
  pure function pg_psibar (p, m, w, psibar) result (ppsibar)
    type(conjspinor) :: ppsibar
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(conjspinor), intent(in) :: psibar
    type(vector) :: vp
    complex(kind=omega_prec), parameter :: one = (1, 0)
    vp = p
    ppsibar = gauss (p*p, m, w) * (f_fv (one, psibar, vp) + m * psibar)
  end function pg_psibar
end module omega_spinor_couplings
