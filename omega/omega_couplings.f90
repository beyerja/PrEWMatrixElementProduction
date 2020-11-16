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
module omega_couplings
  use omega_kinds
  use omega_constants
  use omega_vectors
  use omega_tensors
  implicit none
  private
  public :: g_gg
  public :: x_gg, g_gx
  public :: v_ss, s_vs
  public :: tkv_vv, lkv_vv, tv_kvv, lv_kvv, kg_kgkg
  public :: t5kv_vv, l5kv_vv, t5v_kvv, l5v_kvv, kg5_kgkg, kg_kg5kg
  public :: s_gravs, v_gravv, grav_ss, grav_vv
  public :: t2_vv, v_t2v
  public :: phi_vv, v_phiv
  public :: t2_vv_d5_1, v_t2v_d5_1
  public :: t2_vv_d5_2, v_t2v_d5_2
  public :: t2_vv_d7, v_t2v_d7
  public :: wd_tl
  public :: gauss
  public :: pr_phi, pr_unitarity, pr_feynman, pr_gauge, pr_rxi
  public :: pj_phi, pj_unitarity
  public :: pg_phi, pg_unitarity
  public :: pr_tensor
  integer, parameter, public :: omega_couplings_2003_03_A = 0
contains
  pure function g_gg (g, a1, k1, a2, k2) result (a)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    a = (0, -1) * g * ((k1 - k2) * (a1 * a2) &
                        + ((2*k2 + k1) * a1) * a2 - a1 * ((2*k1 + k2) * a2))
  end function g_gg
  pure function x_gg (g, a1, a2) result (x)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(tensor2odd) :: x
    x = g * (a1 .wedge. a2)
  end function x_gg
  pure function g_gx (g, a1, x) result (a)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1
    type(tensor2odd), intent(in) :: x
    type(vector) :: a
    a = g * (a1 * x)
  end function g_gx
  pure function v_ss (g, phi1, k1, phi2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 - k2) * (g * phi1 * phi2)
  end function v_ss
  pure function s_vs (g, v1, k1, phi2, k2) result (phi)
    complex(kind=omega_prec), intent(in) :: g, phi2
    type(vector), intent(in) :: v1
    type(momentum), intent(in) :: k1, k2
    complex(kind=omega_prec) :: phi
    phi = g * ((k1 + 2*k2) * v1) * phi2
  end function s_vs
  pure function tkv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 - k2) * ((0, 1) * g * (v1*v2))
  end function tkv_vv
  pure function t5kv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 - k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function t5kv_vv
  pure function lkv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 + k2) * ((0, 1) * g * (v1*v2))
  end function lkv_vv
  pure function l5kv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 + k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function l5kv_vv
  pure function tv_kvv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = v2 * ((0, 1) * g * ((2*k2 + k1)*v1))
  end function tv_kvv
  pure function t5v_kvv (g, v1, k1, v2, k2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 + 2*k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function t5v_kvv
  pure function lv_kvv (g, v1, k1, v2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1
    type(vector) :: v
    v = v2 * ((0, -1) * g * (k1*v1))
  end function lv_kvv
  pure function l5v_kvv (g, v1, k1, v2) result (v)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1
    type(vector) :: v
    type(vector) :: k
    k = k1
    v = (0, -1) * g * pseudo_vector (k, v1, v2)
  end function l5v_kvv
  pure function kg_kgkg (g, a1, k1, a2, k2) result (a)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    real(kind=omega_prec) :: k1k1, k2k2, k1k2, kk1, kk2
    complex(kind=omega_prec) :: a1a2, k2a1, ka1, k1a2, ka2
    k1k1 = k1 * k1
    k1k2 = k1 * k2
    k2k2 = k2 * k2
    kk1 = k1k1 + k1k2
    kk2 = k1k2 + k2k2
    k2a1 = k2 * a1
    ka1 = k2a1 + k1 * a1
    k1a2 = k1 * a2
    ka2 = k1a2 + k2 * a2
    a1a2 = a1 * a2
    a = (0, -1) * g * (   (kk2  * k1a2 - k1k2 * ka2 ) * a1 &
                        + (k1k2 * ka1  - kk1  * k2a1) * a2 &
                        + (ka2  * k2a1 - kk2  * a1a2) * k1 &
                        + (kk1  * a1a2 - ka1  * k1a2) * k2 )
  end function kg_kgkg
  pure function kg5_kgkg (g, a1, k1, a2, k2) result (a)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    type(vector) :: kv, k1v, k2v
    kv = - k1 - k2
    k1v = k1
    k2v = k2
    a = (0, -2) * g * (   (k2*A1) * pseudo_vector (kv, k1v, a2 ) &
                        + (k1*A2) * pseudo_vector (kv, A1 , k2v) &
                        - (A1*A2) * pseudo_vector (kv, k1v, k2v) &
                        - (k1*k2) * pseudo_vector (kv, a1 , a2 ) )
  end function kg5_kgkg
  pure function kg_kg5kg (g, a1, k1, a2, k2) result (a)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    type(vector) :: kv, k1v, k2v
    kv = - k1 - k2
    k1v = k1
    k2v = k2
    a = (0, -1) * g * (   (kv*k2v) * pseudo_vector (a2 , k1v, a1) &
                        - (kv*a2 ) * pseudo_vector (k2v, k1v, a1) &
                        -  k2v * pseudo_scalar (kv, a2,  k1v, a1) &
                        +  a2  * pseudo_scalar (kv, k2v, k1v, a1) )
  end function kg_kg5kg
  pure function s_gravs (g, m, k1, k2, t, s) result (phi)
    complex(kind=omega_prec), intent(in) :: g, s
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(tensor), intent(in) :: t
    complex(kind=omega_prec) :: phi, t_tr
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)      
    phi = g * s * (((t*k1)*k2) + ((t*k2)*k1) & 
          - g * (m**2 + (k1*k2))*t_tr)/2.0_omega_prec
  end function s_gravs
  pure function grav_ss (g, m, k1, k2, s1, s2) result (t)
    complex(kind=omega_prec), intent(in) :: g, s1, s2
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t_metric, t 
    t_metric%t = 0
    t_metric%t(0,0) = 1.0_omega_prec
    t_metric%t(1,1) = - 1.0_omega_prec
    t_metric%t(2,2) = - 1.0_omega_prec
    t_metric%t(3,3) = - 1.0_omega_prec
    t = g*s1*s2/2.0_omega_prec * (-(m**2 + (k1*k2)) * t_metric &
         + (k1.tprod.k2) + (k2.tprod.k1))
  end function grav_ss
  pure function v_gravv (g, m, k1, k2, t, v) result (vec)
    complex(kind=omega_prec), intent(in) :: g
    real(kind=omega_prec), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v
    type(tensor), intent(in) :: t
    complex(kind=omega_prec) :: t_tr
    real(kind=omega_prec) :: xi
    type(vector) :: vec
    xi = 1.0_omega_prec
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)      
    vec = (-g)/ 2.0_omega_prec * (((k1*k2) + m**2) * &
         (t*v + v*t - t_tr * v) + t_tr * (k1*v) * k2 &
         - (k1*v) * ((k2*t) + (t*k2)) &
         - ((k1*(t*v)) + (v*(t*k1))) * k2 &
         + ((k1*(t*k2)) + (k2*(t*k1))) * v)
  !!!       Unitarity gauge: xi -> Infinity
  !!!       + (1.0_omega_prec/xi) * (t_tr * ((k1*v)*k2) + &
  !!!       (k2*v)*k2 + (k2*v)*k1 - (k1*(t*v))*k1 + &
  !!!       (k2*v)*(k2*t) - (v*(t*k1))*k1 - (k2*v)*(t*k2)))
  end function v_gravv
  pure function grav_vv (g, m, k1, k2, v1, v2) result (t)
    complex(kind=omega_prec), intent(in) :: g
    type(momentum), intent(in) :: k1, k2
    real(kind=omega_prec), intent(in) :: m
    real(kind=omega_prec) :: xi
    type(vector), intent (in) :: v1, v2
    type(tensor) :: t_metric, t 
    xi = 0.00001_omega_prec
    t_metric%t = 0
    t_metric%t(0,0) = 1.0_omega_prec
    t_metric%t(1,1) = - 1.0_omega_prec
    t_metric%t(2,2) = - 1.0_omega_prec
    t_metric%t(3,3) = - 1.0_omega_prec
    t = (-g)/2.0_omega_prec * ( &
         ((k1*k2) + m**2) * ( &
         (v1.tprod.v2) +  (v2.tprod.v1) - (v1*v2) * t_metric) &
         + (v1*k2)*(v2*k1)*t_metric & 
         - (k2*v1)*((v2.tprod.k1) + (k1.tprod.v2)) &
         - (k1*v2)*((v1.tprod.k2) + (k2.tprod.v1)) &
         + (v1*v2)*((k1.tprod.k2) + (k2.tprod.k1)))
  !!!       Unitarity gauge: xi -> Infinity
  !!!       + (1.0_omega_prec/xi) * ( &
  !!!       ((k1*v1)*(k1*v2) + (k2*v1)*(k2*v2) + (k1*v1)*(k2*v2))* &
  !!!       t_metric) - (k1*v1) * ((k1.tprod.v2) + (v2.tprod.k1)) &
  !!!       - (k2*v2) * ((k2.tprod.v1) + (v1.tprod.k2)))
  end function grav_vv
  pure function phi_vv (g, k1, k2, v1, v2) result (phi)
    complex(kind=omega_prec), intent(in) :: g
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v1, v2
    complex(kind=omega_prec) :: phi
    phi = g * pseudo_scalar (k1, v1, k2, v2)
  end function phi_vv
  pure function v_phiv (g, phi, k1, k2, v) result (w)
    complex(kind=omega_prec), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k1, k2
    type(vector) :: w
    w = g * phi * pseudo_vector (k1, k2, v)
  end function v_phiv
  pure function t2_vv (g, v1, v2) result (t)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(tensor) :: t
    type(tensor) :: tmp
    tmp = v1.tprod.v2
    t%t = g * (tmp%t + transpose (tmp%t))
  end function t2_vv
  pure function v_t2v (g, t, v) result (tv)
    complex(kind=omega_prec), intent(in) :: g
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(vector) :: tv
    type(tensor) :: tmp
    tmp%t = t%t + transpose (t%t)
    tv = g * (tmp * v)
  end function v_t2v
  pure function t2_vv_d5_1 (g, v1, k1, v2, k2) result (t)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (v1 * v2)) * (k1-k2).tprod.(k1-k2)
  end function t2_vv_d5_1
  pure function v_t2v_d5_1 (g, t1, k1, v2, k2) result (tv)
    complex(kind=omega_prec), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    tv = (g * ((k1+2*k2).tprod.(k1+2*k2) * t1)) * v2
  end function v_t2v_d5_1
  pure function t2_vv_d5_2 (g, v1, k1, v2, k2) result (t)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (k2 * v1)) * (k2-k1).tprod.v2
    t%t = t%t + transpose (t%t)
  end function t2_vv_d5_2
  pure function v_t2v_d5_2 (g, t1, k1, v2, k2) result (tv)
    complex(kind=omega_prec), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    type(tensor) :: tmp
    type(momentum) :: k1_k2, k1_2k2
    k1_k2 = k1 + k2
    k1_2k2 = k1_k2 + k2
    tmp%t = t1%t + transpose (t1%t)
    tv = (g * (k1_k2 * v2)) * (k1_2k2 * tmp)
  end function v_t2v_d5_2
  pure function t2_vv_d7 (g, v1, k1, v2, k2) result (t)
    complex(kind=omega_prec), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (k2 * v1) * (k1 * v2)) * (k1-k2).tprod.(k1-k2)
  end function t2_vv_d7
  pure function v_t2v_d7 (g, t1, k1, v2, k2) result (tv)
    complex(kind=omega_prec), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    type(vector) :: k1_k2, k1_2k2
    k1_k2 = k1 + k2
    k1_2k2 = k1_k2 + k2
    tv = (- g * (k1_k2 * v2) * (k1_2k2.tprod.k1_2k2 * t1)) * k2
  end function v_t2v_d7
  pure function wd_tl (p, w) result (width)
    real(kind=omega_prec) :: width
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: w
    if (p*p > 0) then
       width = w
    else
       width = 0
    end if
  end function wd_tl
  pure function gauss (x, mu, w) result (gg)
    real(kind=omega_prec) :: gg
    real(kind=omega_prec), intent(in) :: x, mu, w
    if (w > 0) then
      gg = exp(-(x - mu**2)**2/4.0_omega_prec/mu**2/w**2) * &
           sqrt(sqrt(PI/2)) / w / mu
            else
      gg = 1.0_omega_prec 
    end if
  end function gauss
  pure function pr_phi (p, m, w, phi) result (pphi)
    complex(kind=omega_prec) :: pphi
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    complex(kind=omega_prec), intent(in) :: phi
    pphi = (1 / cmplx (p*p - m**2, m*w, kind=omega_prec)) * phi 
  end function pr_phi
  pure function pj_phi (m, w, phi) result (pphi)
    complex(kind=omega_prec) :: pphi
    real(kind=omega_prec), intent(in) :: m, w
    complex(kind=omega_prec), intent(in) :: phi
    pphi = (0, -1) * sqrt (PI / m / w) * phi 
  end function pj_phi
  pure function pg_phi (p, m, w, phi) result (pphi)
    complex(kind=omega_prec) :: pphi      
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    complex(kind=omega_prec), intent(in) :: phi
    pphi = ((0, 1) * gauss (p*p, m, w)) * phi
  end function pg_phi  
  pure function pr_unitarity (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(vector), intent(in) :: e
    type(vector) :: pv
    pv = p
    pe = - (1 / cmplx (p*p - m**2, m*w, kind=omega_prec)) &
         * (e - (p*e / m**2) * pv)
  end function pr_unitarity
  pure function pj_unitarity (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(vector), intent(in) :: e
    type(vector) :: pv
    pv = p
    pe = (0, 1) * sqrt (PI / m / w) * (e - (p*e / m**2) * pv)
  end function pj_unitarity
  pure function pg_unitarity (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(vector), intent(in) :: e
    type(vector) :: pv
    pv = p
    pe = - gauss (p*p, m, w) &
         * (e - (p*e / m**2) * pv)
  end function pg_unitarity  
  pure function pr_feynman (p, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    type(vector), intent(in) :: e
    pe = - (1 / (p*p)) * e
  end function pr_feynman
  pure function pr_gauge (p, xi, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: xi
    type(vector), intent(in) :: e
    real(kind=omega_prec) :: p2
    type(vector) :: pv
    p2 = p*p
    pv = p
    pe = - (1 / p2) * (e - ((1 - xi) * (p*e) / p2) * pv)
  end function pr_gauge
  pure function pr_rxi (p, m, w, xi, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w, xi
    type(vector), intent(in) :: e
    real(kind=omega_prec) :: p2 
    type(vector) :: pv
    p2 = p*p
    pv = p
    pe = - (1 / cmplx (p2 - m**2, m*w, kind=omega_prec)) &
         * (e - ((1 - xi) * (p*e) / (p2 - xi * m**2)) * pv)
  end function pr_rxi
  pure function pr_tensor (p, m, w, t) result (pt)
    type(tensor) :: pt
    type(momentum), intent(in) :: p
    real(kind=omega_prec), intent(in) :: m, w
    type(tensor), intent(in) :: t
    complex(kind=omega_prec) :: p_dd_t
    real(kind=omega_prec), dimension(0:3,0:3) :: p_uu, p_ud, p_du, p_dd
    integer :: i, j
    p_uu(0,0) = 1 - p%t * p%t / m**2
    p_uu(0,1:3) = - p%t * p%x / m**2
    p_uu(1:3,0) = p_uu(0,1:3)
    do i = 1, 3
       do j = 1, 3
          p_uu(i,j) = - p%x(i) * p%x(j) / m**2
       end do
    end do
    do i = 1, 3
       p_uu(i,i) = - 1 + p_uu(i,i)
    end do
    p_ud(:,0) = p_uu(:,0)
    p_ud(:,1:3) = - p_uu(:,1:3)
    p_du = transpose (p_ud)
    p_dd(:,0) = p_du(:,0)
    p_dd(:,1:3) = - p_du(:,1:3)
    p_dd_t = 0
    do i = 0, 3
       do j = 0, 3
          p_dd_t = p_dd_t + p_dd(i,j) * t%t(i,j)
       end do
    end do
    pt%t = matmul (p_ud, matmul (0.5_omega_prec * (t%t + transpose (t%t)), p_du)) &
         - (p_dd_t / 3.0_omega_prec) * p_uu
    pt%t = pt%t / cmplx (p*p - m**2, m*w, kind=omega_prec)
  end function pr_tensor
end module omega_couplings
