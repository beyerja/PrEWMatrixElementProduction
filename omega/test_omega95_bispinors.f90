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
program test_omega95_bispinors
  use omega_kinds
  use omega95_bispinors
  use omega_vspinor_polarizations
  use omega_testtools
  implicit none
  integer :: i, j
  real(kind=omega_prec) :: m, pabs, qabs, tabs, zabs, w
  real(kind=omega_prec), dimension(4) :: r
  complex(kind=omega_prec) :: one, two
  type(momentum) :: p, q, t, z, p_0
  type(vector) :: vp, vq, vt, vz
  type(vectorspinor) :: testv
  call random_seed ()
  one = 1
  two = 2
  w = 1.4142
  m = 13
  pabs = 42
  qabs = 137
  tabs = 84
  zabs = 3.1415
  p_0%t = m
  p_0%x = 0
  call random_momentum (p, pabs, m)
  call random_momentum (q, qabs, m)
  call random_momentum (t, tabs, m)
  call random_momentum (z, zabs, m)
  call random_number (r)
  do i = 1, 4
     testv%psi(1)%a(i) = (0.0_omega_prec, 0.0_omega_prec)
  end do
  do i = 2, 3
     do j = 1, 4
        testv%psi(i)%a(j) = cmplx (10.0_omega_prec * r(j))
    end do
  end do
  testv%psi(4)%a(1) = (1.0_omega_prec, 0.0_omega_prec)
  testv%psi(4)%a(1) = (0.0_omega_prec, 2.0_omega_prec)
  testv%psi(4)%a(1) = (1.0_omega_prec, 0.0_omega_prec)
  testv%psi(4)%a(1) = (3.0_omega_prec, 0.0_omega_prec)
  vp = p
  vq = q
  vt = t
  vz = z
print *, "*** Checking the equations of motion ***:"
call expect (abs(f_vf(one,vp,u(m,p,+1))-m*u(m,p,+1)), 0, "|[p-m]u(+)|=0")
call expect (abs(f_vf(one,vp,u(m,p,-1))-m*u(m,p,-1)), 0, "|[p-m]u(-)|=0")
call expect (abs(f_vf(one,vp,v(m,p,+1))+m*v(m,p,+1)), 0, "|[p+m]v(+)|=0")
call expect (abs(f_vf(one,vp,v(m,p,-1))+m*v(m,p,-1)), 0, "|[p+m]v(-)|=0")
print *, "*** Checking the normalization ***:"
call expect (s_ff(one,v(m,p,+1),u(m,p,+1)), +2*m, "ubar(+)*u(+)=+2m")
call expect (s_ff(one,v(m,p,-1),u(m,p,-1)), +2*m, "ubar(-)*u(-)=+2m")
call expect (s_ff(one,u(m,p,+1),v(m,p,+1)), -2*m, "vbar(+)*v(+)=-2m")
call expect (s_ff(one,u(m,p,-1),v(m,p,-1)), -2*m, "vbar(-)*v(-)=-2m")
call expect (s_ff(one,v(m,p,+1),v(m,p,+1)),    0, "ubar(+)*v(+)=0  ")
call expect (s_ff(one,v(m,p,-1),v(m,p,-1)),    0, "ubar(-)*v(-)=0  ")
call expect (s_ff(one,u(m,p,+1),u(m,p,+1)),    0, "vbar(+)*u(+)=0  ")
call expect (s_ff(one,u(m,p,-1),u(m,p,-1)),    0, "vbar(-)*u(-)=0  ")
print *, "*** Checking the currents ***:"
call expect (abs(v_ff(one,v(m,p,+1),u(m,p,+1))-2*vp), 0, "ubar(+).V.u(+)=2p")
call expect (abs(v_ff(one,v(m,p,-1),u(m,p,-1))-2*vp), 0, "ubar(-).V.u(-)=2p")
call expect (abs(v_ff(one,u(m,p,+1),v(m,p,+1))-2*vp), 0, "vbar(+).V.v(+)=2p")
call expect (abs(v_ff(one,u(m,p,-1),v(m,p,-1))-2*vp), 0, "vbar(-).V.v(-)=2p")
print *, "*** Checking current conservation ***:"
call expect ((vp-vq)*v_ff(one,v(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).V.u(+))=0")
call expect ((vp-vq)*v_ff(one,v(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).V.u(-))=0")
call expect ((vp-vq)*v_ff(one,u(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).V.v(+))=0")
call expect ((vp-vq)*v_ff(one,u(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).V.v(-))=0")
if (m == 0) then
   print *, "*** Checking axial current conservation ***:"
   call expect ((vp-vq)*a_ff(one,v(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).A.u(+))=0")
   call expect ((vp-vq)*a_ff(one,v(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).A.u(-))=0")
   call expect ((vp-vq)*a_ff(one,u(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).A.v(+))=0")
   call expect ((vp-vq)*a_ff(one,u(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).A.v(-))=0")
end if
print *, "*** Checking polarization vectors: ***"
call expect (conjg(eps(m,p, 1))*eps(m,p, 1), -1, "e( 1).e( 1)=-1")
call expect (conjg(eps(m,p, 1))*eps(m,p,-1),  0, "e( 1).e(-1)= 0")
call expect (conjg(eps(m,p,-1))*eps(m,p, 1),  0, "e(-1).e( 1)= 0")
call expect (conjg(eps(m,p,-1))*eps(m,p,-1), -1, "e(-1).e(-1)=-1")
call expect (                 p*eps(m,p, 1),  0, "    p.e( 1)= 0")
call expect (                 p*eps(m,p,-1),  0, "    p.e(-1)= 0")
if (m > 0) then
   call expect (conjg(eps(m,p, 1))*eps(m,p, 0),  0, "e( 1).e( 0)= 0")
   call expect (conjg(eps(m,p, 0))*eps(m,p, 1),  0, "e( 0).e( 1)= 0")
   call expect (conjg(eps(m,p, 0))*eps(m,p, 0), -1, "e( 0).e( 0)=-1")
   call expect (conjg(eps(m,p, 0))*eps(m,p,-1),  0, "e( 0).e(-1)= 0")
   call expect (conjg(eps(m,p,-1))*eps(m,p, 0),  0, "e(-1).e( 0)= 0")
   call expect (                 p*eps(m,p, 0),  0, "    p.e( 0)= 0")
end if
print *, "*** Checking polarization vectorspinors: ***"
call expect (abs(p * ueps(m, p,  2)),  0, "p.ueps ( 2)= 0")
call expect (abs(p * ueps(m, p,  1)),  0, "p.ueps ( 1)= 0")
call expect (abs(p * ueps(m, p, -1)),  0, "p.ueps (-1)= 0")
call expect (abs(p * ueps(m, p, -2)),  0, "p.ueps (-2)= 0")
call expect (abs(p * veps(m, p,  2)),  0, "p.veps ( 2)= 0")
call expect (abs(p * veps(m, p,  1)),  0, "p.veps ( 1)= 0")
call expect (abs(p * veps(m, p, -1)),  0, "p.veps (-1)= 0")
call expect (abs(p * veps(m, p, -2)),  0, "p.veps (-2)= 0")
print *, "*** in the rest frame ***"
call expect (abs(p_0 * ueps(m, p_0,  2)),  0, "p0.ueps ( 2)= 0")
call expect (abs(p_0 * ueps(m, p_0,  1)),  0, "p0.ueps ( 1)= 0")
call expect (abs(p_0 * ueps(m, p_0, -1)),  0, "p0.ueps (-1)= 0")
call expect (abs(p_0 * ueps(m, p_0, -2)),  0, "p0.ueps (-2)= 0")
call expect (abs(p_0 * veps(m, p_0,  2)),  0, "p0.veps ( 2)= 0")
call expect (abs(p_0 * veps(m, p_0,  1)),  0, "p0.veps ( 1)= 0")
call expect (abs(p_0 * veps(m, p_0, -1)),  0, "p0.veps (-1)= 0")
call expect (abs(p_0 * veps(m, p_0, -2)),  0, "p0.veps (-2)= 0")
print *, "*** Checking the irreducibility condition: ***"
call expect (abs(f_potgr (one, one, ueps(m, p,  2))),  0, "g.ueps ( 2)")
call expect (abs(f_potgr (one, one, ueps(m, p,  1))),  0, "g.ueps ( 1)")
call expect (abs(f_potgr (one, one, ueps(m, p, -1))),  0, "g.ueps (-1)")
call expect (abs(f_potgr (one, one, ueps(m, p, -2))),  0, "g.ueps (-2)")
call expect (abs(f_potgr (one, one, veps(m, p,  2))),  0, "g.veps ( 2)")
call expect (abs(f_potgr (one, one, veps(m, p,  1))),  0, "g.veps ( 1)")
call expect (abs(f_potgr (one, one, veps(m, p, -1))),  0, "g.veps (-1)")
call expect (abs(f_potgr (one, one, veps(m, p, -2))),  0, "g.veps (-2)")
print *, "*** in the rest frame ***"
call expect (abs(f_potgr (one, one, ueps(m, p_0,  2))),  0, "g.ueps ( 2)")
call expect (abs(f_potgr (one, one, ueps(m, p_0,  1))),  0, "g.ueps ( 1)")
call expect (abs(f_potgr (one, one, ueps(m, p_0, -1))),  0, "g.ueps (-1)")
call expect (abs(f_potgr (one, one, ueps(m, p_0, -2))),  0, "g.ueps (-2)")
call expect (abs(f_potgr (one, one, veps(m, p_0,  2))),  0, "g.veps ( 2)")
call expect (abs(f_potgr (one, one, veps(m, p_0,  1))),  0, "g.veps ( 1)")
call expect (abs(f_potgr (one, one, veps(m, p_0, -1))),  0, "g.veps (-1)")
call expect (abs(f_potgr (one, one, veps(m, p_0, -2))),  0, "g.veps (-2)")
print *, "*** Testing vectorspinor normalization ***"
call expect (veps(m,p, 2)*ueps(m,p, 2), -2*m, "ueps( 2).ueps( 2)= -2m")
call expect (veps(m,p, 1)*ueps(m,p, 1), -2*m, "ueps( 1).ueps( 1)= -2m")
call expect (veps(m,p,-1)*ueps(m,p,-1), -2*m, "ueps(-1).ueps(-1)= -2m")
call expect (veps(m,p,-2)*ueps(m,p,-2), -2*m, "ueps(-2).ueps(-2)= -2m")
call expect (ueps(m,p, 2)*veps(m,p, 2),  2*m, "veps( 2).veps( 2)= +2m")
call expect (ueps(m,p, 1)*veps(m,p, 1),  2*m, "veps( 1).veps( 1)= +2m")
call expect (ueps(m,p,-1)*veps(m,p,-1),  2*m, "veps(-1).veps(-1)= +2m")
call expect (ueps(m,p,-2)*veps(m,p,-2),  2*m, "veps(-2).veps(-2)= +2m")
call expect (ueps(m,p, 2)*ueps(m,p, 2),    0, "ueps( 2).veps( 2)=   0")
call expect (ueps(m,p, 1)*ueps(m,p, 1),    0, "ueps( 1).veps( 1)=   0")
call expect (ueps(m,p,-1)*ueps(m,p,-1),    0, "ueps(-1).veps(-1)=   0")
call expect (ueps(m,p,-2)*ueps(m,p,-2),    0, "ueps(-2).veps(-2)=   0")
call expect (veps(m,p, 2)*veps(m,p, 2),    0, "veps( 2).ueps( 2)=   0")
call expect (veps(m,p, 1)*veps(m,p, 1),    0, "veps( 1).ueps( 1)=   0")
call expect (veps(m,p,-1)*veps(m,p,-1),    0, "veps(-1).ueps(-1)=   0")
call expect (veps(m,p,-2)*veps(m,p,-2),    0, "veps(-2).ueps(-2)=   0")
print *, "*** in the rest frame ***"
call expect (veps(m,p_0, 2)*ueps(m,p_0, 2), -2*m, "ueps( 2).ueps( 2)= -2m")
call expect (veps(m,p_0, 1)*ueps(m,p_0, 1), -2*m, "ueps( 1).ueps( 1)= -2m")
call expect (veps(m,p_0,-1)*ueps(m,p_0,-1), -2*m, "ueps(-1).ueps(-1)= -2m")
call expect (veps(m,p_0,-2)*ueps(m,p_0,-2), -2*m, "ueps(-2).ueps(-2)= -2m")
call expect (ueps(m,p_0, 2)*veps(m,p_0, 2),  2*m, "veps( 2).veps( 2)= +2m")
call expect (ueps(m,p_0, 1)*veps(m,p_0, 1),  2*m, "veps( 1).veps( 1)= +2m")
call expect (ueps(m,p_0,-1)*veps(m,p_0,-1),  2*m, "veps(-1).veps(-1)= +2m")
call expect (ueps(m,p_0,-2)*veps(m,p_0,-2),  2*m, "veps(-2).veps(-2)= +2m")
call expect (ueps(m,p_0, 2)*ueps(m,p_0, 2),    0, "ueps( 2).veps( 2)=   0")
call expect (ueps(m,p_0, 1)*ueps(m,p_0, 1),    0, "ueps( 1).veps( 1)=   0")
call expect (ueps(m,p_0,-1)*ueps(m,p_0,-1),    0, "ueps(-1).veps(-1)=   0")
call expect (ueps(m,p_0,-2)*ueps(m,p_0,-2),    0, "ueps(-2).veps(-2)=   0")
call expect (veps(m,p_0, 2)*veps(m,p_0, 2),    0, "veps( 2).ueps( 2)=   0")
call expect (veps(m,p_0, 1)*veps(m,p_0, 1),    0, "veps( 1).ueps( 1)=   0")
call expect (veps(m,p_0,-1)*veps(m,p_0,-1),    0, "veps(-1).ueps(-1)=   0")
call expect (veps(m,p_0,-2)*veps(m,p_0,-2),    0, "veps(-2).ueps(-2)=   0")
print *, "*** Majorana properties of gravitino vertices: ***"
call expect (abs(u (m,q,1) * f_sgr (one, one, ueps(m,p,2), t) + & 
   ueps(m,p,2) * gr_sf(one,one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,-1) * f_sgr (one, one, ueps(m,p,2), t) + & 
!!!    ueps(m,p,2) * gr_sf(one,one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,1) * f_sgr (one, one, ueps(m,p,1), t) + & 
!!!    ueps(m,p,1) * gr_sf(one,one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,-1) * f_sgr (one, one, ueps(m,p,1), t) + & 
!!!    ueps(m,p,1) * gr_sf(one,one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,1) * f_sgr (one, one, ueps(m,p,-1), t) + & 
!!!    ueps(m,p,-1) * gr_sf(one,one,u(m,q,1),t)),  0, "f_sgr   + gr_sf       = 0")
!!! call expect (abs(u (m,q,-1) * f_sgr (one, one, ueps(m,p,-1), t) + & 
!!!    ueps(m,p,-1) * gr_sf(one,one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,1) * f_sgr (one, one, ueps(m,p,-2), t) + & 
!!!    ueps(m,p,-2) * gr_sf(one,one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0")
!!! call expect (abs(u (m,q,-1) * f_sgr (one, one, ueps(m,p,-2), t) + & 
!!!    ueps(m,p,-2) * gr_sf(one,one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0")
call expect (abs(u (m,q,1) * f_slgr (one, one, ueps(m,p,2), t) + & 
   ueps(m,p,2) * gr_slf(one,one,u(m,q,1),t)),  0, "f_slgr    + gr_slf    = 0")
call expect (abs(u (m,q,1) * f_srgr (one, one, ueps(m,p,2), t) + & 
   ueps(m,p,2) * gr_srf(one,one,u(m,q,1),t)),  0, "f_srgr    + gr_srf    = 0")
call expect (abs(u (m,q,1) * f_slrgr (one, two, one, ueps(m,p,2), t) + & 
   ueps(m,p,2) * gr_slrf(one,two,one,u(m,q,1),t)),  0, "f_slrgr   + gr_slrf   = 0")
call expect (abs(u (m,q,1) * f_pgr (one, one, ueps(m,p,2), t) + & 
   ueps(m,p,2) * gr_pf(one,one,u(m,q,1),t)),  0, "f_pgr     + gr_pf     = 0")
call expect (abs(u (m,q,1) * f_vgr (one, vt, ueps(m,p,2), p+q) + & 
   ueps(m,p,2) * gr_vf(one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0")
call expect (abs(u (m,q,1) * f_vlrgr (one, two, vt, ueps(m,p,2), p+q) + & 
   ueps(m,p,2) * gr_vlrf(one,two,vt,u(m,q,1),p+q)),  0, "f_vlrgr   + gr_vlrf   = 0")
!!! call expect (abs(u (m,q,-1) * f_vgr (one, vt, ueps(m,p,2), p+q) + & 
!!!    ueps(m,p,2) * gr_vf(one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(u (m,q,1) * f_vgr (one, vt, ueps(m,p,1), p+q) + & 
!!!    ueps(m,p,1) * gr_vf(one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(u (m,q,-1) * f_vgr (one, vt, ueps(m,p,1), p+q) + & 
!!!    ueps(m,p,1) * gr_vf(one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(u (m,q,1) * f_vgr (one, vt, ueps(m,p,-1), p+q) + & 
!!!    ueps(m,p,-1) * gr_vf(one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(u (m,q,-1) * f_vgr (one, vt, veps(m,p,-1), p+q) + & 
!!!    veps(m,p,-1) * gr_vf(one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(v (m,q,1) * f_vgr (one, vt, ueps(m,p,-2), p+q) + & 
!!!    ueps(m,p,-2) * gr_vf(one,vt,v(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0")
!!! call expect (abs(u (m,q,-1) * f_vgr (one, vt, ueps(m,p,-2), p+q) + & 
!!!    ueps(m,p,-2) * gr_vf(one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0")
call expect (abs(s_grf (one, ueps(m,p,2), u(m,q,1),t) + & 
   s_fgr(one,u(m,q,1),ueps(m,p,2),t)),  0, "s_grf     + s_fgr     = 0")
call expect (abs(sl_grf (one, ueps(m,p,2), u(m,q,1),t) + & 
   sl_fgr(one,u(m,q,1),ueps(m,p,2),t)),  0, "sl_grf    + sl_fgr    = 0")
call expect (abs(sr_grf (one, ueps(m,p,2), u(m,q,1),t) + & 
   sr_fgr(one,u(m,q,1),ueps(m,p,2),t)),  0, "sr_grf    + sr_fgr    = 0")
call expect (abs(slr_grf (one, two, ueps(m,p,2), u(m,q,1),t) + & 
   slr_fgr(one,two,u(m,q,1),ueps(m,p,2),t)),  0, "slr_grf   + slr_fgr   = 0")
call expect (abs(p_grf (one, ueps(m,p,2), u(m,q,1),t) + & 
   p_fgr(one,u(m,q,1),ueps(m,p,2),t)),  0, "p_grf     + p_fgr     = 0")
call expect (abs(v_grf (one, ueps(m,p,2), u(m,q,1),t) + & 
   v_fgr(one,u(m,q,1),ueps(m,p,2),t)),  0, "v_grf     + v_fgr     = 0")
call expect (abs(vlr_grf (one, two, ueps(m,p,2), u(m,q,1),t) + & 
   vlr_fgr(one,two,u(m,q,1),ueps(m,p,2),t)),  0, "vlr_grf   + vlr_fgr   = 0")
call expect (abs(u(m,p,1) * f_potgr (one,one,testv) - testv * gr_potf &
   (one,one,u (m,p,1))), 0, "f_potgr   - gr_potf   = 0")
call expect (abs (pot_fgr (one,u(m,p,1),testv) - pot_grf(one, & 
   testv,u(m,p,1))), 0, "pot_fgr   - pot_grf   = 0")
call expect (abs(u(m,p,1) * f_s2gr (one,one,one,testv) - testv * gr_s2f &
   (one,one,one,u (m,p,1))), 0, "f_s2gr    - gr_s2f    = 0")
call expect (abs (s2_fgr (one,u(m,p,1),one,testv) - s2_grf(one, &
   testv,one,u(m,p,1))), 0, "s2_fgr    - s2_grf    = 0")
call expect (abs(u (m,q,1) * f_svgr (one, one, vt, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_svf(one,one,vt,u(m,q,1))),  0, "f_svgr    + gr_svf    = 0")
call expect (abs(u (m,q,1) * f_slvgr (one, one, vt, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_slvf(one,one,vt,u(m,q,1))),  0, "f_slvgr   + gr_slvf   = 0")
call expect (abs(u (m,q,1) * f_srvgr (one, one, vt, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_srvf(one,one,vt,u(m,q,1))),  0, "f_srvgr   + gr_srvf   = 0")
call expect (abs(u (m,q,1) * f_slrvgr (one, two, one, vt, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_slrvf(one,two,one,vt,u(m,q,1))),  0, "f_slrvgr  + gr_slrvf  = 0")
call expect (abs (sv1_fgr (one,u(m,p,1),vt,ueps(m,q,2)) + sv1_grf(one, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "sv1_fgr   + sv1_grf   = 0")
call expect (abs (sv2_fgr (one,u(m,p,1),one,ueps(m,q,2)) + sv2_grf(one, &
   ueps(m,q,2),one,u(m,p,1))), 0, "sv2_fgr   + sv2_grf   = 0")
call expect (abs (slv1_fgr (one,u(m,p,1),vt,ueps(m,q,2)) + slv1_grf(one, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "slv1_fgr  + slv1_grf  = 0")
call expect (abs (srv2_fgr (one,u(m,p,1),one,ueps(m,q,2)) + srv2_grf(one, &
   ueps(m,q,2),one,u(m,p,1))), 0, "srv2_fgr  + srv2_grf  = 0")
call expect (abs (slrv1_fgr (one,two,u(m,p,1),vt,ueps(m,q,2)) + slrv1_grf(one,two, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "slrv1_fgr + slrv1_grf = 0")
call expect (abs (slrv2_fgr (one,two,u(m,p,1),one,ueps(m,q,2)) + slrv2_grf(one, &
   two,ueps(m,q,2),one,u(m,p,1))), 0, "slrv2_fgr + slrv2_grf = 0")
call expect (abs(u (m,q,1) * f_pvgr (one, one, vt, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_pvf(one,one,vt,u(m,q,1))),  0, "f_pvgr    + gr_pvf    = 0")
call expect (abs (pv1_fgr (one,u(m,p,1),vt,ueps(m,q,2)) + pv1_grf(one, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "pv1_fgr   + pv1_grf   = 0")
call expect (abs (pv2_fgr (one,u(m,p,1),one,ueps(m,q,2)) + pv2_grf(one, &
   ueps(m,q,2),one,u(m,p,1))), 0, "pv2_fgr   + pv2_grf   = 0")
call expect (abs(u (m,q,1) * f_v2gr (one, vt, vz, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_v2f(one,vt,vz,u(m,q,1))),  0, "f_v2gr    + gr_v2f    = 0")
call expect (abs(u (m,q,1) * f_v2lrgr (one, two, vt, vz, ueps(m,p,2)) + & 
   ueps(m,p,2) * gr_v2lrf(one,two,vt,vz,u(m,q,1))),  0, "f_v2lrgr  + gr_v2lrf  = 0")
call expect (abs (v2_fgr (one,u(m,p,1),vt,ueps(m,q,2)) + v2_grf(one, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "v2_fgr    + v2_grf    = 0")
call expect (abs (v2lr_fgr (one,two,u(m,p,1),vt,ueps(m,q,2)) + v2lr_grf(one, two, &
   ueps(m,q,2),vt,u(m,p,1))), 0, "v2lr_fgr  + v2lr_grf  = 0")
print *, "*** Testing the gravitino propagator: ***"
print *, "Transversality:"
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,testv))), 0, "p.pr.test")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,ueps(m,p,2)))),  0, "p.pr.ueps ( 2)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,ueps(m,p,1)))),  0, "p.pr.ueps ( 1)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,ueps(m,p,-1)))), 0, "p.pr.ueps (-1)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,ueps(m,p,-2)))), 0, "p.pr.ueps (-2)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,veps(m,p,2)))),  0, "p.pr.veps ( 2)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,veps(m,p,1)))),  0, "p.pr.veps ( 1)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,veps(m,p,-1)))), 0, "p.pr.veps (-1)")
call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
             pr_grav(p,m,w,veps(m,p,-2)))), 0, "p.pr.veps (-2)")
print *, "Irreducibility:"
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,testv)))), 0, "g.pr.test")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,ueps(m,p,2))))), 0, &
             "g.pr.ueps ( 2)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,ueps(m,p,1))))), 0, &
             "g.pr.ueps ( 1)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,ueps(m,p,-1))))), 0, &
             "g.pr.ueps (-1)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,ueps(m,p,-2))))), 0, &
             "g.pr.ueps (-2)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,veps(m,p,2))))), 0, &
             "g.pr.veps ( 2)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,veps(m,p,1))))), 0, &
             "g.pr.veps ( 1)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,veps(m,p,-1))))), 0, &
             "g.pr.veps (-1)")
call expect (abs(f_potgr (one, one, (cmplx (p*p - m**2, m*w, & 
             kind=omega_prec) * pr_grav(p,m,w,veps(m,p,-2))))), 0, &
             "g.pr.veps (-2)")
end program test_omega95_bispinors
