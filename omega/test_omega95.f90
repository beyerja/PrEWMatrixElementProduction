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
program test_omega95
  use omega_kinds
  use omega95
  use omega_testtools
  implicit none
  real(kind=omega_prec) :: m, pabs, qabs, w
  real(kind=omega_prec), dimension(0:3) :: r
  complex(kind=omega_prec) :: one
  type(momentum) :: p, q
  type(vector) :: vp, vq, vtest
  type(tensor) :: ttest
  integer, dimension(8) :: date_time
  integer :: rsize
  call date_and_time (values = date_time)
  call random_seed (size = rsize)
  call random_seed (put = spread (product (date_time), dim = 1, ncopies = rsize))
  w = 1.4142
  one = 1
  m = 13
  pabs = 42 
  qabs = 137
  call random_number (r)
  vtest%t = cmplx (10.0_omega_prec * r(0))
  vtest%x(1:3) = cmplx (10.0_omega_prec * r(1:3))
  ttest = vtest.tprod.vtest
  call random_momentum (p, pabs, m)
  call random_momentum (q, qabs, m)
  vp = p
  vq = q
  print *, "*** Checking the equations of motion ***:"
  call expect (abs(f_vf(one,vp,u(m,p,+1))-m*u(m,p,+1)), 0, "|[p-m]u(+)|=0")
  call expect (abs(f_vf(one,vp,u(m,p,-1))-m*u(m,p,-1)), 0, "|[p-m]u(-)|=0")
  call expect (abs(f_vf(one,vp,v(m,p,+1))+m*v(m,p,+1)), 0, "|[p+m]v(+)|=0")
  call expect (abs(f_vf(one,vp,v(m,p,-1))+m*v(m,p,-1)), 0, "|[p+m]v(-)|=0")
  call expect (abs(f_fv(one,ubar(m,p,+1),vp)-m*ubar(m,p,+1)), 0, "|ubar(+)[p-m]|=0")
  call expect (abs(f_fv(one,ubar(m,p,-1),vp)-m*ubar(m,p,-1)), 0, "|ubar(-)[p-m]|=0")
  call expect (abs(f_fv(one,vbar(m,p,+1),vp)+m*vbar(m,p,+1)), 0, "|vbar(+)[p+m]|=0")
  call expect (abs(f_fv(one,vbar(m,p,-1),vp)+m*vbar(m,p,-1)), 0, "|vbar(-)[p+m]|=0")
  print *, "*** Checking the normalization ***:"
  call expect (ubar(m,p,+1)*u(m,p,+1), +2*m, "ubar(+)*u(+)=+2m")
  call expect (ubar(m,p,-1)*u(m,p,-1), +2*m, "ubar(-)*u(-)=+2m")
  call expect (vbar(m,p,+1)*v(m,p,+1), -2*m, "vbar(+)*v(+)=-2m")
  call expect (vbar(m,p,-1)*v(m,p,-1), -2*m, "vbar(-)*v(-)=-2m")
  call expect (ubar(m,p,+1)*v(m,p,+1),    0, "ubar(+)*v(+)=0  ")
  call expect (ubar(m,p,-1)*v(m,p,-1),    0, "ubar(-)*v(-)=0  ")
  call expect (vbar(m,p,+1)*u(m,p,+1),    0, "vbar(+)*u(+)=0  ")
  call expect (vbar(m,p,-1)*u(m,p,-1),    0, "vbar(-)*u(-)=0  ")
  print *, "*** Checking the currents ***:"
  call expect (abs(v_ff(one,ubar(m,p,+1),u(m,p,+1))-2*vp), 0, "ubar(+).V.u(+)=2p")
  call expect (abs(v_ff(one,ubar(m,p,-1),u(m,p,-1))-2*vp), 0, "ubar(-).V.u(-)=2p")
  call expect (abs(v_ff(one,vbar(m,p,+1),v(m,p,+1))-2*vp), 0, "vbar(+).V.v(+)=2p")
  call expect (abs(v_ff(one,vbar(m,p,-1),v(m,p,-1))-2*vp), 0, "vbar(-).V.v(-)=2p")
  print *, "*** Checking current conservation ***:"
  call expect ((vp-vq)*v_ff(one,ubar(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).V.u(+))=0")
  call expect ((vp-vq)*v_ff(one,ubar(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).V.u(-))=0")
  call expect ((vp-vq)*v_ff(one,vbar(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).V.v(+))=0")
  call expect ((vp-vq)*v_ff(one,vbar(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).V.v(-))=0")
  if (m == 0) then
     print *, "*** Checking axial current conservation ***:"
     call expect ((vp-vq)*a_ff(one,ubar(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).A.u(+))=0")
     call expect ((vp-vq)*a_ff(one,ubar(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).A.u(-))=0")
     call expect ((vp-vq)*a_ff(one,vbar(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).A.v(+))=0")
     call expect ((vp-vq)*a_ff(one,vbar(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).A.v(-))=0")
  end if
  print *, "*** Checking polarisation vectors: ***"
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
     call expect (                p*eps(m,p, 0),  0, "    p.e( 0)= 0")
  end if
  print *, "*** Checking epsilon tensor: ***"
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,q,1),eps(m,p,1),eps(m,p,0),eps(m,q,0)), "eps(1<->2)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,p,0),eps(m,q,1),eps(m,p,1),eps(m,q,0)), "eps(1<->3)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,q,0),eps(m,q,1),eps(m,p,0),eps(m,p,1)), "eps(1<->4)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,p,1),eps(m,p,0),eps(m,q,1),eps(m,q,0)), "eps(2<->3)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,p,1),eps(m,q,0),eps(m,p,0),eps(m,q,1)), "eps(2<->4)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
               - pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,q,0),eps(m,p,0)), "eps(3<->4)")
  call expect (  pseudo_scalar(eps(m,p,1),eps(m,q,1),eps(m,p,0),eps(m,q,0)), &
                 eps(m,p,1)*pseudo_vector(eps(m,q,1),eps(m,p,0),eps(m,q,0)), "eps'")
  print *, "*** Checking tensors: ***"
  call expect (conjg(p.wedge.q)*(p.wedge.q), (p*p)*(q*q)-(p*q)**2, &
       "[p,q].[q,p]=p.p*q.q-p.q^2")
  call expect (conjg(p.wedge.q)*(q.wedge.p), (p*q)**2-(p*p)*(q*q), &
       "[p,q].[q,p]=p.q^2-p.p*q.q")
  call expect (conjg(p.wedge.eps(m,p, 1))*(p.wedge.eps(m,p, 1)), -p*p, &
       "[p,e( 1)].[p,e( 1)]=-p.p")
  call expect (conjg(p.wedge.eps(m,p, 1))*(p.wedge.eps(m,p,-1)),    0, &
       "[p,e( 1)].[p,e(-1)]=0")
  call expect (conjg(p.wedge.eps(m,p,-1))*(p.wedge.eps(m,p, 1)),    0, &
       "[p,e(-1)].[p,e( 1)]=0")
  call expect (conjg(p.wedge.eps(m,p,-1))*(p.wedge.eps(m,p,-1)), -p*p, &
       "[p,e(-1)].[p,e(-1)]=-p.p")
  if (m > 0) then
     call expect (conjg(p.wedge.eps(m,p, 1))*(p.wedge.eps(m,p, 0)),    0, &
          "[p,e( 1)].[p,e( 0)]=0")
     call expect (conjg(p.wedge.eps(m,p, 0))*(p.wedge.eps(m,p, 1)),    0, &
          "[p,e( 0)].[p,e( 1)]=0")
     call expect (conjg(p.wedge.eps(m,p, 0))*(p.wedge.eps(m,p, 0)), -p*p, &
          "[p,e( 0)].[p,e( 0)]=-p.p")
     call expect (conjg(p.wedge.eps(m,p, 0))*(p.wedge.eps(m,p,-1)),    0, &
          "[p,e( 1)].[p,e(-1)]=0")
     call expect (conjg(p.wedge.eps(m,p,-1))*(p.wedge.eps(m,p, 0)),    0, &
          "[p,e(-1)].[p,e( 0)]=0")
  end if
  call expect (abs ((p.wedge.eps(m,p, 1))*p + (p*p)*eps(m,p, 1)), 0, &
       "[p,e( 1)].p=-p.p*e( 1)]")
  call expect (abs ((p.wedge.eps(m,p, 0))*p + (p*p)*eps(m,p, 0)), 0, &
       "[p,e( 0)].p=-p.p*e( 0)]")
  call expect (abs ((p.wedge.eps(m,p,-1))*p + (p*p)*eps(m,p,-1)), 0, &
       "[p,e(-1)].p=-p.p*e(-1)]")
  call expect (abs (p*(p.wedge.eps(m,p, 1)) - (p*p)*eps(m,p, 1)), 0, &
       "p.[p,e( 1)]=p.p*e( 1)]")
  call expect (abs (p*(p.wedge.eps(m,p, 0)) - (p*p)*eps(m,p, 0)), 0, &
       "p.[p,e( 0)]=p.p*e( 0)]")
  call expect (abs (p*(p.wedge.eps(m,p,-1)) - (p*p)*eps(m,p,-1)), 0, &
       "p.[p,e(-1)]=p.p*e(-1)]")
  print *, "*** Checking polarisation tensors: ***"
  call expect (conjg(eps2(m,p, 2))*eps2(m,p, 2), 1, "e2( 2).e2( 2)=1")
  call expect (conjg(eps2(m,p, 2))*eps2(m,p,-2), 0, "e2( 2).e2(-2)=0")
  call expect (conjg(eps2(m,p,-2))*eps2(m,p, 2), 0, "e2(-2).e2( 2)=0")
  call expect (conjg(eps2(m,p,-2))*eps2(m,p,-2), 1, "e2(-2).e2(-2)=1")
  if (m > 0) then
     call expect (conjg(eps2(m,p, 2))*eps2(m,p, 1), 0, "e2( 2).e2( 1)=0")
     call expect (conjg(eps2(m,p, 2))*eps2(m,p, 0), 0, "e2( 2).e2( 0)=0")
     call expect (conjg(eps2(m,p, 2))*eps2(m,p,-1), 0, "e2( 2).e2(-1)=0")
     call expect (conjg(eps2(m,p, 1))*eps2(m,p, 2), 0, "e2( 1).e2( 2)=0")
     call expect (conjg(eps2(m,p, 1))*eps2(m,p, 1), 1, "e2( 1).e2( 1)=1")
     call expect (conjg(eps2(m,p, 1))*eps2(m,p, 0), 0, "e2( 1).e2( 0)=0")
     call expect (conjg(eps2(m,p, 1))*eps2(m,p,-1), 0, "e2( 1).e2(-1)=0")
     call expect (conjg(eps2(m,p, 1))*eps2(m,p,-2), 0, "e2( 1).e2(-2)=0")
     call expect (conjg(eps2(m,p, 0))*eps2(m,p, 2), 0, "e2( 0).e2( 2)=0")
     call expect (conjg(eps2(m,p, 0))*eps2(m,p, 1), 0, "e2( 0).e2( 1)=0")
     call expect (conjg(eps2(m,p, 0))*eps2(m,p, 0), 1, "e2( 0).e2( 0)=1")
     call expect (conjg(eps2(m,p, 0))*eps2(m,p,-1), 0, "e2( 0).e2(-1)=0")
     call expect (conjg(eps2(m,p, 0))*eps2(m,p,-2), 0, "e2( 0).e2(-2)=0")
     call expect (conjg(eps2(m,p,-1))*eps2(m,p, 2), 0, "e2(-1).e2( 2)=0")
     call expect (conjg(eps2(m,p,-1))*eps2(m,p, 1), 0, "e2(-1).e2( 1)=0")
     call expect (conjg(eps2(m,p,-1))*eps2(m,p, 0), 0, "e2(-1).e2( 0)=0")
     call expect (conjg(eps2(m,p,-1))*eps2(m,p,-1), 1, "e2(-1).e2(-1)=1")
     call expect (conjg(eps2(m,p,-1))*eps2(m,p,-2), 0, "e2(-1).e2(-2)=0")
     call expect (conjg(eps2(m,p,-2))*eps2(m,p, 1), 0, "e2(-2).e2( 1)=0")
     call expect (conjg(eps2(m,p,-2))*eps2(m,p, 0), 0, "e2(-2).e2( 0)=0")
     call expect (conjg(eps2(m,p,-2))*eps2(m,p,-1), 0, "e2(-2).e2(-1)=0")
  end if
  call expect (           abs(p*eps2(m,p, 2)  ), 0, " |p.e2( 2)|  =0")
  call expect (             abs(eps2(m,p, 2)*p), 0, "   |e2( 2).p|=0")
  call expect (           abs(p*eps2(m,p,-2)  ), 0, " |p.e2(-2)|  =0")
  call expect (             abs(eps2(m,p,-2)*p), 0, "   |e2(-2).p|=0")
  if (m > 0) then
     call expect (           abs(p*eps2(m,p, 1)  ), 0, " |p.e2( 1)|  =0")
     call expect (             abs(eps2(m,p, 1)*p), 0, "   |e2( 1).p|=0")
     call expect (           abs(p*eps2(m,p, 0)  ), 0, " |p.e2( 0)|  =0")
     call expect (             abs(eps2(m,p, 0)*p), 0, "   |e2( 0).p|=0")
     call expect (           abs(p*eps2(m,p,-1)  ), 0, " |p.e2(-1)|  =0")
     call expect (             abs(eps2(m,p,-1)*p), 0, "   |e2(-1).p|=0")
  end if
  print *, " *** Checking the graviton propagator:"
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,eps2(m,p,-2)))), 0, "p.pr.e(-2)")
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,eps2(m,p,-1)))), 0, "p.pr.e(-1)")
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,eps2(m,p,0)))), 0, "p.pr.e(0)")
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,eps2(m,p,1)))), 0, "p.pr.e(1)")
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,eps2(m,p,2)))), 0, "p.pr.e(2)")
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=omega_prec) * &
               pr_tensor(p,m,w,ttest))), 0, "p.pr.ttest")
end program test_omega95
