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
module omega_parameters
  use omega_kinds
  use omega_constants
  implicit none
  private
  public :: setup_parameters, print_parameters
  real(kind=omega_prec), dimension(37), save, public :: mass = 0, width = 0
  real(kind=omega_prec), parameter, public :: GeV = 1.0_omega_prec
  real(kind=omega_prec), parameter, public :: MeV = GeV / 1000
  real(kind=omega_prec), parameter, public :: keV = MeV / 1000
  real(kind=omega_prec), parameter, public :: TeV = GeV * 1000
  real(kind=omega_prec), save, public :: &
       alpha = 1.0_omega_prec / 137.0359895_omega_prec, &
       sin2thw = 0.23124_omega_prec
  !!! There is no fundamental reason in defining vev private; 
  !!! moreover it is needed for the K-matrix stuff. We also
  !!! need g, sinthw and costhw for this
  real(kind=omega_prec), save, public :: vev
  real(kind=omega_prec), save, public :: g, sinthw, costhw
  complex(kind=omega_prec), save, public :: &
       qlep = 0, qup = 0, qdwn = 0, gcc = 0, qw = 0, &
       gzww = 0, gwww = 0, ghww = 0, ghhww = 0, ghzz = 0, ghhzz = 0, &
       ghbb = 0, ghtt = 0, ghcc = 0, ghtautau = 0, gh3 = 0, gh4 = 0, &
       ghgaga = 0, ghgaz = 0, ghgg = 0, ghmm = 0, &             
       iqw = 0, igzww = 0, igwww = 0, &
       gw4 = 0, gzzww = 0, gazww = 0, gaaww = 0, &
       ig1a = 0, ig1z = 0, rg5a = 0, rg5z = 0, &
       ig1pkpg4a = 0, ig1pkpg4z = 0, ig1pkmg4a = 0, ig1pkmg4z = 0, &
       ig1mkpg4a = 0, ig1mkpg4z = 0, ig1mkmg4a = 0, ig1mkmg4z = 0, &
       ila = 0, ilz = 0, il5a = 0, il5z = 0, ik5a = 0, ik5z = 0, &
       ialww0 = 0, ialww2 = 0, ialzw0 = 0, ialzw1 = 0, ialzz = 0, &
       alww0 = 0, alww2 = 0, alzw0 = 0, alzw1 = 0, alzz = 0, &
       igdh4 = 0, gdh2w2 = 0, gdh2z2 = 0, gdhw2 = 0, gdhz2 = 0, &
       gs = 0, igs = 0
  complex(kind=omega_prec), save, public :: &
       sinckm12 = 0, sinckm13 = 0, sinckm23 = 0, &
       cosckm12 = 0, cosckm13 = 0, cosckm23 = 0
  complex(kind=omega_prec), save, public :: &
       vckm_11 = 0, vckm_12 = 0, vckm_13 = 0, vckm_21 = 0, &
       vckm_22 = 0, vckm_23 = 0, vckm_31 = 0, vckm_32 = 0, vckm_33 = 0
  complex(kind=omega_prec), save, public :: &
       gccq11 = 0, gccq12 = 0, gccq13 = 0, gccq21 = 0, &
       gccq22 = 0, gccq23 = 0, gccq31 = 0, gccq32 = 0, gccq33 = 0     
  real(kind=omega_prec), save, public :: &
       a4 = 0, a5 = 0, a6 = 0, a7 = 0, a10 = 0
  real(kind=omega_prec), save, public :: &
       g1a = 1, g1z = 1, kappaa = 1, kappaz = 1, lambdaa = 0, lambdaz = 0, &
       g4a = 0, g4z = 0, g5a = 0, g5z = 0, &
       kappa5a = 0, kappa5z = 0, lambda5a = 0, lambda5z = 0, &
       alpha4 = 0, alpha5 = 0, tau4 = 0, tau5 = 0
  real(kind=omega_prec), save, public :: xia = 1, xi0 = 1, xipm = 1
  real(kind=omega_prec), save, public :: kc0 = 0, kp0 = 0, kc1 = 0, &
       kp1 = 0, kc2 = 0, kp2 = 0
  real(kind=omega_prec), save, public :: lam_reg = 0   
  complex(kind=omega_prec), dimension(2), save, public :: &
       gnclep = 0, gncneu = 0, gncup = 0, gncdwn = 0
  complex(kind=omega_prec), save, public :: &
       fudge_o1 = 1, fudge_o2 = 1, fudge_o3 = 1, fudge_o4 = 1
  real(kind=omega_prec), save, public :: &
       fudge_higgs = 1, fudge_km = 1, w_res = 0
  real(kind=omega_prec), dimension(1:5), save, public :: &
       gkm, mkm, wkm
contains
  subroutine setup_parameters ()
    real(kind=omega_prec) :: e, qelep, qeup, qedwn
    mass(1) =   5.0 * MeV
    mass(2) =   3.0 * MeV
    mass(3) = 100.0 * MeV
    mass(4) =   1.2 * GeV
    mass(5) =   4.2 * GeV
    mass(6) = 174.0 * GeV
    width(1:5) = 0
    width(6) =  1.3 * GeV
    mass(11) =    0.51099907 * MeV
    mass(12) =    0
    mass(13) =  105.658389 * MeV
    mass(14) =    0
    mass(15) = 1777.05 * MeV
    mass(16) =    0
    width(11:16) = 0
    mass(21) =  0
    mass(22) =  0
    width(21:22) = 0
    mass(23) =  91.187 * GeV
    width(23) =  2.490 * GeV
    mass(24) =  80.41 * GeV
    width(24) =  2.06 * GeV
    mass(25) = 120.00 * GeV
    width(25) =  5.00 * GeV
    mass(35) = 10000 * GeV
    width(35) =  0
    sinckm12 = 0.0_omega_prec
    sinckm13 = 0.0_omega_prec
    sinckm23 = 0.0_omega_prec
    cosckm12 = sqrt ((1.0_omega_prec - (sinckm12**2)))
    cosckm13 = sqrt ((1.0_omega_prec - (sinckm13**2)))
    cosckm23 = sqrt ((1.0_omega_prec - (sinckm23**2)))
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    e = sqrt (4 * PI * alpha)
    qelep = - 1
    qeup = 2.0_omega_prec / 3.0_omega_prec
    qedwn = - 1.0_omega_prec / 3.0_omega_prec
    sinthw = sqrt (sin2thw)
    costhw = sqrt (1 - sin2thw)
    g = e / sinthw
    gcc = - g / 2 / sqrt (2.0_omega_prec)
    vckm_11 = cosckm12 * cosckm13
    vckm_12 = sinckm12 * cosckm13
    vckm_13 = sinckm13 
    vckm_21 =  - (sinckm12 * cosckm23 +  &
      cosckm12 * sinckm23 * sinckm13)
    vckm_22 = cosckm12 * cosckm23 -  &
      sinckm12 * sinckm23 * sinckm13
    vckm_23 = sinckm23 * cosckm13
    vckm_31 = sinckm12 * sinckm23 -  &
      cosckm12 * cosckm23 * sinckm13
    vckm_32 = - (cosckm12 * sinckm23 +  &
      sinckm12 * cosckm23 * sinckm13)
    vckm_33 = cosckm23 * cosckm13
    gccq11 = gcc * vckm_11
    gccq12 = gcc * vckm_12
    gccq13 = gcc * vckm_13
    gccq21 = gcc * vckm_21
    gccq22 = gcc * vckm_22
    gccq23 = gcc * vckm_23
    gccq31 = gcc * vckm_31
    gccq32 = gcc * vckm_32
    gccq33 = gcc * vckm_33
    gncneu(1) = - g / 2 / costhw * ( + 0.5_omega_prec)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_omega_prec - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_omega_prec - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_omega_prec - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_omega_prec)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_omega_prec)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_omega_prec)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_omega_prec)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qw = e
    iqw = (0,1)*qw
    gzww = g * costhw
    igzww = (0,1)*gzww
    gwww = g
    igwww = (0,1)*gwww
    ghww = mass(24) * g
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! This is for the old SM3:
    !!! ghhww = (0,1) * g / Sqrt(2.0_omega_prec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ghhww = g**2 / 2.0_omega_prec
    ghzz = mass(23) * g / costhw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! This is for the old SM3:
    !!! ghhzz = (0,1) * g / costhw / Sqrt(2.0_omega_prec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ghhzz = g**2 / 2.0_omega_prec / costhw**2
    gw4 = g**2
    gzzww = gzww**2
    gazww = gzww*e
    gaaww = e**2
    vev = 2.0 * mass(24) / g
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    gh3 = - 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! gh4 = mass(25) / vev !!! Old SM3
    ig1a = iqw * g1a
    ig1z = igzww * g1z
    ig1pkpg4a = iqw   * (g1a + kappaa + g4a) / 2
    ig1pkpg4z = igzww * (g1z + kappaz + g4z) / 2
    ig1pkmg4a = iqw   * (g1a + kappaa - g4a) / 2
    ig1pkmg4z = igzww * (g1z + kappaz - g4z) / 2
    ig1mkpg4a = iqw   * (g1a - kappaa + g4a) / 2
    ig1mkpg4z = igzww * (g1z - kappaz + g4z) / 2
    ig1mkmg4a = iqw   * (g1a - kappaa - g4a) / 2
    ig1mkmg4z = igzww * (g1z - kappaz - g4z) / 2
    ila = iqw   * lambdaa / (mass(24)*mass(24))
    ilz = igzww * lambdaz / (mass(24)*mass(24))
    rg5a = qw   * g5a
    rg5z = gzww * g5z
    ik5a = iqw   * kappa5a
    ik5z = igzww * kappa5z
    il5a = iqw   * lambda5a / (mass(24)*mass(24))
    il5z = igzww * lambda5z / (mass(24)*mass(24))
    alww0 = g**4 * (alpha4 + 2 * alpha5)
    alww2 = g**4 * 2 * alpha4
    alzw1 = g**4 / costhw**2 * alpha4
    alzw0 = g**4 / costhw**2 * 2 * alpha5
    alzz = g**4 / costhw**4 * 2 * (alpha4 + alpha5)
    ialww0 = g**2 * sqrt ( - cmplx (alpha4 + 2 * alpha5, kind=omega_prec))
    ialww2 = g**2 * sqrt ( - cmplx (2 * alpha4, kind=omega_prec))
    ialzw1 = g**2 / costhw * sqrt ( - cmplx (alpha4, kind=omega_prec))
    ialzw0 = g**2 / costhw * sqrt ( - cmplx (2 * alpha5, kind=omega_prec))
    ialzz  = g**2 / (costhw*costhw) &
                  * sqrt ( - cmplx (2 * (alpha4 + alpha5), kind=omega_prec))
    gdh2w2 = g * vev * sqrt (cmplx (tau4, kind=omega_prec))
    gdhw2 = g * vev * sqrt (cmplx (tau5 / 2, kind=omega_prec))
    gdh2z2 = g * vev / costhw * sqrt (cmplx (tau4, kind=omega_prec))
    gdhz2 = g * vev / costhw * sqrt (cmplx (tau5 / 2, kind=omega_prec))
    igdh4 = g**2 * sqrt ( - cmplx (8 * (tau4 + tau5), kind=omega_prec))
  end subroutine setup_parameters
  subroutine print_parameters ()
    print *, "Quark masses:"
    print *, mass(2:6:2)
    print *, mass(1:5:2)
    print *, "Lepton masses:"
    print *, mass(12:16:2)
    print *, mass(11:15:2)
    print *, "Quark widths:"
    print *, width(2:6:2)
    print *, width(1:5:2)
    print *, "Lepton widths:"
    print *, width(12:16:2)
    print *, width(11:15:2)
    print *, "SU(2)xU(1) Gauge boson masses/widths:"
    print *, mass(22:24)
    print *, width(22:24)
    print *, "Higgs boson and gluon masses/widths:"
    print *, mass(25), mass(21)
    print *, width(25), width(21) 
    print *, "Neutral current couplings:"
    print *, "U:", gncup
    print *, "D:", gncdwn
    print *, "N:", gncneu
    print *, "L:", gnclep
    print *, "Fermion charges:"
    print *, "U:", qup
    print *, "D:", qdwn
    print *, "L:", qlep
    print *, "TGC:"
    print *, "WWA:", iqw
    print *, "WWZ:", igzww
    print *, "WWW:", igwww
    print *, "WWH:", ghww
    !!! print *, "WWHH:", ghhww**2 !!! Old SM3
    print *, "WWHH:", ghhww
    !!! print *, "ZZHH:", ghhzz**2 !!! Old SM3
    print *, "ZZHH:", ghhzz
  end subroutine print_parameters
end module omega_parameters
