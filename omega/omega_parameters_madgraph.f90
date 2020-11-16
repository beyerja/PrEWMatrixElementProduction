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
module omega_parameters_madgraph
  use omega_kinds
  use omega_parameters
  implicit none
  private
  public :: export_parameters_to_madgraph
  integer, parameter, private :: D = selected_real_kind (14, 100)
  real(kind=D), save, public :: gw = 0, gwwa = 0, gwwz = 0
  real(kind=D), dimension(2), save, public :: gal = 0, gau = 0, gad = 0, gwf = 0
  real(kind=D), dimension(2), save, public :: gzn = 0, gzl = 0, gzu = 0, gzd = 0, g1 = 0
  real(kind=D), save, public :: gwwh = 0, gzzh = 0, ghhh = 0, &
        gwwhh = 0, gzzhh = 0, ghhhh = 0
  complex(kind=D), dimension(2,12), save, public :: gh = 0
  real(kind=D), save, public :: wmass = 0, wwidth = 0, zmass = 0, zwidth = 0
  real(kind=D), save, public :: amass = 0, awidth = 0, hmass = 0, hwidth = 0
  real(kind=D), dimension(12), save, public :: fmass = 0, fwidth = 0
  complex(kind=D), save, public :: fudge_m1 = 1, fudge_m2 = 1, fudge_m3 = 1, fudge_m4 = 1
contains
  subroutine export_parameters_to_madgraph ()
    gal = qlep
    gau = qup
    gad = qdwn
    gzl(1) = gnclep(1) + gnclep(2)
    gzl(2) = gnclep(1) - gnclep(2)
    gzn(1) = gncneu(1) + gncneu(2)
    gzn(2) = gncneu(1) - gncneu(2)
    gzu(1) = gncup(1) + gncup(2)
    gzu(2) = gncup(1) - gncup(2)
    gzd(1) = gncdwn(1) + gncdwn(2)
    gzd(2) = gncdwn(1) - gncdwn(2)
    gwf(1) = 2 * gcc
    gwf(2) = 0
    gwwa = qw
    gwwz = gzww
    gwwh = ghww
    !!! gwwhh = ghhww**2 !!! Old SM3
    gwwhh = ghhww
    gzzh = ghzz
    !!! gzzhh = ghhzz**2 !!! Old SM3
    gzzhh = ghhzz
    ghhh = gh3
    ghhhh = gh4
    ghtt = 0
    ghbb = 0
    ghcc = 0
    ghtautau = 0
    gh3 = 0
    gh4 = 0
    gh(:,1:6) = 0
    gh(:,7) = ghcc
    gh(:,8) = 0
    gh(:,9) = ghtautau
    gh(:,10) = 0
    gh(:,11) = ghtt
    gh(:,12) = ghbb
    fmass(1:2) = mass(11:12)
    fwidth(1:2) = width(11:12)
    fmass(5:6) = mass(13:14)
    fwidth(5:6) = width(13:14)
    fmass(9:10) = mass(15:16)
    fwidth(9:10) = width(15:16)
    fmass(4) = mass(1)
    fwidth(4) = width(1)
    fmass(3) = mass(2)
    fwidth(3) = width(2)
    fmass(8) = mass(3)
    fwidth(8) = width(3)
    fmass(7) = mass(4)
    fwidth(7) = width(4)
    fmass(12) = mass(5)
    fwidth(12) = width(5)
    fmass(11) = mass(6)
    fwidth(11) = width(6)
    amass = mass(22)
    awidth = width(22)
    zmass = mass(23)
    zwidth = width(23)
    wmass = mass(24)
    wwidth = width(24)
    hmass = mass(25)
    hwidth = width(25)
  end subroutine export_parameters_to_madgraph
end module omega_parameters_madgraph
