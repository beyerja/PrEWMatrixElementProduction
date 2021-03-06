subroutine TGCinit (TGCpar)
	use omega_kinds
	use omega_parameters
	
	real(kind=omega_prec), dimension(8), intent(in) :: TGCpar
	real(kind=omega_prec) :: Parg1a, Parg1z, Parg4a, Parg4z, Parka, Parkz, Parla, Parlz
	
	Parg1a = TGCpar(1)
	Parg1z = TGCpar(2)
	Parg4a = TGCpar(3)
	Parg4z = TGCpar(4)
	Parka = TGCpar(5)
	Parkz = TGCpar(6)
	Parla = TGCpar(7)
	Parlz = TGCpar(8)
	call setup_parameters()
	ig1a = iqw * Parg1a
    ig1z = igzww * Parg1z
    ig1pkpg4a = iqw   * (Parg1a + Parka + Parg4a) / 2
    ig1pkpg4z = igzww * (Parg1z + Parkz + Parg4z) / 2
    ig1pkmg4a = iqw   * (Parg1a + Parka - Parg4a) / 2
    ig1pkmg4z = igzww * (Parg1z + Parkz - Parg4z) / 2
    ig1mkpg4a = iqw   * (Parg1a - Parka + Parg4a) / 2
    ig1mkpg4z = igzww * (Parg1z - Parkz + Parg4z) / 2
    ig1mkmg4a = iqw   * (Parg1a - Parka - Parg4a) / 2
    ig1mkmg4z = igzww * (Parg1z - Parkz - Parg4z) / 2
    ila = iqw   * Parla / (mass(24)*mass(24))
    ilz = igzww * Parlz / (mass(24)*mass(24))
end subroutine TGCinit

program grid_ww_sl0muq

use omega_kinds
use omega_parameters
use omega95
use MathOperations
use Rotation
use ww_sl0muq

implicit none
integer, dimension(4,1,16,1) :: zero_ct
integer, dimension(2) :: filestat
integer :: n,i,j,h,Ri
!integer :: fi, fo, f, hi, ho,
integer :: icosthq, iphq, iph, imWp, imWm, icosth, icosthl
integer :: costh_steps, ph_steps, costhq_steps, phq_steps, costhl_steps, phl_steps, mWp_steps, mWm_steps
integer :: costh_average_steps, costhl_average_steps
integer :: num_args, ix

real(kind=omega_prec), dimension(0:3,6) :: p 
real(kind=omega_prec), dimension(4,4) :: polarization
real(kind=omega_prec), dimension(10,8) :: anomalousTGCpar
real(kind=omega_prec), dimension(2,40) :: MatrixElement
real(kind=omega_prec), dimension(0:3) :: p1,p2,p3,p4,p5,p6
!real(kind=omega_prec), dimension(0:3) :: W,Wb
real(kind=omega_prec), dimension(0:3) :: Wpboost, Wmboost
real(kind=omega_prec), dimension(3) :: directionW,directionl,directionq
real(kind=omega_prec), dimension(3) :: p_Wm, p_Wp, p_em
real(kind=omega_prec), dimension(3) :: p_f_Wm_rot, p_f_Wp_rot
real(kind=omega_prec), dimension (3,3) :: Wm_rotation, Wp_rotation
real(kind=omega_prec), dimension(3) :: A
real(kind=omega_prec), dimension(8) :: TGCpar
real(kind=omega_prec), dimension(10) :: anomalousg, anomalousk, anomalousl
real(kind=omega_prec), dimension(4) :: rho_in
real(kind=omega_prec), dimension(16,4,4) :: rho_out
real(kind=omega_prec) :: costh_min, costh_max, costh_width, costh_center
real(kind=omega_prec) :: costhq_min, costhq_max, costhq_width, costhq_center
real(kind=omega_prec) :: costhl_min, costhl_max, costhl_width, costhl_center
real(kind=omega_prec) :: costh, ph, costhq, phq, costhl, phl, mWp, mWm, Wrange, E, ECMS
real(kind=omega_prec) :: currentWmomentum, TGCdev, TGCdevg, TGCdevk, TGCdevl, cw
character(len=12), dimension(:), allocatable :: args
character(len=1024), dimension(2) :: filename
character(len=16) :: energylabel

ECMS = 250._omega_prec
write(energylabel, *) int(ECMS)

TGCdev = 0.0001_8
TGCdevg = TGCdev
TGCdevk = TGCdev
TGCdevl = TGCdev

anomalousg = [ 1._8, 1._8+TGCdevg, 1._8, 1._8, 1._8-TGCdevg, 1._8, 1._8, 1._8+TGCdevg, 1._8, 1._8+TGCdevg ]
anomalousk = [ 1._8, 1._8, 1._8+TGCdevk, 1._8, 1._8, 1._8-TGCdevk, 1._8, 1._8+TGCdevk, 1._8+TGCdevk, 1._8 ]
anomalousl = [ 0._8, 0._8, 0._8, TGCdevl, 0._8, 0._8, -TGCdevl, 0._8, TGCdevl, TGCdevl ]

TGCpar = [ 1._8, 1._8, 0._8, 0._8, 1._8, 1._8, 0._8, 0._8 ]
call TGCinit( TGCpar )

cw = mass(24)/mass(23)
!print *, cw
do i=1,10
	anomalousTGCpar(i,1) = 1._8	
	anomalousTGCpar(i,2) = anomalousg(i)
	anomalousTGCpar(i,3) = 0._8
	anomalousTGCpar(i,4) = 0._8
	anomalousTGCpar(i,5) = anomalousk(i)
	anomalousTGCpar(i,6) = -(anomalousk(i)-1)*((1-cw**2)/cw**2) + anomalousg(i) 
	anomalousTGCpar(i,7) = anomalousl(i)
	anomalousTGCpar(i,8) = anomalousl(i)
end do


!write(*, '(1000F12.3)')( real(anomalousTGCpar(j,2)) ,j=1,10)
!write(*, '(1000F12.3)')( real(anomalousTGCpar(j,5)) ,j=1,10)
!write(*, '(1000F12.3)')( real(anomalousTGCpar(j,6)) ,j=1,10)
!write(*, '(1000F12.3)')( real(anomalousTGCpar(j,7)) ,j=1,10)
!write(*, '(1000F12.3)')( real(anomalousTGCpar(j,8)) ,j=1,10)

!call setup_parameters()
filename(1) = "grid_ww_sl0muq_leptonic"
filename(2) = "grid_ww_sl0muq_hadronic"

filename(1) = trim(filename(1)) // "_" // trim(adjustl(energylabel)) // "GeV"
filename(2) = trim(filename(2)) // "_" // trim(adjustl(energylabel)) // "GeV"

num_args = command_argument_count()
allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
do ix = 1, num_args
	call get_command_argument(ix,args(ix))
	filename(1) = trim(filename(1)) // "_" // trim(args(ix))
	filename(2) = trim(filename(2)) // "_" // trim(args(ix))
	read (args(ix),*) A(ix)
end do
filename(1) = trim(filename(1)) // ".txt"
filename(2) = trim(filename(2)) // ".txt"

polarization(1:4,1) = [ 1._omega_prec, 0._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,2) = [ 0._omega_prec, 1._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,3) = [ 0._omega_prec, 0._omega_prec, 1._omega_prec, 0._omega_prec ]
polarization(1:4,4) = [ 0._omega_prec, 0._omega_prec, 0._omega_prec, 1._omega_prec ]

! incoming e-
p(0:3,1) = [ 0.5_omega_prec * ECMS, 0._omega_prec, 0._omega_prec,  0.5_omega_prec * ECMS] 
p_em = p(1:3,1)
! incoming e+
p(0:3,2) = [ 0.5_omega_prec * ECMS, 0._omega_prec, 0._omega_prec, -0.5_omega_prec * ECMS]

!th_steps = 20
!thl_steps = 10
!thq_steps = 10
costh_steps = 20
costhl_steps = 10
costhq_steps = 10
costh_average_steps = 5
costhl_average_steps = 5
phl_steps = 10
phq_steps = 5
ph_steps = 6
mWp_steps = 1
mWm_steps = 1
Wrange = 2 * width(24)

open(30,file=trim(filename(1)),iostat=filestat(1))
open(31,file=trim(filename(2) ),iostat=filestat(2) )
if(filestat(1)==0 .AND. filestat(2)==0) then

	!th = ( PI / dble(th_steps) ) * A(1)
	!thl = ( PI / dble(thl_steps) ) * A(2)
	!costh =  ( 1.9_omega_prec / dble(costh_steps-1_omega_prec)  ) * (A(1)-1_omega_prec) - 0.95_omega_prec
	!costhl = ( 1.9_omega_prec / dble(costhl_steps-1_omega_prec) ) * (A(2)-1_omega_prec) - 0.95_omega_prec	
	costh_min =  -1_omega_prec + ( 1.5_omega_prec / dble(costh_steps-1) )
	costhl_min = -1_omega_prec + ( 1.5_omega_prec / dble(costhl_steps-1) )
	costh_max =  1_omega_prec - ( 1.5_omega_prec / dble(costh_steps-1) )
	costhl_max = 1_omega_prec - ( 1.5_omega_prec / dble(costhl_steps-1) )
	costh_width = ( costh_max - costh_min) / dble(costh_steps-1)
	costhl_width = ( costhl_max - costhl_min) / dble(costhl_steps-1)
	costh_center =  costh_width * (A(1)-1) + costh_min
	costhl_center = costhl_width * (A(2)-1) + costhl_min
	phl = ( 2_omega_prec * PI / dble(phl_steps-1) ) * (A(3)-1)
	
	!directionl(1) = sin(thl)*cos(phl)
	!directionl(2) = sin(thl)*sin(phl)
	!directionl(3) = cos(thl)
	
	!print *, directionl(1:3)
	do icosthq=1,costhq_steps
		!thq = ( PI / dble(thq_steps) ) * dble( ithq )
		costhq = ( 2_omega_prec / dble(costhq_steps) ) * dble( icosthq ) - 1_omega_prec
	do iphq=1,phq_steps
		phq = ( 2_omega_prec * PI / dble(phq_steps) ) * dble( iphq )
		!directionq(1) = sin(thq)*cos(phq)
		!directionq(2) = sin(thq)*sin(phq)
		!directionq(3) = cos(thq)
		directionq(1) = sqrt(1 - costhq**2)*cos(phq)
		directionq(2) = sqrt(1 - costhq**2)*sin(phq)
		directionq(3) = costhq
		!print *, directionq(1:3)
		MatrixElement(1:2,1:40) = 0
		do iph=1,ph_steps
			ph = ( 2_omega_prec * PI / dble(ph_steps) ) * dble( iph  )
			!directionW(1) = sin(th)*cos(ph)
			!directionW(2) = sin(th)*sin(ph)
			!directionW(3) = cos(th)
			!print *, directionW(1:3)
		do imWp=1,mWp_steps
		mWp = ( ( 2_omega_prec * Wrange / dble(mWp_steps+1) ) * dble( imWp ) ) + mass(24) - Wrange
		do imWm=1,mWm_steps
		mWm = ( ( 2_omega_prec * Wrange / dble(mWm_steps+1) ) * dble( imWm  ) ) + mass(24) - Wrange
		
		do icosth=0,costh_average_steps-1
		do icosthl=0,costhl_average_steps-1
			
			costh = ( costh_width / dble(costh_average_steps) ) * icosth + costh_center - (0.5_omega_prec*costh_width)
			costh = costh + ( 0.5_omega_prec*costh_width / dble(costh_average_steps) )
			costhl = ( costhl_width / dble(costhl_average_steps) ) * icosthl + costhl_center - (0.5_omega_prec*costhl_width)
			costhl = costhl + ( 0.5_omega_prec*costhl_width / dble(costhl_average_steps) )
			
			directionW(1) = sqrt(1 - costh**2)*cos(ph)
			directionW(2) = sqrt(1 - costh**2)*sin(ph)
			directionW(3) = costh
			
			directionl(1) = sqrt(1 - costhl**2)*cos(phl)
			directionl(2) = sqrt(1 - costhl**2)*sin(phl)
			directionl(3) = costhl
			
			!print *, " - " 		
			
      ! Calculate the W 4-momenta (2-body kinematics in ee system)
			E = p(0,1) + p(0,2)
			Wmboost(0) = ( E + ( ( mWm**2 - mWp**2 ) / E ) )  / 2_omega_prec
      Wpboost(0) = E - Wmboost(0)
      
			currentWmomentum = sqrt( Wmboost(0)**2 - mWm**2 )
      p_Wm = - currentWmomentum * directionW(1:3)
      p_Wp = + currentWmomentum * directionW(1:3)
      
			Wmboost(1:3) = p_Wm
      Wpboost(1:3) = -Wmboost(1:3)
      
      ! Rotation matrices to rotate out of the W rest frame coordinates in which W flight is z axis
			Wm_rotation = rotate_out_of(p_em,p_Wm)
			Wp_rotation = rotate_out_of(p_em,p_Wp)
      
			!print *, Wmboost(0:3)
			!print *, Wpboost(0:3)		
      !print *, " - " 		
			
      ! Calculate the fermion momenta in the W rest frame, includes undoing rotation
      ! Here: Case that lepton is in W- system, below when it is in W+
			p5(0) = ( mWm + ( ( mass(13)**2 - mass(14)**2 ) / mWm ) )  / 2_omega_prec
			currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
      p_f_Wm_rot = currentWmomentum * directionl(1:3)
      p_f_Wm_rot = matmul(Wm_rotation, p_f_Wm_rot)
			p5(1:3) = p_f_Wm_rot
			
			!print *, p5(0:3)
			
			p6(0) = mWm - p5(0)
			p6(1:3) = -p5(1:3)

			!print *, p6(0:3)
			
			p4(0) = ( mWp + ( ( mass(1)**2 - mass(2)**2 ) / mWp ) )  / 2_omega_prec
			currentWmomentum = sqrt( p4(0)**2 - mass(1)**2 )
      p_f_Wp_rot = currentWmomentum * directionq(1:3)
      p_f_Wp_rot = matmul(Wp_rotation, p_f_Wp_rot)
			p4(1:3) = p_f_Wp_rot

			!print *, p4(0:3)
			
			p3(0) = mWp - p4(0)
			p3(1:3) = -p4(1:3)
			
			!print *, p3(0:3)
			
			!print *, " - " 		
				
      ! Boost fermion momenta out of W system into detector/lab frame
			call lorentzBOOST(Wpboost, p3, p(0:3,3) )
			call lorentzBOOST(Wpboost, p4, p(0:3,4) )
			call lorentzBOOST(Wmboost, p5, p(0:3,5) )
			call lorentzBOOST(Wmboost, p6, p(0:3,6) )			
			
			!do i=3,6
			!	print *, p(0:3,i)
			!end do
			!print *, " - " 		
			
			zero_ct=0
			n=0
			rho_in(1:4)=1
			h=0
			do Ri=1,10
				TGCpar(1:8) = anomalousTGCpar(Ri,1:8)
				call TGCinit( TGCpar )
				call scatter_diagonal_polarized (p, rho_in, rho_out, zero_ct, n)
				do j=1,4
					h = h + 1
					 do i = 1,16
					   MatrixElement(1,h) = MatrixElement(1,h) + rho_out(i,j,1)
					 end do
				 end do
			end do
			
			!do i=0,9
			!	print *, MatrixElement(1,(4*i)+1:(4*i)+4)-MatrixElement(1,1:4)
			!end do
			!print *, " "
			!print *, MatrixElement(1,5:8)-MatrixElement(1,17:20)
			
			p5(0) = ( mWp + ( ( mass(13)**2 - mass(14)**2 ) / mWp ) )  / 2_omega_prec
			currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
      p_f_Wp_rot = currentWmomentum * directionl(1:3)
      p_f_Wp_rot = matmul(Wp_rotation, p_f_Wp_rot)
			p5(1:3) = p_f_Wp_rot
			
			p6(0) = mWp - p5(0)
			p6(1:3) = -p5(1:3)

			p4(0) = ( mWm + ( ( mass(1)**2 - mass(2)**2 ) / mWm ) )  / 2_omega_prec
			currentWmomentum = sqrt( p4(0)**2 - mass(1)**2 )
      p_f_Wm_rot = currentWmomentum * directionq(1:3)
      p_f_Wm_rot = matmul(Wm_rotation, p_f_Wm_rot)
			p4(1:3) = p_f_Wm_rot
			
			p3(0) = mWm - p4(0)
			p3(1:3) = -p4(1:3)
									
			call lorentzBOOST(Wmboost, p3, p(0:3,3) )
			call lorentzBOOST(Wmboost, p4, p(0:3,4) )
			call lorentzBOOST(Wpboost, p5, p(0:3,5) )
			call lorentzBOOST(Wpboost, p6, p(0:3,6) )		
			
			zero_ct=0
			n=0
			rho_in(1:4)=1
			h=0
			do Ri=1,10
				TGCpar(1:8) = anomalousTGCpar(Ri,1:8)
				call TGCinit( TGCpar )
				call scatter_diagonal_polarized (p, rho_in, rho_out, zero_ct, n)
				do j=1,4
					h = h + 1
					 do i = 1,16
					   MatrixElement(2,h) = MatrixElement(2,h) + rho_out(i,j,4)
					 end do
				 end do	
			 end do
			 
			!do i=0,9
			!	print *, MatrixElement(2,(4*i)+1:(4*i)+4)
			!end do
			!print *, " "
		
		end do
		end do	
		end do
		end do
		end do	
		
		write(30,*) costh_center, costhl_center, phl, costhq, phq, MatrixElement(1,1:40)
		write(31,*) costh_center, costhl_center, phl, costhq, phq, MatrixElement(2,1:40)
	end do	
	end do
else
	write(*,*) "File not opened: Terminating"
end if
close(30)
end program grid_ww_sl0muq










