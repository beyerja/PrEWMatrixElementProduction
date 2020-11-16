program grid_zz_sl0mu_down

use omega_kinds
use omega_parameters
use omega95
use MathOperations
use zz_sl0mu_down

implicit none
integer, dimension(4,1,16,1) :: zero_ct
integer, dimension(2) :: filestat
integer :: n,i,j,h,Ri
!integer :: fi, fo, f, hi, ho,
integer :: ithq, iphq, iph, imZl, imZh
integer :: th_steps, ph_steps, thq_steps, phq_steps, thl_steps, phl_steps, mZl_steps, mZh_steps
integer :: num_args, ix

real(kind=omega_prec), dimension(0:3,6) :: p 
real(kind=omega_prec), dimension(4,4) :: polarization
real(kind=omega_prec), dimension(10,8) :: anomalousTGCpar
real(kind=omega_prec), dimension(2,4) :: MatrixElement
real(kind=omega_prec), dimension(0:3) :: p1,p2,p3,p4,p5,p6
!real(kind=omega_prec), dimension(0:3) :: W,Wb
real(kind=omega_prec), dimension(0:3) :: Zlboost, Zhboost
real(kind=omega_prec), dimension(3) :: directionW,directionl,directionq
real(kind=omega_prec), dimension(3) :: A
real(kind=omega_prec), dimension(8) :: TGCpar
real(kind=omega_prec), dimension(10) :: anomalousg, anomalousk, anomalousl
real(kind=omega_prec), dimension(4,1) :: rho_in
real(kind=omega_prec), dimension(16,4) :: rho_out
real(kind=omega_prec) :: th, ph, thq, phq, thl, phl, mZl, mZh, Zrange
real(kind=omega_prec) :: E, ECMS, currentWmomentum, cw
character(len=12), dimension(:), allocatable :: args
character(len=1024), dimension(1) :: filename


ECMS = 500._omega_prec

call setup_parameters()
filename(1) = "grid_zz_sl0mu_down"
num_args = command_argument_count()
allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
do ix = 1, num_args
	call get_command_argument(ix,args(ix))
	filename(1) = trim(filename(1)) // "_" // trim(args(ix))

	read (args(ix),*) A(ix)
end do
filename(1) = trim(filename(1)) // ".txt"


polarization(1:4,1) = [ 1._omega_prec, 0._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,2) = [ 0._omega_prec, 1._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,3) = [ 0._omega_prec, 0._omega_prec, 1._omega_prec, 0._omega_prec ]
polarization(1:4,4) = [ 0._omega_prec, 0._omega_prec, 0._omega_prec, 1._omega_prec ]

p(0:3,1) = [ 0.5_omega_prec * ECMS, 0._omega_prec, 0._omega_prec,  0.5_omega_prec * ECMS]
p(0:3,2) = [ 0.5_omega_prec * ECMS, 0._omega_prec, 0._omega_prec, -0.5_omega_prec * ECMS]

th_steps = 20
thl_steps = 10
phl_steps = 10
thq_steps = 10
phq_steps = 10
ph_steps = 360
mZl_steps = 1
mZh_steps = 1
Zrange = 2 * width(24)

open(30,file=trim(filename(1)),iostat=filestat(1))

if(filestat(1)==0 .AND. filestat(2)==0) then

	th = ( PI / dble(th_steps) ) * A(1)
	thl = ( PI / dble(thl_steps) ) * A(2)
	phl = ( 2_omega_prec * PI / dble(phl_steps) ) * A(3)
	directionl(1) = sin(thl)*cos(phl)
	directionl(2) = sin(thl)*sin(phl)
	directionl(3) = cos(thl)
	!print *, directionl(1:3)
	do ithq=1,thq_steps
		thq = ( PI / dble(thq_steps) ) * dble( ithq )
	do iphq=1,phq_steps
		phq = ( 2_omega_prec * PI / dble(phq_steps) ) * dble( iphq )
		directionq(1) = sin(thq)*cos(phq)
		directionq(2) = sin(thq)*sin(phq)
		directionq(3) = cos(thq)
		!print *, directionq(1:3)
		MatrixElement(1:2,1:4) = 0
		do iph=1,ph_steps
			ph = ( 2_omega_prec * PI / dble(ph_steps) ) * dble( iph  )
			directionW(1) = sin(th)*cos(ph)
			directionW(2) = sin(th)*sin(ph)
			directionW(3) = cos(th)
			!print *, directionW(1:3)
		do imZl=1,mZl_steps
		mZl = ( ( 2_omega_prec * Zrange / dble(mZl_steps+1) ) * dble( imZl ) ) + mass(24) - Zrange
		do imZh=1,mZh_steps
		mZh = ( ( 2_omega_prec * Zrange / dble(mZh_steps+1) ) * dble( imZh  ) ) + mass(24) - Zrange
			
			!print *, " - " 		
			
			E = p(0,1) + p(0,2)
			Zhboost(0) = ( E + ( ( mZh**2 - mZl**2 ) / E ) )  / 2_omega_prec
			currentWmomentum = sqrt( Zhboost(0)**2 - mZh**2 )
			Zhboost(1:3) = - currentWmomentum * directionW(1:3)
			
			!print *, Zhboost(0:3)
			
			Zlboost(0) = E - Zhboost(0)
			Zlboost(1:3) = -Zhboost(1:3)
			
			!print *, Zlboost(0:3)		
			
			!print *, " - " 		
			
			p5(0) = mZh / 2_omega_prec
			currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
			p5(1:3) = currentWmomentum * directionl(1:3)
			
			!print *, p5(0:3)
			
			p6(0) = mZh - p5(0)
			p6(1:3) = -p5(1:3)

			!print *, p6(0:3)
			
			p3(0) = mZl / 2_omega_prec
			currentWmomentum = sqrt( p3(0)**2 - mass(1)**2 )
			p3(1:3) = currentWmomentum * directionq(1:3)

			!print *, p3(0:3)
			
			p4(0) = mZl - p3(0)
			p4(1:3) = -p3(1:3)
					
				
			call lorentzBOOST(Zlboost, p3, p(0:3,3) )
			call lorentzBOOST(Zlboost, p4, p(0:3,4) )
			call lorentzBOOST(Zhboost, p5, p(0:3,5) )
			call lorentzBOOST(Zhboost, p6, p(0:3,6) )			
				
			
			zero_ct=0
			n=0
			
			do j=1,4
				rho_in(1:4,1) = polarization(1:4,j)
				call scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
				do h =1,3
					do i = 1,16
						MatrixElement(1,j) = MatrixElement(1,j) + rho_out(i,h)
					end do
				end do
			 end do
			
			
		end do
		end do
		end do
		
		
		write(30,*) th, thl, phl, thq, phq, MatrixElement(1,1:4)
	end do	
	end do
else
	write(*,*) "File not opened: Terminating"
end if
close(30)
end program grid_zz_sl0mu_down