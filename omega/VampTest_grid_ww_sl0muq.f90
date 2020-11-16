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

function Binnedelement (x, weights, channel,polarization) result (value)
	use kinds
	use omega_kinds
	use omega_parameters
	use MathOperations
	use ww_sl0muq
	use vamp_grid_type !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	real(kind=omega_prec), dimension(:), intent(in) :: polarization
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	real(kind=omega_prec) :: value
	integer, dimension(4,1,16,1) :: zero_ct
	integer :: n,i,j,h
	real(kind=omega_prec), dimension(0:3,6) :: p 
	real(kind=omega_prec), dimension(0:3) :: p1,p2,p3,p4,p5,p6
	real(kind=omega_prec), dimension(0:3) :: Wpboost, Wmboost
	real(kind=omega_prec), dimension(3) :: directionW,directionl,directionq
	real(kind=omega_prec), dimension(3) :: A
	real(kind=omega_prec), dimension(4,1) :: rho_in
	real(kind=omega_prec), dimension(16,4) :: rho_out
	real(kind=omega_prec) :: E, currentWmomentum, th, ph, thq, phq, thl, phl, mWp, mWm
	th = x(1)
	ph = x(2)
	thq = x(3)
	phq = x(4)
	thl = x(5)
	phl = x(6)
	mWm = mass(24) !x(7)
	mWp = mass(24) !x(8)
	directionl(1) = sin(thl)*cos(phl)
	directionl(2) = sin(thl)*sin(phl)
	directionl(3) = cos(thl)
	directionq(1) = sin(thq)*cos(phq)
	directionq(2) = sin(thq)*sin(phq)
	directionq(3) = cos(thq)
	directionW(1) = sin(th)*cos(ph)
	directionW(2) = sin(th)*sin(ph)
	directionW(3) = cos(th)
	
	!directionW(1) = x(1)
	!directionW(2) = x(2)
	!directionW(3) = sqrt( 1_omega_prec - directionW(1)**2 - directionW(2)**2 )
	!directionq(1) = x(3)
	!directionq(2) = x(4)
	!directionq(3) = sqrt( 1_omega_prec - directionq(1)**2 - directionq(2)**2 )
	!directionl(1) = x(5)
	!directionl(2) = x(6)
	!directionl(3) = sqrt( 1_omega_prec - directionl(1)**2 - directionl(2)**2 )
	 
	
	p(0:3,1) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, 125._omega_prec]
	p(0:3,2) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, -125._omega_prec]
	
	E = p(0,1) + p(0,2)
	Wmboost(0) = ( E + ( ( mWm**2 - mWp**2 ) / E ) )  / 2_omega_prec
	currentWmomentum = sqrt( Wmboost(0)**2 - mWm**2 )
	Wmboost(1:3) = - currentWmomentum * directionW(1:3)
		
	Wpboost(0) = E - Wmboost(0)
	Wpboost(1:3) = -Wmboost(1:3)	
			
	p5(0) = ( mWm + ( ( mass(13)**2 - mass(14)**2 ) / mWm ) )  / 2_omega_prec
	currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
	p5(1:3) = currentWmomentum * directionl(1:3)
			
	p6(0) = mWm - p5(0)
	p6(1:3) = -p5(1:3)
			
	p4(0) = ( mWp + ( ( mass(1)**2 - mass(2)**2 ) / mWp ) )  / 2_omega_prec
	currentWmomentum = sqrt( p4(0)**2 - mass(1)**2 )
	p4(1:3) = currentWmomentum * directionq(1:3)
	
	p3(0) = mWp - p4(0)
	p3(1:3) = -p4(1:3)	
				
	call lorentzBOOST(Wpboost, p3, p(0:3,3) )
	call lorentzBOOST(Wpboost, p4, p(0:3,4) )
	call lorentzBOOST(Wmboost, p5, p(0:3,5) )
	call lorentzBOOST(Wmboost, p6, p(0:3,6) )
	
	zero_ct=0
	n=0
	value = 0
	rho_in(1:4,1) = polarization(1:4)
	call scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
	do i = 1,16
		value = value + rho_out(i,1)
	end do	
end function Binnedelement

function BinnedelementLL (x, data, weights, channel,grids) result (value)
	use kinds
	use omega_kinds
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec), dimension(4) :: polarization
	real(kind=omega_prec) :: value
	interface
		function Binnedelement (x, weights, channel,polarization) result (value)
		use kinds
		use omega_kinds
		use omega_parameters
		use vamp_grid_type !NODEP!
		implicit none 
		real(kind=omega_prec), dimension(:), intent(in) :: x
		real(kind=omega_prec), dimension(:), intent(in) :: polarization
		real(kind=omega_prec), dimension(:), intent(in), optional :: weights
		integer, intent(in), optional :: channel
		real(kind=omega_prec) :: value
		end function Binnedelement
	end interface
	
	polarization(1:4) = [ 1._omega_prec, 0._omega_prec, 0._omega_prec, 0._omega_prec ]
	
	value = Binnedelement(x, weights, channel,polarization)

end function BinnedelementLL

function BinnedelementLR (x, data, weights, channel,grids) result (value)
	use kinds
	use omega_kinds
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec), dimension(4) :: polarization
	real(kind=omega_prec) :: value
	interface
		function Binnedelement (x, weights, channel,polarization) result (value)
		use kinds
		use omega_kinds
		use omega_parameters
		use vamp_grid_type !NODEP!
		implicit none 
		real(kind=omega_prec), dimension(:), intent(in) :: x
		real(kind=omega_prec), dimension(:), intent(in) :: polarization
		real(kind=omega_prec), dimension(:), intent(in), optional :: weights
		integer, intent(in), optional :: channel
		real(kind=omega_prec) :: value
		end function Binnedelement
	end interface
	
	polarization(1:4) = [ 0._omega_prec, 1._omega_prec, 0._omega_prec, 0._omega_prec ]
	
	value = Binnedelement(x, weights, channel,polarization)

end function BinnedelementLR

function BinnedelementRL (x, data, weights, channel,grids) result (value)
	use kinds
	use omega_kinds
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec), dimension(4) :: polarization
	real(kind=omega_prec) :: value
	interface
		function Binnedelement (x, weights, channel,polarization) result (value)
		use kinds
		use omega_kinds
		use omega_parameters
		use vamp_grid_type !NODEP!
		implicit none 
		real(kind=omega_prec), dimension(:), intent(in) :: x
		real(kind=omega_prec), dimension(:), intent(in) :: polarization
		real(kind=omega_prec), dimension(:), intent(in), optional :: weights
		integer, intent(in), optional :: channel
		real(kind=omega_prec) :: value
		end function Binnedelement
	end interface
	
	polarization(1:4) = [ 0._omega_prec, 0._omega_prec, 1._omega_prec, 0._omega_prec ]
	
	value = Binnedelement(x, weights, channel,polarization)

end function BinnedelementRL

function BinnedelementRR (x, data, weights, channel,grids) result (value)
	use kinds
	use omega_kinds
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec), dimension(4) :: polarization
	real(kind=omega_prec) :: value
	interface
		function Binnedelement (x, weights, channel,polarization) result (value)
		use kinds
		use omega_kinds
		use omega_parameters
		use vamp_grid_type !NODEP!
		implicit none 
		real(kind=omega_prec), dimension(:), intent(in) :: x
		real(kind=omega_prec), dimension(:), intent(in) :: polarization
		real(kind=omega_prec), dimension(:), intent(in), optional :: weights
		integer, intent(in), optional :: channel
		real(kind=omega_prec) :: value
		end function Binnedelement
	end interface
	
	polarization(1:4) = [ 0._omega_prec, 0._omega_prec, 0._omega_prec, 1._omega_prec ]
	
	value = Binnedelement(x, weights, channel,polarization)

end function BinnedelementRR



program grid_ww_sl0muq

use omega_kinds
use omega_parameters
use omega95
use MathOperations
use exceptions
use tao_random_numbers
use vamp
use ww_sl0muq
implicit none
integer, dimension(4,1,16,1) :: zero_ct
integer, dimension(2) :: filestat
integer :: n,i,j,h,Ri
!integer :: fi, fo, f, hi, ho,
integer :: ith, ithq, iphq, iph, imWp, imWm
integer :: th_steps, ph_steps, thq_steps, phq_steps, thl_steps, phl_steps, mWp_steps, mWm_steps
integer :: num_args, ix

interface
	function BinnedelementLL (x, data, weights, channel, grids) result (value)
	use kinds
	use omega_kinds
	use omega_parameters
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec) :: value
	end function BinnedelementLL
end interface
interface
	function BinnedelementLR (x, data, weights, channel, grids) result (value)
	use kinds
	use omega_kinds
	use omega_parameters
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec) :: value
	end function BinnedelementLR
end interface
interface
	function BinnedelementRL (x, data, weights, channel, grids) result (value)
	use kinds
	use omega_kinds
	use omega_parameters
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec) :: value
	end function BinnedelementRL
end interface
interface
	function BinnedelementRR (x, data, weights, channel, grids) result (value)
	use kinds
	use omega_kinds
	use omega_parameters
	use vamp !NODEP!
	implicit none 
	real(kind=omega_prec), dimension(:), intent(in) :: x
	class(vamp_data_t), intent(in) :: data
	real(kind=omega_prec), dimension(:), intent(in), optional :: weights
	integer, intent(in), optional :: channel
	type(vamp_grid), dimension(:), intent(in), optional :: grids
	real(kind=omega_prec) :: value
	end function BinnedelementRR
end interface
real(kind=omega_prec), dimension(2) :: weights
integer :: channel
real(kind=omega_prec), dimension(0:3,6) :: p 
real(kind=omega_prec), dimension(4,4) :: polarization
real(kind=omega_prec), dimension(10,8) :: anomalousTGCpar
real(kind=omega_prec), dimension(3,40) :: MatrixElement
real(kind=omega_prec), dimension(0:3) :: p1,p2,p3,p4,p5,p6
!real(kind=omega_prec), dimension(0:3) :: W,Wb
real(kind=omega_prec), dimension(0:3) :: Wpboost, Wmboost
real(kind=omega_prec), dimension(3) :: directionW,directionl,directionq
real(kind=omega_prec), dimension(3) :: A
real(kind=omega_prec), dimension(4) :: pol
real(kind=omega_prec), dimension(8) :: TGCpar,variablex,xwidth
real(kind=omega_prec), dimension(10) :: anomalousg, anomalousk, anomalousl
real(kind=omega_prec), dimension(4,1) :: rho_in
real(kind=omega_prec), dimension(16,4) :: rho_out
real(kind=omega_prec) :: th, ph, thq, phq, thl, phl, mWp, mWm, Wrange, E, currentWmomentum, TGCdev, TGCdevg, TGCdevk, TGCdevl, cw
character(len=12), dimension(:), allocatable :: args
character(len=1024), dimension(2) :: filename
type(exception) :: exc
type(tao_random_state) :: rngLR,rngRL,rng
type(vamp_grid) :: gridLR,gridRL,grid
type(vamp_data_t) :: data
real(kind=default), dimension(2,6) :: domain
real(kind=default) :: integral, error, chi2


TGCdev = 0.001_8
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
num_args = command_argument_count()
allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
do ix = 1, num_args
	call get_command_argument(ix,args(ix))
	!filename(1) = trim(filename(1)) // "_" // trim(args(ix))
	!filename(2) = trim(filename(2)) // "_" // trim(args(ix))
	read (args(ix),*) A(ix)
end do
filename(1) = trim(filename(1)) // ".txt"
filename(2) = trim(filename(2)) // ".txt"

polarization(1:4,1) = [ 1._omega_prec, 0._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,2) = [ 0._omega_prec, 1._omega_prec, 0._omega_prec, 0._omega_prec ]
polarization(1:4,3) = [ 0._omega_prec, 0._omega_prec, 1._omega_prec, 0._omega_prec ]
polarization(1:4,4) = [ 0._omega_prec, 0._omega_prec, 0._omega_prec, 1._omega_prec ]

p(0:3,1) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, 125._omega_prec]
p(0:3,2) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, -125._omega_prec]

th_steps = 1
thl_steps = 1
phl_steps = 1
thq_steps = 1
phq_steps = 1
ph_steps = 1
mWp_steps = 1
mWm_steps = 1
Wrange = 2 * width(24)

open(30,file=trim(filename(1)),iostat=filestat(1))
open(31,file=trim(filename(2) ),iostat=filestat(2) )
if(filestat(1)==0 .AND. filestat(2)==0) then

	th = ( PI / dble(th_steps) ) * A(1)
	thl = ( PI / dble(thl_steps) ) * A(2)
	phl = ( 2_omega_prec * PI / dble(phl_steps) ) * A(3)
	directionl(1) = sin(thl)*cos(phl)
	directionl(2) = sin(thl)*sin(phl)
	directionl(3) = cos(thl)
	!print *, directionl(1:3)
	do ithq=1,1
		thq = ( PI / dble(thq_steps) ) * dble( ithq )
	do iphq=1,1
		phq = ( 2_omega_prec * PI / dble(phq_steps) ) * dble( iphq )
		directionq(1) = sin(thq)*cos(phq)
		directionq(2) = sin(thq)*sin(phq)
		directionq(3) = cos(thq)
		!print *, directionq(1:3)
		MatrixElement(1:3,1:40) = 0
		do iph=1,ph_steps
			ph = ( 2_omega_prec * PI / dble(ph_steps) ) * dble( iph  )
			directionW(1) = sin(th)*cos(ph)
			directionW(2) = sin(th)*sin(ph)
			directionW(3) = cos(th)
			!print *, directionW(1:3)
		do imWp=1,mWp_steps
		mWp = ( ( 2_omega_prec * Wrange / dble(mWp_steps+1) ) * dble( imWp ) ) + mass(24) - Wrange
		do imWm=1,mWm_steps
		mWm = ( ( 2_omega_prec * Wrange / dble(mWm_steps+1) ) * dble( imWm  ) ) + mass(24) - Wrange
			
			!print *, " - " 
			!variablex(1) = th
			!variablex(2) = ph
			!variablex(3) = thq
			!variablex(4) = phq
			!variablex(5) = thl
			!variablex(6) = phl
			!variablex(7) = mWm
			!variablex(8) = mWp
			!xwidth(1) = 0.5_omega_prec * ( PI / dble(th_steps) )
			!xwidth(2) = 0.5_omega_prec * ( 2_omega_prec * PI / dble(ph_steps) ) 
			!xwidth(3) = 0.5_omega_prec * ( PI / dble(thq_steps) )
			!xwidth(4) = 0.5_omega_prec * ( 2_omega_prec * PI / dble(phq_steps) )
			!xwidth(5) = 0.5_omega_prec * ( PI / dble(thl_steps) ) 
			!xwidth(6) = 0.5_omega_prec * ( 2_omega_prec * PI / dble(phl_steps) ) 
			!xwidth(7) = Wrange
			!xwidth(8) = Wrange
						
			variablex(1) = PI / 2_omega_prec
			variablex(2) = PI
			variablex(3) = PI / 2_omega_prec
			variablex(4) = PI
			variablex(5) = PI / 2_omega_prec
			variablex(6) = PI
			xwidth(1) = PI / 2_omega_prec
			xwidth(2) = PI
			xwidth(3) = PI / 2_omega_prec
			xwidth(4) = PI
			xwidth(5) = PI / 2_omega_prec
			xwidth(6) = PI
						
			!variablex(1) = directionW(1)
			!variablex(2) = directionW(2)
			!variablex(3) = directionq(1) 
			!variablex(4) = directionq(2) 
			!variablex(5) = directionl(1)
			!variablex(6) = directionl(2)	
			
					
		
			E = p(0,1) + p(0,2)
			Wmboost(0) = ( E + ( ( mWm**2 - mWp**2 ) / E ) )  / 2_omega_prec
			currentWmomentum = sqrt( Wmboost(0)**2 - mWm**2 )
			Wmboost(1:3) = - currentWmomentum * directionW(1:3)
			
			!print *, Wmboost(0:3)
			
			Wpboost(0) = E - Wmboost(0)
			Wpboost(1:3) = -Wmboost(1:3)
			
			!print *, Wpboost(0:3)		
			
			!print *, " - " 		
			
			p5(0) = ( mWm + ( ( mass(13)**2 - mass(14)**2 ) / mWm ) )  / 2_omega_prec
			currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
			p5(1:3) = currentWmomentum * directionl(1:3)
			
			!print *, p5(0:3)
			
			p6(0) = mWm - p5(0)
			p6(1:3) = -p5(1:3)

			!print *, p6(0:3)
			
			p4(0) = ( mWp + ( ( mass(1)**2 - mass(2)**2 ) / mWp ) )  / 2_omega_prec
			currentWmomentum = sqrt( p4(0)**2 - mass(1)**2 )
			p4(1:3) = currentWmomentum * directionq(1:3)

			!print *, p4(0:3)
			
			p3(0) = mWp - p4(0)
			p3(1:3) = -p4(1:3)
			
			!print *, p3(0:3)
			
			!print *, " - " 		
				
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
			
			h=0
			do Ri=1,10
				TGCpar(1:8) = anomalousTGCpar(Ri,1:8)
				call TGCinit( TGCpar )
				do j=1,4
					h = h + 1
					rho_in(1:4,1) = polarization(1:4,j)
					pol(1:4) = polarization(1:4,j)
					call scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
					 do i = 1,16
					   MatrixElement(1,h) = MatrixElement(1,h) + rho_out(i,1)
					 end do
					!MatrixElement(3,h) = Binnedelement(variablex, weights, channel,pol)
				 end do
			end do
			
			do i=0,9
				print *, MatrixElement(1,(4*i)+1:(4*i)+4)
				!print *, MatrixElement(3,(4*i)+1:(4*i)+4)
				print *, " - "
			end do
			print *, " "
			!print *, MatrixElement(1,5:8)-MatrixElement(1,17:20)
			
			TGCpar(1:8) = anomalousTGCpar(1,1:8)
			call TGCinit( TGCpar )
			
			domain(1,1:6) = variablex(1:6) - xwidth(1:6)
			domain(2,1:6) = variablex(1:6) + xwidth(1:6)
			call tao_random_create (rngLR, seed=0)
			call clear_exception (exc)
			call vamp_create_grid (gridLR, domain, num_calls=100000, exc=exc)
			call handle_exception (exc)
			call clear_exception (exc)
			call vamp_sample_grid (rngLR, gridLR, BinnedelementLR, data, 5,  exc=exc)
			call handle_exception (exc)
			call clear_exception (exc)
			call vamp_discard_integral (gridLR, num_calls=1000000, exc=exc)
			call handle_exception (exc)
			
			!call tao_random_create (rngRL, seed=0)
			!call clear_exception (exc)
			!call vamp_create_grid (gridRL, domain, num_calls=1000, exc=exc)
			!call handle_exception (exc)
			!call clear_exception (exc)
			!call vamp_sample_grid (rngRL, gridRL, BinnedelementRL, data, 3, exc=exc)
			!call handle_exception (exc)
			!call clear_exception (exc)
			!call vamp_discard_integral (gridRL, num_calls=10000, exc=exc)
			!call handle_exception (exc)
			
			
			do Ri=1,5
				TGCpar(1:8) = anomalousTGCpar(Ri,1:8)
				call TGCinit( TGCpar )
				
				!call tao_random_create (rng, seed=0)
				!call clear_exception (exc)
				!call vamp_create_grid (grid, domain, num_calls=100, exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_sample_grid (rng, grid, BinnedelementLL, data, 6,  exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_discard_integral (grid, num_calls=10000, exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_sample_grid (rng, grid, BinnedelementLL, data,4, integral, error, chi2, exc=exc)
				!call handle_exception (exc)
				!print *, "integral = ", integral, "+/-", error, " (chi^2 = ", chi2, ")"
			
				
				call clear_exception (exc)
				call vamp_sample_grid (rngLR, gridLR, BinnedelementLR, data,2, integral, error, chi2, exc=exc)
				call handle_exception (exc)
				print *, "integral = ", integral, "+/-", error, " (chi^2 = ", chi2, ")"

				
				!call clear_exception (exc)
				!call vamp_sample_grid (rngRL, gridRL, BinnedelementRL, data,2, integral, error, chi2, exc=exc)
				!call handle_exception (exc)
				!print *, "integral = ", integral, "+/-", error, " (chi^2 = ", chi2, ")"
			
				!call tao_random_create (rng, seed=0)
				!call clear_exception (exc)
				!call vamp_create_grid (grid, domain, num_calls=1000, exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_sample_grid (rng, grid, BinnedelementRR, data, 6,  exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_discard_integral (grid, num_calls=10000, exc=exc)
				!call handle_exception (exc)
				!call clear_exception (exc)
				!call vamp_sample_grid (rng, grid, BinnedelementRR, data,4, integral, error, chi2, exc=exc)
				!call handle_exception (exc)
				!print *, "integral = ", integral, "+/-", error, " (chi^2 = ", chi2, ")"

				print *, " - "

			end do
			p5(0) = ( mWp + ( ( mass(13)**2 - mass(14)**2 ) / mWp ) )  / 2_omega_prec
			currentWmomentum = sqrt( p5(0)**2 - mass(13)**2 )
			p5(1:3) = currentWmomentum * directionl(1:3)
			
			p6(0) = mWp - p5(0)
			p6(1:3) = -p5(1:3)

			p4(0) = ( mWm + ( ( mass(1)**2 - mass(2)**2 ) / mWm ) )  / 2_omega_prec
			currentWmomentum = sqrt( p4(0)**2 - mass(1)**2 )
			p4(1:3) = currentWmomentum * directionq(1:3)
			
			p3(0) = mWm - p4(0)
			p3(1:3) = -p4(1:3)
									
			call lorentzBOOST(Wmboost, p3, p(0:3,3) )
			call lorentzBOOST(Wmboost, p4, p(0:3,4) )
			call lorentzBOOST(Wpboost, p5, p(0:3,5) )
			call lorentzBOOST(Wpboost, p6, p(0:3,6) )		
			
			zero_ct=0
			n=0
			
			h=0
			do Ri=1,10
				TGCpar(1:8) = anomalousTGCpar(Ri,1:8)
				call TGCinit( TGCpar )
				do j=1,4
					h = h + 1
					rho_in(1:4,1) = polarization(1:4,j)
					call scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
					 do i = 1,16
					   MatrixElement(2,h) = MatrixElement(2,h) + rho_out(i,4)
					 end do
				 end do	
			 end do
			 
			!do i=0,9
			!	print *, MatrixElement(2,(4*i)+1:(4*i)+4)
			!end do
			print *, "End"
			
		end do
		end do
		end do
		
		
		!write(30,*) th, thl, phl, thq, phq, MatrixElement(1,1:40)
		!write(31,*) th, thl, phl, thq, phq, MatrixElement(2,1:40)
	end do	
	end do
else
	write(*,*) "File not opened: Terminating"
end if
close(30)
close(31)
!do ith=1,th_steps
!	th = ( PI / dble(th_steps) ) * dble( ith )
!	print *, th
!end do


end program grid_ww_sl0muq










