subroutine lorentzBOOST( boost, MyVector, BoostedVector) ! Ein Unterprogramm des Typs subroutine	
	! Deklaration der formalen Parameter
	implicit none 
	real(kind=8), dimension(0:3), intent(in)  :: boost, MyVector
	real(kind=8), dimension(0:3), intent(out) :: BoostedVector
	real(kind=8), dimension(0:3,0:3) :: LorentzMatrix
	real(kind=8), dimension(1:3) :: direction
	real(kind=8) :: lg,gb,currentp
	integer :: i,j
	currentp=0
	do i=1,3
		currentp=currentp + boost(i)**2._8
	end do
	lg = boost(0) / sqrt( boost(0)**2._8 - currentp )
	currentp=sqrt( currentp )
	gb = sqrt( lg**2 - 1 )
	print *,'Momentum, Gamma, gammabeta ', currentp, lg, gb
	do i=1,3
		direction(i) = boost(i) / currentp
	end do
	print *,'Direction ', direction
	LorentzMatrix=0
	LorentzMatrix(0,0)=lg
	do i=1,3
		LorentzMatrix(0,i) = -1._8*gb*direction(i)
		LorentzMatrix(i,0) = LorentzMatrix(0,i)
		do j=1,3
			LorentzMatrix(j,i) = ( lg - 1._8) * direction(i) * direction(j)
		end do
		LorentzMatrix(i,i) = LorentzMatrix(i,i) + 1._8
	end do
	BoostedVector=0
	! BoostedVector = boost
	print *, 'Begin Matrix: '
	do i=0,3
	print *, LorentzMatrix(0:3,i)
	do j=0,3
		BoostedVector(i) = BoostedVector(i) + LorentzMatrix(i,j) * MyVector(j)
	end do
	end do
	print *, 'End Matrix'
	! return   
end subroutine lorentzBOOST



program testww_sl0muq
! compile and link with:
!   /usr/bin/gfortran -I. -I../omega-src/bundle/lib/ testww_sl0muq.F90 -L. -lproc -L../omega-src/bundle/lib/ -l omega95
! in processes-src
 
use omega_kinds
use omega_parameters
use omega95
use ww_sl0muq
!  in the module:
!  public :: scatter_nonzero,  scatter_diagonal_nonzero, &
!    scatter_diagonal_colored_nz, scatter_colored_nonzero
!  ie. callable form the outside

implicit none
integer, dimension(4,1,16,1) :: zero_ct
integer :: n,i, fi, fo, f, hi, ho, h, filestat

! 0:3 = four-vector (probablY)
! 8 = total number of in- and out- particle 
real(kind=omega_prec), dimension(0:3,6) :: p 
real(kind=omega_prec), dimension(0:3) :: p1 
real(kind=omega_prec), dimension(0:3) :: p2 
real(kind=omega_prec), dimension(0:3) :: p3 
real(kind=omega_prec), dimension(0:3) :: p4
real(kind=omega_prec), dimension(0:3) :: Wp
real(kind=omega_prec), dimension(0:3) :: Wpboost
real(kind=omega_prec) :: th, ph, thq, phq, thl, phl
type(spinor) :: psi
type(momentum) :: inputmomentum
! 4 = total number of initial helicities
! 16 = ditto, final

real(kind=omega_prec), dimension(4,1) :: rho_in
real(kind=omega_prec), dimension(16,4) :: rho_out

open(30,file="test.txt",iostat=filestat)
if(filestat==0) then
	write(30,*) "Test"
n=6
call setup_parameters()
!call print_parameters ()
p(0:3,1) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, 125._omega_prec]
p(0:3,2) = [ 125._omega_prec , 0._omega_prec, 0._omega_prec, -125._omega_prec]
th = PI / 4.
ph = 0.
thq = PI / 4.
phq = 0.
thl = PI / 4.
phl = 0.
print *, mass(24)
!print *, p1
Wp(0) = 125
Wp(1) = sqrt( (Wp(0)**2._omega_prec) - (mass(24)**2._omega_prec) )*sin(th)*cos(ph)
Wp(2) = sqrt( (Wp(0)**2._omega_prec) - (mass(24)**2._omega_prec) )*sin(th)*sin(ph)
Wp(3) = sqrt( (Wp(0)**2._omega_prec) - (mass(24)**2._omega_prec) )*cos(th)
p1(0) = mass(24)*0.5_omega_prec
p1(1) = p1(0) * sin(th)*cos(ph+PI/2._omega_prec)
p1(2) = p1(0) * sin(th)*sin(ph+PI/2._omega_prec)
p1(3) = p1(0) * cos(th)
! p1 = [ mass(24)*0.5_omega_prec, mass(24)*0.5_omega_prec*sin(th)*cos(ph+PI/2._omega_prec), mass(24)*0.5_omega_prec*sin(th)*sin(ph+PI/2._omega_prec), mass(24)*0.5_omega_prec*cos(th) ]
p2 = p1
do i=1,3
	p2(i) = -p2(i)
end do
!print *, Wp
Wpboost=Wp
call lorentzBOOST(Wpboost, p2, p(0:3,3) )
call lorentzBOOST(Wpboost, p1, p(0:3,4) )
do i=1,3
	Wpboost(i) = -Wpboost(i)
end do
call lorentzBOOST(Wpboost, p1, p(0:3,5) )
call lorentzBOOST(Wpboost, p2, p(0:3,6) )
rho_in = 0
!print *, rho_in(1:4,1)
! rho_in(1:4,1) = [ 0.25_omega_prec, 0.25_omega_prec, 0.25_omega_prec, 0.25_omega_prec ]
! rho_in(1:4,1) = [ 1._omega_prec, 0._omega_prec, 0._omega_prec, 0._omega_prec ]
! rho_in(1:4,1) = [ 0._omega_prec, 1._omega_prec, 0._omega_prec, 0._omega_prec ]
! rho_in(1:4,1) = [ 0._omega_prec, 0._omega_prec, 1._omega_prec, 0._omega_prec ]
! rho_in(1:4,1) = [ 0._omega_prec, 0._omega_prec, 0._omega_prec, 1._omega_prec ]
 rho_in(1:4,1) = [ 0._omega_prec, 0.5_omega_prec, 0.5_omega_prec, 0._omega_prec ]
print *, rho_in(1:4,1)
do hi = 1, 4
do ho = 1, 16
	zero_ct(hi,1,ho,1) = 0
end do
end do
inputmomentum%t = p1(0)
do i=1,3
	inputmomentum%x(i) = p1(i)
end do
psi = u(1._omega_prec,inputmomentum,-1)
!print *, psi
print *, 'Momentums: '
do i=1,6
	print *, p(0:3,i)
end do
 call scatter_diagonal_nonzero (p, rho_in, rho_out, zero_ct, n)
 print *, 'Matrix elements:'
 do i = 1,16
   print *, rho_out(i,1:4)
 end do
!do i=1,16
!	print *, i, zero_ct(1:4,1,i,1)
!end do
else
	write(*,*) "File not opened: Terminating"
end if
close(30)
end program testww_sl0muq
