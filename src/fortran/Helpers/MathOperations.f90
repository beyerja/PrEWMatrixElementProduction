module MathOperations
  public :: lorentzBOOST,SphereCoordinate
contains  
subroutine lorentzBOOST( boost, MyVector, BoostedVector) ! Ein Unterprogramm des Typs subroutine	
  ! Other modules
  use omega_kinds ! Only for O'Mega precision
	! Deklaration der formalen Parameter
  implicit none 
	real(kind=omega_prec), dimension(0:3), intent(in)  :: boost, MyVector
	real(kind=omega_prec), dimension(0:3), intent(out) :: BoostedVector
	real(kind=omega_prec), dimension(0:3,0:3) :: LorentzMatrix
	real(kind=omega_prec), dimension(1:3) :: direction
	real(kind=omega_prec) :: lg,gb,currentp
	integer :: i,j
	currentp=0
	do i=1,3
		currentp=currentp + boost(i)**2._omega_prec
	end do
	lg = boost(0) / sqrt( boost(0)**2._omega_prec - currentp )
	currentp=sqrt( currentp )
	gb = sqrt( lg**2 - 1 )
	! print *,'Momentum, Gamma, gammabeta ', currentp, lg, gb
	do i=1,3
		direction(i) = boost(i) / currentp
	end do
	! print *,'Direction ', direction
	LorentzMatrix=0
	LorentzMatrix(0,0)=lg
	do i=1,3
		LorentzMatrix(0,i) = -1._omega_prec*gb*direction(i)
		LorentzMatrix(i,0) = LorentzMatrix(0,i)
		do j=1,3
			LorentzMatrix(j,i) = ( lg - 1._omega_prec) * direction(i) * direction(j)
		end do
		LorentzMatrix(i,i) = LorentzMatrix(i,i) + 1._omega_prec
	end do
	BoostedVector=0
	! BoostedVector = boost
	! print *, 'Begin Matrix: '
	do i=0,3
	! print *, LorentzMatrix(0:3,i)
	do j=0,3
		BoostedVector(i) = BoostedVector(i) + LorentzMatrix(i,j) * MyVector(j)
	end do
	end do
	! print *, 'End Matrix'
	! return   
end subroutine lorentzBOOST

pure function SphereCoordinate( magnitude, theta, phi ) result (Cartesian)
  ! Other modules
  use omega_kinds ! Only for O'Mega precision
	real(kind=omega_prec), intent(in)  :: magnitude, theta, phi
	real(kind=omega_prec), dimension(3) :: Cartesian
	Cartesian=0
	Cartesian(1) = magnitude*sin(theta)*cos(phi)
	Cartesian(2) = magnitude*sin(theta)*sin(phi)
	Cartesian(3) = magnitude*cos(theta)
	
end function SphereCoordinate
	
end module MathOperations















