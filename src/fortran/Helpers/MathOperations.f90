module MathOperations
  public :: lorentzBOOST,SphereCoordinate
contains  
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
	! print *,'Momentum, Gamma, gammabeta ', currentp, lg, gb
	do i=1,3
		direction(i) = boost(i) / currentp
	end do
	! print *,'Direction ', direction
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
	real(kind=8), intent(in)  :: magnitude, theta, phi
	real(kind=8), dimension(3) :: Cartesian
	Cartesian=0
	Cartesian(1) = magnitude*sin(theta)*cos(phi)
	Cartesian(2) = magnitude*sin(theta)*sin(phi)
	Cartesian(3) = magnitude*cos(theta)
	
end function SphereCoordinate
	
end module MathOperations















