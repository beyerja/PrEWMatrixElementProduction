module Rotation

implicit none

contains

function cross(a, b)
  ! Function calculates and returns the cross product of two vectors
  
  ! Other modules
  use omega_kinds ! Only for O'Mega precision
  
  ! Variable definitions
  real(kind=omega_prec), dimension(3) :: cross
  real(kind=omega_prec), dimension(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

function rotate_out_of(x_in,z_in) result(rotation)
  ! Function calculates the rotation matrix.
  ! z point already in the right direction, x may need to be adjusted.
  ! Rotation matrix can be used to rotate _out_ of the given coordinate system.
  
  ! Other modules
  use omega_kinds ! Only for O'Mega precision
  
  ! Variable definitions
  real(kind=omega_prec), dimension (3,3) :: rotation
  real(kind=omega_prec), dimension(9) :: rotvec
  real(kind=omega_prec), dimension(3), INTENT(IN) :: x_in, z_in
  real(kind=omega_prec), dimension(3) :: x, y, z

  ! Normalise the z axis
  z = 1._omega_prec/norm2(z_in) * z_in
  
  ! Get the y axis -> Cross product of z and x
  y = cross(z,x_in)
  y = 1._omega_prec/norm2(y) * y
  
  ! Now determine correct x axis -> y cross z
  x = cross(y,z)
  x = 1._omega_prec/norm2(x) * x

  ! Build rotation matrix from that
  rotvec(1:3) = x
  rotvec(4:6) = y
  rotvec(7:9) = z
  rotation = reshape(rotvec, [3,3])
end function rotate_out_of

end module Rotation