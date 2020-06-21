! Copyright 2011-20 Max-Planck-Institut für Eisenforschung GmbH
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!---------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!> @brief general quaternion math, not limited to unit quaternions
!> @details w is the real part, (x, y, z) are the imaginary parts.
!> @details https://en.wikipedia.org/wiki/Quaternion
!---------------------------------------------------------------------------------------------------
module quaternions
  use prec

  implicit none
  private


  real(pReal), parameter, public :: P = -1.0_pReal                                                  !< parameter for orientation conversion.

  type, public :: quaternion
    real(pReal), private :: w = 0.0_pReal
    real(pReal), private :: x = 0.0_pReal
    real(pReal), private :: y = 0.0_pReal
    real(pReal), private :: z = 0.0_pReal


  contains
    procedure, private :: add__
    procedure, private :: pos__
    generic,   public  :: operator(+) => add__,pos__

    procedure, private :: sub__
    procedure, private :: neg__
    generic,   public  :: operator(-) => sub__,neg__

    procedure, private :: mul_quat__
    procedure, private :: mul_scal__
    generic,   public  :: operator(*) => mul_quat__, mul_scal__

    procedure, private :: div_quat__
    procedure, private :: div_scal__
    generic,   public  :: operator(/) => div_quat__, div_scal__

    procedure, private :: eq__
    generic,   public  :: operator(==) => eq__

    procedure, private :: neq__
    generic,   public  :: operator(/=) => neq__

    procedure, private :: pow_quat__
    procedure, private :: pow_scal__
    generic,   public  :: operator(**) => pow_quat__, pow_scal__

    procedure, public  :: abs   => abs__
    procedure, public  :: conjg => conjg__
    procedure, public  :: real  => real__
    procedure, public  :: aimag => aimag__

    procedure, public  :: homomorphed
    procedure, public  :: asArray
    procedure, public  :: inverse

  end type

  interface assignment (=)
    module procedure assign_quat__
    module procedure assign_vec__
  end interface assignment (=)

  interface quaternion
    module procedure init__
  end interface quaternion

  interface abs
    procedure abs__
  end interface abs

  interface dot_product
    procedure dot_product__
  end interface dot_product

  interface conjg
    module procedure conjg__
  end interface conjg

  interface exp
    module procedure exp__
  end interface exp

  interface log
    module procedure log__
  end interface log

  interface real
    module procedure real__
  end interface real

  interface aimag
    module procedure aimag__
  end interface aimag

  public :: &
    assignment(=), &
    conjg, aimag, &
    log, exp, &
    abs, dot_product, &
    inverse, &
    real

contains


!---------------------------------------------------------------------------------------------------
!> @brief construct a quaternion from a 4-vector
!---------------------------------------------------------------------------------------------------
type(quaternion) pure function init__(array)

  real(pReal), intent(in), dimension(4) :: array

  init__%w = array(1)
  init__%x = array(2)
  init__%y = array(3)
  init__%z = array(4)

end function init__


!---------------------------------------------------------------------------------------------------
!> @brief assign a quaternion
!---------------------------------------------------------------------------------------------------
elemental pure subroutine assign_quat__(self,other)

  type(quaternion), intent(out) :: self
  type(quaternion), intent(in)  :: other

  self = [other%w,other%x,other%y,other%z]

end subroutine assign_quat__


!---------------------------------------------------------------------------------------------------
!> @brief assign a 4-vector
!---------------------------------------------------------------------------------------------------
pure subroutine assign_vec__(self,other)

  type(quaternion), intent(out)                :: self
  real(pReal),       intent(in), dimension(4)  :: other

  self%w = other(1)
  self%x = other(2)
  self%y = other(3)
  self%z = other(4)

end subroutine assign_vec__


!---------------------------------------------------------------------------------------------------
!> @brief add a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function add__(self,other)

  class(quaternion), intent(in) :: self,other

  add__ = [ self%w,  self%x,  self%y ,self%z] &
        + [other%w, other%x, other%y,other%z]

end function add__


!---------------------------------------------------------------------------------------------------
!> @brief return (unary positive operator)
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pos__(self)

  class(quaternion), intent(in) :: self

  pos__ = self * (+1.0_pReal)

end function pos__


!---------------------------------------------------------------------------------------------------
!> @brief subtract a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function sub__(self,other)

  class(quaternion), intent(in) :: self,other

  sub__ = [ self%w,  self%x,  self%y ,self%z] &
        - [other%w, other%x, other%y,other%z]

end function sub__


!---------------------------------------------------------------------------------------------------
!> @brief negate (unary negative operator)
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function neg__(self)

  class(quaternion), intent(in) :: self

  neg__ = self * (-1.0_pReal)

end function neg__


!---------------------------------------------------------------------------------------------------
!> @brief multiply with a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function mul_quat__(self,other)

  class(quaternion), intent(in) :: self, other

  mul_quat__%w = self%w*other%w - self%x*other%x -      self%y*other%y - self%z*other%z
  mul_quat__%x = self%w*other%x + self%x*other%w + P * (self%y*other%z - self%z*other%y)
  mul_quat__%y = self%w*other%y + self%y*other%w + P * (self%z*other%x - self%x*other%z)
  mul_quat__%z = self%w*other%z + self%z*other%w + P * (self%x*other%y - self%y*other%x)

end function mul_quat__


!---------------------------------------------------------------------------------------------------
!> @brief multiply with a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function mul_scal__(self,scal)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: scal

  mul_scal__ = [self%w,self%x,self%y,self%z]*scal

end function mul_scal__


!---------------------------------------------------------------------------------------------------
!> @brief divide by a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function div_quat__(self,other)

  class(quaternion), intent(in) :: self, other

  div_quat__ = self * (conjg(other)/(abs(other)**2.0_pReal))

end function div_quat__


!---------------------------------------------------------------------------------------------------
!> @brief divide by a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function div_scal__(self,scal)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: scal

  div_scal__ = [self%w,self%x,self%y,self%z]/scal

end function div_scal__


!---------------------------------------------------------------------------------------------------
!> @brief test equality
!---------------------------------------------------------------------------------------------------
logical elemental pure function eq__(self,other)

  class(quaternion), intent(in) :: self,other

  eq__ = all(dEq([ self%w, self%x, self%y, self%z], &
                 [other%w,other%x,other%y,other%z]))

end function eq__


!---------------------------------------------------------------------------------------------------
!> @brief test inequality
!---------------------------------------------------------------------------------------------------
logical elemental pure function neq__(self,other)

  class(quaternion), intent(in) :: self,other

  neq__ = .not. self%eq__(other)

end function neq__


!---------------------------------------------------------------------------------------------------
!> @brief raise to the power of a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pow_quat__(self,expon)

  class(quaternion), intent(in) :: self
  type(quaternion),  intent(in) :: expon

  pow_quat__ = exp(log(self)*expon)

end function pow_quat__


!---------------------------------------------------------------------------------------------------
!> @brief raise to the power of a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pow_scal__(self,expon)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: expon

  pow_scal__ = exp(log(self)*expon)

end function pow_scal__


!---------------------------------------------------------------------------------------------------
!> @brief take exponential
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function exp__(a)

  class(quaternion), intent(in) :: a
  real(pReal)                   :: absImag

  absImag = norm2(aimag(a))

  exp__ = merge(exp(a%w) * [               cos(absImag), &
                             a%x/absImag * sin(absImag), &
                             a%y/absImag * sin(absImag), &
                             a%z/absImag * sin(absImag)], &
                IEEE_value(1.0_pReal,IEEE_SIGNALING_NAN), &
                dNeq0(absImag))

end function exp__


!---------------------------------------------------------------------------------------------------
!> @brief take logarithm
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function log__(a)

  class(quaternion), intent(in) :: a
  real(pReal)                   :: absImag

  absImag = norm2(aimag(a))

  log__ = merge([log(abs(a)), &
                 a%x/absImag * acos(a%w/abs(a)), &
                 a%y/absImag * acos(a%w/abs(a)), &
                 a%z/absImag * acos(a%w/abs(a))], &
                IEEE_value(1.0_pReal,IEEE_SIGNALING_NAN), &
                dNeq0(absImag))

end function log__


!---------------------------------------------------------------------------------------------------
!> @brief return norm
!---------------------------------------------------------------------------------------------------
real(pReal) elemental pure function abs__(self)

  class(quaternion), intent(in) :: self

  abs__ = norm2([self%w,self%x,self%y,self%z])

end function abs__


!---------------------------------------------------------------------------------------------------
!> @brief calculate dot product
!---------------------------------------------------------------------------------------------------
real(pReal) elemental pure function dot_product__(a,b)

  class(quaternion), intent(in) :: a,b

  dot_product__ = a%w*b%w + a%x*b%x + a%y*b%y + a%z*b%z

end function dot_product__


!---------------------------------------------------------------------------------------------------
!> @brief take conjugate complex
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function conjg__(self)

  class(quaternion), intent(in) :: self

  conjg__ = [self%w,-self%x,-self%y,-self%z]

end function conjg__


!---------------------------------------------------------------------------------------------------
!> @brief homomorph
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function homomorphed(self)

  class(quaternion), intent(in) :: self

  homomorphed = - self

end function homomorphed


!---------------------------------------------------------------------------------------------------
!> @brief return as plain array
!---------------------------------------------------------------------------------------------------
pure function asArray(self)

  real(pReal), dimension(4)     :: asArray
  class(quaternion), intent(in) :: self

  asArray = [self%w,self%x,self%y,self%z]

end function asArray


!---------------------------------------------------------------------------------------------------
!> @brief real part (scalar)
!---------------------------------------------------------------------------------------------------
pure function real__(self)

  real(pReal)                   :: real__
  class(quaternion), intent(in) :: self

  real__ = self%w

end function real__


!---------------------------------------------------------------------------------------------------
!> @brief imaginary part (3-vector)
!---------------------------------------------------------------------------------------------------
pure function aimag__(self)

  real(pReal), dimension(3)     :: aimag__
  class(quaternion), intent(in) :: self

  aimag__ = [self%x,self%y,self%z]

end function aimag__


!---------------------------------------------------------------------------------------------------
!> @brief inverse
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function inverse(self)

  class(quaternion), intent(in) :: self

  inverse = conjg(self)/abs(self)**2.0_pReal

end function inverse


end module quaternions
