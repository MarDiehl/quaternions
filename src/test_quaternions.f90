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
!--------------------------------------------------------------------------------------------------
!> @brief check correctness of some quaternions functions
!--------------------------------------------------------------------------------------------------
program test_quaternions
  use prec
  use quaternions

  implicit none
  real(pReal), dimension(4) :: qu
  type(quaternion)          :: q, q_2

  call random_number(qu)
  qu = (qu-0.5_pReal) * 2.0_pReal
  q  = quaternion(qu)

  q_2= qu
  if(any(dNeq(q%asArray(),q_2%asArray())))             error stop 'assign_vec__'

  q_2 = q + q
  if(any(dNeq(q_2%asArray(),2.0_pReal*qu)))            error stop 'add__'

  q_2 = q - q
  if(any(dNeq0(q_2%asArray())))                        error stop 'sub__'

  q_2 = q * 5.0_pReal
  if(any(dNeq(q_2%asArray(),5.0_pReal*qu)))            error stop 'mul__'

  q_2 = q / 0.5_pReal
  if(any(dNeq(q_2%asArray(),2.0_pReal*qu)))            error stop 'div__'

  q_2 = q * 0.3_pReal
  if(dNeq0(abs(q)) .and. q_2 == q)                     error stop 'eq__'

  q_2 = q
  if(q_2 /= q)                                         error stop 'neq__'

  if(dNeq(abs(q),norm2(qu)))                           error stop 'abs__'
  if(dNeq(abs(q)**2.0_pReal, real(q*q%conjg()),1.0e-14_pReal)) &
                                                       error stop 'abs__/*conjg'

  if(any(dNeq(q%asArray(),qu)))                        error stop 'eq__'
  if(dNeq(q%real(),       qu(1)))                      error stop 'real()'
  if(any(dNeq(q%aimag(),  qu(2:4))))                   error stop 'aimag()'

  q_2 = q%homomorphed()
  if(q                 /= q_2*    (-1.0_pReal))        error stop 'homomorphed'
  if(dNeq(q_2%real(),     qu(1)*  (-1.0_pReal)))       error stop 'homomorphed/real'
  if(any(dNeq(q_2%aimag(),qu(2:4)*(-1.0_pReal))))      error stop 'homomorphed/aimag'

  q_2 = conjg(q)
  if(dNeq(abs(q),abs(q_2)))                            error stop 'conjg/abs'
  if(q /= conjg(q_2))                                  error stop 'conjg/involution'
  if(dNeq(q_2%real(),     q%real()))                   error stop 'conjg/real'
  if(any(dNeq(q_2%aimag(),q%aimag()*(-1.0_pReal))))    error stop 'conjg/aimag'

  if(abs(q) > 0.0_pReal) then
    q_2 = q * q%inverse()
    if(     dNeq(real(q_2), 1.0_pReal,1.0e-15_pReal))  error stop 'inverse/real'
    if(any(dNeq0(aimag(q_2),          1.0e-15_pReal))) error stop 'inverse/aimag'

    q_2 = q/abs(q)
    q_2 = conjg(q_2) - inverse(q_2)
    if(any(dNeq0(q_2%asArray(),1.0e-15_pReal)))        error stop 'inverse/conjg'
  endif
  if(dNeq(dot_product(qu,qu),dot_product(q,q)))        error stop 'dot_product'

#if !(defined(__GFORTRAN__) &&  __GNUC__ < 9)
  if (norm2(aimag(q)) > 0.0_pReal) then
    if (dNeq0(abs(q-exp(log(q))),1.0e-13_pReal))       error stop 'exp/log'
    if (dNeq0(abs(q-log(exp(q))),1.0e-13_pReal))       error stop 'log/exp'
  endif
#endif
  print*, 'All fine'
end program test_quaternions
