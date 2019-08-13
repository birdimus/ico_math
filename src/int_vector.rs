// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::float_vector::FloatVector;
use core::arch::x86_64::*;
use core::hash::Hash;
use core::hash::Hasher;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct IntVector {
    pub data: __m128i,
}

impl IntVector {
    /// Returns a new Vector2
    #[inline(always)]
    pub fn new(value: i32) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_set1_epi32(value),
            }
        }
    }
    #[inline(always)]
    pub fn zero() -> IntVector {
        unsafe {
            IntVector {
                data: _mm_setzero_si128(),
            }
        }
    }
    #[inline(always)]
    pub fn value(&self) -> i32 {
        unsafe {
            return _mm_cvtsi128_si32(self.data);
        }
    }

    #[inline(always)]
    pub fn max(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_max_epi32(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn min(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_min_epi32(v1.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn abs(v1: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_abs_epi32(v1.data),
            }
        }
    }
    #[inline(always)]
    pub fn sign(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_sign_epi32(v1.data, v2.data),
            }
        }
    }

    /*
    // TODO: Requires the shifts to be constant - waiting on const generics, or const arguments
    // https://github.com/rust-lang/rfcs/pull/2000
    // could also switch on the shift var. 0-32
    #[inline(always)]
    pub fn shift_right(v1 : IntVector, shift : i32) -> IntVector{
        unsafe{
            IntVector{data : _mm_srli_epi32(v1.data, shift)}
        }
    }

    #[inline(always)]
    pub fn shift_left(v1 : IntVector, shift : i32) -> IntVector{
        unsafe{
            IntVector{data : _mm_slli_epi32(v1.data, shift)}
        }
    }
    */
    #[inline(always)]
    pub fn add(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_add_epi32(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn sub(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_sub_epi32(v1.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn mul(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_mullo_epi32(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn and(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_and_si128(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn or(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_or_si128(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn andnot(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_andnot_si128(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn xor(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_xor_si128(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn equal(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cmpeq_epi32(v1.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn greater(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cmpgt_epi32(v1.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn less(v1: IntVector, v2: IntVector) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cmplt_epi32(v1.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn equals(v1: IntVector, v2: IntVector) -> bool {
        return v1.value() == v2.value();
    }
}

impl From<IntVector> for i32 {
    #[inline(always)]
    fn from(v: IntVector) -> i32 {
        return v.value();
    }
}
impl From<i32> for IntVector {
    #[inline(always)]
    fn from(v: i32) -> IntVector {
        return IntVector::new(v);
    }
}
impl From<FloatVector> for IntVector {
    #[inline(always)]
    fn from(v: FloatVector) -> IntVector {
        unsafe {
            return IntVector {
                data: _mm_cvttps_epi32(v.data),
            };
        }
    }
}

impl<T: Into<IntVector>> core::ops::Add<T> for IntVector {
    type Output = IntVector;
    #[inline(always)]
    fn add(self, _rhs: T) -> IntVector {
        IntVector::add(self, _rhs.into())
    }
}
impl<T: Into<IntVector>> core::ops::AddAssign<T> for IntVector {
    #[inline(always)]
    fn add_assign(&mut self, other: T) {
        *self = IntVector::add(*self, other.into())
    }
}
impl<T: Into<IntVector>> core::ops::Sub<T> for IntVector {
    type Output = IntVector;
    #[inline(always)]
    fn sub(self, _rhs: T) -> IntVector {
        IntVector::sub(self, _rhs.into())
    }
}

impl<T: Into<IntVector>> core::ops::SubAssign<T> for IntVector {
    #[inline(always)]
    fn sub_assign(&mut self, other: T) {
        *self = IntVector::sub(*self, other.into())
    }
}
impl core::ops::Neg for IntVector {
    type Output = IntVector;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return IntVector {
                data: _mm_sub_epi32(_mm_set1_epi32(0), self.data),
            };
        }
    }
}

impl<T: Into<IntVector>> core::ops::Mul<T> for IntVector {
    type Output = IntVector;
    #[inline(always)]
    fn mul(self, _rhs: T) -> IntVector {
        IntVector::mul(self, _rhs.into())
    }
}
impl<T: Into<IntVector>> core::ops::MulAssign<T> for IntVector {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = IntVector::mul(*self, _rhs.into())
    }
}

impl PartialEq for IntVector {
    #[inline(always)]
    fn eq(&self, other: &IntVector) -> bool {
        return IntVector::equals(*self, *other);
    }
}
impl PartialEq<i32> for IntVector {
    #[inline(always)]
    fn eq(&self, other: &i32) -> bool {
        return self.value() == *other;
    }
}
impl PartialEq<IntVector> for i32 {
    #[inline(always)]
    fn eq(&self, other: &IntVector) -> bool {
        return *self == other.value();
    }
}
impl Hash for IntVector {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value().hash(state);
    }
}
