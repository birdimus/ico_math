use crate::int_vector::IntVector;
use crate::sse_extensions::*;
use crate::structure::SIMDVector1;
use crate::structure::SIMDVector2;
use crate::structure::SIMDVector3;
use crate::structure::SIMDVector4;
use core::arch::x86_64::*;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct FloatVector {
    pub data: __m128,
}

impl FloatVector {
    /// Returns a new FloatVector
    #[inline(always)]
    pub fn new(value: f32) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_set1_ps(value),
            }
        }
    }
    #[inline(always)]
    pub fn zero() -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_setzero_ps(),
            }
        }
    }

    #[inline(always)]
    pub fn value(self) -> f32 {
        unsafe {
            return _mm_cvtss_f32(self.data);
        }
    }

    /// Compute the sum of two vectors and return it as a new vector.
    #[inline(always)]
    pub fn add(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_add_ps(self.data, v2.data),
            }
        }
    }

    /// Subtract a vector and return it as a new vector.
    #[inline(always)]
    pub fn sub(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_sub_ps(self.data, v2.data),
            }
        }
    }

    /// Multiply two vectors component-wise, and return it as a new vector.  
    /// (x1 * x2, y1 * y2, z1 * z2, w1 * w2)
    #[inline(always)]
    pub fn mul(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_mul_ps(self.data, v2.data),
            }
        }
    }

    /// Divide two vectors component-wise, and return it as a new vector.  
    /// (x1 / x2, y1 / y2, z1 / z2, w1 / w2)
    #[inline(always)]
    pub fn div(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_div_ps(self.data, v2.data),
            }
        }
    }

    /// Fused Multiply Add.  Result = (a * b) + c.
    #[inline(always)]
    pub fn mul_add(self, v2: FloatVector, v3: FloatVector) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _mm_fmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Fused Multiply Sub.  Result = (a * b) - c.
    #[inline(always)]
    pub fn mul_sub(self, v2: FloatVector, v3: FloatVector) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _mm_fmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Add.  Result = -(a * b) + c.
    #[inline(always)]
    pub fn neg_mul_add(self, v2: FloatVector, v3: FloatVector) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _mm_fnmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Sub.  Result = -(a * b) - c.
    #[inline(always)]
    pub fn neg_mul_sub(self, v2: FloatVector, v3: FloatVector) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _mm_fnmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Compute the bitwise AND of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn and<T: SIMDVector1>(self, v2: T) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_and_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise OR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn or<T: SIMDVector1>(self, v2: T) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_or_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise ANDNOT of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn andnot<T: SIMDVector1>(self, v2: T) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_andnot_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise XOR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn xor<T: SIMDVector1>(self, v2: T) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_xor_ps(self.data, v2.data()),
            }
        }
    }

    /// Equals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn equal(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmpeq_ps(self.data, v2.data))) != 0;
        }
    }
    /// NotEquals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn not_equal(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmpneq_ps(self.data, v2.data))) != 0;
        }
    }
    /// Greater than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater_equal(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmpge_ps(self.data, v2.data))) != 0;
        }
    }
    /// Greater than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmpgt_ps(self.data, v2.data))) != 0;
        }
    }
    /// Less than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less_equal(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmple_ps(self.data, v2.data))) != 0;
        }
    }
    /// Less than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less(self, v2: FloatVector) -> bool {
        unsafe {
            return _mm_movemask_epi8(_mm_castps_si128(_mm_cmplt_ps(self.data, v2.data))) != 0;
        }
    }
    /// Relative and absolute epsilon comparison.  
    /// Uses machine epsilon as absolute, and 4*machine epsilon for relative.
    /// return abs(a - b) <= max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    /// Adapted from Knuth.  
    /// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    #[inline(always)]
    pub fn approx_equal(self, v2: FloatVector) -> bool {
        let delta = (self - v2).abs();
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(FloatVector::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.less_equal(epsilon_bound);
    }

    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_greater(self, v2: FloatVector) -> bool {
        let delta = self.sub(v2);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(FloatVector::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }
    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_less(self, v2: FloatVector) -> bool {
        let delta = v2.sub(self);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(FloatVector::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }
    /// The absolute value of each component of the vector.
    #[inline(always)]
    pub fn abs(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _ico_abs_ps(self.data),
            }
        }
    }

    /// Take the magnitude of the first argument (self), and use the sign of the second argument to produce a new vector
    #[inline(always)]
    pub fn copysign(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _ico_copysign_ps(self.data, v2.data),
            }
        }
    }

    /// Floor function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn floor(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_floor_ps(self.data),
            }
        }
    }

    /// Ceil function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn ceil(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_ceil_ps(self.data),
            }
        }
    }

    /// Round to nearest even function. Returns signed 0 when applicable.
    #[inline(always)]
    pub fn round(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Truncate function.
    #[inline(always)]
    pub fn truncate(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Convert to an integer vector using the floor function.
    #[inline(always)]
    pub fn floor_to_int(self) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cvttps_epi32(self.floor().data),
            }
        }
    }

    /// Convert to an integer vector using the ceil function.
    #[inline(always)]
    pub fn ceil_to_int(self) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cvttps_epi32(self.ceil().data),
            }
        }
    }
    /// Convert to an integer vector using the round function.
    #[inline(always)]
    pub fn round_to_int(self) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cvttps_epi32(self.round().data),
            }
        }
    }

    /// Convert to an integer vector using the truncate function.
    #[inline(always)]
    pub fn truncate_to_int(self) -> IntVector {
        unsafe {
            IntVector {
                data: _mm_cvttps_epi32(self.truncate().data),
            }
        }
    }

    /// Compute the fractional component of each component
    /// Result = X - Floor(x)
    #[inline(always)]
    pub fn frac(self) -> FloatVector {
        return FloatVector::sub(self, FloatVector::floor(self));
    }

    /// Compute the square root of each component
    #[inline(always)]
    pub fn sqrt(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_sqrt_ps(self.data),
            }
        }
    }

    /// Compute the approximate sin of each component
    #[inline(always)]
    pub fn sin(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _ico_sin_ps(self.data),
            }
        }
    }

    /// Compute the approximate cos of each component
    #[inline(always)]
    pub fn cos(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _ico_cos_ps(self.data),
            }
        }
    }

    /// Compute the approximate acos of each component
    #[inline(always)]
    pub fn acos(self) -> FloatVector {
        unsafe {
            FloatVector {
                data: _ico_acos_ps(self.data),
            }
        }
    }

    /// Compute the component-wise max.
    #[inline(always)]
    pub fn max(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_max_ps(self.data, v2.data),
            }
        }
    }

    /// Compute the component-wise min.
    #[inline(always)]
    pub fn min(self, v2: FloatVector) -> FloatVector {
        unsafe {
            FloatVector {
                data: _mm_min_ps(self.data, v2.data),
            }
        }
    }
    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: FloatVector, mask: bool) -> FloatVector {
        let mask_val: u32 = if mask { 0xFFFFFFFF } else { 0 };
        unsafe {
            FloatVector {
                data: _ico_select_ps(
                    self.data,
                    v2.data,
                    _mm_castsi128_ps(_mm_set1_epi32(core::mem::transmute::<u32, i32>(mask_val))),
                ),
            }
        }
    }
    /// Linear interpolate from a to b based on a float.
    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(self, v2: FloatVector, t: T) -> FloatVector {
        unsafe {
            let t_val = t.into().data;
            let tmp = _mm_fnmadd_ps(self.data, t_val, self.data); //a - (a*t)
            FloatVector {
                data: _mm_fmadd_ps(v2.data, t_val, tmp),
            } //b * t + a
        }
    }
}

impl From<FloatVector> for f32 {
    #[inline(always)]
    fn from(v: FloatVector) -> f32 {
        return v.value();
    }
}
impl From<f32> for FloatVector {
    #[inline(always)]
    fn from(v: f32) -> FloatVector {
        return FloatVector::new(v);
    }
}
impl From<IntVector> for FloatVector {
    #[inline(always)]
    fn from(v: IntVector) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _mm_cvtepi32_ps(v.data),
            };
        }
    }
}

impl<T: Into<FloatVector>> core::ops::Add<T> for FloatVector {
    type Output = FloatVector;
    #[inline(always)]
    fn add(self, _rhs: T) -> FloatVector {
        FloatVector::add(self, _rhs.into())
    }
}
impl<T: Into<FloatVector>> core::ops::AddAssign<T> for FloatVector {
    #[inline(always)]
    fn add_assign(&mut self, other: T) {
        *self = FloatVector::add(*self, other.into())
    }
}
impl<T: Into<FloatVector>> core::ops::Sub<T> for FloatVector {
    type Output = FloatVector;
    #[inline(always)]
    fn sub(self, _rhs: T) -> FloatVector {
        FloatVector::sub(self, _rhs.into())
    }
}

impl<T: Into<FloatVector>> core::ops::SubAssign<T> for FloatVector {
    #[inline(always)]
    fn sub_assign(&mut self, other: T) {
        *self = FloatVector::sub(*self, other.into())
    }
}
impl core::ops::Neg for FloatVector {
    type Output = FloatVector;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return FloatVector {
                data: _mm_xor_ps(_ico_signbit_ps(), self.data),
            };
        }
    }
}

impl<T: Into<FloatVector>> core::ops::Mul<T> for FloatVector {
    type Output = FloatVector;
    #[inline(always)]
    fn mul(self, _rhs: T) -> FloatVector {
        FloatVector::mul(self, _rhs.into())
    }
}
impl<T: Into<FloatVector>> core::ops::MulAssign<T> for FloatVector {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = FloatVector::mul(*self, _rhs.into())
    }
}

impl<T: Into<FloatVector>> core::ops::Div<T> for FloatVector {
    type Output = FloatVector;
    #[inline(always)]
    fn div(self, _rhs: T) -> FloatVector {
        FloatVector::div(self, _rhs.into())
    }
}
impl<T: Into<FloatVector>> core::ops::DivAssign<T> for FloatVector {
    #[inline(always)]
    fn div_assign(&mut self, _rhs: T) {
        *self = FloatVector::div(*self, _rhs.into())
    }
}
impl PartialEq for FloatVector {
    #[inline(always)]
    fn eq(&self, other: &FloatVector) -> bool {
        return FloatVector::equal(*self, *other);
    }
}
impl PartialEq<f32> for FloatVector {
    #[inline(always)]
    fn eq(&self, other: &f32) -> bool {
        return self.value() == *other;
    }
}
impl PartialEq<FloatVector> for f32 {
    #[inline(always)]
    fn eq(&self, other: &FloatVector) -> bool {
        return *self == other.value();
    }
}
impl SIMDVector1 for FloatVector {
    fn data(self) -> __m128 {
        return self.data;
    }
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}
impl SIMDVector2 for FloatVector {
    fn data(self) -> __m128 {
        return self.data;
    }
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}
impl SIMDVector3 for FloatVector {
    fn data(self) -> __m128 {
        return self.data;
    }
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}
impl SIMDVector4 for FloatVector {
    fn data(self) -> __m128 {
        return self.data;
    }
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}
// #[cfg(test)]
// mod test;
