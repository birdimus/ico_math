use crate::float_vector::FloatVector;
use crate::sse_extensions::*;
use crate::structure::SIMDVector2;
use crate::vector2_bool::Vector2Bool;
use crate::vector2_int::Vector2Int;
use crate::vector3::Vector3;
use crate::vector4::Vector4;
use core::arch::x86_64::*;

/// A vector of 2 floats (x,y).
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2 {
    pub data: __m128,
}

impl Vector2 {
    /// Construct a new vector from f32 components.
    #[inline(always)]
    pub fn new(x: f32, y: f32) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_set_ps(0.0f32, 0.0f32, y, x),
            }
        }
    }

    /// Set all values of the vector to the same f32 value.
    #[inline(always)]
    pub fn set<T: Into<FloatVector>>(value: T) -> Vector2 {
        return Vector2 {
            data: value.into().data,
        };
    }

    /// Construct a new vector of zeros.
    #[inline(always)]
    pub fn zero() -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_setzero_ps(),
            }
        }
    }

    /// Get the x value of the vector, broadcast to all components as a FloatVector (xxxx).
    #[inline(always)]
    pub fn x(self) -> FloatVector {
        return FloatVector {
            data: self.xxxx().data,
        };
    }

    /// Get the y value of the vector, broadcast to all components as a FloatVector (yyyy).
    #[inline(always)]
    pub fn y(self) -> FloatVector {
        return FloatVector {
            data: self.yyyy().data,
        };
    }

    /// Set the x value of this vector, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_x<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            self.data = _mm_move_ss(self.data, value.into().data);
        }
    }

    /// Set the y value of this vector, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_y<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(value.into().data, self.data);
            self.data = _mm_shuffle_ps(v1, self.data, _ico_shuffle(3, 2, 1, 0));
        }
    }

    /// Compute the 2 element dot-product, and return it as a broadcast FloatVector.
    #[inline(always)]
    pub fn dot(v0: Vector2, v1: Vector2) -> FloatVector {
        unsafe {
            let tmp0 = _mm_mul_ps(v0.data, v1.data);
            let mut tmp1 = _yxzw(tmp0); // _mm_shuffle_ps(tmp0,tmp0, _ico_shuffle(3, 2, 0, 1)); //yxzw

            tmp1 = _mm_add_ps(tmp0, tmp1); //xy,xy,qq,qq
            return FloatVector {
                data: _mm_unpacklo_ps(tmp1, tmp1),
            }; //xy,xy.xy,xy
        }
    }

    /// Rotate the vector by radians
    /// Right handed system, positive rotation is counterclockwise about the axis of rotation.
    #[inline(always)]
    pub fn rotate(v1: Vector2, radians: f64) -> Vector2 {
        let f = radians.sin_cos();

        let sn: f32 = f.0 as f32;
        let cs: f32 = f.1 as f32;

        unsafe {
            // Any values below the epsilon get clamped to zero.  This fixes precision issues around zero.
            let epsilon = _mm_set1_ps(EPSILON_AT_ONE);
            let sncs = _mm_set_ps(-sn, sn, cs, cs);
            let mask = _mm_cmpgt_ps(_ico_abs_ps(sncs), epsilon);

            let masked_sncs = _mm_and_ps(sncs, mask);

            //x1y1x2y2
            let xyxy = _mm_movelh_ps(v1.data, v1.data);
            //x * cs, y * cs, x*sn, y*-sn
            let v2 = _mm_mul_ps(xyxy, masked_sncs);

            ///x1 + y2, y1 + x2,
            //x = x * cs - y * sn;
            //y = x * sn + y * cs;
            return Vector2 {
                data: _mm_add_ps(v2, _wzyx(v2)),
            }; //wzyx
        }
    }

    /// Compute the sum of two vectors and return it as a new vector.
    #[inline(always)]
    pub fn add(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_add_ps(self.data, v2.data),
            }
        }
    }

    /// Subtract a vector and return it as a new vector.
    #[inline(always)]
    pub fn sub(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_sub_ps(self.data, v2.data),
            }
        }
    }

    /// Multiply two vectors component-wise, and return it as a new vector.  
    /// (x1 * x2, y1 * y2, z1 * z2, w1 * w2)
    #[inline(always)]
    pub fn mul(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_mul_ps(self.data, v2.data),
            }
        }
    }

    /// Divide two vectors component-wise, and return it as a new vector.  
    /// (x1 / x2, y1 / y2, z1 / z2, w1 / w2)
    #[inline(always)]
    pub fn div(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_div_ps(self.data, v2.data),
            }
        }
    }

    /// Fused Multiply Add.  Result = (a * b) + c.
    #[inline(always)]
    pub fn mul_add(self, v2: Vector2, v3: Vector2) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _mm_fmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Fused Multiply Sub.  Result = (a * b) - c.
    #[inline(always)]
    pub fn mul_sub(self, v2: Vector2, v3: Vector2) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _mm_fmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Add.  Result = -(a * b) + c.
    #[inline(always)]
    pub fn neg_mul_add(self, v2: Vector2, v3: Vector2) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _mm_fnmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Sub.  Result = -(a * b) - c.
    #[inline(always)]
    pub fn neg_mul_sub(self, v2: Vector2, v3: Vector2) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _mm_fnmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Compute the bitwise AND of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn and<T: SIMDVector2>(self, v2: T) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_and_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise OR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn or<T: SIMDVector2>(self, v2: T) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_or_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise ANDNOT of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn andnot<T: SIMDVector2>(self, v2: T) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_andnot_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise XOR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn xor<T: SIMDVector2>(self, v2: T) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_xor_ps(self.data, v2.data()),
            }
        }
    }

    /// Equals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn equal(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmpeq_ps(self.data, v2.data)),
            }
        }
    }
    /// NotEquals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn not_equal(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmpneq_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater_equal(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmpge_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmpgt_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less_equal(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmple_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less(self, v2: Vector2) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_castps_si128(_mm_cmplt_ps(self.data, v2.data)),
            }
        }
    }

    /// The absolute value of each component of the vector.
    #[inline(always)]
    pub fn abs(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_abs_ps(self.data),
            }
        }
    }

    /// Take the magnitude of the first argument (self), and use the sign of the second argument to produce a new vector
    #[inline(always)]
    pub fn copysign(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_copysign_ps(self.data, v2.data),
            }
        }
    }

    /// Floor function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn floor(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_floor_ps(self.data),
            }
        }
    }

    /// Ceil function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn ceil(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_ceil_ps(self.data),
            }
        }
    }

    /// Round to nearest even function. Returns signed 0 when applicable.
    #[inline(always)]
    pub fn round(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Truncate function.
    #[inline(always)]
    pub fn truncate(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Convert to an integer vector using the floor function.
    #[inline(always)]
    pub fn floor_to_int(self) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_cvttps_epi32(self.floor().data),
            }
        }
    }

    /// Convert to an integer vector using the ceil function.
    #[inline(always)]
    pub fn ceil_to_int(self) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_cvttps_epi32(self.ceil().data),
            }
        }
    }
    /// Convert to an integer vector using the round function.
    #[inline(always)]
    pub fn round_to_int(self) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_cvttps_epi32(self.round().data),
            }
        }
    }

    /// Convert to an integer vector using the truncate function.
    #[inline(always)]
    pub fn truncate_to_int(self) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_cvttps_epi32(self.truncate().data),
            }
        }
    }

    /// Compute the fractional component of each component
    /// Result = X - Floor(x)
    #[inline(always)]
    pub fn frac(self) -> Vector2 {
        return Vector2::sub(self, Vector2::floor(self));
    }

    /// Compute the square root of each component
    #[inline(always)]
    pub fn sqrt(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_sqrt_ps(self.data),
            }
        }
    }

    /// Compute the approximate sin of each component
    #[inline(always)]
    pub fn sin(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_sin_ps(self.data),
            }
        }
    }

    /// Compute the approximate cos of each component
    #[inline(always)]
    pub fn cos(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_cos_ps(self.data),
            }
        }
    }

    /// Compute the approximate acos of each component
    #[inline(always)]
    pub fn acos(self) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_acos_ps(self.data),
            }
        }
    }

    /// Compute the component-wise max.
    #[inline(always)]
    pub fn max(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_max_ps(self.data, v2.data),
            }
        }
    }

    /// Compute the component-wise min.
    #[inline(always)]
    pub fn min(self, v2: Vector2) -> Vector2 {
        unsafe {
            Vector2 {
                data: _mm_min_ps(self.data, v2.data),
            }
        }
    }
    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: Vector2, mask: Vector2Bool) -> Vector2 {
        unsafe {
            Vector2 {
                data: _ico_select_ps(self.data, v2.data, _mm_castsi128_ps(mask.data)),
            }
        }
    }
    /// Linear interpolate from a to b based on a float.
    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(self, v2: Vector2, t: T) -> Vector2 {
        unsafe {
            let t_val = t.into().data;
            let tmp = _mm_fnmadd_ps(self.data, t_val, self.data); //a - (a*t)
            Vector2 {
                data: _mm_fmadd_ps(v2.data, t_val, tmp),
            } //b * t + a
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    #[inline(always)]
    pub fn normalize(self) -> Vector2 {
        let length = FloatVector::sqrt(Vector2::dot(self, self));
        let norm = Vector2::div(self, Vector2::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector2::from(Vector2::dot(norm, norm));
            let mask_less = Vector2::less(
                result_length_sqr,
                Vector2 {
                    data: _ico_two_ps(),
                },
            );
            return Vector2::and(norm, mask_less);
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    /// Also returns the length of the vector prior to normalization.
    #[inline(always)]
    pub fn normalize_length(self) -> (Vector2, FloatVector) {
        let length = FloatVector::sqrt(Vector2::dot(self, self));
        let norm = Vector2::div(self, Vector2::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector2::from(Vector2::dot(norm, norm));
            let mask_less = Vector2::less(
                result_length_sqr,
                Vector2 {
                    data: _ico_two_ps(),
                },
            );
            return (Vector2::and(norm, mask_less), length);
        }
    }

    /// The square magnitude of the vector.  Equal to Dot(self, self).
    #[inline(always)]
    pub fn sqr_magnitude(self) -> FloatVector {
        return Vector2::dot(self, self);
    }

    /// The magnitude of the vector.  Equal to Sqrt(Dot(self, self)).
    #[inline(always)]
    pub fn magnitude(self) -> FloatVector {
        return FloatVector::sqrt(Vector2::dot(self, self));
    }

    #[inline(always)]
    pub fn xxxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyyy(self.data),
            };
        }
    }

    #[inline(always)]
    pub fn xx(self) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _xxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xy(self) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _xyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yx(self) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _yxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yy(self) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _yyzw(self.data),
            };
        }
    }
}

impl From<f32> for Vector2 {
    #[inline(always)]
    fn from(v: f32) -> Vector2 {
        return Vector2::set(v);
    }
}
impl From<FloatVector> for Vector2 {
    #[inline(always)]
    fn from(v: FloatVector) -> Vector2 {
        return Vector2 { data: v.data };
    }
}
impl From<Vector3> for Vector2 {
    #[inline(always)]
    fn from(v: Vector3) -> Vector2 {
        Vector2 { data: v.data }
    }
}
impl From<Vector4> for Vector2 {
    #[inline(always)]
    fn from(v: Vector4) -> Vector2 {
        Vector2 { data: v.data }
    }
}
impl From<Vector2Int> for Vector2 {
    #[inline(always)]
    fn from(v: Vector2Int) -> Vector2 {
        unsafe {
            return Vector2 {
                data: _mm_cvtepi32_ps(v.data),
            };
        }
    }
}

impl core::ops::Add for Vector2 {
    type Output = Vector2;
    #[inline(always)]
    fn add(self, _rhs: Vector2) -> Vector2 {
        Vector2::add(self, _rhs)
    }
}
impl core::ops::AddAssign for Vector2 {
    #[inline(always)]
    fn add_assign(&mut self, other: Vector2) {
        *self = Vector2::add(*self, other)
    }
}
impl core::ops::Sub for Vector2 {
    type Output = Vector2;
    #[inline(always)]
    fn sub(self, _rhs: Vector2) -> Vector2 {
        Vector2::sub(self, _rhs)
    }
}
impl core::ops::SubAssign for Vector2 {
    #[inline(always)]
    fn sub_assign(&mut self, other: Vector2) {
        *self = Vector2::sub(*self, other)
    }
}
impl core::ops::Neg for Vector2 {
    type Output = Vector2;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return Vector2 {
                data: _mm_xor_ps(_ico_signbit_ps(), self.data),
            };
        }
    }
}
impl<T: Into<FloatVector>> core::ops::Mul<T> for Vector2 {
    type Output = Vector2;
    #[inline(always)]
    fn mul(self, _rhs: T) -> Vector2 {
        return Vector2::mul(self, Vector2::from(_rhs.into()));
    }
}
impl<T: Into<FloatVector>> core::ops::MulAssign<T> for Vector2 {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = Vector2::mul(*self, Vector2::from(_rhs.into()));
    }
}
impl core::ops::Mul<Vector2> for FloatVector {
    type Output = Vector2;
    #[inline(always)]
    fn mul(self: FloatVector, _rhs: Vector2) -> Vector2 {
        return Vector2::mul(_rhs, Vector2::from(self));
    }
}

impl<T: Into<FloatVector>> core::ops::Div<T> for Vector2 {
    type Output = Vector2;
    #[inline(always)]
    fn div(self, _rhs: T) -> Vector2 {
        return Vector2::div(self, Vector2::from(_rhs.into()));
    }
}
impl core::ops::Div<Vector2> for FloatVector {
    type Output = Vector2;
    #[inline(always)]
    fn div(self: FloatVector, _rhs: Vector2) -> Vector2 {
        return Vector2::div(Vector2::from(self), _rhs);
    }
}
impl<T: Into<FloatVector>> core::ops::DivAssign<T> for Vector2 {
    #[inline(always)]
    fn div_assign(&mut self, _rhs: T) {
        *self = Vector2::div(*self, Vector2::from(_rhs.into()));
    }
}

impl PartialEq for Vector2 {
    #[inline(always)]
    fn eq(&self, other: &Vector2) -> bool {
        return Vector2::equal(*self, *other).all();
    }
}

impl SIMDVector2 for Vector2 {
    #[inline(always)]
    fn data(self) -> __m128 {
        return self.data;
    }
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}
#[cfg(test)]
mod test;
