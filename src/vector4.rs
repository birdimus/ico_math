use crate::float_vector::FloatVector;
use crate::vector2::Vector2;
use crate::vector3::Vector3;
use core::arch::x86_64::*;

use crate::quaternion::Quaternion;
use crate::sse_extensions::*;
use crate::structure::SIMDVector4;
use crate::vector4_bool::Vector4Bool;
use crate::vector4_int::Vector4Int;
use crate::raw::RawVector_f32;

/// A vector of 4 floats (x,y,z,w).
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4 {
    pub data: __m128,
}

impl Vector4 {
    /// Construct a new vector from f32 components.
    #[inline(always)]
    pub fn new(x: f32, y: f32, z: f32, w: f32) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_set_ps(w, z, y, x),
            }
        }
    }

    /// Set all values of the vector to the same f32 value.
    #[inline(always)]
    pub fn set<T: Into<FloatVector>>(value: T) -> Vector4 {
        return Vector4 {
            data: value.into().data,
        };
    }

    /// Construct a new vector of zeros.
    #[inline(always)]
    pub fn zero() -> Vector4 {
        unsafe {
            Vector4 {
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

    /// Get the z value of the vector, broadcast to all components as a FloatVector (zzzz).
    #[inline(always)]
    pub fn z(self) -> FloatVector {
        return FloatVector {
            data: self.zzzz().data,
        };
    }

    /// Get the w value of the vector, broadcast to all components as a FloatVector (wwww).
    #[inline(always)]
    pub fn w(self) -> FloatVector {
        return FloatVector {
            data: self.wwww().data,
        };
    }

    /// Load a value from aligned memory.
    #[inline(always)]
    pub fn load(raw : &RawVector_f32) -> Vector4{
    	unsafe{
    		// Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
    		return Vector4{data:_mm_load_ps(core::mem::transmute(raw))};
    		
    	}
    }

    /// Store a value to aligned memory.
    #[inline(always)]
    pub fn store(self,dst : &mut RawVector_f32){
    	unsafe{
    		// Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
    		_mm_store_ps(core::mem::transmute(dst), self.data);
    	}
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

    /// Set the z value of this vector, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_z<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(self.data, value.into().data);
            self.data = _mm_shuffle_ps(self.data, v1, _ico_shuffle(3, 0, 1, 0));
        }
    }

    /// Set the w value of this vector, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_w<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(self.data, value.into().data);
            self.data = _mm_shuffle_ps(self.data, v1, _ico_shuffle(0, 2, 1, 0));
        }
    }



    /// Compute the 4 element dot-product, and return it as a broadcast FloatVector.
    #[inline(always)]
    pub fn dot(self, v2: Vector4) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _ico_dp4_ps(self.data, v2.data),
            };
        }
    }

    /// Compute the sum of two vectors and return it as a new vector.
    #[inline(always)]
    pub fn add(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_add_ps(self.data, v2.data),
            }
        }
    }

    /// Subtract a vector and return it as a new vector.
    #[inline(always)]
    pub fn sub(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_sub_ps(self.data, v2.data),
            }
        }
    }

    /// Multiply two vectors component-wise, and return it as a new vector.  
    /// (x1 * x2, y1 * y2, z1 * z2, w1 * w2)
    #[inline(always)]
    pub fn mul(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_mul_ps(self.data, v2.data),
            }
        }
    }

    /// Divide two vectors component-wise, and return it as a new vector.  
    /// (x1 / x2, y1 / y2, z1 / z2, w1 / w2)
    #[inline(always)]
    pub fn div(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_div_ps(self.data, v2.data),
            }
        }
    }

    /// Fused Multiply Add.  Result = (a * b) + c.
    #[inline(always)]
    pub fn mul_add(self, v2: Vector4, v3: Vector4) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_fmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Fused Multiply Sub.  Result = (a * b) - c.
    #[inline(always)]
    pub fn mul_sub(self, v2: Vector4, v3: Vector4) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_fmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Add.  Result = -(a * b) + c.
    #[inline(always)]
    pub fn neg_mul_add(self, v2: Vector4, v3: Vector4) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_fnmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Sub.  Result = -(a * b) - c.
    #[inline(always)]
    pub fn neg_mul_sub(self, v2: Vector4, v3: Vector4) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_fnmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Compute the bitwise AND of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn and<T: SIMDVector4>(self, v2: T) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_and_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise OR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn or<T: SIMDVector4>(self, v2: T) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_or_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise ANDNOT of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn andnot<T: SIMDVector4>(self, v2: T) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_andnot_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise XOR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn xor<T: SIMDVector4>(self, v2: T) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_xor_ps(self.data, v2.data()),
            }
        }
    }

    /// Equals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn equal(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpeq_ps(self.data, v2.data)),
            }
        }
    }

    /// NotEquals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn not_equal(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpneq_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater_equal(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpge_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpgt_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less_equal(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmple_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less(self, v2: Vector4) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmplt_ps(self.data, v2.data)),
            }
        }
    }

    /// Relative and absolute epsilon comparison.  
    /// Uses machine epsilon as absolute, and 4*machine epsilon for relative.
    /// return abs(a - b) <= max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    /// Adapted from Knuth.  
    /// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    #[inline(always)]
    pub fn approx_equal(self, v2: Vector4) -> Vector4Bool {
        let delta = (self - v2).abs();
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector4::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.less_equal(epsilon_bound);
    }

    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_greater(self, v2: Vector4) -> Vector4Bool {
        let delta = self.sub(v2);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector4::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }
    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_less(self, v2: Vector4) -> Vector4Bool {
        let delta = v2.sub(self);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector4::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }

    /// The absolute value of each component of the vector.
    #[inline(always)]
    pub fn abs(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_abs_ps(self.data),
            }
        }
    }

    /// Take the magnitude of the first argument (self), and use the sign of the second argument to produce a new vector
    #[inline(always)]
    pub fn copysign(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_copysign_ps(self.data, v2.data),
            }
        }
    }

    /// Floor function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn floor(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_floor_ps(self.data),
            }
        }
    }

    /// Ceil function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn ceil(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_ceil_ps(self.data),
            }
        }
    }

    /// Round to nearest even function. Returns signed 0 when applicable.
    #[inline(always)]
    pub fn round(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Truncate function.
    #[inline(always)]
    pub fn truncate(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Convert to an integer vector using the floor function.
    #[inline(always)]
    pub fn floor_to_int(self) -> Vector4Int {
        unsafe {
            Vector4Int {
                data: _mm_cvttps_epi32(self.floor().data),
            }
        }
    }

    /// Convert to an integer vector using the ceil function.
    #[inline(always)]
    pub fn ceil_to_int(self) -> Vector4Int {
        unsafe {
            Vector4Int {
                data: _mm_cvttps_epi32(self.ceil().data),
            }
        }
    }
    /// Convert to an integer vector using the round function.
    #[inline(always)]
    pub fn round_to_int(self) -> Vector4Int {
        unsafe {
            Vector4Int {
                data: _mm_cvttps_epi32(self.round().data),
            }
        }
    }

    /// Convert to an integer vector using the truncate function.
    #[inline(always)]
    pub fn truncate_to_int(self) -> Vector4Int {
        unsafe {
            Vector4Int {
                data: _mm_cvttps_epi32(self.truncate().data),
            }
        }
    }

    /// Compute the fractional component of each component
    /// Result = X - Floor(x)
    #[inline(always)]
    pub fn frac(self) -> Vector4 {
        return Vector4::sub(self, Vector4::floor(self));
    }

    /// Compute the square root of each component
    #[inline(always)]
    pub fn sqrt(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_sqrt_ps(self.data),
            }
        }
    }

    /// Compute the approximate sin of each component
    #[inline(always)]
    pub fn sin(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_sin_ps(self.data),
            }
        }
    }

    /// Compute the approximate cos of each component
    #[inline(always)]
    pub fn cos(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_cos_ps(self.data),
            }
        }
    }

    /// Compute the approximate acos of each component
    #[inline(always)]
    pub fn acos(self) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_acos_ps(self.data),
            }
        }
    }

    /// Compute the component-wise max.
    #[inline(always)]
    pub fn max(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_max_ps(self.data, v2.data),
            }
        }
    }

    /// Compute the component-wise min.
    #[inline(always)]
    pub fn min(self, v2: Vector4) -> Vector4 {
        unsafe {
            Vector4 {
                data: _mm_min_ps(self.data, v2.data),
            }
        }
    }
    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: Vector4, mask: Vector4Bool) -> Vector4 {
        unsafe {
            Vector4 {
                data: _ico_select_ps(self.data, v2.data, _mm_castsi128_ps(mask.data)),
            }
        }
    }
    /// Linear interpolate from a to b based on a float.
    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(self, v2: Vector4, t: T) -> Vector4 {
        unsafe {
            let t_val = t.into().data;
            let tmp = _mm_fnmadd_ps(self.data, t_val, self.data); //a - (a*t)
            Vector4 {
                data: _mm_fmadd_ps(v2.data, t_val, tmp),
            } //b * t + a
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    #[inline(always)]
    pub fn normalize(self) -> Vector4 {
        let length = FloatVector::sqrt(Vector4::dot(self, self));
        let norm = Vector4::div(self, Vector4::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector4::from(Vector4::dot(norm, norm));
            let mask_less = Vector4::less(
                result_length_sqr,
                Vector4 {
                    data: _ico_two_ps(),
                },
            );
            return Vector4::and(norm, mask_less);
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    /// Also returns the length of the vector prior to normalization.
    #[inline(always)]
    pub fn normalize_length(self) -> (Vector4, FloatVector) {
        let length = FloatVector::sqrt(Vector4::dot(self, self));
        let norm = Vector4::div(self, Vector4::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector4::from(Vector4::dot(norm, norm));
            let mask_less = Vector4::less(
                result_length_sqr,
                Vector4 {
                    data: _ico_two_ps(),
                },
            );
            return (Vector4::and(norm, mask_less), length);
        }
    }

    /// The square magnitude of the vector.  Equal to Dot(self, self).
    #[inline(always)]
    pub fn sqr_magnitude(self) -> FloatVector {
        return Vector4::dot(self, self);
    }

    /// The magnitude of the vector.  Equal to Sqrt(Dot(self, self)).
    #[inline(always)]
    pub fn magnitude(self) -> FloatVector {
        return FloatVector::sqrt(Vector4::dot(self, self));
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
    pub fn xxxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xxww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xywx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xywx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xywy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xywy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xywz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xywz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xyww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xzww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xwww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _xwww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yxww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyyx(self.data),
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
    pub fn yyyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yywx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yywx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yywy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yywy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yywz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yywz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yyww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _yzww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn ywww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _ywww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zxww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zywx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zywx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zywy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zywy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zywz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zywz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zyww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zwww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zwww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wxww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wxww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wywx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wywx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wywy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wywy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wywz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wywz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wyww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wyww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wzww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wzww(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwxx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwxx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwxy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwxy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwxz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwxz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwxw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwyx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwyx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwyy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwyy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwyz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwyz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwyw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwzx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwzx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwzy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwzy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwzz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwzw(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwwx(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwwx(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwwy(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwwy(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwwz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwwz(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn wwww(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _wwww(self.data),
            };
        }
    }
}

impl From<f32> for Vector4 {
    #[inline(always)]
    fn from(v: f32) -> Vector4 {
        return Vector4::set(v);
    }
}
impl From<FloatVector> for Vector4 {
    #[inline(always)]
    fn from(v: FloatVector) -> Vector4 {
        return Vector4 { data: v.data };
    }
}
impl From<Vector2> for Vector4 {
    #[inline(always)]
    fn from(v: Vector2) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_movelh_ps(v.data, _mm_setzero_ps()),
            };
        }
    }
}
impl From<Vector3> for Vector4 {
    #[inline(always)]
    fn from(v: Vector3) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_castsi128_ps(_mm_srli_si128(
                    _mm_slli_si128(_mm_castps_si128(v.data), 4),
                    4,
                )),
            };
        }
    }
}
impl From<Vector4Int> for Vector4 {
    #[inline(always)]
    fn from(v: Vector4Int) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _mm_cvtepi32_ps(v.data),
            };
        }
    }
}
impl From<Quaternion> for Vector4 {
    #[inline(always)]
    fn from(v: Quaternion) -> Vector4 {
        Vector4 { data: v.data }
    }
}
impl core::ops::Add for Vector4 {
    type Output = Vector4;
    #[inline(always)]
    fn add(self, _rhs: Vector4) -> Vector4 {
        Vector4::add(self, _rhs)
    }
}
impl core::ops::AddAssign for Vector4 {
    #[inline(always)]
    fn add_assign(&mut self, other: Vector4) {
        *self = Vector4::add(*self, other)
    }
}
impl core::ops::Sub for Vector4 {
    type Output = Vector4;
    #[inline(always)]
    fn sub(self, _rhs: Vector4) -> Vector4 {
        Vector4::sub(self, _rhs)
    }
}
impl core::ops::SubAssign for Vector4 {
    #[inline(always)]
    fn sub_assign(&mut self, other: Vector4) {
        *self = Vector4::sub(*self, other)
    }
}
impl core::ops::Neg for Vector4 {
    type Output = Vector4;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return Vector4 {
                data: _mm_xor_ps(_ico_signbit_ps(), self.data),
            };
        }
    }
}
impl<T: Into<FloatVector>> core::ops::Mul<T> for Vector4 {
    type Output = Vector4;
    #[inline(always)]
    fn mul(self, _rhs: T) -> Vector4 {
        return Vector4::mul(self, Vector4::from(_rhs.into()));
    }
}
impl<T: Into<FloatVector>> core::ops::MulAssign<T> for Vector4 {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = Vector4::mul(*self, Vector4::from(_rhs.into()));
    }
}
impl core::ops::Mul<Vector4> for FloatVector {
    type Output = Vector4;
    #[inline(always)]
    fn mul(self: FloatVector, _rhs: Vector4) -> Vector4 {
        return Vector4::mul(_rhs, Vector4::from(self));
    }
}

impl<T: Into<FloatVector>> core::ops::Div<T> for Vector4 {
    type Output = Vector4;
    #[inline(always)]
    fn div(self, _rhs: T) -> Vector4 {
        return Vector4::div(self, Vector4::from(_rhs.into()));
    }
}
impl core::ops::Div<Vector4> for FloatVector {
    type Output = Vector4;
    #[inline(always)]
    fn div(self: FloatVector, _rhs: Vector4) -> Vector4 {
        return Vector4::div(Vector4::from(self), _rhs);
    }
}
impl<T: Into<FloatVector>> core::ops::DivAssign<T> for Vector4 {
    #[inline(always)]
    fn div_assign(&mut self, _rhs: T) {
        *self = Vector4::div(*self, Vector4::from(_rhs.into()));
    }
}

impl PartialEq for Vector4 {
    #[inline(always)]
    fn eq(&self, other: &Vector4) -> bool {
        return Vector4::equal(*self, *other).all();
    }
}

impl SIMDVector4 for Vector4 {
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
