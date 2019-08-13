// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::float_vector::FloatVector;
use crate::raw::RawVector_f32;
use crate::sse_extensions::*;
use crate::structure::SIMDVector3;
use crate::vector2::Vector2;
use crate::vector3_bool::Vector3Bool;
use crate::vector3_int::Vector3Int;
use crate::vector4::Vector4;
use core::arch::x86_64::*;

/// A vector of 3 floats (x,y,z).
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector3 {
    pub data: __m128,
}

impl Vector3 {
    /// Construct a new vector from f32 components.
    ///
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// ```
    #[inline(always)]
    pub fn new(x: f32, y: f32, z: f32) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_set_ps(0.0f32, z, y, x),
            }
        }
    }

    /// Set all values of the vector to the same f32 value.
    ///
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let ones = Vector3::set(1.0);
    /// ```
    #[inline(always)]
    pub fn set<T: Into<FloatVector>>(value: T) -> Vector3 {
        return Vector3 {
            data: value.into().data,
        };
    }

    /// Construct a new vector of zeros.
    ///
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let zeros = Vector3::zero();
    /// ```
    #[inline(always)]
    pub fn zero() -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_setzero_ps(),
            }
        }
    }

    /// Get the x value of the vector, broadcast to all components as a FloatVector (xxxx).
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// let one = one_two_three.x();
    /// ```
    #[inline(always)]
    pub fn x(self) -> FloatVector {
        return FloatVector {
            data: self.xxxx().data,
        };
    }

    /// Get the y value of the vector, broadcast to all components as a FloatVector (yyyy).
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// let two = one_two_three.y();
    /// ```
    #[inline(always)]
    pub fn y(self) -> FloatVector {
        return FloatVector {
            data: self.yyyy().data,
        };
    }

    /// Get the z value of the vector, broadcast to all components as a FloatVector (zzzz).
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// let three = one_two_three.z();
    /// ```
    #[inline(always)]
    pub fn z(self) -> FloatVector {
        return FloatVector {
            data: self.zzzz().data,
        };
    }

    /// Load a value from aligned memory.
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// use ico_math::raw::RawVector_f32;
    /// let raw = RawVector_f32{data:[0.0,1.0,2.0,0.0]};
    /// let one_two_three = Vector3::load(&raw);
    /// ```
    #[inline(always)]
    pub fn load(raw: &RawVector_f32) -> Vector3 {
        unsafe {
            // Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
            return Vector3 {
                data: _mm_load_ps(core::mem::transmute(raw)),
            };
        }
    }

    /// Store a value to aligned memory.
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// use ico_math::raw::RawVector_f32;
    /// let mut raw = RawVector_f32{data:[0.0; 4]};
    /// let one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// one_two_three.store(&mut raw);
    /// ```
    #[inline(always)]
    pub fn store(self, dst: &mut RawVector_f32) {
        unsafe {
            // Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
            _mm_store_ps(core::mem::transmute(dst), self.data);
        }
    }

    /// Set the x value of this vector, leaving the other components unchanged.
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let mut one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// one_two_three.set_x(0.5);
    /// ```
    #[inline(always)]
    pub fn set_x<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            self.data = _mm_move_ss(self.data, value.into().data);
        }
    }

    /// Set the y value of this vector, leaving the other components unchanged.
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let mut one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// one_two_three.set_y(0.5);
    /// ```
    #[inline(always)]
    pub fn set_y<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(value.into().data, self.data);
            self.data = _mm_shuffle_ps(v1, self.data, _ico_shuffle(3, 2, 1, 0));
        }
    }

    /// Set the z value of this vector, leaving the other components unchanged.
    /// This doesn't modify W.
    /// # Examples
    ///
    /// ```
    /// use ico_math::Vector3;
    /// let mut one_two_three = Vector3::new(1.0, 2.0, 3.0);
    /// one_two_three.set_z(0.5);
    /// ```
    #[inline(always)]
    pub fn set_z<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(self.data, value.into().data);
            self.data = _mm_shuffle_ps(self.data, v1, _ico_shuffle(3, 0, 1, 0));
        }
    }

    /// Compute the 3 element dot-product, and return it as a broadcast FloatVector.
    #[inline(always)]
    pub fn dot(self, v2: Vector3) -> FloatVector {
        unsafe {
            let tmp0 = _mm_mul_ps(self.data, v2.data); //xyzw
            let tmp1 = _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(tmp0), 4)); //0xyz
            let tmp2 = _mm_add_ps(tmp0, tmp1); //x xy, yz, wz
            let tmp3 = _mm_moveldup_ps(tmp2); // x x yz yz
            FloatVector {
                data: _mm_add_ps(tmp3, _wzyx(tmp3)),
            }
        }
    }

    /// Right handed cross product.
    #[inline(always)]
    pub fn cross(self, rhs: Vector3) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _ico_cross_ps(self.data, rhs.data),
            };
        }
    }

    /// Compute the sum of two vectors and return it as a new vector.
    #[inline(always)]
    pub fn add(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_add_ps(self.data, v2.data),
            }
        }
    }

    /// Subtract a vector and return it as a new vector.
    #[inline(always)]
    pub fn sub(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_sub_ps(self.data, v2.data),
            }
        }
    }

    /// Multiply two vectors component-wise, and return it as a new vector.  
    /// (x1 * x2, y1 * y2, z1 * z2, w1 * w2)
    #[inline(always)]
    pub fn mul(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_mul_ps(self.data, v2.data),
            }
        }
    }

    /// Divide two vectors component-wise, and return it as a new vector.  
    /// (x1 / x2, y1 / y2, z1 / z2, w1 / w2)
    #[inline(always)]
    pub fn div(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_div_ps(self.data, v2.data),
            }
        }
    }

    /// Fused Multiply Add.  Result = (a * b) + c.
    #[inline(always)]
    pub fn mul_add(self, v2: Vector3, v3: Vector3) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_fmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Fused Multiply Sub.  Result = (a * b) - c.
    #[inline(always)]
    pub fn mul_sub(self, v2: Vector3, v3: Vector3) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_fmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Add.  Result = -(a * b) + c.
    #[inline(always)]
    pub fn neg_mul_add(self, v2: Vector3, v3: Vector3) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_fnmadd_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Negated Fused Multiply Sub.  Result = -(a * b) - c.
    #[inline(always)]
    pub fn neg_mul_sub(self, v2: Vector3, v3: Vector3) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_fnmsub_ps(self.data, v2.data, v3.data),
            };
        }
    }

    /// Compute the bitwise AND of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn and<T: SIMDVector3>(self, v2: T) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_and_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise OR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn or<T: SIMDVector3>(self, v2: T) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_or_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise ANDNOT of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn andnot<T: SIMDVector3>(self, v2: T) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_andnot_ps(self.data, v2.data()),
            }
        }
    }

    /// Compute the bitwise XOR of two vectors.
    /// This function treats inputs as binary data, and doesn't perform any conversions.
    #[inline(always)]
    pub fn xor<T: SIMDVector3>(self, v2: T) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_xor_ps(self.data, v2.data()),
            }
        }
    }

    /// Equals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn equal(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_castps_si128(_mm_cmpeq_ps(self.data, v2.data)),
            }
        }
    }
    /// NotEquals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn not_equal(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_castps_si128(_mm_cmpneq_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater_equal(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_castps_si128(_mm_cmpge_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_castps_si128(_mm_cmpgt_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less_equal(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_castps_si128(_mm_cmple_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less(self, v2: Vector3) -> Vector3Bool {
        unsafe {
            Vector3Bool {
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
    pub fn approx_equal(self, v2: Vector3) -> Vector3Bool {
        let delta = (self - v2).abs();
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector3::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.less_equal(epsilon_bound);
    }

    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_greater(self, v2: Vector3) -> Vector3Bool {
        let delta = self.sub(v2);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector3::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }
    /// Adapted from Knuth with an added absolute epsilon
    /// return (a - b) > max(machine_epsilon, (max( abs(a), abs(b) ) * relative_epsilon);
    #[inline(always)]
    pub fn definitely_less(self, v2: Vector3) -> Vector3Bool {
        let delta = v2.sub(self);
        let abs_a = self.abs();
        let abs_b = v2.abs();
        let epsilon_bound = (abs_a.max(abs_b) * RELATIVE_COMPARISON_EPSILON)
            .max(Vector3::from(ABSOLUTE_COMPARISON_EPSILON));
        return delta.greater(epsilon_bound);
    }

    /// The absolute value of each component of the vector.
    #[inline(always)]
    pub fn abs(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_abs_ps(self.data),
            }
        }
    }

    /// Take the magnitude of the first argument (self), and use the sign of the second argument to produce a new vector
    #[inline(always)]
    pub fn copysign(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_copysign_ps(self.data, v2.data),
            }
        }
    }

    /// Floor function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn floor(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_floor_ps(self.data),
            }
        }
    }

    /// Ceil function.  Returns signed 0 when applicable.
    #[inline(always)]
    pub fn ceil(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_ceil_ps(self.data),
            }
        }
    }

    /// Round to nearest even function. Returns signed 0 when applicable.
    #[inline(always)]
    pub fn round(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Truncate function.
    #[inline(always)]
    pub fn truncate(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_round_ps(self.data, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC),
            }
        }
    }

    /// Convert to an integer vector using the floor function.
    #[inline(always)]
    pub fn floor_to_int(self) -> Vector3Int {
        unsafe {
            Vector3Int {
                data: _mm_cvttps_epi32(self.floor().data),
            }
        }
    }

    /// Convert to an integer vector using the ceil function.
    #[inline(always)]
    pub fn ceil_to_int(self) -> Vector3Int {
        unsafe {
            Vector3Int {
                data: _mm_cvttps_epi32(self.ceil().data),
            }
        }
    }
    /// Convert to an integer vector using the round function.
    #[inline(always)]
    pub fn round_to_int(self) -> Vector3Int {
        unsafe {
            Vector3Int {
                data: _mm_cvttps_epi32(self.round().data),
            }
        }
    }

    /// Convert to an integer vector using the truncate function.
    #[inline(always)]
    pub fn truncate_to_int(self) -> Vector3Int {
        unsafe {
            Vector3Int {
                data: _mm_cvttps_epi32(self.truncate().data),
            }
        }
    }

    /// Compute the fractional component of each component
    /// Result = X - Floor(x)
    #[inline(always)]
    pub fn frac(self) -> Vector3 {
        return Vector3::sub(self, Vector3::floor(self));
    }

    /// Compute the square root of each component
    #[inline(always)]
    pub fn sqrt(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_sqrt_ps(self.data),
            }
        }
    }

    /// Compute the approximate sin of each component
    #[inline(always)]
    pub fn sin(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_sin_ps(self.data),
            }
        }
    }

    /// Compute the approximate cos of each component
    #[inline(always)]
    pub fn cos(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_cos_ps(self.data),
            }
        }
    }
    
    /// Compute the approximate tan of each component
    #[inline(always)]
    pub fn tan(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_tan_ps(self.data),
            }
        }
    }

    /// Compute the approximate acos of each component
    #[inline(always)]
    pub fn acos(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_acos_ps(self.data),
            }
        }
    }

    /// Compute the approximate asin of each component
    #[inline(always)]
    pub fn asin(self) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_asin_ps(self.data),
            }
        }
    }

    /// Compute the component-wise max.
    #[inline(always)]
    pub fn max(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_max_ps(self.data, v2.data),
            }
        }
    }

    /// Compute the component-wise min.
    #[inline(always)]
    pub fn min(self, v2: Vector3) -> Vector3 {
        unsafe {
            Vector3 {
                data: _mm_min_ps(self.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn horizontal_min(self) -> FloatVector {
        let x = self.x();
        let xy = x.min(self.y());
        return xy.min(self.z());
    }

    #[inline(always)]
    pub fn horizontal_max(self) -> FloatVector {
        let x = self.x();
        let xy = x.max(self.y());
        return xy.max(self.z());
    }

    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: Vector3, mask: Vector3Bool) -> Vector3 {
        unsafe {
            Vector3 {
                data: _ico_select_ps(self.data, v2.data, _mm_castsi128_ps(mask.data)),
            }
        }
    }
    /// Linear interpolate from a to b based on a float.
    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(self, v2: Vector3, t: T) -> Vector3 {
        unsafe {
            let t_val = t.into().data;
            let tmp = _mm_fnmadd_ps(self.data, t_val, self.data); //a - (a*t)
            Vector3 {
                data: _mm_fmadd_ps(v2.data, t_val, tmp),
            } //b * t + a
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    #[inline(always)]
    pub fn normalize(self) -> Vector3 {
        let length = FloatVector::sqrt(Vector3::dot(self, self));
        let norm = Vector3::div(self, Vector3::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector3::from(Vector3::dot(norm, norm));
            let mask_less = Vector3::less(
                result_length_sqr,
                Vector3 {
                    data: _ico_two_ps(),
                },
            );
            return Vector3::and(norm, mask_less);
        }
    }

    /// Normalize the vector.  Returns a zero vector if the result would be infinity or NAN.
    /// Also returns the length of the vector prior to normalization.
    #[inline(always)]
    pub fn normalize_length(self) -> (Vector3, FloatVector) {
        let length = FloatVector::sqrt(Vector3::dot(self, self));
        let norm = Vector3::div(self, Vector3::from(length));

        unsafe {
            // This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
            let result_length_sqr = Vector3::from(Vector3::dot(norm, norm));
            let mask_less = Vector3::less(
                result_length_sqr,
                Vector3 {
                    data: _ico_two_ps(),
                },
            );
            return (Vector3::and(norm, mask_less), length);
        }
    }

    /// The square magnitude of the vector.  Equal to Dot(self, self).
    #[inline(always)]
    pub fn sqr_magnitude(self) -> FloatVector {
        return Vector3::dot(self, self);
    }

    /// The magnitude of the vector.  Equal to Sqrt(Dot(self, self)).
    #[inline(always)]
    pub fn magnitude(self) -> FloatVector {
        return FloatVector::sqrt(Vector3::dot(self, self));
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
    pub fn zzzz(self) -> Vector4 {
        unsafe {
            return Vector4 {
                data: _zzzz(self.data),
            };
        }
    }

    #[inline(always)]
    pub fn xxx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xxz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xyz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xzz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _xzzw(self.data),
            };
        }
    }

    #[inline(always)]
    pub fn yxx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yxz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yzz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _yzzw(self.data),
            };
        }
    }

    #[inline(always)]
    pub fn zxx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zxxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zxyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zxz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zxzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zyxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zyyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zyz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zyzw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzx(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zzxw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzy(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zzyw(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn zzz(self) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _zzzw(self.data),
            };
        }
    }
}

impl From<f32> for Vector3 {
    #[inline(always)]
    fn from(v: f32) -> Vector3 {
        return Vector3::set(v);
    }
}

// impl From<Vector3> for [f32;3] {
//     #[inline(always)]
//     fn from(v : Vector3) -> [f32;3] {
//         unsafe{
//         let mut raw= core::mem::MaybeUninit::<RawVector_f32>::uninit().assume_init();//RawVector_f32{data:core::mem::MaybeUninit::<[f32;4]>::uninit().assume_init()};// = RawVector_f32{data:[0.0;4]};
//         v.store(& mut raw);
//         return raw.data_3;
//         }
//     }
// }

impl From<FloatVector> for Vector3 {
    #[inline(always)]
    fn from(v: FloatVector) -> Vector3 {
        return Vector3 { data: v.data };
    }
}

impl From<Vector2> for Vector3 {
    #[inline(always)]
    fn from(v: Vector2) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_movelh_ps(v.data, _mm_setzero_ps()),
            };
        }
    }
}

impl From<Vector4> for Vector3 {
    #[inline(always)]
    fn from(v: Vector4) -> Vector3 {
        Vector3 { data: v.data }
    }
}

impl From<Vector3Int> for Vector3 {
    #[inline(always)]
    fn from(v: Vector3Int) -> Vector3 {
        unsafe {
            return Vector3 {
                data: _mm_cvtepi32_ps(v.data),
            };
        }
    }
}

impl core::ops::Add for Vector3 {
    type Output = Vector3;
    #[inline(always)]
    fn add(self, _rhs: Vector3) -> Vector3 {
        Vector3::add(self, _rhs)
    }
}

impl core::ops::AddAssign for Vector3 {
    #[inline(always)]
    fn add_assign(&mut self, other: Vector3) {
        *self = Vector3::add(*self, other)
    }
}

impl core::ops::Sub for Vector3 {
    type Output = Vector3;
    #[inline(always)]
    fn sub(self, _rhs: Vector3) -> Vector3 {
        Vector3::sub(self, _rhs)
    }
}

impl core::ops::SubAssign for Vector3 {
    #[inline(always)]
    fn sub_assign(&mut self, other: Vector3) {
        *self = Vector3::sub(*self, other)
    }
}

impl core::ops::Neg for Vector3 {
    type Output = Vector3;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return Vector3 {
                data: _mm_xor_ps(_ico_signbit_ps(), self.data),
            };
        }
    }
}

impl<T: Into<FloatVector>> core::ops::Mul<T> for Vector3 {
    type Output = Vector3;
    #[inline(always)]
    fn mul(self, _rhs: T) -> Vector3 {
        return Vector3::mul(self, Vector3::from(_rhs.into()));
    }
}

impl<T: Into<FloatVector>> core::ops::MulAssign<T> for Vector3 {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = Vector3::mul(*self, Vector3::from(_rhs.into()));
    }
}

impl core::ops::Mul<Vector3> for FloatVector {
    type Output = Vector3;
    #[inline(always)]
    fn mul(self: FloatVector, _rhs: Vector3) -> Vector3 {
        return Vector3::mul(_rhs, Vector3::from(self));
    }
}

impl<T: Into<FloatVector>> core::ops::Div<T> for Vector3 {
    type Output = Vector3;
    #[inline(always)]
    fn div(self, _rhs: T) -> Vector3 {
        return Vector3::div(self, Vector3::from(_rhs.into()));
    }
}

impl core::ops::Div<Vector3> for FloatVector {
    type Output = Vector3;
    #[inline(always)]
    fn div(self: FloatVector, _rhs: Vector3) -> Vector3 {
        return Vector3::div(Vector3::from(self), _rhs);
    }
}

impl<T: Into<FloatVector>> core::ops::DivAssign<T> for Vector3 {
    #[inline(always)]
    fn div_assign(&mut self, _rhs: T) {
        *self = Vector3::div(*self, Vector3::from(_rhs.into()));
    }
}

impl PartialEq for Vector3 {
    #[inline(always)]
    fn eq(&self, other: &Vector3) -> bool {
        return Vector3::equal(*self, *other).all();
    }
}

impl SIMDVector3 for Vector3 {
    #[inline(always)]
    fn data(self) -> __m128 {
        return self.data;
    }
    #[inline(always)]
    fn data_i(self) -> __m128i {
        return unsafe { _mm_castps_si128(self.data) };
    }
}

#[cfg(test)]
mod test;
