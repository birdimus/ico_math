// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::float_vector::FloatVector;
use crate::raw::RawVector_f32;
use crate::sse_extensions::*;
use crate::vector3::Vector3;
use crate::vector4::Vector4;
use crate::vector4_bool::Vector4Bool;
use core::arch::x86_64::*;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Quaternion {
    pub data: __m128,
}

pub enum RotationOrder {
    XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX,
}

// The value at which slerp switches to lerp.
const SLERP_EPSILON: f32 = 0.9995;
const EQUALITY_EPSILON_LOWER_BOUND: f32 = 1.0 - 0.000001;
const EQUALITY_EPSILON_UPPER_BOUND: f32 = 1.0 + 0.000001;
impl Quaternion {
    /// Construct a new quaternion from f32 components.  This may result in an invalid quaternion.
    #[inline(always)]
    pub fn new(x: f32, y: f32, z: f32, w: f32) -> Quaternion {
        unsafe {
            Quaternion {
                data: _mm_set_ps(w, z, y, x),
            }
        }
    }

    /// Construct a new identity quaternion.
    #[inline(always)]
    pub fn identity() -> Quaternion {
        unsafe {
            Quaternion {
                data: _mm_set_ps(1.0f32, 0.0f32, 0.0f32, 0.0f32),
            }
        }
    }
    /// Load a value from aligned memory.
    #[inline(always)]
    pub fn load(raw: &RawVector_f32) -> Quaternion {
        unsafe {
            // Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
            return Quaternion {
                data: _mm_load_ps(core::mem::transmute(raw)),
            };
        }
    }

    /// Store a value to aligned memory.
    #[inline(always)]
    pub fn store(self, dst: &mut RawVector_f32) {
        unsafe {
            // Use the sledgehammer cast here.  It's fine because RawVector is aligned and c-like.
            _mm_store_ps(core::mem::transmute(dst), self.data);
        }
    }

    /// Construct a new quaternion from axis-angle format.
    #[inline(always)]
    pub fn angle_axis<T: Into<FloatVector>>(radians: T, axis: Vector3) -> Quaternion {
        let half_angle = radians.into() * 0.5f32;
        let sin_angle = half_angle.sin();
        let cos_angle = half_angle.cos();

        let sin_axis = axis * sin_angle;
        let mut q = Quaternion {
            data: sin_axis.data,
        };
        q.set_w(cos_angle);
        return q;
    }

    /// Get the x value of the quaternion, broadcast to all components as a FloatVector (xxxx).
    #[inline(always)]
    pub fn x(self) -> FloatVector {
        return FloatVector {
            data: Vector4::from(self).xxxx().data,
        };
    }

    /// Get the y value of the quaternion, broadcast to all components as a FloatVector (yyyy).
    #[inline(always)]
    pub fn y(self) -> FloatVector {
        return FloatVector {
            data: Vector4::from(self).yyyy().data,
        };
    }

    /// Get the z value of the quaternion, broadcast to all components as a FloatVector (zzzz).
    #[inline(always)]
    pub fn z(self) -> FloatVector {
        return FloatVector {
            data: Vector4::from(self).zzzz().data,
        };
    }

    /// Get the w value of the quaternion, broadcast to all components as a FloatVector (wwww).
    #[inline(always)]
    pub fn w(self) -> FloatVector {
        return FloatVector {
            data: Vector4::from(self).wwww().data,
        };
    }

    /// Set the x value of this quaternion, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_x<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            self.data = _mm_move_ss(self.data, value.into().data);
        }
    }

    /// Set the y value of this quaternion, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_y<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(value.into().data, self.data);
            self.data = _mm_shuffle_ps(v1, self.data, _ico_shuffle(3, 2, 1, 0));
        }
    }

    /// Set the z value of this quaternion, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_z<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(self.data, value.into().data);
            self.data = _mm_shuffle_ps(self.data, v1, _ico_shuffle(3, 0, 1, 0));
        }
    }

    /// Set the w value of this quaternion, leaving the other components unchanged.
    #[inline(always)]
    pub fn set_w<T: Into<FloatVector>>(&mut self, value: T) {
        unsafe {
            let v1 = _mm_move_ss(self.data, value.into().data);
            self.data = _mm_shuffle_ps(self.data, v1, _ico_shuffle(0, 2, 1, 0));
        }
    }

    #[inline(always)]
    pub fn mul_vec(quat: Quaternion, vec: Vector3) -> Vector3 {
        // t = 2.0 * cross(q.xyz, vec)
        let qv = Vector3 { data: quat.data };
        let t = Vector3::cross(qv, vec);
        let t2 = Vector3::add(t, t);

        // v' = v + q.w * t + cross(q.xyz, t)
        let tmp = Vector3::mul_add(Vector3::from(quat.w()), t2, vec);

        return Vector3::add(tmp, Vector3::cross(qv, t2));
    }

    #[inline(always)]
    pub fn mul(lhs: Quaternion, rhs: Quaternion) -> Quaternion {
        unsafe {
            return Quaternion {
                data: _ico_quat_mul(lhs.data, rhs.data),
            };
        }
    }

    /// Equals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn equal(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpeq_ps(self.data, v2.data)),
            }
        }
    }
    /// NotEquals, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn not_equal(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpneq_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater_equal(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpge_ps(self.data, v2.data)),
            }
        }
    }
    /// Greater than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn greater(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmpgt_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than or equal to, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less_equal(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmple_ps(self.data, v2.data)),
            }
        }
    }
    /// Less than, computed component-wise.  This compares bits, and is exact.
    #[inline(always)]
    pub fn less(self, v2: Quaternion) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_castps_si128(_mm_cmplt_ps(self.data, v2.data)),
            }
        }
    }

    #[inline(always)]
    pub fn dot(self, q2: Quaternion) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _ico_dp4_ps(self.data, q2.data),
            };
        }
    }

    #[inline(always)]
    pub fn sqr_magnitude(self) -> FloatVector {
        return Quaternion::dot(self, self);
    }

    #[inline(always)]
    pub fn magnitude(self) -> FloatVector {
        return FloatVector::sqrt(Quaternion::dot(self, self));
    }

    #[inline(always)]
    pub fn inverse(self) -> Quaternion {
        let mask = Vector4::new(SIGN_BIT, SIGN_BIT, SIGN_BIT, 0.0f32);
        return Quaternion::from(Vector4::xor(Vector4::from(self), mask));
    }

    #[inline(always)]
    pub fn reverse(self) -> Quaternion {
        let mask = FloatVector::new(SIGN_BIT);
        return Quaternion::from(Vector4::xor(Vector4::from(self), mask));
    }

    #[inline(always)]
    pub fn delta(from: Quaternion, to: Quaternion) -> Quaternion {
        return Quaternion::mul(from.inverse(), to);
    }

    #[inline(always)]
    pub fn shortest_delta(self, to: Quaternion) -> Quaternion {
        let dp = Quaternion::dot(self, to);
        let negative = Vector4::less(Vector4::from(dp), Vector4::zero());
        let sign_flip = Vector4::and(Vector4::set(SIGN_BIT), negative);
        let shortest_from = Vector4::xor(Vector4::from(self), sign_flip);
        let inv = Quaternion::inverse(Quaternion::from(shortest_from));
        return Quaternion::from(Vector4::mul(Vector4::from(inv), Vector4::from(to)));
    }

    #[inline(always)]
    pub fn angle(self, b: Quaternion) -> FloatVector {
        let dot = Quaternion::dot(self, b);
        let val = dot.acos();
        return val.add(val);
    }

    /// A fast normalize (q = q / q.magnitude).  When the length is 0, this will divide by 0 - and may produce infinity or NaN.
    /// Only should be used if the Quaternion is known to be non-zero.
    #[inline(always)]
    pub fn renormalize(self) -> Quaternion {
        let len = FloatVector::sqrt(Quaternion::dot(self, self));
        return Quaternion::from(Vector4::from(self) / len);
    }

    /// A safe normalize.  If the quaternion is 0, returns Identity.
    #[inline(always)]
    pub fn normalize(self) -> Quaternion {
        let len = FloatVector::sqrt(Quaternion::dot(self, self));
        let scaled = Vector4::from(self) / len;
        let mask = Vector4::less(scaled.abs(), Vector4::set(core::f32::INFINITY));
        return Quaternion::from(Vector4::select(
            Vector4::from(Quaternion::identity()),
            scaled,
            mask,
        ));
    }
    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: Quaternion, mask: Vector4Bool) -> Quaternion {
        unsafe {
            return Quaternion {
                data: _ico_select_ps(self.data, v2.data, _mm_castsi128_ps(mask.data)),
            };
        }
    }

    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(self, to: Quaternion, t: T) -> Quaternion {
        // If the sign of T is negative, flip this around.
        let t_tmp = t.into();
        let t_negative = Vector4::from(t_tmp).less(Vector4::zero());

        let sign_flipped_to = Quaternion::select(to, to.inverse(), t_negative);
        let t_vec = t_tmp.abs();

        let cos_theta = Quaternion::dot(self, sign_flipped_to);
        // Lerp in the shorter direction.  This maybe should be optional.
        let sign_flip = Vector4::and(Vector4::set(SIGN_BIT), Vector4::from(cos_theta));

        let dest = Vector4::from(sign_flipped_to) * t_vec;
        let inv_t_vec = FloatVector::from(1.0f32) - t_vec;
        let f_vec = Vector4::xor(sign_flip, Vector4::from(inv_t_vec));
        let result = Quaternion::from(Vector4::mul_add(Vector4::from(self), f_vec, dest));

        return Quaternion::normalize(result);
    }

    /// Spherical linear interpolation.  This one in particular is probably slow.  
    /// It uses SIMD approximations for acos and sin, however, it could almost certainly be improved.
    #[inline(always)]
    pub fn slerp<T: Into<FloatVector>>(self, to: Quaternion, t: T) -> Quaternion {
        let t_tmp = t.into();
        let t_negative = Vector4::from(t_tmp).less(Vector4::zero());
        let sign_flipped_to = Quaternion::select(to, to.inverse(), t_negative);
        let t_vec = t_tmp.abs();

        let inv_t_vec = FloatVector::from(1.0f32) - t_vec;
        let cos_theta = Quaternion::dot(self, sign_flipped_to);
        let sign_flip = FloatVector::and(FloatVector::new(SIGN_BIT), cos_theta);
        let abs_cos_theta = FloatVector::xor(sign_flip, cos_theta);

        let lerp =
            t_vec == 0.0 || t_vec == 1.0 || abs_cos_theta.greater(FloatVector::new(SLERP_EPSILON));

        // If we are too close to parallel, switch to lerp.  Also if 1 is 0 or 1 so we can get an exact result.
        if lerp {
            let dest = Vector4::from(sign_flipped_to) * t_vec;
            let f_vec = Vector4::xor(Vector4::from(inv_t_vec), sign_flip);
            let result = Quaternion::from(Vector4::mul_add(Vector4::from(self), f_vec, dest));
            return Quaternion::normalize(result);
        }
        let theta = cos_theta.acos();
        let theta_scale = theta * Vector4::new(1.0f32, t_vec.value(), inv_t_vec.value(), 1.0f32);
        let sins = theta_scale.sin();
        let sin_div = sins / sins.x();

        let a_scale = FloatVector::xor(sign_flip, sin_div.z());
        let dest = Vector4::from(sign_flipped_to) * sin_div.y();
        let result = Vector4::mul_add(Vector4::from(self), Vector4::from(a_scale), dest);
        return Quaternion::from(result);
    }

    #[inline(always)]
    fn euler_xyz_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpackhi_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpackhi_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz);
            let wzyx = _mm_mul_ps(combined, sin_value.x().data);
            let wzyx_2 = _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, 0.0, SIGN_BIT), wzyx);

            let result = _mm_fmadd_ps(combined, cos_value.x().data, _wzyx(wzyx_2));
            return Quaternion { data: result };
        }
    }

    #[inline(always)]
    fn euler_xzy_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpackhi_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz);
            let wzyx = _mm_mul_ps(combined, sin_value.z().data);
            let wzyx_2 = _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, SIGN_BIT, 0.0), wzyx);

            let result = _mm_fmadd_ps(_zyxw(combined), cos_value.z().data, _yzwx(wzyx_2));
            return Quaternion { data: result };
        }
    }

    #[inline(always)]
    fn euler_yxz_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpackhi_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpackhi_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz);
            let wzyx = _mm_mul_ps(combined, sin_value.x().data);
            let wzyx_2 = _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, SIGN_BIT, 0.0), wzyx);

            let result = _mm_fmadd_ps(combined, cos_value.x().data, _wzyx(wzyx_2));
            return Quaternion { data: result };
        }
    }

    #[inline(always)]
    fn euler_yzx_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpackhi_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpackhi_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz);
            let wzyx = _mm_mul_ps(combined, sin_value.x().data);
            let wzyx_2 = _mm_xor_ps(_mm_set_ps(0.0, 0.0, SIGN_BIT, SIGN_BIT), wzyx);

            let result = _mm_fmadd_ps(combined, cos_value.x().data, _wzyx(wzyx_2));
            return Quaternion { data: result };
        }
    }
    #[inline(always)]
    fn euler_zxy_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpackhi_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sz_cz_b = _mm_xor_ps(_mm_set_ps(0.0, 0.0, 0.0, SIGN_BIT), sz_cz);
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpackhi_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz_b);
            let wzyx = _mm_mul_ps(combined, sin_value.x().data);

            let result = _mm_fmadd_ps(combined, cos_value.x().data, _wzyx(wzyx));
            return Quaternion { data: result };
        }
    }
    #[inline(always)]
    fn euler_zyx_normalized(sin_value: Vector3, cos_value: Vector3) -> Quaternion {
        unsafe {
            let sz_cz_0 = _mm_unpackhi_ps(sin_value.data, cos_value.data);
            let sz_cz = _mm_movelh_ps(sz_cz_0, sz_cz_0); //scsc
            let sscc_y_0 = _mm_unpacklo_ps(sin_value.data, cos_value.data); //sxcxsycy
            let sscc_y = _mm_unpacklo_ps(sscc_y_0, sscc_y_0); //sscc

            let combined = _mm_mul_ps(sscc_y, sz_cz);
            let wzyx = _mm_mul_ps(combined, sin_value.y().data);
            let wzyx_2 = _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, SIGN_BIT, 0.0), wzyx);

            let result = _mm_fmadd_ps(_yxzw(combined), cos_value.y().data, _zwyx(wzyx_2));
            return Quaternion { data: result };
        }
    }

    #[inline(always)]
    pub fn euler(radians: Vector3, rotation_order: RotationOrder) -> Quaternion {
        let half_angle = radians * 0.5;
        let sin = half_angle.sin();
        let cos = half_angle.cos();

        match rotation_order {
            RotationOrder::XYZ => return Quaternion::euler_xyz_normalized(sin, cos),
            RotationOrder::XZY => return Quaternion::euler_xzy_normalized(sin, cos),
            RotationOrder::YXZ => return Quaternion::euler_yxz_normalized(sin, cos),
            RotationOrder::YZX => return Quaternion::euler_yzx_normalized(sin, cos),
            RotationOrder::ZXY => return Quaternion::euler_zxy_normalized(sin, cos),
            RotationOrder::ZYX => return Quaternion::euler_zyx_normalized(sin, cos),
        }
    }

    /// Are these quaternions approximately the same, within ~0.163 degrees
    #[inline(always)]
    pub fn approx_equal(self, to: Quaternion) -> bool {
        let v = self.dot(to).abs().value();
        return v > EQUALITY_EPSILON_LOWER_BOUND && v < EQUALITY_EPSILON_UPPER_BOUND;
    }
}

impl From<Vector4> for Quaternion {
    #[inline(always)]
    fn from(v: Vector4) -> Quaternion {
        Quaternion { data: v.data }
    }
}

impl core::ops::Mul<Quaternion> for Quaternion {
    type Output = Quaternion;
    #[inline]
    fn mul(self, _rhs: Quaternion) -> Quaternion {
        return Quaternion::mul(self, _rhs);
    }
}
impl core::ops::MulAssign<Quaternion> for Quaternion {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: Quaternion) {
        *self = Quaternion::mul(*self, _rhs);
    }
}
impl core::ops::Mul<Vector3> for Quaternion {
    type Output = Vector3;
    #[inline]
    fn mul(self, _rhs: Vector3) -> Vector3 {
        return Quaternion::mul_vec(self, _rhs);
    }
}

impl PartialEq for Quaternion {
    #[inline(always)]
    fn eq(&self, other: &Quaternion) -> bool {
        return Quaternion::equal(*self, *other).all();
    }
}

#[cfg(test)]
mod test;
