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

const SLERP_EPSILON: f32 = 0.9995;
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

    // #[inline(always)]
    // pub fn store(self,  dst : &mut RawVector_f32){
    // 	let x : *mut f32 = &mut (dst.data[0]) as *mut f32;
    // 	unsafe{
    // 		_mm_store_ps(x, self.data);
    // 	}
    // }

    // #[inline(always)]
    // pub fn equals(v1 : Quaternion, v2 : Quaternion) -> bool{
    // 	unsafe{
    // 		let d = _mm_cmpeq_ps(v1.data, v2.data);
    // 		return (_mm_movemask_ps(d) ) == 15;
    // 	}
    // }
    #[inline(always)]
    pub fn dot(q1: Quaternion, q2: Quaternion) -> FloatVector {
        unsafe {
            return FloatVector {
                data: _ico_dp4_ps(q1.data, q2.data),
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
    pub fn shortest_delta(from: Quaternion, to: Quaternion) -> Quaternion {
        let dp = Quaternion::dot(from, to);
        let negative = Vector4::less(Vector4::from(dp), Vector4::zero());
        let sign_flip = Vector4::and(Vector4::set(SIGN_BIT), negative);
        let shortest_from = Vector4::xor(Vector4::from(from), sign_flip);
        let inv = Quaternion::inverse(Quaternion::from(shortest_from));
        return Quaternion::from(Vector4::mul(Vector4::from(inv), Vector4::from(to)));
        // unsafe{
        // 	let dp = _ico_dp4_ps(from.data, to.data);
        // 	let negative = _mm_cmplt_ps(dp, _mm_setzero_ps());

        // 	let flip_sign = _mm_and_ps(negative, _ico_signbit_ps());
        // 	let shortest_from = _mm_xor_ps(from.data, flip_sign);
        // 	return Quaternion{data : _mm_mul_ps(_mm_xor_ps(shortest_from, _mm_set_ps(0f32,SIGN_BIT,SIGN_BIT,SIGN_BIT)),to.data)};
        // }
    }

    #[inline(always)]
    pub fn angle(a: Quaternion, b: Quaternion) -> FloatVector {
        let dot = Quaternion::dot(a, b);
        let val = dot.acos();
        return val.add(val);
        // unsafe{
        // 	let w = _mm_cvtss_f32(_ico_dp4_ps(a.data, b.data));
        // 	return 2.0f32 * w.acos();
        // }
    }

    /// A fast normalize (q = q / q.magnitude).  When the length is 0, this will divide by 0 - and may produce infinity or NaN.
    /// Only should be used if the Quaternion is known to be non-zero.
    #[inline(always)]
    pub fn renormalize(q: Quaternion) -> Quaternion {
        let len = FloatVector::sqrt(Quaternion::dot(q, q));
        return Quaternion::from(Vector4::from(q) / len);
        // unsafe{
        // 		let sqr_length = _ico_dp4_ps(q.data, q.data);
        // 	return Quaternion{data : _mm_div_ps(q.data, _mm_sqrt_ps(sqr_length) )};
        // }
    }

    /// A safe normalize.  If the quaternion is 0, returns Identity.
    #[inline(always)]
    pub fn normalize(q: Quaternion) -> Quaternion {
        let len = FloatVector::sqrt(Quaternion::dot(q, q));
        let scaled = Vector4::from(q) / len;
        let mask = Vector4::less(scaled.abs(), Vector4::set(core::f32::INFINITY));
        return Quaternion::from(Vector4::select(
            Vector4::from(Quaternion::identity()),
            scaled,
            mask,
        ));
        // unsafe{
        // 		let sqr_length = _ico_dp4_ps(q.data, q.data);
        // 		//let mask = _mm_cmpgt_ps(sqr_length, _mm_set1_ps(NORMALIZATION_EPSILON));
        // 		let scaled = _mm_div_ps(q.data, _mm_sqrt_ps(sqr_length) );
        // 		// This will return false if the value is infinity, or NaN.
        // 		let mask = _mm_cmplt_ps(_ico_abs_ps(scaled), _mm_set1_ps(core::f32::INFINITY));
        // 	return Quaternion{data : _ico_select_ps(scaled,_mm_set_ps(1.0f32, 0.0f32, 0.0f32, 0.0f32), mask)};
        // }
    }

    #[inline(always)]
    pub fn lerp<T: Into<FloatVector>>(from: Quaternion, to: Quaternion, t: T) -> Quaternion {
        /*  Correct implementation
        float cosT = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
        ct = cosT;
        q2.x = q2.x * t;
        q2.y = q2.y * t;
        q2.z = q2.z * t;
        q2.w = q2.w * t;
        float inv_t = 1.0f - t;
        if(cosT < 0.0f){
            inv_t = -inv_t;
        }
        Quaternion result = q1;
        result.x = result.x * inv_t + q2.x;
        result.y = result.y * inv_t + q2.y;
        result.z = result.z * inv_t + q2.z;
        result.w = result.w * inv_t + q2.w;

        float mag = Mathf.Sqrt(
            result.x * result.x + result.y * result.y + result.z * result.z + result.w*result.w);
        result.x /= mag;
        result.y /= mag;
        result.z /= mag;
        result.w /= mag;
        return result;
        */
        let cos_theta = Quaternion::dot(from, to);
        let sign_flip = Vector4::and(Vector4::set(SIGN_BIT), Vector4::from(cos_theta));

        let t_vec = t.into();
        let dest = Vector4::from(to) * t_vec;
        let inv_t_vec = FloatVector::from(1.0f32) - t_vec;
        let f_vec = Vector4::xor(sign_flip, Vector4::from(inv_t_vec));
        let result = Quaternion::from(Vector4::mul_add(Vector4::from(from), f_vec, dest));
        return Quaternion::normalize(result);
        // unsafe{
        // let cos_theta = _ico_dp4_ps(from.data, to.data);

        // let t_vec = _mm_set1_ps(t);
        // let dest = _mm_mul_ps(to.data, t_vec);

        // 1-t
        // let inv_t_vec = _mm_sub_ps(_ico_one_ps(), t_vec);

        //flip sign based on shortest path
        //if cos is negative, negate the invBlend
        // let sign_flip = _mm_and_ps(_ico_signbit_ps(), cos_theta);
        // let f_vec = _mm_xor_ps(sign_flip, inv_t_vec);

        // let result =   Quaternion{data: _mm_fmadd_ps(from.data, f_vec, dest)};
        // return Quaternion::normalize(result);
    }

    #[inline(always)]
    pub fn slerp<T: Into<f32>>(from: Quaternion, to: Quaternion, t: T) -> Quaternion {
        let t_vec = t.into();
        let inv_t_vec = 1.0f32 - t_vec;
        let cos_theta = Quaternion::dot(from, to);
        let sign_flip = FloatVector::and(FloatVector::new(SIGN_BIT), cos_theta);
        let abs_cos_theta = FloatVector::xor(sign_flip, cos_theta);

        // If we are too close to parallel, switch to lerp.  Also if 1 is 0 or 1 so we can get an exact result.
        if (abs_cos_theta.value() > SLERP_EPSILON || t_vec == 0.0 || t_vec == 1.0 || t_vec == -1.0)
        {
            let dest = Vector4::from(to) * t_vec;
            let f_vec = Vector4::xor(Vector4::from(inv_t_vec), sign_flip);
            let result = Quaternion::from(Vector4::mul_add(Vector4::from(from), f_vec, dest));
            return Quaternion::normalize(result);
        }
        let theta = cos_theta.acos();
        let theta_scale = theta * Vector4::new(1.0f32, t_vec, inv_t_vec, 1.0f32);
        let sins = theta_scale.sin();
        let sin_div = sins / sins.x();

        let a_scale = FloatVector::xor(sign_flip, sin_div.z());
        let dest = Vector4::from(to) * sin_div.y();
        let result = Vector4::mul_add(Vector4::from(from), Vector4::from(a_scale), dest);
        return Quaternion::from(result);
        // unsafe{
        //    	let cos_theta = _ico_dp4_ps(from.data, to.data);
        //    	let sign_flip = _mm_and_ps(_ico_signbit_ps(), cos_theta);
        //    	let abs_cos_theta = _mm_xor_ps(sign_flip, cos_theta);

        //    	// If we are too close to parallel, switch to lerp
        //    	if(_mm_cvtss_f32(abs_cos_theta) > SLERP_EPSILON){
        //    		let t_vec = _mm_set1_ps(t);
        //    		let inv_t_vec = _mm_sub_ps(_ico_one_ps(), t_vec);
        //    		let f_vec = _mm_xor_ps(sign_flip, inv_t_vec);
        // 		let dest = _mm_mul_ps(to.data, t_vec);
        // 		let result =   Quaternion{data: _mm_fmadd_ps(from.data, f_vec, dest)};
        // 		return Quaternion::normalize(result);
        //    	}

        //    	let cos_theta_f = _mm_cvtss_f32(cos_theta);
        //    	let theta : f32 = cos_theta_f.acos();

        //    	//TODO: use sse for fast sines all at once
        // 	let sin_theta : f32 = (1.0 - cos_theta_f*cos_theta_f).sqrt();//theta.sin();
        // 	let a = ((1.0f32-t) * theta).sin() / sin_theta;
        // 	let b = _mm_set1_ps((t * theta).sin() / sin_theta);

        // 	// again this uses the flipped sign.
        // 	let a_scale = _mm_xor_ps(sign_flip, _mm_set1_ps(a));
        // 	let dest = _mm_mul_ps(to.data, b);
        // 	let result =   Quaternion{data: _mm_fmadd_ps(from.data, a_scale, dest)};
        // 	return result;
        // }
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
