use core::arch::x86_64::*;
use crate::Vector3;
use crate::Quaternion;
use crate::RawVec;
use crate::_ico_shuffle;
use crate::_ico_quat_mul;
use crate::_ico_abs_ps;
use crate::_ico_dp4_ps;
use crate::_ico_signbit_ps;
use crate::_ico_select_ps;
use crate::_ico_one_ps;
use crate::NORMALIZATION_EPSILON;
use crate::SIGN_BIT;
const SLERP_EPSILON : f32 = 0.9995;
impl Quaternion{

	#[inline(always)]
	pub fn identity() -> Quaternion {
		unsafe{
			Quaternion { data : _mm_set_ps(1.0f32, 0.0f32, 0.0f32, 0.0f32) }
		}
	}


	#[inline(always)]
	pub fn mul_vec(quat : Quaternion, vec : Vector3) -> Vector3{	
		unsafe{
			// t = 2.0 * cross(q.xyz, vec)
			let qv = Vector3{data: quat.data};
			let t = Vector3::cross(qv,vec);
			let t2 = Vector3::add(t,t);
		
			// v' = v + q.w * t + cross(q.xyz, t)
			let w = Vector3{data:_mm_shuffle_ps(quat.data, quat.data, _ico_shuffle(3, 3, 3, 3))};
			let tmp = Vector3::fmadd(w, t2, vec);
	        
			return Vector3::add(tmp, Vector3::cross(qv, t2));
		}
	}
	#[inline(always)]
	pub fn mul(lhs : Quaternion, rhs : Quaternion) -> Quaternion{	
		unsafe{
			return Quaternion{data : _ico_quat_mul(lhs.data, rhs.data)};
		}
	}	
	#[inline(always)]
	pub fn angle_axis(axis : Vector3, radians : f32) -> Quaternion{	
		let half_angle = radians *0.5f32;

		let sin_cos = half_angle.sin_cos();
		unsafe{
			let mut q = Quaternion{data : _mm_mul_ps(axis.data, _mm_set1_ps(sin_cos.0))};
		
		q.set_w(sin_cos.1);
		return q;
		}
	}



	#[inline(always)]
	pub fn store(self,  dst : &mut RawVec){	
		let x : *mut f32 = &mut (dst.data[0]) as *mut f32;
		unsafe{
			_mm_store_ps(x, self.data);
		}
	}

	//#[inline(always)]
	//pub fn to_angle_axis(self) -> (Vector3, f32){	

	//}

	#[inline(always)]
	pub fn x(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(1, 1, 1, 1)));
		}	
	}

	#[inline(always)]
	pub fn z(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(2, 2, 2, 2)));
		}	
	}

	#[inline(always)]
	pub fn w(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(3, 3, 3, 3)));
		}	
	}

	#[inline(always)]
	pub fn set_x(&mut self, value : f32) {
		unsafe{
			self.data = _mm_move_ss(self.data, _mm_set_ss(value));
		}	
	}

	#[inline(always)]
	pub fn set_y(&mut self, value : f32) {
		unsafe{
			let v1 = _mm_move_ss(_mm_set1_ps(value), self.data);
			self.data = _mm_shuffle_ps(v1,self.data, _ico_shuffle(3, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn set_z(&mut self, value : f32) {
		unsafe{
			self.data = _mm_shuffle_ps(self.data, _mm_set1_ps(value), _ico_shuffle(3, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn set_w(&mut self, value : f32) {
		unsafe{
			let v1 = _mm_move_ss(self.data, _mm_set_ss(value));
			self.data = _mm_shuffle_ps(self.data, v1,  _ico_shuffle(0, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn equals(v1 : Quaternion, v2 : Quaternion) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.data, v2.data);
			return (_mm_movemask_ps(d) ) == 15;
		}
	}
	#[inline(always)]
	pub fn dot(q1 : Quaternion, q2 : Quaternion) -> Quaternion{	
		unsafe{
			return Quaternion{data : _ico_dp4_ps(q1.data,q2.data)};
		}
	}
	#[inline(always)]
	pub fn sqr_magnitude(self) -> f32 {
		return Quaternion::dot(self, self).x();	
	}
	#[inline(always)]
	pub fn magnitude(self) -> f32 {
		unsafe{
		return _mm_cvtss_f32(_mm_sqrt_ps(_ico_dp4_ps(self.data, self.data)));	
		}	
	}

	#[inline(always)]
	pub fn inverse(self) -> Quaternion {
		unsafe{
		return Quaternion{data : _mm_xor_ps(self.data, _mm_set_ps(0f32, SIGN_BIT,SIGN_BIT,SIGN_BIT))};
		}
	}

	#[inline(always)]
	pub fn reverse(self) ->Quaternion{
		unsafe{
			return Quaternion{data : _mm_xor_ps(self.data,_ico_signbit_ps())};
		}
	}
	#[inline(always)]
	pub fn delta(from : Quaternion, to : Quaternion) ->Quaternion{
		return Quaternion::mul(from.inverse(), to);
	}
	#[inline(always)]
	pub fn shortest_delta(from : Quaternion, to : Quaternion) -> Quaternion{

		unsafe{
			let dp = _ico_dp4_ps(from.data, to.data);
			let negative = _mm_cmplt_ps(dp, _mm_setzero_ps());

			let flip_sign = _mm_and_ps(negative, _ico_signbit_ps());
			let shortest_from = _mm_xor_ps(from.data, flip_sign);
			return Quaternion{data : _mm_mul_ps(_mm_xor_ps(shortest_from, _mm_set_ps(0f32,SIGN_BIT,SIGN_BIT,SIGN_BIT)),to.data)};
		}
	
	}
	#[inline(always)]
	pub fn angle( a : Quaternion,  b : Quaternion)->f32{
		unsafe{
			let w = _mm_cvtss_f32(_ico_dp4_ps(a.data, b.data));
			return 2.0f32 * w.acos();
		}
	}

	/// A fast normalize (q = q / q.magnitude).  When the length is 0, this will divide by 0 - and may produce infinity or NaN.
	/// Only should be used if the Quaternion is known to be non-zero.
	#[inline(always)]
	pub fn renormalize( q : Quaternion) ->Quaternion{
		unsafe{
	 		let sqr_length = _ico_dp4_ps(q.data, q.data);
			return Quaternion{data : _mm_div_ps(q.data, _mm_sqrt_ps(sqr_length) )};
		}
	}


	/// A safe normalize.  If the quaternion is 0, returns Identity.
	#[inline(always)]
	pub fn normalize(q : Quaternion)->Quaternion{
		unsafe{
	 		let sqr_length = _ico_dp4_ps(q.data, q.data);
	 		//let mask = _mm_cmpgt_ps(sqr_length, _mm_set1_ps(NORMALIZATION_EPSILON));
	 		let scaled = _mm_div_ps(q.data, _mm_sqrt_ps(sqr_length) );
	 		// This will return false if the value is infinity, or NaN.
	 		let mask = _mm_cmplt_ps(_ico_abs_ps(scaled), _mm_set1_ps(core::f32::INFINITY));
			return Quaternion{data : _ico_select_ps(scaled,_mm_set_ps(1.0f32, 0.0f32, 0.0f32, 0.0f32), mask)};
		}
	}




	#[inline(always)]
	pub fn lerp(from : Quaternion, to : Quaternion, t : f32)->Quaternion{

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
	unsafe{
	    let cos_theta = _ico_dp4_ps(from.data, to.data);

	    let t_vec = _mm_set1_ps(t);
	    let dest = _mm_mul_ps(to.data, t_vec);
	    
	    // 1-t
		let inv_t_vec = _mm_sub_ps(_ico_one_ps(), t_vec);

		//flip sign based on shortest path
		//if cos is negative, negate the invBlend
	    let sign_flip = _mm_and_ps(_ico_signbit_ps(), cos_theta);
	    let f_vec = _mm_xor_ps(sign_flip, inv_t_vec);
	 
	    let result =   Quaternion{data: _mm_fmadd_ps(from.data, f_vec, dest)};
	    return Quaternion::normalize(result);
	}



	#[inline(always)]
	pub fn slerp(from : Quaternion, to : Quaternion, t : f32)->Quaternion{
		unsafe{
	    	let cos_theta = _ico_dp4_ps(from.data, to.data);
	    	let sign_flip = _mm_and_ps(_ico_signbit_ps(), cos_theta);
	    	let abs_cos_theta = _mm_xor_ps(sign_flip, cos_theta);

	    	// If we are too close to parallel, switch to lerp
	    	if(_mm_cvtss_f32(abs_cos_theta) > SLERP_EPSILON){
	    		let t_vec = _mm_set1_ps(t);
	    		let inv_t_vec = _mm_sub_ps(_ico_one_ps(), t_vec);
	    		let f_vec = _mm_xor_ps(sign_flip, inv_t_vec);
				let dest = _mm_mul_ps(to.data, t_vec);
				let result =   Quaternion{data: _mm_fmadd_ps(from.data, f_vec, dest)};
				return Quaternion::normalize(result);
	    	}

	    	let cos_theta_f = _mm_cvtss_f32(cos_theta);
	    	let theta : f32 = cos_theta_f.acos();

	    	//TODO: use sse for fast sines all at once
			let sin_theta : f32 = (1.0 - cos_theta_f*cos_theta_f).sqrt();//theta.sin();
			let a = ((1.0f32-t) * theta).sin() / sin_theta;
			let b = _mm_set1_ps((t * theta).sin() / sin_theta);

			// again this uses the flipped sign.
			let a_scale = _mm_xor_ps(sign_flip, _mm_set1_ps(a));
			let dest = _mm_mul_ps(to.data, b);
			let result =   Quaternion{data: _mm_fmadd_ps(from.data, a_scale, dest)};
			return result;
		}
	}


	//#[inline(always)]
	//pub fn look_rotation(direction : Vector3, up : Vector3)->Quaternion{


	//}
 
}
	
}
	


impl core::ops::Mul<Quaternion> for Quaternion{
	type Output = Quaternion;
	#[inline]
	fn mul(self, _rhs: Quaternion) -> Quaternion{
		return Quaternion::mul(self, _rhs);
	}
}

impl core::ops::Mul<Vector3> for Quaternion{
	type Output = Vector3;
	#[inline]
	fn mul(self, _rhs: Vector3) -> Vector3{
		return Quaternion::mul_vec(self, _rhs);
	}
}
impl PartialEq for Quaternion {
    fn eq(&self, other: &Quaternion) -> bool {
    	return Quaternion::equals(*self, *other);
    }
}