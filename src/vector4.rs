use crate::Vector2;
use crate::Vector3;
use crate::Vector4;
use crate::Vector4Int;

use std::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_floorps_epi32;
use crate::_ico_ceilps_epi32;
use crate::_ico_round_ps;
use crate::_ico_abs_ps;
use crate::_ico_truncate_ps;
use crate::_ico_copysign_ps;
use crate::_ico_dp4_ps;
impl Vector4{
	/// Returns a new Vector4
	#[inline(always)]
	pub fn new(x : f32, y : f32, z : f32, w : f32) -> Vector4{
		unsafe{
			Vector4{data : _mm_set_ps(w, z, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector4 {
		unsafe{
			Vector4 { data : _mm_setzero_ps() }
		}
	}

	#[inline(always)]
	pub fn x(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(1, 1, 1, 1)));
		}	
	}

	#[inline(always)]
	pub fn z(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(2, 2, 2, 2)));
		}	
	}

	#[inline(always)]
	pub fn w(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(3,3,3,3)));
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
	pub fn add(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_add_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_sub_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_mul(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_mul_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_div(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_div_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn fmadd(v1 : Vector4, v2 : Vector4, v3 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_fmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fmsub(v1 : Vector4, v2 : Vector4, v3 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_fmsub_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmadd(v1 : Vector4, v2 : Vector4, v3 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_fnmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmsub(v1 : Vector4, v2 : Vector4, v3 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_fnmsub_ps(v1.data, v2.data, v3.data)}
		}
	}

	#[inline(always)]
	pub fn scale(v1 : Vector4, scalar : f32) -> Vector4{	
		unsafe{
			Vector4{data : _mm_mul_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn div(v1 : Vector4, scalar : f32) -> Vector4{	
		unsafe{
			Vector4{data : _mm_div_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn and(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_and_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn or(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_or_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn andnot(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_andnot_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn xor(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_xor_ps(v1.data, v2.data)}
		}
	}
	
	#[inline(always)]
	pub fn all(v1 : Vector4) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) ) == 15;
		}
	}

	#[inline(always)]
	pub fn any(v1 : Vector4) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) ) != 0;
		}
	}

	#[inline(always)]
	pub fn equals(v1 : Vector4, v2 : Vector4) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.data, v2.data);
			return (_mm_movemask_ps(d) ) == 15;
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmpeq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_not_equal(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmpneq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_greater_equal(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmpge_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmpgt_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less_equal(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmple_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_cmplt_ps(v1.data, v2.data)}
		}
	}


	#[inline(always)]
	pub fn abs(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _ico_abs_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _ico_copysign_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn floor(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _mm_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn ceil(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _mm_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn round(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : Vector4) -> Vector4Int{	
		unsafe{
			Vector4Int{data :  _ico_floorps_epi32(v1.data)}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : Vector4) -> Vector4Int{	
		unsafe{
			Vector4Int{data :  _ico_ceilps_epi32(v1.data)}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _ico_truncate_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn frac(v1 : Vector4) -> Vector4{	
		return Vector4::sub(v1, Vector4::floor(v1));
	}

	#[inline(always)]
	pub fn sqrt(v1 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn max(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			Vector4{data : _mm_min_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn dot(v1 : Vector4, v2 : Vector4) -> Vector4{	
		unsafe{
			return Vector4{data : _ico_dp4_ps(v1.data,v2.data)};
		}
	}
	#[inline(always)]
	pub fn lerp(v1 : Vector4, v2 : Vector4, t : f32) -> Vector4{	
		unsafe{
			let t_val = _mm_set1_ps(t);
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector4{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}

	#[inline(always)]
	pub fn renormalize(v1 : Vector4) -> Vector4{	
		let length = Vector4::sqrt(Vector4::dot(v1,v1));
		return Vector4::component_div(v1, length);
	}

	#[inline(always)]
	pub fn normalize(v1 : Vector4) -> Vector4{	
		let length = Vector4::sqrt(Vector4::dot(v1,v1));
		let norm = Vector4::component_div(v1, length);
		let mask = Vector4::component_less(Vector4::abs(norm), Vector4::from( std::f32::INFINITY));
		return Vector4::and(norm, mask);
	}

	#[inline(always)]
	pub fn normalize_length(v1 : Vector4) -> (Vector4, f32){	
		let length = Vector4::sqrt(Vector4::dot(v1,v1));
		let norm = Vector4::component_div(v1, length);
		let mask = Vector4::component_less(Vector4::abs(norm),Vector4::from( std::f32::INFINITY));
		return (Vector4::and(norm, mask), length.x());
	}

	

	#[inline(always)]
	pub fn sqr_magnitude(&self) -> f32 {
		return Vector4::dot(*self, *self).x();	
	}
	#[inline(always)]
	pub fn magnitude(&self) -> f32 {
		return Vector4::sqrt(Vector4::dot(*self, *self)).x();	
	}

	
}

impl From<f32> for Vector4 {
    fn from(val : f32) -> Vector4 {
       unsafe{
			Vector4{data : _mm_set1_ps(val)}
		}
    }
}
impl From<Vector2> for Vector4 {
    fn from(v : Vector2) -> Vector4 {
        Vector4 { data : v.data }
    }
}
impl From<Vector3> for Vector4 {
    fn from(v : Vector3) -> Vector4 {
        Vector4 { data : v.data }
    }
}
impl From<Vector4Int> for Vector4 {
    fn from(v : Vector4Int) -> Vector4 {
    	unsafe{
        	return Vector4 { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl std::ops::Add for Vector4{
	type Output = Vector4;
	#[inline]
	fn add(self, _rhs: Vector4) -> Vector4{
		Vector4::add(self, _rhs)
	}
}

impl std::ops::Sub for Vector4{
	type Output = Vector4;
	#[inline]
	fn sub(self, _rhs: Vector4) -> Vector4{
		Vector4::sub(self, _rhs)
	}
}

impl std::ops::Mul<f32> for Vector4{
	type Output = Vector4;
	#[inline]
	fn mul(self, _rhs: f32) -> Vector4{
		Vector4::scale(self, _rhs)
	}
}

impl std::ops::Mul<Vector4> for f32{
	type Output = Vector4;
	#[inline]
	fn mul(self, _rhs: Vector4) -> Vector4{
		Vector4::scale(_rhs, self)
	}
}

impl std::ops::Div<f32> for Vector4{
	type Output = Vector4;
	#[inline]
	fn div(self, _rhs: f32) -> Vector4{
		Vector4::div(self, _rhs)
	}
}
	
impl PartialEq for Vector4 {
    fn eq(&self, other: &Vector4) -> bool {
    	return Vector4::equals(*self, *other);
    }
}
