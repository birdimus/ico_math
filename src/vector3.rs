use crate::Vector2;
use crate::Vector3;
use crate::Vector4;
use crate::Vector3Int;

use std::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_floorps_epi32;
use crate::_ico_ceilps_epi32;
use crate::_ico_round_ps;
use crate::_ico_abs_ps;
use crate::_ico_cross_ps;
use crate::_ico_truncate_ps;
use crate::_ico_copysign_ps;
impl Vector3{
	/// Returns a new Vector3
	#[inline(always)]
	pub fn new(x : f32, y : f32, z : f32) -> Vector3{
		unsafe{
			Vector3{data : _mm_set_ps(0.0f32, z, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector3 {
		unsafe{
			Vector3 { data : _mm_setzero_ps() }
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
	pub fn add(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_add_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_sub_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_mul(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_mul_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_div(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_div_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn fmadd(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fmsub(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fmsub_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmadd(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fnmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmsub(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fnmsub_ps(v1.data, v2.data, v3.data)}
		}
	}

	#[inline(always)]
	pub fn scale(v1 : Vector3, scalar : f32) -> Vector3{	
		unsafe{
			Vector3{data : _mm_mul_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn div(v1 : Vector3, scalar : f32) -> Vector3{	
		unsafe{
			Vector3{data : _mm_div_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn and(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_and_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn or(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_or_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn andnot(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_andnot_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn xor(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_xor_ps(v1.data, v2.data)}
		}
	}
	
	#[inline(always)]
	pub fn all(v1 : Vector3) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 7) == 7;
		}
	}

	#[inline(always)]
	pub fn any(v1 : Vector3) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 7) != 7;
		}
	}

	#[inline(always)]
	pub fn equals(v1 : Vector3, v2 : Vector3) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.data, v2.data);
			return (_mm_movemask_ps(d) & 7) == 7;
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpeq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_not_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpneq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_greater_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpge_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpgt_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmple_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmplt_ps(v1.data, v2.data)}
		}
	}


	#[inline(always)]
	pub fn cross(lhs : Vector3, rhs : Vector3) -> Vector3{	
		unsafe{
			return Vector3{data :  _ico_cross_ps(lhs.data, rhs.data)};
		}
	}

	#[inline(always)]
	pub fn abs(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_abs_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_copysign_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn floor(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _mm_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn ceil(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _mm_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn round(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _ico_floorps_epi32(v1.data)}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _ico_ceilps_epi32(v1.data)}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_truncate_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn frac(v1 : Vector3) -> Vector3{	
		return Vector3::sub(v1, Vector3::floor(v1));
	}

	#[inline(always)]
	pub fn sqrt(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn max(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_min_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn dot(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{

			let tmp0 = _mm_mul_ps(v1.data, v2.data); //xyzw
		    let tmp1 = _mm_castsi128_ps(_mm_slli_si128 (_mm_castps_si128(tmp0), 4)); //0xyz
		    let tmp2 = _mm_add_ps(tmp0 , tmp1);//x xy, yz, wz
		    let tmp3 = _mm_moveldup_ps(tmp2); // x x yz yz
			Vector3{data : _mm_add_ps(tmp3, _mm_shuffle_ps(tmp3,tmp3, _ico_shuffle(0, 1, 2, 3)))}
		}
	}
	#[inline(always)]
	pub fn lerp(v1 : Vector3, v2 : Vector3, t : f32) -> Vector3{	
		unsafe{
			let t_val = _mm_set1_ps(t);
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector3{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}

	#[inline(always)]
	pub fn normalize(v1 : Vector3) -> Vector3{	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		return Vector3::component_div(v1, length);
	}

	#[inline(always)]
	pub fn normalize_length(v1 : Vector3) -> (Vector3, f32){	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		return (Vector3::component_div(v1, length), length.x());
	}

	

	#[inline(always)]
	pub fn sqr_magnitude(&self) -> f32 {
		return Vector3::dot(*self, *self).x();	
	}
	#[inline(always)]
	pub fn magnitude(&self) -> f32 {
		return Vector3::sqrt(Vector3::dot(*self, *self)).x();	
	}

	
}

impl From<f32> for Vector3 {
    fn from(val : f32) -> Vector3 {
       unsafe{
			Vector3{data : _mm_set1_ps(val)}
		}
    }
}
impl From<Vector2> for Vector3 {
    fn from(v : Vector2) -> Vector3 {
        Vector3 { data : v.data }
    }
}
impl From<Vector4> for Vector3 {
    fn from(v : Vector4) -> Vector3 {
        Vector3 { data : v.data }
    }
}
impl From<Vector3Int> for Vector3 {
    fn from(v : Vector3Int) -> Vector3 {
    	unsafe{
        	return Vector3 { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl std::ops::Add for Vector3{
	type Output = Vector3;
	#[inline]
	fn add(self, _rhs: Vector3) -> Vector3{
		Vector3::add(self, _rhs)
	}
}

impl std::ops::Sub for Vector3{
	type Output = Vector3;
	#[inline]
	fn sub(self, _rhs: Vector3) -> Vector3{
		Vector3::sub(self, _rhs)
	}
}

impl std::ops::Mul<f32> for Vector3{
	type Output = Vector3;
	#[inline]
	fn mul(self, _rhs: f32) -> Vector3{
		Vector3::scale(self, _rhs)
	}
}

impl std::ops::Mul<Vector3> for f32{
	type Output = Vector3;
	#[inline]
	fn mul(self, _rhs: Vector3) -> Vector3{
		Vector3::scale(_rhs, self)
	}
}

impl std::ops::Div<f32> for Vector3{
	type Output = Vector3;
	#[inline]
	fn div(self, _rhs: f32) -> Vector3{
		Vector3::div(self, _rhs)
	}
}
	
impl PartialEq for Vector3 {
    fn eq(&self, other: &Vector3) -> bool {
    	return Vector3::equals(*self, *other);
    }
}
