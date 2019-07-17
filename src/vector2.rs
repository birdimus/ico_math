use crate::Vector2;
use crate::Vector3;
use crate::Vector4;
use crate::Vector2Int;

use std::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_abs_ps;
use crate::_ico_truncate_ps;
use crate::_ico_copysign_ps;
use crate::_ico_two_ps;
use crate::_ico_signbit_ps;
use crate::*;
impl Vector2{
	/// Returns a new Vector2
	#[inline(always)]
	pub fn new(x : f32, y : f32) -> Vector2{
		unsafe{
			Vector2{data : _mm_set_ps(0.0f32, 0.0f32, y, x)}
		}
	}
	#[inline(always)]
	pub fn set(value : f32) -> Vector2{
		unsafe{
			Vector2{data : _mm_set1_ps(value)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector2 {
		unsafe{
			Vector2 { data : _mm_setzero_ps() }
		}
	}

	#[inline(always)]
	pub fn x(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.yy().data);
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
	pub fn add(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_add_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_sub_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_mul(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_mul_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_div(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_div_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn fmadd(v1 : Vector2, v2 : Vector2, v3 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_fmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fmsub(v1 : Vector2, v2 : Vector2, v3 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_fmsub_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmadd(v1 : Vector2, v2 : Vector2, v3 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_fnmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmsub(v1 : Vector2, v2 : Vector2, v3 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_fnmsub_ps(v1.data, v2.data, v3.data)}
		}
	}

	#[inline(always)]
	pub fn scale(v1 : Vector2, scalar : f32) -> Vector2{	
		unsafe{
			Vector2{data : _mm_mul_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn div(v1 : Vector2, scalar : f32) -> Vector2{	
		unsafe{
			Vector2{data : _mm_div_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn and(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_and_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn or(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_or_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn andnot(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_andnot_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn xor(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_xor_ps(v1.data, v2.data)}
		}
	}
	
	#[inline(always)]
	pub fn all(v1 : Vector2) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 3) == 3;
		}
	}

	#[inline(always)]
	pub fn any(v1 : Vector2) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 3) != 0;
		}
	}

	#[inline(always)]
	pub fn equals(v1 : Vector2, v2 : Vector2) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.data, v2.data);
			return (_mm_movemask_ps(d) & 3) == 3;
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmpeq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_not_equal(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmpneq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_greater_equal(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmpge_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmpgt_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less_equal(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmple_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_cmplt_ps(v1.data, v2.data)}
		}
	}


	#[inline(always)]
	pub fn abs(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  _ico_abs_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  _ico_copysign_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	/// Floor function.  Returns signed 0 when applicable.
	pub fn floor(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  
				_mm_floor_ps(v1.data)}
				//_ico_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Ceil function.  Returns signed 0 when applicable.
	pub fn ceil(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  
				_mm_ceil_ps(v1.data)}
				//_ico_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Round to nearest even function. Returns signed 0 when applicable.
	pub fn round(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : 
			_mm_round_ps(v1.data, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC)}
			// _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : Vector2) -> Vector2Int{	
		unsafe{
			Vector2Int{data :  _mm_cvttps_epi32(_mm_floor_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : Vector2) -> Vector2Int{	
		unsafe{
			Vector2Int{data :  _mm_cvttps_epi32(_mm_ceil_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  
				_mm_round_ps(v1.data, _MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC)}
				//_ico_truncate_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn frac(v1 : Vector2) -> Vector2{	
		return Vector2::sub(v1, Vector2::floor(v1));
	}

	#[inline(always)]
	pub fn sqrt(v1 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn max(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector2, v2 : Vector2) -> Vector2{	
		unsafe{
			Vector2{data : _mm_min_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn dot(v0 : Vector2, v1 : Vector2) -> Vector2{	
		unsafe{
			let tmp0 = _mm_mul_ps(v0.data, v1.data);
			let mut tmp1 = _yxzw(tmp0);// _mm_shuffle_ps(tmp0,tmp0, _ico_shuffle(3, 2, 0, 1)); //yxzw
			
			tmp1 = _mm_add_ps(tmp0, tmp1);//xy,xy,qq,qq
			return Vector2{data : _mm_unpacklo_ps(tmp1,tmp1)};//xy,xy.xy,xy
		}
	}
	#[inline(always)]
	pub fn lerp(v1 : Vector2, v2 : Vector2, t : f32) -> Vector2{	
		unsafe{
			let t_val = _mm_set1_ps(t);
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector2{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}

	#[inline(always)]
	pub fn normalize(v1 : Vector2) -> Vector2{	
		let length = Vector2::sqrt(Vector2::dot(v1,v1));
		let norm = Vector2::component_div(v1, length);
		
		unsafe{
			// This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
			let result_length_sqr = Vector2::dot(norm,norm);
			let mask_less = Vector2::component_less( result_length_sqr, Vector2{data :_ico_two_ps()});
			return Vector2::and(norm, mask_less);
		}
	}

	#[inline(always)]
	pub fn normalize_length(v1 : Vector2) -> (Vector2, f32){	
		let length = Vector2::sqrt(Vector2::dot(v1,v1));
		let norm = Vector2::component_div(v1, length);
		
		unsafe{

			// This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
			let result_length_sqr = Vector2::dot(norm,norm);
			let mask_less = Vector2::component_less( result_length_sqr, Vector2{data :_ico_two_ps()});
			return (Vector2::and(norm, mask_less), length.x());
		}
	}

	

	#[inline(always)]
	pub fn sqr_magnitude(self) -> f32 {
		return Vector2::dot(self, self).x();	
	}
	#[inline(always)]
	pub fn magnitude(self) -> f32 {
		return Vector2::sqrt(Vector2::dot(self, self)).x();	
	}

	#[inline(always)]
	/// Right handed system, positive rotation is counterclockwise about the axis of rotation.
	pub fn rotate(v1 : Vector2, radians : f64) -> Vector2 {
		let f = radians.sin_cos();

		let mut sn : f32 = f.0 as f32;
		let mut cs : f32 = f.1 as f32;

		unsafe{

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
			return Vector2{data: _mm_add_ps(v2,_wzyx(v2))};//wzyx
			
		}
	}

	#[inline(always)] pub fn xxxx(self) -> Vector4 { unsafe{return Vector4{data:_xxxx(self.data)};}}
	#[inline(always)] pub fn yyyy(self) -> Vector4 { unsafe{return Vector4{data:_yyyy(self.data)};}}

	#[inline(always)] pub fn xx(self) -> Vector2 { unsafe{return Vector2{data:_xxzw(self.data)};}}
	#[inline(always)] pub fn xy(self) -> Vector2 { unsafe{return Vector2{data:_xyzw(self.data)};}}
	#[inline(always)] pub fn yx(self) -> Vector2 { unsafe{return Vector2{data:_yxzw(self.data)};}}
	#[inline(always)] pub fn yy(self) -> Vector2 { unsafe{return Vector2{data:_yyzw(self.data)};}}
}


impl From<Vector3> for Vector2 {
    fn from(v : Vector3) -> Vector2 {
        Vector2 { data : v.data }
    }
}
impl From<Vector4> for Vector2 {
    fn from(v : Vector4) -> Vector2 {
        Vector2 { data : v.data }
    }
}
impl From<Vector2Int> for Vector2 {
    fn from(v : Vector2Int) -> Vector2 {
    	unsafe{
        	return Vector2 { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl std::ops::Add for Vector2{
	type Output = Vector2;
	#[inline]
	fn add(self, _rhs: Vector2) -> Vector2{
		Vector2::add(self, _rhs)
	}
}
impl std::ops::AddAssign for Vector2 {
	#[inline(always)]
    fn add_assign(&mut self, other: Vector2) {
        *self = Vector2::add(*self, other)
    }
}
impl std::ops::Sub for Vector2{
	type Output = Vector2;
	#[inline]
	fn sub(self, _rhs: Vector2) -> Vector2{
		Vector2::sub(self, _rhs)
	}
}
impl std::ops::SubAssign for Vector2 {
	#[inline(always)]
    fn sub_assign(&mut self, other: Vector2) {
        *self = Vector2::sub(*self, other)
    }
}
impl std::ops::Neg for Vector2 {
	type Output = Vector2;
	#[inline(always)]
	fn neg(self) -> Self::Output {
		unsafe{
			return Vector2{data:_mm_xor_ps(_ico_signbit_ps(),self.data)};
		}
	}
}
impl std::ops::Mul<f32> for Vector2{
	type Output = Vector2;
	#[inline]
	fn mul(self, _rhs: f32) -> Vector2{
		Vector2::scale(self, _rhs)
	}
}

impl std::ops::Mul<Vector2> for f32{
	type Output = Vector2;
	#[inline]
	fn mul(self, _rhs: Vector2) -> Vector2{
		Vector2::scale(_rhs, self)
	}
}
impl std::ops::MulAssign<f32> for Vector2{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: f32){
		*self = Vector2::scale(*self, _rhs)
	}
}
impl std::ops::Div<f32> for Vector2{
	type Output = Vector2;
	#[inline]
	fn div(self, _rhs: f32) -> Vector2{
		Vector2::div(self, _rhs)
	}
}
impl std::ops::Div<Vector2> for f32{
	type Output = Vector2;
	#[inline(always)]
	fn div(self : f32, _rhs: Vector2) -> Vector2{
		return Vector2::component_div(Vector2::set(self), _rhs);
	}
}
impl std::ops::DivAssign<f32> for Vector2{
	#[inline(always)]
	fn div_assign(&mut self, _rhs: f32){
		*self = Vector2::div(*self, _rhs)
	}
}
impl PartialEq for Vector2 {
    fn eq(&self, other: &Vector2) -> bool {
    	return Vector2::equals(*self, *other);
    }
}
#[cfg(test)]
mod tests;