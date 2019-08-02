use crate::FloatVector;
use crate::Vector2;
use crate::Vector3;
use crate::Vector4;
use crate::Vector3Int;

use core::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_abs_ps;
use crate::_ico_cross_ps;
use crate::_ico_copysign_ps;
use crate::_ico_two_ps;
use crate::_ico_signbit_ps;
use crate::*;
impl Vector3{
	/// Returns a new Vector3
	#[inline(always)]
	pub fn new(x : f32, y : f32, z : f32) -> Vector3{
		unsafe{
			Vector3{data : _mm_set_ps(0.0f32, z, y, x)}
		}
	}
	
	#[inline(always)]
	pub fn set<T : Into<FloatVector>>(value : T) -> Vector3 {
		return Vector3{data : value.into().data};
	}

	#[inline(always)]
	pub fn zero() -> Vector3 {
		unsafe{
			Vector3 { data : _mm_setzero_ps() }
		}
	}

	#[inline(always)]
	pub fn x(self) -> FloatVector {
		return FloatVector{data:self.xxxx().data};
	}

	#[inline(always)]
	pub fn y(self) -> FloatVector {
		return FloatVector{data:self.yyyy().data};
	}

	#[inline(always)]
	pub fn z(self) -> FloatVector {
		return FloatVector{data:self.zzzz().data};
	}

	#[inline(always)]
	pub fn set_x<T : Into<FloatVector>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_move_ss(self.data, value.into().data);
		}	
	}

	#[inline(always)]
	pub fn set_y<T : Into<FloatVector>>(&mut self, value : T) {
		unsafe{
			let v1 = _mm_move_ss(value.into().data, self.data);
			self.data = _mm_shuffle_ps(v1,self.data, _ico_shuffle(3, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn set_z<T : Into<FloatVector>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_shuffle_ps(self.data, value.into().data, _ico_shuffle(3, 2, 1, 0));
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

	/// Multiply v1 and v2, and add v3.
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
	pub fn scale<T : Into<FloatVector>>(v1 : Vector3, scalar : T) -> Vector3{	
		unsafe{
			Vector3{data : _mm_mul_ps(v1.data, scalar.into().data)}
		}
	}

	#[inline(always)]
	pub fn div<T : Into<FloatVector>>(v1 : Vector3, scalar : T) -> Vector3{	
		unsafe{
			Vector3{data : _mm_div_ps(v1.data, scalar.into().data)}
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
			return (_mm_movemask_ps(v1.data) & 7) != 0;
		}
	}

	#[inline(always)]
	pub fn equals(v1 : Vector3, v2 : Vector3) -> bool{	
		let d = Vector3::component_equal(v1,v2);
		return Vector3::all(d);
	}
	#[inline(always)]
	pub fn is_zero(v1 : Vector3) -> bool{	
		let d = Vector3::component_equal(v1,Vector3::zero());
		return Vector3::all(d);
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
	/// Right handed cross product
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
	pub fn copysign(magnitude : Vector3, sign : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_copysign_ps(magnitude.data, sign.data)}
		}
	}
	#[inline(always)]
	/// Floor function.  Returns signed 0 when applicable.
	pub fn floor(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_floor_ps(v1.data)}
				//_ico_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Ceil function.  Returns signed 0 when applicable.
	pub fn ceil(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_ceil_ps(v1.data)}
				//_ico_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Round to nearest even function. Returns signed 0 when applicable.
	pub fn round(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : 
			_mm_round_ps(v1.data, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC)}
			// _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _mm_cvttps_epi32(_mm_floor_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _mm_cvttps_epi32(_mm_ceil_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_round_ps(v1.data, _MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC)}
				//_ico_truncate_ps(v1.data)}
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
	pub fn sin(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_sin_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn cos(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_cos_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn acos(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_acos_ps(v1.data)}
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
	pub fn dot(v1 : Vector3, v2 : Vector3) -> FloatVector{	
		unsafe{

			let tmp0 = _mm_mul_ps(v1.data, v2.data); //xyzw
		    let tmp1 = _mm_castsi128_ps(_mm_slli_si128 (_mm_castps_si128(tmp0), 4)); //0xyz
		    let tmp2 = _mm_add_ps(tmp0 , tmp1);//x xy, yz, wz
		    let tmp3 = _mm_moveldup_ps(tmp2); // x x yz yz
			FloatVector{data : _mm_add_ps(tmp3, _wzyx(tmp3))}
		}
	}
	#[inline(always)]
	pub fn lerp<T : Into<FloatVector>>(v1 : Vector3, v2 : Vector3, t : T) -> Vector3{	
		unsafe{
			let t_val = t.into().data;
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector3{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}



	#[inline(always)]
	pub fn normalize(v1 : Vector3) -> Vector3{	
		let length = FloatVector::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::div(v1, length);
		
		unsafe{
			// This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
			let result_length_sqr = Vector3::from(Vector3::dot(norm,norm));
			let mask_less = Vector3::component_less( result_length_sqr, Vector3{data :_ico_two_ps()});
			return Vector3::and(norm, mask_less);
		}
	}

	#[inline(always)]
	pub fn normalize_length(v1 : Vector3) -> (Vector3, FloatVector){	
		let length = FloatVector::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::div(v1, length);
		
		unsafe{

			// This catches infinity, NAN.  Zero vectors are possible - but that is fine - we failed
			let result_length_sqr = Vector3::from(Vector3::dot(norm,norm));
			let mask_less = Vector3::component_less( result_length_sqr, Vector3{data :_ico_two_ps()});
			return (Vector3::and(norm, mask_less), length);
		}
	}

	

	#[inline(always)]
	pub fn sqr_magnitude(self) -> FloatVector {
		return Vector3::dot(self, self);	
	}

	#[inline(always)]
	pub fn magnitude(self) -> FloatVector {
		return FloatVector::sqrt(Vector3::dot(self, self));	
	}

	#[inline(always)] pub fn xxxx(self) -> Vector4 { unsafe{return Vector4{data:_xxxx(self.data)};}}
	#[inline(always)] pub fn yyyy(self) -> Vector4 { unsafe{return Vector4{data:_yyyy(self.data)};}}
	#[inline(always)] pub fn zzzz(self) -> Vector4 { unsafe{return Vector4{data:_zzzz(self.data)};}}

	#[inline(always)] pub fn xxx(self) -> Vector3 { unsafe{return Vector3{data:_xxxw(self.data)};}}
	#[inline(always)] pub fn xxy(self) -> Vector3 { unsafe{return Vector3{data:_xxyw(self.data)};}}
	#[inline(always)] pub fn xxz(self) -> Vector3 { unsafe{return Vector3{data:_xxzw(self.data)};}}
	#[inline(always)] pub fn xyx(self) -> Vector3 { unsafe{return Vector3{data:_xyxw(self.data)};}}
	#[inline(always)] pub fn xyy(self) -> Vector3 { unsafe{return Vector3{data:_xyyw(self.data)};}}
	#[inline(always)] pub fn xyz(self) -> Vector3 { unsafe{return Vector3{data:_xyzw(self.data)};}}
	#[inline(always)] pub fn xzx(self) -> Vector3 { unsafe{return Vector3{data:_xzxw(self.data)};}}
	#[inline(always)] pub fn xzy(self) -> Vector3 { unsafe{return Vector3{data:_xzyw(self.data)};}}
	#[inline(always)] pub fn xzz(self) -> Vector3 { unsafe{return Vector3{data:_xzzw(self.data)};}}

	#[inline(always)] pub fn yxx(self) -> Vector3 { unsafe{return Vector3{data:_yxxw(self.data)};}}
	#[inline(always)] pub fn yxy(self) -> Vector3 { unsafe{return Vector3{data:_yxyw(self.data)};}}
	#[inline(always)] pub fn yxz(self) -> Vector3 { unsafe{return Vector3{data:_yxzw(self.data)};}}
	#[inline(always)] pub fn yyx(self) -> Vector3 { unsafe{return Vector3{data:_yyxw(self.data)};}}
	#[inline(always)] pub fn yyy(self) -> Vector3 { unsafe{return Vector3{data:_yyyw(self.data)};}}
	#[inline(always)] pub fn yyz(self) -> Vector3 { unsafe{return Vector3{data:_yyzw(self.data)};}}
	#[inline(always)] pub fn yzx(self) -> Vector3 { unsafe{return Vector3{data:_yzxw(self.data)};}}
	#[inline(always)] pub fn yzy(self) -> Vector3 { unsafe{return Vector3{data:_yzyw(self.data)};}}
	#[inline(always)] pub fn yzz(self) -> Vector3 { unsafe{return Vector3{data:_yzzw(self.data)};}}

	#[inline(always)] pub fn zxx(self) -> Vector3 { unsafe{return Vector3{data:_zxxw(self.data)};}}
	#[inline(always)] pub fn zxy(self) -> Vector3 { unsafe{return Vector3{data:_zxyw(self.data)};}}
	#[inline(always)] pub fn zxz(self) -> Vector3 { unsafe{return Vector3{data:_zxzw(self.data)};}}
	#[inline(always)] pub fn zyx(self) -> Vector3 { unsafe{return Vector3{data:_zyxw(self.data)};}}
	#[inline(always)] pub fn zyy(self) -> Vector3 { unsafe{return Vector3{data:_zyyw(self.data)};}}
	#[inline(always)] pub fn zyz(self) -> Vector3 { unsafe{return Vector3{data:_zyzw(self.data)};}}
	#[inline(always)] pub fn zzx(self) -> Vector3 { unsafe{return Vector3{data:_zzxw(self.data)};}}
	#[inline(always)] pub fn zzy(self) -> Vector3 { unsafe{return Vector3{data:_zzyw(self.data)};}}
	#[inline(always)] pub fn zzz(self) -> Vector3 { unsafe{return Vector3{data:_zzzw(self.data)};}}
}


impl From<FloatVector> for Vector3 {
	#[inline(always)]
    fn from(v : FloatVector) -> Vector3 {
        return Vector3 { data :v.data};	
    }
}

impl From<Vector2> for Vector3 {
	#[inline(always)]
    fn from(v : Vector2) -> Vector3 {
    	unsafe{
        return Vector3 { data : _mm_movelh_ps(v.data, _mm_setzero_ps() ) };
    	}
    }
}
impl From<Vector4> for Vector3 {
	#[inline(always)]
    fn from(v : Vector4) -> Vector3 {
        Vector3 { data : v.data }
    }
}
impl From<Vector3Int> for Vector3 {
	#[inline(always)]
    fn from(v : Vector3Int) -> Vector3 {
    	unsafe{
        	return Vector3 { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl core::ops::Add for Vector3{
	type Output = Vector3;
	#[inline(always)]
	fn add(self, _rhs: Vector3) -> Vector3{
		Vector3::add(self, _rhs)
	}
}
impl core::ops::AddAssign for Vector3 {
	#[inline(always)]
    fn add_assign(&mut self, other: Vector3) {
        *self = Vector3::add(*self, other)
    }
}

impl core::ops::Sub for Vector3{
	type Output = Vector3;
	#[inline(always)]
	fn sub(self, _rhs: Vector3) -> Vector3{
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
		unsafe{
			return Vector3{data:_mm_xor_ps(_ico_signbit_ps(),self.data)};
		}
	}
}
impl<T : Into<FloatVector>> core::ops::Mul<T> for Vector3{
	type Output = Vector3;
	#[inline(always)]
	fn mul(self, _rhs: T) -> Vector3{
		Vector3::scale(self, _rhs.into())
	}
}
impl<T : Into<FloatVector>> core::ops::MulAssign<T> for Vector3{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: T){
		*self = Vector3::scale(*self, _rhs.into())
	}
}
impl core::ops::Mul<Vector3> for FloatVector{
	type Output = Vector3;
	#[inline(always)]
	fn mul(self : FloatVector, _rhs: Vector3) -> Vector3{
		Vector3::scale(_rhs, self)
	}
}

impl<T : Into<FloatVector>> core::ops::Div<T> for Vector3{
	type Output = Vector3;
	#[inline(always)]
	fn div(self, _rhs: T) -> Vector3{
		Vector3::div(self, _rhs.into())
	}
}
impl core::ops::Div<Vector3> for FloatVector{
	type Output = Vector3;
	#[inline(always)]
	fn div(self : FloatVector, _rhs: Vector3) -> Vector3{
		return Vector3::component_div(Vector3::from(self), _rhs);
	}
}
impl<T : Into<FloatVector>> core::ops::DivAssign<T> for Vector3{
	#[inline(always)]
	fn div_assign(&mut self, _rhs: T){
		*self = Vector3::div(*self, _rhs.into())
	}
}	
impl PartialEq for Vector3 {
	#[inline(always)]
    fn eq(&self, other: &Vector3) -> bool {
    	return Vector3::equals(*self, *other);
    }
}

#[cfg(test)]
mod test;