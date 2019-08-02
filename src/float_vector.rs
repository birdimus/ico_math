use crate::FloatVector;


use core::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_abs_ps;
use crate::_ico_truncate_ps;
use crate::_ico_copysign_ps;
use crate::_ico_two_ps;
use crate::_ico_signbit_ps;
use crate::*;
impl FloatVector{
	/// Returns a new FloatVector
	#[inline(always)]
	pub fn new(value : f32) -> FloatVector{
		unsafe{
			FloatVector{data : _mm_set1_ps(value)}
		}
	}
	#[inline(always)]
	pub fn zero() -> FloatVector {
		unsafe{
			FloatVector { data : _mm_setzero_ps() }
		}
	}

	#[inline(always)]
	pub fn value(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}





	#[inline(always)]
	pub fn add(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_add_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_sub_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn mul(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_mul_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn div(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_div_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn fmadd(v1 : FloatVector, v2 : FloatVector, v3 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_fmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fmsub(v1 : FloatVector, v2 : FloatVector, v3 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_fmsub_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmadd(v1 : FloatVector, v2 : FloatVector, v3 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_fnmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmsub(v1 : FloatVector, v2 : FloatVector, v3 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_fnmsub_ps(v1.data, v2.data, v3.data)}
		}
	}


	#[inline(always)]
	pub fn and(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_and_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn or(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_or_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn andnot(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_andnot_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn xor(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_xor_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn equals(v1 : FloatVector, v2 : FloatVector) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.data, v2.data);
			return (_mm_movemask_ps(d) ) == 15;
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmpeq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_not_equal(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmpneq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_greater_equal(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmpge_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmpgt_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less_equal(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmple_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_cmplt_ps(v1.data, v2.data)}
		}
	}


	#[inline(always)]
	pub fn abs(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  _ico_abs_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  _ico_copysign_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	/// Floor function.  Returns signed 0 when applicable.
	pub fn floor(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  
				_mm_floor_ps(v1.data)}
				//_ico_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Ceil function.  Returns signed 0 when applicable.
	pub fn ceil(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  
				_mm_ceil_ps(v1.data)}
				//_ico_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Round to nearest even function. Returns signed 0 when applicable.
	pub fn round(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : 
			_mm_round_ps(v1.data, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC)}
			// _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : FloatVector) -> IntVector{	
		unsafe{
			IntVector{data :  _mm_cvttps_epi32(_mm_floor_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : FloatVector) -> IntVector{	
		unsafe{
			IntVector{data :  _mm_cvttps_epi32(_mm_ceil_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  
				_mm_round_ps(v1.data, _MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC)}
				//_ico_truncate_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn frac(v1 : FloatVector) -> FloatVector{	
		return FloatVector::sub(v1, FloatVector::floor(v1));
	}

	#[inline(always)]
	pub fn sqrt(v1 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn max(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : FloatVector, v2 : FloatVector) -> FloatVector{	
		unsafe{
			FloatVector{data : _mm_min_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn lerp(v1 : FloatVector, v2 : FloatVector, t : FloatVector) -> FloatVector{	
		unsafe{
			let tmp = _mm_fnmadd_ps(v1.data, t.data, v1.data); //a - (a*t)
			FloatVector{data : _mm_fmadd_ps(v2.data, t.data, tmp)} //b * t + a
		}
	}


}

impl From<FloatVector> for f32 {
    fn from(v : FloatVector) -> f32 {
    	return v.value();
    }
}
impl From<f32> for FloatVector {
    fn from(v : f32) -> FloatVector {
    	return FloatVector::new(v);
    }
}
impl From<IntVector> for FloatVector {
    fn from(v : IntVector) -> FloatVector {
    	unsafe{
        	return FloatVector { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl core::ops::Add for FloatVector{
	type Output = FloatVector;
	#[inline]
	fn add(self, _rhs: FloatVector) -> FloatVector{
		FloatVector::add(self, _rhs)
	}
}
impl core::ops::AddAssign for FloatVector {
	#[inline(always)]
    fn add_assign(&mut self, other: FloatVector) {
        *self = FloatVector::add(*self, other)
    }
}
impl core::ops::Sub for FloatVector{
	type Output = FloatVector;
	#[inline]
	fn sub(self, _rhs: FloatVector) -> FloatVector{
		FloatVector::sub(self, _rhs)
	}
}
impl core::ops::SubAssign for FloatVector {
	#[inline(always)]
    fn sub_assign(&mut self, other: FloatVector) {
        *self = FloatVector::sub(*self, other)
    }
}
impl core::ops::Neg for FloatVector {
	type Output = FloatVector;
	#[inline(always)]
	fn neg(self) -> Self::Output {
		unsafe{
			return FloatVector{data:_mm_xor_ps(_ico_signbit_ps(),self.data)};
		}
	}
}
impl core::ops::Mul<FloatVector> for FloatVector{
	type Output = FloatVector;
	#[inline]
	fn mul(self, _rhs: FloatVector) -> FloatVector{
		FloatVector::mul(self, _rhs)
	}
}


impl core::ops::MulAssign<FloatVector> for FloatVector{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: FloatVector){
		*self = FloatVector::mul(*self, _rhs)
	}
}
impl core::ops::Div<FloatVector> for FloatVector{
	type Output = FloatVector;
	#[inline]
	fn div(self, _rhs: FloatVector) -> FloatVector{
		FloatVector::div(self, _rhs)
	}
}
impl core::ops::DivAssign<FloatVector> for FloatVector{
	#[inline(always)]
	fn div_assign(&mut self, _rhs: FloatVector){
		*self = FloatVector::div(*self, _rhs)
	}
}
impl PartialEq for FloatVector {
    fn eq(&self, other: &FloatVector) -> bool {
    	return FloatVector::equals(*self, *other);
    }
}
impl PartialEq<f32> for FloatVector {
	fn eq(&self, other: &f32) -> bool {
    	return self.value() == *other;
    }
    
}
impl PartialEq<FloatVector> for f32 {
    fn eq(&self, other: &FloatVector) -> bool {
    	return *self == other.value();
    }
}
// #[cfg(test)]
// mod test;