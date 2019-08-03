use core::arch::x86_64::*;
use core::hash::Hasher;
use core::hash::Hash;
use crate::RawIntVector;
use crate::IntVector;
use crate::Vector2;
use crate::Vector2Int;
use crate::Vector3Int;
use crate::Vector4Int;
use crate::sse_extensions::*;



impl Vector2Int{
	/// Returns a new Vector2
	#[inline(always)]
	pub fn new(x : i32, y : i32) -> Vector2Int{
		unsafe{
			Vector2Int{data : _mm_set_epi32(0, 0, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector2Int {
		unsafe{
			Vector2Int { data : _mm_setzero_si128() }
		}
	}
	#[inline(always)]
	pub fn x(self) -> IntVector {
		return IntVector{data:self.xxxx().data};	
	}

	#[inline(always)]
	pub fn y(self) -> IntVector {
		return IntVector{data:self.yyyy().data};	
	}

	#[inline(always)]
	pub fn set_x<T : Into<i32>>(&self, value : T) {
		unsafe{
			_mm_insert_epi32(self.data, value.into(), 0);
		}	
	}

	#[inline(always)]
	pub fn set_y<T : Into<i32>>(&self, value : T) {
		unsafe{
			_mm_insert_epi32(self.data, value.into(), 1);
		}	
	}

	#[inline(always)]
	pub fn max(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_max_epi32(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_min_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn abs(v1 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_abs_epi32(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_sign_epi32(v1.data, v2.data)}
		}
	}

	/*
	// TODO: Requires the shifts to be constant - waiting on const generics, or const arguments
	// https://github.com/rust-lang/rfcs/pull/2000
	// could also switch on the shift var. 0-32
	#[inline(always)]
	pub fn shift_right(v1 : Vector2Int, shift : i32) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_srli_epi32(v1.data, shift)}
		}
	}

	#[inline(always)]
	pub fn shift_left(v1 : Vector2Int, shift : i32) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_slli_epi32(v1.data, shift)}
		}
	}
	*/
	#[inline(always)]
	pub fn add(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_add_epi32(v1.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_sub_epi32(v1.data, v2.data)};
		}
	}
	#[inline(always)]
	pub fn component_mul(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_mullo_epi32(v1.data, v2.data)};
		}
	}
	#[inline(always)]
	pub fn scale<T : Into<IntVector>>(v1 : Vector2Int, scalar : T) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_mullo_epi32(v1.data, scalar.into().data)};
		}
	}
	#[inline(always)]
	pub fn and(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_and_si128(v1.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn or(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_or_si128(v1.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn andnot(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_andnot_si128(v1.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn xor(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_xor_si128(v1.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_cmpeq_epi32(v1.data, v2.data)};
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_cmpgt_epi32(v1.data, v2.data)};
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			return Vector2Int{data : _mm_cmplt_epi32(v1.data, v2.data)};
		}
	}
	#[inline(always)]
	pub fn all(v1 : Vector2Int) -> bool{	
		unsafe{
			return (_mm_movemask_epi8 (v1.data) & 255 ) == 255;
		}
	}
	#[inline(always)]
	pub fn any(v1 : Vector2Int) -> bool{	
		unsafe{
			return (_mm_movemask_epi8 (v1.data) & 255 ) != 0;
		}
	}
	#[inline(always)]
	pub fn equals(v1 : Vector2Int, v2 : Vector2Int) -> bool{	
		unsafe{
			let d = _mm_cmpeq_epi32(v1.data, v2.data);
			return (_mm_movemask_epi8(d) & 255) == 255;
		}
	}


	#[inline(always)] pub fn xxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxx_i(self.data)};}}
	#[inline(always)] pub fn yyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyy_i(self.data)};}}

	#[inline(always)] pub fn xx(self) -> Vector2Int { unsafe{return Vector2Int{data:_xxzw_i(self.data)};}}
	#[inline(always)] pub fn xy(self) -> Vector2Int { unsafe{return Vector2Int{data:_xyzw_i(self.data)};}}
	#[inline(always)] pub fn yx(self) -> Vector2Int { unsafe{return Vector2Int{data:_yxzw_i(self.data)};}}
	#[inline(always)] pub fn yy(self) -> Vector2Int { unsafe{return Vector2Int{data:_yyzw_i(self.data)};}}
	
}

impl From<IntVector> for Vector2Int {
	#[inline(always)]
    fn from(val : IntVector) -> Vector2Int {
       unsafe{
			return Vector2Int{data : val.data};
		}
    }
}

impl From<Vector2> for Vector2Int {
	#[inline(always)]
    fn from(v : Vector2) -> Vector2Int {
    	unsafe{
        	return Vector2Int { data : _mm_cvttps_epi32(v.data) };
        }
    }
}
impl From<Vector3Int> for Vector2Int {
	#[inline(always)]
    fn from(v : Vector3Int) -> Vector2Int {
        return Vector2Int { data : v.data };
    }
}
impl From<Vector4Int> for Vector2Int {
	#[inline(always)]
    fn from(v : Vector4Int) -> Vector2Int {
        return Vector2Int { data : v.data };
    }
}


impl core::ops::Add for Vector2Int{
	type Output = Vector2Int;
	#[inline(always)]
	fn add(self, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::add(self, _rhs)
	}
}
impl core::ops::AddAssign for Vector2Int {
	#[inline(always)]
    fn add_assign(&mut self, other: Vector2Int) {
        *self = Vector2Int::add(*self, other)
    }
}
impl core::ops::Sub for Vector2Int{
	type Output = Vector2Int;
	#[inline(always)]
	fn sub(self, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::sub(self, _rhs)
	}
}
impl core::ops::SubAssign for Vector2Int {
	#[inline(always)]
    fn sub_assign(&mut self, other: Vector2Int) {
        *self = Vector2Int::sub(*self, other)
    }
}
impl core::ops::Neg for Vector2Int {
	type Output = Vector2Int;
	#[inline(always)]
	fn neg(self) -> Self::Output {
		unsafe{
			return Vector2Int{data:_mm_sub_epi32(_mm_set1_epi32(0), self.data)};
		}
	}
}

impl<T : Into<IntVector>> core::ops::Mul<T> for Vector2Int{
	type Output = Vector2Int;
	#[inline(always)]
	fn mul(self, _rhs: T) -> Vector2Int{
		Vector2Int::scale(self, _rhs.into())
	}
}
impl<T : Into<IntVector>> core::ops::MulAssign<T> for Vector2Int{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: T){
		*self = Vector2Int::scale(*self, _rhs.into())
	}
}
impl core::ops::Mul<Vector2Int> for IntVector{
	type Output = Vector2Int;
	#[inline(always)]
	fn mul(self : IntVector, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::scale(_rhs, self)
	}
}
impl PartialEq for Vector2Int {
	#[inline(always)]
    fn eq(&self, other: &Vector2Int) -> bool {
    	return Vector2Int::equals(*self, *other);
    }
}
impl Hash for Vector2Int {
    fn hash<H: Hasher>(&self, state: &mut H) {
    	let mut dst = RawIntVector{data:[0;4]};
    	unsafe{
    		let x : *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
    		_mm_store_si128(x, self.data);
    	}
    	dst.data[2] = 0;
    	dst.data[3] = 0;
        dst.data.hash(state);
    }
}
