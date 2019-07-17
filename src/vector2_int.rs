use crate::Vector2;
use crate::Vector2Int;
use crate::Vector3Int;
use crate::Vector4Int;
use std::arch::x86_64::*;



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
	pub fn x(&self) -> i32 {
		unsafe{
			return _mm_cvtsi128_si32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(&self) -> i32 {
		unsafe{
			return  _mm_extract_epi32(self.data,1);
		}	
	}


	#[inline(always)]
	pub fn set_x(&self, value : i32) {
		unsafe{
			_mm_insert_epi32(self.data, value, 0);
		}	
	}

	#[inline(always)]
	pub fn set_y(&self, value : i32) {
		unsafe{
			_mm_insert_epi32(self.data, value, 1);
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
			Vector2Int{data : _mm_add_epi32(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_sub_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_mul(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_mullo_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn scale(v1 : Vector2Int, scalar : i32) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_mullo_epi32(v1.data, _mm_set1_epi32(scalar))}
		}
	}
	#[inline(always)]
	pub fn and(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_and_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn or(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_or_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn andnot(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_andnot_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn xor(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_xor_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_cmpeq_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_cmpgt_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector2Int, v2 : Vector2Int) -> Vector2Int{	
		unsafe{
			Vector2Int{data : _mm_cmplt_epi32(v1.data, v2.data)}
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
	
}

impl From<i32> for Vector2Int {
    fn from(val : i32) -> Vector2Int {
       unsafe{
			return Vector2Int{data : _mm_set1_epi32 (val)};
		}
    }
}

impl From<Vector2> for Vector2Int {
    fn from(v : Vector2) -> Vector2Int {
    	unsafe{
        	return Vector2Int { data : _mm_cvttps_epi32(v.data) };
        }
    }
}
impl From<Vector3Int> for Vector2Int {
    fn from(v : Vector3Int) -> Vector2Int {
        return Vector2Int { data : v.data };
    }
}
impl From<Vector4Int> for Vector2Int {
    fn from(v : Vector4Int) -> Vector2Int {
        return Vector2Int { data : v.data };
    }
}


impl std::ops::Add for Vector2Int{
	type Output = Vector2Int;
	#[inline]
	fn add(self, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::add(self, _rhs)
	}
}

impl std::ops::Sub for Vector2Int{
	type Output = Vector2Int;
	#[inline]
	fn sub(self, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::sub(self, _rhs)
	}
}
impl std::ops::Mul<i32> for Vector2Int{
	type Output = Vector2Int;
	#[inline]
	fn mul(self, _rhs: i32) -> Vector2Int{
		Vector2Int::scale(self, _rhs)
	}
}

impl std::ops::Mul<Vector2Int> for i32{
	type Output = Vector2Int;
	#[inline]
	fn mul(self, _rhs: Vector2Int) -> Vector2Int{
		Vector2Int::scale(_rhs, self)
	}
}
impl PartialEq for Vector2Int {
    fn eq(&self, other: &Vector2Int) -> bool {
    	return Vector2Int::equals(*self, *other);
    }
}