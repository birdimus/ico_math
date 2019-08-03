use core::arch::x86_64::*;
use core::hash::Hasher;
use core::hash::Hash;
use crate::Vector3;
use crate::RawIntVector;
use crate::IntVector;
use crate::Vector2Int;
use crate::Vector3Int;
use crate::Vector4Int;
use crate::sse_extensions::*;



impl Vector3Int{
	/// Returns a new Vector3
	#[inline(always)]
	pub fn new(x : i32, y : i32, z : i32) -> Vector3Int{
		unsafe{
			Vector3Int{data : _mm_set_epi32(0, z, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector3Int {
		unsafe{
			Vector3Int { data : _mm_setzero_si128() }
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
	pub fn z(self) -> IntVector {
		return IntVector{data:self.zzzz().data};	
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
	pub fn set_z<T : Into<i32>>(&self, value : T) {
		unsafe{
			_mm_insert_epi32(self.data, value.into(), 2);
		}	
	}


	#[inline(always)]
	pub fn max(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_max_epi32(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_min_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn abs(v1 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_abs_epi32(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_sign_epi32(v1.data, v2.data)}
		}
	}

	/*
	// TODO: Requires the shifts to be constant - waiting on const generics, or const arguments
	// https://github.com/rust-lang/rfcs/pull/2000
	// could also switch on the shift var. 0-32
	#[inline(always)]
	pub fn shift_right(v1 : Vector3Int, shift : i32) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_srli_epi32(v1.data, shift)}
		}
	}

	#[inline(always)]
	pub fn shift_left(v1 : Vector3Int, shift : i32) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_slli_epi32(v1.data, shift)}
		}
	}
	*/
	#[inline(always)]
	pub fn add(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_add_epi32(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_sub_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_mul(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_mullo_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn scale<T : Into<IntVector>>(v1 : Vector3Int, scalar : T) -> Vector3Int{	
		unsafe{
			return Vector3Int{data : _mm_mullo_epi32(v1.data, scalar.into().data)};
		}
	}

	#[inline(always)]
	pub fn and(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_and_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn or(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_or_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn andnot(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_andnot_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn xor(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_xor_si128(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_cmpeq_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_cmpgt_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector3Int, v2 : Vector3Int) -> Vector3Int{	
		unsafe{
			Vector3Int{data : _mm_cmplt_epi32(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn all(v1 : Vector3Int) -> bool{	
		unsafe{
			return (_mm_movemask_epi8 (v1.data) & 4095 ) == 4095;
		}
	}
	#[inline(always)]
	pub fn any(v1 : Vector3Int) -> bool{	
		unsafe{
			return (_mm_movemask_epi8 (v1.data) & 4095 ) != 0;
		}
	}
	#[inline(always)]
	pub fn equals(v1 : Vector3Int, v2 : Vector3Int) -> bool{	
		unsafe{
			let d = _mm_cmpeq_epi32(v1.data, v2.data);
			return (_mm_movemask_epi8(d) & 4095) == 4095;
		}
	}


	#[inline(always)] pub fn xxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxx_i(self.data)};}}
	#[inline(always)] pub fn yyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyy_i(self.data)};}}
	#[inline(always)] pub fn zzzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzzz_i(self.data)};}}

	#[inline(always)] pub fn xxx(self) -> Vector3Int { unsafe{return Vector3Int{data:_xxxw_i(self.data)};}}
	#[inline(always)] pub fn xxy(self) -> Vector3Int { unsafe{return Vector3Int{data:_xxyw_i(self.data)};}}
	#[inline(always)] pub fn xxz(self) -> Vector3Int { unsafe{return Vector3Int{data:_xxzw_i(self.data)};}}
	#[inline(always)] pub fn xyx(self) -> Vector3Int { unsafe{return Vector3Int{data:_xyxw_i(self.data)};}}
	#[inline(always)] pub fn xyy(self) -> Vector3Int { unsafe{return Vector3Int{data:_xyyw_i(self.data)};}}
	#[inline(always)] pub fn xyz(self) -> Vector3Int { unsafe{return Vector3Int{data:_xyzw_i(self.data)};}}
	#[inline(always)] pub fn xzx(self) -> Vector3Int { unsafe{return Vector3Int{data:_xzxw_i(self.data)};}}
	#[inline(always)] pub fn xzy(self) -> Vector3Int { unsafe{return Vector3Int{data:_xzyw_i(self.data)};}}
	#[inline(always)] pub fn xzz(self) -> Vector3Int { unsafe{return Vector3Int{data:_xzzw_i(self.data)};}}

	#[inline(always)] pub fn yxx(self) -> Vector3Int { unsafe{return Vector3Int{data:_yxxw_i(self.data)};}}
	#[inline(always)] pub fn yxy(self) -> Vector3Int { unsafe{return Vector3Int{data:_yxyw_i(self.data)};}}
	#[inline(always)] pub fn yxz(self) -> Vector3Int { unsafe{return Vector3Int{data:_yxzw_i(self.data)};}}
	#[inline(always)] pub fn yyx(self) -> Vector3Int { unsafe{return Vector3Int{data:_yyxw_i(self.data)};}}
	#[inline(always)] pub fn yyy(self) -> Vector3Int { unsafe{return Vector3Int{data:_yyyw_i(self.data)};}}
	#[inline(always)] pub fn yyz(self) -> Vector3Int { unsafe{return Vector3Int{data:_yyzw_i(self.data)};}}
	#[inline(always)] pub fn yzx(self) -> Vector3Int { unsafe{return Vector3Int{data:_yzxw_i(self.data)};}}
	#[inline(always)] pub fn yzy(self) -> Vector3Int { unsafe{return Vector3Int{data:_yzyw_i(self.data)};}}
	#[inline(always)] pub fn yzz(self) -> Vector3Int { unsafe{return Vector3Int{data:_yzzw_i(self.data)};}}

	#[inline(always)] pub fn zxx(self) -> Vector3Int { unsafe{return Vector3Int{data:_zxxw_i(self.data)};}}
	#[inline(always)] pub fn zxy(self) -> Vector3Int { unsafe{return Vector3Int{data:_zxyw_i(self.data)};}}
	#[inline(always)] pub fn zxz(self) -> Vector3Int { unsafe{return Vector3Int{data:_zxzw_i(self.data)};}}
	#[inline(always)] pub fn zyx(self) -> Vector3Int { unsafe{return Vector3Int{data:_zyxw_i(self.data)};}}
	#[inline(always)] pub fn zyy(self) -> Vector3Int { unsafe{return Vector3Int{data:_zyyw_i(self.data)};}}
	#[inline(always)] pub fn zyz(self) -> Vector3Int { unsafe{return Vector3Int{data:_zyzw_i(self.data)};}}
	#[inline(always)] pub fn zzx(self) -> Vector3Int { unsafe{return Vector3Int{data:_zzxw_i(self.data)};}}
	#[inline(always)] pub fn zzy(self) -> Vector3Int { unsafe{return Vector3Int{data:_zzyw_i(self.data)};}}
	#[inline(always)] pub fn zzz(self) -> Vector3Int { unsafe{return Vector3Int{data:_zzzw_i(self.data)};}}
	
}

impl From<IntVector> for Vector3Int {
	#[inline(always)]
    fn from(val : IntVector) -> Vector3Int {
       unsafe{
			return Vector3Int{data : val.data};
		}
    }
}

impl From<Vector3> for Vector3Int {
	#[inline(always)]
    fn from(v : Vector3) -> Vector3Int {
    	unsafe{
        	return Vector3Int { data : _mm_cvttps_epi32(v.data) };
        }
    }
}
impl From<Vector2Int> for Vector3Int {
	#[inline(always)]
    fn from(v : Vector2Int) -> Vector3Int {
        return Vector3Int { data : v.data };
    }
}
impl From<Vector4Int> for Vector3Int {
	#[inline(always)]
    fn from(v : Vector4Int) -> Vector3Int {
        return Vector3Int { data : v.data };
    }
}


impl core::ops::Add for Vector3Int{
	type Output = Vector3Int;
	#[inline(always)]
	fn add(self, _rhs: Vector3Int) -> Vector3Int{
		Vector3Int::add(self, _rhs)
	}
}
impl core::ops::AddAssign for Vector3Int {
	#[inline(always)]
    fn add_assign(&mut self, other: Vector3Int) {
        *self = Vector3Int::add(*self, other)
    }
}
impl core::ops::Sub for Vector3Int{
	type Output = Vector3Int;
	#[inline(always)]
	fn sub(self, _rhs: Vector3Int) -> Vector3Int{
		Vector3Int::sub(self, _rhs)
	}
}
impl core::ops::SubAssign for Vector3Int {
	#[inline(always)]
    fn sub_assign(&mut self, other: Vector3Int) {
        *self = Vector3Int::sub(*self, other)
    }
}
impl core::ops::Neg for Vector3Int {
	type Output = Vector3Int;
	#[inline(always)]
	fn neg(self) -> Self::Output {
		unsafe{
			return Vector3Int{data:_mm_sub_epi32(_mm_set1_epi32(0), self.data)};
		}
	}
}

impl<T : Into<IntVector>> core::ops::Mul<T> for Vector3Int{
	type Output = Vector3Int;
	#[inline(always)]
	fn mul(self, _rhs: T) -> Vector3Int{
		Vector3Int::scale(self, _rhs.into())
	}
}
impl<T : Into<IntVector>> core::ops::MulAssign<T> for Vector3Int{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: T){
		*self = Vector3Int::scale(*self, _rhs.into())
	}
}
impl core::ops::Mul<Vector3Int> for IntVector{
	type Output = Vector3Int;
	#[inline(always)]
	fn mul(self : IntVector, _rhs: Vector3Int) -> Vector3Int{
		Vector3Int::scale(_rhs, self)
	}
}

impl PartialEq for Vector3Int {
	#[inline(always)]
    fn eq(&self, other: &Vector3Int) -> bool {
    	return Vector3Int::equals(*self, *other);
    }
}
impl Hash for Vector3Int {
    fn hash<H: Hasher>(&self, state: &mut H) {
    	let mut dst = RawIntVector{data:[0;4]};
    	unsafe{
    		let x : *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
    		_mm_store_si128(x, self.data);
    	}
    	dst.data[3] = 0;
        dst.data.hash(state);
    }
}