use core::arch::x86_64::*;
use core::hash::Hasher;
use core::hash::Hash;
use crate::structure::SIMDVector4;
use crate::vector4::Vector4;
use crate::vector4_bool::Vector4Bool;
use crate::raw::RawVector_i32;
use crate::int_vector::IntVector;
use crate::vector2_int::Vector2Int;
use crate::vector3_int::Vector3Int;
use crate::sse_extensions::*;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4Int{
	pub data : __m128i,
}


impl Vector4Int{
	/// Returns a new Vector4
	#[inline(always)]
	pub fn new(x : i32, y : i32, z : i32, w : i32) -> Vector4Int{
		unsafe{
			Vector4Int{data : _mm_set_epi32(w, z, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector4Int {
		unsafe{
			Vector4Int { data : _mm_setzero_si128() }
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
	pub fn w(self) -> IntVector {
		return IntVector{data:self.wwww().data};	
	}

	#[inline(always)]
	pub fn set_x<T : Into<i32>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_insert_epi32(self.data, value.into(), 0);
		}	
	}

	#[inline(always)]
	pub fn set_y<T : Into<i32>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_insert_epi32(self.data, value.into(), 1);
		}	
	}

	#[inline(always)]
	pub fn set_z<T : Into<i32>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_insert_epi32(self.data, value.into(), 2);
		}	
	}
	#[inline(always)]
	pub fn set_w<T : Into<i32>>(&mut self, value : T) {
		unsafe{
			self.data = _mm_insert_epi32(self.data, value.into(), 3);
		}	
	}
	
	#[inline(always)]
	pub fn max(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_max_epi32(self.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_min_epi32(self.data, v2.data)}
		}
	}
	/// Choose component wise between A and B based on the mask.  False = A, True = B.
	#[inline(always)]
	pub fn select(self, v2 : Vector4Int, mask : Vector4Bool) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _ico_select_si128(self.data, v2.data, mask.data)}
		}
	}

	#[inline(always)]
	pub fn add(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_add_epi32(self.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_sub_epi32(self.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn mul(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_mullo_epi32(self.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn abs(self) -> Vector4Int{	
		unsafe{
			return Vector4Int{data : _mm_abs_epi32(self.data)};
		}
	}
	#[inline(always)]
	pub fn copysign(self, v2 : Vector4Int) -> Vector4Int{	
		unsafe{
			return Vector4Int{data : _mm_sign_epi32(self.data, v2.data)};
		}
	}

	#[inline(always)]
	pub fn and<T : SIMDVector4>(self, v2 : T) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_and_si128(self.data, v2.data_i())}
		}
	}

	#[inline(always)]
	pub fn or<T : SIMDVector4>(self, v2 : T) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_or_si128(self.data, v2.data_i())}
		}
	}

	#[inline(always)]
	pub fn andnot<T : SIMDVector4>(self, v2 : T) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_andnot_si128(self.data, v2.data_i())}
		}
	}

	#[inline(always)]
	pub fn xor<T : SIMDVector4>(self, v2 : T) -> Vector4Int{	
		unsafe{
			Vector4Int{data : _mm_xor_si128(self.data, v2.data_i())}
		}
	}

	#[inline(always)]
	pub fn equal(self, v2 : Vector4Int) -> Vector4Bool{	
		unsafe{
			Vector4Bool{data : _mm_cmpeq_epi32(self.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn not_equal(self, v2 : Vector4Int) -> Vector4Bool{	
		return Vector4Bool::xor(self.equal(self), self.equal(v2));
	}

	#[inline(always)]
	pub fn greater_equal(self, v2 : Vector4Int) -> Vector4Bool{	
		return Vector4Bool::or(self.equal(v2), self.greater(v2));
	}
	#[inline(always)]
	pub fn greater(self, v2 : Vector4Int) -> Vector4Bool{	
		unsafe{
			Vector4Bool{data : _mm_cmpgt_epi32(self.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn less_equal(self, v2 : Vector4Int) -> Vector4Bool{	
		return Vector4Bool::or(self.equal(v2), self.less(v2));
	}

	#[inline(always)]
	pub fn less(self, v2 : Vector4Int) -> Vector4Bool{	
		unsafe{
			Vector4Bool{data : _mm_cmplt_epi32(self.data, v2.data)}
		}
	}


	
	
	#[inline(always)] pub fn xxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxx_i(self.data)};}}
	#[inline(always)] pub fn xxxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxy_i(self.data)};}}
	#[inline(always)] pub fn xxxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxz_i(self.data)};}}
	#[inline(always)] pub fn xxxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxxw_i(self.data)};}}
	#[inline(always)] pub fn xxyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxyx_i(self.data)};}}
	#[inline(always)] pub fn xxyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxyy_i(self.data)};}}
	#[inline(always)] pub fn xxyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxyz_i(self.data)};}}
	#[inline(always)] pub fn xxyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxyw_i(self.data)};}}
	#[inline(always)] pub fn xxzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxzx_i(self.data)};}}
	#[inline(always)] pub fn xxzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxzy_i(self.data)};}}
	#[inline(always)] pub fn xxzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxzz_i(self.data)};}}
	#[inline(always)] pub fn xxzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxzw_i(self.data)};}}
	#[inline(always)] pub fn xxwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxwx_i(self.data)};}}
	#[inline(always)] pub fn xxwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxwy_i(self.data)};}}
	#[inline(always)] pub fn xxwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxwz_i(self.data)};}}
	#[inline(always)] pub fn xxww(self) -> Vector4Int { unsafe{return Vector4Int{data:_xxww_i(self.data)};}}
	#[inline(always)] pub fn xyxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyxx_i(self.data)};}}
	#[inline(always)] pub fn xyxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyxy_i(self.data)};}}
	#[inline(always)] pub fn xyxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyxz_i(self.data)};}}
	#[inline(always)] pub fn xyxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyxw_i(self.data)};}}
	#[inline(always)] pub fn xyyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyyx_i(self.data)};}}
	#[inline(always)] pub fn xyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyyy_i(self.data)};}}
	#[inline(always)] pub fn xyyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyyz_i(self.data)};}}
	#[inline(always)] pub fn xyyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyyw_i(self.data)};}}
	#[inline(always)] pub fn xyzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyzx_i(self.data)};}}
	#[inline(always)] pub fn xyzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyzy_i(self.data)};}}
	#[inline(always)] pub fn xyzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyzz_i(self.data)};}}
	#[inline(always)] pub fn xyzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyzw_i(self.data)};}}
	#[inline(always)] pub fn xywx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xywx_i(self.data)};}}
	#[inline(always)] pub fn xywy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xywy_i(self.data)};}}
	#[inline(always)] pub fn xywz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xywz_i(self.data)};}}
	#[inline(always)] pub fn xyww(self) -> Vector4Int { unsafe{return Vector4Int{data:_xyww_i(self.data)};}}
	#[inline(always)] pub fn xzxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzxx_i(self.data)};}}
	#[inline(always)] pub fn xzxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzxy_i(self.data)};}}
	#[inline(always)] pub fn xzxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzxz_i(self.data)};}}
	#[inline(always)] pub fn xzxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzxw_i(self.data)};}}
	#[inline(always)] pub fn xzyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzyx_i(self.data)};}}
	#[inline(always)] pub fn xzyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzyy_i(self.data)};}}
	#[inline(always)] pub fn xzyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzyz_i(self.data)};}}
	#[inline(always)] pub fn xzyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzyw_i(self.data)};}}
	#[inline(always)] pub fn xzzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzzx_i(self.data)};}}
	#[inline(always)] pub fn xzzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzzy_i(self.data)};}}
	#[inline(always)] pub fn xzzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzzz_i(self.data)};}}
	#[inline(always)] pub fn xzzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzzw_i(self.data)};}}
	#[inline(always)] pub fn xzwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzwx_i(self.data)};}}
	#[inline(always)] pub fn xzwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzwy_i(self.data)};}}
	#[inline(always)] pub fn xzwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzwz_i(self.data)};}}
	#[inline(always)] pub fn xzww(self) -> Vector4Int { unsafe{return Vector4Int{data:_xzww_i(self.data)};}}
	#[inline(always)] pub fn xwxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwxx_i(self.data)};}}
	#[inline(always)] pub fn xwxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwxy_i(self.data)};}}
	#[inline(always)] pub fn xwxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwxz_i(self.data)};}}
	#[inline(always)] pub fn xwxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwxw_i(self.data)};}}
	#[inline(always)] pub fn xwyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwyx_i(self.data)};}}
	#[inline(always)] pub fn xwyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwyy_i(self.data)};}}
	#[inline(always)] pub fn xwyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwyz_i(self.data)};}}
	#[inline(always)] pub fn xwyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwyw_i(self.data)};}}
	#[inline(always)] pub fn xwzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwzx_i(self.data)};}}
	#[inline(always)] pub fn xwzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwzy_i(self.data)};}}
	#[inline(always)] pub fn xwzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwzz_i(self.data)};}}
	#[inline(always)] pub fn xwzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwzw_i(self.data)};}}
	#[inline(always)] pub fn xwwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwwx_i(self.data)};}}
	#[inline(always)] pub fn xwwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwwy_i(self.data)};}}
	#[inline(always)] pub fn xwwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwwz_i(self.data)};}}
	#[inline(always)] pub fn xwww(self) -> Vector4Int { unsafe{return Vector4Int{data:_xwww_i(self.data)};}}
	#[inline(always)] pub fn yxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxxx_i(self.data)};}}
	#[inline(always)] pub fn yxxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxxy_i(self.data)};}}
	#[inline(always)] pub fn yxxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxxz_i(self.data)};}}
	#[inline(always)] pub fn yxxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxxw_i(self.data)};}}
	#[inline(always)] pub fn yxyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxyx_i(self.data)};}}
	#[inline(always)] pub fn yxyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxyy_i(self.data)};}}
	#[inline(always)] pub fn yxyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxyz_i(self.data)};}}
	#[inline(always)] pub fn yxyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxyw_i(self.data)};}}
	#[inline(always)] pub fn yxzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxzx_i(self.data)};}}
	#[inline(always)] pub fn yxzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxzy_i(self.data)};}}
	#[inline(always)] pub fn yxzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxzz_i(self.data)};}}
	#[inline(always)] pub fn yxzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxzw_i(self.data)};}}
	#[inline(always)] pub fn yxwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxwx_i(self.data)};}}
	#[inline(always)] pub fn yxwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxwy_i(self.data)};}}
	#[inline(always)] pub fn yxwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxwz_i(self.data)};}}
	#[inline(always)] pub fn yxww(self) -> Vector4Int { unsafe{return Vector4Int{data:_yxww_i(self.data)};}}
	#[inline(always)] pub fn yyxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyxx_i(self.data)};}}
	#[inline(always)] pub fn yyxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyxy_i(self.data)};}}
	#[inline(always)] pub fn yyxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyxz_i(self.data)};}}
	#[inline(always)] pub fn yyxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyxw_i(self.data)};}}
	#[inline(always)] pub fn yyyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyx_i(self.data)};}}
	#[inline(always)] pub fn yyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyy_i(self.data)};}}
	#[inline(always)] pub fn yyyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyz_i(self.data)};}}
	#[inline(always)] pub fn yyyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyyw_i(self.data)};}}
	#[inline(always)] pub fn yyzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyzx_i(self.data)};}}
	#[inline(always)] pub fn yyzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyzy_i(self.data)};}}
	#[inline(always)] pub fn yyzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyzz_i(self.data)};}}
	#[inline(always)] pub fn yyzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyzw_i(self.data)};}}
	#[inline(always)] pub fn yywx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yywx_i(self.data)};}}
	#[inline(always)] pub fn yywy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yywy_i(self.data)};}}
	#[inline(always)] pub fn yywz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yywz_i(self.data)};}}
	#[inline(always)] pub fn yyww(self) -> Vector4Int { unsafe{return Vector4Int{data:_yyww_i(self.data)};}}
	#[inline(always)] pub fn yzxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzxx_i(self.data)};}}
	#[inline(always)] pub fn yzxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzxy_i(self.data)};}}
	#[inline(always)] pub fn yzxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzxz_i(self.data)};}}
	#[inline(always)] pub fn yzxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzxw_i(self.data)};}}
	#[inline(always)] pub fn yzyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzyx_i(self.data)};}}
	#[inline(always)] pub fn yzyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzyy_i(self.data)};}}
	#[inline(always)] pub fn yzyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzyz_i(self.data)};}}
	#[inline(always)] pub fn yzyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzyw_i(self.data)};}}
	#[inline(always)] pub fn yzzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzzx_i(self.data)};}}
	#[inline(always)] pub fn yzzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzzy_i(self.data)};}}
	#[inline(always)] pub fn yzzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzzz_i(self.data)};}}
	#[inline(always)] pub fn yzzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzzw_i(self.data)};}}
	#[inline(always)] pub fn yzwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzwx_i(self.data)};}}
	#[inline(always)] pub fn yzwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzwy_i(self.data)};}}
	#[inline(always)] pub fn yzwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzwz_i(self.data)};}}
	#[inline(always)] pub fn yzww(self) -> Vector4Int { unsafe{return Vector4Int{data:_yzww_i(self.data)};}}
	#[inline(always)] pub fn ywxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywxx_i(self.data)};}}
	#[inline(always)] pub fn ywxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywxy_i(self.data)};}}
	#[inline(always)] pub fn ywxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywxz_i(self.data)};}}
	#[inline(always)] pub fn ywxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywxw_i(self.data)};}}
	#[inline(always)] pub fn ywyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywyx_i(self.data)};}}
	#[inline(always)] pub fn ywyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywyy_i(self.data)};}}
	#[inline(always)] pub fn ywyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywyz_i(self.data)};}}
	#[inline(always)] pub fn ywyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywyw_i(self.data)};}}
	#[inline(always)] pub fn ywzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywzx_i(self.data)};}}
	#[inline(always)] pub fn ywzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywzy_i(self.data)};}}
	#[inline(always)] pub fn ywzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywzz_i(self.data)};}}
	#[inline(always)] pub fn ywzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywzw_i(self.data)};}}
	#[inline(always)] pub fn ywwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywwx_i(self.data)};}}
	#[inline(always)] pub fn ywwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywwy_i(self.data)};}}
	#[inline(always)] pub fn ywwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywwz_i(self.data)};}}
	#[inline(always)] pub fn ywww(self) -> Vector4Int { unsafe{return Vector4Int{data:_ywww_i(self.data)};}}
	#[inline(always)] pub fn zxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxxx_i(self.data)};}}
	#[inline(always)] pub fn zxxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxxy_i(self.data)};}}
	#[inline(always)] pub fn zxxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxxz_i(self.data)};}}
	#[inline(always)] pub fn zxxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxxw_i(self.data)};}}
	#[inline(always)] pub fn zxyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxyx_i(self.data)};}}
	#[inline(always)] pub fn zxyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxyy_i(self.data)};}}
	#[inline(always)] pub fn zxyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxyz_i(self.data)};}}
	#[inline(always)] pub fn zxyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxyw_i(self.data)};}}
	#[inline(always)] pub fn zxzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxzx_i(self.data)};}}
	#[inline(always)] pub fn zxzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxzy_i(self.data)};}}
	#[inline(always)] pub fn zxzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxzz_i(self.data)};}}
	#[inline(always)] pub fn zxzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxzw_i(self.data)};}}
	#[inline(always)] pub fn zxwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxwx_i(self.data)};}}
	#[inline(always)] pub fn zxwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxwy_i(self.data)};}}
	#[inline(always)] pub fn zxwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxwz_i(self.data)};}}
	#[inline(always)] pub fn zxww(self) -> Vector4Int { unsafe{return Vector4Int{data:_zxww_i(self.data)};}}
	#[inline(always)] pub fn zyxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyxx_i(self.data)};}}
	#[inline(always)] pub fn zyxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyxy_i(self.data)};}}
	#[inline(always)] pub fn zyxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyxz_i(self.data)};}}
	#[inline(always)] pub fn zyxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyxw_i(self.data)};}}
	#[inline(always)] pub fn zyyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyyx_i(self.data)};}}
	#[inline(always)] pub fn zyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyyy_i(self.data)};}}
	#[inline(always)] pub fn zyyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyyz_i(self.data)};}}
	#[inline(always)] pub fn zyyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyyw_i(self.data)};}}
	#[inline(always)] pub fn zyzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyzx_i(self.data)};}}
	#[inline(always)] pub fn zyzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyzy_i(self.data)};}}
	#[inline(always)] pub fn zyzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyzz_i(self.data)};}}
	#[inline(always)] pub fn zyzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyzw_i(self.data)};}}
	#[inline(always)] pub fn zywx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zywx_i(self.data)};}}
	#[inline(always)] pub fn zywy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zywy_i(self.data)};}}
	#[inline(always)] pub fn zywz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zywz_i(self.data)};}}
	#[inline(always)] pub fn zyww(self) -> Vector4Int { unsafe{return Vector4Int{data:_zyww_i(self.data)};}}
	#[inline(always)] pub fn zzxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzxx_i(self.data)};}}
	#[inline(always)] pub fn zzxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzxy_i(self.data)};}}
	#[inline(always)] pub fn zzxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzxz_i(self.data)};}}
	#[inline(always)] pub fn zzxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzxw_i(self.data)};}}
	#[inline(always)] pub fn zzyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzyx_i(self.data)};}}
	#[inline(always)] pub fn zzyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzyy_i(self.data)};}}
	#[inline(always)] pub fn zzyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzyz_i(self.data)};}}
	#[inline(always)] pub fn zzyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzyw_i(self.data)};}}
	#[inline(always)] pub fn zzzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzzx_i(self.data)};}}
	#[inline(always)] pub fn zzzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzzy_i(self.data)};}}
	#[inline(always)] pub fn zzzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzzz_i(self.data)};}}
	#[inline(always)] pub fn zzzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzzw_i(self.data)};}}
	#[inline(always)] pub fn zzwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzwx_i(self.data)};}}
	#[inline(always)] pub fn zzwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzwy_i(self.data)};}}
	#[inline(always)] pub fn zzwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzwz_i(self.data)};}}
	#[inline(always)] pub fn zzww(self) -> Vector4Int { unsafe{return Vector4Int{data:_zzww_i(self.data)};}}
	#[inline(always)] pub fn zwxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwxx_i(self.data)};}}
	#[inline(always)] pub fn zwxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwxy_i(self.data)};}}
	#[inline(always)] pub fn zwxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwxz_i(self.data)};}}
	#[inline(always)] pub fn zwxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwxw_i(self.data)};}}
	#[inline(always)] pub fn zwyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwyx_i(self.data)};}}
	#[inline(always)] pub fn zwyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwyy_i(self.data)};}}
	#[inline(always)] pub fn zwyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwyz_i(self.data)};}}
	#[inline(always)] pub fn zwyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwyw_i(self.data)};}}
	#[inline(always)] pub fn zwzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwzx_i(self.data)};}}
	#[inline(always)] pub fn zwzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwzy_i(self.data)};}}
	#[inline(always)] pub fn zwzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwzz_i(self.data)};}}
	#[inline(always)] pub fn zwzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwzw_i(self.data)};}}
	#[inline(always)] pub fn zwwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwwx_i(self.data)};}}
	#[inline(always)] pub fn zwwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwwy_i(self.data)};}}
	#[inline(always)] pub fn zwwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwwz_i(self.data)};}}
	#[inline(always)] pub fn zwww(self) -> Vector4Int { unsafe{return Vector4Int{data:_zwww_i(self.data)};}}
	#[inline(always)] pub fn wxxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxxx_i(self.data)};}}
	#[inline(always)] pub fn wxxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxxy_i(self.data)};}}
	#[inline(always)] pub fn wxxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxxz_i(self.data)};}}
	#[inline(always)] pub fn wxxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxxw_i(self.data)};}}
	#[inline(always)] pub fn wxyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxyx_i(self.data)};}}
	#[inline(always)] pub fn wxyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxyy_i(self.data)};}}
	#[inline(always)] pub fn wxyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxyz_i(self.data)};}}
	#[inline(always)] pub fn wxyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxyw_i(self.data)};}}
	#[inline(always)] pub fn wxzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxzx_i(self.data)};}}
	#[inline(always)] pub fn wxzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxzy_i(self.data)};}}
	#[inline(always)] pub fn wxzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxzz_i(self.data)};}}
	#[inline(always)] pub fn wxzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxzw_i(self.data)};}}
	#[inline(always)] pub fn wxwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxwx_i(self.data)};}}
	#[inline(always)] pub fn wxwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxwy_i(self.data)};}}
	#[inline(always)] pub fn wxwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxwz_i(self.data)};}}
	#[inline(always)] pub fn wxww(self) -> Vector4Int { unsafe{return Vector4Int{data:_wxww_i(self.data)};}}
	#[inline(always)] pub fn wyxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyxx_i(self.data)};}}
	#[inline(always)] pub fn wyxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyxy_i(self.data)};}}
	#[inline(always)] pub fn wyxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyxz_i(self.data)};}}
	#[inline(always)] pub fn wyxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyxw_i(self.data)};}}
	#[inline(always)] pub fn wyyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyyx_i(self.data)};}}
	#[inline(always)] pub fn wyyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyyy_i(self.data)};}}
	#[inline(always)] pub fn wyyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyyz_i(self.data)};}}
	#[inline(always)] pub fn wyyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyyw_i(self.data)};}}
	#[inline(always)] pub fn wyzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyzx_i(self.data)};}}
	#[inline(always)] pub fn wyzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyzy_i(self.data)};}}
	#[inline(always)] pub fn wyzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyzz_i(self.data)};}}
	#[inline(always)] pub fn wyzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyzw_i(self.data)};}}
	#[inline(always)] pub fn wywx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wywx_i(self.data)};}}
	#[inline(always)] pub fn wywy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wywy_i(self.data)};}}
	#[inline(always)] pub fn wywz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wywz_i(self.data)};}}
	#[inline(always)] pub fn wyww(self) -> Vector4Int { unsafe{return Vector4Int{data:_wyww_i(self.data)};}}
	#[inline(always)] pub fn wzxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzxx_i(self.data)};}}
	#[inline(always)] pub fn wzxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzxy_i(self.data)};}}
	#[inline(always)] pub fn wzxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzxz_i(self.data)};}}
	#[inline(always)] pub fn wzxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzxw_i(self.data)};}}
	#[inline(always)] pub fn wzyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzyx_i(self.data)};}}
	#[inline(always)] pub fn wzyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzyy_i(self.data)};}}
	#[inline(always)] pub fn wzyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzyz_i(self.data)};}}
	#[inline(always)] pub fn wzyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzyw_i(self.data)};}}
	#[inline(always)] pub fn wzzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzzx_i(self.data)};}}
	#[inline(always)] pub fn wzzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzzy_i(self.data)};}}
	#[inline(always)] pub fn wzzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzzz_i(self.data)};}}
	#[inline(always)] pub fn wzzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzzw_i(self.data)};}}
	#[inline(always)] pub fn wzwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzwx_i(self.data)};}}
	#[inline(always)] pub fn wzwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzwy_i(self.data)};}}
	#[inline(always)] pub fn wzwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzwz_i(self.data)};}}
	#[inline(always)] pub fn wzww(self) -> Vector4Int { unsafe{return Vector4Int{data:_wzww_i(self.data)};}}
	#[inline(always)] pub fn wwxx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwxx_i(self.data)};}}
	#[inline(always)] pub fn wwxy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwxy_i(self.data)};}}
	#[inline(always)] pub fn wwxz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwxz_i(self.data)};}}
	#[inline(always)] pub fn wwxw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwxw_i(self.data)};}}
	#[inline(always)] pub fn wwyx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwyx_i(self.data)};}}
	#[inline(always)] pub fn wwyy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwyy_i(self.data)};}}
	#[inline(always)] pub fn wwyz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwyz_i(self.data)};}}
	#[inline(always)] pub fn wwyw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwyw_i(self.data)};}}
	#[inline(always)] pub fn wwzx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwzx_i(self.data)};}}
	#[inline(always)] pub fn wwzy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwzy_i(self.data)};}}
	#[inline(always)] pub fn wwzz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwzz_i(self.data)};}}
	#[inline(always)] pub fn wwzw(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwzw_i(self.data)};}}
	#[inline(always)] pub fn wwwx(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwwx_i(self.data)};}}
	#[inline(always)] pub fn wwwy(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwwy_i(self.data)};}}
	#[inline(always)] pub fn wwwz(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwwz_i(self.data)};}}
	#[inline(always)] pub fn wwww(self) -> Vector4Int { unsafe{return Vector4Int{data:_wwww_i(self.data)};}}
}


impl From<i32> for Vector4Int {
	#[inline(always)]
    fn from(v : i32) -> Vector4Int {
    	unsafe{
        return Vector4Int { data :_mm_set1_epi32(v)};	
    	}
    }
}
impl From<IntVector> for Vector4Int {
	#[inline(always)]
    fn from(val : IntVector) -> Vector4Int {
       unsafe{
			return Vector4Int{data : val.data};
		}
    }
}
impl From<Vector2Int> for Vector4Int {
	#[inline(always)]
    fn from(v : Vector2Int) -> Vector4Int {
    	unsafe{
        return Vector4Int { data : _mm_castps_si128(_mm_movelh_ps(_mm_castsi128_ps(v.data), _mm_setzero_ps() )) };
    	}
    }
}
impl From<Vector3Int> for Vector4Int {
	#[inline(always)]
    fn from(v : Vector3Int) -> Vector4Int {
    	unsafe{
       	 return Vector4Int { data : _mm_srli_si128 (_mm_slli_si128 (v.data, 4), 4) };
    	}
    }
}
impl From<Vector4> for Vector4Int {
	#[inline(always)]
    fn from(v : Vector4) -> Vector4Int {
    	unsafe{
        	return Vector4Int { data : _mm_cvttps_epi32(v.data) };
        }
    }
}


impl core::ops::Add for Vector4Int{
	type Output = Vector4Int;
	#[inline(always)]
	fn add(self, _rhs: Vector4Int) -> Vector4Int{
		return Vector4Int::add(self, _rhs);
	}
}
impl core::ops::AddAssign for Vector4Int {
	#[inline(always)]
    fn add_assign(&mut self, other: Vector4Int) {
        *self = Vector4Int::add(*self, other);
    }
}
impl core::ops::Sub for Vector4Int{
	type Output = Vector4Int;
	#[inline(always)]
	fn sub(self, _rhs: Vector4Int) -> Vector4Int{
		return Vector4Int::sub(self, _rhs);
	}
}
impl core::ops::SubAssign for Vector4Int {
	#[inline(always)]
    fn sub_assign(&mut self, other: Vector4Int) {
        *self = Vector4Int::sub(*self, other);
    }
}
impl core::ops::Neg for Vector4Int {
	type Output = Vector4Int;
	#[inline(always)]
	fn neg(self) -> Self::Output {
		unsafe{
			return Vector4Int{data:_mm_sub_epi32(_mm_set1_epi32(0), self.data)};
		}
	}
}

impl<T : Into<IntVector>> core::ops::Mul<T> for Vector4Int{
	type Output = Vector4Int;
	#[inline(always)]
	fn mul(self, _rhs: T) -> Vector4Int{
		return Vector4Int::mul(self,Vector4Int::from( _rhs.into()) );
	}
}
impl<T : Into<IntVector>> core::ops::MulAssign<T>  for Vector4Int{
	#[inline(always)]
	fn mul_assign(&mut self, _rhs: T){
		*self = Vector4Int::mul(*self, Vector4Int::from(_rhs.into()));
	}
}
impl core::ops::Mul<Vector4Int> for IntVector{
	type Output = Vector4Int;
	#[inline(always)]
	fn mul(self : IntVector, _rhs: Vector4Int) -> Vector4Int{
		return Vector4Int::mul(_rhs, Vector4Int::from(self));
	}
}

impl PartialEq for Vector4Int {
	#[inline(always)]
    fn eq(&self, other: &Vector4Int) -> bool {
    	return Vector4Int::equal(*self, *other).all();
    }
}

impl Hash for Vector4Int {
    fn hash<H: Hasher>(&self, state: &mut H) {
    	let mut dst = RawVector_i32{data:[0;4]};
    	unsafe{
    		let x : *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
    		_mm_store_si128(x, self.data);
    	}
        dst.data.hash(state);
    }
}
impl SIMDVector4 for Vector4Int{
#[inline(always)]
  fn data(self)->__m128{
  	return unsafe{_mm_castsi128_ps (self.data)};
  }
  fn data_i(self)->__m128i{
  	return unsafe{ self.data};
  }
}