mod vector2;
mod vector3;
mod vector4;
mod vector2_int;
mod vector3_int;
mod vector4_int;
mod quaternion;
mod dual_quaternion;
use std::arch::x86_64::*;
const NORMALIZATION_EPSILON : f32 = 1.0e-10;


#[derive(Copy, Clone)]
#[repr(C, align(16))]
pub struct RawVec{
	data : [f32;4],
}

#[derive(Copy, Clone)]
pub struct Vector2{
	data : __m128,
}

/// A vector 3
#[derive(Copy, Clone)]
pub struct Vector3{
	data : __m128,
}

#[derive(Copy, Clone)]
pub struct Vector4{
	data : __m128,
}

#[derive(Copy, Clone)]
pub struct Vector2Int{
	data : __m128i,
}
#[derive(Copy, Clone)]
pub struct Vector3Int{
	data : __m128i,
}
#[derive(Copy, Clone)]
pub struct Vector4Int{
	data : __m128i,
}

#[derive(Copy, Clone)]
pub struct Quaternion{
	data : __m128,
}

#[derive(Copy, Clone)]
pub struct DualQuaternion{
	real : __m128,
	dual : __m128,
}

#[derive(Copy, Clone)]
pub struct Matrix4x4{
	m : [__m128; 4],
}

/// A replacement for the _MM_SHUFFLE macro.
const fn _ico_shuffle(z : u32, y : u32, x : u32, w : u32) -> u32 {
    (((z) << 6) | ((y) << 4) | ((x) << 2) | (w))
}

#[inline(always)]
unsafe fn _xxxx(m : __m128)->__m128{
	return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,0,0));
}
#[inline(always)]
unsafe fn _yyyy(m : __m128)->__m128{
	return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,1,1));
}
#[inline(always)]
unsafe fn _zzzz(m : __m128)->__m128{
	return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,2,2));
}
#[inline(always)]
unsafe fn _wwww(m : __m128)->__m128{
	return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,3,3));
}
#[inline(always)]
unsafe fn _xzyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 0));
}
#[inline(always)]
unsafe fn _yxzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 1));
}
#[inline(always)]
unsafe fn _yzxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 1));
}
#[inline(always)]
unsafe fn _zxyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 2));
}
#[inline(always)]
unsafe fn _zxwy(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 0, 2));
}
#[inline(always)]
unsafe fn _zyxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 2));
}

#[inline(always)]
unsafe fn _xxyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 0));
}
#[inline(always)]
unsafe fn _xxzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 0));
}
#[inline(always)]
unsafe fn _xyxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 0));
}
#[inline(always)]
unsafe fn _xzxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 0));
}
#[inline(always)]
unsafe fn _yxxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 1));
}
#[inline(always)]
unsafe fn _zxxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 2));
}

#[inline(always)]
unsafe fn _yyxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 1));
}
#[inline(always)]
unsafe fn _yyzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 1));
}
#[inline(always)]
unsafe fn _yxyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 1));
}
#[inline(always)]
unsafe fn _yzyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 1));
}
#[inline(always)]
unsafe fn _xyyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 0));
}
#[inline(always)]
unsafe fn _zyyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 2));
}
#[inline(always)]
unsafe fn _ywxz(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 3, 1));
}

#[inline(always)]
unsafe fn _zzxw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 2));
}
#[inline(always)]
unsafe fn _zzyw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 2));
}
#[inline(always)]
unsafe fn _zxzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 2));
}
#[inline(always)]
unsafe fn _zyzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 2));
}
#[inline(always)]
unsafe fn _xzzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 0));
}
#[inline(always)]
unsafe fn _yzzw(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 1));
}


//Others
#[inline(always)]
unsafe fn _zwxy(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 3, 2));
}
#[inline(always)]
unsafe fn _zwyx(m : __m128)->__m128{
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 3, 2));
}

#[inline(always)]
unsafe fn _yxwz(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 0, 1));
}

#[inline(always)]
unsafe fn _wzxy(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 2, 3));
}
#[inline(always)]
unsafe fn _wzyx(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 2, 3));
}
#[inline(always)]
unsafe fn _wxyz(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 0, 3));
}
#[inline(always)]
unsafe fn _yzwx(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 2, 1));
}
#[inline(always)]
unsafe fn _xyzz(m : __m128)->__m128{
	return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 1, 0));
}

// Generators for common vectors to avoid set.
// https://www.agner.org/optimize/optimizing_assembly.pdf

#[inline(always)]
unsafe fn  _ico_one_epi32() -> __m128i{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return (_mm_srli_epi32(a, 31));
}

#[inline(always)]
unsafe fn  _ico_half_ps()->__m128{
   let mut  a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 26);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
unsafe fn  _ico_one_ps()->__m128{
    //return  _mm_set1_epi32(0x3f800000);
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 25);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
unsafe fn  _ico_two_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 31);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 1));
}
#[inline(always)]
unsafe fn   _ico_nearesttwo_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
unsafe fn   _ico_negtwo_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_slli_epi32(a, 30));
}
#[inline(always)]
unsafe fn _ico_signbit_ps() ->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_slli_epi32(a, 31));
}
#[inline(always)]
unsafe fn  _ico_abs_ps(a : __m128 )->__m128{
	return _mm_andnot_ps(_ico_signbit_ps(), a);
}
#[inline(always)]
unsafe fn _ico_copysign_ps(a : __m128, b : __m128) ->__m128{
    // When we xor the sign bit, both positive or negative = 0.  If they differ, a's sign bit is set to true.
    // But then we mask the sign bit off with the absolute value.
	let tmp = _ico_abs_ps(_mm_xor_ps(a, b));
    // When we xor them back together, if b is negative, it's value comes through. A always 0.
    // 0 ^ 1 = 1, 0 ^ 0 = 0
    
	return _mm_xor_ps(tmp, b);
}
#[inline(always)]
unsafe fn _ico_floorps_epi32(a : __m128)->__m128i{
    ///convert to int using truncate
    let integer = _mm_cvttps_epi32(a);
    // if the truncated value is greater than the old value (negative) subtract 1.
    return _mm_sub_epi32(integer, 
		_mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(_mm_cvtepi32_ps(integer), a)), _ico_one_epi32()));
}
#[inline(always)]
unsafe fn _ico_ceilps_epi32(a : __m128)->__m128i{
    ///convert to int using truncate
    let integer = _mm_cvttps_epi32(a);
    // if the truncated value is greater than the old value (negative) subtract 1.
    return _mm_add_epi32(integer, _mm_and_si128(
		_mm_castps_si128(_mm_cmplt_ps(_mm_cvtepi32_ps(integer), a)), _ico_one_epi32()));
}

#[inline(always)]
unsafe fn _ico_select_ps( a : __m128,   b : __m128,  mask : __m128)->__m128{
    return _mm_blendv_ps(a,b,mask);
	// (((b ^ a) & mask)^a)
 
    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/  
    //return _mm_or_ps(_mm_and_ps(a, mask), _mm_andnot_ps(mask, b));
	//return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}
#[inline(always)]
unsafe fn _ico_select_pd( a : __m128d, b : __m128d,  mask : __m128d)->__m128d{
    return _mm_blendv_pd(a,b,mask);
    // (((b ^ a) & mask)^a)
    
    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
    //return _mm_or_pd(_mm_and_pd(b, mask), _mm_andnot_pd(mask, a));
    //return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}

#[inline(always)]
unsafe fn  _ico_select_si128(a : __m128i, b : __m128i,  mask : __m128i) ->__m128i{
	// (((b ^ a) & mask)^a)
    return _mm_or_si128(_mm_and_si128(b, mask), _mm_andnot_si128(mask, a));
	//return _mm_xor_si128( a, _mm_and_si128( mask, _mm_xor_si128 ( b, a ) ) );
	//_mm_or_si128(_mm_and_si128(a, cond), _mm_andnot_si128(cond, b))
}


///If the float is entirely integer or greater - we return the float
/// Otherwise we cast to int
#[inline(always)]
unsafe fn _ico_truncate_ps(a : __m128) -> __m128{
    

    /// Largest float < integer only is  < 8388608
    /// NaN and int-only floats fail this test.
    let mask = _mm_cmplt_ps(_ico_abs_ps(a), _mm_castsi128_ps(_mm_set1_epi32(0x4b000000))); /// is the float entirely integer?
    let trunc = _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
                                  //if Abs(a) > 8388608
    // select a if greater then 8388608.0f, otherwise select the result of FuncT
    // NaN will fail and get passed through.
    return _ico_select_ps(a,trunc, mask);
}

///Round away from zero, not toward zero to match cmath.
//http://dss.stephanierct.com/DevBlog/?p=8
#[inline(always)]
unsafe fn _ico_round_ps(a : __m128) ->__m128{

    ///convert to int using truncate
    let trunc = _ico_truncate_ps(a);//_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    let remainder = _mm_sub_ps(a, trunc);
    /// If we want to round toward zero, use _ico_nearesttwo_ps()
    let rmd2 = _mm_mul_ps( remainder, _ico_two_ps());
    let rmd2i = _mm_cvttps_epi32(rmd2);    // after being truncated of course
    
    return _mm_add_ps(trunc, _mm_cvtepi32_ps(rmd2i));
}
#[inline(always)]
unsafe fn _ico_cross_ps(lhs : __m128, rhs : __m128) ->__m128{
    // x  <-  a.y*b.z - a.z*b.y
    // y  <-  a.z*b.x - a.x*b.z
    // z  <-  a.x*b.y - a.y*b.x
    // We can save a shuffle by grouping it in this wacky order:
    let mut tmp0 =  _mm_mul_ps(lhs, _mm_shuffle_ps(rhs, rhs, _ico_shuffle(3, 1, 0, 2)));
    //__m128 tmp1 =  _mm_mul_ps(_mm_shuffle_ps(lhs, lhs, _MM_SHUFFLE(3, 1, 0, 2)), rhs);
    //tmp0 = _mm_sub_ps(tmp1, tmp0);
    tmp0 = _mm_fmsub_ps(_mm_shuffle_ps(lhs, lhs, _ico_shuffle(3, 1, 0, 2)), rhs, tmp0);
    
    return _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(3, 1, 0, 2));
}

#[inline(always)]
unsafe fn _ico_dp4_ps( v0 : __m128,  v1 : __m128) ->__m128{
	let mut tmp0 = _mm_mul_ps(v0, v1);
	let mut tmp1 = _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(2, 3, 0, 1)); //yxwz
	tmp0 = _mm_add_ps(tmp0, tmp1);//xy, xy, zw, zw
	tmp1 = _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(0, 1, 2, 3)); //wxyz
	return _mm_add_ps(tmp0, tmp1);
}

#[inline(always)]
unsafe fn _ico_quat_mul(lhs : __m128, rhs : __m128) -> __m128{
    let x = _mm_mul_ps(_wzyx(rhs), _xxxx(lhs));
    let xw = _mm_fmsubadd_ps(_wwww(lhs), rhs, x);
    let y = _mm_xor_ps( _yyyy(lhs), _mm_set_ps(-0.0,-0.0,0.0,0.0));
    let xyw = _mm_fmadd_ps(_zwxy(rhs), y, xw);
    let z = _mm_xor_ps( _zzzz(lhs), _mm_set_ps(-0.0,0.0,0.0,-0.0));
    return _mm_fmadd_ps(_yxwz(rhs),z, xyw);
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {


    	let q = crate::Mesh::new();
        let v = q.clone();
       
       
    	use crate::Vector3;
    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);


		let c = a + b;
        assert_eq!(c.x(), 5.0);

        
    }
}
