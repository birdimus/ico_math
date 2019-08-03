

use core::arch::x86_64::*;
pub const NORMALIZATION_EPSILON : f32 = 1.0e-10;
pub const SIGN_BIT : f32 = -0.0;
pub const EPSILON_AT_ONE : f32 = 0.00000012;
pub const INV_TWO_PI : f32 = 0.159154943091895335768883763372514362034459645740456448747;
pub const TWO_PI : f32 = 6.283185307179586476925286766559005768394338798750211641949;
pub const HALF_PI : f32 = 1.57079632679;
pub const PI : f32 = 3.14159265359;

/// A replacement for the _MM_SHUFFLE macro.
pub const fn _ico_shuffle(z : i32, y : i32, x : i32, w : i32) -> i32 {
    (((z) << 6) | ((y) << 4) | ((x) << 2) | (w))
}


// Generators for common vectors to avoid set.
// https://www.agner.org/optimize/optimizing_assembly.pdf

#[inline(always)]
pub unsafe fn  _ico_one_epi32() -> __m128i{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return (_mm_srli_epi32(a, 31));
}

#[inline(always)]
pub unsafe fn  _ico_half_ps()->__m128{
   let mut  a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 26);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
pub unsafe fn  _ico_one_ps()->__m128{
    //return  _mm_set1_epi32(0x3f800000);
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 25);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
pub unsafe fn  _ico_two_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   a = _mm_slli_epi32(a, 31);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 1));
}
#[inline(always)]
pub unsafe fn   _ico_nearesttwo_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
pub unsafe fn   _ico_negtwo_ps()->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_slli_epi32(a, 30));
}
#[inline(always)]
pub unsafe fn _ico_signbit_ps() ->__m128{
   let mut a = _mm_setzero_si128();
   a = _mm_cmpeq_epi32 (a, a);
   return _mm_castsi128_ps(_mm_slli_epi32(a, 31));
}
#[inline(always)]
pub unsafe fn  _ico_abs_ps(a : __m128 )->__m128{
	return _mm_andnot_ps(_ico_signbit_ps(), a);
}
#[inline(always)]
pub unsafe fn _ico_copysign_ps(a : __m128, b : __m128) ->__m128{
    // When we xor the sign bit, both positive or negative = 0.  If they differ, a's sign bit is set to true.
    // But then we mask the sign bit off with the absolute value.
	let tmp = _ico_abs_ps(_mm_xor_ps(a, b));
    // When we xor them back together, if b is negative, it's value comes through. A always 0.
    // 0 ^ 1 = 1, 0 ^ 0 = 0
    
	return _mm_xor_ps(tmp, b);
}


#[inline(always)]
pub unsafe fn _ico_select_ps( a : __m128,   b : __m128,  mask : __m128)->__m128{
    return _mm_blendv_ps(a,b,mask);
	// (((b ^ a) & mask)^a)
 
    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/  
    //return _mm_or_ps(_mm_and_ps(a, mask), _mm_andnot_ps(mask, b));
	//return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}
#[inline(always)]
pub unsafe fn _ico_select_pd( a : __m128d, b : __m128d,  mask : __m128d)->__m128d{
    return _mm_blendv_pd(a,b,mask);
    // (((b ^ a) & mask)^a)
    
    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
    //return _mm_or_pd(_mm_and_pd(b, mask), _mm_andnot_pd(mask, a));
    //return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}

#[inline(always)]
pub unsafe fn  _ico_select_si128(a : __m128i, b : __m128i,  mask : __m128i) ->__m128i{
	// (((b ^ a) & mask)^a)
    return _mm_or_si128(_mm_and_si128(b, mask), _mm_andnot_si128(mask, a));
	//return _mm_xor_si128( a, _mm_and_si128( mask, _mm_xor_si128 ( b, a ) ) );
	//_mm_or_si128(_mm_and_si128(a, cond), _mm_andnot_si128(cond, b))
}


///If the float is entirely integer or greater - we return the float
/// Otherwise we cast to int
#[inline(always)]
pub unsafe fn _ico_truncate_ps(a : __m128) -> __m128{
    

    /// Largest float < integer only is  < 8388608
    /// NaN and int-only floats fail this test.
    let mask = _mm_cmplt_ps(_ico_abs_ps(a), _mm_castsi128_ps(_mm_set1_epi32(0x4b000000))); /// is the float entirely integer?
    let trunc = _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
                                  //if Abs(a) > 8388608
    // select a if greater then 8388608.0f, otherwise select the result of FuncT
    // NaN will fail and get passed through.
    return _ico_select_ps(a,trunc, mask);
}




#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_floor_ps(a : __m128) -> __m128{
    
    ///convert to int using truncate
    let t = _ico_truncate_ps(a);//_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    // if the truncated value is greater than the old value (negative, non integral) subtract 1.
    return _mm_sub_ps(t, _mm_and_ps(_mm_cmpgt_ps(t, a), _ico_one_ps()));
}

#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_ceil_ps(a : __m128) -> __m128{
    ///convert to int using truncate
    let t = _ico_truncate_ps(a);//_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    // if the truncated value is less than the old value (positive) add 1.
    return _mm_add_ps(t, _mm_and_ps(_mm_cmplt_ps(t, a), _ico_one_ps()));
}

///Round away from zero, not toward zero to match cmath.
/// Returns unsigned 0.
//http://dss.stephanierct.com/DevBlog/?p=8
#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_round_ps(a : __m128) -> __m128{

    ///convert to int using truncate
    let trunc = _ico_truncate_ps(a);//_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    let remainder = _mm_sub_ps(a, trunc);
    /// If we want to round toward zero, use _ico_nearesttwo_ps()
    let rmd2 = _mm_mul_ps( remainder, _ico_two_ps());
    let rmd2i = _mm_cvttps_epi32(rmd2);    // after being truncated of course
    
    return _mm_add_ps(trunc, _mm_cvtepi32_ps(rmd2i));
}

#[inline(always)]
pub unsafe fn _ico_cross_ps(lhs : __m128, rhs : __m128) ->__m128{
    // x  <-  a.y*b.z - a.z*b.y
    // y  <-  a.z*b.x - a.x*b.z
    // z  <-  a.x*b.y - a.y*b.x
    // We can save a shuffle by grouping it in this wacky order:
    let mut tmp0 =  _mm_mul_ps(lhs, _zxyw(rhs));
    //__m128 tmp1 =  _mm_mul_ps(_mm_shuffle_ps(lhs, lhs, _MM_SHUFFLE(3, 1, 0, 2)), rhs);
    //tmp0 = _mm_sub_ps(tmp1, tmp0);
    tmp0 = _mm_fmsub_ps(_zxyw(lhs), rhs, tmp0);
    
    return _zxyw(tmp0);// _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(3, 1, 0, 2));
}

#[inline(always)]
pub unsafe fn _ico_dp4_ps( v0 : __m128,  v1 : __m128) ->__m128{
	let mut tmp0 = _mm_mul_ps(v0, v1);
	let mut tmp1 = _yxwz(tmp0);// _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(2, 3, 0, 1)); //yxwz
	tmp0 = _mm_add_ps(tmp0, tmp1);//xy, xy, zw, zw
	tmp1 = _wzyx(tmp0);//_mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(0, 1, 2, 3)); //wxyz
	return _mm_add_ps(tmp0, tmp1);
}

#[inline(always)]
pub unsafe fn _ico_quat_mul(lhs : __m128, rhs : __m128) -> __m128{
    let x = _mm_mul_ps(_wzyx(rhs), _xxxx(lhs));
    let xw = _mm_fmsubadd_ps(_wwww(lhs), rhs, x);
    let y = _mm_xor_ps( _yyyy(lhs), _mm_set_ps(-0.0,-0.0,0.0,0.0));
    let xyw = _mm_fmadd_ps(_zwxy(rhs), y, xw);
    let z = _mm_xor_ps( _zzzz(lhs), _mm_set_ps(-0.0,0.0,0.0,-0.0));
    return _mm_fmadd_ps(_yxwz(rhs),z, xyw);
}

#[inline(always)]
/// Bounce between 0 and 1
pub unsafe fn _ico_ping_pong(vec : __m128 ) -> __m128{
  // 2 * abs(x + 0.5 - floor(x + 0.5) - 0.5)
  //  2 * abs(x - floor(x + 0.5))
    let shifted = _mm_add_ps(vec, _ico_half_ps());
    let internal = _ico_abs_ps(_mm_sub_ps(vec, _mm_floor_ps(shifted)));
    return _mm_add_ps(internal,internal); // multiply by 2
}


#[inline(always)]
pub unsafe fn _ico_approx_cos01(vec : __m128 ) -> __m128{
    
    ///  2 approximations were used to compute this result
    ///  A = 4x^3 - 6x^2 + 1
    ///  B = -2x^5 + 5x^4 - 5x^2 + 1
    /// The range was cut in half
    ///  A = 2x^3 - 3x^2 + 1
    ///  B = -x^5 + 2.5x^4 - 2.5x^2 + 1
    ///  These were weighted (1-t)*A + (t)*B to minimize error and ensure 0,1,0.25, 0.75, and 0.5 provided exact solutions
    /// Max absolute error on domain [0,1], range [-1,1] is 0.000192f at ~0.25 +- 0.125
    /// Max relative error on domain [0,1]. range [-1,1] approaches 0.00069f at 0.5
    
    
    
    /// 1.115408f Was found empirically as T.
    //// coeff for {x^5, x^4, x^3, x^2}
    
    //// Next, coefficients had to be shifted to account for limitations in floating point representation
    //// COMPUTED SCALARS: (-2.230816f,5.57704f, -0.461632f, -4.884592f );
    ///  x^4, x^3 coefficients were adjusted
    let scalars = _mm_set_ps(-2.230816,5.5770397, -0.4616319, -4.884592 );
    //let scalars = _mm_castsi128_ps(_mm_set_epi32( 0xc00ec5b0, 0x40b2771c, 0xbeec5b04, 0xc09c4e94));
    /// This is a series approximation of cos.  Error < 0.01% in the domain
    let vec2 = _mm_mul_ps(vec, vec);
    let mut result = _mm_fmadd_ps(vec2, _xxxx(scalars), _ico_one_ps());
    
    let vec3 = _mm_mul_ps(vec2, vec);
    result = _mm_fmadd_ps(vec3, _yyyy(scalars), result);
    
    let vec4 = _mm_mul_ps(vec3, vec);
    result = _mm_fmadd_ps(vec4, _zzzz(scalars), result);
    
    let vec5 = _mm_mul_ps(vec4, vec);
    result = _mm_fmadd_ps(vec5, _wwww(scalars), result);

    // As a final step, clamp to between -1 and 1 for safety.
    let one = _ico_one_ps();
    let mask = _mm_cmpgt_ps(_ico_abs_ps(result), one);
    let clamped = _ico_copysign_ps(one, result);
    return _ico_select_ps(result, clamped, mask);
}

#[inline(always)]
pub unsafe fn _ico_cos_ps(vec : __m128 ) -> __m128{
  //let scaled = _mm_div_ps(vec, _mm_set1_ps(TWO_PI));
  let scaled = _mm_mul_ps(vec, _mm_set1_ps(INV_TWO_PI));
  let looped = _ico_ping_pong(scaled);
  return _ico_approx_cos01(looped);
}
#[inline(always)]
pub unsafe fn _ico_cos_deg_ps(vec : __m128 ) -> __m128{
  let scaled = _mm_div_ps(vec, _mm_set1_ps(360.0));
  let looped = _ico_ping_pong(scaled);
  return _ico_approx_cos01(looped);
}
#[inline(always)]
pub unsafe fn _ico_sin_ps(vec : __m128 ) -> __m128{
  //let scaled = _mm_div_ps(vec, _mm_set1_ps(TWO_PI));
  let scaled = _mm_mul_ps(vec, _mm_set1_ps(INV_TWO_PI));

  let shifted = _mm_add_ps(scaled, _ico_half_ps());
  let internal = _mm_sub_ps(scaled, _mm_floor_ps(shifted));
  let sin_shift = _mm_sub_ps(internal, _mm_set1_ps(0.25));

  let looped = _ico_ping_pong(sin_shift);
  return _ico_approx_cos01(looped);
}
#[inline(always)]
pub unsafe fn _ico_sin_deg_ps(vec : __m128 ) -> __m128{
  let scaled = _mm_div_ps(vec, _mm_set1_ps(360.0));

  let shifted = _mm_add_ps(scaled, _ico_half_ps());
  let internal = _mm_sub_ps(scaled, _mm_floor_ps(shifted));
  let sin_shift = _mm_sub_ps(internal, _mm_set1_ps(0.25));

  let looped = _ico_ping_pong(sin_shift);
  return _ico_approx_cos01(looped);
}

#[inline(always)]
pub unsafe fn _ico_tan_ps(vec : __m128 ) -> __m128{
  return _mm_div_ps(_ico_sin_ps(vec), _ico_cos_ps(vec));
}

#[inline(always)]
/// max error 0.000001 radians
/// 4th order approximation
/// acos(x)/sqrt(1-x) goes from pi/2 at 0 to sqrt(2) at 1 ==  0 to -0.15658276442
pub unsafe fn _ico_acos_ps(vec : __m128 ) -> __m128{

  
  let abs_vec = _ico_abs_ps(vec);
  let scalars = _mm_set_ps(-0.2142619,0.08528893,-0.03673658,0.009126827);
  let mut result = _mm_fmadd_ps(_xxxx(scalars), abs_vec, _yyyy(scalars));
  result = _mm_fmadd_ps(result, abs_vec, _zzzz(scalars));
  result = _mm_fmadd_ps(result, abs_vec, _wwww(scalars));
  result = _mm_fmadd_ps(result, abs_vec, _mm_set1_ps(HALF_PI));
  let sqrt_vec = _mm_sqrt_ps(_mm_sub_ps(_ico_one_ps(), abs_vec));
  result = _mm_mul_ps(result, sqrt_vec);

  let inv_mask = _mm_cmplt_ps(vec, _mm_setzero_ps());
  result = _mm_xor_ps(result, _mm_and_ps(_ico_signbit_ps(), inv_mask));
  result = _mm_add_ps(result, _mm_and_ps(inv_mask,  _mm_set1_ps(PI)));
  return result;
}
pub unsafe fn _ico_asin_ps(vec : __m128 ) -> __m128{
  return _mm_sub_ps(_mm_set1_ps(HALF_PI), _ico_acos_ps(vec));
}





#[inline(always)]
pub unsafe fn _xxxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,0,0));
}
#[inline(always)]
pub unsafe fn _xxxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,0,0));
}
#[inline(always)]
pub unsafe fn _xxxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,0,0));
}
#[inline(always)]
pub unsafe fn _xxxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,0,0));
}

#[inline(always)]
pub unsafe fn _xxyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,0,0));
}
#[inline(always)]
pub unsafe fn _xxyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,0,0));
}
#[inline(always)]
pub unsafe fn _xxyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,0,0));
}
#[inline(always)]
pub unsafe fn _xxyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,0,0));
}

#[inline(always)]
pub unsafe fn _xxzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,0,0));
}
#[inline(always)]
pub unsafe fn _xxzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,0,0));
}
#[inline(always)]
pub unsafe fn _xxzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,0,0));
}
#[inline(always)]
pub unsafe fn _xxzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,0,0));
}

#[inline(always)]
pub unsafe fn _xxwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,0,0));
}
#[inline(always)]
pub unsafe fn _xxwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,0,0));
}
#[inline(always)]
pub unsafe fn _xxwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,0,0));
}
#[inline(always)]
pub unsafe fn _xxww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,0,0));
}


#[inline(always)]
pub unsafe fn _xyxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,1,0));
}
#[inline(always)]
pub unsafe fn _xyxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,1,0));
}
#[inline(always)]
pub unsafe fn _xyxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,1,0));
}
#[inline(always)]
pub unsafe fn _xyxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,1,0));
}

#[inline(always)]
pub unsafe fn _xyyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,1,0));
}
#[inline(always)]
pub unsafe fn _xyyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,1,0));
}
#[inline(always)]
pub unsafe fn _xyyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,1,0));
}
#[inline(always)]
pub unsafe fn _xyyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,1,0));
}

#[inline(always)]
pub unsafe fn _xyzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,1,0));
}
#[inline(always)]
pub unsafe fn _xyzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,1,0));
}
#[inline(always)]
pub unsafe fn _xyzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,1,0));
}
#[inline(always)]
pub unsafe fn _xyzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,1,0));
}

#[inline(always)]
pub unsafe fn _xywx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,1,0));
}
#[inline(always)]
pub unsafe fn _xywy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,1,0));
}
#[inline(always)]
pub unsafe fn _xywz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,1,0));
}
#[inline(always)]
pub unsafe fn _xyww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,1,0));
}



#[inline(always)]
pub unsafe fn _xzxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,2,0));
}
#[inline(always)]
pub unsafe fn _xzxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,2,0));
}
#[inline(always)]
pub unsafe fn _xzxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,2,0));
}
#[inline(always)]
pub unsafe fn _xzxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,2,0));
}

#[inline(always)]
pub unsafe fn _xzyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,2,0));
}
#[inline(always)]
pub unsafe fn _xzyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,2,0));
}
#[inline(always)]
pub unsafe fn _xzyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,2,0));
}
#[inline(always)]
pub unsafe fn _xzyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,2,0));
}

#[inline(always)]
pub unsafe fn _xzzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,2,0));
}
#[inline(always)]
pub unsafe fn _xzzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,2,0));
}
#[inline(always)]
pub unsafe fn _xzzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,2,0));
}
#[inline(always)]
pub unsafe fn _xzzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,2,0));
}

#[inline(always)]
pub unsafe fn _xzwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,2,0));
}
#[inline(always)]
pub unsafe fn _xzwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,2,0));
}
#[inline(always)]
pub unsafe fn _xzwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,2,0));
}
#[inline(always)]
pub unsafe fn _xzww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,2,0));
}

#[inline(always)]
pub unsafe fn _xwxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,3,0));
}
#[inline(always)]
pub unsafe fn _xwxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,3,0));
}
#[inline(always)]
pub unsafe fn _xwxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,3,0));
}
#[inline(always)]
pub unsafe fn _xwxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,3,0));
}

#[inline(always)]
pub unsafe fn _xwyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,3,0));
}
#[inline(always)]
pub unsafe fn _xwyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,3,0));
}
#[inline(always)]
pub unsafe fn _xwyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,3,0));
}
#[inline(always)]
pub unsafe fn _xwyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,3,0));
}

#[inline(always)]
pub unsafe fn _xwzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,3,0));
}
#[inline(always)]
pub unsafe fn _xwzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,3,0));
}
#[inline(always)]
pub unsafe fn _xwzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,3,0));
}
#[inline(always)]
pub unsafe fn _xwzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,3,0));
}

#[inline(always)]
pub unsafe fn _xwwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,3,0));
}
#[inline(always)]
pub unsafe fn _xwwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,3,0));
}
#[inline(always)]
pub unsafe fn _xwwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,3,0));
}
#[inline(always)]
pub unsafe fn _xwww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,3,0));
}

#[inline(always)]
pub unsafe fn _yxxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,0,1));
}
#[inline(always)]
pub unsafe fn _yxxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,0,1));
}
#[inline(always)]
pub unsafe fn _yxxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,0,1));
}
#[inline(always)]
pub unsafe fn _yxxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,0,1));
}

#[inline(always)]
pub unsafe fn _yxyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,0,1));
}
#[inline(always)]
pub unsafe fn _yxyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,0,1));
}
#[inline(always)]
pub unsafe fn _yxyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,0,1));
}
#[inline(always)]
pub unsafe fn _yxyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,0,1));
}

#[inline(always)]
pub unsafe fn _yxzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,0,1));
}
#[inline(always)]
pub unsafe fn _yxzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,0,1));
}
#[inline(always)]
pub unsafe fn _yxzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,0,1));
}
#[inline(always)]
pub unsafe fn _yxzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,0,1));
}

#[inline(always)]
pub unsafe fn _yxwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,0,1));
}
#[inline(always)]
pub unsafe fn _yxwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,0,1));
}
#[inline(always)]
pub unsafe fn _yxwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,0,1));
}
#[inline(always)]
pub unsafe fn _yxww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,0,1));
}


#[inline(always)]
pub unsafe fn _yyxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,1,1));
}
#[inline(always)]
pub unsafe fn _yyxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,1,1));
}
#[inline(always)]
pub unsafe fn _yyxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,1,1));
}
#[inline(always)]
pub unsafe fn _yyxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,1,1));
}

#[inline(always)]
pub unsafe fn _yyyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,1,1));
}
#[inline(always)]
pub unsafe fn _yyyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,1,1));
}
#[inline(always)]
pub unsafe fn _yyyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,1,1));
}
#[inline(always)]
pub unsafe fn _yyyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,1,1));
}

#[inline(always)]
pub unsafe fn _yyzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,1,1));
}
#[inline(always)]
pub unsafe fn _yyzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,1,1));
}
#[inline(always)]
pub unsafe fn _yyzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,1,1));
}
#[inline(always)]
pub unsafe fn _yyzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,1,1));
}

#[inline(always)]
pub unsafe fn _yywx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,1,1));
}
#[inline(always)]
pub unsafe fn _yywy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,1,1));
}
#[inline(always)]
pub unsafe fn _yywz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,1,1));
}
#[inline(always)]
pub unsafe fn _yyww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,1,1));
}



#[inline(always)]
pub unsafe fn _yzxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,2,1));
}
#[inline(always)]
pub unsafe fn _yzxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,2,1));
}
#[inline(always)]
pub unsafe fn _yzxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,2,1));
}
#[inline(always)]
pub unsafe fn _yzxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,2,1));
}

#[inline(always)]
pub unsafe fn _yzyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,2,1));
}
#[inline(always)]
pub unsafe fn _yzyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,2,1));
}
#[inline(always)]
pub unsafe fn _yzyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,2,1));
}
#[inline(always)]
pub unsafe fn _yzyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,2,1));
}

#[inline(always)]
pub unsafe fn _yzzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,2,1));
}
#[inline(always)]
pub unsafe fn _yzzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,2,1));
}
#[inline(always)]
pub unsafe fn _yzzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,2,1));
}
#[inline(always)]
pub unsafe fn _yzzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,2,1));
}

#[inline(always)]
pub unsafe fn _yzwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,2,1));
}
#[inline(always)]
pub unsafe fn _yzwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,2,1));
}
#[inline(always)]
pub unsafe fn _yzwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,2,1));
}
#[inline(always)]
pub unsafe fn _yzww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,2,1));
}

#[inline(always)]
pub unsafe fn _ywxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,3,1));
}
#[inline(always)]
pub unsafe fn _ywxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,3,1));
}
#[inline(always)]
pub unsafe fn _ywxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,3,1));
}
#[inline(always)]
pub unsafe fn _ywxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,3,1));
}

#[inline(always)]
pub unsafe fn _ywyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,3,1));
}
#[inline(always)]
pub unsafe fn _ywyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,3,1));
}
#[inline(always)]
pub unsafe fn _ywyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,3,1));
}
#[inline(always)]
pub unsafe fn _ywyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,3,1));
}

#[inline(always)]
pub unsafe fn _ywzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,3,1));
}
#[inline(always)]
pub unsafe fn _ywzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,3,1));
}
#[inline(always)]
pub unsafe fn _ywzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,3,1));
}
#[inline(always)]
pub unsafe fn _ywzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,3,1));
}

#[inline(always)]
pub unsafe fn _ywwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,3,1));
}
#[inline(always)]
pub unsafe fn _ywwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,3,1));
}
#[inline(always)]
pub unsafe fn _ywwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,3,1));
}
#[inline(always)]
pub unsafe fn _ywww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,3,1));
}

#[inline(always)]
pub unsafe fn _zxxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,0,2));
}
#[inline(always)]
pub unsafe fn _zxxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,0,2));
}
#[inline(always)]
pub unsafe fn _zxxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,0,2));
}
#[inline(always)]
pub unsafe fn _zxxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,0,2));
}

#[inline(always)]
pub unsafe fn _zxyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,0,2));
}
#[inline(always)]
pub unsafe fn _zxyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,0,2));
}
#[inline(always)]
pub unsafe fn _zxyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,0,2));
}
#[inline(always)]
pub unsafe fn _zxyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,0,2));
}

#[inline(always)]
pub unsafe fn _zxzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,0,2));
}
#[inline(always)]
pub unsafe fn _zxzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,0,2));
}
#[inline(always)]
pub unsafe fn _zxzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,0,2));
}
#[inline(always)]
pub unsafe fn _zxzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,0,2));
}

#[inline(always)]
pub unsafe fn _zxwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,0,2));
}
#[inline(always)]
pub unsafe fn _zxwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,0,2));
}
#[inline(always)]
pub unsafe fn _zxwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,0,2));
}
#[inline(always)]
pub unsafe fn _zxww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,0,2));
}


#[inline(always)]
pub unsafe fn _zyxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,1,2));
}
#[inline(always)]
pub unsafe fn _zyxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,1,2));
}
#[inline(always)]
pub unsafe fn _zyxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,1,2));
}
#[inline(always)]
pub unsafe fn _zyxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,1,2));
}

#[inline(always)]
pub unsafe fn _zyyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,1,2));
}
#[inline(always)]
pub unsafe fn _zyyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,1,2));
}
#[inline(always)]
pub unsafe fn _zyyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,1,2));
}
#[inline(always)]
pub unsafe fn _zyyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,1,2));
}

#[inline(always)]
pub unsafe fn _zyzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,1,2));
}
#[inline(always)]
pub unsafe fn _zyzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,1,2));
}
#[inline(always)]
pub unsafe fn _zyzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,1,2));
}
#[inline(always)]
pub unsafe fn _zyzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,1,2));
}

#[inline(always)]
pub unsafe fn _zywx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,1,2));
}
#[inline(always)]
pub unsafe fn _zywy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,1,2));
}
#[inline(always)]
pub unsafe fn _zywz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,1,2));
}
#[inline(always)]
pub unsafe fn _zyww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,1,2));
}



#[inline(always)]
pub unsafe fn _zzxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,2,2));
}
#[inline(always)]
pub unsafe fn _zzxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,2,2));
}
#[inline(always)]
pub unsafe fn _zzxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,2,2));
}
#[inline(always)]
pub unsafe fn _zzxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,2,2));
}

#[inline(always)]
pub unsafe fn _zzyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,2,2));
}
#[inline(always)]
pub unsafe fn _zzyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,2,2));
}
#[inline(always)]
pub unsafe fn _zzyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,2,2));
}
#[inline(always)]
pub unsafe fn _zzyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,2,2));
}

#[inline(always)]
pub unsafe fn _zzzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,2,2));
}
#[inline(always)]
pub unsafe fn _zzzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,2,2));
}
#[inline(always)]
pub unsafe fn _zzzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,2,2));
}
#[inline(always)]
pub unsafe fn _zzzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,2,2));
}

#[inline(always)]
pub unsafe fn _zzwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,2,2));
}
#[inline(always)]
pub unsafe fn _zzwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,2,2));
}
#[inline(always)]
pub unsafe fn _zzwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,2,2));
}
#[inline(always)]
pub unsafe fn _zzww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,2,2));
}

#[inline(always)]
pub unsafe fn _zwxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,3,2));
}
#[inline(always)]
pub unsafe fn _zwxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,3,2));
}
#[inline(always)]
pub unsafe fn _zwxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,3,2));
}
#[inline(always)]
pub unsafe fn _zwxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,3,2));
}

#[inline(always)]
pub unsafe fn _zwyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,3,2));
}
#[inline(always)]
pub unsafe fn _zwyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,3,2));
}
#[inline(always)]
pub unsafe fn _zwyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,3,2));
}
#[inline(always)]
pub unsafe fn _zwyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,3,2));
}

#[inline(always)]
pub unsafe fn _zwzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,3,2));
}
#[inline(always)]
pub unsafe fn _zwzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,3,2));
}
#[inline(always)]
pub unsafe fn _zwzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,3,2));
}
#[inline(always)]
pub unsafe fn _zwzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,3,2));
}

#[inline(always)]
pub unsafe fn _zwwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,3,2));
}
#[inline(always)]
pub unsafe fn _zwwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,3,2));
}
#[inline(always)]
pub unsafe fn _zwwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,3,2));
}
#[inline(always)]
pub unsafe fn _zwww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,3,2));
}
#[inline(always)]
pub unsafe fn _wxxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,0,3));
}
#[inline(always)]
pub unsafe fn _wxxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,0,3));
}
#[inline(always)]
pub unsafe fn _wxxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,0,3));
}
#[inline(always)]
pub unsafe fn _wxxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,0,3));
}

#[inline(always)]
pub unsafe fn _wxyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,0,3));
}
#[inline(always)]
pub unsafe fn _wxyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,0,3));
}
#[inline(always)]
pub unsafe fn _wxyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,0,3));
}
#[inline(always)]
pub unsafe fn _wxyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,0,3));
}

#[inline(always)]
pub unsafe fn _wxzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,0,3));
}
#[inline(always)]
pub unsafe fn _wxzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,0,3));
}
#[inline(always)]
pub unsafe fn _wxzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,0,3));
}
#[inline(always)]
pub unsafe fn _wxzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,0,3));
}

#[inline(always)]
pub unsafe fn _wxwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,0,3));
}
#[inline(always)]
pub unsafe fn _wxwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,0,3));
}
#[inline(always)]
pub unsafe fn _wxwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,0,3));
}
#[inline(always)]
pub unsafe fn _wxww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,0,3));
}


#[inline(always)]
pub unsafe fn _wyxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,1,3));
}
#[inline(always)]
pub unsafe fn _wyxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,1,3));
}
#[inline(always)]
pub unsafe fn _wyxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,1,3));
}
#[inline(always)]
pub unsafe fn _wyxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,1,3));
}

#[inline(always)]
pub unsafe fn _wyyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,1,3));
}
#[inline(always)]
pub unsafe fn _wyyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,1,3));
}
#[inline(always)]
pub unsafe fn _wyyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,1,3));
}
#[inline(always)]
pub unsafe fn _wyyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,1,3));
}

#[inline(always)]
pub unsafe fn _wyzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,1,3));
}
#[inline(always)]
pub unsafe fn _wyzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,1,3));
}
#[inline(always)]
pub unsafe fn _wyzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,1,3));
}
#[inline(always)]
pub unsafe fn _wyzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,1,3));
}

#[inline(always)]
pub unsafe fn _wywx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,1,3));
}
#[inline(always)]
pub unsafe fn _wywy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,1,3));
}
#[inline(always)]
pub unsafe fn _wywz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,1,3));
}
#[inline(always)]
pub unsafe fn _wyww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,1,3));
}



#[inline(always)]
pub unsafe fn _wzxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,2,3));
}
#[inline(always)]
pub unsafe fn _wzxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,2,3));
}
#[inline(always)]
pub unsafe fn _wzxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,2,3));
}
#[inline(always)]
pub unsafe fn _wzxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,2,3));
}

#[inline(always)]
pub unsafe fn _wzyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,2,3));
}
#[inline(always)]
pub unsafe fn _wzyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,2,3));
}
#[inline(always)]
pub unsafe fn _wzyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,2,3));
}
#[inline(always)]
pub unsafe fn _wzyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,2,3));
}

#[inline(always)]
pub unsafe fn _wzzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,2,3));
}
#[inline(always)]
pub unsafe fn _wzzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,2,3));
}
#[inline(always)]
pub unsafe fn _wzzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,2,3));
}
#[inline(always)]
pub unsafe fn _wzzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,2,3));
}

#[inline(always)]
pub unsafe fn _wzwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,2,3));
}
#[inline(always)]
pub unsafe fn _wzwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,2,3));
}
#[inline(always)]
pub unsafe fn _wzwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,2,3));
}
#[inline(always)]
pub unsafe fn _wzww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,2,3));
}

#[inline(always)]
pub unsafe fn _wwxx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,0,3,3));
}
#[inline(always)]
pub unsafe fn _wwxy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,0,3,3));
}
#[inline(always)]
pub unsafe fn _wwxz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,0,3,3));
}
#[inline(always)]
pub unsafe fn _wwxw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,0,3,3));
}

#[inline(always)]
pub unsafe fn _wwyx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,1,3,3));
}
#[inline(always)]
pub unsafe fn _wwyy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,1,3,3));
}
#[inline(always)]
pub unsafe fn _wwyz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,1,3,3));
}
#[inline(always)]
pub unsafe fn _wwyw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,1,3,3));
}

#[inline(always)]
pub unsafe fn _wwzx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,2,3,3));
}
#[inline(always)]
pub unsafe fn _wwzy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,2,3,3));
}
#[inline(always)]
pub unsafe fn _wwzz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,2,3,3));
}
#[inline(always)]
pub unsafe fn _wwzw(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,2,3,3));
}

#[inline(always)]
pub unsafe fn _wwwx(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(0,3,3,3));
}
#[inline(always)]
pub unsafe fn _wwwy(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(1,3,3,3));
}
#[inline(always)]
pub unsafe fn _wwwz(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(2,3,3,3));
}
#[inline(always)]
pub unsafe fn _wwww(m : __m128)->__m128{
  return _mm_shuffle_ps(m,m, _ico_shuffle(3,3,3,3));
}
