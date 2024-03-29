// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use core::arch::x86_64::*;

pub const ABSOLUTE_COMPARISON_EPSILON: f32 = 1.00000011920928955078125 - 1.0;
pub const RELATIVE_COMPARISON_EPSILON: f32 = ABSOLUTE_COMPARISON_EPSILON * 4.0; //0.000000238418579e-7
                                                                                // pub const NORMALIZATION_EPSILON: f32 = 1.0e-10;
pub const SIGN_BIT: f32 = -0.0;
// pub const EPSILON_AT_ONE: f32 = 0.00000012;
pub const INV_TWO_PI: f32 = 0.159154943091895335768883763372514362034459645740456448747;
pub const INV_PI: f32 = 1.0 / core::f32::consts::PI;
pub const INV_TWO_PI_64: f64 = 0.5 / core::f64::consts::PI;
pub const INV_360_64: f64 = 1.0 / 360.0;

pub const TWO_PI: f32 = 6.283185307179586476925286766559005768394338798750211641949;
pub const HALF_PI: f32 = 1.57079632679;
pub const PI: f32 = 3.14159265359;
pub const INV_360: f32 = 1.0 / 360.0;
/// A replacement for the _MM_SHUFFLE macro.
#[inline(always)]
pub const fn _ico_shuffle(z: i32, y: i32, x: i32, w: i32) -> i32 {
    (((z) << 6) | ((y) << 4) | ((x) << 2) | (w))
}

// Generators for common vectors to avoid set.
// https://www.agner.org/optimize/optimizing_assembly.pdf

// #[inline(always)]
// pub unsafe fn _ico_one_epi32() -> __m128i {
//     let mut a = _mm_setzero_si128();
//     a = _mm_cmpeq_epi32(a, a);
//     return _mm_srli_epi32(a, 31);
// }

#[inline(always)]
pub unsafe fn _ico_half_ps() -> __m128 {
    return _mm_set1_ps(0.5f32);
    // let mut a = _mm_setzero_si128();
    // a = _mm_cmpeq_epi32(a, a);
    // a = _mm_slli_epi32(a, 26);
    // return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
pub unsafe fn _ico_one_ps() -> __m128 {
    return _mm_set1_ps(1.0f32);
    // let mut a = _mm_setzero_si128();
    // a = _mm_cmpeq_epi32(a, a);
    // a = _mm_slli_epi32(a, 25);
    // return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
}
#[inline(always)]
pub unsafe fn _ico_two_ps() -> __m128 {
    return _mm_set1_ps(2.0f32);
    // let mut a = _mm_setzero_si128();
    // a = _mm_cmpeq_epi32(a, a);
    // a = _mm_slli_epi32(a, 31);
    // return _mm_castsi128_ps(_mm_srli_epi32(a, 1));
}
// #[inline(always)]
// pub unsafe fn _ico_nearesttwo_ps() -> __m128 {
//     let mut a = _mm_setzero_si128();
//     a = _mm_cmpeq_epi32(a, a);
//     return _mm_castsi128_ps(_mm_srli_epi32(a, 2));
// }
// #[inline(always)]
// pub unsafe fn _ico_negtwo_ps() -> __m128 {
//     let mut a = _mm_setzero_si128();
//     a = _mm_cmpeq_epi32(a, a);
//     return _mm_castsi128_ps(_mm_slli_epi32(a, 30));
// }
#[inline(always)]
pub unsafe fn _ico_signbit_ps() -> __m128 {
    return _mm_set1_ps(SIGN_BIT);
    // let mut a = _mm_setzero_si128();
    // a = _mm_cmpeq_epi32(a, a);
    // return _mm_castsi128_ps(_mm_slli_epi32(a, 31));
}
#[inline(always)]
pub unsafe fn _ico_abs_ps(a: __m128) -> __m128 {
    return _mm_andnot_ps(_ico_signbit_ps(), a);
}
#[inline(always)]
pub unsafe fn _ico_copysign_ps(a: __m128, b: __m128) -> __m128 {
    // When we xor the sign bit, both positive or negative = 0.  If they differ, a's sign bit is set to true.
    // But then we mask the sign bit off with the absolute value.
    let tmp = _ico_abs_ps(_mm_xor_ps(a, b));
    // When we xor them back together, if b is negative, it's value comes through. A always 0.
    // 0 ^ 1 = 1, 0 ^ 0 = 0

    return _mm_xor_ps(tmp, b);
}

#[inline(always)]
pub unsafe fn _ico_select_ps(a: __m128, b: __m128, mask: __m128) -> __m128 {
    return _mm_blendv_ps(a, b, mask);
    // (((b ^ a) & mask)^a)

    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
    //return _mm_or_ps(_mm_and_ps(a, mask), _mm_andnot_ps(mask, b));
    //return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}
#[inline(always)]
pub unsafe fn _ico_select_pd(a: __m128d, b: __m128d, mask: __m128d) -> __m128d {
    return _mm_blendv_pd(a, b, mask);
    // (((b ^ a) & mask)^a)

    // https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
    //return _mm_or_pd(_mm_and_pd(b, mask), _mm_andnot_pd(mask, a));
    //return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}

#[inline(always)]
pub unsafe fn _ico_select_si128(a: __m128i, b: __m128i, mask: __m128i) -> __m128i {
    // (((b ^ a) & mask)^a)
    return _mm_or_si128(_mm_and_si128(b, mask), _mm_andnot_si128(mask, a));
    //return _mm_xor_si128( a, _mm_and_si128( mask, _mm_xor_si128 ( b, a ) ) );
    //_mm_or_si128(_mm_and_si128(a, cond), _mm_andnot_si128(cond, b))
}

///If the float is entirely integer or greater - we return the float
/// Otherwise we cast to int
#[inline(always)]
pub unsafe fn _ico_truncate_ps(a: __m128) -> __m128 {
    // Largest float < integer only is  < 8388608
    // NaN and int-only floats fail this test.
    let mask = _mm_cmplt_ps(_ico_abs_ps(a), _mm_castsi128_ps(_mm_set1_epi32(0x4b000000)));
    // is the float entirely integer?
    let trunc = _mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    //if Abs(a) > 8388608
    // select a if greater then 8388608.0f, otherwise select the result of FuncT
    // NaN will fail and get passed through.
    return _ico_select_ps(a, trunc, mask);
}

#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_floor_ps(a: __m128) -> __m128 {
    //convert to int using truncate
    let t = _ico_truncate_ps(a); //_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
                                 // if the truncated value is greater than the old value (negative, non integral) subtract 1.
    return _mm_sub_ps(t, _mm_and_ps(_mm_cmpgt_ps(t, a), _ico_one_ps()));
}

#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_ceil_ps(a: __m128) -> __m128 {
    //convert to int using truncate
    let t = _ico_truncate_ps(a); //_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
                                 // if the truncated value is less than the old value (positive) add 1.
    return _mm_add_ps(t, _mm_and_ps(_mm_cmplt_ps(t, a), _ico_one_ps()));
}

///Round away from zero, not toward zero to match cmath.
/// Returns unsigned 0.
//http://dss.stephanierct.com/DevBlog/?p=8
#[inline(always)]
/// Returns unsigned 0.
pub unsafe fn _ico_round_ps(a: __m128) -> __m128 {
    //convert to int using truncate
    let trunc = _ico_truncate_ps(a); //_mm_cvtepi32_ps(_mm_cvttps_epi32(a));
    let remainder = _mm_sub_ps(a, trunc);
    // If we want to round toward zero, use _ico_nearesttwo_ps()
    let rmd2 = _mm_mul_ps(remainder, _ico_two_ps());
    let rmd2i = _mm_cvttps_epi32(rmd2); // after being truncated of course

    return _mm_add_ps(trunc, _mm_cvtepi32_ps(rmd2i));
}

#[inline(always)]
pub unsafe fn _ico_cross_ps(lhs: __m128, rhs: __m128) -> __m128 {
    // x  <-  a.y*b.z - a.z*b.y
    // y  <-  a.z*b.x - a.x*b.z
    // z  <-  a.x*b.y - a.y*b.x
    // We can save a shuffle by grouping it in this wacky order:
    let mut tmp0 = _mm_mul_ps(lhs, _zxyw(rhs));
    //__m128 tmp1 =  _mm_mul_ps(_mm_shuffle_ps(lhs, lhs, _MM_SHUFFLE(3, 1, 0, 2)), rhs);
    //tmp0 = _mm_sub_ps(tmp1, tmp0);
    tmp0 = _mm_fmsub_ps(_zxyw(lhs), rhs, tmp0);

    return _zxyw(tmp0); // _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(3, 1, 0, 2));
}

#[inline(always)]
pub unsafe fn _ico_dp4_ps(v0: __m128, v1: __m128) -> __m128 {
    let mut tmp0 = _mm_mul_ps(v0, v1);
    let mut tmp1 = _yxwz(tmp0); // _mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(2, 3, 0, 1)); //yxwz
    tmp0 = _mm_add_ps(tmp0, tmp1); //xy, xy, zw, zw
    tmp1 = _wzyx(tmp0); //_mm_shuffle_ps(tmp0, tmp0, _ico_shuffle(0, 1, 2, 3)); //wxyz
    return _mm_add_ps(tmp0, tmp1);
}

#[inline(always)]
pub unsafe fn _ico_quat_mul(lhs: __m128, rhs: __m128) -> __m128 {
    let x = _mm_mul_ps(_wzyx(rhs), _xxxx(lhs));
    let xw = _mm_fmsubadd_ps(_wwww(lhs), rhs, x);
    let y = _mm_xor_ps(_yyyy(lhs), _mm_set_ps(-0.0, -0.0, 0.0, 0.0));
    let xyw = _mm_fmadd_ps(_zwxy(rhs), y, xw);
    let z = _mm_xor_ps(_zzzz(lhs), _mm_set_ps(-0.0, 0.0, 0.0, -0.0));
    return _mm_fmadd_ps(_yxwz(rhs), z, xyw);
}

#[inline(always)]
/// Bounce between 0 and 1
pub unsafe fn _ico_ping_pong(vec: __m128) -> __m128 {
    // 2 * abs(x + 0.5 - floor(x + 0.5) - 0.5)
    //  2 * abs(x - floor(x + 0.5))
    let shifted = _mm_add_ps(vec, _ico_half_ps());
    let internal = _ico_abs_ps(_mm_sub_ps(vec, _mm_floor_ps(shifted)));
    return _mm_add_ps(internal, internal); // multiply by 2
}

/// 1/4 of a cosine curve
/// Max absolute error is < 0.00002
/// Both 0 and 1 should return exact values (1, and 0, respectively)
#[inline(always)]
pub unsafe fn _ico_approx_cos01(vec: __m128) -> __m128 {
    let scalars = _mm_set_ps(-0.05104357, 0.300697, -0.01847417, -1.23117923841858); //-1.231179
    let mut result = _mm_fmadd_ps(vec, _wwww(scalars), _zzzz(scalars));
    result = _mm_fmadd_ps(vec, result, _yyyy(scalars));
    result = _mm_fmadd_ps(vec, result, _xxxx(scalars));
    result = _mm_mul_ps(result, vec);
    result = _mm_fmadd_ps(vec, result, _ico_one_ps());
    return result;
}

// #[inline(always)]
// unsafe fn _ico_do_cos_ps(scaled: __m128) -> __m128 {
//     //reduce range first. - doesn't actually seem to make a difference in precision - probably limited by the mul with inv pi.
//     let ranged = _mm_sub_ps(scaled,_mm_floor_ps(scaled) );
//     let ping_pong = _ico_abs_ps(_mm_sub_ps(ranged, _mm_set1_ps(0.5)));

//     //this contains the sign
//     let sign_driver = _mm_sub_ps(ping_pong, _mm_set1_ps(0.25));
//     //convert the sign ping pong to a 0-1 driver.
//     let driver = _mm_fnmadd_ps(_ico_abs_ps(sign_driver), _mm_set1_ps(4.0), _mm_set1_ps(1.0));

//     return _ico_copysign_ps(_ico_approx_cos01(driver), sign_driver);
// }

#[inline(always)]
unsafe fn _ico_do_cos_pd(scaled: __m256d) -> __m128 {
    let ranged = _mm256_cvtpd_ps(_mm256_sub_pd(scaled, _mm256_floor_pd(scaled)));
    let ping_pong = _ico_abs_ps(_mm_sub_ps(ranged, _mm_set1_ps(0.5)));

    //this contains the sign
    let sign_driver = _mm_sub_ps(ping_pong, _mm_set1_ps(0.25));
    //convert the sign ping pong to a 0-1 driver.
    let driver = _mm_fnmadd_ps(_ico_abs_ps(sign_driver), _mm_set1_ps(4.0), _mm_set1_ps(1.0));

    return _ico_copysign_ps(_ico_approx_cos01(driver), sign_driver);
}

/// Converts to double precision for range reduction, or else we are limited by the precision in PI
#[inline(always)]
pub unsafe fn _ico_cos_ps(vec: __m128) -> __m128 {
    let dbl = _mm256_cvtps_pd(vec);
    let scaled = _mm256_mul_pd(dbl, _mm256_set1_pd(INV_TWO_PI_64));
    return _ico_do_cos_pd(scaled);
}
#[inline(always)]
pub unsafe fn _ico_cos_deg_ps(vec: __m128) -> __m128 {
    let dbl = _mm256_cvtps_pd(vec);
    let scaled = _mm256_mul_pd(dbl, _mm256_set1_pd(INV_360_64));
    return _ico_do_cos_pd(scaled);
}
#[inline(always)]
pub unsafe fn _ico_sin_ps(vec: __m128) -> __m128 {
    //TODO: we could range reduce first, before shifting - but I don't think it matters much
    // SIN range reduction would look like: abs(0.5 - abs(floor(vec) - vec + 0.25)) - 0.25
    let dbl = _mm256_cvtps_pd(vec);
    let scaled = _mm256_fmsub_pd(dbl, _mm256_set1_pd(INV_TWO_PI_64), _mm256_set1_pd(0.25));
    // let scaled = _mm_fmsub_ps(vec, _mm_set1_ps(INV_TWO_PI), _mm_set1_ps(0.25));
    return _ico_do_cos_pd(scaled);
}
#[inline(always)]
pub unsafe fn _ico_sin_deg_ps(vec: __m128) -> __m128 {
    let dbl = _mm256_cvtps_pd(vec);
    let scaled = _mm256_fmsub_pd(dbl, _mm256_set1_pd(INV_360_64), _mm256_set1_pd(0.25));
    // let scaled = _mm_fmsub_ps(vec, _mm_set1_ps(INV_TWO_PI), _mm_set1_ps(0.25));
    return _ico_do_cos_pd(scaled);
}

#[inline(always)]
pub unsafe fn _ico_tan_ps(vec: __m128) -> __m128 {
    let dbl = _mm256_cvtps_pd(vec);
    let scaled = _mm256_mul_pd(dbl, _mm256_set1_pd(INV_TWO_PI_64));
    let cos = _ico_do_cos_pd(scaled);
    let sin = _ico_do_cos_pd(_mm256_sub_pd(scaled, _mm256_set1_pd(0.25)));

    return _mm_div_ps(sin, cos);
}

/// Approximation Valid between -1 and 1
/// 7th order approximation
#[inline(always)]
pub unsafe fn _ico_atan01_ps(value: __m128, offset: __m128) -> __m128 {
    let scalars = _mm_set_ps(-0.03864493, 0.1459276, -0.3210924, 0.9992079);

    //Ax7 + bx5 + cx3 +dx
    //x( x2( x2 (Ax2 + b) + c) +d)
    let x_sqr = _mm_mul_ps(value, value);
    let mut result = _mm_fmadd_ps(x_sqr, _wwww(scalars), _zzzz(scalars));
    result = _mm_fmadd_ps(x_sqr, result, _yyyy(scalars));
    result = _mm_fmadd_ps(x_sqr, result, _xxxx(scalars));
    return _mm_fmadd_ps(value, result, offset);
}

// loosely based on https://www.dsprelated.com/showarticle/1052.php
#[inline(always)]
pub unsafe fn _ico_atan2_ps(y: __m128, x: __m128) -> __m128 {
    let abs_x = _ico_abs_ps(x);
    let abs_y = _ico_abs_ps(y);
    let x_div_y = _mm_div_ps(x, y);
    let y_div_x = _mm_div_ps(y, x);

    //if x >= y
    let x_greater_y = _mm_cmpge_ps(abs_x, abs_y);

    let x_negative = _mm_cmplt_ps(x, _mm_setzero_ps());
    let y_negative = _mm_cmplt_ps(y, _mm_setzero_ps());

    //5 cases
    // if x is positive, 0, otherwise either positive or negative pi offset
    let case_01_offset = _mm_and_ps(
        x_negative,
        _ico_select_ps(
            _mm_set1_ps(core::f32::consts::PI),
            _mm_set1_ps(-core::f32::consts::PI),
            y_negative,
        ),
    );
    // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
    let case_23_offset = _ico_select_ps(
        _mm_set1_ps(-0.5 * core::f32::consts::PI),
        _mm_set1_ps(0.5 * core::f32::consts::PI),
        y_negative,
    );

    let offset = _ico_select_ps(case_23_offset, case_01_offset, x_greater_y);
    let atan_val = _ico_select_ps(x_div_y, y_div_x, x_greater_y);
    let sign_flip = _mm_andnot_ps(x_greater_y, _ico_signbit_ps());

    let atan = _mm_xor_ps(sign_flip, _ico_atan01_ps(atan_val, offset));

    return atan;
}

#[inline(always)]
/// max error 0.000001 radians
/// 4th order approximation
/// acos(x)/sqrt(1-x) goes from pi/2 at 0 to sqrt(2) at 1 ==  0 to -0.15658276442
pub unsafe fn _ico_acos_ps(vec: __m128) -> __m128 {
    let abs_vec = _ico_abs_ps(vec);
    let scalars = _mm_set_ps(-0.2142619, 0.08528893, -0.03673658, 0.009126827);
    let mut result = _mm_fmadd_ps(_xxxx(scalars), abs_vec, _yyyy(scalars));
    result = _mm_fmadd_ps(result, abs_vec, _zzzz(scalars));
    result = _mm_fmadd_ps(result, abs_vec, _wwww(scalars));
    result = _mm_fmadd_ps(result, abs_vec, _mm_set1_ps(HALF_PI));
    let sqrt_vec = _mm_sqrt_ps(_mm_sub_ps(_ico_one_ps(), abs_vec));
    result = _mm_mul_ps(result, sqrt_vec);

    let inv_mask = _mm_cmplt_ps(vec, _mm_setzero_ps());
    result = _mm_xor_ps(result, _mm_and_ps(_ico_signbit_ps(), inv_mask));
    result = _mm_add_ps(result, _mm_and_ps(inv_mask, _mm_set1_ps(PI)));
    return result;
}
pub unsafe fn _ico_asin_ps(vec: __m128) -> __m128 {
    return _mm_sub_ps(_mm_set1_ps(HALF_PI), _ico_acos_ps(vec));
}

#[inline(always)]
pub unsafe fn _xxxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 0, 0));
}

#[inline(always)]
pub unsafe fn _xyxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 0));
}

#[inline(always)]
pub unsafe fn _xyyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 0));
}

#[inline(always)]
pub unsafe fn _xyzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 0));
}

#[inline(always)]
pub unsafe fn _xywx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xywy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xywz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 1, 0));
}

#[inline(always)]
pub unsafe fn _xzxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 2, 0));
}

#[inline(always)]
pub unsafe fn _xwxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 3, 0));
}

#[inline(always)]
pub unsafe fn _yxxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 0, 1));
}

#[inline(always)]
pub unsafe fn _yyxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 1));
}

#[inline(always)]
pub unsafe fn _yyyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 1));
}

#[inline(always)]
pub unsafe fn _yyzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 1));
}

#[inline(always)]
pub unsafe fn _yywx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yywy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yywz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 1, 1));
}

#[inline(always)]
pub unsafe fn _yzxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 2, 1));
}

#[inline(always)]
pub unsafe fn _ywxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 3, 1));
}

#[inline(always)]
pub unsafe fn _zxxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 0, 2));
}

#[inline(always)]
pub unsafe fn _zyxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 2));
}

#[inline(always)]
pub unsafe fn _zyyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 2));
}

#[inline(always)]
pub unsafe fn _zyzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 2));
}

#[inline(always)]
pub unsafe fn _zywx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zywy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zywz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 1, 2));
}

#[inline(always)]
pub unsafe fn _zzxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 2, 2));
}

#[inline(always)]
pub unsafe fn _zwxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _wxxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 0, 3));
}

#[inline(always)]
pub unsafe fn _wyxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 1, 3));
}

#[inline(always)]
pub unsafe fn _wyyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 1, 3));
}

#[inline(always)]
pub unsafe fn _wyzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 1, 3));
}

#[inline(always)]
pub unsafe fn _wywx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wywy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wywz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 1, 3));
}

#[inline(always)]
pub unsafe fn _wzxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 2, 3));
}

#[inline(always)]
pub unsafe fn _wwxx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 0, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwyx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 1, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwzx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzw(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 2, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwwx(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(0, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwwy(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(1, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwwz(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(2, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwww(m: __m128) -> __m128 {
    return _mm_shuffle_ps(m, m, _ico_shuffle(3, 3, 3, 3));
}

/////////. INTEGER SHUFFLES

#[inline(always)]
pub unsafe fn _xxxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 0, 0));
}

#[inline(always)]
pub unsafe fn _xxwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 0, 0));
}
#[inline(always)]
pub unsafe fn _xxww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 0, 0));
}

#[inline(always)]
pub unsafe fn _xyxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 1, 0));
}

#[inline(always)]
pub unsafe fn _xyyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 1, 0));
}

#[inline(always)]
pub unsafe fn _xyzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 1, 0));
}

#[inline(always)]
pub unsafe fn _xywx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xywy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xywz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 1, 0));
}
#[inline(always)]
pub unsafe fn _xyww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 1, 0));
}

#[inline(always)]
pub unsafe fn _xzxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 2, 0));
}

#[inline(always)]
pub unsafe fn _xzwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 2, 0));
}
#[inline(always)]
pub unsafe fn _xzww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 2, 0));
}

#[inline(always)]
pub unsafe fn _xwxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 3, 0));
}

#[inline(always)]
pub unsafe fn _xwwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 3, 0));
}
#[inline(always)]
pub unsafe fn _xwww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 3, 0));
}

#[inline(always)]
pub unsafe fn _yxxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 0, 1));
}

#[inline(always)]
pub unsafe fn _yxwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 0, 1));
}
#[inline(always)]
pub unsafe fn _yxww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 0, 1));
}

#[inline(always)]
pub unsafe fn _yyxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 1, 1));
}

#[inline(always)]
pub unsafe fn _yyyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 1, 1));
}

#[inline(always)]
pub unsafe fn _yyzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 1, 1));
}

#[inline(always)]
pub unsafe fn _yywx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yywy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yywz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 1, 1));
}
#[inline(always)]
pub unsafe fn _yyww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 1, 1));
}

#[inline(always)]
pub unsafe fn _yzxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 2, 1));
}

#[inline(always)]
pub unsafe fn _yzwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 2, 1));
}
#[inline(always)]
pub unsafe fn _yzww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 2, 1));
}

#[inline(always)]
pub unsafe fn _ywxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 3, 1));
}

#[inline(always)]
pub unsafe fn _ywwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 3, 1));
}
#[inline(always)]
pub unsafe fn _ywww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 3, 1));
}

#[inline(always)]
pub unsafe fn _zxxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 0, 2));
}

#[inline(always)]
pub unsafe fn _zxwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 0, 2));
}
#[inline(always)]
pub unsafe fn _zxww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 0, 2));
}

#[inline(always)]
pub unsafe fn _zyxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 1, 2));
}

#[inline(always)]
pub unsafe fn _zyyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 1, 2));
}

#[inline(always)]
pub unsafe fn _zyzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 1, 2));
}

#[inline(always)]
pub unsafe fn _zywx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zywy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zywz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 1, 2));
}
#[inline(always)]
pub unsafe fn _zyww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 1, 2));
}

#[inline(always)]
pub unsafe fn _zzxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 2, 2));
}

#[inline(always)]
pub unsafe fn _zzwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 2, 2));
}
#[inline(always)]
pub unsafe fn _zzww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 2, 2));
}

#[inline(always)]
pub unsafe fn _zwxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 3, 2));
}

#[inline(always)]
pub unsafe fn _zwwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _zwww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 3, 2));
}
#[inline(always)]
pub unsafe fn _wxxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 0, 3));
}

#[inline(always)]
pub unsafe fn _wxwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 0, 3));
}
#[inline(always)]
pub unsafe fn _wxww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 0, 3));
}

#[inline(always)]
pub unsafe fn _wyxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 1, 3));
}

#[inline(always)]
pub unsafe fn _wyyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 1, 3));
}

#[inline(always)]
pub unsafe fn _wyzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 1, 3));
}

#[inline(always)]
pub unsafe fn _wywx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wywy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wywz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 1, 3));
}
#[inline(always)]
pub unsafe fn _wyww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 1, 3));
}

#[inline(always)]
pub unsafe fn _wzxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 2, 3));
}

#[inline(always)]
pub unsafe fn _wzwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 2, 3));
}
#[inline(always)]
pub unsafe fn _wzww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 2, 3));
}

#[inline(always)]
pub unsafe fn _wwxx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 0, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwxw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 0, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwyx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 1, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwyw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 1, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwzx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 2, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwzw_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 2, 3, 3));
}

#[inline(always)]
pub unsafe fn _wwwx_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(0, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwwy_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(1, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwwz_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(2, 3, 3, 3));
}
#[inline(always)]
pub unsafe fn _wwww_i(m: __m128i) -> __m128i {
    return _mm_shuffle_epi32(m, _ico_shuffle(3, 3, 3, 3));
}
