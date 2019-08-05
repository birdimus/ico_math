// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::raw::RawVector_i32;
use crate::structure::SIMDVector4;
use crate::vector2_bool::Vector2Bool;
use crate::vector3_bool::Vector3Bool;
use core::arch::x86_64::*;
use core::hash::Hash;
use core::hash::Hasher;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4Bool {
    pub data: __m128i,
}

impl Vector4Bool {
    #[inline(always)]
    pub fn new(x: bool, y: bool, z: bool, w: bool) -> Vector4Bool {
        let x_val: u32 = if x { 0xFFFFFFFF } else { 0 };
        let y_val: u32 = if y { 0xFFFFFFFF } else { 0 };
        let z_val: u32 = if z { 0xFFFFFFFF } else { 0 };
        let w_val: u32 = if w { 0xFFFFFFFF } else { 0 };
        unsafe {
            return Vector4Bool {
                data: _mm_set_epi32(
                    core::mem::transmute::<u32, i32>(w_val),
                    core::mem::transmute::<u32, i32>(z_val),
                    core::mem::transmute::<u32, i32>(y_val),
                    core::mem::transmute::<u32, i32>(x_val),
                ),
            };
        }
    }
    #[inline(always)]
    pub fn set(value: bool) -> Vector4Bool {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            return Vector4Bool {
                data: _mm_set1_epi32(core::mem::transmute::<u32, i32>(val)),
            };
        }
    }
    #[inline(always)]
    pub fn x(self) -> bool {
        unsafe {
            //1 2 4 8
            return (_mm_movemask_epi8(self.data) & 15) == 15;
        }
    }
    #[inline(always)]
    pub fn y(self) -> bool {
        unsafe {
            //16 32 64 128
            return (_mm_movemask_epi8(self.data) & 240) == 240;
        }
    }
    #[inline(always)]
    pub fn z(self) -> bool {
        unsafe {
            //256 512 1024 2048
            return (_mm_movemask_epi8(self.data) & 3840) == 3840;
        }
    }
    #[inline(always)]
    pub fn w(self) -> bool {
        unsafe {
            // 4096 8192 16384 32768
            return (_mm_movemask_epi8(self.data) & 61440) == 61440;
        }
    }

    #[inline(always)]
    pub fn set_x(&mut self, value: bool) {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 0);
        }
    }
    #[inline(always)]
    pub fn set_y(&mut self, value: bool) {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 1);
        }
    }
    #[inline(always)]
    pub fn set_z(&mut self, value: bool) {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 2);
        }
    }
    #[inline(always)]
    pub fn set_w(&mut self, value: bool) {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 3);
        }
    }

    #[inline(always)]
    pub fn and(self, v2: Vector4Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_and_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn or(self, v2: Vector4Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_or_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn andnot(self, v2: Vector4Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_andnot_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn xor(self, v2: Vector4Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_xor_si128(self.data, v2.data),
            };
        }
    }
    #[inline(always)]
    pub fn equal(self, v2: Vector4Bool) -> Vector4Bool {
        unsafe {
            Vector4Bool {
                data: _mm_cmpeq_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn not_equal(self, v2: Vector4Bool) -> Vector4Bool {
        return Vector4Bool::xor(self.equal(self), self.equal(v2));
    }

    #[inline(always)]
    pub fn all(self: Vector4Bool) -> bool {
        unsafe {
            return (_mm_movemask_epi8(self.data)) == 0x0000FFFF;
        }
    }
    #[inline(always)]
    pub fn any(self: Vector4Bool) -> bool {
        unsafe {
            return _mm_movemask_epi8(self.data) != 0;
        }
    }
}
impl From<bool> for Vector4Bool {
    #[inline(always)]
    fn from(v: bool) -> Vector4Bool {
        return Vector4Bool::set(v);
    }
}
impl From<Vector2Bool> for Vector4Bool {
    #[inline(always)]
    fn from(v: Vector2Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_castps_si128(_mm_movelh_ps(_mm_castsi128_ps(v.data), _mm_setzero_ps())),
            };
        }
    }
}
impl From<Vector3Bool> for Vector4Bool {
    #[inline(always)]
    fn from(v: Vector3Bool) -> Vector4Bool {
        unsafe {
            return Vector4Bool {
                data: _mm_srli_si128(_mm_slli_si128(v.data, 4), 4),
            };
        }
    }
}
impl PartialEq for Vector4Bool {
    #[inline(always)]
    fn eq(&self, other: &Vector4Bool) -> bool {
        return Vector4Bool::equal(*self, *other).all();
    }
}
impl Hash for Vector4Bool {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut dst = RawVector_i32 { data: [0; 4] };
        unsafe {
            let x: *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
            _mm_store_si128(x, self.data);
        }
        dst.data.hash(state);
    }
}
impl SIMDVector4 for Vector4Bool {
    #[inline(always)]
    fn data(self) -> __m128 {
        return unsafe { _mm_castsi128_ps(self.data) };
    }
    fn data_i(self) -> __m128i {
        return self.data;
    }
}
