// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::raw::RawVector_i32;
use crate::structure::SIMDVector3;
use crate::vector2_bool::Vector2Bool;
use crate::vector4_bool::Vector4Bool;
use core::arch::x86_64::*;
use core::hash::Hash;
use core::hash::Hasher;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector3Bool {
    pub data: __m128i,
}

impl Vector3Bool {
    #[inline(always)]
    pub fn new(x: bool, y: bool, z: bool) -> Vector3Bool {
        let x_val: u32 = if x { 0xFFFFFFFF } else { 0 };
        let y_val: u32 = if y { 0xFFFFFFFF } else { 0 };
        let z_val: u32 = if z { 0xFFFFFFFF } else { 0 };

        unsafe {
            return Vector3Bool {
                data: _mm_set_epi32(
                    0,
                    core::mem::transmute::<u32, i32>(z_val),
                    core::mem::transmute::<u32, i32>(y_val),
                    core::mem::transmute::<u32, i32>(x_val),
                ),
            };
        }
    }
    #[inline(always)]
    pub fn set(value: bool) -> Vector3Bool {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            return Vector3Bool {
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
    pub fn and(self, v2: Vector3Bool) -> Vector3Bool {
        unsafe {
            return Vector3Bool {
                data: _mm_and_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn or(self, v2: Vector3Bool) -> Vector3Bool {
        unsafe {
            return Vector3Bool {
                data: _mm_or_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn andnot(self, v2: Vector3Bool) -> Vector3Bool {
        unsafe {
            return Vector3Bool {
                data: _mm_andnot_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn xor(self, v2: Vector3Bool) -> Vector3Bool {
        unsafe {
            return Vector3Bool {
                data: _mm_xor_si128(self.data, v2.data),
            };
        }
    }
    #[inline(always)]
    pub fn equal(self, v2: Vector3Bool) -> Vector3Bool {
        unsafe {
            Vector3Bool {
                data: _mm_cmpeq_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn not_equal(self, v2: Vector3Bool) -> Vector3Bool {
        return Vector3Bool::xor(self.equal(self), self.equal(v2));
    }

    #[inline(always)]
    pub fn all(self: Vector3Bool) -> bool {
        unsafe {
            return (_mm_movemask_epi8(self.data) & 4095) == 4095;
        }
    }
    #[inline(always)]
    pub fn any(self: Vector3Bool) -> bool {
        unsafe {
            return (_mm_movemask_epi8(self.data) & 4095) != 0;
        }
    }
}

impl From<bool> for Vector3Bool {
    #[inline(always)]
    fn from(v: bool) -> Vector3Bool {
        return Vector3Bool::set(v);
    }
}
impl From<Vector2Bool> for Vector3Bool {
    #[inline(always)]
    fn from(v: Vector2Bool) -> Vector3Bool {
        unsafe {
            return Vector3Bool {
                data: _mm_castps_si128(_mm_movelh_ps(_mm_castsi128_ps(v.data), _mm_setzero_ps())),
            };
        }
    }
}
impl From<Vector4Bool> for Vector3Bool {
    #[inline(always)]
    fn from(v: Vector4Bool) -> Vector3Bool {
        Vector3Bool { data: v.data }
    }
}
impl PartialEq for Vector3Bool {
    #[inline(always)]
    fn eq(&self, other: &Vector3Bool) -> bool {
        return Vector3Bool::equal(*self, *other).all();
    }
}
impl Hash for Vector3Bool {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut dst = RawVector_i32 { data: [0; 4] };
        unsafe {
            let x: *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
            _mm_store_si128(x, self.data);
        }
        dst.data[3] = 0;
        dst.data.hash(state);
    }
}
impl SIMDVector3 for Vector3Bool {
    #[inline(always)]
    fn data(self) -> __m128 {
        return unsafe { _mm_castsi128_ps(self.data) };
    }
    #[inline(always)]
    fn data_i(self) -> __m128i {
        return self.data;
    }
}
