use crate::raw::RawVector_i32;
use crate::structure::SIMDVector2;
use crate::vector3_bool::Vector3Bool;
use crate::vector4_bool::Vector4Bool;
use core::arch::x86_64::*;
use core::hash::Hash;
use core::hash::Hasher;
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2Bool {
    pub data: __m128i,
}

impl Vector2Bool {
    #[inline(always)]
    pub fn new(x: bool, y: bool) -> Vector2Bool {
        let x_val: u32 = if x { 0xFFFFFFFF } else { 0 };
        let y_val: u32 = if y { 0xFFFFFFFF } else { 0 };

        unsafe {
            return Vector2Bool {
                data: _mm_set_epi32(
                    0,
                    0,
                    core::mem::transmute::<u32, i32>(y_val),
                    core::mem::transmute::<u32, i32>(x_val),
                ),
            };
        }
    }
    /// Set all values of the vector to the same f32 value.
    #[inline(always)]
    pub fn set(value: bool) -> Vector2Bool {
        let val: u32 = if value { 0xFFFFFFFF } else { 0 };
        unsafe {
            return Vector2Bool {
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
    pub fn and(self, v2: Vector2Bool) -> Vector2Bool {
        unsafe {
            return Vector2Bool {
                data: _mm_and_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn or(self, v2: Vector2Bool) -> Vector2Bool {
        unsafe {
            return Vector2Bool {
                data: _mm_or_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn andnot(self, v2: Vector2Bool) -> Vector2Bool {
        unsafe {
            return Vector2Bool {
                data: _mm_andnot_si128(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn xor(self, v2: Vector2Bool) -> Vector2Bool {
        unsafe {
            return Vector2Bool {
                data: _mm_xor_si128(self.data, v2.data),
            };
        }
    }
    #[inline(always)]
    pub fn equal(self, v2: Vector2Bool) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_cmpeq_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn not_equal(self, v2: Vector2Bool) -> Vector2Bool {
        return Vector2Bool::xor(self.equal(self), self.equal(v2));
    }

    #[inline(always)]
    pub fn all(self: Vector2Bool) -> bool {
        unsafe {
            return (_mm_movemask_epi8(self.data) & 255) == 255;
        }
    }
    #[inline(always)]
    pub fn any(self: Vector2Bool) -> bool {
        unsafe {
            return (_mm_movemask_epi8(self.data) & 255) != 0;
        }
    }
}
impl From<bool> for Vector2Bool {
    #[inline(always)]
    fn from(v: bool) -> Vector2Bool {
        return Vector2Bool::set(v);
    }
}
impl From<Vector3Bool> for Vector2Bool {
    #[inline(always)]
    fn from(v: Vector3Bool) -> Vector2Bool {
        Vector2Bool { data: v.data }
    }
}
impl From<Vector4Bool> for Vector2Bool {
    #[inline(always)]
    fn from(v: Vector4Bool) -> Vector2Bool {
        Vector2Bool { data: v.data }
    }
}
impl PartialEq for Vector2Bool {
    #[inline(always)]
    fn eq(&self, other: &Vector2Bool) -> bool {
        return Vector2Bool::equal(*self, *other).all();
    }
}
impl Hash for Vector2Bool {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut dst = RawVector_i32 { data: [0; 4] };
        unsafe {
            let x: *mut __m128i = &mut (dst.data[0]) as *mut i32 as *mut __m128i;
            _mm_store_si128(x, self.data);
        }
        dst.data[2] = 0;
        dst.data[3] = 0;
        dst.data.hash(state);
    }
}
impl SIMDVector2 for Vector2Bool {
    #[inline(always)]
    fn data(self) -> __m128 {
        return unsafe { _mm_castsi128_ps(self.data) };
    }
    fn data_i(self) -> __m128i {
        return unsafe { self.data };
    }
}
