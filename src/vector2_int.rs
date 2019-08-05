use crate::int_vector::IntVector;
use crate::raw::RawVector_i32;
use crate::sse_extensions::*;
use crate::structure::SIMDVector2;
use crate::vector2::Vector2;
use crate::vector2_bool::Vector2Bool;
use crate::vector3_int::Vector3Int;
use crate::vector4_int::Vector4Int;
use core::arch::x86_64::*;
use core::hash::Hash;
use core::hash::Hasher;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2Int {
    pub data: __m128i,
}

impl Vector2Int {
    /// Returns a new Vector2
    #[inline(always)]
    pub fn new(x: i32, y: i32) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_set_epi32(0, 0, y, x),
            }
        }
    }
    #[inline(always)]
    pub fn set<T: Into<IntVector>>(value: T) -> Vector2Int {
        return Vector2Int {
            data: value.into().data,
        };
    }
    #[inline(always)]
    pub fn zero() -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_setzero_si128(),
            }
        }
    }
    #[inline(always)]
    pub fn x(self) -> IntVector {
        return IntVector {
            data: self.xxxx().data,
        };
    }

    #[inline(always)]
    pub fn y(self) -> IntVector {
        return IntVector {
            data: self.yyyy().data,
        };
    }

    #[inline(always)]
    pub fn set_x<T: Into<i32>>(&mut self, value: T) {
        unsafe {
            self.data = _mm_insert_epi32(self.data, value.into(), 0);
        }
    }

    #[inline(always)]
    pub fn set_y<T: Into<i32>>(&mut self, value: T) {
        unsafe {
            self.data = _mm_insert_epi32(self.data, value.into(), 1);
        }
    }

    #[inline(always)]
    pub fn max(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_max_epi32(self.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn min(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_min_epi32(self.data, v2.data),
            }
        }
    }
    /// Choose component wise between A and B based on the mask.  False = A, True = B.
    #[inline(always)]
    pub fn select(self, v2: Vector2Int, mask: Vector2Bool) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _ico_select_si128(self.data, v2.data, mask.data),
            }
        }
    }

    #[inline(always)]
    pub fn add(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_add_epi32(self.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn sub(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_sub_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn mul(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_mullo_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn abs(self) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _mm_abs_epi32(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn copysign(self, v2: Vector2Int) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _mm_sign_epi32(self.data, v2.data),
            };
        }
    }

    #[inline(always)]
    pub fn and<T: SIMDVector2>(self, v2: T) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_and_si128(self.data, v2.data_i()),
            }
        }
    }

    #[inline(always)]
    pub fn or<T: SIMDVector2>(self, v2: T) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_or_si128(self.data, v2.data_i()),
            }
        }
    }

    #[inline(always)]
    pub fn andnot<T: SIMDVector2>(self, v2: T) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_andnot_si128(self.data, v2.data_i()),
            }
        }
    }

    #[inline(always)]
    pub fn xor<T: SIMDVector2>(self, v2: T) -> Vector2Int {
        unsafe {
            Vector2Int {
                data: _mm_xor_si128(self.data, v2.data_i()),
            }
        }
    }

    #[inline(always)]
    pub fn equal(self, v2: Vector2Int) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_cmpeq_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn not_equal(self, v2: Vector2Int) -> Vector2Bool {
        return Vector2Bool::xor(self.equal(self), self.equal(v2));
    }

    #[inline(always)]
    pub fn greater_equal(self, v2: Vector2Int) -> Vector2Bool {
        return Vector2Bool::or(self.equal(v2), self.greater(v2));
    }
    #[inline(always)]
    pub fn greater(self, v2: Vector2Int) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_cmpgt_epi32(self.data, v2.data),
            }
        }
    }
    #[inline(always)]
    pub fn less_equal(self, v2: Vector2Int) -> Vector2Bool {
        return Vector2Bool::or(self.equal(v2), self.less(v2));
    }

    #[inline(always)]
    pub fn less(self, v2: Vector2Int) -> Vector2Bool {
        unsafe {
            Vector2Bool {
                data: _mm_cmplt_epi32(self.data, v2.data),
            }
        }
    }

    #[inline(always)]
    pub fn xxxx(self) -> Vector4Int {
        unsafe {
            return Vector4Int {
                data: _xxxx_i(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yyyy(self) -> Vector4Int {
        unsafe {
            return Vector4Int {
                data: _yyyy_i(self.data),
            };
        }
    }

    #[inline(always)]
    pub fn xx(self) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _xxzw_i(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn xy(self) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _xyzw_i(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yx(self) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _yxzw_i(self.data),
            };
        }
    }
    #[inline(always)]
    pub fn yy(self) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _yyzw_i(self.data),
            };
        }
    }
}

impl From<i32> for Vector2Int {
    #[inline(always)]
    fn from(v: i32) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _mm_set1_epi32(v),
            };
        }
    }
}
impl From<IntVector> for Vector2Int {
    #[inline(always)]
    fn from(val: IntVector) -> Vector2Int {
        return Vector2Int { data: val.data };
    }
}
impl From<Vector3Int> for Vector2Int {
    #[inline(always)]
    fn from(v: Vector3Int) -> Vector2Int {
        Vector2Int { data: v.data }
    }
}
impl From<Vector4Int> for Vector2Int {
    #[inline(always)]
    fn from(v: Vector4Int) -> Vector2Int {
        Vector2Int { data: v.data }
    }
}
impl From<Vector2> for Vector2Int {
    #[inline(always)]
    fn from(v: Vector2) -> Vector2Int {
        unsafe {
            return Vector2Int {
                data: _mm_cvttps_epi32(v.data),
            };
        }
    }
}

impl core::ops::Add for Vector2Int {
    type Output = Vector2Int;
    #[inline(always)]
    fn add(self, _rhs: Vector2Int) -> Vector2Int {
        return Vector2Int::add(self, _rhs);
    }
}
impl core::ops::AddAssign for Vector2Int {
    #[inline(always)]
    fn add_assign(&mut self, other: Vector2Int) {
        *self = Vector2Int::add(*self, other);
    }
}
impl core::ops::Sub for Vector2Int {
    type Output = Vector2Int;
    #[inline(always)]
    fn sub(self, _rhs: Vector2Int) -> Vector2Int {
        return Vector2Int::sub(self, _rhs);
    }
}
impl core::ops::SubAssign for Vector2Int {
    #[inline(always)]
    fn sub_assign(&mut self, other: Vector2Int) {
        *self = Vector2Int::sub(*self, other);
    }
}
impl core::ops::Neg for Vector2Int {
    type Output = Vector2Int;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        unsafe {
            return Vector2Int {
                data: _mm_sub_epi32(_mm_set1_epi32(0), self.data),
            };
        }
    }
}

impl<T: Into<IntVector>> core::ops::Mul<T> for Vector2Int {
    type Output = Vector2Int;
    #[inline(always)]
    fn mul(self, _rhs: T) -> Vector2Int {
        return Vector2Int::mul(self, Vector2Int::from(_rhs.into()));
    }
}
impl<T: Into<IntVector>> core::ops::MulAssign<T> for Vector2Int {
    #[inline(always)]
    fn mul_assign(&mut self, _rhs: T) {
        *self = Vector2Int::mul(*self, Vector2Int::from(_rhs.into()));
    }
}
impl core::ops::Mul<Vector2Int> for IntVector {
    type Output = Vector2Int;
    #[inline(always)]
    fn mul(self: IntVector, _rhs: Vector2Int) -> Vector2Int {
        return Vector2Int::mul(_rhs, Vector2Int::from(self));
    }
}

impl PartialEq for Vector2Int {
    #[inline(always)]
    fn eq(&self, other: &Vector2Int) -> bool {
        return Vector2Int::equal(*self, *other).all();
    }
}

impl Hash for Vector2Int {
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
impl SIMDVector2 for Vector2Int {
    #[inline(always)]
    fn data(self) -> __m128 {
        return unsafe { _mm_castsi128_ps(self.data) };
    }
    fn data_i(self) -> __m128i {
        return self.data;
    }
}
