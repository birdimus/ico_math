use crate::quaternion::Quaternion;
use crate::sse_extensions::*;
use crate::vector3::Vector3;
use crate::vector4::Vector4;
use core::arch::x86_64::*;

#[derive(Copy, Clone)]
#[repr(C, align(16))]
pub struct Matrix4x4 {
    pub m: [__m128; 4],
}

impl Matrix4x4 {
    #[inline(always)]
    pub fn new(
        m00: f32,
        m01: f32,
        m02: f32,
        m03: f32,
        m10: f32,
        m11: f32,
        m12: f32,
        m13: f32,
        m20: f32,
        m21: f32,
        m22: f32,
        m23: f32,
        m30: f32,
        m31: f32,
        m32: f32,
        m33: f32,
    ) -> Matrix4x4 {
        unsafe {
            return Matrix4x4 {
                m: [
                    _mm_set_ps(m03, m02, m01, m00),
                    _mm_set_ps(m13, m12, m11, m10),
                    _mm_set_ps(m23, m22, m21, m20),
                    _mm_set_ps(m33, m32, m31, m30),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn zero() -> Matrix4x4 {
        unsafe {
            return Matrix4x4 {
                m: [
                    _mm_setzero_ps(),
                    _mm_setzero_ps(),
                    _mm_setzero_ps(),
                    _mm_setzero_ps(),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn identity() -> Matrix4x4 {
        return Matrix4x4::new(
            1.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 1.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 1.0f32,
            0.0f32, 0.0f32, 0.0f32, 0.0f32, 1.0f32,
        );
    }

    #[inline(always)]
    pub fn from_columns(c0: Vector4, c1: Vector4, c2: Vector4, c3: Vector4) -> Matrix4x4 {
        return Matrix4x4 {
            m: [c0.data, c1.data, c2.data, c3.data],
        };
    }

    #[inline(always)]
    pub fn from_rows(r0: Vector4, r1: Vector4, r2: Vector4, r3: Vector4) -> Matrix4x4 {
        unsafe {
            let tmp0 = _mm_unpacklo_ps(r0.data, r1.data);
            let tmp2 = _mm_unpacklo_ps(r2.data, r3.data);
            let tmp1 = _mm_unpackhi_ps(r0.data, r1.data);
            let tmp3 = _mm_unpackhi_ps(r2.data, r3.data);
            return Matrix4x4 {
                m: [
                    _mm_movelh_ps(tmp0, tmp2),
                    _mm_movehl_ps(tmp2, tmp0),
                    _mm_movelh_ps(tmp1, tmp3),
                    _mm_movehl_ps(tmp3, tmp1),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn transpose(&self) -> Matrix4x4 {
        unsafe {
            let tmp0 = _mm_unpacklo_ps(self.m[0], self.m[1]);
            let tmp2 = _mm_unpacklo_ps(self.m[2], self.m[3]);
            let tmp1 = _mm_unpackhi_ps(self.m[0], self.m[1]);
            let tmp3 = _mm_unpackhi_ps(self.m[2], self.m[3]);
            return Matrix4x4 {
                m: [
                    _mm_movelh_ps(tmp0, tmp2),
                    _mm_movehl_ps(tmp2, tmp0),
                    _mm_movelh_ps(tmp1, tmp3),
                    _mm_movehl_ps(tmp3, tmp1),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn add(a: Matrix4x4, b: Matrix4x4) -> Matrix4x4 {
        unsafe {
            return Matrix4x4 {
                m: [
                    _mm_add_ps(a.m[0], b.m[0]),
                    _mm_add_ps(a.m[1], b.m[1]),
                    _mm_add_ps(a.m[2], b.m[2]),
                    _mm_add_ps(a.m[3], b.m[3]),
                ],
            };
        }
    }
    #[inline(always)]
    pub fn sub(a: Matrix4x4, b: Matrix4x4) -> Matrix4x4 {
        unsafe {
            return Matrix4x4 {
                m: [
                    _mm_sub_ps(a.m[0], b.m[0]),
                    _mm_sub_ps(a.m[1], b.m[1]),
                    _mm_sub_ps(a.m[2], b.m[2]),
                    _mm_sub_ps(a.m[3], b.m[3]),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn scale(mat: Matrix4x4, value: f32) -> Matrix4x4 {
        unsafe {
            let tmp0 = _mm_set1_ps(value);
            return Matrix4x4 {
                m: [
                    _mm_mul_ps(mat.m[0], tmp0),
                    _mm_mul_ps(mat.m[1], tmp0),
                    _mm_mul_ps(mat.m[2], tmp0),
                    _mm_mul_ps(mat.m[3], tmp0),
                ],
            };
        }
    }

    #[inline(always)]
    pub fn mul(a: Matrix4x4, b: Matrix4x4) -> Matrix4x4 {
        unsafe {
            let mut col0 = _mm_mul_ps(a.m[0], _xxxx(b.m[0]));
            let mut col1 = _mm_mul_ps(a.m[0], _xxxx(b.m[1]));
            let mut col2 = _mm_mul_ps(a.m[0], _xxxx(b.m[2]));
            let mut col3 = _mm_mul_ps(a.m[0], _xxxx(b.m[3]));

            col0 = _mm_fmadd_ps(a.m[1], _yyyy(b.m[0]), col0);
            col1 = _mm_fmadd_ps(a.m[1], _yyyy(b.m[1]), col1);
            col2 = _mm_fmadd_ps(a.m[1], _yyyy(b.m[2]), col2);
            col3 = _mm_fmadd_ps(a.m[1], _yyyy(b.m[3]), col3);

            col0 = _mm_fmadd_ps(a.m[2], _zzzz(b.m[0]), col0);
            col1 = _mm_fmadd_ps(a.m[2], _zzzz(b.m[1]), col1);
            col2 = _mm_fmadd_ps(a.m[2], _zzzz(b.m[2]), col2);
            col3 = _mm_fmadd_ps(a.m[2], _zzzz(b.m[3]), col3);
            return Matrix4x4 {
                m: [
                    _mm_fmadd_ps(a.m[3], _wwww(b.m[0]), col0),
                    _mm_fmadd_ps(a.m[3], _wwww(b.m[1]), col1),
                    _mm_fmadd_ps(a.m[3], _wwww(b.m[2]), col2),
                    _mm_fmadd_ps(a.m[3], _wwww(b.m[3]), col3),
                ],
            };
        }
    }
    #[inline(always)]
    pub fn mul_vec(mat: Matrix4x4, v1: Vector4) -> Vector4 {
        //x,y,z,w = column[0] * vec.x
        //x,y,z,w = column[1] * vec.y
        //x,y,z,w = column[2] * vec.z
        //x,y,z,w = column[3] * vec.w
        unsafe {
            let mut result = _mm_mul_ps(mat.m[0], _xxxx(v1.data));
            result = _mm_fmadd_ps(mat.m[1], _yyyy(v1.data), result);
            result = _mm_fmadd_ps(mat.m[2], _zzzz(v1.data), result);
            return Vector4 {
                data: _mm_fmadd_ps(mat.m[3], _wwww(v1.data), result),
            };
        }
    }

    #[inline(always)]
    pub fn transform_point(mat: Matrix4x4, v1: Vector3) -> Vector4 {
        //x,y,z = column[0] * vec.x
        //x,y,z = column[1] * vec.y
        //x,y,z = column[2] * vec.z
        // + column[3]
        unsafe {
            let mut result = _mm_fmadd_ps(mat.m[0], _xxxx(v1.data), mat.m[3]);
            result = _mm_fmadd_ps(mat.m[1], _yyyy(v1.data), result);
            result = _mm_fmadd_ps(mat.m[2], _zzzz(v1.data), result);
            return Vector4 { data: result };
        }
    }

    #[inline(always)]
    pub fn transform_vector(mat: Matrix4x4, v1: Vector3) -> Vector4 {
        //x,y,z = column[0] * vec.x
        //x,y,z = column[1] * vec.y
        //x,y,z = column[2] * vec.z
        unsafe {
            let mut result = _mm_mul_ps(mat.m[0], _xxxx(v1.data));
            result = _mm_fmadd_ps(mat.m[1], _yyyy(v1.data), result);
            return Vector4 {
                data: _mm_fmadd_ps(mat.m[2], _zzzz(v1.data), result),
            };
        }
    }

    ///  From Essential Mathematics for games 2nd ed.
    #[inline(always)]
    pub fn from_rotation(q: Quaternion) -> Matrix4x4 {
        unsafe {
            let wzyx = _wzyx(q.data);
            let zwxy = _zwxy(q.data);
            let yxwz = _yxwz(q.data);
            let xyzw = (q.data);

            let m1 = Matrix4x4 {
                m: [
                    // 0 0 - -
                    _mm_xor_ps(_mm_set_ps(SIGN_BIT, SIGN_BIT, 0.0, 0.0), wzyx),
                    // - 0 0 -
                    _mm_xor_ps(_mm_set_ps(SIGN_BIT, 0.0, 0.0, SIGN_BIT), zwxy),
                    // 0,-,0,-
                    _mm_xor_ps(_mm_set_ps(SIGN_BIT, 0.0, SIGN_BIT, 0.0), yxwz),
                    // 0 0 0 0
                    xyzw,
                ],
            };
            let m2 = Matrix4x4 {
                m: [
                    // 0 0 - 0
                    _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, 0.0, 0.0), wzyx),
                    // - 0 0 0
                    _mm_xor_ps(_mm_set_ps(0.0, 0.0, 0.0, SIGN_BIT), zwxy),
                    // 0,-,0,0
                    _mm_xor_ps(_mm_set_ps(0.0, 0.0, SIGN_BIT, 0.0), yxwz),
                    // - - - 0
                    _mm_xor_ps(_mm_set_ps(0.0, SIGN_BIT, SIGN_BIT, SIGN_BIT), xyzw),
                ],
            };

            return Matrix4x4::mul(m2, m1);
        }
    }
    #[inline(always)]
    pub fn from_trs(translation: Vector3, rotation: Quaternion, uniform_scale: f32) -> Matrix4x4 {
        let mut mat = Matrix4x4::from_rotation(rotation);

        unsafe {
            let scale = _mm_set1_ps(uniform_scale);
            mat.m[0] = _mm_mul_ps(mat.m[0], scale);
            mat.m[1] = _mm_mul_ps(mat.m[1], scale);
            mat.m[2] = _mm_mul_ps(mat.m[2], scale);
            //set the position value
            mat.m[3] = _yzwx(_mm_move_ss(_wxyz(translation.data), _mm_set_ss(1.0f32)));
        }
        return mat;
    }

    #[inline(always)]
    pub fn equals(a: Matrix4x4, b: Matrix4x4) -> bool {
        unsafe {
            return _mm_movemask_ps(_mm_cmpeq_ps(a.m[0], b.m[0])) == 15
                && _mm_movemask_ps(_mm_cmpeq_ps(a.m[1], b.m[1])) == 15
                && _mm_movemask_ps(_mm_cmpeq_ps(a.m[2], b.m[2])) == 15
                && _mm_movemask_ps(_mm_cmpeq_ps(a.m[3], b.m[3])) == 15;
        }
    }
}
impl core::ops::Mul<Matrix4x4> for Matrix4x4 {
    type Output = Matrix4x4;
    #[inline]
    fn mul(self, _rhs: Matrix4x4) -> Matrix4x4 {
        return Matrix4x4::mul(self, _rhs);
    }
}
impl PartialEq for Matrix4x4 {
    fn eq(&self, other: &Matrix4x4) -> bool {
        return Matrix4x4::equals(*self, *other);
    }
}
