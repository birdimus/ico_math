use crate::quaternion::Quaternion;
use crate::sse_extensions::*;
use crate::vector3::Vector3;
use core::arch::x86_64::*;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct DualQuaternion {
    pub real: __m128,
    pub dual: __m128,
}

impl DualQuaternion {
    #[inline(always)]
    pub fn new(rotation: Quaternion, translation: Vector3) -> DualQuaternion {
        unsafe {
            let vec = _mm_mul_ps(translation.data, _mm_set_ps(0f32, 0.5f32, 0.5f32, 0.5f32));
            let dual = _ico_quat_mul(vec, rotation.data);
            return DualQuaternion {
                real: rotation.data,
                dual: dual,
            };
        }
    }

    #[inline(always)]
    pub fn rotation(&self) -> Quaternion {
        return Quaternion { data: self.real };
    }
    #[inline(always)]
    pub fn translation(&self) -> Vector3 {
        let inv_rot = self.rotation().inverse();
        unsafe {
            let p = _ico_quat_mul(self.dual, inv_rot.data);
            return Vector3 {
                data: _mm_add_ps(p, p),
            }; //*2
        }
    }

    #[inline(always)]
    pub fn inverse(&self) -> DualQuaternion {
        let real = Quaternion { data: self.real }.inverse();
        let dual = Quaternion { data: self.dual }.inverse();
        return DualQuaternion {
            real: real.data,
            dual: dual.data,
        };
    }

    #[inline(always)]
    pub fn reverse(&self) -> DualQuaternion {
        let real = Quaternion { data: self.real }.reverse();
        let dual = Quaternion { data: self.dual }.reverse();
        return DualQuaternion {
            real: real.data,
            dual: dual.data,
        };
    }
    #[inline(always)]
    pub fn mul(lhs: DualQuaternion, rhs: DualQuaternion) -> DualQuaternion {
        unsafe {
            return DualQuaternion {
                real: _ico_quat_mul(lhs.real, rhs.real),
                dual: _mm_add_ps(
                    _ico_quat_mul(lhs.real, rhs.dual),
                    _ico_quat_mul(lhs.dual, rhs.real),
                ),
            };
        }
    }

    #[inline(always)]
    pub fn renormalize(&self) -> DualQuaternion {
        unsafe {
            let sqr_length = _ico_dp4_ps(self.real, self.real);
            let factor = _mm_sqrt_ps(sqr_length);

            return DualQuaternion {
                real: _mm_div_ps(self.real, factor),
                dual: _mm_div_ps(self.dual, factor),
            };
        }
    }
}

impl PartialEq for DualQuaternion {
    #[inline(always)]
    fn eq(&self, other: &DualQuaternion) -> bool {
        let real_equal = Quaternion::equal(
            Quaternion { data: self.real },
            Quaternion { data: other.real },
        );
        let dual_equal = Quaternion::equal(
            Quaternion { data: self.dual },
            Quaternion { data: other.dual },
        );
        let both_equal = real_equal.and(dual_equal);
        return both_equal.all();
    }
}
