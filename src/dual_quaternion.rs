// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::quaternion::Quaternion;
use crate::vector3::Vector3;
use crate::vector4::Vector4;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct DualQuaternion {
    pub real: Quaternion,
    pub dual: Quaternion,
}

impl DualQuaternion {
    #[inline(always)]
    pub fn identity() -> DualQuaternion {
        return DualQuaternion {
            real: Quaternion::identity(),
            dual: Quaternion::from(Vector4::zero()),
        };
    }

    #[inline(always)]
    pub fn new(rotation: Quaternion, translation: Vector3) -> DualQuaternion {
        let vec = Vector4::from(translation) * 0.5f32;
        let dual = Quaternion::from(vec) * rotation;
        return DualQuaternion {
            real: rotation,
            dual: dual,
        };
    }

    #[inline(always)]
    pub fn rotation(&self) -> Quaternion {
        return self.real;
    }
    #[inline(always)]
    pub fn translation(&self) -> Vector3 {
        let inv_rot = self.real.inverse();
        let p = self.dual * inv_rot;
        let vec = Vector4::from(p);
        return Vector3::from(vec + vec);
    }

    #[inline(always)]
    pub fn inverse(&self) -> DualQuaternion {
        return DualQuaternion {
            real: self.real.inverse(),
            dual: self.dual.inverse(),
        };
    }

    #[inline(always)]
    pub fn reverse(&self) -> DualQuaternion {
        return DualQuaternion {
            real: self.real.reverse(),
            dual: self.dual.reverse(),
        };
    }

    #[inline(always)]
    pub fn mul(lhs: DualQuaternion, rhs: DualQuaternion) -> DualQuaternion {
        let new_dual = Vector4::from(lhs.real * rhs.dual) + Vector4::from(lhs.dual * rhs.real);
        return DualQuaternion {
            real: lhs.real * rhs.real,
            dual: Quaternion::from(new_dual),
        };
    }

    #[inline(always)]
    pub fn renormalize(&self) -> DualQuaternion {
        let length = self.real.magnitude();

        return DualQuaternion {
            real: Quaternion::from(Vector4::from(self.real) / length),
            dual: Quaternion::from(Vector4::from(self.dual) / length),
        };
    }
}

impl PartialEq for DualQuaternion {
    #[inline(always)]
    fn eq(&self, other: &DualQuaternion) -> bool {
        let real_equal = self.real.equal(other.real);
        let dual_equal = self.dual.equal(other.dual);
        let both_equal = real_equal.and(dual_equal);
        return both_equal.all();
    }
}
#[cfg(test)]
mod test;
