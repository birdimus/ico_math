// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#[cfg(test)]
mod test {
    use crate::dual_quaternion::DualQuaternion;
    use crate::quaternion::Quaternion;
    use crate::vector3::Vector3;

    #[test]
    fn identity() {
        let a = DualQuaternion::identity();
        let b = DualQuaternion::new(Quaternion::identity(), Vector3::zero());
        assert_eq!(a.real.x(), b.real.x());
        assert_eq!(a.real.y(), b.real.y());
        assert_eq!(a.real.z(), b.real.z());
        assert_eq!(a.real.w(), b.real.w());
        assert_eq!(a.dual.x(), b.dual.x());
        assert_eq!(a.dual.y(), b.dual.y());
        assert_eq!(a.dual.z(), b.dual.z());
        assert_eq!(a.dual.w(), b.dual.w());
    }

}
