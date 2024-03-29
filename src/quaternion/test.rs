// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#[cfg(test)]
mod test {
    use crate::quaternion::Quaternion;
    use crate::quaternion::RotationOrder;
    use crate::vector3::Vector3;
    use crate::FloatVector;
    #[test]
    fn new() {
        let a = Quaternion::new(1.0, 2.0, 3.0, 4.0);

        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }

    #[test]
    fn identity() {
        let a = Quaternion::identity();

        assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), 0.0);
        assert_eq!(a.w(), 1.0);
    }

    #[test]
    fn set_x() {
        let mut a = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        a.set_x(-11.0);
        assert_eq!(a.x(), -11.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_y() {
        let mut a = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        a.set_y(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), -11.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_z() {
        let mut a = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        a.set_z(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), -11.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_w() {
        let mut a = Quaternion::new(1.0, 2.0, 3.0, 4.0);
        a.set_w(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), -11.0);
    }

    #[test]
    fn normalize() {
        let mut a = Quaternion::new(0.0, 0.0, 0.0, 0.0);
        a = a.normalize();
        assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), 0.0);
        assert_eq!(a.w(), 1.0);
    }
    #[test]
    fn angle_axis() {
        {
            let a = Vector3::new(0.0, 1.0, 0.0);
            let b = 90.0f32.to_radians();

            let c = Quaternion::angle_axis(b, a);
            assert_eq!(c.x(), 0.0);
            assert!(c.y().approx_equal(FloatVector::from(0.70710677)));
            assert_eq!(c.z(), 0.0);
            assert!(c.w().approx_equal(FloatVector::from(0.70710677)));
        }
        {
            let a = Vector3::new(0.0, 1.0, 0.0);
            let b = -90.0f32.to_radians();

            let c = Quaternion::angle_axis(b, a);
            assert_eq!(c.x(), 0.0);
            assert!(c.y().approx_equal(FloatVector::from(-0.70710677)));
            assert_eq!(c.z(), 0.0);
            assert!(c.w().approx_equal(FloatVector::from(0.70710677)));
        }
        {
            let a = Vector3::new(0.0, 1.0, 0.0);
            let b = -270.0f32.to_radians();

            let c = Quaternion::angle_axis(b, a);
            assert_eq!(c.x(), 0.0);
            assert!(c.y().approx_equal(FloatVector::from(-0.70710677)));
            assert_eq!(c.z(), 0.0);
            assert!(c.w().approx_equal(FloatVector::from(-0.70710677)));
        }
    }

    #[test]
    fn lerp() {
        let a = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();

        let to = Quaternion::angle_axis(b, a);
        let from = Quaternion::identity();

        for i in 0..1000 {
            let lerp_val = (i - 500) as f32 / 500.0f32;
            let result = Quaternion::lerp(from, to, lerp_val);
            let ideal = Quaternion::angle_axis(b * lerp_val, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.value() > 0.999 && proximity.value() < 1.001,
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::lerp(from, to, 0.0f32);
            let ideal = from;
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::lerp(from, to, 1.0f32);
            let ideal = to;
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} {}",
                proximity.value(),
                to.sqr_magnitude().value()
            );
        }
        {
            let result = Quaternion::lerp(from, to, -1.0f32);
            let ideal = Quaternion::angle_axis(b * -1.0f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
    }
    #[test]
    fn slerp() {
        let a = Vector3::new(0.0, 1.0, 0.0);
        let b = 170.0f32.to_radians();

        let to = Quaternion::angle_axis(b, a);
        let from = Quaternion::identity();

        //slerp can handle the big rotations.
        for i in 0..1000 {
            let lerp_val = (i - 500) as f32 / 200.0f32;
            let result = Quaternion::slerp(from, to, lerp_val);
            let ideal = Quaternion::angle_axis(b * lerp_val, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.value() > 0.999 && proximity.value() < 1.001,
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::slerp(from, to, 0.5f32);
            let ideal = Quaternion::angle_axis(b * 0.5f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, -0.5f32);
            let ideal = Quaternion::angle_axis(b * -0.5f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, 0.7f32);
            let ideal = Quaternion::angle_axis(b * 0.7f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, -0.7f32);
            let ideal = Quaternion::angle_axis(b * -0.7f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::slerp(from, to, 0.2f32);
            let ideal = Quaternion::angle_axis(b * 0.2f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::slerp(from, to, 0.0f32);
            let ideal = from;
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }

        {
            let result = Quaternion::slerp(from, to, 1.0f32);
            let ideal = to;
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, -1.0f32);
            let ideal = Quaternion::angle_axis(b * -1.0f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, 3.0f32);
            let ideal = Quaternion::angle_axis(b * 3.0f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
        {
            let result = Quaternion::slerp(from, to, -3.5f32);
            let ideal = Quaternion::angle_axis(b * -3.5f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(
                proximity.approx_equal(FloatVector::from(1.0)),
                "{} ",
                proximity.value()
            );
        }
    }

    #[test]
    fn euler_xyz() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::XYZ);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b2, right);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        c *= Quaternion::angle_axis(b, up);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        c *= Quaternion::angle_axis(b3, back);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }
    #[test]
    fn euler_xzy() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::XZY);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b2, right);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        c *= Quaternion::angle_axis(b3, back);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        c *= Quaternion::angle_axis(b, up);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }
    #[test]
    fn euler_yxz() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::YXZ);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b, up);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        c *= Quaternion::angle_axis(b2, right);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        c *= Quaternion::angle_axis(b3, back);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }
    #[test]
    fn euler_yzx() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::YZX);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b, up);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        c *= Quaternion::angle_axis(b3, back);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        c *= Quaternion::angle_axis(b2, right);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }
    #[test]
    fn euler_zxy() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::ZXY);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b3, back);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        c *= Quaternion::angle_axis(b2, right);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        c *= Quaternion::angle_axis(b, up);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }
    #[test]
    fn euler_zyx() {
        let a = Vector3::new(
            45.0f32.to_radians(),
            90.0f32.to_radians(),
            30.0f32.to_radians(),
        );
        //let a = Vector3::new(45.0f32.to_radians(), 0.0,0.0);
        let q = Quaternion::euler(a, RotationOrder::ZYX);

        let back = Vector3::new(0.0, 0.0, 1.0);
        let b3 = 30.0f32.to_radians();
        let mut c = Quaternion::angle_axis(b3, back);

        let up = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();
        c *= Quaternion::angle_axis(b, up);

        let right = Vector3::new(1.0, 0.0, 0.0);
        let b2 = 45.0f32.to_radians();
        c *= Quaternion::angle_axis(b2, right);

        let cmp = c.dot(q).value();
        assert!(cmp > 0.999 && cmp < 1.001, "{} ", cmp);
    }

}
