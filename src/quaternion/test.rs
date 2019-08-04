use super::*;

#[cfg(test)]
mod test {
    use crate::float_vector::FloatVector;
    use crate::quaternion::Quaternion;
    use crate::vector3::Vector3;
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
    fn angle_axis() {
        let a = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();

        let c = Quaternion::angle_axis(b, a);
        assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 0.70710677);
        assert_eq!(c.z(), 0.0);
        assert_eq!(c.w(), 0.70710677);
    }

    #[test]
    fn lerp() {
        let a = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();

        let to = Quaternion::angle_axis(b, a);
        let from = Quaternion::identity();

        {
            let result = Quaternion::lerp(from, to, 0.5f32);
            let ideal = Quaternion::angle_axis(b * 0.5f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999);
        }

        {
            let result = Quaternion::lerp(from, to, 0.7f32);
            let ideal = Quaternion::angle_axis(b * 0.7f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999, "{} ", proximity.value());
        }

        {
            let result = Quaternion::lerp(from, to, 0.2f32);
            let ideal = Quaternion::angle_axis(b * 0.2f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999, "{} ", proximity.value());
        }

        {
            let result = Quaternion::lerp(from, to, 0.0f32);
            let ideal = from;
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() == 1.0, "{} ", proximity.value());
        }

        {
            let result = Quaternion::lerp(from, to, 1.0f32);
            let ideal = to;
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() == 1.0, "{} ", proximity.value());
        }
    }
    #[test]
    fn slerp() {
        let a = Vector3::new(0.0, 1.0, 0.0);
        let b = 90.0f32.to_radians();

        let to = Quaternion::angle_axis(b, a);
        let from = Quaternion::identity();

        {
            let result = Quaternion::slerp(from, to, 0.5f32);
            let ideal = Quaternion::angle_axis(b * 0.5f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999);
        }

        {
            let result = Quaternion::slerp(from, to, 0.7f32);
            let ideal = Quaternion::angle_axis(b * 0.7f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999, "{} ", proximity.value());
        }

        {
            let result = Quaternion::slerp(from, to, 0.2f32);
            let ideal = Quaternion::angle_axis(b * 0.2f32, a);
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() > 0.999, "{} ", proximity.value());
        }

        {
            let result = Quaternion::slerp(from, to, 0.0f32);
            let ideal = from;
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() == 1.0, "{} ", proximity.value());
        }

        {
            let result = Quaternion::slerp(from, to, 1.0f32);
            let ideal = to;
            let proximity = Quaternion::dot(result, ideal);
            assert!(proximity.value() == 1.0, "{} ", proximity.value());
        }
    }

}
