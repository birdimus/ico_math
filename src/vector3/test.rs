#[cfg(test)]
mod test {
    use crate::float_vector::FloatVector;
    use crate::vector3::Vector3;
    use crate::vector3_bool::Vector3Bool;
    use crate::vector4::Vector4;
    use crate::raw::RawVector_f32;
    #[test]
    fn new() {
        let a = Vector3::new(1.0, 2.0, 3.0);

        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
    }
    #[test]
    fn set() {
        let a = Vector3::set(-5.0);

        assert_eq!(a.x(), -5.0);
        assert_eq!(a.y(), -5.0);
        assert_eq!(a.z(), -5.0);
    }
    #[test]
    fn zero() {
        let a = Vector3::zero();

        assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), 0.0);
    }
    #[test]
    fn store() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let mut b = RawVector_f32{data:[0.0;4]};
        a.store(&mut b);
        assert_eq!(b.data[0], 1.0);
        assert_eq!(b.data[1], 2.0);
        assert_eq!(b.data[2], 3.0);
        b.data[0] += 1.0;
        b.data[3] += 1.0;
        let c = Vector3::load(&b);
        let e = c + a;
        e.store(&mut b);
        assert_eq!(b.data[0], 3.0);
        assert_eq!(b.data[1], 4.0);
        assert_eq!(b.data[2], 6.0);
    }
    #[test]
    fn set_x() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a.set_x(-11.0);
        assert_eq!(a.x(), -11.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
    }
    #[test]
    fn set_y() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a.set_y(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), -11.0);
        assert_eq!(a.z(), 3.0);
    }
    #[test]
    fn set_z() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a.set_z(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), -11.0);
    }
    #[test]
    fn add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);

        let c = a + b;
        assert_eq!(c.x(), 5.0);
        assert_eq!(c.y(), 6.0);
        assert_eq!(c.z(), 9.0);
    }
    #[test]
    fn add_assign() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a += Vector3::set(3.0);

        assert_eq!(a.x(), 4.0);
        assert_eq!(a.y(), 5.0);
        assert_eq!(a.z(), 6.0);
    }
    #[test]
    fn sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);

        let c = a - b;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
    }
    #[test]
    fn sub_assign() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a -= Vector3::new(4.0, 4.0, 6.0);

        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -2.0);
        assert_eq!(a.z(), -3.0);
    }
    #[test]
    fn neg() {
        let a = Vector3::new(1.0, -2.0, 3.0);
        let c = -a;

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 2.0);
        assert_eq!(c.z(), -3.0);
    }
    #[test]
    fn component_mul() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);

        let c = Vector3::mul(a, b);
        assert_eq!(c.x(), 4.0);
        assert_eq!(c.y(), 8.0);
        assert_eq!(c.z(), 18.0);
    }
    #[test]
    fn component_div() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);

        let c = Vector3::div(a, b);
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.5);
    }
    #[test]
    fn mul_add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);
        let c = Vector3::new(2.0, -5.0, 10.0);
        let d = Vector3::mul_add(a, b, c);

        assert_eq!(d.x(), 6.0);
        assert_eq!(d.y(), 3.0);
        assert_eq!(d.z(), 28.0);
    }
    #[test]
    fn mul_sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);
        let c = Vector3::new(2.0, -5.0, 10.0);
        let d = Vector3::mul_sub(a, b, c);

        assert_eq!(d.x(), 2.0);
        assert_eq!(d.y(), 13.0);
        assert_eq!(d.z(), 8.0);
    }
    #[test]
    fn neg_mul_add() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);
        let c = Vector3::new(2.0, -5.0, 10.0);
        let d = Vector3::neg_mul_add(a, b, c);

        assert_eq!(d.x(), -2.0);
        assert_eq!(d.y(), -13.0);
        assert_eq!(d.z(), -8.0);
    }
    #[test]
    fn neg_mul_sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 4.0, 6.0);
        let c = Vector3::new(2.0, -5.0, 10.0);
        let d = Vector3::neg_mul_sub(a, b, c);

        assert_eq!(d.x(), -6.0);
        assert_eq!(d.y(), -3.0);
        assert_eq!(d.z(), -28.0);
    }
    #[test]
    fn mul() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = FloatVector::new(-3.0);
        {
            let c = a * b;
            assert_eq!(c.x(), -3.0);
            assert_eq!(c.y(), -6.0);
            assert_eq!(c.z(), -9.0);
        }
        {
            let c = b * a;
            assert_eq!(c.x(), -3.0);
            assert_eq!(c.y(), -6.0);
            assert_eq!(c.z(), -9.0);
        }
    }
    #[test]
    fn mul_assign() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a *= -3.0;

        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -6.0);
        assert_eq!(a.z(), -9.0);
    }
    #[test]
    fn div() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = 4.0;

        let c = a / b;
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.75);
    }
    #[test]
    fn div_inv() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = FloatVector::new(1.0);

        let c = b / a;
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 1.0 / 3.0);
    }
    #[test]
    fn div_assign() {
        let mut a = Vector3::new(1.0, 2.0, 3.0);
        a /= 4.0;

        assert_eq!(a.x(), 0.25);
        assert_eq!(a.y(), 0.5);
        assert_eq!(a.z(), 0.75);
    }
    #[test]
    fn and() {
        let a = Vector3::new(-1.0, -2.0, 3.0);
        let b = Vector3::new(1.0, -2.0, -3.0);

        let c = Vector3::and(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
    }
    #[test]
    fn andnot() {
        let a = Vector3::new(0.0, 0.0, -0.0);
        let b = Vector3::new(1.0, -2.0, -3.0);

        let c = Vector3::andnot(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
    }
    #[test]
    fn or() {
        let a = Vector3::new(-1.0, -2.0, 3.0);
        let b = Vector3::new(1.0, -2.0, -3.0);

        let c = Vector3::or(a, b);
        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
    }
    #[test]
    fn xor() {
        let a = Vector3::new(-1.0, -2.0, 3.0);
        let b = Vector3::new(-0.0, 0.0, -0.0);

        let c = Vector3::xor(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
    }
    #[test]
    fn all() {
        {
            let a = Vector3::new(1.0, 2.0, 3.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), true);
        }
        {
            let a = Vector3::new(0.0, 2.0, 3.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3::new(1.0, 0.0, 3.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3::new(1.0, 2.0, 0.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3::new(0.0, 0.0, 0.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        unsafe {
            use core::arch::x86_64::*;
            let a = Vector3 {
                data: _mm_set_ps(0.0, 1.0, 2.0, 3.0),
            };
            let b = Vector3 {
                data: _mm_set_ps(99.0, 1.0, 2.0, 3.0),
            };
            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::all(c), true);
        }
    }
    #[test]
    fn any() {
        {
            let a = Vector3::new(1.0, 2.0, 3.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3::new(0.0, 0.0, 3.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3::new(1.0, 0.0, 0.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3::new(0.0, 2.0, 0.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3::new(0.0, 0.0, 0.0);
            let b = Vector3::new(1.0, 2.0, 3.0);

            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), false);
        }
        unsafe {
            use core::arch::x86_64::*;
            let a = Vector3 {
                data: _mm_set_ps(0.0, 1.0, 2.0, 3.0),
            };
            let b = Vector3 {
                data: _mm_set_ps(0.0, 99.0, 99.0, 99.0),
            };
            let c = Vector3::equal(a, b);
            assert_eq!(Vector3Bool::any(c), false);
        }
    }
    #[test]
    fn equal() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn not_equal() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::not_equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), -1000.0);
    }
    #[test]
    fn greater_equal() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::greater_equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn greater() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::greater(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn less_equal() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::less_equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
    }
    #[test]
    fn less() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(1.0, 1.0, 4.0);

        let mask = Vector3::less(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
    }
     #[test]
    fn approx_equal() {
        {
        let a = Vector3::new(1.0, 1.0, 1.0);
        let b = Vector3::new(1.0 + 4.0*core::f32::EPSILON, 1.0 - 4.0*core::f32::EPSILON, 1.0 + 5.0*core::f32::EPSILON);

        let mask = Vector3::approx_equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
        }
        {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(core::f32::EPSILON, -core::f32::EPSILON, 2.0 * core::f32::EPSILON);

        let mask = Vector3::approx_equal(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
        }
    }
    #[test]
    fn definitely_greater() {
        {
        let a = Vector3::new(1.0, 1.0, 1.0);
        let b = Vector3::new(1.0 + 4.0*core::f32::EPSILON, 1.0 - 4.0*core::f32::EPSILON, 1.0 + 5.0*core::f32::EPSILON);

        let mask = Vector3::definitely_greater(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), 0.0);
        }
        {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(core::f32::EPSILON, -core::f32::EPSILON, 2.0 * core::f32::EPSILON);

        let mask = Vector3::definitely_greater(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), 0.0);
        }
    }
    #[test]
    fn definitely_less() {
        {
        let a = Vector3::new(1.0, 1.0, 1.0);
        let b = Vector3::new(1.0 + 4.0*core::f32::EPSILON, 1.0 - 4.0*core::f32::EPSILON, 1.0 + 5.0*core::f32::EPSILON);

        let mask = Vector3::definitely_less(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
        }
        {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(core::f32::EPSILON, -core::f32::EPSILON, 2.0 * core::f32::EPSILON);

        let mask = Vector3::definitely_less(a, b);
        let c = Vector3::new(-1000.0, -1000.0, -1000.0);
        let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
        }
    }
    #[test]
    fn cross() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, -1.0, 0.0, 0.0),
            };
            let b = Vector3 {
                data: _mm_set_ps(200.0, 0.0, 0.0, -1.0),
            };
            let c = Vector3::cross(a, b);

            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector3::new(0.0, 0.0, -1.0);
            let b = Vector3::new(-1.0, 0.0, 0.0);
            let c = Vector3::cross(b, a);

            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), -1.0);
            assert_eq!(c.z(), 0.0);
        }
    }

    #[test]
    fn abs() {
        let a = Vector3::new(-1.0, 0.0, -0.0);

        let b = Vector3::abs(a);

        assert_eq!(b.x(), 1.0);
        assert_eq!(b.y(), 0.0);
        assert_eq!(b.z(), 0.0);
        assert_eq!(a.x(), -1.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), -0.0);
    }
    #[test]
    fn copysign() {
        let a = Vector3::new(-1.0, 0.0, -0.0);
        let b = Vector3::new(10.0, -20.0, 5.0);
        let c = Vector3::copysign(b, a);

        assert_eq!(c.x(), -10.0);
        assert_eq!(c.y(), 20.0);
        assert_eq!(c.z(), -5.0);
    }

    #[test]
    fn floor() {
        {
            let a = Vector3::new(-1.1, 0.1, 0.7);
            let c = Vector3::floor(a);

            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector3::floor(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.0));
            let b = Vector3::new(0.0, 0.0, 1.0);
            let c = Vector3::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn ceil() {
        {
            let a = Vector3::new(-1.1, -0.7, 0.7);
            let c = Vector3::ceil(a);

            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 1.0);
        }
        {
            let a = Vector3::ceil(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.0));
            let b = Vector3::new(0.0, 0.0, 1.0);
            let c = Vector3::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn round() {
        {
            let a = Vector3::new(1.5, -0.5, 0.5);
            let c = Vector3::round(a);

            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector3::round(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.0));
            let b = Vector3::new(0.0, 0.0, 1.0);
            let c = Vector3::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn floor_to_int() {
        {
            let a = Vector3::new(-1.1, 0.1, 0.7);
            let c = Vector3::floor_to_int(a);

            assert_eq!(c.x(), -2);
            assert_eq!(c.y(), 0);
            assert_eq!(c.z(), 0);
        }
        {
            let a =
                Vector3::floor_to_int(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.0));

            assert_eq!(a.x(), -2147483648);
            assert_eq!(a.y(), -2147483648);

            assert_eq!(a.z(), 0);
        }
    }

    #[test]
    fn ceil_to_int() {
        {
            let a = Vector3::new(-1.1, 0.1, 0.7);
            let c = Vector3::ceil_to_int(a);

            assert_eq!(c.x(), -1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
        }
        {
            // let a =
            Vector3::ceil_to_int(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.0));

            //assert_eq!(a.x(), 2.0 * 2147483647.0);
            //assert_eq!(a.y(), -2.0 * 2147483648.0);

            //assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn truncate() {
        {
            let a = Vector3::new(-1.6, 0.1, 0.7);
            let c = Vector3::truncate(a);

            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector3::truncate(Vector3::new(2.0 * 2147483647.0, -2.0 * 2147483648.0, -0.1));
            let b = Vector3::new(0.0, 0.0, 1.0);
            let c = Vector3::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn frac() {
        {
            let a = Vector3::frac(Vector3::new(-1.75, 0.1, 0.0));

            assert_eq!(a.x(), 0.25);
            assert_eq!(a.y(), 0.1);
            assert_eq!(a.z(), 0.0);
        }
    }

    #[test]
    fn sqrt() {
        {
            let a = Vector3::sqrt(Vector3::new(4.0, 9.0, 16.0));

            assert_eq!(a.x(), 2.0);
            assert_eq!(a.y(), 3.0);
            assert_eq!(a.z(), 4.0);
        }
    }
    #[test]
    fn sin() {
        {
            let a = Vector3::sin(Vector3::new(
                0.0,
                0.5 * core::f32::consts::PI,
                core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 0.0);
            assert_eq!(a.y(), 1.0);
            assert_eq!(a.z(), 0.0);
        }
        {
            let a = Vector3::sin(Vector3::new(
                0.25 * core::f32::consts::PI,
                -0.25 * core::f32::consts::PI,
                -core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 0.70710678118);
            assert_eq!(a.y(), -0.70710678118);
            assert_eq!(a.z(), 0.0);
        }

        {
            let a = Vector3::sin(Vector3::new(
                100000.25 * core::f32::consts::PI,
                -100000.25 * core::f32::consts::PI,
                -99999.0 * core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 0.70710678118);
            assert_eq!(a.y(), -0.70710678118);
            assert_eq!(a.z(), 0.0);
        }
    }
    #[test]
    fn cos() {
        {
            let a = Vector3::cos(Vector3::new(
                0.0,
                0.5 * core::f32::consts::PI,
                core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 1.0);
            assert_eq!(a.y(), 0.0);
            assert_eq!(a.z(), -1.0);
        }
        {
            let a = Vector3::cos(Vector3::new(
                0.25 * core::f32::consts::PI,
                -0.25 * core::f32::consts::PI,
                -core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 0.70710678118);
            assert_eq!(a.y(), 0.70710678118);
            assert_eq!(a.z(), -1.0);
        }

        {
            let a = Vector3::cos(Vector3::new(
                100000.25 * core::f32::consts::PI,
                -100000.25 * core::f32::consts::PI,
                -99999.0 * core::f32::consts::PI,
            ));

            assert_eq!(a.x(), 0.70710678118);
            assert_eq!(a.y(), 0.70710678118);
            assert_eq!(a.z(), -1.0);
        }
    }
    #[test]
    fn acos() {
        {
            let a = Vector3::acos(Vector3::new(-1.0, 0.0, 1.0));

            assert_eq!(a.x(), core::f32::consts::PI);
            assert_eq!(a.y(), 0.5 * core::f32::consts::PI);
            assert_eq!(a.z(), 0.0);
        }
        /*{
        let a = Vector3::acos(Vector3::new(-0.70710678118, 0.0, 0.70710678118));

        assert_eq!(a.x(), 0.75*core::f32::consts::PI);
        assert_eq!(a.y(), 0.5*core::f32::consts::PI);
        assert_eq!(a.z(), 0.25*core::f32::consts::PI);

        }*/
    }
    #[test]
    fn max() {
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::max(a, b);
            assert_eq!(c.x(), -1.75);
            assert_eq!(c.y(), 0.2);
            assert_eq!(c.z(), 0.0);
        }
    }

    #[test]
    fn min() {
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::min(a, b);
            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.1);
            assert_eq!(c.z(), -0.1);
        }
    }

    #[test]
    fn dot() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, -1.0, 0.0, 0.0),
            };
            let b = Vector3 {
                data: _mm_set_ps(200.0, 0.0, 0.0, -1.0),
            };
            let c = Vector3::dot(a, b);
            let f = Vector4::from(c);
            assert_eq!(f.x(), 0.0);
            assert_eq!(f.y(), 0.0);
            assert_eq!(f.z(), 0.0);
            assert_eq!(f.w(), 0.0);
        }
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, -1.0, 0.0, 5.0),
            };
            let b = Vector3 {
                data: _mm_set_ps(200.0, 1.0, 0.0, 7.0),
            };
            let c = Vector3::dot(a, b);
            let f = Vector4::from(c);
            assert_eq!(f.x(), 34.0);
            assert_eq!(f.y(), 34.0);
            assert_eq!(f.z(), 34.0);
            assert_eq!(f.w(), 34.0);
        }
    }

    #[test]
    fn lerp() {
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::lerp(a, b, 0.0);
            assert_eq!(c.x(), -1.75);
            assert_eq!(c.y(), 0.1);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::lerp(a, b, 1.0);
            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.2);
            assert_eq!(c.z(), -0.1);
        }
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::lerp(a, b, 0.5);
            assert_eq!(c.x(), -1.875);
            assert_eq!(c.y(), 0.15);
            assert_eq!(c.z(), -0.05);
        }
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::lerp(a, b, 2.0);
            assert_eq!(c.x(), -2.25);
            assert_eq!(c.y(), 0.3);
            assert_eq!(c.z(), -0.2);
        }
        {
            let a = Vector3::new(-1.75, 0.1, 0.0);
            let b = Vector3::new(-2.0, 0.2, -0.1);
            let c = Vector3::lerp(a, b, -1.0);
            assert_eq!(c.x(), -1.5);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.1);
        }
    }

    #[test]
    fn normalize() {
        {
            let a = Vector3::new(-2.0, 0.0, 0.0);
            let c = Vector3::normalize(a);
            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector3::new(core::f32::MAX, 1.0, 0.0);
            let c = Vector3::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector3::new(0.0, 0.0, 0.0);
            let c = Vector3::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector3::new(core::f32::NAN, 0.0, 0.0);
            let c = Vector3::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
    }

    #[test]
    fn sqr_magnitude() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, -3.0, 0.0, 0.0),
            };

            let c = Vector3::sqr_magnitude(a);

            assert_eq!(c, 9.0);
        }
    }
    #[test]
    fn magnitude() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, -3.0, 0.0, 0.0),
            };

            let c = Vector3::magnitude(a);

            assert_eq!(c, 3.0);
        }
    }

    #[test]
    fn swizzle() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3 {
                data: _mm_set_ps(100.0, 3.0, 2.0, 1.0),
            };
            {
                let c = a.xxxx();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 1.0);
                assert_eq!(c.w(), 1.0);
            }
            {
                let c = a.yyyy();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 2.0);
                assert_eq!(c.w(), 2.0);
            }
            {
                let c = a.zzzz();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 3.0);
                assert_eq!(c.w(), 3.0);
            }
            {
                let c = a.xxx();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.xxy();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.xxz();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.xyx();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.xyy();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.xyz();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.xzx();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.xzy();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.xzz();
                assert_eq!(c.x(), 1.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 3.0);
            }

            {
                let c = a.yxx();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.yxy();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.yxz();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.yyx();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.yyy();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.yyz();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.yzx();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.yzy();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.yzz();
                assert_eq!(c.x(), 2.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 3.0);
            }

            {
                let c = a.zxx();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.zxy();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.zxz();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 1.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.zyx();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.zyy();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.zyz();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 2.0);
                assert_eq!(c.z(), 3.0);
            }
            {
                let c = a.zzx();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 1.0);
            }
            {
                let c = a.zzy();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 2.0);
            }
            {
                let c = a.zzz();
                assert_eq!(c.x(), 3.0);
                assert_eq!(c.y(), 3.0);
                assert_eq!(c.z(), 3.0);
            }
        }
    }
}
