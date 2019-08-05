// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#[cfg(test)]
mod test {
    use crate::float_vector::FloatVector;
    use crate::raw::RawVector_f32;
    use crate::vector4::Vector4;
    use crate::vector4_bool::Vector4Bool;
    #[test]
    fn new() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);

        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set() {
        let a = Vector4::set(-5.0);

        assert_eq!(a.x(), -5.0);
        assert_eq!(a.y(), -5.0);
        assert_eq!(a.z(), -5.0);
        assert_eq!(a.w(), -5.0);
    }
    #[test]
    fn zero() {
        let a = Vector4::zero();

        assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), 0.0);
        assert_eq!(a.w(), 0.0);
    }
    #[test]
    fn store() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let mut b = RawVector_f32 { data: [0.0; 4] };
        a.store(&mut b);
        assert_eq!(b.data[0], 1.0);
        assert_eq!(b.data[1], 2.0);
        assert_eq!(b.data[2], 3.0);
        assert_eq!(b.data[3], 4.0);
        b.data[0] += 1.0;
        b.data[3] += 1.0;
        let c = Vector4::load(&b);
        let e = c + a;
        e.store(&mut b);
        assert_eq!(b.data[0], 3.0);
        assert_eq!(b.data[1], 4.0);
        assert_eq!(b.data[2], 6.0);
        assert_eq!(b.data[3], 9.0);
    }

    #[test]
    fn set_x() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a.set_x(-11.0);
        assert_eq!(a.x(), -11.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_y() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a.set_y(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), -11.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_z() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a.set_z(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), -11.0);
        assert_eq!(a.w(), 4.0);
    }
    #[test]
    fn set_w() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a.set_w(-11.0);
        assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        assert_eq!(a.w(), -11.0);
    }
    #[test]
    fn add() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);

        let c = a + b;
        assert_eq!(c.x(), 5.0);
        assert_eq!(c.y(), 6.0);
        assert_eq!(c.z(), 9.0);
        assert_eq!(c.w(), 8.0);
    }
    #[test]
    fn add_assign() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a += Vector4::set(3.0);

        assert_eq!(a.x(), 4.0);
        assert_eq!(a.y(), 5.0);
        assert_eq!(a.z(), 6.0);
        assert_eq!(a.w(), 7.0);
    }
    #[test]
    fn sub() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);

        let c = a - b;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        assert_eq!(c.w(), 0.0);
    }
    #[test]
    fn sub_assign() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a -= Vector4::new(4.0, 4.0, 6.0, 4.0);

        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -2.0);
        assert_eq!(a.z(), -3.0);
        assert_eq!(a.w(), 0.0);
    }
    #[test]
    fn neg() {
        let a = Vector4::new(1.0, -2.0, 3.0, 4.0);
        let c = -a;

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 2.0);
        assert_eq!(c.z(), -3.0);
        assert_eq!(c.w(), -4.0);
    }
    #[test]
    fn component_mul() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);

        let c = Vector4::mul(a, b);
        assert_eq!(c.x(), 4.0);
        assert_eq!(c.y(), 8.0);
        assert_eq!(c.z(), 18.0);
        assert_eq!(c.w(), 16.0);
    }
    #[test]
    fn component_div() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);

        let c = Vector4::div(a, b);
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.5);
        assert_eq!(c.w(), 1.0);
    }
    #[test]
    fn mul_add() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);
        let c = Vector4::new(2.0, -5.0, 10.0, 4.0);
        let d = Vector4::mul_add(a, b, c);

        assert_eq!(d.x(), 6.0);
        assert_eq!(d.y(), 3.0);
        assert_eq!(d.z(), 28.0);
        assert_eq!(d.w(), 20.0);
    }
    #[test]
    fn mul_sub() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);
        let c = Vector4::new(2.0, -5.0, 10.0, 4.0);
        let d = Vector4::mul_sub(a, b, c);

        assert_eq!(d.x(), 2.0);
        assert_eq!(d.y(), 13.0);
        assert_eq!(d.z(), 8.0);
        assert_eq!(d.w(), 12.0);
    }
    #[test]
    fn neg_mul_add() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);
        let c = Vector4::new(2.0, -5.0, 10.0, 4.0);
        let d = Vector4::neg_mul_add(a, b, c);

        assert_eq!(d.x(), -2.0);
        assert_eq!(d.y(), -13.0);
        assert_eq!(d.z(), -8.0);
        assert_eq!(d.w(), -12.0);
    }
    #[test]
    fn neg_mul_sub() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(4.0, 4.0, 6.0, 4.0);
        let c = Vector4::new(2.0, -5.0, 10.0, 4.0);
        let d = Vector4::neg_mul_sub(a, b, c);

        assert_eq!(d.x(), -6.0);
        assert_eq!(d.y(), -3.0);
        assert_eq!(d.z(), -28.0);
        assert_eq!(d.w(), -20.0);
    }
    #[test]
    fn mul() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = FloatVector::new(-3.0);
        {
            let c = a * b;
            assert_eq!(c.x(), -3.0);
            assert_eq!(c.y(), -6.0);
            assert_eq!(c.z(), -9.0);
            assert_eq!(c.w(), -12.0);
        }
        {
            let c = b * a;
            assert_eq!(c.x(), -3.0);
            assert_eq!(c.y(), -6.0);
            assert_eq!(c.z(), -9.0);
            assert_eq!(c.w(), -12.0);
        }
    }
    #[test]
    fn mul_assign() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a *= -3.0;

        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -6.0);
        assert_eq!(a.z(), -9.0);
        assert_eq!(a.w(), -12.0);
    }
    #[test]
    fn div() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = 4.0;

        let c = a / b;
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.75);
        assert_eq!(c.w(), 1.0);
    }
    #[test]
    fn div_inv() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = FloatVector::new(1.0);

        let c = b / a;
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 1.0 / 3.0);
        assert_eq!(c.w(), 0.25);
    }
    #[test]
    fn div_assign() {
        let mut a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        a /= 4.0;

        assert_eq!(a.x(), 0.25);
        assert_eq!(a.y(), 0.5);
        assert_eq!(a.z(), 0.75);
        assert_eq!(a.w(), 1.0);
    }
    #[test]
    fn and() {
        let a = Vector4::new(-1.0, -2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, -2.0, -3.0, 4.0);

        let c = Vector4::and(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
        assert_eq!(c.w(), 4.0);
    }
    #[test]
    fn andnot() {
        let a = Vector4::new(0.0, 0.0, -0.0, 4.0);
        let b = Vector4::new(1.0, -2.0, -3.0, 4.0);

        let c = Vector4::andnot(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
        assert_eq!(c.w(), 0.0);
    }
    #[test]
    fn or() {
        let a = Vector4::new(-1.0, -2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, -2.0, -3.0, 4.0);

        let c = Vector4::or(a, b);
        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        assert_eq!(c.w(), 4.0);
    }
    #[test]
    fn xor() {
        let a = Vector4::new(-1.0, -2.0, 3.0, 4.0);
        let b = Vector4::new(-0.0, 0.0, -0.0, 4.0);

        let c = Vector4::xor(a, b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        assert_eq!(c.w(), 0.0);
    }
    #[test]
    fn all() {
        {
            let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), true);
        }
        {
            let a = Vector4::new(0.0, 2.0, 3.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4::new(1.0, 0.0, 3.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4::new(1.0, 2.0, 0.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4::new(1.0, 2.0, 3.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
    }
    #[test]
    fn any() {
        {
            let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4::new(0.0, 0.0, 3.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4::new(1.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4::new(0.0, 2.0, 0.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 4.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(1.0, 2.0, 3.0, 4.0);

            let c = Vector4::equal(a, b);
            assert_eq!(Vector4Bool::any(c), false);
        }
    }
    #[test]
    fn equal() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::equal(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn not_equal() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::not_equal(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), -1000.0);
    }
    #[test]
    fn greater_equal() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::greater_equal(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn greater() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::greater(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
    }
    #[test]
    fn less_equal() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::less_equal(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
    }
    #[test]
    fn less() {
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::new(1.0, 1.0, 4.0, 4.0);

        let mask = Vector4::less(a, b);
        let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
        let d = Vector4::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
    }
    #[test]
    fn approx_equal() {
        {
            let a = Vector4::new(1.0, 1.0, 1.0, 1.0);
            let b = Vector4::new(
                1.0 + 4.0 * core::f32::EPSILON,
                1.0 - 4.0 * core::f32::EPSILON,
                1.0 + 5.0 * core::f32::EPSILON,
                1.0 - 5.0 * core::f32::EPSILON,
            );

            let mask = Vector4::approx_equal(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), -1000.0);
            assert_eq!(d.y(), -1000.0);
            assert_eq!(d.z(), 0.0);
            assert_eq!(d.w(), 0.0);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(
                core::f32::EPSILON,
                -core::f32::EPSILON,
                2.0 * core::f32::EPSILON,
                -2.0 * core::f32::EPSILON,
            );

            let mask = Vector4::approx_equal(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), -1000.0);
            assert_eq!(d.y(), -1000.0);
            assert_eq!(d.z(), 0.0);
            assert_eq!(d.w(), 0.0);
        }
    }
    #[test]
    fn definitely_greater() {
        {
            let a = Vector4::new(1.0, 1.0, 1.0, 1.0);
            let b = Vector4::new(
                1.0 + 4.0 * core::f32::EPSILON,
                1.0 - 4.0 * core::f32::EPSILON,
                1.0 + 5.0 * core::f32::EPSILON,
                1.0 - 5.0 * core::f32::EPSILON,
            );

            let mask = Vector4::definitely_greater(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), 0.0);
            assert_eq!(d.y(), 0.0);
            assert_eq!(d.z(), 0.0);
            assert_eq!(d.w(), -1000.0);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(
                core::f32::EPSILON,
                -core::f32::EPSILON,
                2.0 * core::f32::EPSILON,
                -2.0 * core::f32::EPSILON,
            );

            let mask = Vector4::definitely_greater(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), 0.0);
            assert_eq!(d.y(), 0.0);
            assert_eq!(d.z(), 0.0);
            assert_eq!(d.w(), -1000.0);
        }
    }
    #[test]
    fn definitely_less() {
        {
            let a = Vector4::new(1.0, 1.0, 1.0, 1.0);
            let b = Vector4::new(
                1.0 + 4.0 * core::f32::EPSILON,
                1.0 - 4.0 * core::f32::EPSILON,
                1.0 + 5.0 * core::f32::EPSILON,
                1.0 - 5.0 * core::f32::EPSILON,
            );

            let mask = Vector4::definitely_less(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), 0.0);
            assert_eq!(d.y(), 0.0);
            assert_eq!(d.z(), -1000.0);
            assert_eq!(d.w(), 0.0);
        }
        {
            let a = Vector4::new(0.0, 0.0, 0.0, 0.0);
            let b = Vector4::new(
                core::f32::EPSILON,
                -core::f32::EPSILON,
                2.0 * core::f32::EPSILON,
                -2.0 * core::f32::EPSILON,
            );

            let mask = Vector4::definitely_less(a, b);
            let c = Vector4::new(-1000.0, -1000.0, -1000.0, -1000.0);
            let d = Vector4::and(c, mask);

            assert_eq!(d.x(), 0.0);
            assert_eq!(d.y(), 0.0);
            assert_eq!(d.z(), -1000.0);
            assert_eq!(d.w(), 0.0);
        }
    }
    #[test]
    fn abs() {
        let a = Vector4::new(-1.0, 0.0, -0.0, 4.0);

        let b = Vector4::abs(a);

        assert_eq!(b.x(), 1.0);
        assert_eq!(b.y(), 0.0);
        assert_eq!(b.z(), 0.0);
        assert_eq!(a.x(), -1.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), -0.0);
    }
    #[test]
    fn copysign() {
        let a = Vector4::new(-1.0, 0.0, -0.0, 4.0);
        let b = Vector4::new(10.0, -20.0, 5.0, 4.0);
        let c = Vector4::copysign(b, a);

        assert_eq!(c.x(), -10.0);
        assert_eq!(c.y(), 20.0);
        assert_eq!(c.z(), -5.0);
    }

    #[test]
    fn floor() {
        {
            let a = Vector4::new(-1.1, 0.1, 0.7, 4.0);
            let c = Vector4::floor(a);

            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector4::floor(Vector4::new(
                2.0 * 2147483647.0,
                -2.0 * 2147483648.0,
                -0.0,
                1.0,
            ));
            let b = Vector4::new(0.0, 0.0, 1.0, -1.0);
            let c = Vector4::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn ceil() {
        {
            let a = Vector4::new(-1.1, -0.7, 0.7, 4.0);
            let c = Vector4::ceil(a);

            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 1.0);
        }
        {
            let a = Vector4::ceil(Vector4::new(
                2.0 * 2147483647.0,
                -2.0 * 2147483648.0,
                -0.0,
                1.0,
            ));
            let b = Vector4::new(0.0, 0.0, 1.0, 1.0);
            let c = Vector4::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn round() {
        {
            let a = Vector4::new(1.5, -0.5, 0.5, 4.0);
            let c = Vector4::round(a);

            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector4::round(Vector4::new(
                2.0 * 2147483647.0,
                -2.0 * 2147483648.0,
                -0.0,
                1.0,
            ));
            let b = Vector4::new(0.0, 0.0, 1.0, 1.0);
            let c = Vector4::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn floor_to_int() {
        {
            let a = Vector4::new(-1.1, 0.1, 0.7, 4.0);
            let c = Vector4::floor_to_int(a);

            assert_eq!(c.x(), -2);
            assert_eq!(c.y(), 0);
            assert_eq!(c.z(), 0);
        }
        {
            let a = Vector4::floor_to_int(Vector4::new(
                2.0 * 2147483647.0,
                -2.0 * 2147483648.0,
                -0.0,
                1.0,
            ));

            assert_eq!(a.x(), -2147483648);
            assert_eq!(a.y(), -2147483648);

            assert_eq!(a.z(), 0);
        }
    }

    #[test]
    fn ceil_to_int() {
        {
            let a = Vector4::new(-1.1, 0.1, 0.7, 4.0);
            let c = Vector4::ceil_to_int(a);

            assert_eq!(c.x(), -1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
        }
        {
            // let a = Vector4::ceil_to_int(Vector4::new(
            //     2.0 * 2147483647.0,
            //     -2.0 * 2147483648.0,
            //     -0.0,
            //     1.0,
            // ));

            //assert_eq!(a.x(), 2.0 * 2147483647.0);
            //assert_eq!(a.y(), -2.0 * 2147483648.0);

            //assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn truncate() {
        {
            let a = Vector4::new(-1.6, 0.1, 0.7, 4.0);
            let c = Vector4::truncate(a);

            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector4::truncate(Vector4::new(
                2.0 * 2147483647.0,
                -2.0 * 2147483648.0,
                -0.1,
                1.0,
            ));
            let b = Vector4::new(0.0, 0.0, 1.0, 1.0);
            let c = Vector4::copysign(b, a);
            assert_eq!(a.x(), 2.0 * 2147483647.0);
            assert_eq!(a.y(), -2.0 * 2147483648.0);

            assert_eq!(c.z(), -1.0);
        }
    }

    #[test]
    fn frac() {
        {
            let a = Vector4::frac(Vector4::new(-1.75, 0.1, 0.0, 4.0));

            assert_eq!(a.x(), 0.25);
            assert_eq!(a.y(), 0.1);
            assert_eq!(a.z(), 0.0);
        }
    }

    #[test]
    fn sqrt() {
        {
            let a = Vector4::sqrt(Vector4::new(4.0, 9.0, 16.0, 4.0));

            assert_eq!(a.x(), 2.0);
            assert_eq!(a.y(), 3.0);
            assert_eq!(a.z(), 4.0);
        }
    }

    #[test]
    fn max() {
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::max(a, b);
            assert_eq!(c.x(), -1.75);
            assert_eq!(c.y(), 0.2);
            assert_eq!(c.z(), 0.0);
        }
    }

    #[test]
    fn min() {
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::min(a, b);
            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.1);
            assert_eq!(c.z(), -0.1);
        }
    }

    #[test]
    fn dot() {
        {
            // make sure this works with garbage in w.
            let a = Vector4::new(100.0, -1.0, 2.0, 1.0);
            let b = Vector4::new(10.0, 1.0, 6.0, -1.0);
            let c = Vector4::from(Vector4::dot(a, b));

            assert_eq!(c.x(), 1010.0);
            assert_eq!(c.y(), 1010.0);
            assert_eq!(c.z(), 1010.0);
            assert_eq!(c.w(), 1010.0);
        }
        {
            // make sure this works with garbage in w.
            let a = Vector4::new(3.0, -1.0, -3.0, 5.0);
            let b = Vector4::new(2.0, 1.0, 1.0, 7.0);
            let c = Vector4::from(Vector4::dot(a, b));

            assert_eq!(c.x(), 37.0);
            assert_eq!(c.y(), 37.0);
            assert_eq!(c.z(), 37.0);
            assert_eq!(c.w(), 37.0);
        }
    }

    #[test]
    fn lerp() {
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::lerp(a, b, 0.0);
            assert_eq!(c.x(), -1.75);
            assert_eq!(c.y(), 0.1);
            assert_eq!(c.z(), 0.0);
        }
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::lerp(a, b, 1.0);
            assert_eq!(c.x(), -2.0);
            assert_eq!(c.y(), 0.2);
            assert_eq!(c.z(), -0.1);
        }
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::lerp(a, b, 0.5);
            assert_eq!(c.x(), -1.875);
            assert_eq!(c.y(), 0.15);
            assert_eq!(c.z(), -0.05);
        }
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::lerp(a, b, 2.0);
            assert_eq!(c.x(), -2.25);
            assert_eq!(c.y(), 0.3);
            assert_eq!(c.z(), -0.2);
        }
        {
            let a = Vector4::new(-1.75, 0.1, 0.0, 4.0);
            let b = Vector4::new(-2.0, 0.2, -0.1, 4.0);
            let c = Vector4::lerp(a, b, -1.0);
            assert_eq!(c.x(), -1.5);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.1);
        }
    }

    #[test]
    fn normalize() {
        {
            let a = Vector4::new(-2.0, 0.0, 0.0, 0.0);
            let c = Vector4::normalize(a);
            assert_eq!(c.x(), -1.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector4::new(core::f32::MAX, 1.0, 0.0, 4.0);
            let c = Vector4::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector4::new(0.0, 0.0, 0.0, 4.0);
            let c = Vector4::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
        {
            //SHOULD NOT BE NAN
            let a = Vector4::new(core::f32::NAN, 0.0, 0.0, 4.0);
            let c = Vector4::normalize(a);
            assert_eq!(c.x(), 0.0);
            assert_eq!(c.y(), 0.0);
            assert_eq!(c.z(), 0.0);
        }
    }

    #[test]
    fn sqr_magnitude() {
        let a = Vector4::new(0.0, 0.0, 0.0, -3.0);

        let c = Vector4::sqr_magnitude(a);

        assert_eq!(c, 9.0);
    }
    #[test]
    fn magnitude() {
        let a = Vector4::new(2.0, 2.0, 2.0, 2.0);

        let c = Vector4::magnitude(a);

        assert_eq!(c, 4.0);
    }

    #[test]
    fn swizzle() {
        //OMG total swizzle coverage sux
        let a = Vector4::new(1.0, 2.0, 3.0, 4.0);
        {
            let c = a.xxxx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xxxy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xxxz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xxxw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xxyx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xxyy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xxyz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xxyw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xxzx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xxzy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xxzz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xxzw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xxwx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xxwy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xxwz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xxww();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xyxx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xyxy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xyxz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xyxw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xyyx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xyyy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xyyz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xyyw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xyzx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xyzy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xyzz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xyzw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xywx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xywy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xywz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xyww();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xzxx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xzxy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xzxz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xzxw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xzyx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xzyy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xzyz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xzyw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xzzx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xzzy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xzzz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xzzw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xzwx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xzwy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xzwz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xzww();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xwxx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xwxy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xwxz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xwxw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xwyx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xwyy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xwyz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xwyw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xwzx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xwzy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xwzz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xwzw();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.xwwx();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.xwwy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.xwwz();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.xwww();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yxxx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yxxy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yxxz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yxxw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yxyx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yxyy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yxyz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yxyw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yxzx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yxzy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yxzz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yxzw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yxwx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yxwy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yxwz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yxww();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yyxx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yyxy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yyxz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yyxw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yyyx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
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
            let c = a.yyyz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yyyw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yyzx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yyzy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yyzz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yyzw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yywx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yywy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yywz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yyww();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yzxx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yzxy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yzxz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yzxw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yzyx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yzyy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yzyz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yzyw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yzzx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yzzy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yzzz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yzzw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.yzwx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.yzwy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.yzwz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.yzww();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.ywxx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.ywxy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.ywxz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.ywxw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.ywyx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.ywyy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.ywyz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.ywyw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.ywzx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.ywzy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.ywzz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.ywzw();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.ywwx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.ywwy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.ywwz();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.ywww();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zxxx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zxxy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zxxz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zxxw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zxyx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zxyy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zxyz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zxyw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zxzx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zxzy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zxzz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zxzw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zxwx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zxwy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zxwz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zxww();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zyxx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zyxy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zyxz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zyxw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zyyx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zyyy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zyyz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zyyw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zyzx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zyzy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zyzz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zyzw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zywx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zywy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zywz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zyww();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zzxx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zzxy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zzxz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zzxw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zzyx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zzyy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zzyz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zzyw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zzzx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zzzy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
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
            let c = a.zzzw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zzwx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zzwy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zzwz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zzww();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zwxx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zwxy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zwxz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zwxw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zwyx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zwyy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zwyz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zwyw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zwzx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zwzy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zwzz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zwzw();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.zwwx();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.zwwy();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.zwwz();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.zwww();
            assert_eq!(c.x(), 3.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wxxx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wxxy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wxxz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wxxw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wxyx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wxyy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wxyz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wxyw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wxzx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wxzy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wxzz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wxzw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wxwx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wxwy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wxwz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wxww();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 1.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wyxx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wyxy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wyxz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wyxw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wyyx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wyyy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wyyz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wyyw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wyzx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wyzy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wyzz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wyzw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wywx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wywy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wywz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wyww();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 2.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wzxx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wzxy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wzxz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wzxw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wzyx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wzyy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wzyz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wzyw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wzzx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wzzy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wzzz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wzzw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wzwx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wzwy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wzwz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wzww();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 3.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wwxx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wwxy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wwxz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wwxw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 1.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wwyx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wwyy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wwyz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wwyw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 2.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wwzx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wwzy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wwzz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wwzw();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 3.0);
            assert_eq!(c.w(), 4.0);
        }
        {
            let c = a.wwwx();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 1.0);
        }
        {
            let c = a.wwwy();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 2.0);
        }
        {
            let c = a.wwwz();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 3.0);
        }
        {
            let c = a.wwww();
            assert_eq!(c.x(), 4.0);
            assert_eq!(c.y(), 4.0);
            assert_eq!(c.z(), 4.0);
            assert_eq!(c.w(), 4.0);
        }
    }
}
