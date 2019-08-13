// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2/.
//

#[cfg(test)]
mod test {
    use crate::int_vector::IntVector;
    use crate::raw::RawVector_i32;
    use crate::vector4_int::Vector4Int;
    use crate::vector4_bool::Vector4Bool;
    #[test]
    fn new() {
        let a = Vector4Int::new(1, 2, 3, 4);

        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), 3);
        assert_eq!(a.w(), 4);
    }
    #[test]
    fn set() {
        let a = Vector4Int::set(-5);

        assert_eq!(a.x(), -5);
        assert_eq!(a.y(), -5);
        assert_eq!(a.z(), -5);
        assert_eq!(a.w(), -5);
    }
    #[test]
    fn zero() {
        let a = Vector4Int::zero();

        assert_eq!(a.x(), 0);
        assert_eq!(a.y(), 0);
        assert_eq!(a.z(), 0);
        assert_eq!(a.w(), 0);
    }
    #[test]
    fn store() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let mut b = RawVector_i32 { data: [0; 4] };
        a.store(&mut b);
        assert_eq!(b.data[0], 1);
        assert_eq!(b.data[1], 2);
        assert_eq!(b.data[2], 3);
        assert_eq!(b.data[3], 4);
        b.data[0] += 1;
        b.data[3] += 1;
        let c = Vector4Int::load(&b);
        let e = c + a;
        e.store(&mut b);
        assert_eq!(b.data[0], 3);
        assert_eq!(b.data[1], 4);
        assert_eq!(b.data[2], 6);
        assert_eq!(b.data[3], 9);
    }

    #[test]
    fn set_x() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a.set_x(-11);
        assert_eq!(a.x(), -11);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), 3);
        assert_eq!(a.w(), 4);
    }
    #[test]
    fn set_y() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a.set_y(-11);
        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), -11);
        assert_eq!(a.z(), 3);
        assert_eq!(a.w(), 4);
    }
    #[test]
    fn set_z() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a.set_z(-11);
        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), -11);
        assert_eq!(a.w(), 4);
    }
    #[test]
    fn set_w() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a.set_w(-11);
        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), 3);
        assert_eq!(a.w(), -11);
    }
    #[test]
    fn add() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(4, 4, 6, 4);

        let c = a + b;
        assert_eq!(c.x(), 5);
        assert_eq!(c.y(), 6);
        assert_eq!(c.z(), 9);
        assert_eq!(c.w(), 8);
    }
    #[test]
    fn add_assign() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a += Vector4Int::set(3);

        assert_eq!(a.x(), 4);
        assert_eq!(a.y(), 5);
        assert_eq!(a.z(), 6);
        assert_eq!(a.w(), 7);
    }
    #[test]
    fn sub() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(4, 4, 6, 4);

        let c = a - b;
        assert_eq!(c.x(), -3);
        assert_eq!(c.y(), -2);
        assert_eq!(c.z(), -3);
        assert_eq!(c.w(), 0);
    }
    #[test]
    fn sub_assign() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a -= Vector4Int::new(4, 4, 6, 4);

        assert_eq!(a.x(), -3);
        assert_eq!(a.y(), -2);
        assert_eq!(a.z(), -3);
        assert_eq!(a.w(), 0);
    }
    #[test]
    fn neg() {
        let a = Vector4Int::new(1, -2, 3, 4);
        let c = -a;

        assert_eq!(c.x(), -1);
        assert_eq!(c.y(), 2);
        assert_eq!(c.z(), -3);
        assert_eq!(c.w(), -4);
    }
    #[test]
    fn component_mul() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(4, 4, 6, 4);

        let c = Vector4Int::mul(a, b);
        assert_eq!(c.x(), 4);
        assert_eq!(c.y(), 8);
        assert_eq!(c.z(), 18);
        assert_eq!(c.w(), 16);
    }


    #[test]
    fn mul() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = IntVector::new(-3);
        {
            let c = a * b;
            assert_eq!(c.x(), -3);
            assert_eq!(c.y(), -6);
            assert_eq!(c.z(), -9);
            assert_eq!(c.w(), -12);
        }
        {
            let c = b * a;
            assert_eq!(c.x(), -3);
            assert_eq!(c.y(), -6);
            assert_eq!(c.z(), -9);
            assert_eq!(c.w(), -12);
        }
    }
    #[test]
    fn mul_assign() {
        let mut a = Vector4Int::new(1, 2, 3, 4);
        a *= -3;

        assert_eq!(a.x(), -3);
        assert_eq!(a.y(), -6);
        assert_eq!(a.z(), -9);
        assert_eq!(a.w(), -12);
    }

    #[test]
    fn and() {
        let a = Vector4Int::new(1, -2, 0, 4);
        let b = Vector4Int::new(3, -2, -3, 0);

        let c = Vector4Int::and(a, b);
        assert_eq!(c.x(), 1);
        assert_eq!(c.y(), -2);
        assert_eq!(c.z(), 0);
        assert_eq!(c.w(), 0);
    }
    #[test]
    fn andnot() {
        let a = Vector4Int::new(0, 1, 0, 4);
        let b = Vector4Int::new(1, 3, -3, 0);

        let c = Vector4Int::andnot(a, b);
        assert_eq!(c.x(), 1);
        assert_eq!(c.y(), 2);
        assert_eq!(c.z(), -3);
        assert_eq!(c.w(), 0);
    }
    #[test]
    fn or() {
        let a = Vector4Int::new(1, 2, 3, 0);
        let b = Vector4Int::new(2, 1, 3, 0);

        let c = Vector4Int::or(a, b);
        assert_eq!(c.x(), 3);
        assert_eq!(c.y(), 3);
        assert_eq!(c.z(), 3);
        assert_eq!(c.w(), 0);
    }
    #[test]
    fn xor() {
        let a = Vector4Int::new(1, 0, 3, 0);
        let b = Vector4Int::new(3, 3, 5, 0);

        let c = Vector4Int::xor(a, b);
        assert_eq!(c.x(), 2);
        assert_eq!(c.y(), 3);
        assert_eq!(c.z(), 6);
        assert_eq!(c.w(), 0);
    }
    #[test]
    fn all() {
        {
            let a = Vector4Int::new(1, 2, 3, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), true);
        }
        {
            let a = Vector4Int::new(0, 2, 3, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4Int::new(1, 0, 3, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4Int::new(1, 2, 0, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4Int::new(1, 2, 3, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
        {
            let a = Vector4Int::new(0, 0, 0, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::all(c), false);
        }
    }
    #[test]
    fn any() {
        {
            let a = Vector4Int::new(1, 2, 3, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4Int::new(0, 0, 3, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4Int::new(1, 0, 0, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4Int::new(0, 2, 0, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4Int::new(0, 0, 0, 4);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), true);
        }
        {
            let a = Vector4Int::new(0, 0, 0, 0);
            let b = Vector4Int::new(1, 2, 3, 4);

            let c = Vector4Int::equal(a, b);
            assert_eq!(Vector4Bool::any(c), false);
        }
    }
    #[test]
    fn equal() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::equal(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn not_equal() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::not_equal(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), -1000);
    }
    #[test]
    fn greater_equal() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::greater_equal(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn greater() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::greater(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn less_equal() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::less_equal(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), -1000);
    }
    #[test]
    fn less() {
        let a = Vector4Int::new(1, 2, 3, 4);
        let b = Vector4Int::new(1, 1, 4, 4);

        let mask = Vector4Int::less(a, b);
        let c = Vector4Int::new(-1000, -1000, -1000, -1000);
        let d = Vector4Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), -1000);
    }
    
    #[test]
    fn abs() {
        let a = Vector4Int::new(-1, 0, -0, 4);

        let b = Vector4Int::abs(a);

        assert_eq!(b.x(), 1);
        assert_eq!(b.y(), 0);
        assert_eq!(b.z(), 0);
        assert_eq!(a.x(), -1);
        assert_eq!(a.y(), 0);
        assert_eq!(a.z(), -0);
    }
    #[test]
    fn sign() {
        let a = Vector4Int::new(10, -20, 5, -4);
        let b = Vector4Int::new(-1, 1, -0, -4);
        
        let c = Vector4Int::sign(a, b);

        assert_eq!(c.x(), -10);
        assert_eq!(c.y(), -20);
        assert_eq!(c.z(), 0);
        assert_eq!(c.w(), 4);
    }




    #[test]
    fn max() {
        {
            let a = Vector4Int::new(-175, 10, 0, 4);
            let b = Vector4Int::new(-200, 20, -10, 4);
            let c = Vector4Int::max(a, b);
            assert_eq!(c.x(), -175);
            assert_eq!(c.y(), 20);
            assert_eq!(c.z(), 0);
        }
    }

    #[test]
    fn min() {
        {
            let a = Vector4Int::new(-175, 10, 0, 4);
            let b = Vector4Int::new(-200, 20, -10, 4);
            let c = Vector4Int::min(a, b);
            assert_eq!(c.x(), -200);
            assert_eq!(c.y(), 10);
            assert_eq!(c.z(), -10);
        }
    }




    #[test]
    fn swizzle() {
        //OMG total swizzle coverage sux this should be a macro.  TODO: investigate a better way.
        let a = Vector4Int::new(1, 2, 3, 4);
        {
            let c = a.xxxx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xxxy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xxxz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xxxw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xxyx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xxyy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xxyz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xxyw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xxzx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xxzy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xxzz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xxzw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xxwx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xxwy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xxwz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xxww();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xyxx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xyxy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xyxz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xyxw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xyyx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xyyy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xyyz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xyyw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xyzx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xyzy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xyzz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xyzw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xywx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xywy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xywz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xyww();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xzxx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xzxy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xzxz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xzxw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xzyx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xzyy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xzyz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xzyw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xzzx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xzzy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xzzz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xzzw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xzwx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xzwy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xzwz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xzww();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xwxx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xwxy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xwxz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xwxw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xwyx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xwyy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xwyz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xwyw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xwzx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xwzy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xwzz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xwzw();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.xwwx();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.xwwy();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.xwwz();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.xwww();
            assert_eq!(c.x(), 1);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yxxx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yxxy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yxxz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yxxw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yxyx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yxyy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yxyz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yxyw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yxzx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yxzy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yxzz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yxzw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yxwx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yxwy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yxwz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yxww();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yyxx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yyxy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yyxz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yyxw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yyyx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yyyy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yyyz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yyyw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yyzx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yyzy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yyzz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yyzw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yywx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yywy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yywz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yyww();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yzxx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yzxy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yzxz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yzxw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yzyx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yzyy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yzyz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yzyw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yzzx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yzzy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yzzz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yzzw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.yzwx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.yzwy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.yzwz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.yzww();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.ywxx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.ywxy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.ywxz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.ywxw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.ywyx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.ywyy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.ywyz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.ywyw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.ywzx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.ywzy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.ywzz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.ywzw();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.ywwx();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.ywwy();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.ywwz();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.ywww();
            assert_eq!(c.x(), 2);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zxxx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zxxy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zxxz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zxxw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zxyx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zxyy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zxyz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zxyw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zxzx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zxzy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zxzz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zxzw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zxwx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zxwy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zxwz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zxww();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zyxx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zyxy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zyxz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zyxw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zyyx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zyyy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zyyz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zyyw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zyzx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zyzy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zyzz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zyzw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zywx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zywy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zywz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zyww();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zzxx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zzxy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zzxz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zzxw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zzyx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zzyy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zzyz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zzyw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zzzx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zzzy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zzzz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zzzw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zzwx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zzwy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zzwz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zzww();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zwxx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zwxy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zwxz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zwxw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zwyx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zwyy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zwyz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zwyw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zwzx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zwzy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zwzz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zwzw();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.zwwx();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.zwwy();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.zwwz();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.zwww();
            assert_eq!(c.x(), 3);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wxxx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wxxy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wxxz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wxxw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wxyx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wxyy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wxyz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wxyw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wxzx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wxzy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wxzz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wxzw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wxwx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wxwy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wxwz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wxww();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 1);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wyxx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wyxy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wyxz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wyxw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wyyx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wyyy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wyyz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wyyw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wyzx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wyzy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wyzz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wyzw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wywx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wywy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wywz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wyww();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 2);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wzxx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wzxy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wzxz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wzxw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wzyx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wzyy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wzyz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wzyw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wzzx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wzzy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wzzz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wzzw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wzwx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wzwy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wzwz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wzww();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 3);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wwxx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wwxy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wwxz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wwxw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 1);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wwyx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wwyy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wwyz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wwyw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 2);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wwzx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wwzy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wwzz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wwzw();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 3);
            assert_eq!(c.w(), 4);
        }
        {
            let c = a.wwwx();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 1);
        }
        {
            let c = a.wwwy();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 2);
        }
        {
            let c = a.wwwz();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 3);
        }
        {
            let c = a.wwww();
            assert_eq!(c.x(), 4);
            assert_eq!(c.y(), 4);
            assert_eq!(c.z(), 4);
            assert_eq!(c.w(), 4);
        }
    }
}
