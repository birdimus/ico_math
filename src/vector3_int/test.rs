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
    use crate::vector3_bool::Vector3Bool;
    use crate::vector3_int::Vector3Int;
    #[test]
    fn new() {
        let a = Vector3Int::new(1, 2, 3);

        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), 3);
    }
    #[test]
    fn set() {
        let a = Vector3Int::set(-5);

        assert_eq!(a.x(), -5);
        assert_eq!(a.y(), -5);
        assert_eq!(a.z(), -5);
    }
    #[test]
    fn zero() {
        let a = Vector3Int::zero();

        assert_eq!(a.x(), 0);
        assert_eq!(a.y(), 0);
        assert_eq!(a.z(), 0);
    }
    #[test]
    fn store() {
        let a = Vector3Int::new(1, 2, 3);
        let mut b = RawVector_i32 { data: [0; 4] };
        a.store(&mut b);
        assert_eq!(b.data[0], 1);
        assert_eq!(b.data[1], 2);
        assert_eq!(b.data[2], 3);
        b.data[0] += 1;
        b.data[3] += 1;
        let c = Vector3Int::load(&b);
        let e = c + a;
        e.store(&mut b);
        assert_eq!(b.data[0], 3);
        assert_eq!(b.data[1], 4);
        assert_eq!(b.data[2], 6);
    }

    #[test]
    fn set_x() {
        let mut a = Vector3Int::new(1, 2, 3);
        a.set_x(-11);
        assert_eq!(a.x(), -11);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), 3);
    }
    #[test]
    fn set_y() {
        let mut a = Vector3Int::new(1, 2, 3);
        a.set_y(-11);
        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), -11);
        assert_eq!(a.z(), 3);
    }
    #[test]
    fn set_z() {
        let mut a = Vector3Int::new(1, 2, 3);
        a.set_z(-11);
        assert_eq!(a.x(), 1);
        assert_eq!(a.y(), 2);
        assert_eq!(a.z(), -11);
    }

    #[test]
    fn add() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(4, 4, 6);

        let c = a + b;
        assert_eq!(c.x(), 5);
        assert_eq!(c.y(), 6);
        assert_eq!(c.z(), 9);
    }
    #[test]
    fn add_assign() {
        let mut a = Vector3Int::new(1, 2, 3);
        a += Vector3Int::set(3);

        assert_eq!(a.x(), 4);
        assert_eq!(a.y(), 5);
        assert_eq!(a.z(), 6);
    }
    #[test]
    fn sub() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(4, 4, 6);

        let c = a - b;
        assert_eq!(c.x(), -3);
        assert_eq!(c.y(), -2);
        assert_eq!(c.z(), -3);
    }
    #[test]
    fn sub_assign() {
        let mut a = Vector3Int::new(1, 2, 3);
        a -= Vector3Int::new(4, 4, 6);

        assert_eq!(a.x(), -3);
        assert_eq!(a.y(), -2);
        assert_eq!(a.z(), -3);
    }
    #[test]
    fn neg() {
        let a = Vector3Int::new(1, -2, 3);
        let c = -a;

        assert_eq!(c.x(), -1);
        assert_eq!(c.y(), 2);
        assert_eq!(c.z(), -3);
    }
    #[test]
    fn component_mul() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(4, 4, 6);

        let c = Vector3Int::mul(a, b);
        assert_eq!(c.x(), 4);
        assert_eq!(c.y(), 8);
        assert_eq!(c.z(), 18);
    }

    #[test]
    fn mul() {
        let a = Vector3Int::new(1, 2, 3);
        let b = IntVector::new(-3);
        {
            let c = a * b;
            assert_eq!(c.x(), -3);
            assert_eq!(c.y(), -6);
            assert_eq!(c.z(), -9);
        }
        {
            let c = b * a;
            assert_eq!(c.x(), -3);
            assert_eq!(c.y(), -6);
            assert_eq!(c.z(), -9);
        }
    }
    #[test]
    fn mul_assign() {
        let mut a = Vector3Int::new(1, 2, 3);
        a *= -3;

        assert_eq!(a.x(), -3);
        assert_eq!(a.y(), -6);
        assert_eq!(a.z(), -9);
    }

    #[test]
    fn and() {
        let a = Vector3Int::new(1, -2, 0);
        let b = Vector3Int::new(3, -2, -3);

        let c = Vector3Int::and(a, b);
        assert_eq!(c.x(), 1);
        assert_eq!(c.y(), -2);
        assert_eq!(c.z(), 0);
    }
    #[test]
    fn andnot() {
        let a = Vector3Int::new(0, 1, 0);
        let b = Vector3Int::new(1, 3, -3);

        let c = Vector3Int::andnot(a, b);
        assert_eq!(c.x(), 1);
        assert_eq!(c.y(), 2);
        assert_eq!(c.z(), -3);
    }
    #[test]
    fn or() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(2, 1, 3);

        let c = Vector3Int::or(a, b);
        assert_eq!(c.x(), 3);
        assert_eq!(c.y(), 3);
        assert_eq!(c.z(), 3);
    }
    #[test]
    fn xor() {
        let a = Vector3Int::new(1, 0, 3);
        let b = Vector3Int::new(3, 3, 5);

        let c = Vector3Int::xor(a, b);
        assert_eq!(c.x(), 2);
        assert_eq!(c.y(), 3);
        assert_eq!(c.z(), 6);
    }
    #[test]
    fn all() {
        {
            let a = Vector3Int::new(1, 2, 3);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), true);
        }
        {
            let a = Vector3Int::new(0, 2, 3);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3Int::new(1, 0, 3);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3Int::new(1, 2, 0);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        {
            let a = Vector3Int::new(0, 0, 0);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), false);
        }
        unsafe {
            use core::arch::x86_64::*;
            let a = Vector3Int {
                data: _mm_set_epi32(0, 1, 2, 3),
            };
            let b = Vector3Int {
                data: _mm_set_epi32(99, 1, 2, 3),
            };
            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::all(c), true);
        }
    }
    #[test]
    fn any() {
        {
            let a = Vector3Int::new(1, 2, 3);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3Int::new(0, 0, 3);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3Int::new(1, 0, 0);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3Int::new(0, 2, 0);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), true);
        }
        {
            let a = Vector3Int::new(0, 0, 0);
            let b = Vector3Int::new(1, 2, 3);

            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), false);
        }
        unsafe {
            use core::arch::x86_64::*;
            let a = Vector3Int {
                data: _mm_set_epi32(0, 1, 2, 3),
            };
            let b = Vector3Int {
                data: _mm_set_epi32(0, 99, 99, 99),
            };
            let c = Vector3Int::equal(a, b);
            assert_eq!(Vector3Bool::any(c), false);
        }
    }
    #[test]
    fn equal() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::equal(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn not_equal() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::not_equal(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), -1000);
    }
    #[test]
    fn greater_equal() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::greater_equal(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn greater() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::greater(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), -1000);
        assert_eq!(d.z(), 0);
    }
    #[test]
    fn less_equal() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::less_equal(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), -1000);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), -1000);
    }
    #[test]
    fn less() {
        let a = Vector3Int::new(1, 2, 3);
        let b = Vector3Int::new(1, 1, 4);

        let mask = Vector3Int::less(a, b);
        let c = Vector3Int::new(-1000, -1000, -1000);
        let d = Vector3Int::and(c, mask);

        assert_eq!(d.x(), 0);
        assert_eq!(d.y(), 0);
        assert_eq!(d.z(), -1000);
    }

    #[test]
    fn abs() {
        let a = Vector3Int::new(-1, 0, -0);

        let b = Vector3Int::abs(a);

        assert_eq!(b.x(), 1);
        assert_eq!(b.y(), 0);
        assert_eq!(b.z(), 0);
        assert_eq!(a.x(), -1);
        assert_eq!(a.y(), 0);
        assert_eq!(a.z(), -0);
    }
    #[test]
    fn sign() {
        let a = Vector3Int::new(10, -20, 5);
        let b = Vector3Int::new(-1, 1, -0);

        let c = Vector3Int::sign(a, b);

        assert_eq!(c.x(), -10);
        assert_eq!(c.y(), -20);
        assert_eq!(c.z(), 0);
    }

    #[test]
    fn max() {
        {
            let a = Vector3Int::new(-175, 10, 0);
            let b = Vector3Int::new(-200, 20, -10);
            let c = Vector3Int::max(a, b);
            assert_eq!(c.x(), -175);
            assert_eq!(c.y(), 20);
            assert_eq!(c.z(), 0);
        }
    }

    #[test]
    fn min() {
        {
            let a = Vector3Int::new(-175, 10, 0);
            let b = Vector3Int::new(-200, 20, -10);
            let c = Vector3Int::min(a, b);
            assert_eq!(c.x(), -200);
            assert_eq!(c.y(), 10);
            assert_eq!(c.z(), -10);
        }
    }

    #[test]
    fn swizzle() {
        use core::arch::x86_64::*;
        unsafe {
            // make sure this works with garbage in w.
            let a = Vector3Int {
                data: _mm_set_epi32(100, 3, 2, 1),
            };
            {
                let c = a.xxxx();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 1);
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
                let c = a.zzzz();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 3);
                assert_eq!(c.w(), 3);
            }
            {
                let c = a.xxx();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.xxy();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.xxz();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.xyx();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.xyy();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.xyz();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.xzx();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.xzy();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.xzz();
                assert_eq!(c.x(), 1);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 3);
            }

            {
                let c = a.yxx();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.yxy();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.yxz();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.yyx();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.yyy();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.yyz();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.yzx();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.yzy();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.yzz();
                assert_eq!(c.x(), 2);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 3);
            }

            {
                let c = a.zxx();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.zxy();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.zxz();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 1);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.zyx();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.zyy();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.zyz();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 2);
                assert_eq!(c.z(), 3);
            }
            {
                let c = a.zzx();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 1);
            }
            {
                let c = a.zzy();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 2);
            }
            {
                let c = a.zzz();
                assert_eq!(c.x(), 3);
                assert_eq!(c.y(), 3);
                assert_eq!(c.z(), 3);
            }
        }
    }
}
