use super::*;


#[cfg(test)]
mod test {
	use crate::float_vector::FloatVector;
	use crate::vector2::Vector2;
	use crate::vector2_bool::Vector2Bool;

	#[test]
    fn new() {

    	let a = Vector2::new(1.0,2.0);

	    assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        
    }
    #[test]
    fn set() {

        let a = Vector2::set(-5.0);

        assert_eq!(a.x(), -5.0);
        assert_eq!(a.y(), -5.0);
    }
    #[test]
    fn zero() {

    	let a = Vector2::zero();

	    assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        
    }
    #[test]
    fn set_x() {

    	let mut a = Vector2::new(1.0,2.0);
    	a.set_x(-11.0);
	    assert_eq!(a.x(), -11.0);
        assert_eq!(a.y(), 2.0);
        
    }
    #[test]
    fn set_y() {

    	let mut a = Vector2::new(1.0,2.0);
    	a.set_y(-11.0);
	    assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), -11.0);
        
    }

    #[test]
    fn add() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);

		let c = a + b;
        assert_eq!(c.x(), 5.0);
        assert_eq!(c.y(), 6.0);
        
    }
    #[test]
    fn add_assign() {

    	let mut a = Vector2::new(1.0,2.0);
    	a += Vector2::set(3.0);

        assert_eq!(a.x(), 4.0);
        assert_eq!(a.y(), 5.0);
        
    }
    #[test]
    fn sub() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);

		let c = a - b;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -2.0);
        
    }
    #[test]
    fn sub_assign() {

    	let mut a = Vector2::new(1.0,2.0);
		a -= Vector2::new(4.0,4.0);

		
        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -2.0);
        
    }
    #[test]
    fn neg() {

    	let a = Vector2::new(1.0,-2.0);
		let c = -a;

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 2.0);
        
    }
    #[test]
    fn component_mul() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);

		let c = Vector2::mul(a,b);
        assert_eq!(c.x(), 4.0);
        assert_eq!(c.y(), 8.0);
        
    }
    #[test]
    fn component_div() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);

		let c = Vector2::div(a,b);
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        
    }
    #[test]
    fn mul_add() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);
		let c = Vector2::new(2.0,-5.0);
		let d = Vector2::mul_add(a,b,c);

        assert_eq!(d.x(), 6.0);
        assert_eq!(d.y(), 3.0);
        
    }
    #[test]
    fn mul_sub() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);
		let c = Vector2::new(2.0,-5.0);
		let d = Vector2::mul_sub(a,b,c);

        assert_eq!(d.x(), 2.0);
        assert_eq!(d.y(), 13.0);
        
    }
    #[test]
    fn neg_mul_add() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);
		let c = Vector2::new(2.0,-5.0);
		let d = Vector2::neg_mul_add(a,b,c);

        assert_eq!(d.x(), -2.0);
        assert_eq!(d.y(), -13.0);
        
    }
    #[test]
    fn neg_mul_sub() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(4.0,4.0);
		let c = Vector2::new(2.0,-5.0);
		let d = Vector2::neg_mul_sub(a,b,c);

        assert_eq!(d.x(), -6.0);
        assert_eq!(d.y(), -3.0);
        
    }
    #[test]
    fn mul() {

    	let a = Vector2::new(1.0,2.0);
		let b = FloatVector::new(-3.0);
		{
		let c = a * b;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -6.0);
    	}
        {
		let c = b * a;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -6.0);
    	}
    }
    #[test]
    fn mul_assign() {

    	let mut a = Vector2::new(1.0,2.0);
		a *= -3.0;

        assert_eq!(a.x(), -3.0);
        assert_eq!(a.y(), -6.0);
    	
    }
    #[test]
    fn div() {

    	let a = Vector2::new(1.0,2.0);
		let b = 4.0;

		let c = a/b;
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        
    }
    #[test]
    fn div_inv() {

    	let a = Vector2::new(1.0,2.0);
		let b = FloatVector::new(1.0);

		let c = b / a;
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), 0.5);
        
    }
    #[test]
    fn div_assign() {

    	let mut a = Vector2::new(1.0,2.0);
		a /= 4.0;

        assert_eq!(a.x(), 0.25);
        assert_eq!(a.y(), 0.5);
        
    }
    #[test]
    fn and() {

    	let a = Vector2::new(-1.0,-2.0);
		let b = Vector2::new(1.0,-2.0);

		let c = Vector2::and(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        
    }
    #[test]
    fn andnot() {

    	let a = Vector2::new(0.0,0.0);
		let b = Vector2::new(1.0,-2.0);

		let c = Vector2::andnot(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        
    }
    #[test]
    fn or() {

    	let a = Vector2::new(-1.0,-2.0);
		let b = Vector2::new(1.0,-2.0);

		let c = Vector2::or(a,b);
        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), -2.0);
        
    }
    #[test]
    fn xor() {

    	let a = Vector2::new(-1.0,-2.0);
		let b = Vector2::new(-0.0,0.0);

		let c = Vector2::xor(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        
    }
    #[test]
    fn all() {
    	{
	    	let a = Vector2::new(1.0,2.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::all(c), true);
		}	
		{
	    	let a = Vector2::new(0.0,2.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::all(c), false);
		}	
		{
	    	let a = Vector2::new(1.0,0.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::all(c), false);
		}

		{
	    	let a = Vector2::new(0.0,0.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::all(c), false);
		}	
		unsafe{
			use core::arch::x86_64::*;
			let a = Vector2 { data : _mm_set_ps(0.0, 1.0,2.0, 0.0)};
			let b = Vector2 { data : _mm_set_ps(99.0, 1.0,2.0, 0.0)};
			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::all(c), true);
		}

    }
    #[test]
    fn any() {
    	{
	    	let a = Vector2::new(1.0,2.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::any(c), true);
		}		
		{
	    	let a = Vector2::new(1.0,0.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::any(c), true);
		}
		{
	    	let a = Vector2::new(0.0,2.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::any(c), true);
		}	
		{
	    	let a = Vector2::new(0.0,0.0);
			let b = Vector2::new(1.0,2.0);

			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::any(c), false);
		}	
		unsafe{
			use core::arch::x86_64::*;
			let a = Vector2 { data : _mm_set_ps(0.0, 1.0,2.0,3.0)};
			let b = Vector2 { data : _mm_set_ps(0.0, 99.0,99.0,99.0)};
			let c = Vector2::equal(a,b);
			assert_eq!(Vector2Bool::any(c), false);
		}

    }
    #[test]
    fn equal() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::equal(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        
    }
    #[test]
    fn not_equal() {

    	let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::not_equal(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        
    }
    #[test]
	fn greater_equal() {
		let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::greater_equal(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
	}
	#[test]
	fn greater() {
		let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::greater(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
	}
	#[test]
	fn less_equal() {
		let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::less_equal(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
	}
	#[test]
	fn less() {
		let a = Vector2::new(1.0,2.0);
		let b = Vector2::new(1.0,1.0);

		let mask = Vector2::less(a,b);
		let c = Vector2::new(-1000.0,-1000.0);
		let d = Vector2::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
	}


	#[test]
	fn abs(){	
		let a = Vector2::new(-1.0,0.0);

		let b = Vector2::abs(a);

        assert_eq!(b.x(), 1.0);
        assert_eq!(b.y(), 0.0);

	}
	#[test]
	fn copysign(){	
		let a = Vector2::new(-1.0,0.0);
		let b = Vector2::new(10.0,-20.0);
		let c = Vector2::copysign(b,a);

        assert_eq!(c.x(), -10.0);
        assert_eq!(c.y(), 20.0);

	}
	
	#[test]
	fn floor(){	
		{
		let a = Vector2::new(-1.1,0.7);
		let c = Vector2::floor(a);

        assert_eq!(c.x(), -2.0);
        assert_eq!(c.y(), 0.0);
        
    	}	
    	{
		let a = Vector2::floor(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		let b = Vector2::new(0.0,0.0);
		let c = Vector2::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

    	}
	}


	#[test]
	fn ceil(){	
		{
		let a = Vector2::new(-1.1,-0.7);
		let c = Vector2::ceil(a);

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
    	}
        {
		let a = Vector2::ceil(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		let b = Vector2::new(0.0,0.0);
		let c = Vector2::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

    	}
	}

	#[test]
	fn round(){	
		{
		let a = Vector2::new(1.5,-0.5);
		let c = Vector2::round(a);

        assert_eq!(c.x(), 2.0);
        assert_eq!(c.y(), 0.0);

    	}
    	{
		let a = Vector2::round(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		let b = Vector2::new(0.0,0.0);
		let c = Vector2::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

    	}
	}

	
	#[test]
	fn floor_to_int(){	
		{
		let a = Vector2::new(-1.1,0.7);
		let c = Vector2::floor_to_int(a);

        assert_eq!(c.x(), -2);
        assert_eq!(c.y(), 0);
        
    	}	
    	{
		let a = Vector2::floor_to_int(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		
        assert_eq!(a.x(), -2147483648);
        assert_eq!(a.y(), -2147483648);

    	}
	}

	#[test]
	fn ceil_to_int(){	
		{
		let a = Vector2::new(-1.1,0.7);
		let c = Vector2::ceil_to_int(a);

        assert_eq!(c.x(), -1);
        assert_eq!(c.y(), 1);
        
    	}	
    	{
		let a = Vector2::ceil_to_int(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		
        //assert_eq!(a.x(), 2.0 * 2147483647.0);
        //assert_eq!(a.y(), -2.0 * 2147483648.0);

        //assert_eq!(c.z(), -1.0);
    	}
	}

	#[test]
	fn truncate(){	
		{
		let a = Vector2::new(-1.6,0.7);
		let c = Vector2::truncate(a);

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
        
    	}	
    	{
		let a = Vector2::truncate(Vector2::new(2.0 * 2147483647.0,-2.0 * 2147483648.0));
		let b = Vector2::new(0.0,0.0);
		let c = Vector2::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

    	}
	}
	
	#[test]
	fn frac(){	
		{
		let a = Vector2::frac(Vector2::new(-1.75,0.1));

        assert_eq!(a.x(), 0.25);
        assert_eq!(a.y(), 0.1);
        
    	}	
	}
	
	#[test]
	fn sqrt(){	
		{
		let a = Vector2::sqrt(Vector2::new(4.0,9.0));

        assert_eq!(a.x(), 2.0);
        assert_eq!(a.y(), 3.0);
        
    	}	
	}

	#[test]
	fn max(){	
		{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::max(a,b);
        assert_eq!(c.x(), -1.75);
        assert_eq!(c.y(), 0.2);
        
    	}	
	}

	#[test]
	fn min(){	
		{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::min(a,b);
        assert_eq!(c.x(), -2.0);
        assert_eq!(c.y(), 0.1);
        
    	}	
	}
	
	#[test]
	fn dot(){	
		use core::arch::x86_64::*;
		unsafe{
			// make sure this works with garbage in zw.
		let a = Vector2{data : _mm_set_ps(100.0, -1.0, 0.0, 0.0)};
		let b = Vector2{data : _mm_set_ps(200.0, 3.0, 0.0, -1.0)};
		let c = Vector2::from(Vector2::dot(a,b));

		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 0.0);
    	}
    	unsafe{
			// make sure this works with garbage in zw.
		let a = Vector2{data : _mm_set_ps(100.0, -1.0, -1.0, 5.0)};
		let b = Vector2{data : _mm_set_ps(200.0, 1.0, 1.0, 7.0)};
		let c = Vector2::from(Vector2::dot(a,b));

		assert_eq!(c.x(), 34.0);
        assert_eq!(c.y(), 34.0);
    	}
	}

	
	#[test]
	fn lerp(){	
		{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::lerp(a,b, 0.0);
		assert_eq!(c.x(), -1.75);
        assert_eq!(c.y(), 0.1);
    	}
    	{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::lerp(a,b, 1.0);
		assert_eq!(c.x(), -2.0);
        assert_eq!(c.y(), 0.2);
    	}
    	{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::lerp(a,b, 0.5);
		assert_eq!(c.x(), -1.875);
        assert_eq!(c.y(), 0.15);
    	}
    	{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::lerp(a,b, 2.0);
		assert_eq!(c.x(), -2.25);
        assert_eq!(c.y(), 0.3);
    	}
    	{
		let a = Vector2::new(-1.75,0.1);
		let b = Vector2::new(-2.0,0.2);
		let c = Vector2::lerp(a,b, -1.0);
		assert_eq!(c.x(), -1.5);
        assert_eq!(c.y(), 0.0);
    	}
	}
	

	#[test]
	fn normalize(){	
		{
		let a = Vector2::new(-2.0,0.0);
		let c = Vector2::normalize(a);
		assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
    	}
    	{
    	//SHOULD NOT BE NAN
		let a = Vector2::new( core::f32::MAX,1.0);
		let c = Vector2::normalize(a);
		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 0.0);
    	}
    	{
    	//SHOULD NOT BE NAN
		let a = Vector2::new(0.0,0.0);
		let c = Vector2::normalize(a);
		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 0.0);
    	}
    	{
    	//SHOULD NOT BE NAN
		let a = Vector2::new(core::f32::NAN,0.0);
		let c = Vector2::normalize(a);
		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 0.0);
    	}
    	
	}

	

	

	#[test]
	fn sqr_magnitude()  {
		
		use core::arch::x86_64::*;
		unsafe{
			// make sure this works with garbage in zw.
		let a = Vector2{data : _mm_set_ps(100.0, -3.0, -3.0, 0.0)};

		let c = Vector2::sqr_magnitude(a);

		assert_eq!(c, 9.0);
    	}	
	}
	#[test]
	fn magnitude() {
		use core::arch::x86_64::*;
		unsafe{
			// make sure this works with garbage in zw.
		let a = Vector2{data : _mm_set_ps(100.0, -3.0, -3.0, 0.0)};

		let c = Vector2::magnitude(a);

		assert_eq!(c, 3.0);
    	}	
	
	}
	#[test]
    fn rotate() {
        use core::arch::x86_64::*;
        unsafe{
            // make sure this works with garbage in zw.
        let a = Vector2{data : _mm_set_ps(100.0, -3.0, 1.0, 0.0)};

        let c = Vector2::rotate(a, 0.5* core::f64::consts::PI);

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
        }   
        unsafe{
            // make sure this works with garbage in zw.
        let a = Vector2{data : _mm_set_ps(100.0, -3.0, 0.0, 1.0)};

        let c = Vector2::rotate(a, 0.5* core::f64::consts::PI);

        assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 1.0);
        } 
        unsafe{
            // make sure this works with garbage in zw.
        let a = Vector2{data : _mm_set_ps(100.0, -3.0, 0.0, -1.0)};

        let c = Vector2::rotate(a, 0.5* core::f64::consts::PI);

        assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), -1.0);
        } 
        unsafe{
            // make sure this works with garbage in zw.
        let a = Vector2{data : _mm_set_ps(100.0, -3.0, 1.0, 0.0)};

        let c = Vector2::rotate(a, core::f64::consts::PI);

        assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), -1.0);
        }  
    }

	#[test]
	fn swizzle() {
		use core::arch::x86_64::*;
		unsafe{
			// make sure this works with garbage in zw.
		let a = Vector2{data : _mm_set_ps(100.0, 3.0, 2.0, 1.0)};
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
			let c = a.xx();
			assert_eq!(c.x(), 1.0);
			assert_eq!(c.y(), 1.0);
		}
        {
            let c = a.xy();
            assert_eq!(c.x(), 1.0);
            assert_eq!(c.y(), 2.0);
        }
        {
            let c = a.yx();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 1.0);
        }
        {
            let c = a.yy();
            assert_eq!(c.x(), 2.0);
            assert_eq!(c.y(), 2.0);
        }

    	}	
	
	}
}