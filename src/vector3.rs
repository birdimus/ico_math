use crate::Vector2;
use crate::Vector3;
use crate::Vector4;
use crate::Vector3Int;

use std::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::_ico_abs_ps;
use crate::_ico_cross_ps;
use crate::_ico_copysign_ps;

impl Vector3{
	/// Returns a new Vector3
	#[inline(always)]
	pub fn new(x : f32, y : f32, z : f32) -> Vector3{
		unsafe{
			Vector3{data : _mm_set_ps(0.0f32, z, y, x)}
		}
	}
	#[inline(always)]
	pub fn zero() -> Vector3 {
		unsafe{
			Vector3 { data : _mm_setzero_ps() }
		}
	}

	#[inline(always)]
	pub fn x(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(1, 1, 1, 1)));
		}	
	}

	#[inline(always)]
	pub fn z(self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(2, 2, 2, 2)));
		}	
	}

	#[inline(always)]
	pub fn set_x(&mut self, value : f32) {
		unsafe{
			self.data = _mm_move_ss(self.data, _mm_set_ss(value));
		}	
	}

	#[inline(always)]
	pub fn set_y(&mut self, value : f32) {
		unsafe{
			let v1 = _mm_move_ss(_mm_set1_ps(value), self.data);
			self.data = _mm_shuffle_ps(v1,self.data, _ico_shuffle(3, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn set_z(&mut self, value : f32) {
		unsafe{
			self.data = _mm_shuffle_ps(self.data, _mm_set1_ps(value), _ico_shuffle(3, 2, 1, 0));
		}	
	}

	#[inline(always)]
	pub fn add(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_add_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn sub(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_sub_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_mul(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_mul_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_div(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_div_ps(v1.data, v2.data)}
		}
	}

	/// Multiply v1 and v2, and add v3.
	#[inline(always)]
	pub fn fmadd(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fmsub(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fmsub_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmadd(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fnmadd_ps(v1.data, v2.data, v3.data)}
		}
	}
	#[inline(always)]
	pub fn fnmsub(v1 : Vector3, v2 : Vector3, v3 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_fnmsub_ps(v1.data, v2.data, v3.data)}
		}
	}

	#[inline(always)]
	pub fn scale(v1 : Vector3, scalar : f32) -> Vector3{	
		unsafe{
			Vector3{data : _mm_mul_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn div(v1 : Vector3, scalar : f32) -> Vector3{	
		unsafe{
			Vector3{data : _mm_div_ps(v1.data, _mm_set1_ps(scalar))}
		}
	}

	#[inline(always)]
	pub fn and(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_and_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn or(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_or_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn andnot(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_andnot_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn xor(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_xor_ps(v1.data, v2.data)}
		}
	}
	
	#[inline(always)]
	pub fn all(v1 : Vector3) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 7) == 7;
		}
	}

	#[inline(always)]
	pub fn any(v1 : Vector3) -> bool{	
		unsafe{
			return (_mm_movemask_ps(v1.data) & 7) != 0;
		}
	}

	#[inline(always)]
	pub fn equals(v1 : Vector3, v2 : Vector3) -> bool{	
		let d = Vector3::component_equal(v1,v2);
		return Vector3::all(d);
	}

	#[inline(always)]
	pub fn component_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpeq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_not_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpneq_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn component_greater_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpge_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_greater(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmpgt_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less_equal(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmple_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn component_less(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_cmplt_ps(v1.data, v2.data)}
		}
	}


	#[inline(always)]
	/// Right handed cross product
	pub fn cross(lhs : Vector3, rhs : Vector3) -> Vector3{	
		unsafe{
			return Vector3{data :  _ico_cross_ps(lhs.data, rhs.data)};
		}
	}

	#[inline(always)]
	pub fn abs(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_abs_ps(v1.data)}
		}
	}
	#[inline(always)]
	pub fn copysign(magnitude : Vector3, sign : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _ico_copysign_ps(magnitude.data, sign.data)}
		}
	}
	#[inline(always)]
	/// Floor function.  Returns signed 0 when applicable.
	pub fn floor(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_floor_ps(v1.data)}
				//_ico_floor_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Ceil function.  Returns signed 0 when applicable.
	pub fn ceil(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_ceil_ps(v1.data)}
				//_ico_ceil_ps(v1.data)}
		}
	}

	#[inline(always)]
	/// Round to nearest even function. Returns signed 0 when applicable.
	pub fn round(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : 
			_mm_round_ps(v1.data, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC)}
			// _ico_round_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn floor_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _mm_cvttps_epi32(_mm_floor_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn ceil_to_int(v1 : Vector3) -> Vector3Int{	
		unsafe{
			Vector3Int{data :  _mm_cvttps_epi32(_mm_ceil_ps(v1.data))}
		}
	}

	#[inline(always)]
	pub fn truncate(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  
				_mm_round_ps(v1.data, _MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC)}
				//_ico_truncate_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn frac(v1 : Vector3) -> Vector3{	
		return Vector3::sub(v1, Vector3::floor(v1));
	}

	#[inline(always)]
	pub fn sqrt(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[inline(always)]
	pub fn max(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[inline(always)]
	pub fn min(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_min_ps(v1.data, v2.data)}
		}
	}
	#[inline(always)]
	pub fn dot(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{

			let tmp0 = _mm_mul_ps(v1.data, v2.data); //xyzw
		    let tmp1 = _mm_castsi128_ps(_mm_slli_si128 (_mm_castps_si128(tmp0), 4)); //0xyz
		    let tmp2 = _mm_add_ps(tmp0 , tmp1);//x xy, yz, wz
		    let tmp3 = _mm_moveldup_ps(tmp2); // x x yz yz
			Vector3{data : _mm_add_ps(tmp3, _mm_shuffle_ps(tmp3,tmp3, _ico_shuffle(0, 1, 2, 3)))}
		}
	}
	#[inline(always)]
	pub fn lerp(v1 : Vector3, v2 : Vector3, t : f32) -> Vector3{	
		unsafe{
			let t_val = _mm_set1_ps(t);
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector3{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}

	#[inline(always)]
	pub fn renormalize(v1 : Vector3) -> Vector3{	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		return Vector3::component_div(v1, length);
	}

	#[inline(always)]
	pub fn normalize(v1 : Vector3) -> Vector3{	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::component_div(v1, length);
		let mask = Vector3::component_less(Vector3::abs(norm), Vector3::from( std::f32::INFINITY));
		return Vector3::and(norm, mask);
	}

	#[inline(always)]
	pub fn normalize_length(v1 : Vector3) -> (Vector3, f32){	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::component_div(v1, length);
		let mask = Vector3::component_less(Vector3::abs(norm),Vector3::from( std::f32::INFINITY));
		return (Vector3::and(norm, mask), length.x());
	}

	

	#[inline(always)]
	pub fn sqr_magnitude(self) -> f32 {
		return Vector3::dot(self, self).x();	
	}
	#[inline(always)]
	pub fn magnitude(self) -> f32 {
		return Vector3::sqrt(Vector3::dot(self, self)).x();	
	}

	
}

impl From<f32> for Vector3 {
    fn from(val : f32) -> Vector3 {
       unsafe{
			Vector3{data : _mm_set1_ps(val)}
		}
    }
}
impl From<Vector2> for Vector3 {
    fn from(v : Vector2) -> Vector3 {
        Vector3 { data : v.data }
    }
}
impl From<Vector4> for Vector3 {
    fn from(v : Vector4) -> Vector3 {
        Vector3 { data : v.data }
    }
}
impl From<Vector3Int> for Vector3 {
    fn from(v : Vector3Int) -> Vector3 {
    	unsafe{
        	return Vector3 { data : _mm_cvtepi32_ps(v.data) };
        }
    }
}

impl std::ops::Add for Vector3{
	type Output = Vector3;
	#[inline]
	fn add(self, _rhs: Vector3) -> Vector3{
		Vector3::add(self, _rhs)
	}
}

impl std::ops::Sub for Vector3{
	type Output = Vector3;
	#[inline]
	fn sub(self, _rhs: Vector3) -> Vector3{
		Vector3::sub(self, _rhs)
	}
}

impl std::ops::Mul<f32> for Vector3{
	type Output = Vector3;
	#[inline]
	fn mul(self, _rhs: f32) -> Vector3{
		Vector3::scale(self, _rhs)
	}
}

impl std::ops::Mul<Vector3> for f32{
	type Output = Vector3;
	#[inline]
	fn mul(self, _rhs: Vector3) -> Vector3{
		Vector3::scale(_rhs, self)
	}
}

impl std::ops::Div<f32> for Vector3{
	type Output = Vector3;
	#[inline]
	fn div(self, _rhs: f32) -> Vector3{
		Vector3::div(self, _rhs)
	}
}
	
impl PartialEq for Vector3 {
    fn eq(&self, other: &Vector3) -> bool {
    	return Vector3::equals(*self, *other);
    }
}

#[cfg(test)]
mod tests {
	use crate::Vector3;
	#[test]
    fn new() {

    	let a = Vector3::new(1.0,2.0,3.0);

	    assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        
    }
    #[test]
    fn zero() {

    	let a = Vector3::zero();

	    assert_eq!(a.x(), 0.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), 0.0);
        
    }
    #[test]
    fn set_x() {

    	let mut a = Vector3::new(1.0,2.0,3.0);
    	a.set_x(-11.0);
	    assert_eq!(a.x(), -11.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), 3.0);
        
    }
    #[test]
    fn set_y() {

    	let mut a = Vector3::new(1.0,2.0,3.0);
    	a.set_y(-11.0);
	    assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), -11.0);
        assert_eq!(a.z(), 3.0);
        
    }
    #[test]
    fn set_z() {

    	let mut a = Vector3::new(1.0,2.0,3.0);
    	a.set_z(-11.0);
	    assert_eq!(a.x(), 1.0);
        assert_eq!(a.y(), 2.0);
        assert_eq!(a.z(), -11.0);
        
    }
    #[test]
    fn add() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);

		let c = a + b;
        assert_eq!(c.x(), 5.0);
        assert_eq!(c.y(), 6.0);
        assert_eq!(c.z(), 9.0);
        
    }
    #[test]
    fn sub() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);

		let c = a - b;
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        
    }
    #[test]
    fn component_mul() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);

		let c = Vector3::component_mul(a,b);
        assert_eq!(c.x(), 4.0);
        assert_eq!(c.y(), 8.0);
        assert_eq!(c.z(), 18.0);
        
    }
    #[test]
    fn component_div() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);

		let c = Vector3::component_div(a,b);
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.5);
        
    }
    #[test]
    fn fmadd() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);
		let c = Vector3::new(2.0,-5.0,10.0);
		let d = Vector3::fmadd(a,b,c);

        assert_eq!(d.x(), 6.0);
        assert_eq!(d.y(), 3.0);
        assert_eq!(d.z(), 28.0);
        
    }
    #[test]
    fn fmsub() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);
		let c = Vector3::new(2.0,-5.0,10.0);
		let d = Vector3::fmsub(a,b,c);

        assert_eq!(d.x(), 2.0);
        assert_eq!(d.y(), 13.0);
        assert_eq!(d.z(), 8.0);
        
    }
    #[test]
    fn fnmadd() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);
		let c = Vector3::new(2.0,-5.0,10.0);
		let d = Vector3::fnmadd(a,b,c);

        assert_eq!(d.x(), -2.0);
        assert_eq!(d.y(), -13.0);
        assert_eq!(d.z(), -8.0);
        
    }
    #[test]
    fn fnmsub() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(4.0,4.0,6.0);
		let c = Vector3::new(2.0,-5.0,10.0);
		let d = Vector3::fnmsub(a,b,c);

        assert_eq!(d.x(), -6.0);
        assert_eq!(d.y(), -3.0);
        assert_eq!(d.z(), -28.0);
        
    }
    #[test]
    fn scale() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = -3.0;

		let c = Vector3::scale(a,b);
        assert_eq!(c.x(), -3.0);
        assert_eq!(c.y(), -6.0);
        assert_eq!(c.z(), -9.0);
        
    }
    #[test]
    fn div() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = 4.0;

		let c = Vector3::div(a,b);
        assert_eq!(c.x(), 0.25);
        assert_eq!(c.y(), 0.5);
        assert_eq!(c.z(), 0.75);
        
    }
    #[test]
    fn and() {

    	let a = Vector3::new(-1.0,-2.0,3.0);
		let b = Vector3::new(1.0,-2.0,-3.0);

		let c = Vector3::and(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
        
    }
    #[test]
    fn andnot() {

    	let a = Vector3::new(0.0,0.0,-0.0);
		let b = Vector3::new(1.0,-2.0,-3.0);

		let c = Vector3::andnot(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), 3.0);
        
    }
    #[test]
    fn or() {

    	let a = Vector3::new(-1.0,-2.0,3.0);
		let b = Vector3::new(1.0,-2.0,-3.0);

		let c = Vector3::or(a,b);
        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        
    }
    #[test]
    fn xor() {

    	let a = Vector3::new(-1.0,-2.0,3.0);
		let b = Vector3::new(-0.0,0.0,-0.0);

		let c = Vector3::xor(a,b);
        assert_eq!(c.x(), 1.0);
        assert_eq!(c.y(), -2.0);
        assert_eq!(c.z(), -3.0);
        
    }
    #[test]
    fn all() {
    	{
	    	let a = Vector3::new(1.0,2.0,3.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), true);
		}	
		{
	    	let a = Vector3::new(0.0,2.0,3.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), false);
		}	
		{
	    	let a = Vector3::new(1.0,0.0,3.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), false);
		}
		{
	    	let a = Vector3::new(1.0,2.0,0.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), false);
		}	
		{
	    	let a = Vector3::new(0.0,0.0,0.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), false);
		}	
		unsafe{
			use std::arch::x86_64::*;
			let a = Vector3 { data : _mm_set_ps(0.0, 1.0,2.0,3.0)};
			let b = Vector3 { data : _mm_set_ps(99.0, 1.0,2.0,3.0)};
			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::all(c), true);
		}

    }
    #[test]
    fn any() {
    	{
	    	let a = Vector3::new(1.0,2.0,3.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), true);
		}	
		{
	    	let a = Vector3::new(0.0,0.0,3.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), true);
		}	
		{
	    	let a = Vector3::new(1.0,0.0,0.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), true);
		}
		{
	    	let a = Vector3::new(0.0,2.0,0.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), true);
		}	
		{
	    	let a = Vector3::new(0.0,0.0,0.0);
			let b = Vector3::new(1.0,2.0,3.0);

			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), false);
		}	
		unsafe{
			use std::arch::x86_64::*;
			let a = Vector3 { data : _mm_set_ps(0.0, 1.0,2.0,3.0)};
			let b = Vector3 { data : _mm_set_ps(0.0, 99.0,99.0,99.0)};
			let c = Vector3::component_equal(a,b);
			assert_eq!(Vector3::any(c), false);
		}

    }
    #[test]
    fn component_equal() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_equal(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), 0.0);
        
    }
    #[test]
    fn component_not_equal() {

    	let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_not_equal(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), -1000.0);
        
    }
    #[test]
	fn component_greater_equal() {
		let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_greater_equal(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
	}
	#[test]
	fn component_greater() {
		let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_greater(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), -1000.0);
        assert_eq!(d.z(), 0.0);
	}
	#[test]
	fn component_less_equal() {
		let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_less_equal(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), -1000.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
	}
	#[test]
	fn component_less() {
		let a = Vector3::new(1.0,2.0,3.0);
		let b = Vector3::new(1.0,1.0,4.0);

		let mask = Vector3::component_less(a,b);
		let c = Vector3::new(-1000.0,-1000.0,-1000.0);
		let d = Vector3::and(c, mask);

        assert_eq!(d.x(), 0.0);
        assert_eq!(d.y(), 0.0);
        assert_eq!(d.z(), -1000.0);
	}
	#[test]
	fn cross(){	

		{
		let a = Vector3::new(0.0,0.0,-1.0);
		let b = Vector3::new(-1.0,0.0,0.0);
		let c = Vector3::cross(a,b);

		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), 1.0);
        assert_eq!(c.z(), 0.0);
    	}
    	{
		let a = Vector3::new(0.0,0.0,-1.0);
		let b = Vector3::new(-1.0,0.0,0.0);
		let c = Vector3::cross(b,a);

		assert_eq!(c.x(), 0.0);
        assert_eq!(c.y(), -1.0);
        assert_eq!(c.z(), 0.0);
    	}
	}

	#[test]
	fn abs(){	
		let a = Vector3::new(-1.0,0.0,-0.0);

		let b = Vector3::abs(a);

        assert_eq!(b.x(), 1.0);
        assert_eq!(b.y(), 0.0);
        assert_eq!(b.z(), 0.0);
        assert_eq!(a.x(), -1.0);
        assert_eq!(a.y(), 0.0);
        assert_eq!(a.z(), -0.0);
	}
	#[test]
	fn copysign(){	
		let a = Vector3::new(-1.0,0.0,-0.0);
		let b = Vector3::new(10.0,-20.0, 5.0);
		let c = Vector3::copysign(b,a);

        assert_eq!(c.x(), -10.0);
        assert_eq!(c.y(), 20.0);
        assert_eq!(c.z(), -5.0);

	}
	
	#[test]
	fn floor(){	
		{
		let a = Vector3::new(-1.1,0.1,0.7);
		let c = Vector3::floor(a);

        assert_eq!(c.x(), -2.0);
        assert_eq!(c.y(), 0.0);
        assert_eq!(c.z(), 0.0);
        
    	}	
    	{
		let a = Vector3::floor(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.0));
		let b = Vector3::new(0.0,0.0,1.0);
		let c = Vector3::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

        assert_eq!(c.z(), -1.0);
    	}
	}


	#[test]
	fn ceil(){	
		{
		let a = Vector3::new(-1.1,-0.7,0.7);
		let c = Vector3::ceil(a);

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
        assert_eq!(c.z(), 1.0);
    	}
        {
		let a = Vector3::ceil(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.0));
		let b = Vector3::new(0.0,0.0,1.0);
		let c = Vector3::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

        assert_eq!(c.z(), -1.0);
    	}
	}

	#[test]
	fn round(){	
		{
		let a = Vector3::new(1.5,-0.5,0.5);
		let c = Vector3::round(a);

        assert_eq!(c.x(), 2.0);
        assert_eq!(c.y(), 0.0);
        assert_eq!(c.z(), 0.0);

    	}
    	{
		let a = Vector3::round(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.0));
		let b = Vector3::new(0.0,0.0,1.0);
		let c = Vector3::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

        assert_eq!(c.z(), -1.0);
    	}
	}

	
	#[test]
	fn floor_to_int(){	
		{
		let a = Vector3::new(-1.1,0.1,0.7);
		let c = Vector3::floor_to_int(a);

        assert_eq!(c.x(), -2);
        assert_eq!(c.y(), 0);
        assert_eq!(c.z(), 0);
        
    	}	
    	{
		let a = Vector3::floor_to_int(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.0));
		
        assert_eq!(a.x(), -2147483648);
        assert_eq!(a.y(), -2147483648);

        assert_eq!(a.z(), 0);
    	}
	}

	#[test]
	fn ceil_to_int(){	
		{
		let a = Vector3::new(-1.1,0.1,0.7);
		let c = Vector3::ceil_to_int(a);

        assert_eq!(c.x(), -1);
        assert_eq!(c.y(), 1);
        assert_eq!(c.z(), 1);
        
    	}	
    	{
		let a = Vector3::ceil_to_int(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.0));
		
        //assert_eq!(a.x(), 2.0 * 2147483647.0);
        //assert_eq!(a.y(), -2.0 * 2147483648.0);

        //assert_eq!(c.z(), -1.0);
    	}
	}

	#[test]
	fn truncate(){	
		{
		let a = Vector3::new(-1.6,0.1,0.7);
		let c = Vector3::truncate(a);

        assert_eq!(c.x(), -1.0);
        assert_eq!(c.y(), 0.0);
        assert_eq!(c.z(), 0.0);
        
    	}	
    	{
		let a = Vector3::truncate(Vector3::new(2.0 * 2147483647.0,-2.0 * 2147483648.0, -0.1));
		let b = Vector3::new(0.0,0.0,1.0);
		let c = Vector3::copysign(b,a);
        assert_eq!(a.x(), 2.0 * 2147483647.0);
        assert_eq!(a.y(), -2.0 * 2147483648.0);

        assert_eq!(c.z(), -1.0);
    	}
	}
	
	#[test]
	fn frac(){	
		{
		let a = Vector3::frac(Vector3::new(-1.75,0.1,0.0));

        assert_eq!(a.x(), 0.25);
        assert_eq!(a.y(), 0.1);
        assert_eq!(a.z(), 0.0);
        
    	}	
	}
	/*
	#[test]
	fn sqrt(v1 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data :  _mm_sqrt_ps(v1.data)}
		}
	}

	#[test]
	fn max(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_max_ps(v1.data, v2.data)}
		}
	}

	#[test]
	fn min(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_min_ps(v1.data, v2.data)}
		}
	}
	#[test]
	fn dot(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{

			let tmp0 = _mm_mul_ps(v1.data, v2.data); //xyzw
		    let tmp1 = _mm_castsi128_ps(_mm_slli_si128 (_mm_castps_si128(tmp0), 4)); //0xyz
		    let tmp2 = _mm_add_ps(tmp0 , tmp1);//x xy, yz, wz
		    let tmp3 = _mm_moveldup_ps(tmp2); // x x yz yz
			Vector3{data : _mm_add_ps(tmp3, _mm_shuffle_ps(tmp3,tmp3, _ico_shuffle(0, 1, 2, 3)))}
		}
	}
	#[test]
	fn lerp(v1 : Vector3, v2 : Vector3, t : f32) -> Vector3{	
		unsafe{
			let t_val = _mm_set1_ps(t);
			let tmp = _mm_fnmadd_ps(v1.data, t_val, v1.data); //a - (a*t)
			Vector3{data : _mm_fmadd_ps(v2.data, t_val, tmp)} //b * t + a
		}
	}

	#[test]
	fn renormalize(v1 : Vector3) -> Vector3{	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		return Vector3::component_div(v1, length);
	}

	#[test]
	fn normalize(v1 : Vector3) -> Vector3{	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::component_div(v1, length);
		let mask = Vector3::component_less(Vector3::abs(norm), Vector3::from( std::f32::INFINITY));
		return Vector3::and(norm, mask);
	}

	#[test]
	fn normalize_length(v1 : Vector3) -> (Vector3, f32){	
		let length = Vector3::sqrt(Vector3::dot(v1,v1));
		let norm = Vector3::component_div(v1, length);
		let mask = Vector3::component_less(Vector3::abs(norm),Vector3::from( std::f32::INFINITY));
		return (Vector3::and(norm, mask), length.x());
	}

	

	#[test]
	fn sqr_magnitude(self) -> f32 {
		return Vector3::dot(self, self).x();	
	}
	#[test]
	fn magnitude(self) -> f32 {
		return Vector3::sqrt(Vector3::dot(self, self)).x();	
	}
	*/
}