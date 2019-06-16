use std::arch::x86_64::*;
use crate::_ico_shuffle;
use crate::Vector3;
use crate::Vector4;
use crate::Quaternion;
use crate::DualQuaternion;
use crate::_ico_quat_mul;
use crate::_ico_dp4_ps;
impl DualQuaternion{

	#[inline(always)]
	pub fn new(rotation : Quaternion, translation : Vector3) -> DualQuaternion{
		unsafe{

			let vec = _mm_mul_ps(translation.data, _mm_set_ps(0f32,0.5f32,0.5f32,0.5f32));
			let dual = _ico_quat_mul(vec,rotation.data);
			return DualQuaternion{real : rotation.data, dual : dual };
		}
	}


	#[inline(always)]
	pub fn rotation(&self) ->Quaternion{
		return Quaternion { data : self.real };
	}
	#[inline(always)]
	pub fn translation(&self) ->Vector3{
		let inv_rot = self.rotation().inverse();
		unsafe{
			let p = _ico_quat_mul(self.dual, inv_rot.data);
			return Vector3 { data : _mm_add_ps(p,p)}; //*2
		}
	}

	#[inline(always)]
	pub fn inverse(&self) ->DualQuaternion{
		let real = Quaternion { data : self.real }.inverse();
		let dual = Quaternion { data : self.dual }.inverse();
		return DualQuaternion{real : real.data, dual : dual.data };
	}

	#[inline(always)]
	pub fn reverse(&self) ->DualQuaternion{
		let real = Quaternion { data : self.real }.reverse();
		let dual = Quaternion { data : self.dual }.reverse();
		return DualQuaternion{real : real.data, dual : dual.data };
	}
	#[inline(always)]
	pub fn mul(lhs : DualQuaternion, rhs : DualQuaternion) -> DualQuaternion{	
		unsafe{
			return DualQuaternion{real : _ico_quat_mul(lhs.real, rhs.real), 
			dual : _mm_add_ps( 	_ico_quat_mul(lhs.real, rhs.dual),
								_ico_quat_mul(lhs.dual, rhs.real))};
		}
	}


	#[inline(always)]
	pub fn renormalize(&self) ->DualQuaternion{
		unsafe{
			let sqr_length = _ico_dp4_ps(self.real, self.real);
			let factor = _mm_sqrt_ps(sqr_length);

			return DualQuaternion{real : _mm_div_ps(self.real, factor), 
								 dual : _mm_div_ps(self.dual, factor)};

		}

	}

	#[inline(always)]
	pub fn equals(v1 : DualQuaternion, v2 : DualQuaternion) -> bool{	
		unsafe{
			let d = _mm_cmpeq_ps(v1.real, v2.real);
			let d2 = _mm_cmpeq_ps(v1.dual, v2.dual);
			let d3 = _mm_and_ps(d, d2);
			return (_mm_movemask_ps(d) ) == 15;
		}
	}

	
}

impl PartialEq for DualQuaternion {
    fn eq(&self, other: &DualQuaternion) -> bool {
    	return DualQuaternion::equals(*self, *other);
    }
}