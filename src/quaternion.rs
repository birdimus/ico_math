use std::arch::x86_64::*;
use crate::Vector3;
use crate::Quaternion;
use crate::RawVec;
use crate::_ico_shuffle;

impl Quaternion{

	#[inline(always)]
	pub fn identity() -> Quaternion {
		unsafe{
			Quaternion { data : _mm_set_ps(1.0f32, 0.0f32, 0.0f32, 0.0f32) }
		}
	}


	#[inline(always)]
	pub fn mul(quat : Quaternion, vec : Vector3) -> Vector3{	
		unsafe{
			// t = 2.0 * cross(q.xyz, vec)
			let qv = Vector3{data: quat.data};
			let t = Vector3::cross(qv,vec);
			let t2 = Vector3::add(t,t);
		
			// v' = v + q.w * t + cross(q.xyz, t)
			let w = Vector3{data:_mm_shuffle_ps(quat.data, quat.data, _ico_shuffle(3, 3, 3, 3))};
			let tmp = Vector3::fmadd(w, t2, vec);
	        
			return Vector3::add(tmp, Vector3::cross(qv, t2));
		}
	}

	#[inline(always)]
	pub fn angle_axis(axis : Vector3, radians : f32) -> Quaternion{	
		let half_angle = radians *0.5f32;

		let sin_cos = half_angle.sin_cos();
		unsafe{
			let mut q = Quaternion{data : _mm_mul_ps(axis.data, _mm_set1_ps(sin_cos.0))};
		
		q.set_w(sin_cos.1);
		return q;
		}
	}
	#[inline(always)]
	pub fn store(&self,  dst : &mut RawVec){	
		let x : *mut f32 = &mut (dst.data[0]) as *mut f32;
		unsafe{
			_mm_store_ps(x, self.data);
		}
	}

	//#[inline(always)]
	//pub fn to_angle_axis(&self) -> (Vector3, f32){	

	//}

	#[inline(always)]
	pub fn x(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(self.data);
		}	
	}

	#[inline(always)]
	pub fn y(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(1, 1, 1, 1)));
		}	
	}

	#[inline(always)]
	pub fn z(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(2, 2, 2, 2)));
		}	
	}

	#[inline(always)]
	pub fn w(&self) -> f32 {
		unsafe{
			return _mm_cvtss_f32(_mm_shuffle_ps(self.data, self.data, _ico_shuffle(3, 3, 3, 3)));
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
	pub fn set_w(&mut self, value : f32) {
		unsafe{
			let v1 = _mm_move_ss(self.data, _mm_set_ss(value));
			self.data = _mm_shuffle_ps(self.data, v1,  _ico_shuffle(0, 2, 1, 0));
		}	
	}
}