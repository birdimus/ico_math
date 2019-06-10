use std::arch::x86_64::*;
use crate::Vector3;
use crate::Quaternion;
use crate::_ico_shuffle;

impl Quaternion{


	#[inline(always)]
	pub fn add2(v1 : Vector3, v2 : Vector3) -> Vector3{	
		unsafe{
			Vector3{data : _mm_add_ps(v1.data, v2.data)}
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
}