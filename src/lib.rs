/*

Major TODOs
HASH needs a custom implementation for float, vecs - int should probably store & mask.
Quaternion needs the euler conversion methods ported.
SinCos vec implementations should be used.
Load / Store stuff
More tests

*/



#![no_std]
mod sse_extensions;
mod float_vector;
mod int_vector;
mod vector2;
mod vector3;
mod vector4;
mod vector2_int;
mod vector3_int;
mod vector4_int;
mod vector4_bool;
mod quaternion;
mod dual_quaternion;
mod matrix4x4;

use core::arch::x86_64::*;


pub trait SIMDVector1 {
  fn data(self)->__m128;
}
pub trait SIMDVector2 {
  fn data(self)->__m128;
}
pub trait SIMDVector3 {
  fn data(self)->__m128;
}
pub trait SIMDVector4 {
  fn data(self)->__m128;
}


#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct RawFloatVector{
  pub data : [f32;4],
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct RawIntVector{
  data : [i32;4],
}



#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct FloatVector{
  data : __m128,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct IntVector{
  data : __m128i,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2{
	data : __m128,
}

/// A vector 3
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector3{
	data : __m128,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4{
	data : __m128,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2Int{
	data : __m128i,
}
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector3Int{
	data : __m128i,
}
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4Int{
	data : __m128i,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2Bool{
  data : __m128i,
}
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector3Bool{
  data : __m128i,
}
#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector4Bool{
  data : __m128i,
}


#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Quaternion{
	data : __m128,
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct DualQuaternion{
	real : __m128,
	dual : __m128,
}

#[derive(Copy, Clone)]
#[repr(C, align(16))]
pub struct Matrix4x4{
	m : [__m128; 4],
}




