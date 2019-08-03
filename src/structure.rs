use core::arch::x86_64::*;



pub trait SIMDVector1 {
  fn data(self)->__m128;
  fn data_i(self)->__m128i;
}
pub trait SIMDVector2 {
  fn data(self)->__m128;
  fn data_i(self)->__m128i;
}
pub trait SIMDVector3 {
  fn data(self)->__m128;
  fn data_i(self)->__m128i;
}
pub trait SIMDVector4 {
  fn data(self)->__m128;
  fn data_i(self)->__m128i;
}




