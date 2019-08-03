use core::arch::x86_64::*;
use crate::Vector4Bool;
use crate::SIMDVector4;

impl SIMDVector4 for Vector4Bool{
#[inline(always)]
  fn data(self)->__m128{
    return unsafe{_mm_castsi128_ps (self.data)};
  }
}

impl Vector4Bool{
#[inline(always)]
  pub fn new(x : bool, y : bool, z : bool, w : bool) -> Vector4Bool{
    let x_val : u32 = if x {0xFFFFFFFF} else{0};
    let y_val : u32  = if y {0xFFFFFFFF} else{0};
    let z_val : u32  = if z {0xFFFFFFFF} else{0};
    let w_val : u32  = if w {0xFFFFFFFF} else{0};
    unsafe{
      return Vector4Bool{data : _mm_set_epi32(
        core::mem::transmute::<u32, i32>(w_val),
        core::mem::transmute::<u32, i32>(z_val),
        core::mem::transmute::<u32, i32>(y_val), 
        core::mem::transmute::<u32, i32>(x_val))};
    }
  }
  #[inline(always)]
  pub fn x(self) -> bool {
    unsafe{//1 2 4 8
      return (_mm_movemask_epi8 (self.data) & 15 ) == 15;
    }
  }
  #[inline(always)]
  pub fn y(self) -> bool {
    unsafe{ //16 32 64 128
      return (_mm_movemask_epi8 (self.data) & 240 ) == 240;
    }
  }
  #[inline(always)]
  pub fn z(self) -> bool {
    unsafe{ //256 512 1024 2048
      return (_mm_movemask_epi8 (self.data) & 3840 ) == 3840;
    }
  }
  #[inline(always)]
  pub fn w(self) -> bool {
    unsafe{ // 4096 8192 16384 32768
      return (_mm_movemask_epi8 (self.data) & 61440 ) == 61440;
    }
  }


  #[inline(always)]
  pub fn set_x(&mut self, value : bool) {
    let val : u32 = if value {0xFFFFFFFF} else{0};
    unsafe{
       self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 0);
    } 
  }
  #[inline(always)]
  pub fn set_y(&mut self, value : bool) {
    let val : u32 = if value {0xFFFFFFFF} else{0};
    unsafe{
       self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 1);
    } 
  }
  #[inline(always)]
  pub fn set_z(&mut self, value : bool) {
    let val : u32 = if value {0xFFFFFFFF} else{0};
    unsafe{
       self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 2);
    } 
  }
  #[inline(always)]
  pub fn set_w(&mut self, value : bool) {
    let val : u32 = if value {0xFFFFFFFF} else{0};
    unsafe{
       self.data = _mm_insert_epi32(self.data, core::mem::transmute::<u32, i32>(val), 3);
    } 
  }

  #[inline(always)]
  pub fn all(self : Vector4Bool) -> bool{  
    unsafe{
      return (_mm_movemask_epi8 (self.data)  ) == 0x0000FFFF;
    }
  }
  #[inline(always)]
  pub fn any(self : Vector4Bool) -> bool{  
    unsafe{
      return _mm_movemask_epi8 (self.data) != 0;
    }
  }

}

