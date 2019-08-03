use core::arch::x86_64::*;
use crate::structure::SIMDVector2;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Vector2Bool{
  pub data : __m128i,
}

impl SIMDVector2 for Vector2Bool{
#[inline(always)]
  fn data(self)->__m128{
    return unsafe{_mm_castsi128_ps (self.data)};
  }
  fn data_i(self)->__m128i{
    return unsafe{ self.data};
  }
}

impl Vector2Bool{
#[inline(always)]
  pub fn new(x : bool, y : bool) -> Vector2Bool{
    let x_val : u32 = if x {0xFFFFFFFF} else{0};
    let y_val : u32  = if y {0xFFFFFFFF} else{0};

    unsafe{
      return Vector2Bool{data : _mm_set_epi32(
        0,
        0,
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
  pub fn all(self : Vector2Bool) -> bool{  
    unsafe{
      return (_mm_movemask_epi8 (self.data) & 255 ) == 255;
    }
  }
  #[inline(always)]
  pub fn any(self : Vector2Bool) -> bool{  
    unsafe{
      return (_mm_movemask_epi8 (self.data) & 255 ) != 0;
    }
  }

}

