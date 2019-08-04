#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct RawVector_f32 {
    pub data: [f32; 4],
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct RawVector_i32 {
    pub data: [i32; 4],
}
