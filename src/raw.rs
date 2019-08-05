// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

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
