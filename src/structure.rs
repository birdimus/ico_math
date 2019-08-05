// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use core::arch::x86_64::*;

pub trait SIMDVector1 {
    fn data(self) -> __m128;
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector2 {
    fn data(self) -> __m128;
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector3 {
    fn data(self) -> __m128;
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector4 {
    fn data(self) -> __m128;
    fn data_i(self) -> __m128i;
}
