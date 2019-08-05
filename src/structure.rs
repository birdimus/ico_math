// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use core::arch::x86_64::*;

pub trait SIMDVector1 {
	#[inline(always)]
    fn data(self) -> __m128;
    #[inline(always)]
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector2 {
	#[inline(always)]
    fn data(self) -> __m128;
    #[inline(always)]
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector3 {
	#[inline(always)]
    fn data(self) -> __m128;
    #[inline(always)]
    fn data_i(self) -> __m128i;
}
pub trait SIMDVector4 {
	#[inline(always)]
    fn data(self) -> __m128;
    #[inline(always)]
    fn data_i(self) -> __m128i;
}
