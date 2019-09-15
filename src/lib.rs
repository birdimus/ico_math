// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

#![cfg_attr(not(any(test, feature = "std")), no_std)]

mod dual_quaternion;
mod float_vector;
mod int_vector;
mod matrix4x4;
pub mod physics;
mod quaternion;
pub mod raw;
pub mod sse_extensions;
pub mod structure;
mod vector2;
mod vector2_bool;
mod vector2_int;
mod vector3;
mod vector3_bool;
mod vector3_int;
mod vector4;
mod vector4_bool;
mod vector4_int;

pub use dual_quaternion::DualQuaternion;
pub use float_vector::FloatVector;
pub use int_vector::IntVector;
pub use matrix4x4::Matrix4x4;
pub use quaternion::Quaternion;
pub use quaternion::RotationOrder;
pub use vector2::Vector2;
pub use vector2_bool::Vector2Bool;
pub use vector2_int::Vector2Int;
pub use vector3::Vector3;
pub use vector3_bool::Vector3Bool;
pub use vector3_int::Vector3Int;
pub use vector4::Vector4;
pub use vector4_bool::Vector4Bool;
pub use vector4_int::Vector4Int;

#[cfg(not(any(test, feature = "std")))]
#[panic_handler]
fn my_panic(_info: &core::panic::PanicInfo) -> ! {
    loop {}
}

//Is this required?
//RUSTFLAGS='-C target-features=+avx2,+fma'
