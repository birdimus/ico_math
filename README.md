# ico_math
An opinionated SIMD Math Library for games and graphics in Rust.

The goal of this library is to facilitate the development of low-level game engine functionality, as well as provide a user-friendly interface to common math operations without sacrificing performance.

This library is still under development and may change.  Unit test coverage and documentation is not complete.  Where possible, I've tried to document error bounds (such as on SSE sin, cos, acos approximations).

## features
Ico Math supports the standard primitive types (Vector2, Vector3, Vector4, Quaternion, Matrix4x4, DualQuaternion), as well as native SIMD integer types (Vector2Int, Vector3Int, Vector4Int), boolean masks (Vector2Bool, Vector3Bool, Vector4Bool), and broadcast types (FloatVector, IntVector) where all elements are the same.

The library also provides aligned raw vector types for load / store operations.

Additionally, the library provides a complete set of CG style swizzle operators, such as vec.zwxy(), or vec.xxyy().

## requirements and dependencies

- Only uses the core::arch::x86_64 module.
- No external dependencies.
- #![no_std]
- Requires FMA and SSE 4.1

## conventions

- Right Handed Coordinate System
- Column Major
- Radians unless specified
- Equality is bitwise unless specified (for instance. Vector3.approx_equals() is approximate).

This library takes a middle-ground approach to type safety.  Vector2, Vector3, Vector4, and Quaternion are concrete types.  For perfomance reasons there are separate integer based Vector2Int, Vector3Int and Vector4Int types.  However, there are no 'point' types, nor are there 'radian' or 'degree' types.

Most methods require explicit casts to the correct type, from conversions between most types are provided.  The only exception are bitwise operators, which allow operations between any same sized Vector type (int, bool, float).  Vector3Int XOR Vector3 is allowed, but Vector2 XOR Vector3 is not allowed.

Casts upward (from Vector2 to Vector4) always zero-initialize the additional channels (zw, in this case). Casts downward (Vector4 to Vector2) do not zero channels. 

## future work

- Encoded load / store operations (16 uint, etc).
- Additional matrix methods and types.
- Methods for shifts which take const arguments (without const generics).
- More unit test coverage (specifically bool, and matrix).
- More documentation.
- Perfomance tests.
- Basic physics.
- Removing code duplication through use of impl trait, subject to analysis of inlining behavior.

## feedback

I'm new to Rust, and there are certainly things I can do better.
I'm open to feedback, and love bug reports.  

There is presently a lot of code duplication I intend to remove once I can validate inlining perfomance.  I've resisted macros in favor of readability.

## ico?
Short for Icosahedra, the plural form of icosahedron, a 20 sided polyhedron.  Also the name of my LLC.

Copyright 2019 Icosahedra, LLC

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
