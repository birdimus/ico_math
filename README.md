# ico_math
An opinionated Rust SIMD Math Library for games and graphics.

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

## future work

- Encoded load / store operations (16 uint, etc).
- Additional matrix methods and types.
- More unit test coverage (specifically bool, int, and matrix).
- More documentation.
- Perfomance tests.
- Basic physics.
- Removing code duplication through use of impl trait, subject to analysis of inlining behavior.

## feedback

This is an opinionated library (as close to c as possible), but I'm also new to Rust, and there are probably things I can do better.
I'm open to suggestions, bug reports, or feature requests.

## ico?
Short for Icosahedra, the plural form of icosahedron, a 20 sided polyhedron.  Also the name of my LLC.
