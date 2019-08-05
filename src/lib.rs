#![no_std]
pub mod dual_quaternion;
pub mod float_vector;
pub mod int_vector;
pub mod matrix4x4;
pub mod quaternion;
pub mod raw;
pub mod sse_extensions;
mod structure;
pub mod vector2;
pub mod vector2_bool;
pub mod vector2_int;
pub mod vector3;
pub mod vector3_bool;
pub mod vector3_int;
pub mod vector4;
pub mod vector4_bool;
pub mod vector4_int;

/*

Major TODOs
HASH needs a custom implementation for float, vecs - int should probably store & mask.
Quaternion needs the euler conversion methods ported.
SinCos vec implementations should be used.
Load / Store stuff
More tests

*/

//from Knuth
//  bool approximatelyEqual(float a, float b, float epsilon)
// {
//     return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// bool essentiallyEqual(float a, float b, float epsilon)
// {
//     return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// bool definitelyGreaterThan(float a, float b, float epsilon)
// {
//     return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// bool definitelyLessThan(float a, float b, float epsilon)
// {
//     return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }
