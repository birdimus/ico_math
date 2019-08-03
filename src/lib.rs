#![no_std]
mod sse_extensions;
mod structure;
mod raw;
mod float_vector;
mod int_vector;
mod vector2;
mod vector2_bool;
mod vector2_int;
mod vector3;
mod vector3_bool;
mod vector3_int;
mod vector4;
mod vector4_bool;
mod vector4_int;
mod quaternion;
mod dual_quaternion;
mod matrix4x4;

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






















