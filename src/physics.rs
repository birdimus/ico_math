// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::float_vector::FloatVector;
use crate::vector3::Vector3;
use crate::vector4::Vector4;


#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Segment {
    pub start: Vector3,
    pub end: Vector3,
}

impl Segment {
    #[inline(always)]
    pub fn closest_point(self, point : Vector3) -> Vector3 {

        let ab = self.end - self.start;
        // Clamp 0-1
        let t =  (point - self.start).dot(ab).min(FloatVector::new(1.0)).max(FloatVector::zero());
    
        let t_scaled = Vector3::from(t / ab.sqr_magnitude());
        return Vector3::mul_add(ab, t_scaled, self.start);
    }
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Ray {
    pub origin: Vector3,
    pub normalized_direction: Vector3,
}
impl Ray {


    #[inline(always)]
    pub fn closest_point(self, point : Vector3) -> Vector3 {
        let mut ab_cos = Vector3::from((point - self.origin).dot(self.normalized_direction));
        ab_cos = ab_cos.max(Vector3::zero());
        return Vector3::mul_add(ab_cos, self.normalized_direction, self.origin);
    }

}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Line {
    pub origin: Vector3,
    pub normalized_direction: Vector3,
}
impl Line {

    #[inline(always)]
    pub fn new(point_a : Vector3, point_b : Vector3) -> Line {
        return Line{origin:point_a, normalized_direction:(point_b - point_a ).normalize()};
    }

    #[inline(always)]
    pub fn closest_point(self, point : Vector3) -> Vector3 {
        let ab_cos = Vector3::from((point - self.origin).dot(self.normalized_direction));

        return Vector3::mul_add(ab_cos, self.normalized_direction, self.origin);
    }

}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Sphere {
    pub data: Vector4,
}

impl Sphere {

    
    #[inline(always)]
    pub fn new<T: Into<FloatVector>>(center : Vector3, radius : T) -> Sphere {
        let mut internal = Vector4{data:center.data};
        internal.set_w(radius);
        return Sphere{data:internal};
    }
    #[inline(always)]
    pub fn center(self) -> Vector3 {
        return Vector3::from(self.data);
    }

    #[inline(always)]
    pub fn radius(self) -> FloatVector {
        return self.data.w();
    }

    #[inline(always)]
    pub fn test_point(self, point : Vector3) -> bool {
        let sqr_dist = Vector4::from((point - self.center()).sqr_magnitude());
        //we only care about w, but we can skip the shuffle.
        let self_sqr = self.data.mul(self.data);  
        let inside = sqr_dist.less_equal(self_sqr);
        return inside.w();
    }

    #[inline(always)]
    pub fn test_line(self, line : Line) -> bool {
        let p = line.closest_point(self.center());
        return self.test_point(p);
    }

    #[inline(always)]
    pub fn test_ray(self, ray : Ray) -> bool {
        let p = ray.closest_point(self.center());
        return self.test_point(p);
    }

    #[inline(always)]
    pub fn test_segment(self, segment : Segment) -> bool {
        let p = segment.closest_point(self.center());
        return self.test_point(p);
    }
}
// #[cfg(test)]
// mod test;
