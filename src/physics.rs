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
pub struct AABB {
    pub center: Vector3, 
    pub extents: Vector3, 
}
impl AABB {
    #[inline(always)]
    pub fn new(center : Vector3, extents : Vector3) -> AABB {
        return AABB{center:center, extents:extents};
    }
    #[inline(always)]
    pub fn min(self) -> Vector3 {
        return self.center - self.extents;
    }
    #[inline(always)]
    pub fn max(self) -> Vector3 {
        return self.center + self.extents;
    }
    #[inline(always)]
    pub fn center(self) -> Vector3 {
        return self.center;
    }
    #[inline(always)]
    pub fn extents(self) -> Vector3 {
        return self.extents;
    }
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Plane {
    //Points (P) on plane satisfy: distance = dot(normal, P)
    pub data: Vector4, //normal = xyz, w=distance
}

impl Plane {

    #[inline(always)]
    pub fn new(point : Vector3, normal : Vector3) -> Plane {
        let dist = point.dot(normal);
        let mut plane = Plane{data:Vector4::from(normal)};
        plane.data.set_w(dist);
        return plane;
    }

    #[inline(always)]
    pub fn compute_ccw(point_a : Vector3, point_b : Vector3, point_c : Vector3) -> Plane {
        let first = point_b - point_a;
        let second = point_c - point_a;
        let norm = first.cross(second).normalize();
        return Plane::new(point_a, norm);
    }
    

    #[inline(always)]
    pub fn normal(self) -> Vector3 {
        return Vector3::from(self.data);
    }

    #[inline(always)]
    pub fn distance(self) -> FloatVector {
        return self.data.w();
    }

    #[inline(always)]
    pub fn distance_point(self, point : Vector3) -> FloatVector {
        return (point.dot(self.normal()) - self.distance()).abs();
    }
}

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
    pub data: Vector4, //center = xyz, radius = w
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
        let radius_sqr = self.data.mul(self.data);  
        let inside = sqr_dist.less_equal(radius_sqr);
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

    /// Is the sphere intersecting the plane.
    #[inline(always)]
    pub fn test_plane(self, plane : Plane) -> bool {
        let dist = Vector4::from(plane.distance_point(self.center()));
        let inside = dist.less_equal(self.data);
        return inside.w();
    }

    /// Is the sphere fully behind the plane.
    #[inline(always)]
    pub fn test_inside_plane(self, plane : Plane) -> bool {
       let distance = Vector4::from(self.center().dot(plane.normal()) - plane.distance());
       let behind = distance.less(-self.data);
       return behind.w();
       
    }
    /// Is the sphere intersecting or behind the plane.
    #[inline(always)]
    pub fn test_halfspace_plane(self, plane : Plane) -> bool {
       let distance = Vector4::from(self.center().dot(plane.normal()) - plane.distance());
       let behind = distance.less_equal(self.data);
       return behind.w();
       
    }

    #[inline(always)]
    pub fn test_sphere(self, other : Sphere) -> bool {
        let vec = self.center() - other.center();
        let sqr_dist = Vector4::from(vec.sqr_magnitude());
        let total_radius = self.data + other.data;
        let total_radius_sqr = total_radius.mul(total_radius);
        let intersecting = sqr_dist.less_equal(total_radius_sqr);
        return intersecting.w();
    }
}
// #[cfg(test)]
// mod test;
