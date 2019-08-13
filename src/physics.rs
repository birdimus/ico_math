// Copyright 2019 Icosahedra, LLC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//

use crate::float_vector::FloatVector;
use crate::quaternion::Quaternion;
use crate::structure::SIMDVector3;
use crate::vector3::Vector3;
use crate::vector4::Vector4;

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct Box {
    pub center: Vector3,
    pub extents: Vector3,
    pub rotation: Quaternion, //Unorthodox
}
impl Box {
    #[inline(always)]
    pub fn closest_point(self, point: Vector3) -> Vector3 {
        let shifted = point - self.center;
        let local = self.rotation.inverse() * shifted;
        let clamped = local.max(-self.extents).min(self.extents);

        return (self.rotation * clamped) + self.center;
    }

    /// Get the closest point on the box to the center of the sphere, and then do a sqr-distance check.
    #[inline(always)]
    pub fn test_sphere(self, sphere: Sphere) -> bool {
        let closest_point = self.closest_point(sphere.center());
        return sphere.test_point(closest_point);
    }
}

#[derive(Copy, Clone, Debug)]
#[repr(C, align(16))]
pub struct AABB {
    pub center: Vector3,
    pub extents: Vector3,
}
impl AABB {
    #[inline(always)]
    pub fn new(center: Vector3, extents: Vector3) -> AABB {
        return AABB {
            center: center,
            extents: extents,
        };
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

    #[inline(always)]
    pub fn closest_point(self, point: Vector3) -> Vector3 {
        let min_clamp = point.max(self.min());
        return min_clamp.min(self.max());
    }

    #[inline(always)]
    pub fn test_point(self, point: Vector3) -> bool {
        let abs_delta = (self.center - point).abs();
        let outside = abs_delta.less_equal(self.extents);
        return outside.all();
    }

    #[inline(always)]
    pub fn test_ray(self, ray: Ray) -> bool {
        let ray_to_bounds = self.center - ray.origin;
        let origin_inside_bounds = ray_to_bounds.abs().less_equal(self.extents);
        if origin_inside_bounds.all() {
            return true;
        }

        let inverse = Vector3::from(1.0).div(ray.normalized_direction);
        let infinity = Vector3::set(core::f32::INFINITY);
        let infinite = inverse.abs().equal(infinity);
        // If an inverse is infinite AND originOutsideBounds, we fail.
        // Otherwise, all infinte vectors should be ignored.
        // This should catch our NaNs and early out.
        if origin_inside_bounds.andnot(infinite).any() {
            return false;
        }

        let t1 = (ray_to_bounds - self.extents).mul(inverse);
        let t2 = (ray_to_bounds + self.extents).mul(inverse);

        let mut t_min = t1.min(t2);
        let mut t_max = t1.max(t2);

        // Set unused components to 0 on the min vector.
        t_min = Vector3 {
            data: infinite.data(),
        }
        .andnot(t_min);

        t_max = Vector3::select(t_max, infinity, infinite);

        let t_near = t_min.horizontal_max();
        let t_far = t_max.horizontal_min();

        // Point = ray.origin + ray.direction * tnear
        return t_far.greater_equal(t_near);
    }

    #[inline(always)]
    pub fn test_aabb(self, other: AABB) -> bool {
        let abs_delta = (self.center - other.center).abs();
        let extents_sum = self.extents + other.extents;
        let overlap = abs_delta.less_equal(extents_sum);
        return overlap.all();
    }

    #[inline(always)]
    pub fn contains_segment(self, segment: Segment) -> bool {
        return self.test_point(segment.start) && self.test_point(segment.end);
    }

    #[inline(always)]
    pub fn contains_bounds(self, other: AABB) -> bool {
        let min = self.min().less_equal(other.min());
        let max = self.max().greater_equal(other.max());
        return min.and(max).all();
    }

    #[inline(always)]
    pub fn contains_sphere(self, sphere: Sphere) -> bool {
        // Absolute distance vector between sphere and bounds
        let abs_dir = (sphere.center() - self.center).abs() + Vector3::from(sphere.radius());
        // Check that each component is less than the bounds extents
        return abs_dir.less_equal(self.extents).all();
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
    pub fn new(point: Vector3, normal: Vector3) -> Plane {
        let dist = point.dot(normal);
        let mut plane = Plane {
            data: Vector4::from(normal),
        };
        plane.data.set_w(dist);
        return plane;
    }

    #[inline(always)]
    pub fn compute_ccw(point_a: Vector3, point_b: Vector3, point_c: Vector3) -> Plane {
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
    pub fn distance_point(self, point: Vector3) -> FloatVector {
        return (point.dot(self.normal()) - self.distance()).abs();
    }

    #[inline(always)]
    pub fn closest_point(self, point: Vector3) -> Vector3 {
        let distance = self.distance_point(point);
        // point - normal*distance
        return Vector3::neg_mul_add(self.normal(), Vector3::from(distance), point);
    }
    #[inline(always)]
    pub fn closest_point_distance(self, point: Vector3) -> (Vector3, FloatVector) {
        let distance = self.distance_point(point);
        // point - normal*distance
        return (
            Vector3::neg_mul_add(self.normal(), Vector3::from(distance), point),
            distance,
        );
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
    pub fn closest_point(self, point: Vector3) -> Vector3 {
        let ab = self.end - self.start;
        // Clamp 0-1
        let t = (point - self.start)
            .dot(ab)
            .min(FloatVector::new(1.0))
            .max(FloatVector::zero());

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
    pub fn closest_point(self, point: Vector3) -> Vector3 {
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
    pub fn new(point_a: Vector3, point_b: Vector3) -> Line {
        return Line {
            origin: point_a,
            normalized_direction: (point_b - point_a).normalize(),
        };
    }

    #[inline(always)]
    pub fn closest_point(self, point: Vector3) -> Vector3 {
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
    pub fn new<T: Into<FloatVector>>(center: Vector3, radius: T) -> Sphere {
        let mut internal = Vector4 { data: center.data };
        internal.set_w(radius);
        return Sphere { data: internal };
    }
    #[inline(always)]
    pub fn center(self) -> Vector3 {
        return Vector3::from(self.data);
    }

    #[inline(always)]
    pub fn radius(self) -> FloatVector {
        return self.data.w();
    }

    /// Is the point inside or on the sphere.
    #[inline(always)]
    pub fn test_point(self, point: Vector3) -> bool {
        let sqr_dist = Vector4::from((point - self.center()).sqr_magnitude());
        //we only care about w, but we can skip the shuffle.
        let radius_sqr = self.data.mul(self.data);
        let inside = sqr_dist.less_equal(radius_sqr);
        return inside.w();
    }

    /// Does the line intersect or touch the sphere.
    #[inline(always)]
    pub fn test_line(self, line: Line) -> bool {
        let p = line.closest_point(self.center());
        return self.test_point(p);
    }

    /// Does the ray intersect or touch the sphere.
    #[inline(always)]
    pub fn test_ray(self, ray: Ray) -> bool {
        let p = ray.closest_point(self.center());
        return self.test_point(p);
    }

    /// Does the segment intersect or touch the sphere.
    #[inline(always)]
    pub fn test_segment(self, segment: Segment) -> bool {
        let p = segment.closest_point(self.center());
        return self.test_point(p);
    }

    /// Is the sphere intersecting or touching the plane.
    #[inline(always)]
    pub fn test_plane(self, plane: Plane) -> bool {
        let dist = Vector4::from(plane.distance_point(self.center()));
        let inside = dist.less_equal(self.data);
        return inside.w();
    }

    /// Is the sphere fully behind the plane.
    #[inline(always)]
    pub fn test_inside_plane(self, plane: Plane) -> bool {
        let distance = Vector4::from(self.center().dot(plane.normal()) - plane.distance());
        let behind = distance.less(-self.data);
        return behind.w();
    }
    /// Is the sphere intersecting or behind the plane.
    #[inline(always)]
    pub fn test_halfspace_plane(self, plane: Plane) -> bool {
        let distance = Vector4::from(self.center().dot(plane.normal()) - plane.distance());
        let behind = distance.less_equal(self.data);
        return behind.w();
    }

    /// Are the two spheres touching or intersecting.
    #[inline(always)]
    pub fn test_sphere(self, other: Sphere) -> bool {
        let vec = self.center() - other.center();
        let sqr_dist = Vector4::from(vec.sqr_magnitude());
        let total_radius = self.data + other.data;
        let total_radius_sqr = total_radius.mul(total_radius);
        let intersecting = sqr_dist.less_equal(total_radius_sqr);
        return intersecting.w();
    }

    #[inline(always)]
    pub fn contains_segment(self, segment: Segment) -> bool {
        let sqr_dist_start = Vector4::from((segment.start - self.center()).sqr_magnitude());
        let sqr_dist_end = Vector4::from((segment.end - self.center()).sqr_magnitude());
        let radius_sqr = self.data.mul(self.data);
        let contained = sqr_dist_start
            .less_equal(radius_sqr)
            .and(sqr_dist_end.less_equal(radius_sqr));
        return contained.w();
    }

    #[inline(always)]
    pub fn contains_bounds(self, bounds: AABB) -> bool {
        let abs_dir =
            Vector4::from(((self.center() - bounds.center).abs() + bounds.extents).sqr_magnitude());
        let radius_sqr = self.data.mul(self.data);
        let contained = abs_dir.less_equal(radius_sqr);
        return contained.w();
    }

    #[inline(always)]
    pub fn contains_sphere(self, inner: Sphere) -> bool {
        let sqr_dist = Vector4::from((self.center() - inner.center()).sqr_magnitude());
        let radius = self.data - inner.data;
        let r_sqr = radius.mul(radius);
        let contained = sqr_dist.less_equal(r_sqr);
        return contained.w();
    }
}
// #[cfg(test)]
// mod test;
