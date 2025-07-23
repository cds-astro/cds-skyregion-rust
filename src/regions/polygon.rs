//! Polygon on the unit sphere.  
//! This code contains a re-write (sometime very close to the copy/paste)
//! of codes I originally write for [cdshealpix](https://github.com/cds-astro/cds-healpix-rust).
//!
use std::{cmp::Ordering, collections::BTreeSet, ops::Range};

use cdshealpix::{
  nested::{bmoc::BMOC, custom_polygon_coverage},
  sph_geom::ContainsSouthPoleMethod,
};

use crate::{
  common::{
    error::SkyRegionError,
    lat_deg2rad, lon_deg2rad,
    math::{
      cross_product, lonlat_to_xyz, minus, norm, normalized, scalar_product, time, Coo,
      RefToLocalRotMatrix, UnitVec3, Vec3, HALF_PI, PI, TWICE_PI,
    },
  },
  SkyRegion,
};

#[derive(Debug)]
enum PolygonType {
  /// Convex polygon: https://en.wikipedia.org/wiki/Convex_polygon
  Convex,
  /// Concave polygon: https://en.wikipedia.org/wiki/Simple_polygon
  Concave,
  // see polygon tesselation (or split into simple polygons)
  // see: https://en.wikipedia.org/wiki/Boolean_operations_on_polygons
  // see: https://en.wikipedia.org/wiki/Greiner%E2%80%93Hormann_clipping_algorithm
  // see: https://stackoverflow.com/questions/4876065/is-there-an-easy-and-fast-way-of-checking-if-a-polygon-is-self-intersecting
  /// Self intersecting polygon: https://en.wikipedia.org/wiki/List_of_self-intersecting_polygons
  SelfIntersecting,
}

impl PolygonType {
  fn area(&self, poly: &Polygon) -> f64 {
    match self {
      PolygonType::Convex => {
        // https://en.wikipedia.org/wiki/Spherical_trigonometry#Area_and_spherical_excess
        // somme of the N internal angles - (N-2) * pi
        // N = number of side
        // We can get the angle from the segment plane normals: |sin(angle)| = ||n1 x n2||
        // the problem is that the angle returned by an arcsin is in [0, pi/2] instead of [0, pi].
        // To compute the cosine, we use vector differences:
        //   d1 = v1 - cos(dist(v1, v2))v2 and d2 = v3 - cos(dist(v2, v3))v2
        //   we need the cos(angular distance) so that the diff vector is perpendicular to the common
        //   line defined by the common vertex v2
        // then
        //   cos(angle) = d1.d2 / ||d1||*||d2||
        // If v1 (or v3) and v2 are very close, the norm of d1 may be very close to 0
        // (but it should not be problematic).
        // And we use:
        //   angle = arctan2(|sin(angle)|, cos(angle))
        let sum = poly
          .segments
          .iter()
          .zip(poly.segments.iter().cycle().skip(1))
          .map(|(lhseg, rhseg)| {
            assert_eq!(lhseg.v2_lon, rhseg.v1_lon); // l.v2 = r.v1 = common vertex
                                                    // We normalized normals to avoid dividing by the product of two possibly small norms
            let sin_ang = norm(&cross_product(&normalized(&lhseg.n), &normalized(&rhseg.n)));
            // We get already normalized normals to avoid dividing by the product of two possibly small norms
            let cos_ang = scalar_product(
              &lhseg.normal_to_v2_in_v1v2_plane(),
              &rhseg.normal_to_v1_in_v1v2_plane(),
            );
            sin_ang.atan2(cos_ang)
          })
          .sum::<f64>();
        sum - (poly.segments.len() - 2) as f64 * PI
      }
      PolygonType::Concave => {
        // We use the formula provided in https://en.wikipedia.org/wiki/Spherical_trigonometry
        // * s = 1/2 (a + b + c)
        // *  tan S/2 = [tan(a/2) tan(b/2) sin(C)] / [1 + tan(a/2) tan(b/2) cos(C)]
        // Wikpedia: "Because some triangles are badly characterized by their edges
        // (e.g., if a = b ≈ 1/2 c {\displaystyle a=b\approx {\frac {1}{2}}c} a = b \approx \frac12c),
        // it is often better to use the formula for the excess in terms of two edges and their included angle"
        // Here we use
        // * C = Delta(lon)
        // * a = pi/2 - lat1
        // * b = pi/2 - lat2
        let area = poly
          .segments
          .iter()
          .map(|seg| {
            let dlon = (seg.v2_lon - seg.v1_lon).abs();
            let (dlon, cross_primary_meridian) = if dlon > PI {
              (dlon - PI, true)
            } else {
              (dlon, false)
            };
            let (sin_dlon, cos_dlon) = dlon.sin_cos();
            let tan_ha_tanhb =
              (0.5 * (HALF_PI - seg.v1_lat)).tan() * (0.5 * (HALF_PI - seg.v2_lat)).tan();
            let triangle_area =
              2.0 * ((tan_ha_tanhb * sin_dlon) / (1.0 + tan_ha_tanhb * cos_dlon)).atan();
            let v1_is_left = (seg.v1_lon < seg.v2_lon) != cross_primary_meridian;
            if v1_is_left {
              triangle_area
            } else {
              -triangle_area
            }
          })
          .sum::<f64>()
          .abs();
        if poly.contains_south_pole {
          4.0 * PI - area
        } else {
          area
        }
      }
      PolygonType::SelfIntersecting => {
        eprintln!("WARNING: the area returned for self-intersecting polygons is wrong!");
        eprintln!("         The polygon MUST BE split into simple polygons first!");
        PolygonType::Concave.area(poly)
      }
    }
  }
}

/// A polygon divides the unit sphere into two complementary areas.
/// * In case one of the area contains both poles, we choose by default the
///   area that do not contains any pole.
/// * In case both area contains one pole, we chose the area containing the pole
///   in the hemisphere containing the polygon gravity center.
///
/// You can revert this behaviour by setting the `complement` flag.  
///
/// Another approach is to provide a control point which MUST be inside the
/// polygon (the control point MUST NOT be a vertex or along an edge).
///
/// Another possible approach is to consider that the opposite of the gravity
/// center of the polygon is out of the polygon (we consider the opposite because
/// for self intersecting polygons the gravity center may also easily be out
/// of the polygon).
///
/// Another approach may be to compute the area of the polygon and choose by
/// default the smallest area. Unfortunately, computing the area is far from
/// trivial for self-intersecting polygons.
/// * convex polygon: can't be larger than an hemisphere, so cannot contains
///   both poles (except if poles are vertices or on edges).
///   So if the polygon (with complement = false) is convex and its surface area
///   is larger than 2PI, we must take the complement by default.
/// * concave polygons: the computed area is always the one that do not contains
///   both poles. Again, if the computed area is larger than 2PI,
///   we must take the complement by default.
/// * self-intersecting polygons: to compute the area, the polygon must be split
///   in simple polygons, but the algorithm (on the sphere) if not trivial and
///   requires more time to be investigated (to be adapted from Bentley-Ottmann).
///   In this particular case, the default may not be necessarily the polygon
///   having the smallest surface area.
///
// https://www.mathopenref.com/coordpolygonarea2.html
// https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm
//
// May be helpful to compute area of a zone
//  * Area of a lune = 2*lune_angle: https://www.sjsu.edu/faculty/watkins/sphere.htm
#[derive(Debug)]
pub struct Polygon {
  /// Vertices list. For each vertices, coordinates:
  /// * `lon` must be in `[0, 2\pi[`
  /// * `lat` must be in `[-\pi/2, \pi/2]`
  vertices: Vec<Coo>,
  /// One entry per segment. A segment here is made of:
  /// * the LHS vertex longitude
  /// * the RHS vertex longitude
  /// * both vertices (not normalized) vectorial product pointing to in the north hemisphere.
  ///
  /// It means that we make a copy of each vertex longitude,
  /// but all necessary data is aligned in memory.
  segments: Vec<Segment>,
  /// Polygon type
  ptype: PolygonType,
  /// `true` if the south pole is in the polygon
  contains_south_pole: bool,
}

impl Polygon {
  pub fn from_deg(vertices_deg: Vec<(f64, f64)>, complement: bool) -> Result<Self, SkyRegionError> {
    let mut vertices = vertices_deg;
    for (l, b) in &mut vertices {
      *l = lon_deg2rad(*l)?;
      *b = lat_deg2rad(*b)?;
    }
    Ok(Self::new(vertices, complement))
  }

  /// For all (both simple and self-intersecting) polygons:
  /// * if at least one meridian do not intersects the polygon, both poles are
  ///   considered outside the polygon.
  /// * else (i.e. all meridian intersect the polygon), the pole in the hemisphere
  ///   containing the polygon gravity center is considered as being inside the polygon.
  ///
  /// For a different approach, provided a control point that must be inside
  /// the polygon, see [new_with_control_point](#new_with_control_point).
  /// # Args:
  /// * `complement`: by default
  ///   If `complement` is set to `true`, then the inside of the polygon is the
  ///   complement surface of the default one.
  /// # Panics
  /// * if `vertices.len() < 3`
  /// * if, for each vertex:
  ///     + if `lon` not in `[0, 2\pi[`
  ///     + if `lat` not in `[-\pi/2, \pi/2]`
  pub fn new(vertices: Vec<Coo>, complement: bool) -> Self {
    assert!(vertices.len() > 2);
    let vertices_xyz = vertices
      .iter()
      .map(|(lon, lat)| {
        assert!(0.0 <= *lon && *lon < TWICE_PI);
        assert!(-HALF_PI <= *lat && *lat <= HALF_PI);
        lonlat_to_xyz(*lon, *lat)
      })
      .collect::<Vec<UnitVec3>>();
    let vertices_lon = vertices.iter().map(|&(lon, _)| lon).collect::<Vec<f64>>();
    let vertices_lonlat_xyz = vertices
      .iter()
      .cloned()
      .zip(vertices_xyz.iter().cloned())
      .collect::<Vec<((f64, f64), UnitVec3)>>();
    // Compute segments
    let mut normals: Vec<Vec3> = Vec::with_capacity(vertices.len() + 1);
    let segments = vertices_lonlat_xyz
      .iter()
      .zip(vertices_lonlat_xyz.iter().cycle().skip(1))
      .enumerate()
      .map(
        |(i, (((v1_lon, v1_lat), v1_xyz), ((v2_lon, v2_lat), v2_xyz)))| {
          let n = cross_product(v1_xyz, v2_xyz);
          normals.push(n);
          let (x, y, z) = n;
          if z < 0.0 {
            Segment::new(
              i as u32,
              *v1_lon,
              *v1_lat,
              *v2_lon,
              *v2_lat,
              *v1_xyz,
              *v2_xyz,
              (-x, -y, -z),
            )
          } else {
            Segment::new(
              i as u32,
              *v1_lon,
              *v1_lat,
              *v2_lon,
              *v2_lat,
              *v1_xyz,
              *v2_xyz,
              (x, y, z),
            )
          }
        },
      )
      .collect::<Vec<Segment>>();
    // Compute the direction of the gravity center
    let g: UnitVec3 = gravity_center_direction(&vertices_xyz[..]);
    let is_g_dot_n_positive = scalar_product(&g, &normals[0]) >= 0.0;
    // Test if convex or not
    // Remark: we intentionally left aside the last-to-1st segment comparison
    //         because the number of negative scalar product should be equals to
    //         or greater than 2 if the polygon is not convex (if not, then
    //         a '.cylce()' must be added.
    let is_simple = is_simple(&segments[..]);
    let is_convex = is_simple
      && normals
        .iter()
        .any(|n| is_g_dot_n_positive != (scalar_product(&g, n) >= 0.0));
    let ptype = if is_convex {
      PolygonType::Convex
    } else if is_simple {
      PolygonType::Concave
    } else {
      PolygonType::SelfIntersecting
    };
    // Computes mean latitude
    /*let lat_mean = vertices.iter()
    .fold(0.0_f64, |lat_sum, (_, lat)| lat_sum + lat)
    / vertices.len() as f64;*/
    // Remark: at the end, sum_dlon should be = 0 or -2PI or 2PI
    // (we simply test if its absolute value is smaller than PI
    let contains_south_pole = vertices_lon
      .iter()
      .zip(vertices_lon.iter().cycle().skip(1))
      .fold(0.0, |sum_dlon, (lon1, lon2)| {
        let dlon = lon2 - lon1;
        let abs_dlon = dlon.abs();
        if abs_dlon <= PI {
          sum_dlon + dlon
        } else if dlon > 0.0 {
          sum_dlon - (TWICE_PI - abs_dlon)
        } else {
          sum_dlon + (TWICE_PI - abs_dlon)
        }
      })
      .abs()
      > PI
      && g.2 < 0.0;
    // compute surface area?
    Self {
      vertices,
      segments,
      ptype,
      contains_south_pole: if complement {
        !contains_south_pole
      } else {
        contains_south_pole
      },
    }
  }

  /*pub fn new_with_control_point(vertices: Vec<Coo>, pos_in_polygone: Coo) -> Polygon {
    let mut polygon = Polygon::new(vertices, false);
    if !polygon.contains(pos_in_polygone.0, pos_in_polygone.1) {
      polygon.contains_south_pole = !polygon.contains_south_pole;
    }
    polygon
  }*/

  pub fn from_box_deg(
    lon_deg: f64,
    lat_deg: f64,
    a_deg: f64,
    b_deg: f64,
    pa_deg: f64,
  ) -> Result<Self, SkyRegionError> {
    let lon = lon_deg2rad(lon_deg)?;
    let lat = lat_deg2rad(lat_deg)?;
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let pa = pa_deg.to_radians();
    if a <= 0.0 || HALF_PI <= a {
      Err(SkyRegionError::UnexpectedBoxA { a_deg })
    } else if b <= 0.0 || a <= b {
      Err(SkyRegionError::UnexpectedBoxB { b_deg })
    } else if pa <= 0.0 || HALF_PI <= pa {
      Err(SkyRegionError::UnexpectedBoxPA { pa_deg })
    } else {
      Ok(Self::from_box(lon, lat, a, b, pa))
    }
  }

  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2]`
  /// * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
  /// * if `b` not in `]0, a[`
  /// * if `pa` not in `[0, \pi[`
  pub fn from_box(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> Self {
    assert!((0.0..TWICE_PI).contains(&lon));
    assert!((-HALF_PI..=HALF_PI).contains(&lat));
    assert!(0.0 < a && a <= HALF_PI);
    assert!(0.0 < b && b < a);
    assert!((0.0..PI).contains(&pa));
    // Compute spherical coordinates
    let frame_rotation = RefToLocalRotMatrix::from_center(lon, lat);
    // By application of the Thales theorem, the new point has the property:
    //   sin(new_lat) / sin(dlat) = (cos(dlon) * cos(new_lat)) / cos(dlat)
    // With (imagine looking a the triangle in the (x, z) plane:
    // * old_x = cos(lat)
    // * new_x = cos(dlon) * cos(new_lat)
    // * old_z = sin(dlat)
    // * new_z = sin(new_lat)
    // Leading to:
    //   tan(new_lat) = cos(dlon) * tan(dlat)
    let lon = a;
    let (sin_lon, cos_lon) = lon.sin_cos();
    let lat = (cos_lon * b.tan()).atan();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_pa, cos_pa) = pa.sin_cos();
    // Rotation by the position angle
    // - upper right (before rotation by PA)
    let (x1, y1, z1) = (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat);
    // - apply rotation (sin and cos are revere since theta = pi/2 - pa)
    let (y2, z2) = (y1 * sin_pa - z1 * cos_pa, y1 * cos_pa + z1 * sin_pa);
    let mut vertices = Vec::with_capacity(4);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower right (before rotation by PA) = (y1, -z1)
    let (y2, z2) = (y1 * sin_pa + z1 * cos_pa, y1 * cos_pa - z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower left (before rotation by PA) = (-y1, -z1)
    let (y2, z2) = (-y1 * sin_pa + z1 * cos_pa, -y1 * cos_pa - z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - upper left (before rotation by PA) = (-y1, z1)
    let (y2, z2) = (-y1 * sin_pa - z1 * cos_pa, -y1 * cos_pa + z1 * sin_pa);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    Self::new(vertices, false)
  }

  /// # Input:
  /// * `lon` longitude of the position to be tested.
  /// * `x, y, z`: Euclidean coordinates of the position to be tested
  ///
  /// # Remark:
  /// * `lon` is assumed to fit with `x, y, z`, i.e. `lon = atan2(y, x)`
  ///   and `lon` in `[0, 2\pi[` radians (this is **not** checked!)  
  fn odd_num_intersect_going_south(&self, lon: f64, x: f64, y: f64, z: f64) -> bool {
    self.segments.iter().fold(
      false,
      |c,
       &Segment {
         i: _,
         v1_lon,
         v1_lat: _,
         v2_lon,
         v2_lat: _,
         v1: _,
         v2: _,
         n,
       }| {
        if is_in_lon_range(lon, v1_lon, v2_lon) && cross_plane_going_south(x, y, z, n.0, n.1, n.2) {
          !c
        } else {
          c
        }
      },
    )
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    custom_polygon_coverage(
      depth,
      &self.vertices,
      if self.contains_south_pole {
        &ContainsSouthPoleMethod::ContainsSouthPole
      } else {
        &ContainsSouthPoleMethod::DoNotContainsSouthPole
      },
      true,
    )
  }
}

impl SkyRegion for Polygon {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    let (x, y, z) = lonlat_to_xyz(lon, lat);
    self.contains_south_pole != self.odd_num_intersect_going_south(lon, x, y, z)
  }

  fn area(&self) -> f64 {
    self.ptype.area(self)
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    self.to_bmoc(depth).to_flagged_ranges()
  }
}

/// Returns `true` if the given point `p` longitude is between the given vertices
/// `v1` and `v2` longitude range.
fn is_in_lon_range(lon: f64, v1_lon: f64, v2_lon: f64) -> bool {
  // Lets note
  //   - lonA = v1_lon
  //   - lonB = v2_lon
  //   - lon0 = lon
  // When (lonB - lonA).abs() <= PI
  //   => lonB > lon0 != lonA > lon0  like in PNPOLY
  //   A    B    lonA <= lon0 && lon0 < lonB
  // --[++++[--
  //   B    A    lonB <= lon0 && lon0 < lonA
  //
  // But when (lonB - lonA).abs() > PI, then the test must be
  //  =>   lonA >= lon0 == lonB >= lon0
  // <=> !(lonA >= lon0 != lonB >= lon0)
  //    A  |  B    (lon0 < lonB) || (lonA <= lon0)
  //  --[++|++[--
  //    B  |  A    (lon0 < lonA) || (lonB <= lon0)
  let dlon = v2_lon - v1_lon;
  if dlon < 0.0 {
    (dlon >= -PI) == (v2_lon <= lon && lon < v1_lon)
  } else {
    (dlon <= PI) == (v1_lon <= lon && lon < v2_lon)
  }
}

/// Returns `true` if the line at constant `(x, y)` and decreasing `z` going
/// from the given point toward south intersects the plane of given normal vector.
/// The normal vector must have a positive z coordinate (=> is in the north hemisphere)
fn cross_plane_going_south(x: f64, y: f64, z: f64, nx: f64, ny: f64, nz: f64) -> bool {
  // If the scalar product is positive, the point is above the plane, and
  // thus intersect it going south
  scalar_prod(x, y, z, nx, ny, nz) > 0.0
}

/// Reminder: for unit vector, the scalar product is the cosine of the angle
/// between both input vector
fn scalar_prod(x: f64, y: f64, z: f64, nx: f64, ny: f64, nz: f64) -> f64 {
  x * nx + y * ny + z * nz
}

/// Compute the direction of the gravity center of the input list of vertices.
/// If the center of gravity is the center of the sphere, then we remove
/// the last vertex (recursively).
fn gravity_center_direction(vertices: &[UnitVec3]) -> UnitVec3 {
  let (gx, gy, gz) = vertices
    .iter()
    .fold((0.0f64, 0.0f64, 0.0f64), |(gx, gy, gz), &(x, y, z)| {
      (gx + x, gy + y, gz + z)
    });
  let norm = (gx * gx + gy * gy + gz * gz).sqrt();
  if norm > 0.0 {
    (gx / norm, gy / norm, gz / norm)
  } else {
    gravity_center_direction(&vertices[..vertices.len() - 1])
  }
}

#[derive(Debug)]
struct Segment {
  /// index of the segment (for uniq identification)
  i: u32,
  /// Longitude of the first vertex of the segment
  v1_lon: f64,
  /// Latitude of the first vertex of the segment
  v1_lat: f64,
  /// Longitude of the second vertex of the segment
  v2_lon: f64,
  /// Latitude of the second vertex of the segment
  v2_lat: f64,
  /// First vertex Cartesian coordinates
  v1: UnitVec3,
  /// Second vertex Cartesian coordinates
  v2: UnitVec3,
  /// Normal vector to the plane defined by both vertices
  n: Vec3,
}

impl Ord for Segment {
  fn cmp(&self, other: &Self) -> Ordering {
    self.i.cmp(&other.i)
  }
}

impl PartialOrd for Segment {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl Eq for Segment {}

impl PartialEq for Segment {
  fn eq(&self, other: &Self) -> bool {
    self.i == other.i
  }
}

impl Segment {
  #[allow(clippy::too_many_arguments)]
  fn new(
    i: u32,
    v1_lon: f64,
    v1_lat: f64,
    v2_lon: f64,
    v2_lat: f64,
    v1: UnitVec3,
    v2: UnitVec3,
    n: Vec3,
  ) -> Self {
    Self {
      i,
      v1_lon,
      v1_lat,
      v2_lon,
      v2_lat,
      v1,
      v2,
      n,
    }
  }

  /*fn arc_length(&self) -> f64 {
    let cos = scalar_product(&self.v1, &self.v2);
    let sin = norm(&cross_product(&self.v1, &self.v2));
    // Here sin is an absolute value, it's ok since we want a positive angle.
    sin.atan2(cos)
  }*/

  fn normal_to_v1_in_v1v2_plane(&self) -> UnitVec3 {
    // * Project v1 on v2 (warning: test if scalar prod is not already = 0)
    // * n = v1 - projected v1
    let cos_arc_length = scalar_product(&self.v1, &self.v1);
    if cos_arc_length == 0.0 {
      self.v2
    } else {
      let v2_proj_on_v1 = time(cos_arc_length, &self.v1);
      normalized(&minus(&self.v2, &v2_proj_on_v1))
    }
  }

  fn normal_to_v2_in_v1v2_plane(&self) -> UnitVec3 {
    // * Project v1 on v2 (warning: test if scalar prod is not already = 0)
    // * n = v1 - projected v1
    let cos_arc_length = scalar_product(&self.v1, &self.v1);
    if cos_arc_length == 0.0 {
      // Case pi/2
      self.v1
    } else {
      let v1_proj_on_v2 = time(cos_arc_length, &self.v2);
      normalized(&minus(&self.v1, &v1_proj_on_v2))
    }
  }

  fn cross_primary_meridian(&self) -> bool {
    (self.v2_lon - self.v1_lon).abs() > PI
  }

  /// WARNING: also returns true is the segments intersects the opposite segment
  /// (the segment made by the opposite of both vertices)
  fn intersects(&self, rhs: &Segment) -> bool {
    // great-circle arcs intersection = (n1 x n2)
    // ip1 = n1 x (n1 x n2)
    //     = (n1.n2)n1 - (n1.n1)n2
    //     = (n1.n2)n1 - n2
    // ip2 = n2 x (n1 x n2)
    //     = (n2.n2)n1 - (n2.n1)n2
    //     = n1 - (n2.n1)n2
    let n1_dot_n2 = scalar_product(&self.n, &rhs.n);
    let ip = (
      n1_dot_n2 * self.n.0 - rhs.n.0,
      n1_dot_n2 * self.n.1 - rhs.n.1,
      n1_dot_n2 * self.n.2 - rhs.n.2,
    );
    // let i: Vec3 = cross_product(&self.n, &rhs.n);
    // We could have used "is_in_lon_range" to check if i or its opposite is the correct vector
    // let ip = cross_product(&i, &self.n);
    let in_self = (scalar_product(&self.v1, &ip) > 0.0) != (scalar_product(&self.v2, &ip) > 0.0);
    let ip = (
      self.n.0 - n1_dot_n2 * rhs.n.0,
      self.n.1 - n1_dot_n2 * rhs.n.1,
      self.n.2 - n1_dot_n2 * rhs.n.2,
    );
    //let ip = cross_product(&i, &rhs.n);
    let in_rhs = (scalar_product(&rhs.v1, &ip) > 0.0) != (scalar_product(&rhs.v2, &ip) > 0.0);
    // Remark: here we actually test if vec_i or its opposite is inside seg
    //                           and if vec_i or its opposite is inside segment
    // So we can e.g. have vec_i inside seg and its opposite inside segment.
    // but we have checked (externally!!) that longitudes ranges are compatible.
    in_self && in_rhs
  }
}

/// Event triggered when the sweep line moves.  
/// One event per segment vertex.
struct Event<'a> {
  /// Longitude of the vertex triggering the event.
  lon: f64,
  /// Tells if the vertex is the left vertex (i.e. smallest longitude when
  /// `cross_pm` is false, the largest longitude is `cross_pm` is true).
  is_left: bool,
  /// `true` if the segments crosses the Prime Meridian
  cross_pm: bool,
  /// Reference pointing to the full segment
  segment_ref: &'a Segment,
}

/// Contains all segments intersecting the sweep line at the current position.
struct SweepLine<'a> {
  /// Current segments crossing the sweep line
  crossing_segments: BTreeSet<&'a Segment>,
}

impl<'a> SweepLine<'a> {
  /// Init the sweep line with segments crossing the prime meridian.
  fn new(events: &[Event<'a>]) -> SweepLine<'a> {
    SweepLine {
      crossing_segments: events
        .iter()
        .filter(|&event| event.cross_pm && !event.is_left)
        .map(|event| event.segment_ref)
        .collect(),
    }
  }

  /*fn size(&self) -> usize {
    self.crossing_segments.len()
  }*/

  fn insert(&mut self, segment_ref: &'a Segment) {
    assert!(self.crossing_segments.insert(segment_ref));
  }

  fn remove(&mut self, segment_ref: &Segment) {
    //-> Option<&'a Segment> {
    // let pos = self.crossing_segments.iter().position(|&seg_ref| seg_ref == segment_ref)?;
    // Some(self.crossing_segments.remove(pos))
    assert!(self.crossing_segments.remove(segment_ref));
  }

  fn contains_segment_crossing(&self, segment: &Segment) -> bool {
    self
      .crossing_segments
      .iter()
      .any(|&seg| segment.intersects(seg))
  }
}

// Complex calculations (like Bentley-Ottmann Algorithm (sweep line)

/// Implementation of the [Shamos-Hoey](https://www.webcitation.org/6ahkPQIsN)
/// Algorithm adapted for polygons on the unit sphere.
/// * In spherical, we cannot order along `lat` because a great circle arc
///   may contains a point with a latitude larger that both vertices latitude
///   (e.g. considering a vertices having similar latitudes and
///   `Delta Lon \approx \pi`). It is not the case with longitudes.
/// * Thus, I am not sure we can use the above-below relation.
///   E.g. the largest latitude of a great-circle arc with a Delta(Lon) about
///   PI can be well above both vertices latitude. This large great circle arc (L)
///   may intersect a small great circle arc (S) around the largest latitude
///   while two other small great circle arcs below (B) and above (A) S
///   have latitudes larger than the large L latitudes.
///   To cope with this problem, one should probably compute that latitude
///   intersecting the SL at each event (and re-order the SL stack accordingly).
/// * With respect to the Euclidean case, we have to be careful with segments
///   crossing the primary meridian (the sweep line must start containing such
///   segments). In the Bentley-Ottmann Algorithm, the intersection
///   must be added only if its longitude is below PI, if it is above PI, it
///   will be added at the end of the sweep line, when the segment crossing
///   the primary meridian are inserted again in the sweep line.
fn is_simple(segments: &[Segment]) -> bool {
  let mut events = build_events(segments);
  events.sort_by(|lhs, rhs| {
    // Unwrap is safe here since we already checked that 0 <= lon < 2pi in Polygon::new()
    let order = lhs.lon.partial_cmp(&rhs.lon).unwrap();
    match &order {
      // The ordering is chosen such that removal happens before insertion
      Ordering::Equal => match (lhs.is_left, rhs.is_left) {
        (false, true) => Ordering::Less,
        (true, false) => Ordering::Greater,
        _ => Ordering::Equal,
      },
      _ => order,
    }
  });
  // Init the sweep line
  let mut sweep_line = SweepLine::new(&events[..]);
  for event in events {
    if event.is_left {
      // Regular algo: test intersection only with above and bellow segments
      // Here: test with all great circle arcs
      if sweep_line.contains_segment_crossing(event.segment_ref) {
        return true;
      }
      sweep_line.insert(event.segment_ref)
    } else {
      // Regular algo: take above -- take bellow and check if intersects
      sweep_line.remove(event.segment_ref)
    }
  }
  false
}

fn build_events(segments: &[Segment]) -> Vec<Event> {
  segments
    .iter()
    .flat_map(|seg| {
      let cross_pm = seg.cross_primary_meridian();
      let is_v1_left = (seg.v1_lon < seg.v2_lon) != cross_pm;
      vec![
        Event {
          lon: seg.v1_lon,
          is_left: is_v1_left,
          cross_pm,
          segment_ref: seg,
        },
        Event {
          lon: seg.v2_lon,
          is_left: !is_v1_left,
          cross_pm,
          segment_ref: seg,
        },
      ]
    })
    .collect()
}
