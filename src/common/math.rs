pub const HALF_PI: f64 = 0.5 * std::f64::consts::PI;
pub const PI: f64 = std::f64::consts::PI;
pub const TWICE_PI: f64 = 2.0 * std::f64::consts::PI;

/// Equatorial coordinate `(lon, lat)` with
/// * `lon`: longitude, in `[0, 2\pi[` radians
/// * `lat`: latitude, in `[-\pi/2, \pi/2]` radians
pub type Coo = (f64, f64);
/// Vector in a 3-dimensional Euclidean space.
pub type Vec3 = (f64, f64, f64);
/// Unit vector in a 3-dimensional Euclidean space, i.e. the norm = 1
/// (i.e. position on the unit sphere).
pub type UnitVec3 = (f64, f64, f64);

/*
/// Both the equatorial coordinate and its associated unit vector
/// (packed together). Purpose: use the version best fitting a need
/// without having to convert from one to another (high CPU cost due to
/// trigonometric functions).
type CooAndUnitVec3 = (Coo, Vec3);
*/

pub trait Customf64 {
  fn pow2(self) -> f64;
  fn twice(self) -> f64;
  fn half(self) -> f64;
}

impl Customf64 for f64 {
  /// Returns x^2
  fn pow2(self) -> f64 {
    self * self
  }
  /// Returns 2 * x
  fn twice(self) -> f64 {
    2.0 * self // self + self (I hope the compiler know the simple shift bit to be used for x2)
  }
  /// Returns x / 2
  fn half(self) -> f64 {
    0.5 * self
  }
}

pub fn haversine_dist(p1_lon: f64, p1_lat: f64, p2_lon: f64, p2_lat: f64) -> f64 {
  let shs = squared_half_segment(p2_lon - p1_lon, p2_lat - p1_lat, p1_lat.cos(), p2_lat.cos());
  sphe_dist(shs)
}

/// Returns the angular distance corresponding to the given squared half great-circle arc segment.
pub fn sphe_dist(squared_half_segment: f64) -> f64 {
  squared_half_segment.sqrt().asin().twice()
}

/// Returns `(s/2)^2` with `s` the segment (i.e. the Euclidean distance) between
/// the two given points  `P1` and `P2` on the unit-sphere.
/// We recall that `s = 2 sin(ad/2)` with `ad` the angular distance between the two points.
/// # Input
/// - `dlon` the longitude difference, i.e. (P2.lon - P1.lon), in radians
/// - `dlat` the latitude difference, i.e. (P2.lat - P1.lat), in radians
/// - `cos_lat1` cosine of the latitude of the first point
/// - `cos_lat2` cosine of the latitude of the second point
pub fn squared_half_segment(dlon: f64, dlat: f64, cos_lat1: f64, cos_lat2: f64) -> f64 {
  dlat.half().sin().pow2() + cos_lat1 * cos_lat2 * dlon.half().sin().pow2()
}

/// Returns `(s/2)^2` with `s` the segment (i.e. the Euclidean distance) between
/// the two given points  `P1` and `P2` on the unit-sphere,
/// from the given angular distance.
/// We recall that `s = 2 sin(ad/2)` with `ad` the angular distance.
/// `(s/2)^2 = (2 sin(ad/2) / 2)^2 = sin(ad/2)^2`
pub fn squared_half_segment_from_dist(sphe_dist: f64) -> f64 {
  sphe_dist.half().asin().pow2()
}

// TODO: remove this code and use the one in cdshealpix!!!
/// Components of the 3x3 rotation matrix transforming a vector into the
/// reference frame to a vector into the local frame (i.e. the frame in
/// which the position of the projection center is (1, 0, 0).
/// Remark:
/// * `r22 =  cos(lon)`
/// * `r21 = -sin(lon)`
/// * `r33 =  cos(lat)`
/// * `r13 =  sin(lat)`
#[derive(Debug)]
pub struct RefToLocalRotMatrix {
  r11: f64,
  r12: f64,
  r13: f64,
  r21: f64,
  r22: f64,
  r23: f64,
  r31: f64,
  r32: f64,
  r33: f64,
}

impl RefToLocalRotMatrix {
  /// Compute the reference to local rotation matrix from a projection center
  pub fn from_center(lon: f64, lat: f64) -> RefToLocalRotMatrix {
    let (sa, ca) = lon.sin_cos(); // ca, sa stands for cos(alpha), sin(alpha)
    let (sd, cd) = lat.sin_cos(); // cd, sd stands for cos(delta), sin(delta)
    RefToLocalRotMatrix {
      r11: ca * cd,
      r12: sa * cd,
      r13: sd,
      r21: -sa,
      r22: ca,
      r23: 0.0,
      r31: -ca * sd,
      r32: -sa * sd,
      r33: cd,
    }
  }

  /// Transform local to global (or reference) coordinates, by
  /// rotating the input (local) vector using the transpose of the rotation matrix.
  pub fn to_global_xyz(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    (
      self.r11 * x + self.r21 * y + self.r31 * z,
      self.r12 * x + self.r22 * y + self.r32 * z,
      self.r13 * x + self.r23 * y + self.r33 * z,
    )
  }

  /// Transform local to global (or reference) coordinates, by
  /// rotating the input (local) vector using the transpose of the rotation matrix.
  pub fn to_global_coo(&self, x: f64, y: f64, z: f64) -> (f64, f64) {
    let (x, y, z) = self.to_global_xyz(x, y, z);
    xyz_to_lonlat(x, y, z)
  }
}

pub fn lonlat_to_xyz(lon: f64, lat: f64) -> (f64, f64, f64) {
  let (sa, ca) = lon.sin_cos(); // ca, sa stands for cos(alpha), sin(alpha)
  let (sd, cd) = lat.sin_cos(); // cd, sd stands for cos(delta), sin(delta)
  (ca * cd, sa * cd, sd)
}

pub fn xyz_to_lonlat(x: f64, y: f64, z: f64) -> (f64, f64) {
  // Length of the projection on the xy plane
  let r2 = x.pow2() + y.pow2();
  // Latitude in [-pi/2, pi/2] (ok, since cos always positive here)
  let lat = z.atan2(r2.sqrt());
  // Compute the longitude in [-pi, pi]
  let r2 = y.atan2(x);
  // Conforms to convention: Longitude in [0, 2*PI]
  let lon = if r2 < 0.0 { TWICE_PI + r2 } else { r2 };
  (lon, lat)
}

pub fn squared_norm(v: &Vec3) -> f64 {
  v.0 * v.0 + v.1 * v.1 + v.2 * v.2
}

pub fn time(cte: f64, v: &Vec3) -> Vec3 {
  (cte * v.0, cte * v.1, cte * v.2)
}

pub fn minus(lhs: &Vec3, rhs: &Vec3) -> Vec3 {
  (lhs.0 - rhs.0, lhs.1 - rhs.1, lhs.2 - rhs.2)
}
pub fn norm(v: &Vec3) -> f64 {
  squared_norm(v).sqrt()
}

pub fn normalized(v: &Vec3) -> UnitVec3 {
  time(1.0 / norm(v), v)
}

pub fn scalar_product(lhs: &(f64, f64, f64), rhs: &(f64, f64, f64)) -> f64 {
  lhs.0 * rhs.0 + lhs.1 * rhs.1 + lhs.2 * rhs.2
}

pub fn cross_product(lhs: &(f64, f64, f64), rhs: &(f64, f64, f64)) -> (f64, f64, f64) {
  (
    lhs.1 * rhs.2 - lhs.2 * rhs.1,
    lhs.2 * rhs.0 - lhs.0 * rhs.2,
    lhs.0 * rhs.1 - lhs.1 * rhs.0,
  )
}
