use std::ops::Range;

use cdshealpix::nested::{bmoc::BMOC, elliptical_cone_coverage_custom};

use crate::{
  common::{
    error::SkyRegionError,
    lat_deg2rad, lon_deg2rad,
    math::{
      sphe_dist, squared_half_segment, xyz_to_lonlat, RefToLocalRotMatrix, HALF_PI, PI, TWICE_PI,
    },
  },
  SkyRegion,
};

#[derive(Debug)]
pub struct EllipticalCone {
  /// Center longitude, in `[0, 2\pi[` radians
  lon: f64,
  /// Center latitude in `[-\pi/2, \pi/2]` radians
  lat: f64,
  /// Semi-major axis, in radians
  a: f64,
  /// Semi-minor axis, in radians
  b: f64,
  /// Positional angle (east-of-north), in `[0, \pi[` radians
  pa: f64,
  // Derived quantities to speed up computations
  // Global to local frame rotation matrix
  // to_local_matrix: RefToLocalRotMatrix,
  // Euclidean coordinates of the first focus
  // x_f0: f64, y_f0: f64, z_f0: f64,
  /// Equatorial coordinates of the first focus
  lon_f0: f64,
  lat_f0: f64,
  /// To speed up computatons
  cos_lat_f0: f64,
  // Euclidean coordinates of the first focus
  // x_f1: f64, y_f1: f64, z_f1: f64,
  /// Equatorial coordinates of the first focus
  lon_f1: f64,
  lat_f1: f64,
  /// To speed up computatons
  cos_lat_f1: f64,
}

impl EllipticalCone {
  /// Create a new elliptical cone from parameters expressed in degrees.
  /// # params
  /// * `lon_deg`: longitude of the center of the elliptical cone, in degrees
  /// * `lat_deg`: latitude of the center of the elliptical cone, in degrees
  /// * `a_deg`: semi-major axis, in degrees
  /// * `b_deg`: semi-minor axis, in degrees
  /// * `pa_deg`: position angle, in degrees
  pub fn from_deg(
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
      Err(SkyRegionError::UnexpectedEllipseA { a_deg })
    } else if b <= 0.0 || a <= b {
      Err(SkyRegionError::UnexpectedEllipseB { b_deg })
    } else if pa <= 0.0 || PI <= pa {
      Err(SkyRegionError::UnexpectedEllipsePA { pa_deg })
    } else {
      Ok(Self::new(lon, lat, a, b, pa))
    }
  }

  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2]`
  /// * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
  /// * if `b` not in `]0, a[`
  /// * if `pa` not in `[0, \pi[`
  pub fn new(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> Self {
    assert!((0.0..TWICE_PI).contains(&lon));
    assert!((-HALF_PI..=HALF_PI).contains(&lat));
    assert!(0.0 < a && a <= HALF_PI);
    assert!(0.0 < b && b < a);
    assert!(0.0 < pa && pa < PI);
    let (sin_pa, cos_pa) = pa.sin_cos();
    let (sina, cosa) = a.sin_cos();
    let (sinb, cosb) = b.sin_cos();
    // Compute foci coordinates.
    // * a: semi-major angular distance
    // * b: semi-minor angular distance
    // * g: semi-foci distance
    //         B           -
    //         |           | b
    // A---F---O---F---A'  -
    //     a     g
    // |-------|---|
    //   OA = a; OB = b; FO = g
    // Let's note f0 = OF and f1 = OF'
    //   f0 + f1 = 2g; FM + MF' = 2a
    // We consider the spherical triangle FOB
    // Using the spherical pythagorean theorem (considering FM + MF' = 2a => FB = a):
    //    cos(a) = cos(g)cos(b)    (in Euclidean: g^2 = a^2 - b^2).
    // => cos(g) = cos(a) / cos(b)
    // => sin(g) = sqrt(1-cos^2(g)) = sqrt(1 - (cos(a)/cos(b))^2)
    // The foci cordinates are (X=+-sin(g), Y=0) in the canonical frame
    // (frame in which the major axis is the x-axis)
    // Taking the cos of f0 + f1 = 2a, we also get
    //   cos(f0)cos(f1) - sin(f0)sin(f1) = cos(2a)
    let cosg = cosa / cosb;
    // For small angles, the following formula is more precise than sqrt(1 - sin^2)
    // sing = sqrt(1 - cosg^2) = sqrt((cob^2 - cosa^2) / cosb^2) = sqrt(sina^2 - sinb^2) / cosb
    //      = sqrt( (sina + sinb) * (sina - sinb) ) / cosb
    let sing = ((sina + sinb) * (sina - sinb)).sqrt() / cosb;
    // F0 in the local frame
    //   x = sqrt(1 - (y^2 + z^2) = cosg
    //   Rotation of (y=sing, z=0):
    //   y = y * cos(pi/2 - pa) - z * sin(pi/2 - pa) = sin(g) * sin(pa)
    //   z = y * sin(pi/2 - pa) + z * cos(pi/2 - pa) = sin(g) * cos(pa)
    let (x_f0, y_f0, z_f0) = (cosg, sing * sin_pa, sing * cos_pa);
    // F1 in the local frame
    //   x = sqrt(1 - (y^2 + z^2) = cosg
    //   (y, z) = rotateEllipse(new double[]{-sing, 0}
    let (x_f1, y_f1, z_f1) = (cosg, -y_f0, -z_f0);
    // Compute reference to local frame rotation matrix
    let rotation = RefToLocalRotMatrix::from_center(lon, lat);
    // Rotate to global frame
    let (x_f0, y_f0, z_f0) = rotation.to_global_xyz(x_f0, y_f0, z_f0);
    let (x_f1, y_f1, z_f1) = rotation.to_global_xyz(x_f1, y_f1, z_f1);
    let (lon_f0, lat_f0) = xyz_to_lonlat(x_f0, y_f0, z_f0);
    let (lon_f1, lat_f1) = xyz_to_lonlat(x_f1, y_f1, z_f1);
    Self {
      lon,
      lat,
      a,
      b,
      pa,
      // to_local_matrix: rotation,
      // x_f0, y_f0, z_f0,
      lon_f0,
      lat_f0,
      cos_lat_f0: lat_f0.cos(),
      // x_f1, y_f1, z_f1,
      lon_f1,
      lat_f1,
      cos_lat_f1: lat_f0.cos(),
    }
  }

  fn dist_to_f0_and_f1(&self, lon: f64, lat: f64) -> (f64, f64) {
    let cos_lat = lat.cos();
    (
      sphe_dist(squared_half_segment(
        lon - self.lon_f0,
        lat - self.lat_f0,
        self.cos_lat_f0,
        cos_lat,
      )),
      sphe_dist(squared_half_segment(
        lon - self.lon_f1,
        lat - self.lat_f1,
        self.cos_lat_f1,
        cos_lat,
      )),
    )
  }

  fn half_sum_dist_to_foci(&self, lon: f64, lat: f64) -> f64 {
    let (d_f1, d_f2) = self.dist_to_f0_and_f1(lon, lat);
    0.5 * (d_f1 + d_f2)
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    elliptical_cone_coverage_custom(depth, 3, self.lon, self.lat, self.a, self.b, self.pa)
  }
}

impl SkyRegion for EllipticalCone {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    (self.lon - lon).abs() < self.a && self.half_sum_dist_to_foci(lon, lat) < self.a
  }

  fn area(&self) -> f64 {
    if self.a > 1.0e-4 {
      // => ~20 arcsec
      // integral sqrt((cos^2(x)*tan^2(b) - tan^2(b)/tan^2(a) * sin^2(x))/(1 + (cos^2(x)*tan^2(b) - tan^2(b)/tan^2(a) * sin^2(x)))) from 0 to a
      eprintln!("WARNING: the surface area of the elliptical cone is an approximation!");
      PI * self.a.sin() * self.b.sin()
    } else {
      // Euclidean approximation
      PI * self.a * self.b
    }
  }

  /*fn characteristic_depth(&self) -> u8 {
    // The largest between b or a / 8.
    let r = self.b.max(0.1 * self.a);
    if has_best_starting_depth(r) {
      (best_starting_depth(r) + 2).min(DEPTH_MAX)
    } else {
      0
    }
  }*/

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    self.to_bmoc(depth).to_flagged_ranges()
  }
}
