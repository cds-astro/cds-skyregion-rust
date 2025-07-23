use std::ops::Range;

use cdshealpix::{
  best_starting_depth, has_best_starting_depth,
  nested::{bmoc::BMOC, cone_coverage_approx_custom},
  DEPTH_MAX,
};

use crate::{
  common::{
    error::SkyRegionError,
    lat_deg2rad, lon_deg2rad,
    math::{
      squared_half_segment, squared_half_segment_from_dist, Customf64, HALF_PI, PI, TWICE_PI,
    },
  },
  SkyRegion,
};

#[derive(Debug)]
pub struct Cone {
  /// Center longitude, in `[0, 2\pi[` radians
  lon: f64,
  /// Center latitude in `[-\pi/2, \pi/2]` radians
  lat: f64,
  /// Radius, in radians
  r: f64,
  /// Derived quantity to speed up computations: cosine of the latitude.
  cos_lat: f64,
  /// Derived quantity to speed up computation: radius expressed as a Squared Half Segment.
  shs: f64,
}

impl Cone {
  /// Create a new cone from parameters expressed in degrees.
  /// # params
  /// * `lon_deg`: longitude of the center of the cone, in degrees
  /// * `lat_deg`: latitude of the center of the cone, in degrees
  /// * `r_deg`: radius of the cone, in degrees
  pub fn from_deg(lon_deg: f64, lat_deg: f64, r_deg: f64) -> Result<Self, SkyRegionError> {
    let lon = lon_deg2rad(lon_deg)?;
    let lat = lat_deg2rad(lat_deg)?;
    let r = r_deg.to_radians();
    if r <= 0.0 || PI <= r {
      Err(SkyRegionError::UnexpectedConeRadius { r_deg })
    } else {
      Ok(Self::new(lon, lat, r))
    }
  }

  /// Quick creation of a new cone, input parameters are in radians, and **must** have
  /// already been checked to avoid a `panic!`.
  ///
  /// # params
  /// * `lon`: longitude of the center of the cone, in radians
  /// * `lat`: latitude of the center of the cone, in radians
  /// * `r`: radius of the cone, in radians
  ///
  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2[`
  /// * if `r` not in `]0, \pi[`
  pub fn new(lon: f64, lat: f64, r: f64) -> Self {
    assert!((0.0..TWICE_PI).contains(&lon));
    assert!((-HALF_PI..=HALF_PI).contains(&lat));
    assert!(0.0 < r && r < PI);
    Self {
      lon,
      lat,
      r,
      cos_lat: lat.cos(),
      shs: squared_half_segment_from_dist(r),
    }
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    cone_coverage_approx_custom(depth, 3, self.lon, self.lat, self.r)
  }
}

impl SkyRegion for Cone {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    (self.lat - lat).abs() < self.r
      && squared_half_segment(lon - self.lon, lat - self.lat, self.cos_lat, lat.cos()) < self.shs
  }

  fn area(&self) -> f64 {
    if self.r > 1.0e-4 {
      // => ~20 arcsec
      // Area of a cone: S = 2 * pi * (1 - cos(theta))
      TWICE_PI * (1.0 - self.r.cos())
    } else {
      // Eucliean distance between center and border of the cone
      //   2 * sin(theta / 2)
      // => S = 2 * pi * [ (2*sin(theta/2))^2 - sin(theta)^2 ]
      // =>   = 2 * pi  * (2*sin(theta/2) + sin(theta)) * (2*sin(theta/2) - sin(theta))
      // If we need faster computation, we could use the sin Taylor series.
      /*let sin_t = self.r.sin();
      let twice_sin_half_t = 2.0 * (*self.r * 0.5).sin();
      2 * pi * (twice_sin_half_t + sin_t) * (twice_sin_half_t - sin_t)*/
      // Euclidean approximation is ok
      PI * self.r.pow2()
    }
  }

  fn characteristic_depth(&self) -> u8 {
    if has_best_starting_depth(self.r) {
      (best_starting_depth(self.r) + 2).min(DEPTH_MAX)
    } else {
      2
    }
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    self.to_bmoc(depth).to_flagged_ranges()
  }
}

#[cfg(test)]
mod tests {

  use crate::SkyRegion;
  use super::*;

  #[test]
  fn test_cone() {
    let lon_deg = 0.035323_f64;
    let lat_deg = 36.585958_f64;
    let r_deg: f64 = 0.005;
    let cone = Cone::new(
      lon_deg.to_radians(),
      lat_deg.to_radians(),
      r_deg.to_radians(),
    );
    println!("SortedHpxRanges: {:?}", cone.sorted_hpx_ranges(6));
    assert!(true);
  }

  #[test]
  fn test_cone_2() {
    let lon_deg = 8.401978_f64;
    let lat_deg = 84.675171_f64;
    let r_deg: f64 = 0.0008;
    let cone = Cone::new(
      lon_deg.to_radians(),
      lat_deg.to_radians(),
      r_deg.to_radians(),
    );
    println!("SortedHpxRanges: {:?}", cone.sorted_hpx_ranges(6));

    let r_deg: f64 = 0.0004;
    let cone = Cone::new(
      lon_deg.to_radians(),
      lat_deg.to_radians(),
      r_deg.to_radians(),
    );
    println!("SortedHpxRanges: {:?}", cone.sorted_hpx_ranges(6));

    assert!(true);
  }

  #[test]
  fn test_cone_3() {
    let lon_deg = 314.990401_f64;
    let lat_deg = -0.008369_f64;
    let r_deg: f64 = 0.008;
    let cone = Cone::new(
      lon_deg.to_radians(),
      lat_deg.to_radians(),
      r_deg.to_radians(),
    );
    println!("SortedHpxRanges: {:?}", cone.sorted_hpx_ranges(6));

    assert!(true);
  }
}
