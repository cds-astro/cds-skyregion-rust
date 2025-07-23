use std::ops::Range;

use cdshealpix::{
  best_starting_depth, has_best_starting_depth, nested::ring_coverage_approx_custom, DEPTH_MAX,
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
pub struct Ring {
  /// Center longitude, in `[0, 2\pi[` radians
  lon: f64,
  /// Center latitude in `[-\pi/2, \pi/2]` radians
  lat: f64,
  /// Internal radius, in radians
  r_min: f64,
  /// External radius, in radians
  r_max: f64,
  /// Derived quantity to speed up computations
  cos_lat: f64,
  /// Dereived quantity to speed up computation.
  /// Square Half Segment
  shs_min: f64,
  /// Dereived quantity to speed up computation.
  /// Square Half Segment
  shs_max: f64,
}

impl Ring {
  /// Create a new ring from parameters expressed in degrees.
  /// # params
  /// * `lon_deg`: longitude of the center of the ring, in degrees
  /// * `lat_deg`: latitude of the center of the ring, in degrees
  /// * `r_min_deg`: internal radius of the ring, in degrees
  /// * `r_max_deg`: external radius of the ring, in degrees
  pub fn from_deg(
    lon_deg: f64,
    lat_deg: f64,
    r_min_deg: f64,
    r_max_deg: f64,
  ) -> Result<Self, SkyRegionError> {
    let lon = lon_deg2rad(lon_deg)?;
    let lat = lat_deg2rad(lat_deg)?;
    let r_min = r_min_deg.to_radians();
    let r_max = r_max_deg.to_radians();
    if r_min <= 0.0 || r_max <= 0.0 || PI <= r_min || PI <= r_max || r_min >= r_max {
      Err(SkyRegionError::UnexpectedRingRadii {
        r_max_deg,
        r_min_deg,
      })
    } else {
      Ok(Self::new(lon, lat, r_min, r_max))
    }
  }

  /// # Panics
  /// * if `lon` not in `[0, 2\pi[`
  /// * if `lat` not in `[-\pi/2, \pi/2[`
  /// * if `r` not in `]0, \pi[`
  pub fn new(lon: f64, lat: f64, r_min: f64, r_max: f64) -> Self {
    assert!((0.0..TWICE_PI).contains(&lon));
    assert!((-HALF_PI..=HALF_PI).contains(&lat));
    assert!(0.0 < r_min && r_min < r_max && r_max < PI);
    Self {
      lon,
      lat,
      r_min,
      r_max,
      cos_lat: lat.cos(),
      shs_min: squared_half_segment_from_dist(r_min),
      shs_max: squared_half_segment_from_dist(r_max),
    }
  }

  pub fn shs_is_ok(&self, shs: f64) -> bool {
    self.shs_min < shs && shs <= self.shs_max
  }
}

impl SkyRegion for Ring {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    (self.lat - lat).abs() < self.r_max
      && self.shs_is_ok(squared_half_segment(
        lon - self.lon,
        lat - self.lat,
        self.cos_lat,
        lat.cos(),
      ))
  }

  fn area(&self) -> f64 {
    if self.r_max > 1.0e-4 {
      // => ~20 arcsec
      // Area of a cone: S = 2 * pi * (1 - cos(theta))
      // 2 * pi * (1 - cos(theta_max)) - 2 * pi * (1 - cos(theta_min))
      // =>  2 * pi (cos(theta_min) - cos(theta_max))
      TWICE_PI * (self.r_min.cos() - self.r_max.cos())
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
      PI * (self.r_max.pow2() - self.r_min.pow2())
    }
  }

  fn characteristic_depth(&self) -> u8 {
    let r_mean = 0.5 * (self.r_max + self.r_min);
    if has_best_starting_depth(r_mean) {
      (best_starting_depth(r_mean) + 2).min(DEPTH_MAX)
    } else {
      2
    }
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    ring_coverage_approx_custom(depth, 3, self.lon, self.lat, self.r_min, self.r_max)
      .to_flagged_ranges()
  }
}
