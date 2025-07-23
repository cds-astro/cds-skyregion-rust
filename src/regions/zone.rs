use std::ops::Range;

use astrometry::jname::{JNameZone, JZone};
use cdshealpix::nested::{bmoc::BMOC, zone_coverage};

use crate::{
  common::{
    error::SkyRegionError,
    lat_deg2rad, lon_deg2rad,
    math::{HALF_PI, PI, TWICE_PI},
  },
  SkyRegion,
};

#[derive(Debug)]
pub struct Zone {
  /// Minimal longitude (inclusive), in `[0, 2\pi[` radians
  lon_min: f64,
  /// Minimal latitude (inclusive if positive, else exclusive) in `[-\pi/2, \pi/2[` radians
  lat_min: f64,
  /// Maximal longitude (exclusive), in `[0, 2\pi]` radians
  lon_max: f64,
  /// Maximal latitude (exlusive if positive, else inclusive) in `]-\pi/2, \pi/2]` radians
  lat_max: f64,
  /// Tells if the zone cross the primary meridian (in this case, lon_min > lon_max)
  cross_primary_meridian: bool,
}

impl Zone {
  /// Create a new zone from parameters expressed in degrees.
  /// # params
  /// * `lon_min_deg`: longitude of the bottom-left corner, in degrees
  /// * `lat_min_deg`: latitude of the bottom-left corner, in degrees
  /// * `lon_max_deg`: longitude of the top-right corner, in degrees
  /// * `lat_max_deg`: latitude of the top-right corner, in degrees
  pub fn from_deg(
    lon_min_deg: f64,
    lat_min_deg: f64,
    lon_max_deg: f64,
    lat_max_deg: f64,
  ) -> Result<Self, SkyRegionError> {
    let lon_min = lon_deg2rad(lon_min_deg)?;
    let lat_min = lat_deg2rad(lat_min_deg)?;
    let lon_max = lon_deg2rad(lon_max_deg)?;
    let lat_max = lat_deg2rad(lat_max_deg)?;
    Ok(Self::new(lon_min, lat_min, lon_max, lat_max))
  }

  /// # Remarks
  /// * If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
  /// * The north pole is included only if `lon_min == 0 && lat_max == pi/2`
  /// # Panics
  /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
  /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
  /// * `lat_min >= lat_max`
  pub fn new(lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> Zone {
    assert!((0.0..TWICE_PI).contains(&lon_min) && 0.0 < lon_max && lon_max <= TWICE_PI);
    assert!((-HALF_PI..HALF_PI).contains(&lat_min) && -HALF_PI < lat_max && lat_max <= HALF_PI);
    assert!(lat_min <= lat_max);
    Zone {
      lon_min,
      lat_min,
      lon_max,
      lat_max,
      cross_primary_meridian: lon_min > lon_max,
    }
  }

  /// # Input
  /// * `jname`: a jname in one (or a mix of) the following formats
  ///     + `(ACRONYM )JHH+DD(a)`
  ///     + `(ACRONYM )JHHh+DDd(a)`
  ///     + `(ACRONYM )JHHMM+DDMM(a)`
  ///     + `(ACRONYM )JHHMMm+DDMMm(a)`
  ///     + `(ACRONYM )JHHMMSS+DDMMSS(a)`
  ///     + `(ACRONYM )JHHMMSS.S+DDMMSS.S(a)`
  ///     + `(ACRONYM )JDDD.d+DD.dd(a)`
  /// * `rounded`: Should be `false`, set to `true` only if the input JNAME has been
  ///   wrongly computed rounding the coordinates instead of truncating them.
  /// * `epsilon`: an extra size to be added to the zone. Should be `None` except if the JNAME has
  ///   been computed truncating coordinated (e.g. 10^-4 deg) previously rounded (e.g. at 10^-6 deg).
  ///   E.g., if VizieR formated output at 10-6 deg coordinates have been used to compute the JNAME,
  ///   this factor should be 0.5e-6. Example:
  ///     + `RA  =   0.00049996 deg`: RA vlaue
  ///     + `RAF = 000.000500 deg`: output formatting `RA` using `%10.6f`
  ///     + `JNAME JDDD.dddd+...` from `RA` : `J000.0004`
  ///     + `JNAME JDDD.dddd+...` from `RAF`: `J000.0005`
  ///     + `RA range` from `RA `: `[0.0004, 0.0005]`
  ///     + `RA range` from `RAF`: `[0.0005, 0.0006]`
  ///     + `RA range` from `RAF` with 1e-6 correction: `[0.000499, 0.000601]`
  /// # Remark
  /// * We take into account `epsilon` only if `rounded = false`.
  pub fn from_jname(jname: &str, rounded: bool, epsilon: &Option<f64>) -> Result<Zone, String> {
    let jz = jname.parse::<JNameZone>().map_err(|e| e.to_string())?;
    let JZone {
      ra_min_deg,
      ra_max_deg,
      dec_min_deg,
      dec_max_deg,
    } = match (rounded, epsilon) {
      (false, None) => jz.to_zone(),
      (false, Some(eps)) => jz.to_zone_eps(*eps),
      (true, None) => jz.to_zone_round(),
      (true, Some(eps)) => jz.to_zone_round_eps(*eps),
    };
    Ok(Self::new(
      ra_min_deg.to_radians(),
      dec_min_deg.to_radians(),
      ra_max_deg.to_radians(),
      dec_max_deg.to_radians(),
    ))
  }

  pub fn to_bmoc(&self, depth: u8) -> BMOC {
    zone_coverage(
      depth,
      self.lon_min,
      self.lat_min,
      self.lon_max,
      self.lat_max,
    )
  }
}

impl SkyRegion for Zone {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    let mut b = if self.lat_min >= 0.0 || self.lat_min == -HALF_PI {
      self.lat_min <= lat
    } else {
      self.lat_min < lat
    };
    b &= if self.lat_max <= 0.0 || self.lat_max == HALF_PI {
      lat <= self.lat_max
    } else {
      lat < self.lat_max
    };
    b && if self.cross_primary_meridian {
      self.lon_min <= lon || lon < self.lon_max
    } else {
      self.lon_min <= lon && lon < self.lon_max
    }
  }

  fn area(&self) -> f64 {
    // Surface area of a zone in lat = 2*pi*h = 2*pi*(h2-h1) (h = the height of the zone)
    // We find it from the surface area of a cap: 2*pi*h = 2*pi*(1-sin(lat)).
    // We then us the cross-mutliplication:
    //   S = 2*pi*(sin(lat_max) - sin(lat_min)) * (lon_max - lon_min)/2*pi
    //     = (sin(lat_max) - sin(lat_min)) * (lon_max - lon_min)
    let dlon = if self.cross_primary_meridian {
      PI - (self.lon_min - self.lon_max)
    } else {
      self.lon_max - self.lon_min
    };
    (self.lat_max.sin() - self.lat_min.sin()) * dlon
  }

  /*fn characteristic_depth(&self) -> u8 {
    let dlon = self.lon_max - self.lon_min;
    let dlat = self.lat_max - self.lat_min;
    let (a, b) = if dlon > dlat {
      (dlon, dlat)
    }  else {
      (dlat, dlon)
    };
    // The largest between b or a / 8.
    let r = b.max(0.1 * a);
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

#[cfg(test)]
mod tests {
  use super::{super::SkyArea, Zone};

  #[test]
  fn test_jname_to_zone() {
    let jname = "J054214.40-702748.7";
    let lon_deg = 85.560010_f64;
    let lat_deg = -70.463551_f64;
    let zone = Zone::from_jname(jname, false, &None).unwrap();
    println!("zone: {:?}", zone);
    println!("{:?}", zone.sorted_hpx_ranges(10));
    println!(
      "idx10: {}",
      cdshealpix::nested::hash(10, lon_deg.to_radians(), lat_deg.to_radians())
    );
    assert!(zone.contains(lon_deg.to_radians(), lat_deg.to_radians()));
  }
}
