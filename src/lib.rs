use std::ops::Range;

use cdshealpix::{best_starting_depth, has_best_starting_depth, DEPTH_MAX};

pub mod common;
pub mod regions;

use self::common::math::{PI, TWICE_PI};

/// Define a SkyRegion with a method to test if a given point is inside or not.
pub trait SkyRegion: Send + Sync {
  /// Returns `true` if the given position is inside the region.
  /// # Params
  /// * `lon`: longitude, in `[0, 2\pi[` radians
  /// * `lat`: latitude in `[-\pi/2, \pi/2]` radians
  fn contains(&self, lon: f64, lat: f64) -> bool;

  /// Returns the surface area, may be an approximation.
  fn area(&self) -> f64;

  /// Provides the characteristic HEALPix depth of the cell overlapping the region.
  /// The returned value is a compromise between the total number of cells and out-of-area coverage.
  fn characteristic_depth(&self) -> u8 {
    let area = self.area();
    let eq_cone_radius = if area < 3.0e-8 {
      // => cone of radius ~20.1 arcsec
      // Euclidean approximation
      (area / PI).sqrt()
    } else {
      // Exact spherical computation
      (1.0 - (area / TWICE_PI)).acos()
    };
    if has_best_starting_depth(eq_cone_radius) {
      (best_starting_depth(eq_cone_radius) + 2).min(DEPTH_MAX) // if a cone => min 16 cells
    } else {
      2
    }
  }

  /// Ordered list of HEALPix cell ranges of given `depth` overlapping the region.
  /// The associated boolean is `true` if all cells in the range are for sure inside the region.
  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)>;
  // Returns an IntoIter ? (add a parameter type)
  // Add a sorted_contains(&mut self, ) ? // => increasing healpix order

  // For xmatches!
  // fn sorted_extended_hpx_ranges(&self, depth: u8) -> List<Range>;
}

impl<S: SkyRegion + ?Sized> SkyRegion for Box<S> {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    (**self).contains(lon, lat)
  }
  fn area(&self) -> f64 {
    (**self).area()
  }
  fn characteristic_depth(&self) -> u8 {
    (**self).characteristic_depth()
  }
  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    (**self).sorted_hpx_ranges(depth)
  }
}

/// Defines a process depending on a `SkyRegion`.
/// This trait is made to allow for monomorphisation (at the cost of higher compilation time and binaries).
/// If the binary is too large, or the compilation time too long, try `SkyRegionDynProcess` instead.
pub trait SkyRegionProcess {
  type Output;
  type Error;

  fn exec<S: SkyRegion>(self, region: S) -> Result<Self::Output, Self::Error>;
}

/// Defines an action depending on a dyn `SkyRegion`.
/// This trait is made to avoid monomorphisation.
pub trait SkyRegionDynProcess {
  type Output;
  type Error;

  fn exec_dyn(self, region: Box<dyn SkyRegion>) -> Result<Self::Output, Self::Error>;
}

#[cfg(test)]
mod tests {
  use super::{
    SkyRegion,
    regions::cone::Cone
  };

  #[test]
  fn test_cone() {
    let cone = Cone::new(0.0, 0.0, 0.5);
    assert!(cone.contains(0.01, 0.01))
  }
}
