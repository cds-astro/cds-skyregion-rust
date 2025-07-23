use std::ops::Range;

use cdshealpix::{n_hash, nested::hash};

use crate::{common::math::PI, SkyRegion};

#[derive(Debug)]
pub struct HpxCell {
  /// HealpixDepth
  depth: u8,
  /// HEALPix cell hash value (i.e. cell number at the given depth)
  hash: u64,
}

impl HpxCell {
  /// # Panics
  /// * if `depth` not in `[0, 30[`
  /// * if `hash` not in `[0, 12 * 2^(2*depth)[`
  pub fn new(depth: u8, hash: u64) -> Self {
    assert!(hash < n_hash(depth));
    Self { depth, hash }
  }
}

impl SkyRegion for HpxCell {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    self.hash == hash(self.depth, lon, lat)
  }

  fn area(&self) -> f64 {
    (4.0 * PI) / (n_hash(self.depth) as f64)
  }

  fn characteristic_depth(&self) -> u8 {
    self.depth
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    if depth == self.depth {
      vec![(
        Range {
          start: self.hash,
          end: self.hash + 1,
        },
        true,
      )]
    } else if depth < self.depth {
      let start = self.hash >> ((self.depth - depth) << 1);
      vec![(
        Range {
          start,
          end: start + 1,
        },
        false,
      )]
    } else {
      let twice_dd = (depth - self.depth) << 1;
      let start = self.hash << twice_dd;
      let end = (self.hash + 1) << twice_dd;
      vec![(Range { start, end }, true)]
    }
  }
}
