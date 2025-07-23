use std::ops::Range;

use cdshealpix::{n_hash, nested::hash};

use crate::{common::math::PI, SkyRegion};

#[derive(Debug)]
pub struct HpxRange {
  /// HealpixDepth
  depth: u8,
  /// HEALPix cell hash value (i.e. cell number at the given depth)
  range: Range<u64>,
}

impl HpxRange {
  /// # Panics
  /// * if `depth` not in `[0, 30[`
  /// * if `hash` not in `[0, 12 * 2^(2*depth)[`
  pub fn new(depth: u8, range: Range<u64>) -> Self {
    assert!(range.end < n_hash(depth));
    Self { depth, range }
  }
}

impl SkyRegion for HpxRange {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    let h = hash(self.depth, lon, lat);
    self.range.start <= h && h < self.range.end
  }

  fn area(&self) -> f64 {
    (4.0 * PI) * ((self.range.end - self.range.start) as f64) / (n_hash(self.depth) as f64)
  }

  fn characteristic_depth(&self) -> u8 {
    self.depth
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    if depth == self.depth {
      vec![(self.range.clone(), true)]
    } else if depth < self.depth {
      let twice_dd = (self.depth - depth) << 1;
      let start = self.range.start >> twice_dd;
      let end = self.range.end >> twice_dd;
      if (end << twice_dd) == self.range.end {
        vec![(
          Range { start, end },
          (start << twice_dd) == self.range.start,
        )]
      } else {
        vec![(
          Range {
            start,
            end: end + 1,
          },
          false,
        )]
      }
    } else {
      let twice_dd = (depth - self.depth) << 1;
      let start = self.range.start << twice_dd;
      let end = self.range.end << twice_dd;
      vec![(Range { start, end }, true)]
    }
  }
}
