use std::{ops::Range, ptr::slice_from_raw_parts};

use cdshealpix::{n_hash, nested::external_edge, nested::hash};

use crate::{common::math::PI, SkyRegion};

// I did this implementation for the sake of genericity.
// But ideally, iterating over the (ordered) ranges (on HEALpix ordered sources),
// I could have tested (in the contains() method) only the 'current' range.

#[derive(Debug)]
pub struct HpxRanges {
  /// HealpixDepth
  depth: u8,
  /// HEALPix cell hash value (i.e. cell number at the given depth)
  ranges: Vec<Range<u64>>,
}

impl HpxRanges {
  /// # Panics
  /// * if `depth` not in `[0, 30[`
  /// * if greatest value not in `[0, 12 * 2^(2*depth)]`
  /// * if input ranges are not ordered or are overlapping
  pub fn new(depth: u8, ranges: Vec<Range<u64>>) -> Self {
    check_sorted_and_non_overlapping(&ranges);
    assert!(
      ranges.last().map(|r| r.end).unwrap_or(0) <= n_hash(depth),
      "{} > {}",
      ranges.last().map(|r| r.end).unwrap_or(0),
      n_hash(depth)
    );
    Self { depth, ranges }
  }

  /// Returns ranges from an Healpix cell plus its external border of given depth.
  pub fn from_cell_with_border(
    cell_depth: u8,
    cell_hash: u64,
    external_border_depth: u8,
  ) -> HpxRanges {
    assert!(external_border_depth > cell_depth);
    let delta_depth = external_border_depth - cell_depth;
    // We can add them without merging them since we are sure the HEALPix index of cells
    // on external borders are distinct (Z-order curve + delte_depth > 0)

    /*let mut ranges: Vec<Range<u64>> = external_edge_sorted(cell_depth, cell_hash, delta_depth)
      .into_iter()
      .map(|h| *h..*h + 1)
      .collect();
    let shift = delta_depth << 1;
    let cell_range = cell_hash << shift..(cell_hash + 1) << shift;
    let i = match ranges.binary_search_by_key(&cell_range.start, |range| range.end) {
      Ok(i) => {
        ranges[i].end = cell_range.end;
        i
      }
      Err(i) => {
        ranges.insert(i, cell_range);
        i
      }
    };
    Add it, then re-merge!!!
    if i + 1 < ranges.len() && ranges[i].end == ranges[i + 1].start {
      ranges[i].end = ranges[i + 1].end;
      ranges.remove(i + 1);
    }*/

    let mut ranges: Vec<Range<u64>> = external_edge(cell_depth, cell_hash, delta_depth)
      .iter()
      .map(|h| *h..*h + 1)
      .collect();
    let shift = delta_depth << 1;
    ranges.push(cell_hash << shift..(cell_hash + 1) << shift);
    ranges.sort_by_key(|range| range.start);

    let mut merged_ranges: Vec<Range<u64>> = Vec::with_capacity(ranges.len());
    for new_range in ranges {
      if let Some(Range { start: _, end }) = merged_ranges.last_mut() {
        // merge overlapping ranges
        if new_range.start <= *end {
          if *end < new_range.end {
            *end = new_range.end
          }
        } else {
          merged_ranges.push(new_range);
        }
      } else {
        merged_ranges.push(new_range);
      }
    }

    Self::new(external_border_depth, merged_ranges)
  }
}

fn check_sorted_and_non_overlapping(ranges: &[Range<u64>]) {
  for (r1, r2) in ranges.iter().zip(ranges.iter().skip(1)) {
    assert!(r1.end < r2.start)
  }
}

impl SkyRegion for HpxRanges {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    // By accepting a (&mut self), we could cache the last hash with its results
    // and/or keep track of the progression in the ranges to avoid having to look
    // on the full table (but we have to be sure then the interrogation is made
    // using HEALPix ordering).
    let h = hash(self.depth, lon, lat);
    // See the vec of ranges as a vec of elements
    let len = self.ranges.len() << 1;
    let ptr = self.ranges.as_ptr();
    let result: &[u64] = unsafe { &*slice_from_raw_parts(ptr as *const u64, len) };
    // Performs a binary search in it
    match result.binary_search(&h) {
      Ok(i) => i & 1 == 0,  // index must be even (lower bound of a range)
      Err(i) => i & 1 == 1, // index must be odd (inside a range)
    }
  }

  fn area(&self) -> f64 {
    (4.0 * PI) * self.ranges.iter().map(|r| r.end - r.start).sum::<u64>() as f64
      / (n_hash(self.depth) as f64)
  }

  fn characteristic_depth(&self) -> u8 {
    self.depth
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    if depth == self.depth {
      self.ranges.iter().map(|r| (r.clone(), true)).collect()
    } else if depth < self.depth {
      let twice_dd = (self.depth - depth) << 1;
      let (mut range_vec, last_range) = self
        .ranges
        .iter()
        .map(|range| {
          // degrade the range
          let start = range.start >> twice_dd;
          let end = range.end >> twice_dd;
          if (end << twice_dd) == range.end {
            (Range { start, end }, (start << twice_dd) == range.start)
          } else {
            (
              Range {
                start,
                end: end + 1,
              },
              false,
            )
          }
        })
        .fold(
          (Vec::with_capacity(self.ranges.len()), None),
          |(mut ranges, acc_range_flag_opt), (cur_range, cur_flag)| {
            // Possible fusion of ranges
            match acc_range_flag_opt {
              Some::<(Range<u64>, bool)>((mut acc_range, acc_flag)) => {
                if acc_range.end <= cur_range.start {
                  // can be =, not supposed to be <
                  // Range fusion (=> flag = false)
                  acc_range.end = cur_range.end;
                  (ranges, Some((acc_range, false)))
                } else {
                  ranges.push((acc_range, acc_flag));
                  (ranges, Some((cur_range, cur_flag)))
                }
              }
              None => (ranges, Some((cur_range, cur_flag))),
            }
          },
        );
      if let Some(last_range) = last_range {
        range_vec.push(last_range);
      }
      range_vec
    } else {
      let twice_dd = (depth - self.depth) << 1;
      self
        .ranges
        .iter()
        .map(|range| {
          let start = range.start << twice_dd;
          let end = range.end << twice_dd;
          (Range { start, end }, true)
        })
        .collect()
    }
  }
}
