// MULTICONE
// Attach to each cell the list of cone (of same radius) overlapping the cell.
// The main interest is not to return duplicates in case of overlapping cones.
//
// We could generalize to "multi area", but handling the difficulty would be to
// deal with the heterogeneity of HEALPix depth.

// Make the list of cells/ranges
// and a list range/cell <--> List of cone ?

use std::ops::Range;

#[cfg(feature = "rayon")]
use rayon::{prelude::*, ThreadPool, ThreadPoolBuilder};

use cdshealpix::{
  best_starting_depth, has_best_starting_depth,
  nested::{bmoc::BMOCBuilderFixedDepth, cone_coverage_approx_custom, hash},
  DEPTH_MAX,
};
use moc::{
  moc::{
    range::{CellSelection, RangeMOC},
    RangeMOCIntoIterator,
  },
  qty::{Hpx, MocQty},
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

// Assume ASCII position = 10 + 1 (separator) + 10 + 1 (\n) characters
// 100 MB = 4.7 million positions
//  40 MB = 2.0 million positions
// HDD if 3.00 ms per pos => 100 minutes !!
// SSD if 0.15 ms per pos =>   5 minutes !!
/// Maximum number of positions that can be processed at the same time.
/// We limit because of the MOC intersection computation time.
/// Multicone is not a xmatch!!
const MAX_N_POSITIONS: u32 = 1_000_000; // => ~20MB if a row = 22 charaters
                                        // =>  50 min with HDD (if one seek per row)
                                        // => 2.5 min with SSD (if one seek per row)
#[derive(Debug, Clone)]
struct LonLatCosLat {
  lon: f64,
  lat: f64,
  cos_lat: f64,
}

/// MultiCone is primary made for a reasonable number (a few thousand) of possibly large
/// overlapping cones.
/// The typical use case is retrieving the sources around XMM Newton observation centers.
///
#[derive(Debug)]
pub struct MultiCone {
  /// Cone centres (lon, lat) in radians
  positions: Vec<(f64, f64)>,
  /// Radius, in radians
  r: f64,
  /// Allow for multi-threading
  n_threads: u16,
  // Derived quantity to speed up computation.
  /// Square Half Segment
  shs: f64,
  /// Depth of the hp hash used to build the hash map
  hash_map_cells_depth: u8,
  /// Mask used to compute the hash function returning an index in the hash map
  hash_map_mask: u32,
  /// Home made Hash Map implementation
  hash_map: Vec<Box<[LonLatCosLat]>>,
}

impl MultiCone {
  /// Create a new multicone from parameters expressed in degrees.
  /// # params
  /// * `positions_deg`: (longitude, latitudes) tuples, in degrees
  /// * `r_deg`: radius of the cone, in degrees
  /// * `n_threads`: with rayon feature, none returns the maximum of threads available, else returns 1.
  pub fn from_deg(
    positions_deg: Vec<(f64, f64)>,
    r_deg: f64,
    n_threads: Option<u16>,
  ) -> Result<Self, SkyRegionError> {
    let mut positions = positions_deg;
    // Convert position from deg to radians
    for (l, b) in &mut positions {
      *l = lon_deg2rad(*l)?;
      *b = lat_deg2rad(*b)?;
    }
    let r = r_deg.to_radians();
    if r <= 0.0 || PI <= r {
      Err(SkyRegionError::UnexpectedConeRadius { r_deg })
    } else {
      Self::new(positions, r, n_threads)
    }
  }

  /// # Panics
  /// * if one `lon` not in `[0, 2\pi[`
  /// * if one `lat` not in `[-\pi/2, \pi/2[`
  /// * if `r` not in `]0, \pi[`
  /// * `n_threads`: with rayon feature, none returns the maximum of threads available, else returns 1.
  pub fn new(
    positions: Vec<(f64, f64)>,
    r: f64,
    n_threads: Option<u16>,
  ) -> Result<Self, SkyRegionError> {
    if positions.len() > MAX_N_POSITIONS as usize {
      return Err(SkyRegionError::TooManyMulticoneCones {
        max: MAX_N_POSITIONS as usize,
        found: positions.len(),
      });
    }
    assert!(0.0 < r && r < PI);
    let hash_map_cells_depth = if has_best_starting_depth(r) {
      // * best_starting_depth => max 9 cells
      // * best_starting_depth - 1 => max 4 cells
      best_starting_depth(r).max(1) - 1
    } else {
      0
    };
    let n_pos = positions.len() as u32;
    // look at the previous power of 2 (to be used as hash map size)
    let hash_map_size = n_pos.next_power_of_two();
    let hash_map_mask = hash_map_size - 1;
    // The idea here is to duplicate the LonLatCosLat elements to avoid a memory indirection.
    // So our HashMap directly store arrays of LonLatCosLat (so for a same hash, object
    // are consecutive in memory).
    // Performances to be compared with a classical HashMap
    #[cfg(feature = "rayon")]
    let n_threads = n_threads.unwrap_or_else(|| num_cpus::get() as u16);
    #[cfg(not(feature = "rayon"))]
    let n_threads = n_threads.unwrap_or(1);
    let mut hash_map_tmp: Vec<(u32, LonLatCosLat)> = if n_threads <= 1 {
      compute_hash_map_tmp_one_thread(&positions, r, hash_map_cells_depth, hash_map_mask)
    } else {
      #[cfg(not(feature = "rayon"))]
      let hash_map_tmp =
        compute_hash_map_tmp_one_thread(&positions, r, hash_map_cells_depth, hash_map_mask);
      #[cfg(feature = "rayon")]
      let hash_map_tmp = compute_hash_map_tmp_multi_thread(
        &positions,
        r,
        hash_map_cells_depth,
        hash_map_mask,
        n_threads,
      );
      hash_map_tmp
    };
    // Build the final hash map (a simple array of vec)
    let mut hash_map: Vec<Box<[LonLatCosLat]>> = Vec::with_capacity(hash_map_size as usize);
    let mut prev_index = 0_u32;
    let mut entries: Vec<LonLatCosLat> = Default::default();
    for (index, entry) in hash_map_tmp.drain(..) {
      if index != prev_index {
        hash_map.push(std::mem::take(&mut entries).into_boxed_slice());
        for _ in prev_index + 1..index {
          hash_map.push(Vec::with_capacity(0).into_boxed_slice())
        }
        prev_index = index;
        assert_eq!(hash_map.len() as u32, index);
      }
      entries.push(entry);
    }
    hash_map.push(entries.into_boxed_slice());
    for _ in hash_map.len()..hash_map_size as usize {
      hash_map.push(Vec::with_capacity(0).into_boxed_slice())
    }
    assert_eq!(hash_map.len() as u32, hash_map_size);
    Ok(Self {
      positions,
      r,
      n_threads,
      shs: squared_half_segment_from_dist(r),
      hash_map_cells_depth,
      hash_map_mask,
      hash_map,
    })
  }

  pub fn is_multithreaded(&self) -> bool {
    self.n_threads > 1
  }
}

#[cfg(feature = "rayon")]
pub fn get_thread_pool(n_threads: u16) -> ThreadPool {
  ThreadPoolBuilder::new()
    .num_threads((n_threads as usize).min(num_cpus::get()))
    .build()
    .expect("Error initializing the thread pool")
}

fn compute_hash_map_tmp_one_thread(
  positions: &[(f64, f64)],
  r: f64,
  hash_map_cells_depth: u8,
  hash_map_mask: u32,
) -> Vec<(u32, LonLatCosLat)> {
  let mut hash_map_tmp: Vec<(u32, LonLatCosLat)> = positions
    .iter()
    .flat_map(|(lon, lat)| {
      assert!(0.0 <= *lon && *lon < TWICE_PI);
      assert!(-HALF_PI <= *lat && *lat <= HALF_PI);
      let cos_lat = lat.cos();
      let lon_lat_coslat = LonLatCosLat {
        lon: *lon,
        lat: *lat,
        cos_lat,
      };
      // we do not store the original index, the fast rejection (based on latitude diff)
      // will quickly reject most cases
      cone_coverage_approx_custom(hash_map_cells_depth, 2, *lon, *lat, r)
        .into_flat_iter()
        .map(move |idx| (idx as u32 & hash_map_mask, lon_lat_coslat.clone()))
    })
    .collect();
  hash_map_tmp.sort_by(|(idx1, _), (idx2, _)| idx1.cmp(idx2));
  hash_map_tmp
}

#[cfg(feature = "rayon")]
fn compute_hash_map_tmp_multi_thread(
  positions: &Vec<(f64, f64)>,
  r: f64,
  hash_map_cells_depth: u8,
  hash_map_mask: u32,
  n_threads: u16,
) -> Vec<(u32, LonLatCosLat)> {
  get_thread_pool(n_threads).install(|| {
    let mut hash_map_tmp = positions
      .par_iter()
      .flat_map_iter(|(lon, lat)| {
        assert!(0.0 <= *lon && *lon < TWICE_PI);
        assert!(-HALF_PI <= *lat && *lat <= HALF_PI);
        let cos_lat = lat.cos();
        let lon_lat_coslat = LonLatCosLat {
          lon: *lon,
          lat: *lat,
          cos_lat,
        };
        // we do not store the original index, the fast rejection (based on latitude diff)
        // will quickly reject most cases
        cone_coverage_approx_custom(hash_map_cells_depth, 2, *lon, *lat, r)
          .into_flat_iter()
          .map(move |idx| (idx as u32 & hash_map_mask, lon_lat_coslat.clone()))
      })
      .collect::<Vec<(u32, LonLatCosLat)>>();
    hash_map_tmp.par_sort_by(|(idx1, _), (idx2, _)| idx1.cmp(idx2));
    hash_map_tmp
  })
}

impl SkyRegion for MultiCone {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    let hash_map_index = hash(self.hash_map_cells_depth, lon, lat) as u32 & self.hash_map_mask;
    for llc in self.hash_map[hash_map_index as usize].as_ref() {
      if (llc.lat - lat).abs() < self.r
        && squared_half_segment(lon - llc.lon, lat - llc.lat, llc.cos_lat, lat.cos()) < self.shs
      {
        return true;
      }
    }
    false
  }

  fn area(&self) -> f64 {
    // Use MOCs ?!
    // Compute the precise n_cells of one cone by moc
    // Compute the n_cells of n MOCS
    // => compute the correction factor
    // let moc = RangeMOC::from_small_cones(depth, 2, stdin.lock().lines().filter_map(line2cone), None)
    // RangeMOC::from_large_cones(depth, 2, stdin.lock().lines().filter_map(line2cone))
    let depth = if has_best_starting_depth(self.r) {
      (best_starting_depth(self.r) + 4).min(DEPTH_MAX)
    } else {
      5
    };

    // Compute a correction factor
    const DELTA_DEPTH: u8 = 2;
    #[cfg(feature = "rayon")]
    let tot_n_cells: u64 = if self.is_multithreaded() {
      get_thread_pool(self.n_threads).install(|| {
        self
          .positions
          .par_iter()
          .map(|(lon, lat)| {
            RangeMOC::<u64, Hpx<u64>>::from_cone(
              *lon,
              *lat,
              self.r,
              depth,
              DELTA_DEPTH,
              CellSelection::All,
            )
            .n_depth_max_cells()
          })
          .sum()
      })
    } else {
      self
        .positions
        .iter()
        .map(|(lon, lat)| {
          RangeMOC::<u64, Hpx<u64>>::from_cone(
            *lon,
            *lat,
            self.r,
            depth,
            DELTA_DEPTH,
            CellSelection::All,
          )
          .n_depth_max_cells()
        })
        .sum()
    };
    #[cfg(not(feature = "rayon"))]
    let tot_n_cells: u64 = self
      .positions
      .iter()
      .map(|(lon, lat)| {
        RangeMOC::<u64, Hpx<u64>>::from_cone(
          *lon,
          *lat,
          self.r,
          depth,
          DELTA_DEPTH,
          CellSelection::All,
        )
        .n_depth_max_cells()
      })
      .sum();

    let union_n_cells = if depth <= (self.hash_map_cells_depth + 2) {
      // It means that the MOC of a single cone is made of a few cells
      // (so few chances to have a cell fully inside the cone).
      #[cfg(feature = "rayon")]
      let res = if self.is_multithreaded() {
        get_thread_pool(self.n_threads).install(|| {
          let mut cells: Vec<u64> = self
            .positions
            .par_iter()
            .flat_map_iter(|(lon, lat)| {
              cone_coverage_approx_custom(depth, DELTA_DEPTH, *lon, *lat, self.r).into_flat_iter()
            })
            .collect();
          cells.par_sort();
          RangeMOC::<u64, Hpx<u64>>::from_fixed_depth_cells(depth, cells.into_iter(), None)
        })
      } else {
        RangeMOC::<u64, Hpx<u64>>::from_small_cones(
          depth,
          2,
          self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
          None,
        )
      };
      #[cfg(not(feature = "rayon"))]
      let res = RangeMOC::<u64, Hpx<u64>>::from_small_cones(
        depth,
        2,
        self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
        None,
      );
      res
    } else {
      // Too many cones to compute the full union with BMOC code;
      // put aside the boolean telling that a cell is fully inside a cone
      #[cfg(feature = "rayon")]
      let res = if self.is_multithreaded() {
        get_thread_pool(self.n_threads).install(|| {
          self
            .positions
            .par_iter()
            .map(|(lon, lat)| {
              RangeMOC::<u64, Hpx<u64>>::from_cone(
                *lon,
                *lat,
                self.r,
                depth,
                DELTA_DEPTH,
                CellSelection::All,
              )
            })
            .reduce(
              || RangeMOC::<u64, Hpx<u64>>::new_empty(depth),
              |a, b| a.or(&b),
            )
        })
      } else {
        RangeMOC::<u64, Hpx<u64>>::from_large_cones(
          depth,
          2,
          CellSelection::All,
          self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
        )
      };
      #[cfg(not(feature = "rayon"))]
      let res = RangeMOC::<u64, Hpx<u64>>::from_large_cones(
        depth,
        2,
        CellSelection::All,
        self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
      );
      res
    }
    .n_depth_max_cells();
    let correction_factor = union_n_cells as f64 / tot_n_cells as f64;
    // End correction factor

    let one_cone_area = if self.r > 1.0e-4 {
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
    };
    correction_factor * self.positions.len() as f64 * one_cone_area
  }

  fn characteristic_depth(&self) -> u8 {
    if has_best_starting_depth(self.r) {
      (best_starting_depth(self.r) + 2).min(DEPTH_MAX)
    } else {
      2
    }
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    const DELTA_DEPTH: u8 = 2;
    const MAX_POS_BMOC: usize = 1000;
    if depth <= (self.hash_map_cells_depth + 2) {
      // It means that the MOC of a single cone is made of a few cells
      // (so few chances to have a cell fully inside the cone).
      #[cfg(feature = "rayon")]
      let moc = if self.is_multithreaded() {
        get_thread_pool(self.n_threads).install(|| {
          let mut cells: Vec<u64> = self
            .positions
            .par_iter()
            .flat_map_iter(|(lon, lat)| {
              cone_coverage_approx_custom(depth, DELTA_DEPTH, *lon, *lat, self.r).into_flat_iter()
            })
            .collect();
          cells.par_sort();
          RangeMOC::<u64, Hpx<u64>>::from_fixed_depth_cells(depth, cells.into_iter(), None)
        })
      } else {
        RangeMOC::<u64, Hpx<u64>>::from_small_cones(
          depth,
          2,
          self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
          None,
        )
      };
      #[cfg(not(feature = "rayon"))]
      let moc = RangeMOC::<u64, Hpx<u64>>::from_small_cones(
        depth,
        2,
        self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
        None,
      );
      let shift = Hpx::<u64>::shift_from_depth_max(depth);
      moc
        .into_range_moc_iter()
        .map(|range| (range.start >> shift..range.end >> shift, false))
        .collect()
    } else if self.positions.len() > MAX_POS_BMOC {
      // Too many cones to compute the full union with BMOC code;
      // put aside the boolean telling that a cell is fully inside a cone
      #[cfg(feature = "rayon")]
      let moc = if self.is_multithreaded() {
        get_thread_pool(self.n_threads).install(|| {
          self
            .positions
            .par_iter()
            .map(|(lon, lat)| {
              RangeMOC::<u64, Hpx<u64>>::from_cone(
                *lon,
                *lat,
                self.r,
                depth,
                DELTA_DEPTH,
                CellSelection::All,
              )
            })
            .reduce(
              || RangeMOC::<u64, Hpx<u64>>::new_empty(depth),
              |a, b| a.or(&b),
            )
        })
      } else {
        RangeMOC::<u64, Hpx<u64>>::from_large_cones(
          depth,
          2,
          CellSelection::All,
          self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
        )
      };
      #[cfg(not(feature = "rayon"))]
      let moc = RangeMOC::<u64, Hpx<u64>>::from_large_cones(
        depth,
        2,
        CellSelection::All,
        self.positions.iter().map(|(lon, lat)| (*lon, *lat, self.r)),
      );

      let shift = Hpx::<u64>::shift_from_depth_max(depth);
      moc
        .into_range_moc_iter()
        .map(|range| (range.start >> shift..range.end >> shift, false))
        .collect()
    } else {
      #[cfg(feature = "rayon")]
      let bmoc = if self.is_multithreaded() {
        get_thread_pool(self.n_threads).install(|| {
          self
            .positions
            .par_iter()
            .map(|(lon, lat)| cone_coverage_approx_custom(depth, DELTA_DEPTH, *lon, *lat, self.r))
            .reduce(
              || BMOCBuilderFixedDepth::new(depth, false).to_bmoc().unwrap(),
              |a, b| a.or(&b),
            )
        })
      } else {
        self
          .positions
          .iter()
          .map(|(lon, lat)| cone_coverage_approx_custom(depth, DELTA_DEPTH, *lon, *lat, self.r))
          .reduce(|a, b| a.or(&b))
          .unwrap_or_else(|| BMOCBuilderFixedDepth::new(depth, false).to_bmoc().unwrap())
      };
      #[cfg(not(feature = "rayon"))]
      let bmoc = self
        .positions
        .iter()
        .map(|(lon, lat)| cone_coverage_approx_custom(depth, DELTA_DEPTH, *lon, *lat, self.r))
        .reduce(|a, b| a.or(&b))
        .unwrap_or_else(|| BMOCBuilderFixedDepth::new(depth, false).to_bmoc().unwrap());

      bmoc.to_flagged_ranges()
    }
  }
}

#[cfg(test)]
mod tests {

  use super::*;

  #[test]
  fn test_multicone() {
    let pos_deg: Vec<(f64, f64)> = vec![
      (0.035323, 36.585958),
      (0.053309, 38.304050),
      (0.092231, -49.107945),
      (0.112292, -16.697009),
      (0.139568, 16.669049),
      (0.174446, 55.722460),
      (0.201418, 30.395890),
      (0.235772, 35.316720),
    ];
    let r_deg: f64 = 0.005;
    match MultiCone::new(
      pos_deg
        .into_iter()
        .map(|(ra, dec)| (ra.to_radians(), dec.to_radians()))
        .collect(),
      r_deg.to_radians(),
      Some(1),
    ) {
      Ok(multicone) => {
        println!("SortedHpxRanges: {:?}", multicone.sorted_hpx_ranges(6));
        assert!(true);
      }
      Err(e) => {
        println!("Error: {:?}", e);
        assert!(false);
      }
    }
  }
}
