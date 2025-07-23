// STC-S visitor!!

use std::ops::Range;

use nom::{
  error::{convert_error, VerboseError},
  Err,
};
use thiserror::Error;

use cdshealpix::nested::bmoc::BMOC;
use stc_s::{
  space::{
    common::{
      region::{BoxParams, CircleParams, ConvexParams, EllipseParams, PolygonParams},
      FillFrameRefposFlavor, Flavor, Frame, FromPosToVelocity, SpaceUnit,
    },
    position::Position,
    positioninterval::PositionInterval,
  },
  visitor::{impls::donothing::VoidVisitor, CompoundVisitor, SpaceVisitor, StcVisitResult},
  Stc,
};

use crate::{
  common::{
    error::SkyRegionError,
    math::{HALF_PI, PI, TWICE_PI},
  },
  regions::{cone::Cone, ellipse::EllipticalCone, polygon::Polygon, zone::Zone},
  SkyRegion,
};

#[derive(Debug)]
pub struct Stcs {
  stc: StcsElem,
}
impl Stcs {
  /// Create new StcsQuery from the given STC-S string.
  ///
  /// # WARNING
  /// * `DIFFERENCE` is interpreted as a symmetrical difference (it is a `MINUS` in the STC standard)
  /// * `Polygon` do not follow the STC-S standard: here self-intersecting polygons are supported
  /// * No implicit conversion: the STC-S will be rejected if
  ///     + the frame is different from `ICRS`
  ///     + the flavor is different from `Spher2`
  ///     + the units are different from `degrees`
  /// * Time, Spectral and Redshift sub-phrases are ignored
  ///
  /// # Params
  /// * `ascii_stcs`: lthe STC-S string
  ///
  /// # Output
  /// - The new StcsQuery (or an error)
  pub fn new(stcs: &str) -> Result<Self, SkyRegionError> {
    match Stc::parse::<VerboseError<&str>>(stcs.trim()) {
      Ok((rem, stcs)) => {
        if !rem.is_empty() {
          return Err(SkyRegionError::Stcs(Stc2StcQueryError::ParseHasRemaining {
            rem: rem.to_string(),
          }));
        }
        let StcVisitResult { space, .. } =
          stcs.accept(VoidVisitor, Stc2StcQuery, VoidVisitor, VoidVisitor);
        match space {
          None => Err(Stc2StcQueryError::NoSpaceFound),
          Some(space_res) => space_res.map(|stc| Self { stc }),
        }
      }
      Err(err) => Err(match err {
        Err::Incomplete(_) => Stc2StcQueryError::ParseIncomplete {
          msg: String::from("Incomplete parsing."),
        },
        Err::Error(e) => Stc2StcQueryError::ParseIncomplete {
          msg: convert_error(stcs, e),
        },
        Err::Failure(e) => Stc2StcQueryError::ParseIncomplete {
          msg: convert_error(stcs, e),
        },
      }),
    }
    .map_err(SkyRegionError::Stcs)
  }
}
impl SkyRegion for Stcs {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    self.stc.contains(lon, lat)
  }

  fn area(&self) -> f64 {
    self.stc.area()
  }

  fn characteristic_depth(&self) -> u8 {
    self.stc.characteristic_depth()
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    self.stc.to_bmoc(depth).to_flagged_ranges()
  }
}

#[derive(Debug)]
enum StcsElem {
  AllSky,
  Cone(Cone),
  Ellipse(EllipticalCone),
  Zone(Zone),
  Polygon(Polygon), // Box => Polygon
  Not {
    elem: Box<StcsElem>,
  },
  Union {
    elems: Vec<StcsElem>,
  },
  Intersection {
    elems: Vec<StcsElem>,
  },
  Difference {
    left: Box<StcsElem>,
    right: Box<StcsElem>,
  },
}
impl StcsElem {
  fn to_bmoc(&self, depth: u8) -> BMOC {
    match self {
      Self::AllSky => BMOC::new_allsky(depth),
      Self::Cone(elem) => elem.to_bmoc(depth),
      Self::Ellipse(elem) => elem.to_bmoc(depth),
      Self::Zone(elem) => elem.to_bmoc(depth),
      Self::Polygon(elem) => elem.to_bmoc(depth),
      Self::Not { elem } => elem.to_bmoc(depth).not(),
      Self::Union { elems } => elems
        .iter()
        .map(|elem| elem.to_bmoc(depth))
        .reduce(|acc, curr| acc.or(&curr))
        .unwrap_or(BMOC::new_empty(depth)),
      Self::Intersection { elems } => elems
        .iter()
        .map(|elem| elem.to_bmoc(depth))
        .reduce(|acc, curr| acc.and(&curr))
        .unwrap_or(BMOC::new_empty(depth)),
      Self::Difference { left, right } => left.to_bmoc(depth).xor(&right.to_bmoc(depth)),
    }
  }
}

impl SkyRegion for StcsElem {
  fn contains(&self, lon: f64, lat: f64) -> bool {
    match self {
      Self::AllSky => true,
      Self::Cone(elem) => elem.contains(lon, lat),
      Self::Ellipse(elem) => elem.contains(lon, lat),
      Self::Zone(elem) => elem.contains(lon, lat),
      Self::Polygon(elem) => elem.contains(lon, lat),
      Self::Not { elem } => !elem.contains(lon, lat),
      Self::Union { elems } => {
        // 'true' if at least one element contains the point
        for e in elems {
          if e.contains(lon, lat) {
            return true;
          }
        }
        false
      }
      Self::Intersection { elems } => {
        // 'true' if all elements contain the point
        for e in elems {
          if !e.contains(lon, lat) {
            return false;
          }
        }
        true
      }
      Self::Difference { left, right } => left.contains(lon, lat) ^ right.contains(lon, lat),
    }
  }

  // Approximate result!!!
  fn area(&self) -> f64 {
    match self {
      Self::AllSky => 4.0 * PI,
      Self::Cone(elem) => elem.area(),
      Self::Ellipse(elem) => elem.area(),
      Self::Zone(elem) => elem.area(),
      Self::Polygon(elem) => elem.area(),
      Self::Not { elem } => (4.0 * PI) - elem.area(),
      Self::Union { elems } => {
        // Consider disjoint areas
        elems.iter().map(|elem| elem.area()).sum::<f64>()
      }
      Self::Intersection { elems } => {
        // Consider overlapping areas
        elems
          .iter()
          .map(|elem| elem.area())
          .reduce(|min, curr| min.min(curr))
          .unwrap_or(0.0)
      }
      Self::Difference { left, right } => {
        // Consider disjoint areas and symmetric difference
        left.area() + right.area()
      }
    }
  }

  fn characteristic_depth(&self) -> u8 {
    match self {
      Self::AllSky => 0,
      Self::Cone(elem) => elem.characteristic_depth(),
      Self::Ellipse(elem) => elem.characteristic_depth(),
      Self::Zone(elem) => elem.characteristic_depth(),
      Self::Polygon(elem) => elem.characteristic_depth(),
      Self::Not { elem } => elem.characteristic_depth(),
      Self::Union { elems } => {
        // Consider disjoint areas
        elems
          .iter()
          .map(|elem| elem.characteristic_depth())
          .reduce(|max, curr| max.max(curr))
          .unwrap_or(0)
      }
      Self::Intersection { elems } => {
        // Consider overlapping areas
        elems
          .iter()
          .map(|elem| elem.characteristic_depth())
          .reduce(|max, curr| max.max(curr))
          .unwrap_or(0)
      }
      Self::Difference { left, right } => {
        // Consider disjoint areas and symmetric difference
        left
          .characteristic_depth()
          .max(right.characteristic_depth())
      }
    }
  }

  fn sorted_hpx_ranges(&self, depth: u8) -> Vec<(Range<u64>, bool)> {
    self.to_bmoc(depth).to_flagged_ranges()
  }
}

#[derive(Error, Debug)]
pub enum Stc2StcQueryError {
  #[error("Frame other than ICRS not supported (yet). Found: {found:?}")]
  FrameIsNotICRS { found: Frame },
  #[error("Flavor other than Spher2 not supported (yet). Found: {found:?}")]
  FlavorIsNotSpher2 { found: Flavor },
  #[error("Units ther than 'deg' not (yet?!) supported. Found: {found:?}")]
  UnitsNotSupported { found: Vec<SpaceUnit> },
  #[error("Convex shape not (yet?!) supported.")]
  ConvexNotSupported,
  #[error("Simple position not supported.")]
  SimplePositionNotSupported,
  #[error("Empty position range not supported.")]
  EmptyPositionRangeNotSupported,
  #[error("Position interval not supported.")]
  PositionIntervalNotSupported,
  #[error("invalid header (expected {expected}, found {found})")]
  WrongNumberOfParams { expected: u8, found: u8 },
  #[error("Longitude value out of bounds. Expected: [0, 360[. Actual: {value}")]
  WrongLongitude { value: f64 },
  #[error("Latitude value out of bounds. Expected: [-90, 90[. Actual: {value}")]
  WrongLatitude { value: f64 },
  #[error("STC-S string parsing not complete. Remaining: {rem}")]
  ParseHasRemaining { rem: String },
  #[error("STC-S string parsing incomplete: {msg}")]
  ParseIncomplete { msg: String },
  #[error("STC-S string parsing error: {msg}")]
  ParseFailure { msg: String },
  #[error("STC-S string parsing failure: {msg}")]
  ParseError { msg: String },
  #[error("No space sub-phrase found in STC-S string")]
  NoSpaceFound,
  #[error("Custom error: {msg}")]
  Custom { msg: String },
}

#[derive(Debug, Clone)]
struct Stc2StcQuery;

impl CompoundVisitor for Stc2StcQuery {
  type Value = StcsElem;
  type Error = Stc2StcQueryError;

  fn visit_allsky(&mut self) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::AllSky)
  }

  fn visit_circle(&mut self, circle: &CircleParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = circle
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty circle longitude"),
      })?;
    let lat_deg = circle
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty circle latitude"),
      })?;
    let radius_deg = circle.radius();
    // Convert params
    let lon = lon_deg2rad(*lon_deg)?;
    let lat = lat_deg2rad(*lat_deg)?;
    let r = radius_deg.to_radians();
    if r <= 0.0 || PI <= r {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Radius out of bounds in  CIRCLE. Expected: ]0, 180[. Actual: {}.",
          r
        ),
      })
    } else {
      Ok(StcsElem::Cone(Cone::new(lon, lat, r)))
    }
  }

  fn visit_ellipse(&mut self, ellipse: &EllipseParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = ellipse
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse longitude"),
      })?;
    let lat_deg = ellipse
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse latitude"),
      })?;
    let a_deg = ellipse.radius_a();
    let b_deg = ellipse.radius_b();
    let pa_deg = ellipse.pos_angle();
    // Convert params
    let lon = lon_deg2rad(*lon_deg)?;
    let lat = lat_deg2rad(*lat_deg)?;
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let pa = pa_deg.to_radians();
    if a <= 0.0 || HALF_PI <= a {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Semi-major axis out of bounds. Expected: ]0, 90[. Actual: {}.",
          a_deg
        ),
      })
    } else if b <= 0.0 || a <= b {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Semi-minor axis out of bounds. Expected: ]0, {}[. Actual: {}.",
          a_deg, b_deg
        ),
      })
    } else if pa <= 0.0 || PI <= pa {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Position angle out of bounds. Expected: [0, 180[. Actual: {}.",
          pa_deg
        ),
      })
    } else {
      Ok(StcsElem::Ellipse(EllipticalCone::new(lon, lat, a, b, pa)))
    }
  }

  fn visit_box(&mut self, skybox: &BoxParams) -> Result<Self::Value, Self::Error> {
    // Get params
    let lon_deg = skybox
      .center()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse longitude"),
      })?;
    let lat_deg = skybox
      .center()
      .get(1)
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty ellipse latitude"),
      })?;
    let mut a_deg = skybox
      .bsize()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty bsize on latitude"),
      })?;
    let mut b_deg = skybox
      .bsize()
      .first()
      .ok_or_else(|| Stc2StcQueryError::Custom {
        msg: String::from("Empty bsize on longitude"),
      })?;
    let mut pa_deg = skybox.bsize().first().copied().unwrap_or(90.0);
    if a_deg < b_deg {
      std::mem::swap(&mut b_deg, &mut a_deg);
      pa_deg = 90.0 - pa_deg;
    }
    // Convert params
    let lon = lon_deg2rad(*lon_deg)?;
    let lat = lat_deg2rad(*lat_deg)?;
    let a = a_deg.to_radians();
    let b = b_deg.to_radians();
    let pa = pa_deg.to_radians();
    if a <= 0.0 || HALF_PI <= a {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Box semi-major axis out of bounds. Expected: ]0, 90[. Actual: {}.",
          a_deg
        ),
      })
    } else if b <= 0.0 || a <= b {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Box semi-minor axis out of bounds. Expected: ]0, {}[. Actual: {}.",
          a_deg, b_deg
        ),
      })
    } else if !(0.0..PI).contains(&pa) {
      Err(Stc2StcQueryError::Custom {
        msg: format!(
          "Position angle out of bounds. Expected: [0, 180[. Actual: {}.",
          pa_deg
        ),
      })
    } else {
      Ok(StcsElem::Polygon(Polygon::from_box(lon, lat, a, b, pa)))
    }
  }

  fn visit_polygon(&mut self, polygon: &PolygonParams) -> Result<Self::Value, Self::Error> {
    let vertices_deg = polygon.vertices();
    let vertices = vertices_deg
      .iter()
      .step_by(2)
      .zip(vertices_deg.iter().skip(1).step_by(2))
      .map(|(lon_deg, lat_deg)| {
        let lon = lon_deg2rad(*lon_deg)?;
        let lat = lat_deg2rad(*lat_deg)?;
        Ok((lon, lat))
      })
      .collect::<Result<Vec<(f64, f64)>, Stc2StcQueryError>>()?;
    // Do something to use the right convention (depending on vertex order)!!
    Ok(StcsElem::Polygon(Polygon::new(vertices, false)))
  }

  fn visit_convex(&mut self, _convex: &ConvexParams) -> Result<Self::Value, Self::Error> {
    Err(Stc2StcQueryError::ConvexNotSupported)
  }

  fn visit_not(&mut self, elem: Self::Value) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Not {
      elem: Box::new(elem),
    })
  }

  fn visit_union(&mut self, elems: Vec<Self::Value>) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Union { elems })
  }

  fn visit_intersection(&mut self, elems: Vec<Self::Value>) -> Result<Self::Value, Self::Error> {
    Ok(StcsElem::Intersection { elems })
  }

  fn visit_difference(
    &mut self,
    left: Self::Value,
    right: Self::Value,
  ) -> Result<Self::Value, Self::Error> {
    // Warning: we interpret 'difference' as being a 'symmetrical difference', i.e. xor (not minus)
    Ok(StcsElem::Difference {
      left: Box::new(left),
      right: Box::new(right),
    })
  }
}

impl SpaceVisitor for Stc2StcQuery {
  type Value = StcsElem;
  type Error = Stc2StcQueryError;
  type C = Self;

  fn new_compound_visitor(
    &self,
    fill_frame_refpos_flavor: &FillFrameRefposFlavor,
    from_pos_to_velocity: &FromPosToVelocity,
  ) -> Result<Self, Self::Error> {
    // Check ICRS frame
    let frame = fill_frame_refpos_flavor.frame();
    if frame != Frame::ICRS {
      return Err(Stc2StcQueryError::FrameIsNotICRS { found: frame });
    }
    // Check SPHER2 flavor
    let flavor = fill_frame_refpos_flavor.flavor();
    if let Some(flavor) = flavor {
      if flavor != Flavor::Spher2 {
        return Err(Stc2StcQueryError::FlavorIsNotSpher2 { found: flavor });
      }
    }
    // Check units
    let opt_units = from_pos_to_velocity.unit().cloned();
    if let Some(units) = opt_units {
      for unit in units.iter().cloned() {
        if unit != SpaceUnit::Deg {
          return Err(Stc2StcQueryError::UnitsNotSupported { found: units });
        }
      }
    }
    Ok(self.clone())
  }

  fn visit_position_simple(self, _: &Position) -> Result<Self::Value, Self::Error> {
    Err(Stc2StcQueryError::SimplePositionNotSupported)
  }

  fn visit_position_interval(
    self,
    interval: &PositionInterval,
  ) -> Result<Self::Value, Self::Error> {
    // We use compound visitor only to check interval parameters
    self.new_compound_visitor(&interval.pre, &interval.post)?;
    let corners = interval
      .lo_hi_limits
      .iter()
      .step_by(2)
      .zip(interval.lo_hi_limits.iter().skip(1).step_by(2))
      .map(|(lon_deg, lat_deg)| {
        let lon = lon_deg2rad(*lon_deg)?;
        let lat = lat_deg2rad(*lat_deg)?;
        Ok((lon, lat))
      })
      .collect::<Result<Vec<(f64, f64)>, Stc2StcQueryError>>()?;
    let corners_it = corners
      .iter()
      .cloned()
      .step_by(2)
      .zip(corners.iter().cloned().skip(1).step_by(2));
    let mut zones = corners_it
      .map(|((ra_min, dec_min), (ra_max, dec_max))| {
        StcsElem::Zone(Zone::new(ra_min, dec_min, ra_max, dec_max))
      })
      .collect::<Vec<StcsElem>>();
    match zones.len() {
      0 => Err(Stc2StcQueryError::EmptyPositionRangeNotSupported),
      1 => Ok(zones.drain(..).next_back().unwrap()), // Unwrap ok since we tested the number of elements
      _ => Ok(StcsElem::Union { elems: zones }),
    }
  }

  fn visit_allsky(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_circle(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_ellipse(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_box(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_polygon(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_convex(self, _: StcsElem) -> Result<Self::Value, Self::Error> {
    unreachable!() // because an error is raised before calling this
  }

  fn visit_not(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_union(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_intersection(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }

  fn visit_difference(self, elem: StcsElem) -> Result<Self::Value, Self::Error> {
    Ok(elem)
  }
}

fn lon_deg2rad(lon_deg: f64) -> Result<f64, Stc2StcQueryError> {
  let mut lon = lon_deg.to_radians();
  if lon == TWICE_PI {
    lon = 0.0;
  }
  if !(0.0..TWICE_PI).contains(&lon) {
    Err(Stc2StcQueryError::WrongLongitude { value: lon_deg })
  } else {
    Ok(lon)
  }
}

fn lat_deg2rad(lat_deg: f64) -> Result<f64, Stc2StcQueryError> {
  let lat = lat_deg.to_radians();
  if !(-HALF_PI..=HALF_PI).contains(&lat) {
    Err(Stc2StcQueryError::WrongLatitude { value: lat_deg })
  } else {
    Ok(lat)
  }
}
