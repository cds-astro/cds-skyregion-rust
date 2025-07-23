use thiserror::Error;

use crate::regions::stcs::Stc2StcQueryError;

#[derive(Error, Debug)]
pub enum SkyRegionError {
  #[error("Unexpected SkyRegion longitude. Expected: [0, 360[. Actual: {lon_deg}.")]
  UnexpectedLongitude { lon_deg: f64 },
  #[error("Unexpected SkyRegion latitude. Expected: [-90, 90[. Actual: {lat_deg}.")]
  UnexpectedLatitude { lat_deg: f64 },
  #[error("Unexpected cone radius. Expected: ]0, 180[. Actual: {r_deg}.")]
  UnexpectedConeRadius { r_deg: f64 },
  #[error("Unexpected ring radii. Expected: (rmax in ]0, 180[, rmin in ]0, rmax[). Actual: (rmax = {r_max_deg}, rmin = {r_min_deg}).")]
  UnexpectedRingRadii { r_max_deg: f64, r_min_deg: f64 },
  #[error("Unexpected ellipse semi-major axis. Expected: ]0, 90]. Actual: {a_deg}.")]
  UnexpectedEllipseA { a_deg: f64 },
  #[error("Unexpected ellipse semi-minor axis. Expected: ]0, a[. Actual: {b_deg}.")]
  UnexpectedEllipseB { b_deg: f64 },
  #[error("Unexpected ellipse position angle. Expected: [0, 180[. Actual: {pa_deg}.")]
  UnexpectedEllipsePA { pa_deg: f64 },
  #[error("Unexpected box semi-major axis. Expected: ]0, 90]. Actual: {a_deg}.")]
  UnexpectedBoxA { a_deg: f64 },
  #[error("Unexpected box semi-minor axis. Expected: ]0, a[. Actual: {b_deg}.")]
  UnexpectedBoxB { b_deg: f64 },
  #[error("Unexpected box position angle. Expected: [0, 180[. Actual: {pa_deg}.")]
  UnexpectedBoxPA { pa_deg: f64 },
  #[error("Too many cones in multicone. Expected max: {max}. Actual: {found}.")]
  TooManyMulticoneCones { max: usize, found: usize },
  #[error("STC-S error: {0}.")]
  Stcs(Stc2StcQueryError),
  #[error("Custom error: {msg:?}")]
  Custom { msg: String },
}
