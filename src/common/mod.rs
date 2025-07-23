pub mod error;
pub mod math;

use self::{
  error::SkyRegionError,
  math::{HALF_PI, TWICE_PI},
};

pub fn lon_deg2rad(lon_deg: f64) -> Result<f64, SkyRegionError> {
  let mut lon = lon_deg.to_radians();
  if lon == TWICE_PI {
    lon = 0.0;
  }
  if !(0.0..TWICE_PI).contains(&lon) {
    Err(SkyRegionError::UnexpectedLongitude { lon_deg })
  } else {
    Ok(lon)
  }
}

pub fn lat_deg2rad(lat_deg: f64) -> Result<f64, SkyRegionError> {
  let lat = lat_deg.to_radians();
  if !(-HALF_PI..=HALF_PI).contains(&lat) {
    Err(SkyRegionError::UnexpectedLatitude { lat_deg })
  } else {
    Ok(lat)
  }
}
