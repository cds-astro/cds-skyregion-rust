<meta charset="utf-8"/>

# `skyregion`

Defines possibly complex sky region (including STC-S regions) and test if a point on the unit-sphere is inside or not.
Also provide B-MOC like coverage to quickly retrieve (and possibly post-filter) indexed sources?  

## About

This Rust library was originally part of the QAT2S code but is now independent to be used elsewhere.

## Supported regions

* Cone
* Ring
* Elliptical Cone
* Multi-Cone
* Polygon
* Box (see one of the Polygon constructors)
* Zone
* JName (see one of the Zone constructors)
* STC-C (including complex STC-S with multiple union, intersection, ...)
* HEALPix:
  + single index
  + index range
  + set of (sorted and non-overlapping) ranges

## Features

You may activate the `rayon` feature to benefit from parallelism in the `multicone` region.
To do so, add in you `Cargo.toml` file:
```toml
skyregion = { version = 0.1.0, features = ["rayon"] }
```

## License

Like most projects in Rust, this project is licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or
  http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or
  http://opensource.org/licenses/MIT)

at your option.


## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this project by you, as defined in the Apache-2.0 license,
shall be dual licensed as above, without any additional terms or conditions.

### Warning

The code is formatted using 2 tab spaces instead of the regular 4:

```bash
cargo fmt -- --config tab_spaces=2
```

