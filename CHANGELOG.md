# Changelog

## [2.28.1] - 2024-03-25

### Added
- New `--pam-bands` command option to calculate PAM along band structure paths
- New `--int-pamdos` command option to integrate PAM-resolved DOS and calculate difference between positive and negative states
- Temperature control for PAM calculations via `--pam-temperature` (default: 0)
- Frequency range control for DOS integration using `--fmin` and `--fmax`

### Changed
- Renamed `-ldos` flag to `--pamdos` for better clarity and consistency
- Updated plotting functionality to show PAM along band structure and PAM-resolved DOS side by side
- Improved graph proportions when saving with `-s` flag

### Fixed
- Integration functionality for PAM-resolved DOS
- Temperature-dependent PAM calculations
- Frequency range filtering for DOS integration

## [2.28.0] - Initial Release
- Initial implementation of PAM and PAM-resolved DOS features 