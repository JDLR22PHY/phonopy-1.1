# Changelog

## [Unreleased] - 2025-04-09

### Added
- New script `plot_phonon_dispersion_color.py` in Postprocessing_scripts for colored phonon dispersion plotting
- Command line argument `--pam-cmap` to specify custom colormap for PAM bands plotting
- Command line argument `--with-pam-bands` to plot PAMDOS alongside PAM projected phonon dispersion
- Fixed colorbar range for PAM bands plot to [-1, 1]

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