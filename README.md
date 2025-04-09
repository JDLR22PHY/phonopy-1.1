# Phonopy with PAM

This version incorporates capabilities to compute the phonon angular momentum (PAM) and the corresponding density of states (DOS) resolved by the sign of the angular momentum. In this version, you can compute not only the total phonon DOS but also a DOS that is resolved into contributions from positive or negative PAM. In addition, the integration of both the total DOS and the PAM-resolved DOS is now available.

> **Important:**  
> The temperature-dependent calculations are now supported through the `--pam-temperature` option (default: 0). The frequency range for DOS integration can be controlled using `--fmin` and `--fmax` options. **The symmetry is turned off in this version.**

## New Features

- **Phonon Angular Momentum Calculator:**  
  A new calculator that computes the phonon angular momentum for a given q-mesh. This feature allows one to analyze the angular momentum content of phonon modes.

- **PAM-Resolved DOS:**  
  In addition to the total DOS, you can now compute a density of states resolved by phonon angular momentum. That is, the integration can be performed separately for positive and negative phonon angular momentum contributions.

- **DOS Integration:**  
  The integrated DOS is available both for the total phonon DOS and for the PAM-resolved DOS. The integration can be performed over a specified frequency range using `--fmin` and `--fmax` options.

- **Band Structure PAM:**  
  New capability to calculate PAM along band structure paths using the `--pam-bands` option. This requires a pre-existing `band.yaml` file.

- **Colored Phonon Dispersion Plotting:**  
  A new script `plot_phonon_dispersion_color.py` has been added to the Postprocessing_scripts directory. This script allows for visualization of phonon dispersion with PAM projections. See the script's README for detailed usage instructions.

## Installation

You can compile Phonopy from the source code as explained in the [official installation instructions](https://phonopy.github.io/phonopy/install.html#installation-from-source-code).

Here are the commands:
```bash
#These commands are for creating an environment
python -m venv phonopy_test
source phonopy_test/bin/activate
git clone https://github.com/brayanb1701/phonopy-1.git
cd phonopy-1
git checkout pam_ldos_integration
git tag -a v2.28.1 -m "version 2.28.1"
python -m pip install --upgrade pip
pip install numpy wheel setuptools
sudo apt-get install build-essential python3-dev
sudo apt-get update
rm -rf build/ dist/ *.egg-info/
pip install . -vvv || true
pip install . -vvv
```

> **Disclaimer:**  
> A less recommended (but possible) approach for testing these new features is to copy the `phonopy` folder (the folder inside the repository that bears the same name) directly into the installation directory of the Phonopy library in your environment. **Warning:** This method is risky since it bypasses a full reinstallation/compilation of Phonopy and may lead to unexpected behavior if there are mismatches between the modified code and other installed modules.

## Usage

Once installed, you can run the new features using the `phonopy-load` command. For example, the following command initializes a mesh sampling calculation and then computes and plots the PAM-resolved DOS:

```bash
phonopy-load --mesh 31 31 31 -p --pamdos --fmin 0 -s
```

To calculate PAM along band structure paths:
```bash
phonopy-load --pam-bands -p -s
```

### Explanation of the Command Options

- `--mesh 61 61 61`  
  Specifies the q-mesh grid. In this example, a 61×61×61 mesh is used for sampling the phonon dispersion and DOS.

- `-p`  
  Tells Phonopy to plot the results. This flag enables the graphical output.

- `--pamdos`  
  Activates the feature for calculating the **phonon angular momentum–resolved DOS (PAM-DOS)**. When this flag is set, Phonopy will run the PAM DOS routines.

- `--pam-bands`  
  Calculates PAM along band structure paths. Requires a pre-existing `band.yaml` file.

- `--pam-temperature`  
  Sets the temperature for PAM calculations (default: 0).

- `--fmin` and `--fmax`  
  Sets the frequency range for DOS integration.

- `--int-pamdos`  
  Integrates the PAM-resolved DOS and calculates the difference between positive and negative states.

- `-s`  
  Instructs Phonopy to save the plotted graph. This is recommended for better graph proportions.

- `--pam-cmap`: Specify a custom colormap for PAM bands plotting (e.g., 'viridis', 'plasma', 'inferno')

- `--with-pam-bands`: When used with PAMDOS calculation, plots the PAMDOS alongside the PAM projected phonon dispersion. Note that the frequency range is automatically scaled to match the PAMDOS range.

## Example Workflow

1. **Standard PAM DOS Calculation:**  
   Run the command below to calculate the PAM-resolved DOS on a 61×61×61 mesh:
   ```bash
   phonopy-load --mesh 31 31 31 -p --pamdos --fmin 0 -s
   ```
   In this case, Phonopy will:
   - Initialize a 31×31×31 q-mesh.
   - Run the phonon angular momentum calculations.
   - Compute and integrate the PAM-resolved DOS.
   - Plot and save the PAM DOS.

2. **PAM DOS with Custom Colormap and Bands:**  
   To calculate PAM DOS with a custom colormap and display the bands:
   ```bash
   phonopy-load --mesh 31 31 31 -p --pamdos --pam-cmap viridis --with-pam-bands -s
   ```
   This will:
   - Calculate the PAM DOS as before
   - Use the 'viridis' colormap for visualization
   - Display the PAM projected phonon dispersion alongside the DOS
   - Save the plot with proper proportions

3. **Band Structure PAM Calculation:**  
   To calculate PAM along band structure paths:
   ```bash
   phonopy-load --pam-bands -p -s
   ```
   This will:
   - Read the existing band structure from `band.yaml`
   - Calculate PAM along the band structure paths
   - Plot and save the results

4. **PAM DOS Integration:**  
   To integrate the PAM-resolved DOS and calculate the difference between positive and negative states:
   ```bash
   phonopy-load --mesh 61 61 61 --pamdos --int-pamdos --fmin 0 --fmax 10
   ```

## Notes

- The integration of the DOS is currently performed with 600 frequency points (fixed resolution).
- Temperature-dependent calculations are supported through the `--pam-temperature` option.
- The frequency range for DOS integration can be controlled using `--fmin` and `--fmax`.
- For best visualization results, always use the `-s` flag to save the plots rather than displaying them interactively.
- When using `--pam-bands`, ensure you have a valid `band.yaml` file in your working directory.
