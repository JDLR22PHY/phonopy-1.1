# Phonopy with PAM

This version incorporates capabilities to compute the phonon angular momentum (PAM) and the corresponding density of states (DOS) resolved by the sign of the angular momentum. In this version, you can compute not only the total phonon DOS but also a DOS that is resolved into contributions from positive or negative PAM. In addition, the integration of both the total DOS and the PAM-resolved DOS is now available.

> **Important:**  
> In this release, the Bose–Einstein distribution is **not** included in the DOS integration, the temperature is fixed (i.e. no temperature-dependent occupation), and the number of frequency points (the resolution of the DOS) is currently hardcoded to 600.

## New Features

- **Phonon Angular Momentum Calculator:**  
  A new calculator that computes the phonon angular momentum for a given q-mesh. This feature allows one to analyze the angular momentum content of phonon modes.

- **PAM-Resolved DOS:**  
  In addition to the total DOS, you can now compute a density of states resolved by phonon angular momentum. That is, the integration can be performed separately for positive and negative phonon angular momentum contributions.

- **DOS Integration:**  
  The integrated DOS is available both for the total phonon DOS and for the PAM-resolved DOS. Note that at this stage, the integration does not account for the Bose–Einstein distribution, the temperature is fixed, and the frequency resolution is set to 600 points.

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
python -m pip install --upgrade pip
pip install numpy wheel setuptools
sudo apt-get install build-essential python3-dev
sudo apt-get update
rm -rf build/ dist/ *.egg-info/
pip install . -vvv
```

> **Disclaimer:**  
> A less recommended (but possible) approach for testing these new features is to copy the `phonopy` folder (the folder inside the repository that bears the same name) directly into the installation directory of the Phonopy library in your environment. **Warning:** This method is risky since it bypasses a full reinstallation/compilation of Phonopy and may lead to unexpected behavior if there are mismatches between the modified code and other installed modules.

## Usage

Once installed, you can run the new features using the `phonopy-load` command. For example, the following command initializes a mesh sampling calculation and then computes and plots the PAM-resolved DOS:

```bash
phonopy-load --mesh 61 61 61 -p -ldos --fmin 0 -s
```

Alternatively, running the `phonopy` command without the `-ldos` option will perform a standard calculation (e.g., band structure or mesh mode) without calculating the PAM-resolved DOS.

### Explanation of the Command Options

- `--mesh 61 61 61`  
  Specifies the q-mesh grid. In this example, a 61×61×61 mesh is used for sampling the phonon dispersion and DOS.

- `-p`  
  Tells Phonopy to plot the results. This flag enables the graphical output. (Note: The plotting mode is only active in modes that support it—here, it is used with mesh mode.)

- `-ldos`  
  Activates the new feature for calculating the **phonon angular momentum–resolved DOS (PAM-DOS)**. When this flag is set, Phonopy will run the new PAM DOS routines in place of the traditional DOS calculation.

- `--fmin 0`  
  Sets the minimum frequency for the DOS calculation. Only phonon frequencies above this threshold will be considered when integrating the DOS.

- `-s`  
  Instructs Phonopy to save the plotted graph (as a file, for example in PDF format) instead of just displaying it interactively.

## Example Workflow

1. **Standard PAM DOS Calculation:**  
   Run the command below to calculate the PAM-resolved DOS on a 61×61×61 mesh:
   ```bash
   phonopy-load --mesh 61 61 61 -p -ldos --fmin 0 -s
   ```
   In this case, Phonopy will:
   - Initialize a 61×61×61 q-mesh.
   - Run the new phonon angular momentum calculations.
   - Compute and integrate the PAM-resolved DOS (for both total and positive/negative PAM contributions).
   - Plot and save the PAM DOS (as indicated by `-p` and `-s`).

2. **Standard Calculation without PAM DOS:**  
   If you omit the `-ldos` flag, Phonopy will execute its standard mesh or band mode calculations:
   ```bash
   phonopy-load --mesh 61 61 61 -p --fmin 0 -s
   ```
   This command runs a mesh sampling calculation and plots the total DOS without resolving the angular momentum contributions.

## Notes

- The integration of the DOS is currently performed with 600 frequency points (fixed resolution).
- The temperature used in the DOS integration is fixed (no temperature dependence), and the Bose–Einstein distribution is not applied at this time.
- These new features are experimental and subject to further improvement in future releases.
