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
