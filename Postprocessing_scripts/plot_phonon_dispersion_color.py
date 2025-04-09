#!/usr/bin/env python
import sys
import subprocess
import importlib.util
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

logging.info("Starting script execution")

def check_and_install_packages():
    """Check if required packages are installed and install them if needed."""
    logging.info("Checking required packages")
    required_packages = {
        'numpy': 'numpy',
        'matplotlib': 'matplotlib',
        'yaml': 'pyyaml'
    }
    
    missing_packages = []
    for module, package in required_packages.items():
        if importlib.util.find_spec(module) is None:
            missing_packages.append(package)
    
    if missing_packages:
        logging.info(f"Installing missing packages: {', '.join(missing_packages)}")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing_packages)
            logging.info("Packages installed successfully.")
        except subprocess.CalledProcessError:
            logging.error("Failed to install packages. Please install them manually.")
            print(f"pip install {' '.join(missing_packages)}")
            sys.exit(1)
    else:
        logging.info("All required packages are already installed")

# Check and install required packages
check_and_install_packages()

import numpy as np
import matplotlib.pyplot as plt
import yaml
import argparse
import re
from matplotlib.collections import LineCollection
from collections import defaultdict


def read_band_yaml(filename):
    """Read the band.yaml file and extract relevant data."""
    logging.info(f"Reading band.yaml file: {filename}")
    try:
        # Get file size to report progress
        file_size = os.path.getsize(filename)
        logging.info(f"YAML file size: {file_size / (1024*1024):.2f} MB")
        
        if file_size > 50 * 1024 * 1024:  # Larger than 50MB
            logging.warning("Large YAML file detected. Loading may take significant time and memory.")
            
        logging.info("Processing YAML file (this may take a while for large files)")
        
        # Try using LibYAML for faster parsing if available
        try:
            logging.info("Attempting to use faster C-based YAML loader")
            with open(filename, 'r') as f:
                data = yaml.load(f, Loader=yaml.CSafeLoader)
        except (AttributeError, ImportError):
            logging.info("Falling back to standard YAML loader")
            with open(filename, 'r') as f:
                data = yaml.safe_load(f)
                
        logging.info("Successfully loaded band.yaml file")
        return data
    except MemoryError:
        logging.error("Memory error when loading YAML file. Try using --max-qpoints to process a subset of data.")
        raise
    except Exception as e:
        logging.error(f"Error reading band.yaml file: {str(e)}")
        raise


def read_vesta_colors(filename):
    """Read the VESTA colors from the ini file."""
    logging.info(f"Reading VESTA colors file: {filename}")
    element_colors = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 7:  # Format: Z Element ... R G B
                    element = parts[1]
                    r, g, b = float(parts[5]), float(parts[6]), float(parts[7])
                    element_colors[element] = (r, g, b)
        
        logging.info(f"Loaded colors for {len(element_colors)} elements")
        return element_colors
    except Exception as e:
        logging.error(f"Error reading VESTA colors file: {str(e)}")
        raise


def calculate_weights(eigenvectors, points):
    """Calculate the weights per element from eigenvectors."""
    natom = len(eigenvectors)
    
    # Group atoms by element type
    element_indices = defaultdict(list)
    for i, point in enumerate(points):
        element = point['symbol']
        element_indices[element].append(i)
    
    # Calculate magnitudes for each atom
    magnitudes = np.zeros(natom)
    for atom_idx in range(natom):
        # Each atom has x, y, z components; each component has real and imaginary parts
        # Sum squared magnitudes for all directions
        total_sq_mag = 0
        for direction in range(3):  # x, y, z
            real_part = eigenvectors[atom_idx][direction][0]
            imag_part = eigenvectors[atom_idx][direction][1]
            total_sq_mag += real_part**2 + imag_part**2
        
        magnitudes[atom_idx] = total_sq_mag
    
    # Normalize to ensure sum of squares equals 1
    if np.sum(magnitudes) > 0:
        magnitudes = magnitudes / np.sum(magnitudes)
    
    # Sum up weights by element
    element_weights = {}
    for element, indices in element_indices.items():
        element_weights[element] = np.sum(magnitudes[indices])
    
    return element_weights


def mix_colors(element_weights, element_colors):
    """Mix colors based on element weights using RGB method."""
    r, g, b = 0, 0, 0
    
    for element, weight in element_weights.items():
        if element in element_colors:
            r += element_colors[element][0] * weight
            g += element_colors[element][1] * weight
            b += element_colors[element][2] * weight
    
    # Ensure values are in range [0, 1]
    r = max(0, min(1, r))
    g = max(0, min(1, g))
    b = max(0, min(1, b))
    
    return (r, g, b)


def calculate_distances(band_data):
    """Calculate the cumulative distances for the q-points."""
    logging.info("Calculating q-point distances")
    return [phonon['distance'] for phonon in band_data['phonon']]


def get_band_structure(band_data):
    """Extract band structure data from the YAML file."""
    logging.info("Extracting band structure data")
    q_distances = calculate_distances(band_data)
    
    # Get number of bands (modes) per q-point
    n_bands = len(band_data['phonon'][0]['band'])
    logging.info(f"Found {n_bands} bands at each q-point")
    
    # Get frequencies for each band across all q-points
    n_qpoints = len(band_data['phonon'])
    logging.info(f"Processing data for {n_qpoints} q-points")
    frequencies = np.zeros((n_qpoints, n_bands))
    
    # Array to store element weights for each band at each q-point
    element_weights_data = []
    
    logging.info("Calculating element weights for each band and q-point (this may take a while)")
    for i, phonon in enumerate(band_data['phonon']):
        if i % 100 == 0 and i > 0:
            logging.info(f"Processed {i}/{n_qpoints} q-points")
            
        for j, band in enumerate(phonon['band']):
            frequencies[i, j] = band['frequency']
            
            # Calculate weights for each element based on eigenvectors
            element_weights = calculate_weights(band['eigenvector'], band_data['points'])
            element_weights_data.append((i, j, element_weights))
    
    logging.info("Finished calculating element weights")
    return q_distances, frequencies, element_weights_data


def plot_phonon_bands(q_distances, frequencies, element_weights_data, element_colors, method='mixed', band_data=None, figure_size=(10, 7)):
    """Plot phonon bands with element-projected colors."""
    logging.info(f"Plotting phonon bands using {method} method")
    n_qpoints = len(q_distances)
    n_bands = frequencies.shape[1]
    
    plt.figure(figsize=figure_size)
    
    # Get unique elements for legend
    unique_elements = set()
    for _, _, element_weights in element_weights_data:
        unique_elements.update(element_weights.keys())
    
    logging.info(f"Found elements in the structure: {', '.join(sorted(unique_elements))}")
    
    # Plot band by band
    if method == 'mixed':
        # Mixed RGB approach
        logging.info("Using mixed RGB coloring method")
        for band_idx in range(n_bands):
            if band_idx % 5 == 0:
                logging.info(f"Plotting band {band_idx+1}/{n_bands}")
                
            points = np.array([q_distances, frequencies[:, band_idx]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            colors = []
            for i in range(len(segments)):
                # Get the element weights for this q-point and band
                for q_idx, b_idx, element_weights in element_weights_data:
                    if q_idx == i and b_idx == band_idx:
                        mixed_color = mix_colors(element_weights, element_colors)
                        colors.append(mixed_color)
                        break
            
            lc = LineCollection(segments, colors=colors, linewidths=2.5)
            plt.gca().add_collection(lc)
            
    elif method == 'scatter':
        # Scatter plot approach with color intensity based on element weight
        logging.info("Using scatter point coloring method")
        for element in unique_elements:
            # Skip if element color not found
            if element not in element_colors:
                logging.warning(f"No color found for element {element}, skipping")
                continue
                
            element_color = element_colors[element]
            logging.info(f"Plotting points for element {element}")
            
            points_plotted = 0
            for q_idx, band_idx, element_weights in element_weights_data:
                weight = element_weights.get(element, 0)
                if weight > 0.01:  # Only plot points with significant contribution
                    plt.scatter(q_distances[q_idx], frequencies[q_idx, band_idx], 
                                color=element_color, s=weight*20, alpha=weight, zorder=3)
                    points_plotted += 1
            
            logging.info(f"Plotted {points_plotted} points for element {element}")
    
    logging.info("Adding high-symmetry points and labels")
    # Plot high-symmetry points and labels
    special_points = []
    labels = []
    
    if band_data is not None:
        # Extract special points from segment_nqpoint
        idx = 0
        for i, nqpoint in enumerate(band_data['segment_nqpoint']):
            special_points.append(q_distances[idx])
            # Extract label from the segment info
            start_label = band_data['labels'][i][0]
            labels.append(start_label)
            
            idx += nqpoint - 1
            
        # Add the final point
        special_points.append(q_distances[-1])
        labels.append(band_data['labels'][-1][1])
    else:
        # If band_data not provided, just use min and max q_distances
        special_points = [q_distances[0], q_distances[-1]]
        labels = ['', '']
    
    # Plot vertical lines at high-symmetry points
    for x in special_points:
        plt.axvline(x=x, color='k', linestyle='-', linewidth=0.5)
    
    # Set the x-axis ticks and labels
    plt.xticks(special_points, labels)
    
    # Set plot limits and labels
    plt.xlim(q_distances[0], q_distances[-1])
    frequency_max = np.max(frequencies)
    plt.ylim(min(0, np.min(frequencies)), max(frequency_max, 18))
    plt.ylabel('Frequency (THz)')

    # Create legend for elements
    legend_elements = []
    for element in sorted(unique_elements):
        if element in element_colors:
            from matplotlib.patches import Patch
            legend_elements.append(Patch(facecolor=element_colors[element], label=element))
    
    plt.legend(handles=legend_elements, loc='upper right')
    

    
    # Create legend for elements (only for mixed method, scatter method handles its own legend)
    if method == 'mixed':
        legend_elements = []
        for element in sorted(unique_elements):
            if element in element_colors:
                from matplotlib.patches import Patch
                legend_elements.append(Patch(facecolor=element_colors[element], label=element))
        
        plt.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    logging.info("Finished creating plot")
    return plt.gcf()


if __name__ == "__main__":
    logging.info("Parsing command line arguments")
    parser = argparse.ArgumentParser(description='Plot phonon dispersion with element-projected colors')
    parser.add_argument('--band', default='band.yaml', help='Path to band.yaml file')
    parser.add_argument('--colors', default='elements-VESTA.ini.txt', help='Path to VESTA colors file')
    parser.add_argument('--method', choices=['mixed', 'scatter'], default='mixed', 
                        help='Coloring method: mix RGB values or use scatter points')
    parser.add_argument('--output', default='phonon_dispersion_colored.png', help='Output filename')
    parser.add_argument('--max-qpoints', type=int, default=None, 
                        help='Maximum number of q-points to process (for testing/debugging)')
    parser.add_argument('--downsample', type=int, default=1,
                        help='Downsample factor for q-points (e.g., 2 means use every 2nd point)')
    parser.add_argument('--no-display', action='store_true',
                        help='Do not display the plot window (useful for headless environments)')
    parser.add_argument('--dpi', type=int, default=300, 
                        help='DPI for the output figure (default: 300)')
    parser.add_argument('--figure-size', type=str, default='10,7',
                        help='Figure size in inches as width,height (default: 10,7)')
    args = parser.parse_args()
    
    logging.info(f"Arguments: band={args.band}, colors={args.colors}, method={args.method}, output={args.output}, "
                f"max_qpoints={args.max_qpoints}, downsample={args.downsample}, no_display={args.no_display}, "
                f"dpi={args.dpi}, figure_size={args.figure_size}")
    
    # Parse figure size
    try:
        figure_width, figure_height = map(float, args.figure_size.split(','))
    except:
        logging.warning(f"Invalid figure size format: {args.figure_size}. Using default 10,7.")
        figure_width, figure_height = 10, 7
    
    # Set the backend to non-interactive if requested
    if args.no_display:
        logging.info("Using non-interactive Agg backend")
        import matplotlib
        matplotlib.use('Agg')
    
    try:
        # Read the data
        band_data = read_band_yaml(args.band)
        
        # Process data based on max_qpoints and downsample options
        if args.max_qpoints is not None and args.max_qpoints > 0:
            logging.info(f"Limiting to {args.max_qpoints} q-points")
            # Keep only a subset of q-points
            band_data['phonon'] = band_data['phonon'][:args.max_qpoints]
            # When limiting q-points, we need to disable high-symmetry point labels
            logging.info("Disabling high-symmetry point labels due to q-point limiting")
            band_data_for_plot = None
        elif args.downsample > 1:
            logging.info(f"Downsampling q-points by factor of {args.downsample}")
            
            # Store original segment info for later use
            original_segment_nqpoint = band_data['segment_nqpoint'].copy()
            original_labels = band_data['labels'].copy()
            
            # Keep every nth q-point, but ensure we keep special points
            downsampled_phonons = []
            special_indices = []
            
            # Get all special point indices
            current_index = 0
            for nqpoint in original_segment_nqpoint:
                special_indices.append(current_index)  # Start of segment
                special_indices.append(current_index + nqpoint - 1)  # End of segment
                current_index += nqpoint
            
            # Downsample while preserving special points
            for i, phonon in enumerate(band_data['phonon']):
                if i % args.downsample == 0 or i in special_indices:
                    downsampled_phonons.append(phonon)
            
            band_data['phonon'] = downsampled_phonons
            
            # Calculate new segment_nqpoint values based on downsampling
            new_segment_nqpoint = []
            current_index = 0
            
            for i, nqpoint in enumerate(original_segment_nqpoint):
                # Calculate how many points will be in this segment after downsampling
                if i == 0:
                    # For first segment, count explicitly
                    points_before = 0
                    points_after = nqpoint - 1
                else:
                    points_before = sum(original_segment_nqpoint[:i])
                    points_after = points_before + nqpoint - 1
                
                # Count how many points we'll have after downsampling (including both endpoints)
                new_nqpoint = sum(1 for j in range(points_before, points_after + 1) 
                                if j % args.downsample == 0 or j in special_indices)
                
                # Ensure we have at least 2 points per segment (start and end)
                new_nqpoint = max(2, new_nqpoint)
                new_segment_nqpoint.append(new_nqpoint)
                current_index += new_nqpoint
            
            # Update segment info in band_data
            band_data['segment_nqpoint'] = new_segment_nqpoint
            
            # Since we preserved special points, we can use the band_data for plotting
            band_data_for_plot = band_data
            logging.info(f"Downsampled from {len(band_data['phonon'])} to {len(downsampled_phonons)} points while preserving special points")
        else:
            band_data_for_plot = band_data
        
        element_colors = read_vesta_colors(args.colors)
        
        # Get band structure data
        q_distances, frequencies, element_weights_data = get_band_structure(band_data)
        
        # Plot the phonon dispersion
        fig = plot_phonon_bands(q_distances, frequencies, element_weights_data, element_colors, 
                             method=args.method, band_data=band_data_for_plot,
                             figure_size=(figure_width, figure_height))
        
        # Save the figure with specified DPI
        logging.info(f"Saving figure to {args.output} with DPI={args.dpi}")
        fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        logging.info(f"Phonon dispersion plot saved as {args.output}")
        
        # Show the plot if display is enabled
        if not args.no_display:
            logging.info("Displaying plot")
            plt.show()
        
        logging.info("Script completed successfully")
    except Exception as e:
        logging.error(f"Error during script execution: {str(e)}", exc_info=True)
        sys.exit(1) 