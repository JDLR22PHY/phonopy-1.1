"""Calculation of Phonon Angular Momentum (PAM)"""

# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the phonopy project nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.

import warnings
import numpy as np

# --- Existing Function ---

def phonon_angular_momentum(freq: np.ndarray, 
                           polar_vec: np.ndarray, 
                           temp: float = 0.0) -> np.ndarray:
    """
    Calculate phonon angular momentum (PAM) in units of ℏ.
    
    Parameters
    ----------
    freq : np.ndarray
        Phonon frequencies in THz.
    polar_vec : np.ndarray
        Phonon eigenvectors. Can be:
        - shape=(nqpts, nbnds, natoms, 3) for directly using components
        - shape=(nqpts, nbnds, 3*natoms) which will be reshaped 
    temp : float, optional
        Temperature in Kelvin (default: 0.0).
        
    Returns
    -------
    np.ndarray
        Phonon angular momentum in units of ℏ; shape=(3, nqpts, nbnds).
    """
    KB = 8.617330337217213e-05  # eV/K
    THZ_TO_EV = 0.00413566733   # Conversion factor from THz to eV

    # Handle temperature effects
    if np.isclose(temp, 0.0):
        nbose = 0.5
    else:
        nbose = 0.5 + 1.0 / (np.exp(freq * THZ_TO_EV / (KB * temp)) - 1.0)

    # Determine if we need to reshape the eigenvectors
    if len(polar_vec.shape) == 3:  # shape=(nqpts, nbnds, 3*natoms)
        # Determine number of atoms
        num_atoms = polar_vec.shape[2] // 3
        if polar_vec.shape[2] != 3 * num_atoms:
            raise ValueError(f"Eigenvectors third dimension size {polar_vec.shape[2]} is not divisible by 3.")
        # Reshape to separate atomic and Cartesian components
        E = polar_vec.reshape(polar_vec.shape[0], polar_vec.shape[1], num_atoms, 3)
    else:  # Already shape=(nqpts, nbnds, natoms, 3)
        E = polar_vec

    # Define the component pairs for angular momentum calculation (Jx, Jy, Jz)
    ixyz = [[1, 2], [2, 0], [0, 1]]  # For Jx, Jy, Jz respectively

    # Initialize Jxyz array
    Jxyz = np.zeros((3, freq.shape[0], freq.shape[1]), dtype=np.float64)

    for ii in range(3):
        e1 = E[:, :, :, ixyz[ii][0]]  # e_y or e_z for Jx, etc.
        e2 = E[:, :, :, ixyz[ii][1]]  # e_z or e_x for Jx, etc.
        # 2 * sum_over_atoms(Im(e1 * conj(e2)))
        Jxyz[ii] = 2.0 * np.sum((e1.conj() * e2).imag, axis=2)

    return Jxyz * nbose  # shape=(3, nqpts, nbnds)

# --- Helper Function to Read YAML Data ---

def read_band_yaml(filename: str):
    """
    Read phonon data from a YAML file.

    Parameters
    ----------
    filename : str
        Path to the YAML file (e.g. band.yaml).

    Returns
    -------
    tuple
        (Bcell, distances, frequencies, qpoints, segment_nqpoint, labels, eigenvectors)
        where 'eigenvectors' is an array with shape (nqpts, nbnds, 2, 3*num_atoms)
        and the last axis contains the real and imaginary parts.
    """
    import os
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader
    import yaml

    def open_file(file_path):
        _, ext = os.path.splitext(file_path)
        if ext in ('.xz', '.lzma'):
            import lzma
            return lzma.open(file_path)
        elif ext == '.gz':
            import gzip
            return gzip.open(file_path)
        else:
            return open(file_path, 'r')

    with open_file(filename) as f:
        data = yaml.load(f, Loader=Loader)

    freqs, dists, qpoints, labels, eigvec = [], [], [], [], []
    Bcell = np.array(data['reciprocal_lattice'])
    for v in data['phonon']:
        labels.append(v.get('label', None))
        freqs.append([f['frequency'] for f in v['band']])
        if 'eigenvector' in v['band'][0]:
            eigvec.append([np.array(f['eigenvector']) for f in v['band']])
        qpoints.append(v['q-position'])
        dists.append(v['distance'])
    # If all labels are None, try to use top-level "labels".
    if all(x is None for x in labels):
        if 'labels' in data:
            ss = np.array(data['labels'])
            labels = list(ss[0])
            for ii, f in enumerate(ss[:-1, 1] == ss[1:, 0]):
                if not f:
                    labels[-1] += r'|' + ss[ii+1, 0]
                labels.append(ss[ii+1, 1])
        else:
            labels = []
    segment_nqpoint = data.get('segment_nqpoint', [])
    return (Bcell, np.array(dists), np.array(freqs),
            np.array(qpoints), segment_nqpoint, labels, np.array(eigvec))

# --- New Function: Calculate PAM (with Optional YAML Input) ---

def calculate_pam(yaml_file: str = None,
                  freq: np.ndarray = None,
                  polar_vec: np.ndarray = None,
                  temp: float = 0.0):
    """
    Calculate phonon angular momentum (PAM) either by reading data from a YAML file
    or by using the provided frequency and eigenvector arrays.

    Parameters
    ----------
    yaml_file : str, optional
        Path to the YAML file containing phonon data (e.g., "band.yaml").
    freq : np.ndarray, optional
        Phonon frequencies in THz; shape=(nqpts, nbnds). Used if yaml_file is None.
    polar_vec : np.ndarray, optional
        Eigenvector data; shape=(nqpts, nbnds, 3*num_atoms). Used if yaml_file is None.
        (If the eigenvectors are stored as separate real and imaginary parts,
         they must be combined prior to calling this function.)
    temp : float, optional
        Temperature in Kelvin (default: 0.0).

    Returns
    -------
    If yaml_file is provided:
        tuple: (Jxyz, distances, frequencies, segment_nqpoint, labels)
    Otherwise:
        np.ndarray: Jxyz (angular momentum array with shape (3, nqpts, nbnds))
    """
    if yaml_file is not None:
        # Read phonon data from YAML
        Bcell, distances, frequencies, qpoints, segment_nqpoint, labels, eigvec = read_band_yaml(yaml_file)
        if eigvec.size == 0:
            raise ValueError("PHONON EIGENVECTORS MUST NOT BE EMPTY!")
        # Combine real and imaginary parts: assume eigvec shape is (nqpts, nbnds, 2, 3*num_atoms)
        E_complex = eigvec[..., 0] + 1j * eigvec[..., 1]
        Jxyz = phonon_angular_momentum(frequencies, E_complex, temp=temp)
        return Jxyz, distances, frequencies, segment_nqpoint, labels
    else:
        if freq is None or polar_vec is None:
            raise ValueError("Either yaml_file or both freq and polar_vec must be provided.")
        Jxyz = phonon_angular_momentum(freq, polar_vec, temp=temp)
        return Jxyz

# --- Existing Functions for Writing and Plotting (as provided earlier) ---

def write_pam_data(output_file: str,
                   distances: np.ndarray,
                   frequencies: np.ndarray,
                   Jxyz: np.ndarray,
                   segment_nqpoint) -> None:
    """
    Write phonon angular momentum data to a file.
    
    The file will contain data organized by band with columns:
    Distance, Frequency, Jx, Jy, Jz.
    
    Parameters
    ----------
    output_file : str
        Path to the output file.
    distances : np.ndarray
        1D array of distances along the phonon path.
    frequencies : np.ndarray
        2D array of phonon frequencies with shape (nqpts, nbnds).
    Jxyz : np.ndarray
        3D array of angular momentum components with shape (3, nqpts, nbnds).
    segment_nqpoint : list or array
        List of numbers of q-points in each segment.
    """
    with open(output_file, 'w') as f:
        f.write("# Data organized by band\n")
        f.write("# Each band is separated by a blank line\n")
        f.write("# Columns: Distance\tFrequency\tJx\tJy\tJz\n\n")
        nqpts, nbnds = frequencies.shape
        for band in range(nbnds):
            f.write(f"# Band {band+1}\n")
            data = np.column_stack((distances,
                                    frequencies[:, band],
                                    Jxyz[0, :, band],
                                    Jxyz[1, :, band],
                                    Jxyz[2, :, band]))
            np.savetxt(f, data, fmt='%.4f', delimiter='\t')
            f.write("\n")

def plot_pam_bands_from_data(
    Jxyz: np.ndarray,
    distances: np.ndarray,
    frequencies: np.ndarray,
    segment_nqpoint: list,
    labels: list = None,
    direction: str = 'x',
    plt_type: str = 'scatter',
    normalization: str = 'per_direction',
    ax = None,
    figsize: tuple = None
):
    """
    Plot phonon angular momentum (PAM) projection on band structure using pre-calculated data.
    
    Parameters
    ----------
    Jxyz : np.ndarray
        3D array of angular momentum components; shape=(3, nqpts, nbnds).
    distances : np.ndarray
        1D array of distances along the phonon path.
    frequencies : np.ndarray
        2D array of phonon frequencies; shape=(nqpts, nbnds).
    segment_nqpoint : list
        List of q-point counts per segment.
    labels : list, optional
        High-symmetry point labels for x-axis ticks.
    direction : str, optional
        Direction to plot: 'x', 'y', or 'z' (default: 'x').
    plt_type : str, optional
        Plot type: 'scatter' (default) or 'colormap'.
    normalization : str, optional
        Color normalization method: 'per_direction' (default) or 'all'.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new figure and axes are created.
    figsize : tuple, optional
        Figure size when creating a new figure.
        
    Returns
    -------
    tuple
        (fig, ax, color_mapper) - The figure, axes, and color mapper (for colorbar creation)
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.collections import LineCollection
    
    # Dictionary to map direction strings to indices
    dir_to_idx = {'x': 0, 'y': 1, 'z': 2}
    if direction not in dir_to_idx:
        raise ValueError(f"Direction must be one of 'x', 'y', 'z', got {direction}")
    
    comp = dir_to_idx[direction]
    
    # Create figure if needed
    if ax is None:
        if figsize is None:
            figsize = (8, 6)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Define colormap
    cmap_list = ["PiYG", "PuOr", "seismic"]
    
    # Set up normalization
    if normalization == 'per_direction':
        norm = mpl.colors.Normalize(vmin=Jxyz[comp, :, :].min(), vmax=Jxyz[comp, :, :].max())
    else:
        norm = mpl.colors.Normalize(vmin=Jxyz.min(), vmax=Jxyz.max())
    
    # Create scalar mappable for colorbar
    s_m = mpl.cm.ScalarMappable(cmap=cmap_list[comp], norm=norm)
    s_m.set_array([])
    
    # Plot each band
    for band in range(frequencies.shape[1]):
        # Plot in segments
        start = 0
        for nseg in segment_nqpoint:
            x = distances[start:start+nseg]
            y = frequencies[start:start+nseg, band]
            z = Jxyz[comp, start:start+nseg, band]
            
            # Plot gray line first
            line_width = 2.0 if plt_type == 'colormap' else 0.5
            ax.plot(x, y, lw=line_width, color='gray')
            
            # For colormap, add colored line on top
            if plt_type == 'colormap':
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments,
                                   colors=[s_m.to_rgba(ww)
                                          for ww in (z[1:] + z[:-1])/2.])
                lc.set_linewidth(2.0)
                lc.set_alpha(0.8)
                ax.add_collection(lc)
            
            start += nseg
    
    # Add scatter points if using scatter plot type
    if plt_type == 'scatter':
        for band in range(frequencies.shape[1]):
            ax.scatter(distances, frequencies[:, band],
                       s=np.abs(Jxyz[comp, :, band]) * 20,
                       c=Jxyz[comp, :, band],
                       cmap=cmap_list[comp],
                       norm=norm)
    
    # Add vertical lines at segment boundaries
    for bnd in np.cumsum(segment_nqpoint)[:-1]:
        ax.axvline(x=distances[int(bnd)], ls='--', color='gray', alpha=0.8, lw=0.5)
    
    # Add grid and set axis properties
    ax.axhline(y=0, ls='-', lw=0.5, color='k', alpha=0.6)
    ax.grid(True, ls='--', lw=0.5, color='gray', alpha=0.5)
    ax.set_xlim(distances.min(), distances.max())
    
    # Set x-ticks at high-symmetry points
    tick_indices = np.r_[[0], np.cumsum(segment_nqpoint) - 1]
    ax.set_xticks(distances[tick_indices])
    if labels:
        ax.set_xticklabels(labels)
    
    ax.set_ylabel('Frequency (THz)')
    ax.set_title(f'J{direction} Band Structure')
    
    return fig, ax, s_m

def plot_pam_data(distances: np.ndarray,
                 frequencies: np.ndarray,
                 Jxyz: np.ndarray,
                 segment_nqpoint,
                 labels=None,
                 direction: str = 'a',
                 plt_type: str = 'scatter',
                 normalization: str = 'per_direction',
                 layout: str = 'v',
                 figsize=None,
                 ax = None,
                 band_index: int = None):
    """
    Plot phonon angular momentum (PAM) along the phonon dispersion.

    Parameters
    ----------
    distances : np.ndarray
        1D array of distances along the phonon path.
    frequencies : np.ndarray
        2D array of phonon frequencies; shape=(nqpts, nbnds).
    Jxyz : np.ndarray
        3D array of angular momentum components; shape=(3, nqpts, nbnds).
    segment_nqpoint : list or array
        List of q-point counts per segment.
    labels : list, optional
        High-symmetry point labels for x-axis ticks.
    direction : str, optional
        'a' to plot all three components, or 'x', 'y', 'z' for a single component (default: 'a').
    plt_type : str, optional
        Plot type: 'scatter' (default) or 'colormap'.
    normalization : str, optional
        Color normalization method: 'per_direction' (default) or 'all'.
    layout : str, optional
        Layout of subplots: 'v' (vertical, default) or 'h' (horizontal).
    figsize : tuple, optional
        Figure size; if None, defaults will be chosen based on layout.
    figname : str, optional
        Filename to save the figure.
    band_index : int, optional
        If provided, only the specified band (index) is plotted.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the plot.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.collections import LineCollection

    nsub = 3 if direction == 'a' else 1
    if figsize is None:
        figsize = (6.4, 9.0) if layout == 'v' and nsub == 3 else (6.4, 4.8) if nsub == 1 else (12.0, 4.8)

    if layout == 'v':
        fig, axes = plt.subplots(nsub, 1, figsize=figsize, dpi=300)
    else:
        fig, axes = plt.subplots(1, nsub, figsize=figsize, dpi=300)
    if nsub == 1:
        axes = [axes]
    else:
        axes = list(axes)

    cmap_list = ["PiYG", "PuOr", "seismic"]
    if normalization == 'per_direction':
        norm_all = {comp: mpl.colors.Normalize(vmin=Jxyz[comp, :, :].min(), vmax=Jxyz[comp, :, :].max())
                    for comp in range(3)}
    else:
        overall_norm = mpl.colors.Normalize(vmin=Jxyz.min(), vmax=Jxyz.max())

    bands_to_plot = [band_index] if band_index is not None else range(frequencies.shape[1])
    for i, ax in enumerate(axes):
        comp = i if direction == 'a' else 'xyz'.index(direction)
        
        # Set up scalarmappable for colorbar
        norm_val = norm_all[comp] if normalization == 'per_direction' else overall_norm
        s_m = mpl.cm.ScalarMappable(cmap=cmap_list[comp], norm=norm_val)
        s_m.set_array([Jxyz[comp]])
        
        # Always plot gray lines for each band segment
        for band in bands_to_plot:
            start = 0
            for nseg in segment_nqpoint:
                x = distances[start:start+nseg]
                y = frequencies[start:start+nseg, band]
                z = Jxyz[comp, start:start+nseg, band]
                
                # Plot gray line first (thickness depends on plot type)
                line_width = 2.0 if plt_type == 'colormap' else 0.5
                ax.plot(x, y, lw=line_width, color='gray')
                
                # For colormap, add colored line on top
                if plt_type == 'colormap':
                    points = np.array([x, y]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
                    lc = LineCollection(segments,
                                       colors=[s_m.to_rgba(ww)
                                              for ww in (z[1:] + z[:-1])/2.])
                    lc.set_linewidth(2.0)
                    lc.set_alpha(0.8)
                    ax.add_collection(lc)
                
                start += nseg
                
        # Add scatter points if using scatter plot type
        if plt_type == 'scatter':
            for band in bands_to_plot:
                sc = ax.scatter(distances, frequencies[:, band],
                                s=np.abs(Jxyz[comp, :, band]) * 20,
                                c=Jxyz[comp, :, band],
                                cmap=cmap_list[comp],
                                norm=norm_val)
        
        # Add vertical lines at segment boundaries
        for bnd in np.cumsum(segment_nqpoint)[:-1]:
            ax.axvline(x=distances[int(bnd)], ls='--', color='gray', alpha=0.8, lw=0.5)
        
        ax.axhline(y=0, ls='-', lw=0.5, color='k', alpha=0.6)
        ax.grid(True, ls='--', lw=0.5, color='gray', alpha=0.5)
        ax.set_xlim(distances.min(), distances.max())
        
        segment_nqpoint_np = np.array(segment_nqpoint)
        # Create indices array: first point (0) and last point of each segment
        tick_indices = np.r_[[0], np.cumsum(segment_nqpoint_np) - 1]
        # Set the ticks at these distance positions
        ax.set_xticks(distances[tick_indices])
        if labels:
            ax.set_xticklabels(labels)
        
        ax.set_ylabel('Frequency (THz)')
        
        # Add colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right' if layout == 'v' else 'top',
                                  size='2%' if layout == 'v' else '4%',
                                  pad=0.04 if layout == 'v' else 0.03)
        
        cbar = plt.colorbar(s_m, cax=cax, extend='both', shrink=0.5,
                            orientation='vertical' if layout == 'v' else 'horizontal')
        cbar.ax.tick_params(which='both', labelsize='small')
        
        if layout == 'h':
            cbar.ax.xaxis.set_ticks_position('top')
            cbar.ax.xaxis.set_label_position('top')
            cbar.ax.set_xlabel(r'$J_{}\,/\,\hbar$'.format('xyz'[comp]))
        else:
            cbar.ax.text(1.60, 1.02, r'$J_{}\,/\,\hbar$'.format('xyz'[comp]),
                         ha="left", va="bottom", fontweight='bold',
                         transform=cbar.ax.transAxes,
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.5, lw=0.5))
    
    plt.tight_layout(pad=1.0)
    return fig
        