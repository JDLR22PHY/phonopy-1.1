"""Calculation of density of states."""

# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import warnings

import numpy as np

from phonopy.phonon.mesh import Mesh
from phonopy.phonon.tetrahedron_mesh import TetrahedronMesh
from phonopy.structure.tetrahedron_method import TetrahedronMethod
from phonopy.phonon.pam import phonon_angular_momentum

class NormalDistribution:
    """Class to represent normal distribution."""

    def __init__(self, sigma):
        """Init method."""
        self._sigma = sigma

    def calc(self, x):
        """Return normal distribution."""
        return (
            1.0
            / np.sqrt(2 * np.pi)
            / self._sigma
            * np.exp(-(x**2) / 2.0 / self._sigma**2)
        )


class CauchyDistribution:
    """Class to represent Cauchy distribution."""

    def __init__(self, gamma):
        """Init method."""
        self._gamma = gamma

    def calc(self, x):
        """Return Cauchy distribution."""
        return self._gamma / np.pi / (x**2 + self._gamma**2)


class Dos:
    """Base class to calculate density of states."""

    def __init__(self, mesh_object: Mesh, sigma=None, use_tetrahedron_method=False):
        """Init method."""
        self._mesh_object = mesh_object
        self._frequencies = mesh_object.frequencies
        self._weights = mesh_object.weights
        self._tetrahedron_mesh = None
        if use_tetrahedron_method and sigma is None:
            self._tetrahedron_mesh = TetrahedronMesh(
                mesh_object.dynamical_matrix.primitive,
                self._frequencies,
                mesh_object.mesh_numbers,
                np.array(mesh_object.grid_address, dtype="int_"),
                np.array(mesh_object.grid_mapping_table, dtype="int_"),
                mesh_object.ir_grid_points,
            )
        self._frequency_points = None
        self._sigma = sigma
        self.set_draw_area()
        self.set_smearing_function("Normal")

    @property
    def frequency_points(self):
        """Return frequency points."""
        return self._frequency_points

    def set_smearing_function(self, function_name):
        """Set function form for smearing method.

        Parameters
        ----------
        function_name : str
            'Normal': smearing is done by normal distribution.
            'Cauchy': smearing is done by Cauchy distribution.

        """
        if function_name == "Cauchy":
            self._smearing_function = CauchyDistribution(self._sigma)
        else:
            self._smearing_function = NormalDistribution(self._sigma)

    def integrate_states_in_frequency_range(self, freq_min, freq_max, states, temperature):
        """Integrate the given states over a specified frequency range,
        multiplying by the Bose-Einstein distribution.

        Parameters
        ----------
        freq_min : float
            The lower bound of the frequency range.
        freq_max : float
            The upper bound of the frequency range.
        states : ndarray
            The phonon states to integrate, must be an array of the same length
            as self.frequency_points or compatible shape.
        temperature : float
            Temperature in Kelvin.

        Returns
        -------
        integral : float or ndarray
            The integral of the states over the specified frequency range.
            If states is multidimensional, the integration is performed over
            the frequency axis.

        """
        # Ensure that frequency points and states have compatible shapes
        if states.shape[-1] != len(self._frequency_points):
            raise ValueError("The last dimension of states must match the length of frequency_points.")

        # Create a boolean mask for the frequency range
        freq_mask = (self._frequency_points >= freq_min) & (self._frequency_points <= freq_max)

        # Select the frequencies and states within the specified range
        freq_subset = self._frequency_points[freq_mask]
        states_subset = states[..., freq_mask]

        # Constants for unit conversions and Boltzmann constant
        KB = 8.617330337217213e-05  # Boltzmann constant in eV/K
        THZ_TO_EV = 0.00413566733   # Conversion factor from THz to eV

        # Calculate the Bose-Einstein distribution for the frequencies
        if np.isclose(temperature, 0.0):
            # At zero temperature, the Bose-Einstein distribution is zero
            bose_einstein = np.zeros_like(freq_subset)
        else:
            exponent = freq_subset * THZ_TO_EV / (KB * temperature)
            # To prevent overflow in exp(), cap the exponent to a reasonable value
            exponent = np.clip(exponent, None, 700)
            bose_einstein = 1.0 / (np.exp(exponent) - 1.0)

        # Multiply states by the Bose-Einstein distribution
        #weighted_states = states_subset * bose_einstein
        weighted_states = states_subset * 1
        
        # Perform the integration using the trapezoidal rule over the frequency axis
        integral = np.trapz(weighted_states, freq_subset, axis=-1)

        return integral

    def set_sigma(self, sigma):
        """Set sigma."""
        self._sigma = sigma

    def set_draw_area(self, freq_min=None, freq_max=None, freq_pitch=None):
        """Set frequency points."""
        f_min = self._frequencies.min()
        f_max = self._frequencies.max()

        if self._sigma is None:
            self._sigma = (f_max - f_min) / 100.0

        if freq_min is None:
            f_min -= self._sigma * 10
        else:
            f_min = freq_min

        if freq_max is None:
            f_max += self._sigma * 10
        else:
            f_max = freq_max

        if freq_pitch is None:
            f_delta = (f_max - f_min) / 600.0
        else:
            f_delta = freq_pitch
        self._frequency_points = np.arange(f_min, f_max + f_delta * 0.1, f_delta)


class TotalDos(Dos):
    """Class to calculate total DOS."""

    def __init__(self, mesh_object: Mesh, sigma=None, use_tetrahedron_method=False, 
                 freq_min=None, freq_max=None, temperature=0.1):
        """Init method."""
        super().__init__(
            mesh_object,
            sigma=sigma,
            use_tetrahedron_method=use_tetrahedron_method,
        )
        self._dos = None
        self._freq_Debye = None
        self._Debye_fit_coef = None
        self._openmp_thm = True

        self.freq_min = freq_min
        self.freq_max = freq_max
        self.temperature = temperature
        self.integrated_dos = None  # To store the result of integration

    def run(self):
        """Calculate total DOS."""
        if self._tetrahedron_mesh is None:
            print("Smearing")
            self._dos = np.array(
                [self._get_density_of_states_at_freq(f) for f in self._frequency_points]
            )
        else:
            print("Tetra")
            if self._openmp_thm:
                self._run_tetrahedron_method_dos()
            else:
                self._dos = np.zeros_like(self._frequency_points)
                thm = self._tetrahedron_mesh
                thm.set(value="I", frequency_points=self._frequency_points)
                for i, iw in enumerate(thm):
                    self._dos += np.sum(iw * self._weights[i], axis=1)
        if self.freq_min is None:
            self.freq_min = 0.1
        if self.freq_max is None:
            self.freq_max = 15
        # Perform integration over frequency range with Bose-Einstein distribution
        if self.freq_min is not None and self.freq_max is not None:
            self.integrated_dos = self.integrate_states_in_frequency_range(
                self.freq_min, self.freq_max, self._dos, self.temperature
            )
            print(
                f"Integrated DOS from {self.freq_min} to {self.freq_max} THz at {self.temperature} K: {self.integrated_dos}"
            )

    @property
    def dos(self):
        """Return total DOS."""
        return self._dos

    def get_dos(self):
        """Return frequency points and total DOS.

        Returns
        -------
        tuple
            (frequency_points, total_dos)

        """
        warnings.warn(
            "TotalDos.get_dos() is deprecated. "
            "Use frequency_points and dos attributes instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._frequency_points, self._dos

    def get_Debye_frequency(self):
        """Return a kind of Debye frequency."""
        return self._freq_Debye

    def set_Debye_frequency(self, num_atoms, freq_max_fit=None):
        """Calculate a kind of Debye frequency."""
        try:
            from scipy.optimize import curve_fit
        except ImportError as exc:
            raise ModuleNotFoundError("You need to install python-scipy.") from exc

        def Debye_dos(freq, a):
            return a * freq**2

        freq_min = self._frequency_points.min()
        freq_max = self._frequency_points.max()

        if freq_max_fit is None:
            N_fit = int(len(self._frequency_points) / 4.0)  # Hard coded
        else:
            N_fit = int(
                freq_max_fit / (freq_max - freq_min) * len(self._frequency_points)
            )
        popt, pcov = curve_fit(
            Debye_dos, self._frequency_points[0:N_fit], self._dos[0:N_fit]
        )
        a2 = popt[0]
        self._freq_Debye = (3 * 3 * num_atoms / a2) ** (1.0 / 3)
        self._Debye_fit_coef = a2

    def plot(self, ax, xlabel=None, ylabel=None, draw_grid=True, flip_xy=False):
        """Plot total DOS."""
        if flip_xy:
            _xlabel = "Density of states"
            _ylabel = "Frequency"
        else:
            _xlabel = "Frequency"
            _ylabel = "Density of states"

        if xlabel is not None:
            _xlabel = xlabel
        if ylabel is not None:
            _ylabel = ylabel

        plot_total_dos(
            ax,
            self._frequency_points,
            self._dos,
            freq_Debye=self._freq_Debye,
            Debye_fit_coef=self._Debye_fit_coef,
            xlabel=_xlabel,
            ylabel=_ylabel,
            draw_grid=draw_grid,
            flip_xy=flip_xy,
        )

    def write(self, filename="total_dos.dat"):
        """Write total DOS to total_dos.dat."""
        if self._tetrahedron_mesh is None:
            comment = "Sigma = %f" % self._sigma
        else:
            comment = "Tetrahedron method"

        write_total_dos(
            self._frequency_points, self._dos, comment=comment, filename=filename
        )

    def _run_tetrahedron_method_dos(self):
        mesh_numbers = self._mesh_object.mesh_numbers
        cell = self._mesh_object.dynamical_matrix.primitive
        reciprocal_lattice = np.linalg.inv(cell.cell)
        tm = TetrahedronMethod(reciprocal_lattice, mesh=mesh_numbers)
        self._dos = run_tetrahedron_method_dos(
            mesh_numbers,
            self._frequency_points,
            self._frequencies,
            self._mesh_object.grid_address,
            self._mesh_object.grid_mapping_table,
            tm.tetrahedra,
        )

    def _get_density_of_states_at_freq(self, f):
        return np.sum(
            np.dot(self._weights, self._smearing_function.calc(self._frequencies - f))
        ) / np.sum(self._weights)


class ProjectedDos(Dos):
    """Class to calculate projected DOS.

    Attributes
    ----------
    projected_dos : ndarray
        Projected DOS.
        shape=(pdos, frequency_points), dtype="double"
        The first dimension depends on init argument of `xyz_projection`.
        With it True, the length is the number of atoms times 3, otherwise
        the number of atoms.

    """

    def __init__(
        self,
        mesh_object: Mesh,
        sigma=None,
        use_tetrahedron_method=False,
        direction=None,
        xyz_projection=False,
    ):
        """Init method."""
        super().__init__(
            mesh_object,
            sigma=sigma,
            use_tetrahedron_method=use_tetrahedron_method,
        )
        self._eigenvectors = self._mesh_object.eigenvectors
        self._projected_dos = None

        if xyz_projection:
            self._eigvecs2 = np.abs(self._eigenvectors) ** 2
        else:
            num_atom = self._frequencies.shape[1] // 3
            i_x = np.arange(num_atom, dtype="int") * 3
            i_y = np.arange(num_atom, dtype="int") * 3 + 1
            i_z = np.arange(num_atom, dtype="int") * 3 + 2
            if direction is None:
                self._eigvecs2 = np.abs(self._eigenvectors[:, i_x, :]) ** 2
                self._eigvecs2 += np.abs(self._eigenvectors[:, i_y, :]) ** 2
                self._eigvecs2 += np.abs(self._eigenvectors[:, i_z, :]) ** 2
            else:
                d = np.array(direction, dtype="double")
                d /= np.linalg.norm(direction)
                proj_eigvecs = self._eigenvectors[:, i_x, :] * d[0]
                proj_eigvecs += self._eigenvectors[:, i_y, :] * d[1]
                proj_eigvecs += self._eigenvectors[:, i_z, :] * d[2]
                self._eigvecs2 = np.abs(proj_eigvecs) ** 2

        self._openmp_thm = True

    @property
    def partial_dos(self):
        """Return partial DOS."""
        warnings.warn(
            "PartialDos.partial_dos attribute is deprecated. "
            "Use projected_dos attribute instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._projected_dos

    @property
    def projected_dos(self):
        """Return projected DOS."""
        return self._projected_dos

    def run(self):
        """Calculate projected DOS."""
        if self._tetrahedron_mesh is None:
            self._run_smearing_method()
        else:
            if self._openmp_thm:
                self._run_tetrahedron_method_dos()
            else:
                self._run_tetrahedron_method()

    def get_partial_dos(self):
        """Return partial DOS.

        Returns
        -------
        tuple
            frequency_points: Sampling frequencies
            projected_dos: [atom_index, frequency_points_index]

        """
        warnings.warn(
            "ProjectedDos.get_partial_dos() is deprecated. "
            "Use frequency_points and projected_dos attributes instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._frequency_points, self._projected_dos

    def plot(
        self,
        ax,
        indices=None,
        legend=None,
        legend_prop=None,
        legend_frameon=True,
        xlabel=None,
        ylabel=None,
        draw_grid=True,
        flip_xy=False,
    ):
        """Plot projected DOS."""
        if flip_xy:
            _xlabel = "Partial density of states"
            _ylabel = "Frequency"
        else:
            _xlabel = "Frequency"
            _ylabel = "Partial density of states"

        if xlabel is not None:
            _xlabel = xlabel
        if ylabel is not None:
            _ylabel = ylabel

        plot_projected_dos(
            ax,
            self._frequency_points,
            self._projected_dos,
            indices=indices,
            legend=legend,
            legend_prop=legend_prop,
            legend_frameon=legend_frameon,
            xlabel=_xlabel,
            ylabel=_ylabel,
            draw_grid=draw_grid,
            flip_xy=flip_xy,
        )

    def write(self, filename="projected_dos.dat"):
        """Write projected DOS to projected_dos.dat."""
        if self._tetrahedron_mesh is None:
            comment = "Sigma = %f" % self._sigma
        else:
            comment = "Tetrahedron method"

        write_projected_dos(
            self._frequency_points,
            self._projected_dos,
            comment=comment,
            filename=filename,
        )

    def _run_smearing_method(self):
        num_pdos = self._eigvecs2.shape[1]
        num_freqs = len(self._frequency_points)
        self._projected_dos = np.zeros((num_pdos, num_freqs), dtype="double")
        weights = self._weights / float(np.sum(self._weights))
        for i, freq in enumerate(self._frequency_points):
            amplitudes = self._smearing_function.calc(self._frequencies - freq)
            for j in range(self._projected_dos.shape[0]):
                self._projected_dos[j, i] = np.dot(
                    weights, self._eigvecs2[:, j, :] * amplitudes
                ).sum()

    def _run_tetrahedron_method(self):
        num_pdos = self._eigvecs2.shape[1]
        num_freqs = len(self._frequency_points)
        self._projected_dos = np.zeros((num_pdos, num_freqs), dtype="double")
        thm = self._tetrahedron_mesh
        thm.set(value="I", frequency_points=self._frequency_points)
        for i, iw in enumerate(thm):
            w = self._weights[i]
            self._projected_dos += np.dot(iw * w, self._eigvecs2[i].T).T

    def _run_tetrahedron_method_dos(self):
        mesh_numbers = self._mesh_object.mesh_numbers
        cell = self._mesh_object.dynamical_matrix.primitive
        reciprocal_lattice = np.linalg.inv(cell.cell)
        tm = TetrahedronMethod(reciprocal_lattice, mesh=mesh_numbers)
        pdos = run_tetrahedron_method_dos(
            mesh_numbers,
            self._frequency_points,
            self._frequencies,
            self._mesh_object.grid_address,
            self._mesh_object.grid_mapping_table,
            tm.tetrahedra,
            coef=self._eigvecs2,
        )
        self._projected_dos = pdos.T


class PartialDos(ProjectedDos):
    """Class to calculate partial DOS."""

    def __init__(
        self,
        mesh_object: Mesh,
        sigma=None,
        use_tetrahedron_method=False,
        direction=None,
        xyz_projection=False,
    ):
        """Init method."""
        warnings.warn(
            "PartialDos class is deprecated. Use ProjectedDOS instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(
            mesh_object,
            sigma=sigma,
            use_tetrahedron_method=use_tetrahedron_method,
            direction=direction,
            xyz_projection=xyz_projection,
        )


def get_pdos_indices(symmetry):
    """Return atomic indieces grouped by symmetry."""
    mapping = symmetry.get_map_atoms()
    return [list(np.where(mapping == i)[0]) for i in symmetry.get_independent_atoms()]


def write_total_dos(
    frequency_points, total_dos, comment=None, filename="total_dos.dat"
):
    """Write total_dos.dat."""
    with open(filename, "w") as fp:
        if comment is not None:
            fp.write("# %s\n" % comment)

        for freq, dos in zip(frequency_points, total_dos):
            fp.write("%20.10f%20.10f\n" % (freq, dos))


def write_partial_dos(
    frequency_points, partial_dos, comment=None, filename="partial_dos.dat"
):
    """Write partial_dos.dat."""
    warnings.warn(
        "write_partial_dos() is deprecated. Use write_projected_dos() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    write_projected_dos(
        frequency_points, partial_dos, comment=comment, filename=filename
    )


def write_projected_dos(
    frequency_points, projected_dos, comment=None, filename="projected_dos.dat"
):
    """Write projected_dos.dat."""
    with open(filename, "w") as fp:
        if comment is not None:
            fp.write("# %s\n" % comment)

        for freq, pdos in zip(frequency_points, projected_dos.T):
            fp.write("%20.10f" % freq)
            fp.write(("%20.10f" * len(pdos)) % tuple(pdos))
            fp.write("\n")


def plot_total_dos(
    ax,
    frequency_points,
    total_dos,
    freq_Debye=None,
    Debye_fit_coef=None,
    xlabel=None,
    ylabel=None,
    draw_grid=True,
    flip_xy=False,
):
    """Plot total DOS."""
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_tick_params(which="both", direction="in")
    ax.yaxis.set_tick_params(which="both", direction="in")

    if freq_Debye is not None:
        freq_pitch = frequency_points[1] - frequency_points[0]
        num_points = int(freq_Debye / freq_pitch)
        freqs = np.linspace(0, freq_Debye, num_points + 1)

    if flip_xy:
        ax.plot(total_dos, frequency_points, "r-", linewidth=1)
        if freq_Debye:
            ax.plot(
                np.append(Debye_fit_coef * freqs**2, 0),
                np.append(freqs, freq_Debye),
                "b-",
                linewidth=1,
            )
    else:
        ax.plot(frequency_points, total_dos, "r-", linewidth=1)
        if freq_Debye:
            ax.plot(
                np.append(freqs, freq_Debye),
                np.append(Debye_fit_coef * freqs**2, 0),
                "b-",
                linewidth=1,
            )

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.grid(draw_grid)


def plot_partial_dos(
    ax,
    frequency_points,
    partial_dos,
    indices=None,
    legend=None,
    legend_prop=None,
    legend_frameon=True,
    xlabel=None,
    ylabel=None,
    draw_grid=True,
    flip_xy=False,
):
    """Plot partial DOS."""
    warnings.warn(
        "plot_partial_dos() is deprecated. Use plot_projected_dos() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    plot_projected_dos(
        ax,
        frequency_points,
        partial_dos,
        indices=indices,
        legend=legend,
        legend_prop=legend_prop,
        legend_frameon=legend_frameon,
        xlabel=xlabel,
        ylabel=ylabel,
        draw_grid=draw_grid,
        flip_xy=flip_xy,
    )


def plot_projected_dos(
    ax,
    frequency_points,
    projected_dos,
    indices=None,
    legend=None,
    legend_prop=None,
    legend_frameon=True,
    xlabel=None,
    ylabel=None,
    draw_grid=True,
    flip_xy=False,
):
    """Plot projected DOS."""
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_tick_params(which="both", direction="in")
    ax.yaxis.set_tick_params(which="both", direction="in")

    plots = []
    num_pdos = len(projected_dos)

    if indices is None:
        indices = []
        for i in range(num_pdos):
            indices.append([i])

    for set_for_sum in indices:
        pdos_sum = np.zeros_like(frequency_points)
        for i in set_for_sum:
            if i > num_pdos - 1:
                print("Index number '%d' is specified," % (i + 1))
                print("but it is not allowed to be larger than the number of " "atoms.")
                raise ValueError
            if i < 0:
                print(
                    "Index number '%d' is specified, but it must be "
                    "positive." % (i + 1)
                )
                raise ValueError
            pdos_sum += projected_dos[i]
        if flip_xy:
            plots.append(ax.plot(pdos_sum, frequency_points, linewidth=1))
        else:
            plots.append(ax.plot(frequency_points, pdos_sum, linewidth=1))

    if legend is not None:
        ax.legend(legend, prop=legend_prop, frameon=legend_frameon)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.grid(draw_grid)

# --------------------- New Additions for PAM Integration ---------------------

class PAMDos(Dos):
    """Class to calculate phonon DOS separated by positive and negative PAM for each axis."""

    def __init__(
        self, mesh_object: Mesh, sigma=None, use_tetrahedron_method=False, temperature=0.0, threshold=0.000000001, 
        freq_min=None,
        freq_max=None,
    ):
        """Initialize PAMDos."""
        super().__init__(
            mesh_object, sigma=sigma, use_tetrahedron_method=use_tetrahedron_method
        )
        self.temperature = temperature
        self.threshold = threshold
        self.Jxyz = None  # To store PAM
        self.dos_positive = None  # shape (3, frequency_points)
        self.dos_negative = None  # shape (3, frequency_points)
        self.axis_labels = ["x", "y", "z"]
        self._openmp_thm = True


        self.freq_min = freq_min
        self.freq_max = freq_max
        self.integrated_positive_dos = None  # To store integrated positive PAM DOS
        self.integrated_negative_dos = None  # To store integrated negative PAM DOS

    def run(self, int_pamdos=False):
        """Calculate PAM-resolved DOS."""
        # Calculate PAM
        self.Jxyz = phonon_angular_momentum(
            self._frequencies, self._mesh_object.eigenvectors, self.temperature
        )  # shape=(3, n_ir_grid_points, num_band)

        if self._tetrahedron_mesh is None:
            # Use smearing method
            self._run_smearing_method()
        else:
            # Use tetrahedron method
            self._run_tetrahedron_method_dos()


    def save_integration_results(self, filename="pam_dos_integration_results.txt"):
        """
        Computes the integrated positive and negative PAM DOS over the frequency range
        [freq_min, freq_max] at the given temperature, and saves the results to a text file.
        
        The output file contains:
        - Frequency range and temperature information.
        - The integrated positive PAM DOS.
        - The integrated negative PAM DOS.
        - The difference (positive minus negative).
        
        Parameters
        ----------
        filename : str
            The name of the output file. Defaults to "pam_dos_integration_results.txt".
        
        Raises
        ------
        ValueError
            If either freq_min or freq_max is not defined.
        """
        if self.freq_min is None:
            self.freq_min = 0

        if self.freq_max is None:
            self.freq_max = self._frequency_points.max()

        if int_pamdos:
            # Perform integration over frequency range
            if self.freq_min is not None and self.freq_max is not None:
                self.integrated_positive_dos = self.integrate_states_in_frequency_range(
                    self.freq_min, self.freq_max, self.dos_positive, self.temperature
                )
                self.integrated_negative_dos = self.integrate_states_in_frequency_range(
                    self.freq_min, self.freq_max, self.dos_negative, self.temperature
                )
                print(
                    f"Integrated Positive PAM DOS from {self.freq_min} to {self.freq_max} THz at {self.temperature} K:"
                )
                print(self.integrated_positive_dos)
                print(
                    f"Integrated Negative PAM DOS from {self.freq_min} to {self.freq_max} THz at {self.temperature} K:"
                )
                print(self.integrated_negative_dos)
                diff = self.integrated_positive_dos - self.integrated_negative_dos
                print(
                    f"Difference between Integrated Positive and Negative PAM DOS from {self.freq_min} to {self.freq_max} THz at {self.temperature} K:"
                )
                print(diff)

                # Save the integration results to a file
                filename = "int_pamdos_results.dat"
                with open(filename, "w") as f:
                    f.write(f"Integrated Positive PAM DOS: {self.integrated_positive_dos}\n")
                    f.write(f"Integrated Negative PAM DOS: {self.integrated_negative_dos}\n")
                    f.write(f"Difference: {diff}\n")
                    f.write(f"Frequency range: {self.freq_min} - {self.freq_max} THz\n")
                    f.write(f"Temperature: {self.temperature} K\n")
                print(f"Data saved successfully in {filename}")

    def _run_smearing_method(self):
        """Calculate PAM-resolved DOS using the smearing method."""
        # Initialize DOS arrays for positive and negative PAM per axis
        num_axes = 3
        self.dos_positive = np.zeros(
            (num_axes, len(self._frequency_points)), dtype="double"
        )
        self.dos_negative = np.zeros(
            (num_axes, len(self._frequency_points)), dtype="double"
        )

        # Expand weights to shape (nqpts, 1) for broadcasting
        expanded_weights = self._weights[:, np.newaxis]  # shape=(nqpts, 1)

        # Broadcast weights to shape (nqpts, nbnds)
        broadcast_weights = np.repeat(
            expanded_weights, self._frequencies.shape[1], axis=1
        )  # shape=(nqpts, nbnds)

        # For each axis, separate phonon modes by PAM sign and calculate DOS
        for axis in range(num_axes):
            # Positive PAM
            positive_mask = self.Jxyz[axis] > self.threshold  # shape=(nqpts, nbnds)
            positive_freqs = self._frequencies[positive_mask]  # shape=(num_positive,)
            positive_j = self.Jxyz[axis][positive_mask]
            # Multiply the original weight by |PAM|
            positive_weights = broadcast_weights[positive_mask] * np.abs(positive_j)
            self.dos_positive[axis] = self._calculate_dos(
                positive_freqs, positive_weights
            )

            # Negative PAM
            negative_mask = self.Jxyz[axis] < -self.threshold  # shape=(nqpts, nbnds)
            negative_freqs = self._frequencies[negative_mask]  # shape=(num_negative,)
            negative_j = self.Jxyz[axis][negative_mask]
            negative_weights = broadcast_weights[negative_mask] * np.abs(negative_j)
            self.dos_negative[axis] = self._calculate_dos(
                negative_freqs, negative_weights
            )

    def _run_tetrahedron_method_dos(self):
        """Calculate PAM-resolved DOS using the tetrahedron method."""
        mesh_numbers = self._mesh_object.mesh_numbers
        cell = self._mesh_object.dynamical_matrix.primitive
        reciprocal_lattice = np.linalg.inv(cell.cell)
        tm = TetrahedronMethod(reciprocal_lattice, mesh=mesh_numbers)

        # Initialize DOS arrays
        num_axes = 3
        self.dos_positive = np.zeros(
            (num_axes, len(self._frequency_points)), dtype="double"
        )
        self.dos_negative = np.zeros(
            (num_axes, len(self._frequency_points)), dtype="double"
        )

        num_grid_points = self._frequencies.shape[0]
        num_bands = self._frequencies.shape[1]

        for axis in range(num_axes):
            # Positive PAM
            positive_coef = np.where(self.Jxyz[axis] > self.threshold,
                           np.abs(self.Jxyz[axis]),
                           0)
            coef_positive = positive_coef.reshape(num_grid_points, 1, num_bands)
            # Negative PAM
            negative_coef = np.where(self.Jxyz[axis] < -self.threshold,
                           np.abs(self.Jxyz[axis]),
                           0)
            coef_negative = negative_coef.reshape(num_grid_points, 1, num_bands)
            # Calculate DOS for positive PAM modes
            dos_pos = run_tetrahedron_method_dos(
                mesh_numbers,
                self._frequency_points,
                self._frequencies,
                self._mesh_object.grid_address,
                self._mesh_object.grid_mapping_table,
                tm.tetrahedra,
                coef=coef_positive,
            )
            # dos_pos has shape (len(frequency_points), num_pdos)
            self.dos_positive[axis] = dos_pos[:, 0]

            # Calculate DOS for negative PAM modes
            dos_neg = run_tetrahedron_method_dos(
                mesh_numbers,
                self._frequency_points,
                self._frequencies,
                self._mesh_object.grid_address,
                self._mesh_object.grid_mapping_table,
                tm.tetrahedra,
                coef=coef_negative,
            )
            self.dos_negative[axis] = dos_neg[:, 0]

    def _calculate_dos(self, freqs_subset, weights_subset):
        """Calculate DOS for a subset of frequencies and weights using the smearing function."""
        # Vectorized calculation using broadcasting
        delta = self._frequency_points[:, np.newaxis] - freqs_subset[np.newaxis, :]
        dos = np.sum(weights_subset * self._smearing_function.calc(delta), axis=1)
        dos /= np.sum(self._weights)
        return dos

    def get_pam_dos(self):
        """Return frequency points and PAM-resolved DOS.

        Returns
        -------
        tuple
            (frequency_points, dos_positive, dos_negative)
            Each dos_positive and dos_negative has shape (3, frequency_points)
        """
        return self._frequency_points, self.dos_positive, self.dos_negative

    def plot_pam_dos(self, axes=None, figsize=(12, 9), xlabel="Frequency (THz)", ylabel="Density of States (DOS)", draw_grid=True):
        """
        Plot PAM-resolved DOS with three subplots, one for each PAM component.
        Positive and negative PAM DOS are plotted on opposite sides of the y-axis.
        
        Parameters
        ----------
        figsize : tuple, optional
            Size of the figure in inches (width, height).
        xlabel : str, optional
            Label for the x-axis.
        ylabel : str, optional
            Label for the y-axis.
        draw_grid : bool, optional
            Whether to draw grid lines on the plots.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Define labels and colors for each PAM component
        pam_labels = ['Jx', 'Jy', 'Jz']
        colors = ['r', 'g', 'b']  # Red for x, Green for y, Blue for z
        
        # Determine the maximum absolute DOS value for consistent y-axis limits
        max_dos = max(
            np.max(self.dos_positive),
            np.max(self.dos_negative)
        )
        if axes is None:
            # Create a figure with 3 vertically stacked subplots
            fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
        
        for i, ax in enumerate(axes):
            # Plot positive PAM DOS
            ax.plot(
                self._frequency_points,
                self.dos_positive[i],
                color=colors[i],
                linestyle='-',
                label=f'{pam_labels[i]} PAM > 0'
            )
            
            # Plot negative PAM DOS (mirrored below the x-axis)
            ax.plot(
                self._frequency_points,
                -self.dos_negative[i],
                color=colors[i],
                linestyle='--',
                label=f'{pam_labels[i]} PAM < 0'
            )
            
            # Set subplot title
            ax.set_title(f'Phonon DOS with {pam_labels[i]} Angular Momentum')
            
            # Set y-axis label
            ax.set_ylabel(ylabel)
            
            # Set y-axis limits for symmetry
            ax.set_ylim(-max_dos * 1.1, max_dos * 1.1)
            
            # Add legend
            ax.legend(loc='upper right')
            
            # Add grid if requested
            if draw_grid:
                ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
        
        # Set common x-axis label
        axes[-1].set_xlabel(xlabel)
    def write_pam_dos(self, filename="pam_dos.dat"):
        """
        Write PAM-resolved DOS to a file.

        The file has columns:
        Frequency, Jx+, Jx-, Jy+, Jy-, Jz+, Jz-

        Parameters
        ----------
        filename : str
            Output filename.
        """
        with open(filename, "w") as fp:
            fp.write("# Frequency\tJx+\tJx-\tJy+\tJy-\tJz+\tJz-\n")
            for i, freq in enumerate(self._frequency_points):
                line = f"{freq:.10f}"
                for axis in range(3):
                    line += f"\t{self.dos_positive[axis][i]:.10f}\t{self.dos_negative[axis][i]:.10f}"
                line += "\n"
                fp.write(line)


# --------------------- End of PAM Integration Additions ---------------------


def run_tetrahedron_method_dos(
    mesh,
    frequency_points,
    frequencies,
    grid_address,
    grid_mapping_table,
    relative_grid_address,
    coef=None,
):
    """Return (P)DOS calculated by tetrahedron method in C."""
    try:
        import phonopy._phonopy as phonoc
    except ImportError as exc:
        raise RuntimeError("Phonopy C-extension has to be built properly.") from exc

    if coef is None:
        _coef = np.ones((frequencies.shape[0], 1, frequencies.shape[1]), dtype="double")
    else:
        _coef = np.array(coef, dtype="double", order="C")
    arr_shape = frequencies.shape + (len(frequency_points), _coef.shape[1])
    dos = np.zeros(arr_shape, dtype="double")

    phonoc.tetrahedron_method_dos(
        dos,
        np.array(mesh, dtype="int_"),
        frequency_points,
        frequencies,
        _coef,
        np.array(grid_address, dtype="int_", order="C"),
        np.array(grid_mapping_table, dtype="int_", order="C"),
        relative_grid_address,
    )
    if coef is None:
        return dos[:, :, :, 0].sum(axis=0).sum(axis=0) / np.prod(mesh)
    else:
        return dos.sum(axis=0).sum(axis=0) / np.prod(mesh)


def get_dos_frequency_range(freqs, dos):
    """Return reasonable frequency range."""
    i_min = 0
    i_max = 1000

    for i, (_, d) in enumerate(zip(freqs, dos)):
        if d > 1e-5:
            i_min = i
            break

    for i, (_, d) in enumerate(zip(freqs[::-1], dos[::-1])):
        if d > 1e-5:
            i_max = len(freqs) - 1 - i
            break

    f_min = freqs[i_min]
    if f_min > 0:
        f_min = 0

    f_max = freqs[i_max]
    f_max += (f_max - f_min) * 0.05

    return f_min, f_max

