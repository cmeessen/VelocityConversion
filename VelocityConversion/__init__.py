# -*- coding: utf-8 -*-
################################################################################
#                     Copyright (C) 2017 by Christian Meeßen                   #
#                                                                              #
#                    This file is part of VelocityConversion                   #
#                                                                              #
#   VelocityConversion is free software: you can redistribute it and/or modify #
#     it under the terms of the GNU General Public License as published by     #
#            the Free Software Foundation version 3 of the License.            #
#                                                                              #
#   VelocityConversion is distributed in the hope that it will be useful, but  #
#          WITHOUT ANY WARRANTY; without even the implied warranty of          #
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU       #
#                   General Public License for more details.                   #
#                                                                              #
#      You should have received a copy of the GNU General Public License       #
#   along with VelocityConversion. If not, see <http://www.gnu.org/licenses/>. #
################################################################################
"""
This code is a python implementation of the p- and s-wave velocity to density
conversion approach after Goes et al. (2000)

To obtain help on the usage run
    Conversion.py --help

Recommended citation
Meeßen, Christian (2017): VelocityConversion. V. v1.0.1. GFZ Data Services.
http://doi.org/10.5880/GFZ.6.1.2017.001

For questions or suggestions please contact
Christian Meeßen
Section 6.1 Basin Modelling
GFZ Potsdam
14473 Potsdam, Germany
christian.meessen@gfz-potsdam.de
"""

import numpy as np
import sys
import math
import inspect
import os
import platform
from warnings import warn
import __main__


class MantleConversion:
    """
    The class that provides the code for the velocity conversion.
    """

    def __init__(self):

        # Define folder separator depending on the OS
        self.OS = platform.system()
        if self.OS == 'Windows':
            self.FolderSep = '\\'
        else:
            self.FolderSep = '/'
        self.PyFile = inspect.getfile(inspect.currentframe())
        self.PyPath = os.path.dirname(os.path.abspath(self.PyFile))

        # Other program related parameters
        self.FileIn = None
        self.FileOut = None
        self.Verbose = False
        self.SimpleP = False
        self.NN = False
        self.UseAlpha = 'Alpha'
        self.VelType = None
        self.scaleV = 1.

        # Physical constants
        self.g = 9.81          # gravitational acceleration
        self.TKelvin = 273.15  # offset temperature degC - Kelvin
        self.T0 = 300.         # room temperature
        self.P0 = 101300.      # atmospheric pressure in Pa
        self.R = 8.31446       # universal gas constant

        # Rock properties
        self.Composition = None
        self.XFe = None

        # Initialise mineral database
        self.MinDB = self.LoadMineralProperties()

        # AK135 Earth Model (Kennet et al., 1995)
        self.AK135Model = np.array([[0, 1.02],
                                    [3., 1.02],
                                    [3., 2.],
                                    [3.3, 2.],
                                    [3.3, 2.6],
                                    [10., 2.6],
                                    [10., 2.92],
                                    [18., 2.92],
                                    [18., 3.641],
                                    [43., 3.5801],
                                    [80., 3.502],
                                    [80., 3.502],
                                    [120., 3.4268],
                                    [120., 3.4268],
                                    [165., 3.3711],
                                    [210., 3.3243],
                                    [210., 3.3243],
                                    [260., 3.3663],
                                    [310., 3.4110]])
        self.MaxDepth = -1000*self.AK135Model[-1, 0]

        # Settings for temperature-density lookup tables
        # 3000 K top from Saxena and Shen (1992)
        self.Tmin = 300.
        self.Tmax = 3000.
        self.dT = 1.      # temperature increment
        self.nT = int((self.Tmax - self.Tmin)/self.dT + 1)

        # Attenuation properties
        self.Qmode = 'Sobolev et al. (1996)'
        self.A = 0.148
        self.f = 0.02
        self.omega = 2.*math.pi*self.f
        self.a = 0.15
        self.H = 500000.
        self.V = 0.000020

        # Check whether run as terminal program or loaded as module
        self.as_module = not hasattr(__main__, '__file__')

    def __raise__(self, error_type, msg):
        if self.as_module:
            if error_type is Warning:
                warn(msg)
            else:
                raise error_type(msg)
        else:
            if error_type is Warning:
                msg = "WARNING: " + msg
            else:
                msg = "ERROR: " + msg
            print(msg)
            if error_type is not Warning:
                sys.exit()

    def __check_vmin__(self):
        Vmin = np.amin(self.DataRaw[:, 3])
        if Vmin < 1000:
            msg = "Minimum velocity is " + str(Vmin) + "!\n"
            msg += "It appears that the velocities are provided in km/s but"
            msg += "they should be provided in m/s."
            if not self.as_module:
                msg += "\nUse -scaleV to modify the values on the fly."
            self.__raise__(Warning, msg)

    def __check_mineralogy__(self):
        """
        Check mineralogy information for consistency.
        """
        # Check if last line contains iron content information
        # Does not overwrite XFe if already defined by command line argument
        if self.Mineralogy['phase'][-1] == 'XFe' or \
           self.Mineralogy['phase'][-1] == 'xfe':
            if self.XFe is None:
                self.XFe = float(self.Mineralogy['fraction'][-1])
            else:
                msg = 'XFe already defined by command line argument'
                self.__raise__(Warning, msg)
            self.Mineralogy = np.delete(self.Mineralogy, -1)
        else:
            self.XFe = 0.0
            msg = 'XFe undefined!'
            self.__raise__(Warning, msg)
        print('> XFe assigned to', self.XFe)

        # Check if all mineral phases are known
        for i in self.Mineralogy['phase']:
            if i not in self.MinDB['phase']:
                msg = 'Unknown mineral phase: ' + str(i)
                self.__raise__(ValueError, msg)
        print('> All mineral phases are known')

        # Check if sum of mineral fractions equals 1
        frac_sum = np.round((np.sum(self.Mineralogy['fraction'])), decimals=10)
        if frac_sum != 1.0:
            msg = 'Sum of mineral phase fractions is ' + str(frac_sum)
            self.__raise__(ValueError, msg)
        else:
            print('> Sum of mineral phases equals 1.0')

        self.n_Phases = len(self.Mineralogy['phase'])

    def AK135(self, depth, simple=False):
        """
        Return Pressure at given depth based on AK135 model.
        Depth range: 18 to 310 km
        AK135: http://rses.anu.edu.au/seismology/ak135/ak135f.html

        Paramters:

        * depth : float
            Depth in m below sea level (positive!)

        * simple : boolean
            If `true`, return a lithostatic pressure based on an average
            density of 3000 kg/m3

        Returns:

        * Pressure : float
            Pressure in Pa at depth
        """
        if simple is False:
            DepthCalculated = False
            i = 0
            Pressure = 0.
            depth = np.absolute(depth)
            while DepthCalculated is False:
                if i > self.AK135Model.shape[0]:
                    print('Limit of AK135 model exceeded. End loop.')
                    break
                z1, Rho1 = self.AK135Model[i]*1000.
                z2, Rho2 = self.AK135Model[i+1]*1000.
                if depth < self.AK135Model[i+1, 0]*1000.:
                    Rho2 = Rho1 + (depth - z1)/(z2 - z1)*(Rho2 - Rho1)
                    dz = depth - z1
                    DepthCalculated = True
                else:
                    dz = z2 - z1
                Pressure = Pressure + dz*(Rho1 + (Rho2 - Rho1)/2)*self.g
                i += 1
            return Pressure
        else:
            self.SimpleP = True
            return 3000.*self.g*np.absolute(depth)

    def Alpha(self, P, T, Mineral, p=4.0):
        """
        Returns alpha for a specific mineral depending on the initially defined
        dependency. `self.UseAlpha` determines how the thermal conductivity is
        being calculated:

            * `self.UseAlpha = "Alpha"`
                Neither P nor T dependent, returns alpha_0 from `MinDB.csv`.

            * `self.UseAlpha = "AlphaT"`
                T-dependent conductivity based on Eqn. 5 in Saxena and Shen
                (1992)

            * `self.UseAlpha = "AlphaPT"`
                P- and T dependent conductivity which is caluclated by an
                inverse-distance weighed interpolation between the four closest
                grid points of `self.AlphaPT=AlphaDB.csv`.

        Parameters:

        * P : float
            Pressure in Pa

        * T : float
            Temperature in K

        * Mineral : int
            Mineral index in self.AlphaPT table

        * p : float
            Power parameter. Greater values of p assign greater influence to
            values closest to the interpolated point, with the result turning
            into a mosaic of tiles. For 2D interpolation, p <= 2 cause
            interpolated values to be dominated by points far away (wikipedia).

        Returns:

        * AlphaCalc : float
            Interpolated Alpha value for given parameters
        """

        P = P/1E9  # Transformation to GPa

        if self.UseAlpha == 'Alpha':
            # Alpha neither P nor T dependent
            for i in range(self.n_Phases):
                AlphaCalc = self.alpha[0, Mineral]

        elif self.UseAlpha == 'AlphaT':
            # Alpha only T dependent as used by Goes et al (2000)
            for i in range(self.n_Phases):
                AlphaCalc = self.alpha[0, Mineral] + self.alpha[1, Mineral]*T\
                    + self.alpha[2, Mineral]/T + self.alpha[3, Mineral]/T/T

        elif self.UseAlpha == 'AlphaPT':
            # Alpha P,T dependent, deduced from Hacker and Abers (2004) tables
            # 1 - Get row indices from P and T tables
            try:
                idxP = max(np.where(self.AlphaP <= P)[0])
            except ValueError:
                msg = "Pressure " + str(P) + " out of range!"
                self.__raise__(RuntimeError, msg)

            try:
                idxT = max(np.where(self.AlphaT <= T)[0])
            except ValueError:
                msg = "Temperature " + str(T) + " out of range!"
                self.__raise__(RuntimeError, msg)

            # 2 - Inverse-distance interpolation to get Alpha(P,T)
            Alpha_Tmin = min(self.AlphaT)
            Alpha_Tmax = max(self.AlphaT)
            Alpha_Pmin = min(self.AlphaP)
            Alpha_Pmax = max(self.AlphaP)

            # Analyse cases for T and define indices to lookup in T table
            if T == Alpha_Tmin or T == Alpha_Tmax:
                idxT0 = idxT1 = idxT
            else:
                idxT0 = idxT
                idxT1 = idxT + 1

            # Analyse cases for P and define indices to look up in P table
            if P == Alpha_Pmin or P == Alpha_Pmax:
                idxP0 = idxP1 = idxP
            else:
                idxP0 = idxP
                idxP1 = idxP + 1

            T0 = self.AlphaT[idxT0]
            T1 = self.AlphaT[idxT1]
            P0 = self.AlphaP[idxP0]
            P1 = self.AlphaP[idxP1]

            # Get the 4 surrounding alpha values
            Alphas = np.empty([2, 2])
            Alphas[0, 0] = self.AlphaPT[idxP0, idxT0, Mineral]
            Alphas[0, 1] = self.AlphaPT[idxP0, idxT1, Mineral]
            Alphas[1, 0] = self.AlphaPT[idxP1, idxT0, Mineral]
            Alphas[1, 1] = self.AlphaPT[idxP1, idxT1, Mineral]

            # Determine distances
            Dist = np.empty([2, 2])
            Dist[0, 0] = math.sqrt((P - P0)**2 + (T - T0)**2)
            Dist[0, 1] = math.sqrt((P - P0)**2 + (T - T1)**2)
            Dist[1, 0] = math.sqrt((P - P1)**2 + (T - T0)**2)
            Dist[1, 1] = math.sqrt((P - P1)**2 + (T - T1)**2)

            # Sum up
            Sum1 = Sum2 = 0
            for Pi in range(2):
                for Ti in range(2):
                    wi = Dist[Pi, Ti]**p
                    Sum1 += Alphas[Pi, Ti] / wi
                    Sum2 += 1/wi
            AlphaCalc = Sum1/Sum2

        return AlphaCalc

    def AssignPhases(self):
        """
        Create arrays for mineral properties Rho0, K, dKdT, etc. for phases
        defined in self.Mineralogy.
        """

        # Create empty arrays for the properties
        self.Rho0 = np.empty([self.n_Phases])
        self.dRhodX = self.Rho0.copy()
        self.K = self.Rho0.copy()
        self.dKdX = self.Rho0.copy()
        self.dKdP = self.Rho0.copy()
        self.dKdPdX = self.Rho0.copy()
        self.dKdT = self.Rho0.copy()
        self.mu = self.Rho0.copy()
        self.dmudX = self.Rho0.copy()
        self.dmudP = self.Rho0.copy()
        self.dmudT = self.Rho0.copy()
        self.alpha = np.empty([4, self.n_Phases])
        self.Composition = np.empty([self.n_Phases])

        # Transfer composition
        self.Composition = self.Mineralogy['fraction']

        # Clear all unused mineral phases in MinDB
        DelRows = []
        j = 0
        for i in self.MinDB['phase']:
            if i not in self.Mineralogy['phase']:
                DelRows.append(j)
            j += 1
        self.MinDB = np.delete(self.MinDB, DelRows)

        # Assign properties
        j = 0
        for i in self.Mineralogy['phase']:
            # Get index of phase i in MinDB
            IDX_phase = np.argwhere(self.MinDB['phase'] == i)[0][0]
            # Assign properties to the minerals
            self.Rho0[j] = self.MinDB['rho0'][IDX_phase]
            self.dRhodX[j] = self.MinDB['dRhodX'][IDX_phase]
            self.K[j] = self.MinDB['K0'][IDX_phase]
            self.dKdX[j] = self.MinDB['dKdX'][IDX_phase]
            self.dKdP[j] = self.MinDB['dKdP'][IDX_phase]
            self.dKdPdX[j] = self.MinDB['dKdPdX'][IDX_phase]
            self.dKdT[j] = self.MinDB['dKdT'][IDX_phase]
            self.mu[j] = self.MinDB['mu0'][IDX_phase]
            self.dmudX[j] = self.MinDB['dmudX'][IDX_phase]
            self.dmudP[j] = self.MinDB['dmudP'][IDX_phase]
            self.dmudT[j] = self.MinDB['dmudT'][IDX_phase]
            self.alpha[0, j] = self.MinDB['alpha0'][IDX_phase]
            self.alpha[1, j] = self.MinDB['alpha1'][IDX_phase]
            self.alpha[2, j] = self.MinDB['alpha2'][IDX_phase]
            self.alpha[3, j] = self.MinDB['alpha3'][IDX_phase]
            j += 1

    def CalcPT(self):
        """
        Estimate temperatures of input velocities. Synthetic data is stored in
        `self.SynPTVRho`. `SynPTVRho` dictionary structure::

            [
                { Depth_i   : np.array[T, V, Rho],
                  Depth_i+1 : np.array[T, V, Rho],
                  ...
                }
            ]

        The observed velocity is compared with synthetic velocities and then
        linearly interpolated between the neares calculated temperatures.

        Returns:

        * Nothing
            Updates Result_T and Result_Rho
        """
        print('Starting temperature estimation')
        self.Result_T = np.zeros([self.DataRaw.shape[0]])
        self.Result_Rho = np.zeros([self.DataRaw.shape[0]])
        j = self.ShowProgress()
        j_max = len(self.DataRaw)
        for val in self.DataRaw:
            Obs_z = val[2]
            Obs_v = val[3]
            i = 0
            i_max = self.SynPTVRho[Obs_z].shape[0] - 1
            # Search until observed temperature is greater than calculated
            while Obs_v < self.SynPTVRho[Obs_z][i, 1]:
                if (i + 1) <= i_max:
                    i += 1
                else:
                    i = i_max + 1
                    break
            # If found a match linearly interpolate temperature and density
            if i <= i_max:
                Syn_T1, Syn_v1, Syn_Rho1 = self.SynPTVRho[Obs_z][i]
                Syn_T2, Syn_v2, Syn_Rho2 = self.SynPTVRho[Obs_z][i-1]
                T_term1 = (Syn_T2 - Syn_T1)/(Syn_v2 - Syn_v1)
                self.Result_T[j] = T_term1 * (Obs_v - Syn_v1) + Syn_T1
                Rho_term1 = (Syn_Rho2 - Syn_Rho1)/(Syn_v2 - Syn_v1)
                self.Result_Rho[j] = Rho_term1*(Obs_v - Syn_v1) + Syn_Rho1
            # Return -1 for T or Rho if no match could be found
            else:
                self.Result_T[j] = self.T0 - 1
                self.Result_Rho[j] = -1
            j = self.ShowProgress(j, j_max)
        print('> Done!')

    def CalcVRho(self, z, T):
        """
        Calculate velocities based on depth and temperature

        Parameters:

        * z : double
            Depth in km

        * T : numpy array
            Array containing temperature values

        Returns:

        * TableVRho : numpy array
            Array contining V and Rho according for specified z and T
        """
        # Calculate P
        P = self.AK135(z)
        if self.Verbose:
            print('> z='+str(z)+'km; P='+str(P)+'Pa')

        # Calculate mu minerals
        if self.Verbose:
            print('> Calculate mu minerals')
        mu_minerals = np.zeros([self.n_Phases, self.nT])
        for i in range(self.n_Phases):
            mu_minerals[i, :] = self.mu[i] + (T - self.T0)*self.dmudT[i] +\
                (P - self.P0)*self.dmudP[i] + self.XFe*self.dmudX[i]

        # Calculate mu rock
        if self.Verbose:
            print('> Calculate mu rock')
        mu_voigt = np.zeros([self.nT])
        mu_reuss = np.zeros([self.nT])
        for i in range(self.n_Phases):
            mu_voigt = mu_voigt + self.Composition[i]*mu_minerals[i, :]
            mu_reuss = mu_reuss + self.Composition[i]/mu_minerals[i, :]
        mu_reuss = 1./mu_reuss
        mu_rock = (mu_voigt + mu_reuss)/2.

        # K minerals
        if self.Verbose:
            print('> Calculate K minerals')
        K_minerals = np.zeros([self.n_Phases, self.nT])
        for i in range(self.n_Phases):
            K_minerals[i, :] = self.K[i] + (T - self.T0)*self.dKdT[i] +\
                (P - self.P0)*(self.dKdP[i] + self.XFe*self.dKdPdX[i]) + \
                self.XFe*self.dKdX[i]

        # For P-waves K_rock is required
        if self.VelType == 'P':
            if self.Verbose:
                print('> Calculate K rock')
            K_voigt = np.zeros([self.nT])
            K_reuss = np.zeros([self.nT])
            for i in range(self.n_Phases):
                K_voigt = K_voigt + self.Composition[i]*K_minerals[i, :]
                K_reuss = K_reuss + self.Composition[i]/K_minerals[i, :]
            K_reuss = 1./K_reuss
            K_rock = (K_voigt + K_reuss)/2.

        # Calculate alpha_minerals
        if self.Verbose:
            print('> Calculate alpha minerals')
        alpha_minerals = np.zeros([self.n_Phases, self.nT])
        for i in range(self.n_Phases):
            Tidx = 0
            for Ti in T:
                alpha_minerals[i, Tidx] = self.Alpha(P, Ti, i)
                Tidx += 1

        # Calculate rho minerals
        if self.Verbose:
            print('> Calculate rho minerals')
        rho_minerals = np.empty([self.n_Phases, self.nT])
        rho_XFe = self.Rho0 + self.XFe*self.dRhodX
        for i in range(self.n_Phases):
            rho_minerals[i, :] = (rho_XFe[i]*(1 - alpha_minerals[i, :]*(T -
                                  self.T0) + (P - self.P0)/K_minerals[i, :]))

        # Calculate rho rock
        if self.Verbose:
            print('> Calculate rho rock')
        rho_rock = np.zeros([self.nT])
        for i in range(self.n_Phases):
            rho_rock = rho_rock + self.Composition[i]*rho_minerals[i, :]

        # Calculate attenuation
        if self.Verbose:
            print('> Calculate attenuation')
        Q = np.zeros([self.nT])
        E = self.H + P*self.V
        Q = self.A*(self.omega**self.a)*np.exp((self.a*E)/self.R/T)
        if self.VelType == 'P':
            L = 4.*mu_rock/(3.*K_rock + 4.*mu_rock)
            Q = Q/L
        elif self.VelType != 'S':
            raise ValueError('Wrong velocity type.')

        # Calculate v_syn
        if self.Verbose:
            print('> Calculate v_syn')
        if self.VelType == 'S':
            v_anh = (mu_rock/rho_rock)**0.5
        elif self.VelType == 'P':
            v_anh = ((K_rock + 4./3.*mu_rock)/rho_rock)**0.5
        else:
            raise ValueError('Wrong velocity type.')
        v_anel = 1. - 2./Q/math.tan(self.a*math.pi/2.)

        TableVRho = np.empty([self.nT, 2])
        TableVRho[:, 0] = v_anh*v_anel
        TableVRho[:, 1] = rho_rock

        return TableVRho

    def FillTables(self):
        """
        Create arrays of V and T for the depths in FileIn. Also loads tables
        for Alpha(P,T) if necessary. The arrays are stored in the dictionary
        `self.SynPTVRho`::

            [
                { Depth_i   : np.array[T, V, Rho],
                  Depth_i+1 : np.array[T, V, Rho],
                  ...}
            ]
        """

        if self.UseAlpha == 'AlphaPT':
            self.LoadAlphaTables()

        print('Filling tables')
        nDepths = self.Depths.shape[0]
        arr_T = np.linspace(self.Tmin, self.Tmax, self.nT)
        self.SynPTVRho = {}

        if self.SimpleP:
            CalcPType = 'Simple'
        else:
            CalcPType = 'AK135'

        print('> Number of depth values:', nDepths)
        print('> Number of temperatures:', self.nT)
        print('> Pressure calculation  :', CalcPType)
        print('> T range:', self.Tmin, 'to', self.Tmax, 'steps', self.dT)

        for i in range(nDepths):
            # arr_TVRho contains T and corresponding V
            arr_TVRho = np.zeros([self.nT, 3])
            arr_TVRho[:, 0] = arr_T
            arr_TVRho[:, 1:] = self.CalcVRho(self.Depths[i], arr_T)
            depth_i = {self.Depths[i]: arr_TVRho}
            self.SynPTVRho.update(depth_i)

    def ShowHelp(self):
        """
        Display function help
        """
        print()
        print("Inverts seismic velocities for temperature and density based")
        print("on input composition. By default the expansion coefficient is")
        print("treated as pressure- and temperature-independent.")
        print()
        print("Usage of ", sys.argv[0] + ":")
        print()
        print("Minimum requirements")
        print("--------------------")
        print()
        print(sys.argv[0], "FileIn -type <P|S> -comp <Filename>")
        print("    FileIn")
        print("        Input file name and path")
        print("    -type <P|S>")
        print("        Defines wave type P or S")
        print("    -comp <Filename>")
        print("        Comma-separated file containing mantle rock composition")
        print("    Example")
        print("       ", sys.argv[0], "Input.dat -type S -comp pyrolite.csv")
        print("        The output file will be InputOut.dat")
        print()
        print("Input file requirements")
        print("-----------------------")
        print()
        print("    Input file order and units must be: X Y Z V")
        print("    X Y - Coordinates in any units")
        print("    Z   - Depth in meters below or above sea level")
        print("    V   - Seismic velocity in m/s")
        print("    Can contain header lines if they are marked, i.e. with #")
        print()
        print("Optional arguments")
        print("------------------")
        print()
        print("    -AlphaT")
        print("        Calculate Alpha based on Temperature after Saxena and")
        print("        Shen (1992). Default: P/T-independent.")
        print()
        print("    -AlphaPT")
        print("        Calculate P/T-dependent Alpha based on excel worksheet")
        print("        from Hacker and Abers (2004). Default: P/T-independent.")
        print()
        print("    -dT")
        print("        Changes the temperature increment in the P-T tables.")
        print("        Default increment is 1K.")
        print()
        print("    -NN")
        print("        Output file header information reduced to \# of points")
        print()
        print("    -out <FileOut>")
        print("        Define the output path and/or file name")
        print("        Example: -out ../Output.dat")
        print()
        print("    -scaleV <value>")
        print("        Scale the velocity with the given value")
        print()
        print("    -setQ <1|2>")
        print("        Define anelasticity parameters after 1 - Sobolev et al.")
        print("        (1996) or 2 - Berckhemer et al. (1982). Default: 1.")
        print("        Example: -setQ 2")
        print()
        print("    -v | -verbose")
        print("        Displays debugging messages.")
        print()
        print("    -XFe <XFe> | -xfe <XFe>")
        print("        Define iron content XFe in mole fractions. Default")
        print("        value XFe = 0.0")
        print("        Example: -XFe 0.1")
        print()
        print("For questions or bug reports contact")
        print("christian.meessen@gfz-potsdam.de")
        sys.exit()

    def LoadAlphaTables(self, FileIn='AlphaDB.csv'):
        """
        Load tables for P- and T-dependent thermal conductivity, alpha(P,T),
        that were previously extracted from the worksheet by Hacker and Abers
        (2004, see `./AlphaPT`).

        Parameters:

        * FileIn : string
            File name of table containing Alpha values. Can be commented with
            '#' and must have a row containing names of the columns. Minimum
            requried columns: P, T, Mineral1, .. IMPORTANT: the data must be
            sorted 1st after T, 2nd after P (descending).

        Output:

        * self.AlphaPT : numpy array
            Three-dimensional array where 1st dimension is P, 2nd dimension T
            and 3rd dimension the mineral

        * self.AlphaP : numpy array
            Pressures values corresponding to indices in AlphaPT

        * self.AlphaT : numpy array
            Temperature values corresponding to indices in AlphaPT
        """
        print("Importing tables for Alpha(P,T)")
        FileIn = self.PyPath + self.FolderSep + FileIn
        AlphaTables = np.genfromtxt(FileIn, names=True, dtype=None,
                                    skip_header=2, delimiter=';')
        nP = len(np.unique(AlphaTables['P']))
        nT = len(np.unique(AlphaTables['T']))
        nMinerals = len(self.Mineralogy['phase'])
        # Transfer data to 2D arrays
        self.AlphaPT = np.empty([nP, nT, nMinerals])
        self.AlphaP = np.unique(AlphaTables['P'])
        self.AlphaT = np.unique(AlphaTables['T'])
        for i in range(nMinerals):
            i_phase = self.Mineralogy['phase'][i]
            self.AlphaPT[:, :, i] = np.reshape(AlphaTables[i_phase], (nP, nT))

    def LoadMineralogy(self, FileIn=None):
        """
        Import mineralogy of the mantle rock from an input file and perform
        checks on the input file.

        Parameters:

        * FileIn : string
            Input file name and path
        """

        # Load file
        if FileIn is None:
            msg = "No mineralogy defined. Using default mineralogy."
            self.__raise__(Warning, msg)
            self.DefaultMineralogy()
        else:
            print('Loading mineralogy', FileIn)
            RawData = np.genfromtxt(FileIn, names=True, dtype=None,
                                    delimiter=',', encoding='utf8')

        # Check if all dtypes known
        ValidDtypes = ['phase', 'fraction']
        for i in RawData.dtype.names:
            if i not in ValidDtypes:
                msg = 'Invalid data type \"' + str(i) + '\" in ' + FileIn
                self.__raise__(ValueError, msg)
        print('> All datatypes known')

        # Sort alphabetically
        self.Mineralogy = np.sort(RawData, order='phase')
        self.__check_mineralogy__()
        self.AssignPhases()

    def LoadMineralProperties(self, FileIn=None, delim=';'):
        """
        Import Mineral properties from database file. The input file must be a
        delimited text file with the following columnar structure::

            Full phase name - Self explanatory
            phase  - Mineral phase
            rho0   - Density at T0, P0 / kg/m3
            dRhodX - Density change with XFe
            K0     - Isothermal bulk modulus / Pa
            dKdX   - Change of K0 with XFe / Pa
            dKdP   - Pressure derivative of K0 / Pa/Pa
            dKdPdX - Change of dKdP with XFe
            dKdT   - Change of K0 with T / Pa/K
            mu0    - Shear modulus / Pa
            dmudX  - Change of mu0 with XFe
            dmudP  - Pressure derivative of mu0 / Pa/Pa
            dmudT  - Temperature derivative of mu0 / Pa/K
            alpha0 - Parameterisation for alpha(T) from Saxena and Shen (1992)
            alpha1 - Parameterisation for alpha(T) from Saxena and Shen (1992)
            alpha2 - Parameterisation for alpha(T) from Saxena and Shen (1992)
            alpha3 - Parameterisation for alpha(T) from Saxena and Shen (1992)

        Parameters:

        * FileIn : String
            Input file name.

        * delim : String
            String delimiter.

        Returns:

        * DBase : Numpy Array
            Structured numpy array containing database information.
        """
        if FileIn is None:
            FileIn = self.PyPath + self.FolderSep + 'MinDB.csv'
            if self.Verbose:
                print("Importing mineral properties from", FileIn)
        DBase = np.genfromtxt(FileIn, names=True, dtype=None, delimiter=delim,
                              encoding='utf8')
        return DBase

    def LoadFile(self):
        """
        Load the input file given by `self.FileIn`. The structure of input file
        must be `X Y Z Vel`.
        """
        print("Loading input file:", self.FileIn)
        print("> Velocity scale factor is", self.scaleV)
        self.DataRaw = np.loadtxt(self.FileIn)
        # Remove rows below maximum depth, defined by self.MaxDepth
        DelRows = np.where(np.abs(self.DataRaw[:, 2]) >= np.abs(self.MaxDepth))
        if len(DelRows[0]) > 0:
            print("> Removing depth values beyond", self.MaxDepth, "m.")
        self.DataRaw = np.delete(self.DataRaw, DelRows, axis=0)
        self.Depths = np.unique(self.DataRaw[:, 2])
        self.DataRaw[:, 3] *= self.scaleV
        self.__check_vmin__()

    def ReadArgs(self):
        """
        Read input arguments from console and define an output filename if none
        was specified.
        """
        if len(sys.argv) == 1:
            self.Usage()

        if sys.argv[1] == '-h' or sys.argv[1] == '--help':
            self.ShowHelp()
        else:
            self.FileIn = sys.argv[1]
        n_args = len(sys.argv)

        i = 2
        while i < n_args:
            if sys.argv[i] == '-AlphaT':
                self.UseAlpha = 'AlphaT'
            elif sys.argv[i] == '-AlphaPT':
                self.UseAlpha = 'AlphaPT'
            elif sys.argv[i] == '-comp':
                self.LoadMineralogy(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-dT':
                self.dT = float(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-NN':
                self.NN = True
            elif sys.argv[i] == '-out':
                self.FileOut = sys.argv[i+1]
                i += 1
            elif sys.argv[i] == '-scaleV':
                self.scaleV = float(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-setQ':
                self.SetQMode(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-type':
                self.SetVelType(sys.argv[i+1])
                i += 1
            elif sys.argv[i] == '-verbose' or sys.argv[i] == '-v':
                self.Verbose = True
            elif sys.argv[i] == '-xfe' or sys.argv[i] == '-XFe':
                self.SetXFe(XFe=sys.argv[i+1])
                i += 1
            else:
                print("Unknown argument", sys.argv[i])
                self.Usage()
            i += 1

        if self.FileOut is None:
            self.FileOut = self.FileIn[:-4] + '_RhoT.dat'

        self.TestArgs()

    def SaveFile(self):
        """
        Save results to an output file.
        """
        StrComment = "# "
        Output = np.empty([self.DataRaw.shape[0], 6])
        Output[:, 0:4] = self.DataRaw[:, 0:4]
        Output[:, 4] = self.Result_T - self.TKelvin
        Output[:, 5] = self.Result_Rho
        Output_h = "Temperature output\n"
        if self.FileIn:
            Output_h += "Input file: " + str(self.FileIn) + "\n"
        Output_h += "Velocity scale factor: " + str(self.scaleV) + "\n"
        Output_h += "Mantle composition:\n"
        for i in self.Mineralogy:
            Output_h += str(i['phase']) + " - " + str(i['fraction']) + "\n"
        Output_h += "XFe - " + str(self.XFe) + "\n"
        if self.SimpleP:
            Output_h += "Pressure calculation: Simplified\n"
        else:
            Output_h += "Pressure calculation: AK135\n"
        Output_h += "Wave frequency (Omega) / Hz: " + str(self.f) + "\n"
        Output_h += "Anelasticity parameters: " + self.Qmode + "\n"
        if self.UseAlpha == 'Alpha':
            Output_h += "Alpha depending on: Nothing\n"
        elif self.UseAlpha == 'AlphaT':
            Output_h += "Alpha depending on: Temperature\n"
        elif self.UseAlpha == 'AlphaPT':
            Output_h += "Alpha depending on: Temperature, Pressure\n"
        Output_h += "Columns:\n"
        Output_h += "1 - X\n"
        Output_h += "2 - Y\n"
        Output_h += "3 - Z / masl\n"
        Output_h += "4 - V_" + self.VelType + " / m/s\n"
        Output_h += "5 - T_syn / degC\n"
        Output_h += "6 - Rho / kg/m3"
        if self.NN:
            Output_h = str(self.DataRaw.shape[0])
            StrComment = ''
        fmtstring = '%f %f %f %.2f %.1f %.1f'
        print('Saving results to', self.FileOut)
        np.savetxt(self.FileOut, Output, header=Output_h, fmt=fmtstring,
                   comments=StrComment)
        print('> Done!')

    def SetAlpha(self, mode):
        """Set the dependency of the thermal expansion coefficient.

        If `const`, will use a pressure- and temperature-independent expansion
        coefficient. The value used is `alpha0` (see `MinDB.csv`). If `T`, the
        temperature dependency after Saxena and Shen (1992) will be applied. If
        `PT`, will use a pressure- and temperature-dependent thermal expansion
        coefficient that was obtained from Hacker and Abers (2004) using the
        visual basic script in the repository folder `AlphaPT`.

        Parameters
        ----------
        mode : str
            Valid modes are `const`, `T` and `PT`
        """
        valid_modes = {
            'const': 'Alpha',
            'T': 'AlphaT',
            'PT': 'AlphaPT'
        }
        self.UseAlpha = valid_modes[mode]

    def SetMineralogy(self, assemblage):
        """
        Manually define the mineralogy

        Parameters
        ----------
        assemblage : dir
            dir with keys that correspond to the `phase` column in MinDB.csv
        """
        n_Phases = len(assemblage)
        phases = list(assemblage.keys())

        names = ['phase', 'fraction']
        formats = ['U10', 'f8']
        dtypes = {'names':names, 'formats':formats}

        self.Mineralogy = np.array(list(assemblage.items()), dtype=dtypes)
        self.__check_mineralogy__()
        self.AssignPhases()

    def SetQMode(self, mode):
        """Define the attenuation model

        Defines the QMode which is being used for attenuation.

        Parameters
        ----------
        mode : int, str
            `1` will use the model of Sobolev et al. (1996), `2` will use
            Berckhemer et al. (1982).
        """
        mode = str(mode)
        if mode == '1':
            self.a = 0.15
            self.A = 0.148
            self.H = 500000.
            self.V = 0.000020
            self.Qmode = 'Sobolev et al. (1996)'
        elif mode == '2':
            self.a = 0.25
            self.A = 0.0002
            self.H = 584000.
            self.V = 0.000021
            self.Qmode = 'Berckhemer et al. (1982)'
        else:
            msg = "Unknown attenuation mode " + str(mode)
            self.__raise__(ValueError, msg)
        print("Using Q after", self.Qmode)
        print("    a =", self.a)
        print("    A =", self.A)
        print("    H =", self.H, "J/mol")
        print("    V =", self.V, "m3/mol")

    def SetVelType(self, TypeString):
        """
        Define velocity type and frequency.

        Parameters:

        * TypeString : string
            S or P for velocity type
        """
        if TypeString == 'S' or TypeString == 's':
            self.VelType = 'S'
            self.omega = 2*math.pi*0.02
        elif TypeString == 'P' or TypeString == 'p':
            self.VelType = 'P'
            self.omega = 2*math.pi*1
        else:
            msg = "Unknown velocity type " + TypeString
            self.__raise__(ValueError, msg)

    def SetXFe(self, XFe=None):
        """
        Define the iron content of the rock.
        """
        if XFe is None:
            self.XFe = 0
        else:
            self.XFe = float(XFe)

    def ShowProgress(self, i_step=None, i_max=None):
        """
        Show progress of a calculation. To initialise use i=ShowProgress()
        where i is the loop counter. To display the progress during the loop
        and to add counts to i use i=ShowProgress(i, i_max) during the loop.

        Parameters:

        * i_step : int
            Calculation step

        * i_max : int
            Maximum calculation step

        Returns:

        * i : int
            Initialises i or returns i_step + 1
        """
        if i_step is None:
            sys.stdout.write("> Progress: 0%")
            i = 0
        elif i_step < i_max - 1:
            progress = np.round((float(i_step)/float(i_max)) * 100, 2)
            progress_last = np.round((float(i_step - 1)/float(i_max)) * 100, 2)
            if progress > progress_last:
                sys.stdout.write("\r> Progress: %d%%" % progress)
            i = i_step + 1
        else:
            sys.stdout.write("\r> Progress: 100%\n")
            i = None
        sys.stdout.flush()
        return i

    def TestArgs(self):
        """
        Test whether minimum requirements for code execution are given.
        """

        if self.FileIn is None:
            msg = "No input file name given.")
            self.__raise__(ValueError, msg)
        if self.VelType is None:
            msg = "Error: No velocity type given.")
            self.__raise__(ValueError, msg)
        if self.UseAlpha == 'Alpha':
            print("Using constant expansion coefficient")
        elif self.UseAlpha == 'AlphaT':
            print("Using T-dependent expansion coefficient")
        elif self.UseAlpha == 'AlphaPT':
            print("Using P/T-dependent expansion coefficient")

    def Usage(self):
        """
        Print short script usage
        """
        print()
        print("Usage: Conversion.py FileIn -type <P|S> [optional args]")
        print("    Optional arguments:")
        print("        -AlphaT")
        print("        -AlphaPT")
        print("        -dT <val>")
        print("        -comp <Filename>")
        print("        -h | --help")
        print("        -NN")
        print("        -out <FileOut>")
        print("        -scaleV <value>")
        print("        -setQ <1|2>")
        print("        -v | -verbose")
        print("        -XFe <val>")
        print("    For more help:", sys.argv[0], "--help")
        print()
        sys.exit()


if __name__ == "__main__":
    Instance = MantleConversion()
    Instance.ReadArgs()
    Instance.LoadFile()
    Instance.FillTables()
    Instance.CalcPT()
    Instance.SaveFile()
