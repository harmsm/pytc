__description__ = \
"""
Classes for loading experimental ITC data and associating those data with a
model.

Units:
    Volumes are in microliters
    Temperatures are in Kelvin
    Concentrations are in molar
    Energy is `units`, where `units` is specified when instantiating the
    ITCExperiment class.  It must be a in the AVAIL_UNITS dictionary.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

from ..base import PytcExperiment

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

import random, string, os

class BaseITCExperiment(PytcExperiment):
    """
    Class that holds an experimental ITC measurement and a model that describes it.
    """

    AVAIL_UNITS = {"cal/mol":1.9872036,
                   "kcal/mol":0.0019872036,
                   "J/mol":8.3144598,
                   "kJ/mol":0.0083144598}

    def __init__(self,data_file,model,shot_start=1,units="cal/mol",
                 uncertainty=0.1,**model_kwargs):
        """

        Parameters
        ----------

        data_file: string
            integrated heats file written out by origin software.
        model: ITCModel subclass instance
            ITCModel subclass to use for modeling
        shot_start: int
            what shot to use as the first real point.  Shots start at 0, so
            default=1 discards first point.
        units : string
            file units ("cal/mol","kcal/mol","J/mol","kJ/mol")
        uncertainty : float > 0.0
            uncertainty in integrated heats (set to same for all shots, unless
            specified in something like NITPIC output file).

        **model_kwargs: any keyword arguments to pass to the model.  Any
                        keywords passed here will override whatever is
                        stored in the data_file.
        """


        # ITC specific information
        self._shot_start = shot_start

        # Run parent init
        super().__init__(data_file,model,units,uncertainty,**model_kwargs)


    def _initialize_model(self,**model_kwargs):

        # Initialize model using information read from heats file
        self._model = self._model_class(S_cell=self.stationary_cell_conc,
                                        T_syringe=self.titrant_syringe_conc,
                                        cell_volume=self.cell_volume,
                                        shot_volumes=self._shots,
                                        **model_kwargs)

    @property
    def predicted(self):
        """
        Return heats calculated by the model with parameters defined in params
        dictionary.
        """

        if len(self._model.predicted) == 0:
            return np.array(())

        return self._model.predicted[self._shot_start:]

    @property
    def dilution_heats(self):
        """
        Return dilution heats calculated by the model with parameters defined
        in params dictionary.
        """

        if len(self._model.dilution_heats) == 0:
            return np.array(())

        return self._model.dilution_heats[self._shot_start:]


    @property
    def shot_start(self):
        """
        Starting shot to use.
        """

        return self._shot_start

    @shot_start.setter
    def shot_start(self,value):
        """
        Change starting shot.
        """

        self._shot_start = value

    @property
    def obs(self):
        """
        Return experimental heats.
        """
        return self._heats[self._shot_start:]

    @obs.setter
    def obs(self,obs):
        """
        Set the heats.
        """

        self._heats[self._shot_start:] = obs[:]

    @property
    def obs_stdev(self):
        """
        Standard deviation on the uncertainty of the heat.
        """

        return self._heats_stdev[self._shot_start:]

    @obs_stdev.setter
    def obs_stdev(self,obs_stdev):
        """
        Set the standard deviation on the uncertainty of the heat.
        """

        self._heats_stdev[self._shot_start:] = obs_stdev[:]

    @property
    def mol_injected(self):
        """
        Return the mols injected over shots.
        """

        # uL * mol/L * L/1e6 uL -> mol
        return self._shots[self._shot_start:]*self.titrant_syringe_conc*1e-6

    @property
    def mole_ratio(self):
        """
        Return the mole ratio of titrant to stationary.
        """
        return self._model.mole_ratio[self._shot_start:]


    def plot(self,
             fig=None,ax=None,
             correct_molar_ratio=False,
             subtract_dilution=False,
             normalize_heat_to_shot=False,
             draw_fit=False,
             draw_expt=True,
             color="black",
             data_symbol="o",
             markersize=8,
             linewidth=1.5,
             alpha=1.0):
        """
        Plot the experimental data and fit results.

        Parameters
        ----------
        correct_molar_ratio : bool
            correct the molar ratio using fx_competent
        subtract_dilution : bool
            subtract the heat of dilution
        normalize_heat_to_shot : bool
            divide heats by mols titrant injected
        color_list : list of things matplotlib can interpret as colors
            color of each series
        data_symol : character
            symbol to use to plot data
        linewidth : float
            width of line for fits
        num_samples : int
            number of samples to draw when drawing fits like Bayesian fits with
            multiple fits.

        Returns matplotlib Figure and AxesSubplot instances that can be further
        manipulated by the user of the API.
        """

        if fig is None and ax is None:

            # Make graph of appropraite size
            fig = plt.figure(figsize=(5.5,6))

            # Create two panel graph
            gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
            ax = []
            ax.append(fig.add_subplot(gs[0]))
            ax.append(fig.add_subplot(gs[1],sharex=ax[0]))

            # Clean up graphs
            for i in range(2):
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)

                ax[i].yaxis.set_ticks_position('left')
                ax[i].xaxis.set_ticks_position('bottom')

            # Add labels to top plot and remove x-axis
            u = self.units

            if normalize_heat_to_shot:
                ax[0].set_ylabel("heat per mol titrant ({})".format(u))
            else:
                new_u = u.split("/")[0]
                ax[0].set_ylabel("observed heat ({})".format(new_u))

            plt.setp(ax[0].get_xticklabels(), visible=False)

            # Add labels to the residuals plot
            m = self.mole_ratio
            ax[1].plot([np.min(m),np.max(m)],[0,0],"--",linewidth=1.0,color="gray")
            ax[1].set_xlabel("molar ratio (titrant/stationary)")
            ax[1].set_ylabel("residual")

        else:
            if fig is None or ax is None:
                err = "either both fig and ax must be specified or neither \n"
                err += "can be specified.\n"
                raise ValueError(err)

        # Extract fit info for this experiment
        mr = self.mole_ratio
        obs = self.obs
        obs_stdev = self.obs_stdev
        calc = self.predicted

        if len(calc) > 0:

            # Try to correct molar ratio for competent fraction
            if correct_molar_ratio:
                try:
                    mr = mr/self.param_values["fx_competent"]
                except KeyError:
                    pass

            # Subtract dilution is requested
            if subtract_dilution:
                obs = obs - self.dilution_heats
                calc = calc - self.dilution_heats

            if normalize_heat_to_shot:
                obs = obs/self.mol_injected
                calc = calc/self.mol_injected

        # Draw fit lines and residuals
        if draw_fit and len(self.predicted) > 0:

            marker_style = dict(color=color,
                                linestyle='-',
                                linewidth=linewidth,
                                alpha=alpha,
                                marker=data_symbol,
                                markersize=markersize*2,
                                markerfacecoloralt=color,
                                fillstyle="none")

            ax[0].plot(mr,calc,**marker_style)
            ax[1].plot(mr,(calc-obs),data_symbol,color=color,alpha=alpha,markersize=markersize)

        # If this is the last sample, plot the experimental data
        if draw_expt:
            ax[0].errorbar(mr,obs,obs_stdev,fmt=data_symbol,color=color,
                           markersize=markersize,linestyle='none')

        fig.set_tight_layout(True)

        return fig, ax
