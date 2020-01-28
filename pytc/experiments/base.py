__description__ = \
"""
Classes for loading experimental data and associating those data with a
model.

Units:
    Volumes are in microliters
    Temperatures are in Kelvin
    Concentrations are in molar
    Energy is `units`, where `units` is specified when instantiating the
    PytcExperiment class.  It must be a in the AVAIL_UNITS dictionary.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-20"

import numpy as np
import pandas as pd

import random, string, os, re

class PytcExperiment:
    """
    Class that holds an experimental measurement and a model that describes it.
    """

    AVAIL_UNITS = {"cal/mol":1.9872036,
                   "kcal/mol":0.0019872036,
                   "J/mol":8.3144598,
                   "kJ/mol":0.0083144598}

    def __init__(self,
                 data_file,
                 model,
                 units="cal/mol",
                 uncertainty=0.1,
                 **model_kwargs):
        """

        Parameters
        ----------

        data_file: string
            file containing experimental results
        model: PytcModel subclass
            PytcModel subclass to use for modeling
        units : string
            file units ("cal/mol","kcal/mol","J/mol","kJ/mol")
        uncertainty : float > 0.0
            uncertainty in integrated heats (set to same for all shots, unless
            specified in the data file).

        **model_kwargs: any keyword arguments to pass to the model.  Any
                        keywords passed here will override whatever is
                        stored in the data_file.
        """

        # record the data file to load and fit against
        self.data_file = data_file

        # model function (will be initialized later)
        self._model_class = model

        # Deal with units
        self._units = units
        try:
            self._R = self.AVAIL_UNITS[self._units]
        except KeyError:
            err = "units must be one of:\n"
            for k in self.AVAIL_UNITS.keys():
                err += "    {}\n".format(k)
            err += "\n"

            raise ValueError(err)

        # For numerical reasons, there should always be *some* uncertainty
        self._uncertainty = uncertainty
        if self._uncertainty == 0.0:
            self._uncertainty = 1e-12

        # Create unique id for this experiment
        r = "".join([random.choice(string.ascii_letters) for i in range(20)])
        self._experiment_id = "{}_{}".format(self.data_file,r)

        # Read file and initialize model. (These can be re-defined in
        # subclasses)
        self._read_file()
        self._initialize_model(**model_kwargs)

    def _read_file(self):
        """
        Read the data file.
        """

        self._read_df(self.data_file)

    def _initialize_model(self,**model_kwargs):
        """
        Initialize the model describing the experiment.
        """

        self._model_class(**model_kwargs)

    def _read_df(self,df):
        """
        Read a pandas data frame, sticking column names in as attributes to the
        class.
        """

        # Load in data frame
        if type(df) is not pd.DataFrame:
            if type(df) is str:
                if os.path.isfile(df):
                    if df[-4:].lower() == ".csv":
                        df = pd.read_csv(df)
                    elif df[-4:].lower in ["xlsx",".xls"]:
                        df = pd.read_excel(df)
                    else:
                        err = "file type for {} not recognized\n".format(df)
                        raise ValueError(err)
                else:
                    err = "file {} not found.\n".format(df)
                    raise FileNotFoundError(err)
            else:
                err = "df should either be a dataframe or string pointing to a file\n"
                raise ValueError(err)

        # ---------------------------
        # Sanitize data frame columns
        # ---------------------------

        # Stuff we're going to want to remove
        p = re.compile("[: .]")

        # Go through each columns
        new_columns = []
        times_column_seen = {}
        for c in df.columns:

            # Split on stuff to remove
            splits = [s for s in p.split(c) if s != ""]
            new_column = "_".join(splits)

            # If we ended up creating a newly duplicated column, append a
            # number to it.,
            try:
                times_column_seen[new_column] += 1
                new_column = "{}_{}".format(new_column,times_column_seen[new_column]-1)
            except KeyError:
                times_column_seen[new_column] = 1

            new_columns.append(new_column)

        # Rename columns
        df.columns = new_columns

        # Grab the data frame
        self._df = df.copy()

        # Make columns accessible as attributes to the class
        for k in self._df.columns:
            try:
                self.__dict__[k]
                err = "input data frame uses a reserved name ({})]\n".format(k)
                raise ValueError(err)
            except KeyError:
                self.__dict__[k] = self._df[k]

    def plot(self,
             fig=None,ax=None,
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

    @property
    def obs(self):
        """
        """

        return np.array([])

    @property
    def param_values(self):
        """
        Values of fit parameters.
        """

        return self._model.param_values

    @property
    def param_stdevs(self):
        """
        Standard deviations on fit parameters.
        """

        return self._model.param_stdevs

    @property
    def param_ninetyfives(self):
        """
        95% confidence intervals on fit parmeters.
        """

        return self._model.param_ninetyfives

    @property
    def model(self):
        """
        Fitting model.
        """

        return self._model

    @property
    def experiment_id(self):
        """
        Return a unique experimental id.
        """

        return self._experiment_id

    @property
    def R(self):
        """
        Experiment gas constant.
        """

        return self._R

    @property
    def units(self):
        """
        Units for file.
        """

        return self._units

    @units.setter
    def units(self,units):
        """
        Change the units.
        """

        # Deal with units
        self._units = units
        try:
            self._R = self.AVAIL_UNITS[self._units]
        except KeyError:
            err = "units must be one of:\n"
            for k in self.AVAIL_UNITS.keys():
                err += "    {}\n".format(k)
            err += "\n"

            raise ValueError(err)
