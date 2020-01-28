
from .base import PytcExperiment

import numpy as np


class GenericExperiment(PytcExperiment):

    def __init__(self,
                 data_file,
                 model,
                 obs_column="obs",
                 obs_stdev_column=None,
                 uncertainty=0.1,
                 units="cal/mol",
                 **model_kwargs):

        super().__init__(data_file,model,units,uncertainty,**model_kwargs)

        self._obs_column = obs_column
        self._obs_stdev_column = obs_stdev_column

        # Read in observations
        try:
            self._obs = self._df[obs_column]
        except KeyError:
            err = "The column '{}' is not in the spreadsheet.\n".format(self._obs_column)
            raise ValueError(err)

        # Read in observation standar deviations
        if self._obs_stdev_column is not None:

            try:
                self._obs_stdev = self._df[obs_stdev_column]
            except KeyError:
                err = "The column '{}' is not in the spreadsheet.\n".format(self._obs_stdev_column)
                raise ValueError(err)

        # Just assign uncertainty
        else:
            self._obs_stdev = self._uncertainty*np.ones(len(self._obs))

    def _initialize_model(self,**model_kwargs):

        self._model = self._model_class(df=self._df,**model_kwargs)

    @property
    def predicted(self):
        return self._model.predicted

    @property
    def obs(self):
        return self._obs

    @property
    def obs_stdev(self):
        return self._obs_stdev
