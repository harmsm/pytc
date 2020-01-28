

def plot(global_fit,
         fig=None,ax=None,
         correct_molar_ratio=False,subtract_dilution=False,
         normalize_heat_to_shot=False,color_list=None,
         data_symbol="o",linewidth=1.5,num_samples=100,
         markersize=8,
         each_expt_own_plot=False):
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

    # Make list of colors
    if color_list is None:
        N = len(self._expt_list_stable_order)
        color_list = [plt.cm.brg(i/N) for i in range(N)]

    # Sanity check on colors
    if len(color_list) < len(self._expt_list_stable_order):
        err = "Number of specified colors is less than number of experiments.\n"
        raise ValueError(err)

    try:
        # If there are samples:
        if len(self._fitter.samples) > 0:
            s = self._fitter.samples
            these_samples = s[np.random.randint(len(s),size=num_samples)]
        else:
            these_samples = [self._fitter.estimate]
    except AttributeError:

        # If fit has not been done, create dummy version
        self._prep_fit()
        these_samples = [np.array(self._flat_param)]

    # If there are multiple samples, assign them partial transparency
    if len(these_samples) == 1:
        alpha = 1.0
    else:
        alpha = 0.1


    # plots holds [fig,ax] pairs for each plot we're making, map_expt_to_plot
    # indexes each epxeriment to its appropraite plot.
    plots = []
    map_expt_to_plot = []
    _map_type_to_plot = {}
    _plot_counter = 0
    for j, expt_name in enumerate(self._expt_list_stable_order):

        # If every experiment gets its own plot, do so
        if each_expt_own_plot:
            plots.append([None,None])
            map_expt_to_plot.append(_plot_counter)
            _plot_counter += 1

        # Otherwise, group experiments of the same type together
        else:
            expt_type = type(self._expt_dict[expt_name])

            # If we've already seen this experiment type, find the plot
            # we should use for this experiment
            try:
                index = _map_type_to_plot[expt_type]
                map_expt_to_plot.append(index)

            # If we haven't already seen the expermient type, make a new
            # plot to hold it.
            except KeyError:
                _map_type_to_plot[expt_type] = _plot_counter
                map_expt_to_plot.append(_plot_counter)
                plots.append([None,None])
                _plot_counter += 1

    for i, s in enumerate(these_samples):

        # Update calculation for this sample
        self._y_calc(s)
        for j, expt_name in enumerate(self._expt_list_stable_order):

            fig, ax = plots[map_expt_to_plot[j]]

            # Extract fit info for this experiment
            e = self._expt_dict[expt_name]
            fig, ax = e.plot(fig,ax,
                             color=color_list[j],alpha=alpha,
                             draw_fit=True,draw_expt=False)

            # If this is the last sample, plot the experimental data
            if i == len(these_samples) - 1:
                fig, ax = e.plot(fig,ax,
                                 color=color_list[j],
                                 draw_fit=False,draw_expt=True)

            # update fig and ax
            plots[map_expt_to_plot[j]] = [fig, ax]

    for i in range(len(plots)):
        plots[i][0].set_tight_layout(True)

    return plots


def corner_plot(global_fit,filter_params=("competent","dilution","intercept","heat")):
    """
    Create a "corner plot" that shows distributions of values for each
    parameter, as well as cross-correlations between parameters.

    Parameters
    ----------
    filter_params: do not plot parameters that have these values somewhere
                   in their names.
    """

    try:
        return global_fit.fitter.corner_plot(filter_params)
    except AttributeError:
        # If the fit has not been done, return an empty plot
        dummy_fig = plt.figure(figsize=(5.5,6))
        return dummy_fig
