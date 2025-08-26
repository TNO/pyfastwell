import numpy as np
from matplotlib import pyplot as plt
from typing import Any, Optional, Callable, List, Dict

def get_random_value_double_triangular(minimum: float, maximum: float, most_likely: float) -> float:
    """
    Calculate random value using a double-triangular distribution

    Parameters
    ----------

    minumum: Minimum value

    maximum: Maximum value

    most_likely: Most    likely value

    Returns
    -------
     Random value between min and max
    """
    random_value = np.random.random()

    if random_value <= 0.5:
        random_value = np.sqrt(2 * random_value)
    else:
        random_value = (4 - np.sqrt(8 - (8 * random_value))) / 2

    # Scale to input ranges
    if random_value < 1:
        random_value = minimum + (most_likely - minimum) * (random_value)
    else:
        random_value = most_likely + (maximum - most_likely) * (random_value - 1)
    return random_value

class Stochasticmodel(object):
    """
      Stochastic model for parameter sampling and forward modeling.

      This class manages model parameters, samples them according to specified distributions,
      and applies them to a forward model for ensemble simulations.


      Notes
      -----
      - Locked parameters remain constant across the ensemble.
      - Free parameters are sampled according to their specified distributions.
      - The class supports normal, uniform, and triangular distributions for sampling.
    """

    O_TIME = 0  # column number of time in the forward model output

    LOCKED = 'locked'
    FREEPARAM = 'free_param'
    DISTPARAM = 'free_dist'



    def __init__(self, argdict: Optional[Dict[str, Any]] = None) -> None:
        """
        Initialize the Stochasticmodel object.

        Parameters:
        ----------
        argdict: dictionary with keys 'locked', 'free_param', and 'free dist' containing the model parameters

        Notes
        ______
        The model parameters are subdivided in 'locked' parameters, which typically include the time array for which simulation is performed
        and additional 'locked' parameters which are not variable for the ensemble. The 'free_parameters' contains a list of parameter names which are allowed to vary
        and are sampled for the Stochasticmodel corresponding to the 'free_dist' distribition specifications

        """
        if argdict is None:
            self.paramdict = {}
        self.argdict = argdict
        # check if dict contains 'locked' and 'free_param' keys
        if Stochasticmodel.LOCKED not in self.argdict:
            self.argdict[Stochasticmodel.LOCKED] = {}
        if Stochasticmodel.FREEPARAM not in self.argdict:
            self.argdict[Stochasticmodel.FREEPARAM] = []
        if Stochasticmodel.DISTPARAM not in self.argdict:
            self.argdict[Stochasticmodel.DISTPARAM] = []

    def get_from_dict(self, name: str, param: np.ndarray) -> Any:
        """
        Get a parameter from the dictionary based on its name, either from 'locked' or 'free_param'
        If the name is not found in 'free_param', it will look in 'locked'.
        If the name is not found in either, it will raise an exception.
        This method is used to retrieve the parameter values for the ensemble realisations.

        Parameters
        ----------
        name: name of the parameter to be retrieved from the dictionary (either 'locked' or 'free_param'), 'free_param' takes prevalence

        param: ensembles realisation corresponding to 'free_param' sampled values

        Returns
        -------
        value of the parameter corresponding to the name, either from 'locked' or from 'free_param'
        """
        try:
            param_ind = -1
            for index, k in enumerate(self.argdict[Stochasticmodel.FREEPARAM]):
                if k == name:
                    param_ind = index
                    break
            if param_ind == -1:
                raise KeyError
            else:
                ret = np.array(param)[param_ind]
        except KeyError:
            try:
                g = self.argdict[Stochasticmodel.LOCKED]
                ret = self.argdict[Stochasticmodel.LOCKED][name]
            except KeyError :
                print('Stochasticmodel.get_from_dict()  name not found, name = ', name)
                raise Exception("An error occurred in get_from_dict")
        return ret

    def getmedianparams(self)-> np.ndarray:
        """
        Get the median values of the free parameters

        Returns
        -------
        median values of the free parameters in one dimensional array
        """

        freeparam = self.argdict[Stochasticmodel.FREEPARAM]
        distparam = self.argdict[Stochasticmodel.DISTPARAM]
        median = []
        for param in freeparam:
            if distparam[param]['dist'] == 'normal':
                median.append(distparam[param]['values'][0])
            elif distparam[param]['dist'] == 'uniform':
                median.append((distparam[param]['values'][0] + distparam[param]['values'][1]) / 2.0)
            elif distparam[param]['dist'] == 'triangular':
                median.append(distparam[param]['values'][2])
        return np.array(median)

    def generate_ensemble(self, nsamples: Optional[int] = None) -> List[np.ndarray]:
        """
        Generate an ensemble of model parameters based on the distribution specifications in the dictionary

        Parameters
        ----------
        nsamples: number of samples to generate

        Returns
        -------
        list of parameter arrays (dimension equal to number of free parameters), list size equal to nsamples

        """
        if nsamples is None:
            nsamples = self.nsamples
        else:
            self.nsamples = nsamples

        freeparam = self.argdict[Stochasticmodel.FREEPARAM]
        distparam = self.argdict[Stochasticmodel.DISTPARAM]
        m = []
        for i in range(nsamples):
            samples = []
            #    nsamples):
            #for i in len(freeparam)
            #    # variables in forward model are assigned a value (priors are created)(mean, stdev, shape of array)
            for param in freeparam:
                if distparam[param]['dist'] == 'normal':
                    sample = np.random.normal(distparam[param]['values'][0], distparam[param]['values'][1], 1)[0]
                elif distparam[param]['dist'] == 'uniform':
                    sample = np.random.uniform(distparam[param]['values'][0], distparam[param]['values'][1], 1)[0]
                elif distparam[param]['dist'] == 'triangular':
                    sample = get_random_value_double_triangular(distparam[param]['values'][0],distparam[param]['values'][1],distparam[param]['values'][2])
                samples.append(sample)
            mr = np.array(samples)  # array of all values per parameter in prior
            m.append(mr)  # list of prior array

        return m

    def run_ensemble(self, m: List[np.ndarray], forward: Callable[[np.ndarray], Any]) -> List[np.ndarray]:
        """
        Run the forward model for the ensemble of parameters

        Parameters
        ----------
        m: ensemble of parameters, corresponding to the free parameters and created by generate_ensemble
        forward: forward model function

        Returns
        -------
        results of the forward model for the ensemble
        """

        results = []
        for i in range(len(m)):
            results.append(forward(m[i]))
        return results



    def plot_distribution(self, results:np.ndarray, name:str):
        """
        Plot the distribution of a parameter forward model results

        Parameters:
        ----------
        results: results of the forward model for the ensemble ndarray

        name: name of the plot

        Returns
        -------
        None
        """

        plt.figure(figsize=(10, 6))
        plt.hist(results, bins=30, density=True, alpha=0.6, color='g')
        plt.xlabel(name)
        plt.ylabel('Density')
        plt.title('Histogram of ' + name)
        plt.grid(True)
        plt.show()

    def expectation_plot(
            self,
            results: np.ndarray,
            name: str,
            percent: bool = True,
            expectation: bool = True,
            pvals: Optional[List[float]] = None,
    ) -> None:

        """
        Plot the expectation of a parameter forward model results

        Parameters
        ----------
        results: results of the forward model for the ensemble ndarray

        name: name of the plot

        percent: if True, the y-axis is in percent, otherwise in fraction

        expectation: if True, the y-axis is reversed interpreting the results as expectation

        pvals: pvalues to be plotted (list)


        Returns
        -------
        None
        """


        # Sort the sample values
        sorted_samples = np.sort(results)

        ymax =1
        if percent:
            ymax =100

        if expectation:
            sorted_samples = sorted_samples[::-1]

        # Compute the cumulative density
        cdf = np.arange(1, len(sorted_samples) + 1) / len(sorted_samples)
        cdf *= ymax
        plt.figure(figsize=(10, 6))
        plt.plot(sorted_samples, cdf, marker='none', linestyle='-')
        plt.xlabel(name)
        plt.ylabel('Cumulative Density')
        plt.ylim(0,ymax)
        plt.title('Expectation plot of ' + name)
        plt.grid(True)


        if (pvals!=None):
            for pval in pvals:
                v = np.interp(pval, cdf, sorted_samples)
                plt.plot(v, pval, 'ro')
                label =  'p' + str(pval) + ' = ' + "{:.2f}".format(v)
                plt.annotate(label, (v, pval), textcoords="offset points", xytext=(5, 5), ha='center')



        plt.show()
