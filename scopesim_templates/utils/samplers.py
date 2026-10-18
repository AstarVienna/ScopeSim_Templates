"""A collection of samplers for generating distributions."""

import numpy as np

from .general_utils import get_rng


def randomvariate(pdf, n, xmin, xmax, seed=None):
    """
    Create random numbers according to an arbitrary PDF.

    Uses the rejection method for generating random numbers derived from an
    arbitrary probability distribution. For reference, see Bevington's book,
    page 84. Based on rejection*.py.
    Usage:
    randomvariate(P,N,xmin,xmax)
    where
    P : probability distribution function from which you want to generate
    random numbers
    N : desired number of random values
    xmin,xmax : range of random numbers desired
    Returns:
    the sequence (ran,ntrials) where
    ran : array of shape N with the random variates that follow the input P
    ntrials : number of trials the code needed to achieve N
    Here is the algorithm:
    - generate x' in the desired range
    - generate y' between Pmin and Pmax (Pmax is the maximal value of your pdf)
    - if y'<P(x') accept x', otherwise reject
    - repeat until desired number is achieved
    Rodrigo Nemmen
    Nov. 2011
    seed : optional seed or np.random.Generator
    """
    rng = get_rng(seed)
    x = np.linspace(xmin, xmax, 10000)
    y = pdf(x)
    pmin = y.min()
    pmax = y.max()
    naccept = 0
    ntrial = 0
    ran = []
    while naccept < n:
        x = rng.uniform(xmin, xmax)
        y = rng.uniform(pmin, pmax)
        ntrial = ntrial + 1
        if y <= pdf(x):
            ran.append(x)
            naccept = naccept + 1

    ran = np.asarray(ran)
    return ran, ntrial


def np_rvs(pdf, n, xmin, xmax, sampling=10000, seed=None):
    """
    Numpy adaptation of randomvariate.

    This works in only one dimension
    This version is 100 to 1000 times faster than randomvariate for reasonable
    sized samples

    Parameters
    ----------
    pdf : TYPE
        Probability density function to sample.
    n : TYPE
        Number of returned random values.
    xmin, xmax : TYPE
        Range of random numbers desired.
    sampling : TYPE, optional
        Initial sampling of the parameter space. Increase for complex
        distributions. The default is 10000.
    seed : int, np.random.Generator, optional
        Seed or generator for reproducible samples. Default None.

    Returns
    -------
    TYPE
        np.array with the values.

    """
    x = np.linspace(xmin, xmax, sampling)
    y = pdf(x)
    pmin = y.min()
    pmax = y.max()

    rng = get_rng(seed)
    rel_area = (xmax - xmin) * (pmax - pmin) / np.trapezoid(y, x)
    m = int(n * rel_area) * 2
    ran = np.array([])

    while ran.size < n:
        X = rng.uniform(xmin, xmax, m)
        Y = rng.uniform(pmin, pmax, m)
        ran = np.append(ran, X[Y <= pdf(X)])  # the evaluation pdf(X) is the slow part for complex distributions

    return ran[:n]


def _metropolis_sampler(pdf, nsamples, xmin, xmax, rng):
    """
    Taken from here:
    https://stackoverflow.com/questions/51050658/how-to-generate-random-numbers-with-predefined-probability-distribution

    Usage
        s1 = Sersic1D(amplitude=1, r_eff=5, n=4)
        p = lambda r: s1(r)
        samples = list(metropolis_sampler(p, N))
    """
    x = rng.uniform(xmin, xmax)  # start somewhere

    for i in range(nsamples):
        trial = rng.uniform(xmin, xmax)  # random neighbour from the proposal distribution
        acceptance = pdf(trial) / pdf(x)

        if rng.uniform() < acceptance:
            x = trial

        yield x


def metropolis_sampler(pdf, nsamples, xmin, xmax, seed=None):
    """
    This is a version of a metropolis sampler.

    Depending of the characteristics of the distribution it can be slower or
    comparable in speed to rvs
    """
    result = list(_metropolis_sampler(pdf=pdf, nsamples=nsamples,
                                      xmin=xmin, xmax=xmax,
                                      rng=get_rng(seed)))
    return np.array(result)
