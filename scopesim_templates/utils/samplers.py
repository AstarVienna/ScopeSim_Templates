"""A collection of samplers for generating distributions."""

import numpy as np


def randomvariate(pdf, n, xmin, xmax):
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
    TODO: Do optimization, e.g. using yield, etc
    """
    x = np.linspace(xmin, xmax, 10000)
    y = pdf(x)
    pmin = y.min()
    pmax = y.max()
    naccept = 0
    ntrial = 0
    ran = []
    while naccept < n:
        x = np.random.uniform(xmin, xmax)
        y = np.random.uniform(pmin, pmax)
        if y <= pdf(x):
            ran.append(x)
            naccept = naccept + 1
            #  print pdf(x)
            ntrial = ntrial + 1

    ran = np.asarray(ran)
    return ran, ntrial


def np_rvs(pdf, n, xmin, xmax, sampling=10000):
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

    Returns
    -------
    TYPE
        np.array with the values.

    """
    x = np.linspace(xmin, xmax, sampling)
    y = pdf(x)
    pmin = y.min()
    pmax = y.max()

    rel_area = (xmax - xmin) * (pmax - pmin) / np.trapz(y, x)
    m = int(n * rel_area) * 2
    ran = np.array([])

    while ran.size < n:
        X = np.random.uniform(xmin, xmax, m)
        Y = np.random.uniform(pmin, pmax, m)
        ran = np.append(ran, X[Y <= pdf(X)])  # the evaluation pdf(X) is the slow part for complex distributions

    return ran[:n]


def _metropolis_sampler(pdf, nsamples, xmin, xmax):
    """
    Taken from here:
    https://stackoverflow.com/questions/51050658/how-to-generate-random-numbers-with-predefined-probability-distribution

    Usage
        s1 = Sersic1D(amplitude=1, r_eff=5, n=4)
        p = lambda r: s1(r)
        samples = list(metropolis_sampler(p, N))
    """
    x = np.random.uniform(xmin, xmax)  # start somewhere

    for i in range(nsamples):
        trial = np.random.uniform(xmin, xmax)  # random neighbour from the proposal distribution
        acceptance = pdf(trial) / pdf(x)

        if np.random.uniform() < acceptance:
            x = trial

        yield x


def metropolis_sampler(pdf, nsamples, xmin, xmax):
    """
    This is a version of a metropolis sampler.

    Depending of the characteristics of the distribution it can be slower or
    comparable in speed to rvs
    """
    result = list(_metropolis_sampler(pdf=pdf, nsamples=nsamples,
                                      xmin=xmin, xmax=xmax))
    return np.array(result)
