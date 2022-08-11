'''
Generate multivariate von Mises Fisher samples.
This solution originally appears here:
http://stats.stackexchange.com/questions/156729/sampling-from-von-mises-fisher-distribution-in-python
Also see:
Sampling from vMF on S^2:
    https://www.mitsuba-renderer.org/~wenzel/files/vmf.pdf
    http://www.stat.pitt.edu/sungkyu/software/randvonMisesFisher3.pdf

This code was taken from the following project:
https://github.com/clara-labs/spherecluster
'''
import numpy as np


__all__ = ['sample_vMF']

def sample_vMF(mu, kappa, num_samples=1):
    """Generate N-dimensional samples from von Mises Fisher
    distribution around center mu \in R^N with concentration kappa.
    mu and kappa may be vectors,
    mu should have shape (dim, N), kappa should be vector of length N.
    num_samples is only used if mu is of shape (3, 1) and kappa is scalar.
    """

    if isinstance(kappa, np.ndarray):
        dim = mu.shape[0]
        assert mu.shape[1] == kappa.size
    else:
        dim = mu.size
        mu = np.vstack(num_samples*(mu,)).T
        kappa = np.full(num_samples, kappa)

    # sample offset from center (on sphere) with spread kappa
    w = _sample_weight(kappa, dim)

    # sample a point v on the unit sphere that's orthogonal to mu
    v = _sample_orthonormal_to(mu)

    # compute new point
    result = v * np.sqrt(1. - w**2) + w * mu
    #returns result.T to be backwards compatible
    return result.T

def _sample_weight(kappa, dim):
    """Rejection sampling scheme for sampling distance from center on
    surface of the sphere.
    """
    dim = dim - 1  # since S^{n-1}
    try:
        size = kappa.size
    except AttributeError:
        size = 1
    b = dim / (np.sqrt(4. * kappa**2 + dim**2) + 2 * kappa)
    x = (1. - b) / (1. + b)
    c = kappa * x + dim * np.log(1 - x**2)
    w = np.zeros(size)
    idx = np.zeros(size, dtype=int)
    while True:
        #TODO make more efficient
        #samples for every kappa and only uses the according sample
        #only if it is still needed, else discarded
        where_zero = np.nonzero(idx == 0)
        done = np.all(idx)
        if done:
            return w
        z = np.random.beta(dim / 2., dim / 2., size=size)
        _w = (1. - (1. + b) * z) / (1. - (1. - b) * z)
        u = np.random.uniform(low=0, high=1, size=size)
        _idx = np.nonzero(kappa * _w + dim * np.log(1. - x * _w) - c >= np.log(u))
        if _idx[0].size == 0:
            continue
        else:
            #sort into correct slots of w
            w[where_zero] = _w[where_zero]
            idx[_idx] = 1

def _sample_orthonormal_to(mu):
    """Sample point on sphere orthogonal to mu."""
    v = np.random.randn(*reversed(mu.shape))
    proj_mu_v =  mu * np.diag(np.matmul(v, mu)) / np.linalg.norm(mu, axis=0)
    orthto = v - proj_mu_v.T
    return orthto.T / np.linalg.norm(orthto.T, axis=0)

