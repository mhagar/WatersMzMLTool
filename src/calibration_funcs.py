"""
Functions that can be fit to calibration scans
"""
from numpy import log10

# Intended for use with scipy.optimize.curve_fit
def parabola(X, a, b, c):
    """
    X: array of mz values
    """
    return a + (b * X) + (c * X**2)


def fifth_order_polynomial(X, a, b, c, d, e, f):
    """
    X: array of mz values
    """
    return a + (b * X) + (c * X**2) + (d * X**3) + + (e * X**4) + (f * X**5)


def parabola_w_intensity_term(X, a, b, c, d):
    """
    X: array of mz values and intsy values
    """
    mz = X[:, 0]
    # intsy = X[:, 1]
    return a + (b * mz) + (c * mz**2)
    # return a + (b * mz) + (c * mz**2) + (d * log10(intsy))