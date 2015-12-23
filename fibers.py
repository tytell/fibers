import numpy as np
from scipy import signal, linalg, argsort


def ndautocorr(V, insz, window=None, beta=1, normalize=True):
    """
    Autocorrelation of a volume V with a smaller subvolume of itself

    :param V: Volume
    :param insz: Size of the smaller region.  Will be centered
    :param window: None, 'hamming', or 'kaiser'
    :param beta: Kaiser beta parameter
    :param normalize: If True, normalizes the volume by subtracting the mean and dividing by the standard deviation.
    :return: Returns the autocorrelation matrix
    """
    outsz = V.shape

    if normalize:
        V = (V-np.mean(V))/np.std(V)

    if window is not None:
        if window == 'hamming':
            windowfcn = np.hamming
        elif window == 'kaiser':
            windowfcn = lambda x: np.kaiser(x,beta=beta)
        else:
            raise ValueError('Unrecognized window function')

        W = np.ones((outsz[0]))
        for dim1,outsz1 in enumerate(outsz):
            # make the window on one axis
            w1 = windowfcn(outsz1)
            # set up a slice to broadcast along the current dimension
            ss = [np.newaxis for i in outsz]
            ss[dim1] = slice(None)
            # and multiply to build up the ND window
            W = W * w1[ss]

        V = V*W

    inrng = [slice(out1//2-in1//2, out1//2-in1//2 + in1)
            for (in1,out1) in zip(insz,outsz)]
    outrng = [slice(None,None,-1) for in1 in insz]

    C = signal.fftconvolve(V[outrng],V[inrng],mode='valid')

    return C


def hessian(C, ctr=None, di=1):
    """
    Returns the Hessian matrix of C, evaluated at the point ctr, given a spacing di.
    Uses 5th order central difference approximations to the derivative.

    :param C: 2 or 3D ndarray
    :param ctr: Point at which to evaluate the Hessian
    :param di: Spacing between points in C
    :return: Returns the Hessian matrix
    """

    # 5th order derivatives
    # first derivative
    a1 = np.array([-1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0])/60
    # second derivative
    a2 = np.array([2, -27, 270, -490, 270, -27, 2])/180.0

    # round down
    n = len(a1)//2
    if ctr is None:
        # round up if lengths are odd
        ctr = (np.array(C.shape)+1)//2

    ic = np.arange(-n,n+1)

    nd = C.ndim
    assert((nd == 2) or (nd == 3))

    H = np.zeros((nd,nd))

    # get second derivatives along the diagonal
    for i in xrange(nd):
        r = list(ctr)
        r[i] += ic
        H[i,i] = sum(C[tuple(r)] * a2) / di**2

    # cross derivatives off diagonal
    ic1 = ic[:,np.newaxis]
    ic2 = ic[np.newaxis,:]

    for i in xrange(nd):
        for j in xrange(i+1,nd):
            r = list(ctr)
            r[i] += ic1
            r[j] += ic2
            H[i,j] = np.sum( np.sum( C[tuple(r)] * a1[:,np.newaxis], axis=0) * a1) / di**2

            # mirror the cross derivatives
            H[j,i] = H[i,j]

    return H


def fiber_angle(V, insz, returncorr=False, **kw):
    """
    Estimate the angle of fibers through a volume.
    First takes the autocorrelation of the volume, then takes the eigenvalues and eigenvectors of the Hessian
    matrix at the center of the autocorrelation volume.  The eigenvector corresponding to the smallest magnitude
    (negative) eigenvalue is the dominant angle

    Extra keyword parameters are passed on to the ndautocorr function.

    :param V: Volume.  2 or 3D array
    :param insz: Size of the smaller region to use for the autocorrelation
    :param returncorr: True to return the autocorrelation matrix
    :return: vr, w, C: Eigenvectors and eigenvalues.  Sorted from smallest magnitude to largest.  C is the
    autocorrelation matrix, if returncorr is True
    """

    # get the autocorrelation
    C = ndautocorr(V, insz, **kw)
    ctr = (np.array(V.shape) - np.array(insz)+1)//2

    # and the hessian in the center
    H = hessian(C, ctr=ctr)

    # then eigenvalues and eigenvectors of the Hessian matrix
    w,vr = linalg.eig(H)

    # look for the smallest magnitude eigenvalue
    order = argsort(np.abs(w))
    vr = vr[:, order]
    w = w[order]

    if returncorr:
        return vr, w, C
    else:
        return vr, w