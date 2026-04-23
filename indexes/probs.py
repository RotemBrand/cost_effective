r"""
Define probability classes to make a single api for working with
different probabilities type(Float\Poly\Array)
"""

from abc import ABC, abstractmethod
from copy import copy
import numpy as np
from numpy.polynomial.polynomial import polydiv, polyval, polypow, polyroots
import scipy.special as sp

def polysub(p1, p2):
    if p1.size == p2.size:
        return p1 - p2
    res = np.zeros(max([p1.size, p2.size]))
    res[:p1.size] += p1
    res[:p2.size] -= p2
    return res

def polyadd(p1, p2):
    if p1.size == p2.size:
        return p1 + p2
    res = np.zeros(max([p1.size, p2.size]))
    res[:p1.size] += p1
    res[:p2.size] += p2
    return res

def polymul(p1, p2):
    return np.convolve(p1, p2)

class Prob(ABC):
    """
    The abstract class used to define the api of a probability.
    Contain basic actions such as adding, multiply etc...
    """

    def __init__(self, prob):
        self.prob = prob

    @abstractmethod
    def __add__(self, other):
        pass

    @abstractmethod
    def __mul__(self, other):
        pass

    @abstractmethod
    def __sub__(self, other):
        pass

    @abstractmethod
    def __truediv__(self, other):
        pass

    def __floordiv__(self, other):
        return self.__truediv__(other)

    def __radd__(self, other):
        return self.__add__(other)
    
    def __rsub__(self, other):
        return self.__sub__(other)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    

    @classmethod
    def fast_prod(cls, fast_obj1, fast_obj2):
        """
        The function allows to make multiplications of the self.prob objects without create an
        instance of Prob to speed the calculation.
        """
        return fast_obj1 * fast_obj2

    @classmethod
    def is_zero(cls, prob) -> bool:
        """
        Return true if the probability is zero
        """
        if isinstance(prob, Prob):
            return np.all(prob.prob == 0)
        return np.all(prob == 0)

    @abstractmethod
    def __str__(self):
        pass

    def __repr__(self):
        return self.__str__()

    def is_close(self, other, tol: float) -> bool:
        return np.max(np.abs((self - other).prob)) < tol

    def __neg__(self) -> "Prob":
        return self.__class__(0) - self

    def __eq__(self, other) -> bool:
        return np.all(self.prob == other.prob)
    
    def __lt__(self, other) -> bool:
        return self._compare(other, op="lt")

    def __le__(self, other) -> bool:
        return self._compare(other, op="le")

    def __gt__(self, other) -> bool:
        return self._compare(other, op="gt")

    def __ge__(self, other) -> bool:
        return self._compare(other, op="ge")

    def _compare(self, other, op: str) -> bool:
        if isinstance(self.prob, np.ndarray) and isinstance(other.prob, np.ndarray):
            # Pad arrays with zeros to match size
            max_len = max(len(self.prob), len(other.prob))
            prob1 = np.pad(self.prob, (0, max_len - len(self.prob)), constant_values=0).tolist()
            prob2 = np.pad(other.prob, (0, max_len - len(other.prob)), constant_values=0).tolist()
        elif not isinstance(self.prob, np.ndarray) and not isinstance(other.prob, np.ndarray):
            # Direct comparison for non-array probabilities
            prob1 = self.prob
            prob2 = other.prob
        else:
            # Raise error for mismatched types
            raise TypeError("Both probabilities must be either arrays or non-arrays for comparison.")

        # Perform the comparison
        if op == "lt":
            return prob1 < prob2
        elif op == "le":
            return prob1 <= prob2
        elif op == "gt":
            return prob1 > prob2
        elif op == "ge":
            return prob1 >= prob2
        else:
            raise ValueError("Invalid comparison operator.")
    
    need_to_check_zero_events = False


class Float(Prob):
    """
    Represent a float
    """

    def __init__(self, prob):
        self.prob = float(prob)

    def __add__(self, other):
        return Float(self.prob + float(other))

    def __mul__(self, other):
        return Float(self.prob * float(other))

    def __sub__(self, other):
        return Float(self.prob - float(other))

    def __truediv__(self, other):
        return Float(self.prob / float(other))

    def __float__(self):
        return self.prob

    def __str__(self):
        return f"Float({self.prob})"

    def __call__(self, x):
        raise ValueError("__call__ not implemented for Float")


class Array(Prob):
    """
    Represent an array of floats
    """

    def __init__(self, prob):
        self.prob = np.array(prob)

    def __add__(self, other) -> "Array" :
        if isinstance(other, Prob):
            return Array(self.prob + other.prob)
        return Array(self.prob + np.array(other))

    def __mul__(self, other) -> "Array":
        if isinstance(other, Prob):
            return Array(self.prob * other.prob)
        return Array(self.prob * np.array(other))

    def __sub__(self, other) -> "Array":
        if isinstance(other, Prob):
            return Array(self.prob - other.prob)
        return Array(self.prob - np.array(other))

    def __truediv__(self, other) -> "Array":
        # TODO: what to do with zero div?
        if isinstance(other, Prob):
            return Array(self.prob / other.prob)
        return Array(self.prob / np.array(other))

    def __iter__(self):
        return iter(self.prob)

    def __str__(self) -> str:
        return f"Array({self.prob})"


class Poly(Prob):
    """
    Represent a polynomial by its coefficients array
    """

    def __init__(self, prob):
        if isinstance(prob, Poly):
            self.prob = np.trim_zeros(np.array(prob.prob), "b")
        elif not hasattr(prob, "__len__"):
            self.prob = np.array([prob])
        else:
            self.prob = np.trim_zeros(np.array(prob), "b")
        if len(self.prob) == 0:
            self.prob = np.array([0])
        if self.prob.ndim == 0:
            self.prob = self.prob.reshape(1)

    def __add__(self, other):
        if hasattr(other, "prob"):
            return self.__class__(polyadd(self.prob, other.prob))
        return self.__class__(polyadd(self.prob, other))

    def __mul__(self, other):
        if hasattr(other, "prob"):
            return self.__class__(polymul(self.prob, other.prob))
        return self.__class__(polymul(self.prob, other))

    def __sub__(self, other):
        if hasattr(other, "prob"):
            return self.__class__(polysub(self.prob, other.prob))
        return self.__class__(polysub(self.prob, other))

    def __pow__(self, exponent):
        return self.__class__(polypow(self.prob, exponent))

    def roots(self):
        return polyroots(self.prob)

    def __truediv__(self, other):
        # TODO: what to do with zero div?
        if isinstance(other, Prob):
            return self.__class__(polydiv(self.prob, other.prob)[0])
        return self.__class__(polydiv(self.prob, other)[0])

    def __getitem__(self, index):
        if 0 <= index < len(self.prob):
            return self.prob[index]
        else:
            return 0

    @classmethod
    def fast_prod(cls, fast_obj1, fast_obj2):
        return np.convolve(fast_obj1, fast_obj2)

    def __str__(self):
        return f"Poly({self.prob})"

    def eval(self, x):
        """
        Evaluate the value of the polynomial at point x
        """
        return polyval(x, self.prob)

    def __call__(self, x):
        """
        Evaluate the value of the polynomial at point x
        """
        return self.eval(x)

    @property
    def degree(self):
        return len(self.prob) - 1

    def to_binompoly(self, degree=None) -> "BinomPoly":
        if degree is None:
            degree = self.degree
        return BinomPoly(BinomPoly.from_poly(self.prob, degree=degree), degree=degree)

def array_to_size(arr: np.array, s: int) -> np.array:
    """## Resize arr to length s. if len(arr) < s, return only
    the first s elements. Otherwise, pad with zeros

    ### Args:
        - `arr (np.array)`: array
        - `s (int)`: new size for the array

    ### Returns:
        - `np.array`: New array from length s
    """
    if not hasattr(arr, "__len__"):
        input = np.array([arr])
    else:
        input = np.array(arr)
    if s <= len(input):
        return input[:s]
    res = np.zeros(s)
    res[: len(input)] = input
    return res


def create_poly_class(degree: int) -> ABC:
    """Create a class of polynomials with fixed degree

    Parameters
    ----------
    degree : int
        The maximal degree of the polynomials

    Returns
    -------
    ABC
        A class that inherited from Poly but has a fixed degree

    Raises
    ------
    TypeError
        degree is a non int or negative
    """

    # type error
    if not isinstance(degree, int) or (degree <= 0):
        raise TypeError("The type of degree should be int > 0")

    # create class
    class FixedDegreePolynomProb(Poly):
        def __init__(self, prob):
            if isinstance(prob, Poly):
                input_prob = prob.prob
            elif not hasattr(prob, "__len__"):
                input_prob = np.array([prob])
            else:
                input_prob = np.array(prob)
            self.prob = array_to_size(input_prob, degree + 1)
            if self.prob.ndim == 0:
                self.prob = self.prob.reshape(1)

        def __add__(self, other):
            if isinstance(other, Poly):
                return FixedDegreePolynomProb(self.prob + other.prob)
            return FixedDegreePolynomProb(polyadd(self.prob, other))

        def __mul__(self, other):
            if isinstance(other, Poly):
                return FixedDegreePolynomProb(polymul(self.prob, other.prob))
            return FixedDegreePolynomProb(polymul(self.prob, other))

        def __sub__(self, other):
            if isinstance(other, Poly):
                return FixedDegreePolynomProb(polysub(self.prob, other.prob))
            return FixedDegreePolynomProb(polysub(self.prob, other))

        def __truediv__(self, other):
            if isinstance(other, Poly):
                return FixedDegreePolynomProb(polydiv(self.prob, other.prob)[0])
            return FixedDegreePolynomProb(polydiv(self.prob, other)[0])

        def __floordiv__(self, other):
            if isinstance(other, Poly):
                return FixedDegreePolynomProb(polydiv(self.prob, other.prob)[0])
            return FixedDegreePolynomProb(polydiv(self.prob, other)[0])

        @classmethod
        def fast_prod(cls, fast_obj1, fast_obj2):
            return np.convolve(fast_obj1, fast_obj2)[: degree + 1]

        def to_fast_object(self):
            return copy(self.prob)

        def __str__(self):
            trim_poly = np.trim_zeros(self.prob[: degree + 1], "b")
            if len(trim_poly) == 0:
                trim_poly = np.array([0])
            return f"PolyDeg{degree}({trim_poly})"

    # Dynamically set the class name
    FixedDegreePolynomProb.__name__ = f"PolyDeg{degree}"
    FixedDegreePolynomProb.need_to_check_zero_events = True
    return FixedDegreePolynomProb


#### tries ######


class PolySlow(Poly):
    def to_fast_object(self):
        return PolySlow(self.prob)

    def __str__(self):
        return f"PolySlow({self.prob})"


class BinomPoly(Prob):
    def __init__(self, prob, degree: int):
        if degree < 0:
            raise TypeError("degree should be non-negative int")
        if isinstance(prob, Poly):
            self.prob = BinomPoly.from_poly(prob.prob[: degree + 1], degree=degree)
        elif isinstance(prob, BinomPoly):
            if prob.degree == self.degree:
                self.prob = copy(prob.prob)
            else:
                self.prob = BinomPoly.from_poly(
                    prob.to_poly()[: degree + 1], degree=degree
                )
        else:
            if hasattr(prob, "__len__"):
                self.prob = np.array(prob)[: degree + 1]
            else:
                self.prob = np.array([prob])
        self._degree = degree

    @classmethod
    def from_poly(cls, poly, degree=None) -> np.array:
        if isinstance(poly, Poly):
            prob = poly.prob
        else:
            prob = np.array(poly)
        m = len(prob) - 1 if degree is None else degree
        res_p = np.zeros(m + 1)
        for i in range(0, m + 1):
            res_p[i] = prob[i] if i < len(prob) else 0
            for j in range(0, i):
                res_p[i] -= (-1) ** (i - j) * sp.comb(m - j, i - j) * res_p[j]
        return res_p

    def to_poly(self) -> Poly:
        m = self.degree
        res = [0] * (m + 1)
        for i in range(m + 1):
            for j in range(i + 1):
                res[i] += (-1) ** (i - j) * sp.comb(m - j, i - j) * self.prob[j]
        return Poly(res)

    def __getitem__(self, index):
        if 0 <= index < len(self.prob):
            return self.prob[index]
        else:
            return 0

    def __add__(self, other):
        raise ValueError(
            "__add__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    def __mul__(self, other):
        raise ValueError(
            "__mul__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    def __sub__(self, other):
        raise ValueError(
            "__sub__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    def __truediv__(self, other):
        raise ValueError(
            "__truediv__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    @classmethod
    def fast_prod(cls, fast_obj1, fast_obj2):
        raise ValueError(
            "__fast_prod__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    def to_fast_object(self):
        raise ValueError(
            "__to_fast_object__ not implemented for BinomPoly. Please use create_binompoly_class(degree) insted"
        )

    def __str__(self):
        return f"BinomPoly({self.prob},)"

    def eval(self, x):
        m = self.degree
        return np.sum(
            [c * (x**d) * ((1 - x) ** (m - d)) for d, c in enumerate(self.prob)], axis=0
        )

    def __call__(self, x):
        return self.eval(x)

    @property
    def degree(self):
        return self._degree


def create_binompoly_class(degree: int):
    # type error
    if not isinstance(degree, int) or (degree <= 0):
        raise TypeError("The type of degree should be int > 0")

    # create class
    class BinomPolyFixedDegree(BinomPoly):
        def __init__(self, prob):
            super().__init__(prob, degree=degree)
            self.prob = array_to_size(self.prob, degree + 1)

        def __add__(self, other):
            if isinstance(other, BinomPolyFixedDegree):
                return self.__class__(self.prob + other.prob)
            return self.__class__(polyadd(self.prob, other))

        def __mul__(self, other):
            if isinstance(other, BinomPolyFixedDegree):
                return self.__class__(np.convolve(self.prob, other.prob))
            return self.__class__(np.convolve(self.prob, other))

        def __sub__(self, other):
            if isinstance(other, BinomPolyFixedDegree):
                return self.__class__(self.prob - other.prob)
            return self.__class__(polysub(self.prob, other))

        def __truediv__(self, other):
            if isinstance(other, BinomPolyFixedDegree):
                return self.__class__(polydiv(self.prob, other.prob)[0])
            return self.__class__(polydiv(self.prob, other)[0])

        def __floordiv__(self, other):
            if isinstance(other, BinomPolyFixedDegree):
                return self.__class__(polydiv(self.prob, other.prob)[0])
            return self.__class__(polydiv(self.prob, other)[0])

        @classmethod
        def fast_prod(cls, fast_obj1, fast_obj2):
            return np.convolve(fast_obj1, fast_obj2)[: degree + 1]

        def to_fast_object(self):
            return copy(self.prob)

        def __str__(self):
            trimed_poly = np.trim_zeros(self.prob, "b")
            if len(trimed_poly) == 0:
                trimed_poly = np.array([0])
            if np.all(np.round(trimed_poly) == trimed_poly):
                trimed_poly = trimed_poly.astype(int)
            return f"BinomPolyDeg{degree}({trimed_poly})"

    # Dynamically set the class name
    BinomPolyFixedDegree.__name__ = f"BinomPolyDeg{degree}"
    BinomPolyFixedDegree.need_to_check_zero_events = True
    return BinomPolyFixedDegree
