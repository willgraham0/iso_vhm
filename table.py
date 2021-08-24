from abc import ABC, abstractmethod


class Table(ABC):
    """Abstract base class for lookup tables for parameters a and alpha."""
    _a_values = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    _alpha_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    @property
    @abstractmethod
    def values(self):
        pass

    def _get_a_index(self, value: float) -> int:
        return self._a_values.index(value)

    def _get_alpha_index(self, value: float) -> int:
        return self._alpha_values.index(value)

    def lookup(self, a: float, alpha: float) -> float:
        """Return the table value for the given values of a and alpha."""
        a_index = self._get_a_index(a)
        alpha_index = self._get_alpha_index(alpha)
        return self.values[alpha_index][a_index]


class FvQvt(Table):
    """Concrete class for retrieving values of (Fv/Qv)t for various values of a and alpha."""
    values = (
        (0.354, 0.334, 0.308, 0.293, 0.276, 0.238, 0.200),
        (0.387, 0.373, 0.354, 0.343, 0.331, 0.300, 0.265),
        (0.418, 0.408, 0.396, 0.388, 0.379, 0.357, 0.329),
        (0.447, 0.441, 0.433, 0.428, 0.423, 0.409, 0.390),
        (0.474, 0.471, 0.468, 0.465, 0.463, 0.457, 0.448),
        (0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500),
    )


class MAlpha(Table):
    """Concrete class for retrieving values of m alpha for various values of a and alpha."""
    values = (
        (1.172, 1.200, 1.239, 1.264, 1.295, 1.378, 1.500),
        (0.902, 0.917, 0.937, 0.950, 0.965, 1.006, 1.067),
        (0.653, 0.661, 0.670, 0.676, 0.683, 0.701, 0.729),
        (0.422, 0.425, 0.429, 0.431, 0.434, 0.440, 0.450),
        (0.205, 0.206, 0.207, 0.207, 0.208, 0.209, 0.211),
        (0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000),
    )
