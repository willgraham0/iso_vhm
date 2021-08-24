from abc import ABC, abstractmethod

import numpy as np

from table import FvQvt, MAlpha


class Envelope(ABC):
    """Abstract base class for generating VHM envelopes."""

    def __init__(self, qv, qh, qm):
        self.qv = qv
        self.qh = qh
        self.qm = qm

    @abstractmethod
    def get_fh(self, fv, fm):
        pass

    @abstractmethod
    def get_fm(self, fv, fh):
        pass

    def _get_fh(self, qm, fv, fm, a=0):
        return np.nan_to_num(
            self.qh * np.sqrt(
                4 * a * (fv / self.qv) * (1 - fv / self.qv)
                + 16 * (1 - a) * (fv / self.qv) ** 2
                * (1 - fv / self.qv) ** 2 - (fm / qm) ** 2
            )
        )

    def _get_fm(self, qm, fv, fh, a=0):
        return np.nan_to_num(
            qm * np.sqrt(
                4 * a * (fv / self.qv) * (1 - fv / self.qv)
                + 16 * (1 - a) * (fv / self.qv) ** 2
                * (1 - fv / self.qv) ** 2 - (fh / self.qh) ** 2
            )
        )


class Sand(Envelope):
    """Concrete class for generating VHM envelope for Sand."""

    def __init__(self, qv, qh, qm, bm_b):
        super().__init__(qv, qh, qm)
        self.bm_b = bm_b

    def get_fh(self, fv, fm):
        """Return Fh values for the given Fv and Fm values."""
        qmp = self._get_qmp(fv)
        fh_unmodified = self._get_fh(self.qm, fv, fm)
        fh_modified = self._get_fh(qmp, fv, fm)
        return np.maximum(fh_unmodified, fh_modified)

    def get_fm(self, fv, fh):
        """Return Fm values for the given Fv and Fh values."""
        qmp = self._get_qmp(fv)
        fm_unmodified = self._get_fm(self.qm, fv, fh)
        fm_modified = self._get_fm(qmp, fv, fh)
        return np.maximum(fm_unmodified, fm_modified)

    def _get_qmp(self, fv):
        try:
            length = len(fv)
        except TypeError:
            length = 1
        qmps = self.qm * self.bm_b ** 3 * np.ones(length)
        qmpv = 0.15 * ((self.qm * 0.12) / (0.075 * self.qh)) * fv
        return np.minimum(qmps, qmpv)


class Clay(Envelope):
    """Concrete class for generating VHM envelope for Clay."""

    def __init__(self, qv, qh, qm, a, alpha, suction=True):
        super().__init__(qv, qh, qm)
        self.a = a
        self.alpha = alpha
        self.suction = suction
        self.fv_qvts = FvQvt()
        self.m_alphas = MAlpha()
        self.adhesive = self.a > 0.5

    def get_fh(self, fv, fm):
        """Return Fh values for the given Fv and Fm values."""
        if not self.adhesive:
            return self._get_fh(self.qm, fv, fm, a=self.a)
        else:
            fv_qvt = self.fv_qvts.lookup(self.a, self.alpha)
            m_alpha = self.m_alphas.lookup(self.a, self.alpha)

            if self.suction:
                pass
            else:
                pass

    def get_fm(self, fv, fh):
        """Return Fm values for the given Fv and Fh values."""
        if not self.adhesive:
            return self._get_fm(self.qm, fv, fh, a=self.a)
        else:
            fv_qvt = self.fv_qvts.lookup(self.a, self.alpha)
            m_alpha = self.m_alphas.lookup(self.a, self.alpha)

            if self.suction:
                pass
            else:
                pass
