"""
Created:  2018-08-08
Modified: 2023-07-11
Author:   Christopher Albert <albert@tugraz.at>
"""
#%%
from numpy import sin, cos, zeros_like, ones_like, array, sqrt, arctan2

B0 = 1.0  # magnetic field modulus normalization
iota0 = 1.0  # constant part of rotational transform
a = 0.5  # (equivalent) minor radius
R0 = 1.0  # (equivalent) major radius
Z0 = 0.0  # (equivalent) vertical position


class field:
    """Model tokamak field with B ~ B0*(1-r/R0*cos(th))"""

    def evaluate(self, r, th, ph):
        cth = cos(th)
        sth = sin(th)
        zer = zeros_like(r)

        self.Ar = zer
        self.hr = zer

        self.Ath = B0 * (r**2 / 2.0 - r**3 / (3.0 * R0) * cth)
        self.dAth = array(
            (B0 * (r - r**2 / R0 * cth), B0 * r**3 * sth / (3.0 * R0), zer)
        )
        self.d2Ath = array(
            (
                B0 * (1.0 - 2.0 * r / R0 * cth),
                B0 * r**3 * cth / (3.0 * R0),
                zer,
                B0 * r**2 / R0 * sth,
                zer,
                zer,
            )
        )

        self.Aph = -B0 * iota0 * (r**2 / 2.0 - r**4 / (4.0 * a**2))
        self.dAph = array((-B0 * iota0 * (r - r**3 / a**2), zer, zer))
        self.d2Aph = array(
            (
                -B0 * iota0 * (1.0 - 3.0 * r**2 / a**2),
                zer,
                zer,
                zer,
                zer,
                zer,
            )
        )

        self.hth = iota0 * (1.0 - r**2 / a**2) * r**2 / R0
        self.dhth = array(
            (
                (2.0 * iota0 * r * (a**2 - 2.0 * r**2)) / (a**2 * R0),
                zer,
                zer,
            )
        )
        self.d2hth = array(
            (
                (2.0 * iota0 * (a**2 - 6.0 * r**2)) / (a**2 * R0),
                zer,
                zer,
                zer,
                zer,
                zer,
            )
        )

        self.hph = R0 + r * cth
        self.dhph = array((cth, -r * sth, zer))
        self.d2hph = array((zer, -r * cth, zer, -sth, zer))

        self.B = B0 * (1.0 - r / R0 * cth)
        self.dB = array((-B0 / R0 * cth, B0 * r / R0 * sth, zer))
        self.d2B = array((zer, B0 * r / R0 * cth, zer, B0 / R0 * sth, zer, zer))

        self.Phie = zer
        self.dPhie = array((zer, zer, zer))
        self.d2Phie = array((zer, zer, zer, zer, zer, zer))


def jac_tor_cyl(R, PH, Z):
    r = sqrt((R - R0) ** 2 + (Z - Z0) ** 2)
    return array(
        [
            [(R - R0) / r, zeros_like(R), (Z - Z0) / r],
            [-(Z - Z0) / r**2, zeros_like(R), (R - R0) / r**2],
            [zeros_like(R), ones_like(R), zeros_like(R)],
        ]
    )


class field_cyl:
    def __init__(self):
        self.field = field()

    def evaluate(self, R, PH, Z):
        r = sqrt((R - R0) ** 2 + (Z - Z0) ** 2)
        th = arctan2(Z - Z0, R - R0)

        self.field.evaluate(r, th, PH)

        J = jac_tor_cyl(R, PH, Z)

        self.AR = (
            self.field.Ar * J[0, 0]
            + self.field.Ath * J[1, 0]
            + self.field.Aph * J[2, 0]
        )

        self.APH = (
            self.field.Ar * J[0, 1]
            + self.field.Ath * J[1, 1]
            + self.field.Aph * J[2, 1]
        )

        self.AZ = (
            self.field.Ar * J[0, 2]
            + self.field.Ath * J[1, 2]
            + self.field.Aph * J[2, 2]
        )

        self.hR = (
            self.field.hr * J[0, 0]
            + self.field.hth * J[1, 0]
            + self.field.hph * J[2, 0]
        )

        self.hPH = (
            self.field.hr * J[0, 1]
            + self.field.hth * J[1, 1]
            + self.field.hph * J[2, 1]
        )

        self.hZ = (
            self.field.hr * J[0, 2]
            + self.field.hth * J[1, 2]
            + self.field.hph * J[2, 2]
        )

f = field_cyl()

f.evaluate(R0 + 0.1, 0.0, Z0)
print([f.AR, f.APH, f.AZ])
print([f.hR, f.hPH, f.hZ])

R = R0 + 0.1
Z = Z0
PHstart = 0.0

def dPHdx1(x):
    f.evaluate(R, PH, Z)
    return

# %%
