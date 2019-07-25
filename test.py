import unittest
import numpy as np
from VelocityConversion import MantleConversion

def assemblage():
    a = {
        "ol": 0.617,
        "cpx": 0.133,
        "opx": 0.052,
        "gnt": 0.153,
        "jd": 0.045,
        "XFe": 0.11
    }
    return a

class TestVelocityConversion(unittest.TestCase):

    def test_vs_AlphaConst(self):
        MC = MantleConversion()
        MC.LoadArray(np.array([[0, 0, -50e3, 4287.65]]))
        MC.SetVelType('S')
        MC.SetMineralogy(assemblage())
        MC.FillTables()
        MC.CalcPT()
        self.assertAlmostEqual(MC.Result_T[0], 1301.05, 1)
        self.assertAlmostEqual(MC.Result_Rho[0], 3335.2, 1)

    def test_vs_AlphaT(self):
        MC = MantleConversion()
        MC.LoadArray(np.array([[0, 0, -50e3, 4287.65]]))
        MC.SetAlpha('T')
        MC.SetVelType('S')
        MC.SetMineralogy(assemblage())
        MC.FillTables()
        MC.CalcPT()
        self.assertAlmostEqual(MC.Result_T[0], 1321.95, 1)
        self.assertAlmostEqual(MC.Result_Rho[0], 3280.2, 1)

    def test_vs_AlphaPT(self):
        MC = MantleConversion()
        MC.LoadArray(np.array([[0, 0, -50e3, 4287.65]]))
        MC.SetAlpha('PT')
        MC.SetVelType('S')
        MC.SetMineralogy(assemblage())
        MC.FillTables()
        MC.CalcPT()
        self.assertAlmostEqual(MC.Result_T[0], 1319.45, 1)
        self.assertAlmostEqual(MC.Result_Rho[0], 3287.0, 1)

    def test_vp_AlphaConst(self):
        MC = MantleConversion()
        MC.LoadArray(np.array([[0, 0, -50e3, 8e3]]))
        MC.SetVelType('P')
        MC.SetMineralogy(assemblage())
        MC.FillTables()
        MC.CalcPT()
        self.assertAlmostEqual(MC.Result_T[0], 1141.9, 1)
        self.assertAlmostEqual(MC.Result_Rho[0], 3347.7, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2, buffer=True)
