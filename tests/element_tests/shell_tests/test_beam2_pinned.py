import unittest

import numpy as np

from tacs import TACS, constitutive, elements


class ElementTest(unittest.TestCase):
    def setUp(self):
        num_nodes = 2
        vars_per_node = 6
        num_vars = num_nodes * vars_per_node

        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-10
        else:
            self.dh = 1e-5
            self.rtol = 1e-2
        self.dtype = TACS.dtype

        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)
        self.print_level = 0
        self.elem_index = 0
        self.time = 0.0

        self.xpts = np.array([0.0, 0.0, 0.0, 1.7, 0.0, 0.0], dtype=self.dtype)
        np.random.seed(30)
        self.vars = np.random.rand(num_vars).astype(self.dtype)
        self.dvars = self.vars.copy()
        self.ddvars = self.vars.copy()

        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        self.props = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            alpha=cte,
            kappa=kappa,
        )

        ref_axis = np.array([0.0, 1.0, 1.0], dtype=self.dtype)
        self.transforms = [elements.BeamRefAxisTransform(ref_axis)]

        self.con_objects = [
            constitutive.IsoTubeBeamConstitutive(
                self.props, t=0.1, tNum=1, d=1.0, dNum=0
            ),
            constitutive.CompositeTubeBeamConstitutive(
                E11=240.0e3,
                E22=20.0e3,
                G12=15.0e3,
                nu12=0.25,
                rho=1600.0,
                X_c=900.0,
                X_t=1200.0,
                layup_angles=[0.0, 45.0, -45.0, 0.0],
                d=1.0,
                tw=0.1,
                dNum=0,
                twNum=1,
            ),
        ]

        self.matrix_types = [
            TACS.STIFFNESS_MATRIX,
            TACS.MASS_MATRIX,
            TACS.GEOMETRIC_STIFFNESS_MATRIX,
        ]

        elements.SeedRandomGenerator(0)

    def _make_element(self, transform, con):
        return elements.Beam2Pinned(transform, con)

    def test_element_residual(self):
        dh = 1e-5
        rtol = 1e-2
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        fail = elements.TestElementResidual(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            dh,
                            self.print_level,
                            self.atol,
                            rtol,
                        )
                        self.assertFalse(fail)

    def test_element_jacobian(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        fail = elements.TestElementJacobian(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            -1,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_adj_res_product(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        dvs = element.getDesignVars(self.elem_index)
                        fail = elements.TestAdjResProduct(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            dvs,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_adj_res_xpt_product(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        fail = elements.TestAdjResXptProduct(
                            element,
                            self.elem_index,
                            self.time,
                            self.xpts,
                            self.vars,
                            self.dvars,
                            self.ddvars,
                            self.dh,
                            self.print_level,
                            self.atol,
                            self.rtol,
                        )
                        self.assertFalse(fail)

    def test_element_mat_dv_sens(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        dvs = element.getDesignVars(self.elem_index)
                        for matrix_type in self.matrix_types:
                            with self.subTest(matrix_type=matrix_type):
                                fail = elements.TestElementMatDVSens(
                                    element,
                                    matrix_type,
                                    self.elem_index,
                                    self.time,
                                    self.xpts,
                                    self.vars,
                                    dvs,
                                    self.dh,
                                    self.print_level,
                                    self.atol,
                                    self.rtol,
                                )
                                self.assertFalse(fail)

    def test_element_mat_xpt_sens(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        for matrix_type in self.matrix_types:
                            with self.subTest(matrix_type=matrix_type):
                                fail = elements.TestElementMatXptSens(
                                    element,
                                    matrix_type,
                                    self.elem_index,
                                    self.time,
                                    self.xpts,
                                    self.vars,
                                    self.dh,
                                    self.print_level,
                                    self.atol,
                                    self.rtol,
                                )
                                self.assertFalse(fail)

    def test_element_mat_sv_sens(self):
        for transform in self.transforms:
            with self.subTest(transform=transform):
                for con in self.con_objects:
                    with self.subTest(con=con.getObjectName()):
                        element = self._make_element(transform, con)
                        for matrix_type in self.matrix_types:
                            with self.subTest(matrix_type=matrix_type):
                                fail = elements.TestElementMatSVSens(
                                    element,
                                    matrix_type,
                                    self.elem_index,
                                    self.time,
                                    self.xpts,
                                    self.vars,
                                    self.dh,
                                    self.print_level,
                                    self.atol,
                                    self.rtol,
                                )
                                self.assertFalse(fail)


if __name__ == "__main__":
    unittest.main()
