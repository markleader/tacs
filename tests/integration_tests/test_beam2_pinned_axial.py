import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements, functions


class Beam2PinnedAxialTest(unittest.TestCase):
    def setUp(self):
        self.comm = MPI.COMM_WORLD
        self.dtype = TACS.dtype

        self.length = 2.0
        self.force = 125.0
        self.inner_diameter = 0.08
        self.wall_thickness = 0.01
        self.rho = 2700.0
        self.E = 70.0e9
        self.nu = 0.3
        self.ys = 270.0e6

        self.area = 0.25 * np.pi * (
            (self.inner_diameter + self.wall_thickness) ** 2
            - self.inner_diameter**2
        )
        self.axial_stiffness = self.E * self.area / self.length
        self.mass = self.rho * self.area * self.length

        self.ref_axis = np.array([0.0, 1.0, 0.0], dtype=self.dtype)

    def _make_element(self):
        props = constitutive.MaterialProperties(
            rho=self.rho, E=self.E, nu=self.nu, ys=self.ys
        )
        con = constitutive.IsoTubeBeamConstitutive(
            props,
            d=self.inner_diameter,
            dNum=0,
            t=self.wall_thickness,
            tNum=1,
            Kb=1.0,
            Lb=self.length,
        )
        transform = elements.BeamRefAxisTransform(self.ref_axis)
        return elements.Beam2Pinned(transform, con)

    def _create_assembler(self, constrain_for_axial_solve):
        creator = TACS.Creator(self.comm, 6)

        if self.comm.rank == 0:
            ptr = np.array([0, 2], dtype=np.intc)
            conn = np.array([0, 1], dtype=np.intc)
            comp_ids = np.array([0], dtype=np.intc)
            creator.setGlobalConnectivity(2, ptr, conn, comp_ids)

            if constrain_for_axial_solve:
                bcnodes = np.array([0, 1], dtype=np.intc)
                bcvars = np.array(
                    [0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5], dtype=np.intc
                )
                bcptr = np.array([0, 6, 11], dtype=np.intc)
                creator.setBoundaryConditions(bcnodes, bcvars, bcptr)

            xpts = np.array(
                [0.0, 0.0, 0.0, self.length, 0.0, 0.0], dtype=self.dtype
            )
            creator.setNodes(xpts)

        creator.setElements([self._make_element()])
        return creator.createTACS()

    def _dense_matrix(self, assembler, alpha, gamma):
        assembler.zeroVariables()
        res = assembler.createVec()
        mat = assembler.createSchurMat()
        assembler.assembleJacobian(alpha, 0.0, gamma, res, mat)
        return mat.getDenseMatrix()

    @unittest.skipIf(MPI.COMM_WORLD.size != 1, "Dense matrix inspection is serial-only")
    def test_stiffness_and_mass_structure(self):
        assembler = self._create_assembler(constrain_for_axial_solve=False)

        K = self._dense_matrix(assembler, alpha=1.0, gamma=0.0)
        M = self._dense_matrix(assembler, alpha=0.0, gamma=1.0)

        rot = np.array([3, 4, 5, 9, 10, 11], dtype=np.intc)
        transverse = np.array([1, 2, 7, 8], dtype=np.intc)

        np.testing.assert_allclose(K[np.ix_(rot, np.arange(12))], 0.0, atol=1e-12)
        np.testing.assert_allclose(K[np.ix_(np.arange(12), rot)], 0.0, atol=1e-12)
        np.testing.assert_allclose(
            K[np.ix_(transverse, np.arange(12))], 0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            K[np.ix_(np.arange(12), transverse)], 0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            K[np.ix_([0, 6], [0, 6])],
            self.axial_stiffness * np.array([[1.0, -1.0], [-1.0, 1.0]]),
            rtol=1e-12,
            atol=1e-12,
        )

        for dof_pair in ([0, 6], [1, 7], [2, 8]):
            np.testing.assert_allclose(
                M[np.ix_(dof_pair, dof_pair)],
                (self.mass / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]]),
                rtol=1e-12,
                atol=1e-12,
            )
        np.testing.assert_allclose(M[np.ix_(rot, np.arange(12))], 0.0, atol=1e-12)
        np.testing.assert_allclose(M[np.ix_(np.arange(12), rot)], 0.0, atol=1e-12)

    def test_exact_axial_response_and_function_smoke(self):
        assembler = self._create_assembler(constrain_for_axial_solve=True)

        res = assembler.createVec()
        ans = assembler.createVec()
        mat = assembler.createSchurMat()
        force = assembler.createVec()

        assembler.zeroVariables()
        assembler.assembleJacobian(1.0, 0.0, 0.0, res, mat)

        force_arr = force.getArray()
        force_arr[:] = 0.0
        force_arr[6] = self.force

        assembler.applyBCs(force)
        res.axpy(-1.0, force)

        pc = TACS.Pc(mat)
        pc.factor()
        gmres = TACS.KSM(mat, pc, 20, 2)
        gmres.setTolerances(1e-12, 1e-30)
        gmres.solve(res, ans)
        ans.scale(-1.0)

        state = ans.getArray().copy()
        np.testing.assert_allclose(
            state[6], self.force / self.axial_stiffness, rtol=1e-12, atol=1e-12
        )

        assembler.setVariables(ans)
        func_list = [
            functions.StructuralMass(assembler),
            functions.KSFailure(assembler, ksWeight=10.0),
        ]
        func_vals = assembler.evalFunctions(func_list)

        np.testing.assert_allclose(func_vals[0], self.mass, rtol=1e-12, atol=1e-12)
        self.assertTrue(np.isfinite(np.real(func_vals[1])))


if __name__ == "__main__":
    unittest.main()
