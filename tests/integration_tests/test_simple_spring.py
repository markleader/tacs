import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
A 6 DOF spring element fixed at node 1 with a uniform force of 100 N applied in each direction. 

The stiffness values for the spring are given below:
kx = 1 N/m
ky = 2 N/m
kz = 3 N/m
krx = 4 N/m
kry = 5 N/m
krz = 6 N/m
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/simple_spring.bdf")

FUNC_REFS = {
    "uniform_force_mass": 0.0,
    "uniform_force_compliance": 24500.000000000007,
    "uniform_force_x_disp": 99.93068528194402,
    "uniform_force_y_disp": 49.930685281944015,
    "uniform_force_z_disp": 33.264018615277344,
}

# Force to apply
f = np.ones(6) * 100.0

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 1  # this is how many MPI processes to use for this TestCase.

    def setup_pytacs(self, comm, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Overwrite default check values
        if dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            # Set one thickness dv for every component
            k = np.arange(6) + 1
            con = constitutive.DOFSpringConstitutive(k=k)
            transform = elements.SpringRefAxisTransform([0.0, 1.0, 0.0])
            elem = elements.SpringElement(transform, con)
            return elem

        # Set up constitutive objects and elements
        fea_assembler.initialize(elemCallBack)

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # This problem has no dvs

        # Define perturbation array that moves all nodes on shell
        xpts = fea_assembler.getOrigNodes()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_assembler, problems):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        for problem in problems:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[1.0, 0.0, 0.0],
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 1.0, 0.0],
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 1.0],
            )
        func_list = ["mass", "compliance", "x_disp", "y_disp", "z_disp"]
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Create tacs transient problem
        problem = fea_assembler.createStaticProblem("uniform_force")
        problem.addLoadToNodes(1, f, nastranOrdering=False)
        return [problem]
