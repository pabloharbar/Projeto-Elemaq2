class Material:
    def __init__(
        self,
        elasticity_module: float,
        poisson_coef: float,
        yield_stress: float,
        ultimate_stress: float,
        bending_stress_strength: float,
        contact_stress_strength: float,
    ) -> None:
        self.elasticity_module = elasticity_module
        self.poisson_coef = poisson_coef
        self.stiffness_module = self.elasticity_module / (2 * (1 + self.poisson_coef))
        self.yield_stress = yield_stress
        self.ultimate_stress = ultimate_stress
        self.bending_stress_strength = bending_stress_strength
        self.contact_stress_strength = contact_stress_strength
