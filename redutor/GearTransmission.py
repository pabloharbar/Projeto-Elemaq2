import math

from redutor import Gear


class GearTransmission:
    def __init__(
        self, gear1: Gear, gear2: Gear, seconds_of_use: float, position: float
    ) -> None:
        self.gear1 = gear1
        self.gear2 = gear2
        self.seconds_of_use = seconds_of_use
        self.position = position

    def calculate_forces(
        self, input_power: float, input_velocity: float, roller_efficiency: float
    ):
        self.w1 = input_velocity * math.pi / 30
        self.p1 = input_power
        self.T1 = input_power / self.w1

        Z1 = self.gear1.number_of_teeths
        Z2 = self.gear2.number_of_teeths
        self.gearing_efficiency = 1 - 0.5 * (Z2 - Z1) / (Z1 * Z2)

        self.w2 = self.w1 * self.gear1.primitive_diam / self.gear2.primitive_diam
        print(self.gear1.primitive_diam, self.gear2.primitive_diam)
        self.p2 = input_power * roller_efficiency * self.gearing_efficiency
        self.T2 = self.p2 / self.w2

        self.Ft1 = self.T1 * 2 / self.gear1.primitive_diam
        self.Ft2 = self.T2 * 2 / self.gear2.primitive_diam
        self.Fn1 = self.Ft1 / math.cos(self.gear1.pressure_angle)
        self.Fn2 = self.Ft2 / math.cos(self.gear2.pressure_angle)
        self.Fr1 = self.Fn1 * math.sin(self.gear1.pressure_angle)
        self.Fr2 = self.Fn2 * math.sin(self.gear2.pressure_angle)

    def _calculate_bending_stress(self):
        Kb = 1
        Ka = 1
        Km = 1.6
        Ks = 1
        Ki = 1
        Qv = 7
        B = (12 - Qv) ** (2 / 3) / 4
        A = 50 + 56 * (1 - B)
        Vt = self.w1 * self.gear1.primitive_diam / 2
        Kv = A / (A + (0.2 * Vt) ** 0.5)
        K_total = Ka * Km * Ks * Kb * Ki
        aux1 = self.gear1.thickness * self.gear1.modulo * self.gear1.J * Kv
        aux2 = self.gear2.thickness * self.gear2.modulo * self.gear2.J * Kv
        self.sigma_b1 = self.Ft1 * K_total / aux1
        self.sigma_b2 = self.Ft2 * K_total / aux2

    def _calculate_contact_stress(self):
        Kf = 1
        Ka = 1
        Km = 1.6
        Ks = 1
        Qv = 7
        B = (12 - Qv) ** (2 / 3) / 4
        A = 50 + 56 * (1 - B)
        Vt = self.w1 * self.gear1.primitive_diam / 2
        Kv = A / (A + (0.2 * Vt) ** 0.5)

        aux_Cp = 2 * math.pi * (1 - self.gear1.material.poisson_coef**2)
        Cp = (1 / (aux_Cp / self.gear1.material.elasticity_module)) ** 0.5
        rp = self.gear1.primitive_diam / 2
        rho_p = (
            (rp + self.gear1.modulo) ** 2
            - (rp * math.cos(self.gear1.pressure_angle)) ** 2
        ) ** 0.5 - math.pi * self.gear1.modulo * math.cos(self.gear1.pressure_angle)
        rho_g = (
            math.sin(self.gear1.pressure_angle)
            * (self.gear1.primitive_diam + self.gear2.primitive_diam)
            / 2
            - rho_p
        )
        I = math.cos(self.gear1.pressure_angle) / (
            (1 / rho_p + 1 / rho_g) * self.gear1.primitive_diam
        )
        self.sigma_c1 = (
            Cp
            * (
                self.Ft1
                * Ka
                * Km
                * Ks
                * Kf
                / (self.gear1.thickness * I * self.gear1.primitive_diam * Kv)
            )
            ** 0.5
        )
        self.sigma_c2 = (
            Cp
            * (
                self.Ft2
                * Ka
                * Km
                * Ks
                * Kf
                / (self.gear2.thickness * I * self.gear2.primitive_diam * Kv)
            )
            ** 0.5
        )

    def _calculate_bending_fatigue(self):
        Kt = 1
        Kr = 1
        N1 = self.w1 * self.seconds_of_use / (2 * math.pi)
        N2 = self.w2 * self.seconds_of_use / (2 * math.pi)
        Kl1 = 1.3558 * N1 ** (-0.0178)
        Kl2 = 1.3558 * N2 ** (-0.0178)
        Sfb1 = Kl1 * self.gear1.material.bending_stress_strength / (Kt * Kr)
        Sfb2 = Kl2 * self.gear2.material.bending_stress_strength / (Kt * Kr)
        self.CSb1 = Sfb1 / self.sigma_b1
        self.CSb2 = Sfb2 / self.sigma_b2

    def _calculate_contact_fatigue(self):
        Ct = 1
        Cr = 1
        Ch = 1
        N1 = self.w1 * self.seconds_of_use / (2 * math.pi)
        N2 = self.w2 * self.seconds_of_use / (2 * math.pi)
        Cl1 = 1.4488 * N1 ** (-0.023)
        Cl2 = 1.4488 * N2 ** (-0.023)

        Sfc1 = Cl1 * self.gear1.material.contact_stress_strength * Ch / (Ct * Cr)
        Sfc2 = Cl2 * self.gear2.material.contact_stress_strength * Ch / (Ct * Cr)

        self.CSc1 = (Sfc1 / self.sigma_c1) ** 2
        self.CSc2 = (Sfc2 / self.sigma_c2) ** 2

    def _print_report(self):
        print(f"Tensão de flexão no pinhao: {self.sigma_b1}")
        print(f"Tensão de flexão no coroa: {self.sigma_b2}")
        print(f"Tensão de contato no pinhao: {self.sigma_c1}")
        print(f"Tensão de contato no coroa: {self.sigma_c2}")
        print(f"Coeficiente de segurança de contato pinhao: {self.CSc1}")
        print(f"Coeficiente de segurança de contato coroa: {self.CSc2}")
        print(f"Coeficiente de segurança de flexao pinhao: {self.CSb1}")
        print(f"Coeficiente de segurança de flexao coroa: {self.CSb2}")

    def calculate_stress(self):
        self._calculate_bending_stress()
        self._calculate_contact_stress()
        self._calculate_bending_fatigue()
        self._calculate_contact_fatigue()
        self._print_report()
