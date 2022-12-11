import math

from redutor import Pulley


class PulleyTransmission:
    def __init__(
        self, polia1: Pulley, polia2: Pulley, power: float, position: float
    ) -> None:
        self.P1 = polia1
        self.P2 = polia2
        self.power = power
        self.position = position

    def _calculate_geometry(self):
        # Dimensionamento da transmissao pela correia
        D = self.P2.primitive_diameter
        d = self.P1.primitive_diameter
        i = D / d
        Cc1 = D + d
        self.Ld1 = 2 * Cc1 + math.pi * (D + d) / 2 + (D - d) ** 2 / (4 * Cc1)
        self.comprimento_correia = 1.29  # Primeiro valor maior que Ld1 do catálogo SKF
        a = 2 * self.comprimento_correia - math.pi * (D + d)
        self.distancia_centros = (a + (a**2 - 8 * (D - d) ** 2) ** 0.5) / 8

    def _input_constants(self):
        self.theta = 161.4508026
        self.mi = 0.3
        self.phi = math.radians(34)
        self.Kb = 576 * 4.44822 * 0.0254  # lbf in para N m
        self.Kc = 0.965 * 4.44822 / 0.3048**2  # lbf.s2/ft2 para N.s2/m2
        self.K = 1193 * 4.44822  # lbf para N
        self.b = 10.926

    def _calculate_forces(self):
        V = self.P1.angular_velocity * self.P1.primitive_diameter / 2
        self.Fc = self.Kc * (V / 1000) ** 2
        e = math.exp(self.mi * self.theta / math.sin(self.phi / 2))
        T = self.power / self.P1.angular_velocity
        self.F2 = (self.Fc - 2 * T / self.P1.primitive_diameter - self.Fc * e) / (1 - e)
        self.F1 = self.F2 + 2 * T / self.P1.primitive_diameter

    def _calculate_durability(self):
        V = self.P1.angular_velocity * self.P1.primitive_diameter / 2
        T1 = self.F1 + self.Kb / self.P1.primitive_diameter
        T2 = self.F1 + self.Kb / self.P2.primitive_diameter
        Np = (self.K / T1) ** (-self.b) + (self.K / T2) ** (-self.b)  # ** (-1)
        Np = Np if Np < 10e9 else 10e9
        self.t = Np * (math.pi * self.P1.primitive_diameter) / (720 * V)

    def _print_report(self):
        print(f"Polia\nDistancia entre centros calculada: {self.distancia_centros}")
        print(f"Comprimento da correia: {self.Ld1}->{self.comprimento_correia}")
        print(f"Forças calculadas\nF1: {self.F1}\nF2: {self.F2}\nFc: {self.Fc}")
        print(f"Vida da correia: {self.t}")

    def calculate_transmission(self):
        self._calculate_geometry()
        self._input_constants()
        self._calculate_forces()
        self._calculate_durability()
        self._print_report()
