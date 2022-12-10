import math

from redutor import Material


class Gear:
    def __init__(
        self,
        number_of_teeths: float,
        pressure_angle: float,
        modulo: float,
        thickness_factor: float,
        J_bending_stress: float,
        material: Material,
    ) -> None:
        self.primitive_diam = number_of_teeths * modulo
        self.pressure_angle = math.radians(pressure_angle)
        self.modulo = modulo
        self.thickness = modulo * thickness_factor
        self.number_of_teeths = number_of_teeths
        self.J = J_bending_stress
        self.material = material
