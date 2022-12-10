import math


class Pulley:
    def __init__(self, primitive_diameter: float, input_velocity: float) -> None:
        self.primitive_diameter = primitive_diameter
        self.angular_velocity = input_velocity * math.pi / 30
