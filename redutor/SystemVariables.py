class SystemVariables:
    def __init__(
        self,
        input_power: float,
        input_velocity: float,
        seconds_of_use: float,
        roller_efficiency: float,
        belt_efficiency: float,
    ) -> None:
        self.input_power = input_power
        self.input_velocity = input_velocity
        self.seconds_of_use = seconds_of_use
        self.roller_efficiency = roller_efficiency
        self.belt_efficiency = belt_efficiency
