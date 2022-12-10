import matplotlib.pyplot as plt

from redutor import (
    Gear,
    GearTransmission,
    Material,
    Pulley,
    PulleyTransmission,
    Shaft,
    SystemVariables,
)


def main():
    systemVar = SystemVariables(
        input_power=49370,
        input_velocity=2335,
        seconds_of_use=365 * (2 / 7) * 8 * 60 * 60,
        roller_efficiency=0.96,
    )
    steel = Material(
        elasticity_module=190e9,
        poisson_coef=0.3,
        yield_stress=580e6,
        ultimate_stress=690e6,
        bending_stress_strength=230e6,
        contact_stress_strength=700e6,
    )

    # Transmissao polia
    Polia1 = Pulley(primitive_diameter=0.118, input_velocity=systemVar.input_velocity)
    Polia2 = Pulley(
        primitive_diameter=0.236, input_velocity=systemVar.input_velocity / 2
    )
    Transm_Polia = PulleyTransmission(
        polia1=Polia1, polia2=Polia2, power=systemVar.input_power
    )
    Transm_Polia.calculate_transmission()

    # Transmissao engrenagem
    Engrenagem1B = Gear(
        number_of_teeths=19,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.24,
        material=steel,
    )
    Engrenagem2A = Gear(
        number_of_teeths=57,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.28,
        material=steel,
    )
    Engrenagem2B = Gear(
        number_of_teeths=19,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.24,
        material=steel,
    )
    Engrenagem3A = Gear(
        number_of_teeths=95,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.29,
        material=steel,
    )
    Engrenagem3B = Gear(
        number_of_teeths=43,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.27,
        material=steel,
    )
    Engrenagem4A = Gear(
        number_of_teeths=86,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=12,
        J_bending_stress=0.29,
        material=steel,
    )
    Transmissao1_2 = GearTransmission(
        gear1=Engrenagem1B, gear2=Engrenagem2A, seconds_of_use=systemVar.seconds_of_use
    )
    Transmissao1_2.calculate_forces(
        input_power=systemVar.input_power * systemVar.roller_efficiency,
        input_velocity=systemVar.input_velocity / 2,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao1_2.calculate_stress()
    Transmissao2_3 = GearTransmission(
        gear1=Engrenagem2B, gear2=Engrenagem3A, seconds_of_use=systemVar.seconds_of_use
    )
    Transmissao2_3.calculate_forces(
        input_power=systemVar.input_power * systemVar.roller_efficiency,
        input_velocity=systemVar.input_velocity / 2,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao2_3.calculate_stress()
    Transmissao3_4 = GearTransmission(
        gear1=Engrenagem3B, gear2=Engrenagem4A, seconds_of_use=systemVar.seconds_of_use
    )
    Transmissao3_4.calculate_forces(
        input_power=systemVar.input_power * systemVar.roller_efficiency,
        input_velocity=systemVar.input_velocity / 2,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao3_4.calculate_stress()

    Eixo1 = Shaft(
        length=1,
        resolution=10,
        material=steel,
        sections=[
            [0, 0.06],
            [0.025, 0.09],
            [0.155, 0.065],
            [0.190, 0.050],
            [0.300, 0.050],
        ],
        acting_forces={1: (0, 10, 0), 3: (0, 20, 0), 5: (0, -50, 0), 7: (0, 20, 0)},
        label="Eixo 1",
    )
    Eixo1.calculate_acting_forces()
    Eixo1.export_plots()


if __name__ == "__main__":
    main()
