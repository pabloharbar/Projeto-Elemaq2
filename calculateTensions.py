import math

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
        input_power=2937,
        input_velocity=2335,
        seconds_of_use=365 * (2 / 7) * 8 * 60 * 60,
        roller_efficiency=0.96,
        belt_efficiency=0.96,
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
        polia1=Polia1, polia2=Polia2, power=systemVar.input_power, position=0.195
    )
    Transm_Polia.calculate_transmission()

    # Transmissao engrenagem
    Engrenagem1B = Gear(
        number_of_teeths=19,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.24,
        material=steel,
    )
    Engrenagem2A = Gear(
        number_of_teeths=57,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.28,
        material=steel,
    )
    Engrenagem2B = Gear(
        number_of_teeths=19,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.24,
        material=steel,
    )
    Engrenagem3A = Gear(
        number_of_teeths=95,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.29,
        material=steel,
    )
    Engrenagem3B = Gear(
        number_of_teeths=43,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.27,
        material=steel,
    )
    Engrenagem4A = Gear(
        number_of_teeths=86,
        pressure_angle=20,
        modulo=0.003,
        thickness_factor=14,
        J_bending_stress=0.29,
        material=steel,
    )
    Transmissao1_2 = GearTransmission(
        gear1=Engrenagem1B,
        gear2=Engrenagem2A,
        seconds_of_use=systemVar.seconds_of_use,
        position=0.0975,
    )
    Transmissao1_2.calculate_forces(
        input_power=systemVar.input_power
        * systemVar.roller_efficiency
        * systemVar.belt_efficiency,
        input_velocity=systemVar.input_velocity / 2,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao1_2.calculate_stress()
    Transmissao2_3 = GearTransmission(
        gear1=Engrenagem2B,
        gear2=Engrenagem3A,
        seconds_of_use=systemVar.seconds_of_use,
        position=0.0475,
    )
    Transmissao2_3.calculate_forces(
        input_power=Transmissao1_2.p2,
        input_velocity=systemVar.input_velocity / 6,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao2_3.calculate_stress()
    Transmissao3_4 = GearTransmission(
        gear1=Engrenagem3B,
        gear2=Engrenagem4A,
        seconds_of_use=systemVar.seconds_of_use,
        position=0.0975,
    )
    Transmissao3_4.calculate_forces(
        input_power=Transmissao2_3.p2,
        input_velocity=systemVar.input_velocity / 30,
        roller_efficiency=systemVar.roller_efficiency,
    )
    Transmissao3_4.calculate_stress()

    # Definição das reações
    Ra_distance = 0.133
    Rb_distance = 0.014

    # Eixo 1
    proj = math.cos(math.radians(9.38))
    Ty = -(Transm_Polia.F1 * proj + Transm_Polia.F2 * proj)

    Ra1y = (
        Transmissao1_2.Fn1 * (Transmissao1_2.position - Rb_distance)
        - Ty * (Transm_Polia.position - Rb_distance)
    ) / (Ra_distance - Rb_distance)
    Rb1y = Transmissao1_2.Fn1 - Ra1y - Ty

    Ra1x = -(
        Transmissao1_2.Ft1
        * (Transmissao1_2.position - Rb_distance)
        / (Ra_distance - Rb_distance)
    )
    Rb1x = -Transmissao1_2.Ft1 - Ra1x

    # Eixo 2
    Ra2y = (
        -Transmissao1_2.Fn2 * (Transmissao1_2.position - Rb_distance)
        + Transmissao2_3.Ft1 * (Transmissao2_3.position - Rb_distance)
    ) / (Ra_distance - Rb_distance)
    Rb2y = Transmissao2_3.Ft1 - Ra2y - Transmissao1_2.Fn2

    Ra2x = (
        Transmissao1_2.Ft2 * (Transmissao1_2.position - Rb_distance)
        - Transmissao2_3.Fn1 * (Transmissao2_3.position - Rb_distance)
    ) / (Ra_distance - Rb_distance)
    Rb2x = Transmissao1_2.Ft2 - Ra2x - Transmissao2_3.Fn1

    # Eixo 3
    Ra3y = -(
        Transmissao2_3.Ft2 * (Transmissao2_3.position - Rb_distance)
        + Transmissao3_4.Ft1 * (Transmissao3_4.position - Rb_distance)
    ) / (Ra_distance - Rb_distance)
    Rb3y = -Ra3y - Transmissao2_3.Ft2 - Transmissao3_4.Ft1

    Ra3x = (
        Transmissao2_3.Fn2 * (Transmissao2_3.position - Rb_distance)
        - Transmissao3_4.Fn1 * (Transmissao3_4.position - Rb_distance)
    ) / (Ra_distance - Rb_distance)
    Rb3x = -Ra3x + Transmissao2_3.Fn2 - Transmissao3_4.Fn1

    # Eixo 4
    Ra4y = (
        Transmissao3_4.Ft2
        * (Transmissao3_4.position - Rb_distance)
        / (Ra_distance - Rb_distance)
    )
    Rb4y = Transmissao3_4.Ft2 - Ra4y

    Ra4x = (
        Transmissao3_4.Fn2
        * (Transmissao3_4.position - Rb_distance)
        / (Ra_distance - Rb_distance)
    )
    Rb4x = Transmissao3_4.Fn2 - Ra4x

    print(f"Ra1y :{Ra1y}, Rb1y :{Rb1y}, T1y + T2y : {Ty}, W21y : {-Transmissao1_2.Fn1}")
    print(f"Ra1x : {Ra1x}, Rb1x: {Rb1x}, W21x : {Transmissao1_2.Ft1}")
    print(
        f"Ra2y : {Ra2y}, Rb2y : {Rb2y}, W12y : {Transmissao1_2.Fn2}, W32y : {-Transmissao2_3.Ft1}"
    )
    print(
        f"Ra2x : {Ra2x}, Rb2x : {Rb2x}, W12x : {-Transmissao1_2.Ft2}, W32x : {Transmissao2_3.Fn1}"
    )
    print(
        f"Ra3y : {Ra3y}, Rb3y : {Rb3y}, W23y : {Transmissao2_3.Ft2}, W43y : {Transmissao3_4.Ft1}"
    )
    print(
        f"Ra3x : {Ra3x}, R3bx : {Rb3x}, W23x : {-Transmissao2_3.Fn2}, W43x : {Transmissao3_4.Fn1}"
    )
    print(f"Ra4y : {Ra4y}, Rb4y : {Rb4y}, W34y : {-Transmissao3_4.Ft2}")
    print(f"Ra4x : {Ra4x}, Rb4x : {Rb4x}, W34y : {-Transmissao3_4.Fn2}")

    Eixo1 = Shaft(
        length=0.21,
        resolution=1000,
        material=steel,
        sections=[
            [0, 0.017],
            [0.02, 0.02],
            [0.14, 0.023],
            [0.185, 0.019],
        ],
        acting_forces={
            Rb_distance: (Rb1x, Rb1y, 0),  # R1B
            Transmissao1_2.position: (Transmissao1_2.Ft1, -Transmissao1_2.Fn1, 0),
            Ra_distance: (Ra1x, Ra1y, 0),
            Transm_Polia.position: (0, Ty, 0),
        },
        label="Eixo 1",
        correction_points=[Rb_distance, Ra_distance],
        Torque=Transmissao1_2.T1,
        stress_focus=[[0.0895, 0.1055], [0.19, 0.21]],
    )
    Eixo1.calculate_acting_forces()
    Eixo1.calculate_stress()
    Eixo1.export_plots()

    Eixo2 = Shaft(
        length=0.16,
        resolution=1000,
        material=steel,
        sections=[
            [0, 0.017],
            [0.02, 0.02],
            [0.14, 0.023],
        ],
        acting_forces={
            Rb_distance: (Rb2x, Rb2y, 0),  # R1B
            Transmissao1_2.position: (-Transmissao1_2.Ft2, Transmissao1_2.Fn2, 0),
            Ra_distance: (Ra2x, Ra2y, 0),
            Transmissao2_3.position: (Transmissao2_3.Fn1, -Transmissao2_3.Ft1, 0),
        },
        label="Eixo 2",
        correction_points=[Rb_distance, Ra_distance],
        Torque=Transmissao2_3.T1,
        stress_focus=[[0.0395, 0.0555], [0.0895, 0.1055]],
    )
    Eixo2.calculate_acting_forces()
    Eixo2.calculate_stress()
    Eixo2.export_plots()

    Eixo3 = Shaft(
        length=0.16,
        resolution=1000,
        material=steel,
        sections=[
            [0, 0.017],
            [0.02, 0.02],
            [0.14, 0.023],
        ],
        acting_forces={
            Rb_distance: (Rb3x, Rb3y, 0),  # R1B
            Transmissao2_3.position: (-Transmissao2_3.Fn2, Transmissao2_3.Ft2, 0),
            Ra_distance: (Ra3x, Ra3y, 0),
            Transmissao3_4.position: (Transmissao3_4.Fn1, Transmissao3_4.Ft1, 0),
        },
        label="Eixo 3",
        correction_points=[Rb_distance, Ra_distance],
        Torque=Transmissao3_4.T1,
        stress_focus=[[0.0395, 0.0555], [0.0895, 0.1055]],
    )
    Eixo3.calculate_acting_forces()
    Eixo3.calculate_stress()
    Eixo3.export_plots()

    Eixo4 = Shaft(
        length=0.16,
        resolution=1000,
        material=steel,
        sections=[
            [0, 0.017],
            [0.02, 0.02],
            [0.14, 0.023],
        ],
        acting_forces={
            Rb_distance: (Rb4x, Rb4y, 0),  # R1B
            Ra_distance: (Ra4x, Ra4y, 0),
            Transmissao3_4.position: (-Transmissao3_4.Fn2, -Transmissao3_4.Ft2, 0),
        },
        label="Eixo 4",
        correction_points=[Rb_distance, Ra_distance],
        Torque=Transmissao3_4.T2,
        stress_focus=[[0.0895, 0.1055]],
    )
    Eixo4.calculate_acting_forces()
    Eixo4.calculate_stress()
    Eixo4.export_plots()

    print(
        f"Transmissão:\nw0={systemVar.input_velocity * math.pi / 30}, P0={systemVar.input_power}"
    )
    print(f"w1={Transmissao1_2.w1}, P1={Transmissao1_2.p1}, T1={Transmissao1_2.T1}")
    print(f"w2={Transmissao1_2.w2}, P2={Transmissao1_2.p2}, T2={Transmissao1_2.T2}")
    print(f"w3={Transmissao2_3.w2}, P3={Transmissao3_4.p1}, T3={Transmissao2_3.T2}")
    print(f"w4={Transmissao3_4.w2}, P4={Transmissao3_4.p2}, T4={Transmissao3_4.T2}")


if __name__ == "__main__":
    main()
