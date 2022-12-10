import collections
import math
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from redutor import Material


class Shaft:
    def __init__(
        self,
        length: float,
        resolution: int,
        material: Material,
        sections: List[Tuple[float, float]],
        acting_forces: Dict[float, Tuple[float, float, float]],
        label: str,
    ) -> None:
        self.label = label
        self.x = np.array(
            [i for i in np.arange(0, length + length / resolution, 1 / resolution)]
        )
        d = []
        for pos in self.x:
            diameter = 0
            for i, sec in enumerate(sections):
                if i == (len(sections) - 1):
                    continue
                diameter = (
                    sec[1] if pos < sections[i + 1][0] and pos >= sec[0] else diameter
                )
            d.append(diameter)
        self.Diam = np.array(d)
        self.material = material
        self.J = np.array([math.pi * diam**4 for diam in self.Diam])
        self.I = np.array([math.pi * diam**4 for diam in self.Diam])
        print(self.I)
        self.acting_forces = collections.OrderedDict(sorted(acting_forces.items()))

    def _calculate_deflection(self):
        def _mod_trap_int(x, y) -> np.ndarray:
            result = [0]
            for i in range(0, len(x) - 1):
                value = result[i] + (y[i] + y[i + 1]) * (x[i + 1] - x[i]) / 2
                # print(value, y[i] + y[i + 1])
                result.append(value)
            # print(result)
            return np.array(result)

        Ix = self.material.elasticity_module * self.I * 10e12
        M_Ei_y = [my / ix for my, ix in zip(self.Mx, Ix)]
        M_Ei_z = [mz / ix for mz, ix in zip(self.Mx, Ix)]
        self.def_ang_y = _mod_trap_int(self.x, M_Ei_y)
        self.def_y = _mod_trap_int(self.x, self.def_ang_y)
        self.def_ang_z = _mod_trap_int(self.x, M_Ei_z)
        self.def_z = _mod_trap_int(self.x, self.def_ang_z)
        self.def_ang = [
            (dy**2 + dz**2) ** 0.5 for dy, dz in zip(self.def_ang_y, self.def_ang_z)
        ]
        print(M_Ei_y)

    def calculate_acting_forces(self):
        def _macaulay(x: float, n: float, direction: int) -> float:
            result = 0
            for a, force in self.acting_forces.items():
                if x >= a:
                    result += force[direction] * (x - a) ** n
            return result

        self.Vx = [_macaulay(x=x, n=0, direction=0) for x in self.x]
        self.Mx = [_macaulay(x=x, n=1, direction=0) for x in self.x]

        self.Vy = [_macaulay(x=x, n=0, direction=1) for x in self.x]
        self.My = [_macaulay(x=x, n=1, direction=1) for x in self.x]

        self.Vz = [_macaulay(x=x, n=0, direction=2) for x in self.x]
        self.Mz = [_macaulay(x=x, n=1, direction=2) for x in self.x]

        self.V = [
            (vx**2 + vy**2 + vz**2) ** 0.5
            for vx, vy, vz in zip(self.Vx, self.Vy, self.Vz)
        ]
        self.M = [
            (mx**2 + my**2 + mz**2) ** 0.5
            for mx, my, mz in zip(self.Mx, self.My, self.Mz)
        ]
        self._calculate_deflection()

    def export_plots(self):
        x = self.x * 1000
        my = [my * 1000 for my in self.My]
        mz = [mz * 1000 for mz in self.Mz]
        m = [m * 1000 for m in self.M]

        plt.plot(x, self.Vz)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Esforço cortante V [N] no plano x-z")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Vz.png")
        plt.close()

        plt.plot(x, self.Vy)
        plt.ylabel("Esforço cortante V [N] no plano x-y")
        plt.xlabel("Distância z [mm]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Vy.png")
        plt.close()

        plt.plot(x, self.V)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Magnitude total do esforço cortante [N]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_V.png")
        plt.close()

        plt.plot(x, my)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Momento na direção z (N.mm)")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_My.png")
        plt.close()

        plt.plot(x, mz)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Momento na direção y (N.mm)")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Mz.png")
        plt.close()

        plt.plot(x, m)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Magnitude total do momento [N.mm]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_M.png")
        plt.close()

        plt.plot(x, self.def_ang, color="blue")
        plt.plot(x, self.def_ang_z, color="red")
        plt.plot(x, self.def_ang_y, color="black")
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Deflexão não corrigida")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Def_nao_cor.png")
        plt.close()

    def calculate_stress(self):
        return
