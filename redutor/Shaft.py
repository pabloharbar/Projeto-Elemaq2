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
        correction_points: Tuple[float, float],
        Torque: float,
        stress_focus: List[Tuple[float, float]],
    ) -> None:
        self.stress_focus = stress_focus
        self.Torque = Torque
        self.correction_points = correction_points
        self.label = label
        self.z = np.array(
            [i for i in np.arange(0, length + length / resolution, 1 / resolution)]
        )
        d = []
        for pos in self.z:
            diameter = 0
            for i, sec in enumerate(sections):
                if i == (len(sections) - 1):
                    diameter = sec[1] if pos >= sec[0] else diameter
                else:
                    diameter = (
                        sec[1]
                        if pos < sections[i + 1][0] and pos >= sec[0]
                        else diameter
                    )
            d.append(diameter)
        self.Diam = np.array(d)
        self.material = material
        self.J = np.array([math.pi * diam**4 / 32 for diam in self.Diam])
        self.I = np.array([math.pi * diam**4 / 64 for diam in self.Diam])
        self.acting_forces = collections.OrderedDict(sorted(acting_forces.items()))

    def _correct_deflection(self):
        index0 = np.where(self.z == self.correction_points[0])[0]
        index1 = np.where(self.z == self.correction_points[1])[0]
        C3x = math.atan(
            (self.def_x[index1] - self.def_x[index0])
            / (self.correction_points[1] - self.correction_points[0])
        )
        C3y = math.atan(
            (self.def_y[index1] - self.def_y[index0])
            / (self.correction_points[1] - self.correction_points[0])
        )
        C3z = math.atan(
            (self.def_z[index1] - self.def_z[index0])
            / (self.correction_points[1] - self.correction_points[0])
        )
        self.def_tot_cor_x = self.def_x - C3x * self.def_x
        self.def_tot_cor_y = self.def_y - C3y * self.def_y
        self.def_tot_cor_z = self.def_z - C3z * self.def_z
        self.def_tot_cor = (
            self.def_tot_cor_x**2 + self.def_tot_cor_y**2 + self.def_tot_cor_z**2
        ) ** 0.5

        self.def_ang_cor_x = self.def_ang_x - C3x
        self.def_ang_cor_y = self.def_ang_y - C3y
        self.def_ang_cor_z = self.def_ang_z - C3z
        self.def_ang_cor = (
            self.def_ang_cor_x**2 + self.def_ang_cor_z**2 + self.def_ang_cor_y**2
        ) ** 0.5

    def _calculate_deflection(self):
        def _mod_trap_int(x, y) -> np.ndarray:
            result = [0]
            for i in range(0, len(x) - 1):
                value = result[i] + (y[i] + y[i + 1]) * (x[i + 1] - x[i]) / 2
                result.append(value)
            # print(result)
            return np.array(result)

        Iz = self.material.elasticity_module * self.I
        M_Ei_x = [mx / iz for mx, iz in zip(self.Mx, Iz)]
        M_Ei_y = [my / iz for my, iz in zip(self.My, Iz)]
        M_Ei_z = [mz / iz for mz, iz in zip(self.Mz, Iz)]
        self.def_ang_x = _mod_trap_int(self.z, M_Ei_x)
        self.def_x = _mod_trap_int(self.z, self.def_ang_x)
        self.def_ang_y = _mod_trap_int(self.z, M_Ei_y)
        self.def_y = _mod_trap_int(self.z, self.def_ang_y)
        self.def_ang_z = _mod_trap_int(self.z, M_Ei_z)
        self.def_z = _mod_trap_int(self.z, self.def_ang_z)
        self.def_ang = [
            (dx**2 + dy**2 + dz**2) ** 0.5
            for dx, dy, dz in zip(self.def_ang_x, self.def_ang_y, self.def_ang_z)
        ]
        self._correct_deflection()

    def calculate_acting_forces(self):
        def _macaulay(x: float, n: float, direction: int) -> float:
            result = 0
            for a, force in self.acting_forces.items():
                if x >= a:
                    result += force[direction] * (x - a) ** n
            return result

        self.Vx = [_macaulay(x=x, n=0, direction=0) for x in self.z]
        self.Mx = [_macaulay(x=x, n=1, direction=0) for x in self.z]

        self.Vy = [_macaulay(x=x, n=0, direction=1) for x in self.z]
        self.My = [_macaulay(x=x, n=1, direction=1) for x in self.z]

        self.Vz = [_macaulay(x=x, n=0, direction=2) for x in self.z]
        self.Mz = [_macaulay(x=x, n=1, direction=2) for x in self.z]

        self.V = [
            (vx**2 + vy**2 + vz**2) ** 0.5
            for vx, vy, vz in zip(self.Vx, self.Vy, self.Vz)
        ]
        self.M = [
            (mx**2 + my**2 + mz**2) ** 0.5
            for mx, my, mz in zip(self.Mx, self.My, self.Mz)
        ]
        self._calculate_deflection()

    def _evaluate_fatigue(self):
        self.tensao_alt = self.sigma
        self.tensao_med = (self.tau_xy * 3) ** 0.5
        Cs = 1.189 * (self.Diam * 1000) ** (-0.097)
        Ce = 1
        Cf = min(1, 4.51 * (self.material.ultimate_stress / 10e6) ** (-0.265))
        Ct = 1
        Cr = 0.868
        tensao_fad = 0.5 * self.material.ultimate_stress * Cs * Ce * Cf * Ct * Cr
        print(np.mean(Cs), Ce, Cf, Ct, Cr)
        # Linha de carga 3
        self.N_fad = (
            tensao_fad
            * self.material.ultimate_stress
            / (
                self.tensao_alt * self.material.ultimate_stress
                + self.tensao_med * tensao_fad
            )
        )

    def _evaluate_stress_focus(self):
        q = 1 / (1 + self.material.neuber_constant / ((self.Diam / 2) ** 0.5))
        qs = 1 / (1 + self.material.neuber_constant_shear / ((self.Diam / 2) ** 0.5))

        # Chavetas
        Ktc = 2.1
        Ktsc = 3.0
        Kt = []
        Kts = []
        for pos in self.z:
            kt = 0
            kts = 0
            for focus in self.stress_focus:
                kt = Ktc if pos >= focus[0] and pos <= focus[1] else kt
                kts = Ktsc if pos >= focus[0] and pos <= focus[1] else kts

            Kt.append(kt)
            Kts.append(kts)

        Kf = 1 + q * (np.array(Kt) - 1)
        Kfs = 1 + qs * (np.array(Kts) - 1)
        self.sigma_x = self.sigma_x * Kf
        self.sigma_y = self.sigma_y * Kf
        self.sigma_z = self.sigma_z * Kf
        self.sigma = self.sigma * Kf
        self.tau_xy = self.tau_xy * Kfs

        # Criterio de falha estatico
        self.sigma_eq = (
            (
                (self.sigma_x - self.sigma_y) ** 2
                + self.sigma_x**2
                + self.sigma_y**2
                + 6 * self.tau_xy**2
            )
            / 2
        ) ** 0.5
        self.N_est = self.material.yield_stress / self.sigma_eq

        self._evaluate_fatigue()

    def calculate_stress(self):
        self.sigma_x = self.Mx * self.Diam / (2 * self.I)
        self.sigma_y = self.My * self.Diam / (2 * self.I)
        self.sigma_z = self.Mz * self.Diam / (2 * self.I)
        self.sigma = self.M * self.Diam / (2 * self.I)
        self.tau_xy = self.Torque * self.Diam / (2 * self.J)

        self._evaluate_stress_focus()

    def export_plots(self):
        x = self.z * 1000
        mx = [mx * 1000 for mx in self.Mx]
        my = [my * 1000 for my in self.My]
        mz = [mz * 1000 for mz in self.Mz]
        m = [m * 1000 for m in self.M]

        plt.plot(x, self.Diam)
        plt.ylim(0, max(self.Diam) * 1.1)
        plt.savefig(f"./output/{self.label}_Geometry.png")
        plt.close()

        plt.plot(x, self.Vz)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Esforço cortante V [N] no plano z-z")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Vz.png")
        plt.close()

        plt.plot(x, self.Vy)
        plt.ylabel("Esforço cortante V [N] no plano z-y")
        plt.xlabel("Distância z [mm]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Vy.png")
        plt.close()

        plt.plot(x, self.Vx)
        plt.ylabel("Esforço cortante V [N] no plano z-x")
        plt.xlabel("Distância z [mm]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Vx.png")
        plt.close()

        plt.plot(x, self.V)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Magnitude total do esforço cortante [N]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_V.png")
        plt.close()

        plt.plot(x, mx)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Momento na direção y (N.mm)")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Mx.png")
        plt.close()

        plt.plot(x, my)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Momento na direção x (N.mm)")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_My.png")
        plt.close()

        plt.plot(x, mz)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Momento na direção z (N.mm)")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_Mz.png")
        plt.close()

        plt.plot(x, m)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Magnitude total do momento [N.mm]")
        plt.xlim(0, max(x))
        plt.savefig(f"./output/{self.label}_M.png")
        plt.close()

        plt.plot(x, self.def_x, color="blue", label="def_x")
        plt.plot(x, self.def_y, color="red", label="def_y")
        plt.plot(x, self.def_z, color="black", label="def_z")
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Deflexão não corrigida")
        plt.xlim(0, max(x))
        plt.legend()
        plt.savefig(f"./output/{self.label}_Def_nao_cor.png")
        plt.close()

        plt.plot(x, self.def_tot_cor, color="blue", label="def_tot_cor")
        plt.plot(x, self.def_tot_cor_x, color="yellow", label="def_x_cor")
        plt.plot(x, self.def_tot_cor_y, color="red", label="def_y_cor")
        plt.plot(x, self.def_tot_cor_z, color="black", label="def_z_cor")
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Deflexão total corrigida")
        plt.xlim(0, max(x))
        plt.legend()
        plt.savefig(f"./output/{self.label}_Def_cor.png")
        plt.close()

        plt.plot(x, self.def_ang, color="blue", label="def_ang")
        plt.plot(x, self.def_ang_z, color="red", label="def_ang_z")
        plt.plot(x, self.def_ang_y, color="black", label="def_ang_y")
        plt.plot(x, self.def_ang_x, color="yellow", label="def_ang_x")
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Deflexão angular não corrigida")
        plt.xlim(0, max(x))
        plt.legend()
        plt.savefig(f"./output/{self.label}_Def_ang_nao_cor.png")
        plt.close()

        plt.plot(x, self.def_ang_cor, color="blue", label="def_ang_cor")
        plt.plot(x, self.def_ang_cor_x, color="yellow", label="def_ang_cor_x")
        plt.plot(x, self.def_ang_cor_y, color="red", label="def_ang_cor_y")
        plt.plot(x, self.def_ang_cor_z, color="black", label="def_ang_cor_z")
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Deflexão angular corrigida")
        plt.xlim(0, max(x))
        plt.legend()
        plt.savefig(f"./output/{self.label}_Def_ang_cor.png")
        plt.close()

        plt.plot(x, self.N_est)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Coeficiente de segurança estático")
        plt.xlim(0, max(x))
        plt.legend()
        plt.savefig(f"./output/{self.label}_N-estatico.png")
        plt.close()

        plt.plot(x, self.N_fad)
        plt.xlabel("Distância z [mm]")
        plt.ylabel("Coeficiente de segurança de fadiga")
        plt.xlim(0, max(x))
        plt.ylim(0, 50)
        plt.legend()
        plt.savefig(f"./output/{self.label}_N-fadiga.png")
        plt.close()
