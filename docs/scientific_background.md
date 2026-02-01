# Scientific Background

This document provides the scientific context for the `pyrene-dimer-analyzer` package, focusing on the photophysics of pyrene and the geometric requirements for excimer formation.

## Pyrene Photophysics: Monomers and Excimers

Pyrene is a polycyclic aromatic hydrocarbon (PAH) well-known for its unique fluorescent properties. When a pyrene molecule absorbs a photon, it is promoted to an excited singlet state (S1). From this state, it can return to the ground state (S0) by emitting a photon, a process known as fluorescence.

-   **Monomer Emission**: An isolated, excited pyrene molecule exhibits a characteristic, highly structured fluorescence spectrum in the violet region (approximately 375-400 nm). This is referred to as **monomer emission**.

-   **Excimer Emission**: If an excited pyrene molecule encounters a ground-state pyrene molecule in close proximity and with a suitable orientation, they can form a transient, excited-state dimer called an **excimer** (a portmanteau of "excited dimer"). This excimer is only stable in the excited state and dissociates upon returning to the ground state. The excimer has its own distinct, broad, and structureless fluorescence spectrum shifted to a longer wavelength (blue/green region, approximately 460-500 nm). This is known as **excimer emission**.

The ratio of excimer to monomer emission intensity (Ie/Im) is highly sensitive to the local concentration and mobility of pyrene molecules, making pyrene a popular fluorescent probe for studying microenvironments, such as in cell membranes or polymers.

## Geometric Requirements for Excimer Formation

The formation of a stable pyrene excimer is not arbitrary; it requires a very specific geometric arrangement between the two pyrene rings. The key parameters, which this package is designed to calculate, are the plane-plane angle (θ), the inter-plane distance (d), and the degree of π-π overlap.

| Property                | Strong Excimer | Weak Excimer | Monomer Only |
| ----------------------- | -------------- | ------------ | ------------ |
| **θ (Plane-plane angle)** | < 20°          | 20-60°       | > 60°        |
| **d (Inter-plane distance)** | 3.3-3.7 Å      | 3.7-4.5 Å    | > 4.5 Å      |
| **π-overlap**           | > 70%          | 30-70%       | < 30%        |
| **Emission Wavelength** | 460-500 nm     | Mixed        | 375-400 nm   |

### 1. Plane-Plane Angle (θ)

The ideal geometry for excimer formation is a parallel, face-to-face arrangement of the two pyrene planes. As the angle between the planes increases, the efficiency of the electronic coupling required for excimer formation decreases dramatically.

> According to Birks (1970), an eclipsed geometry where the two pyrene molecules are perfectly aligned is optimal for excimer formation [1].

### 2. Inter-Plane Distance (d)

The distance between the two pyrene planes is critical. The optimal distance for π-stacking interactions that lead to excimer formation is in the range of 3.3 to 3.7 Å.

> Stevens (1968) determined the equilibrium distance for the pyrene crystal excimer to be 3.34 Å, providing a key benchmark for the ideal stacking distance [2].

If the planes are too far apart, the necessary orbital overlap cannot be achieved. If they are too close, repulsive forces dominate.

### 3. π-Overlap

Even with the correct angle and distance, significant overlap between the π-orbitals of the two pyrene rings is necessary. The percentage of overlap is a direct measure of the electronic communication between the two moieties.

> Recent studies have emphasized that the degree of π-overlap can be even more critical than the inter-plane distance for achieving effective excimer formation [4].

## Application to Covalently-Linked Dimers

In the context of this project, we are analyzing pyrene dimers that are covalently linked by a bridge. The chemical nature and steric bulk of the substituent groups (R-groups) on this bridge can impose conformational constraints that control the relative orientation of the two pyrene units.

**Hypothesis**: Bulkier R-groups (e.g., Ethyl < Isopropyl < Cyclohexyl < tert-Butyl) will introduce steric hindrance that forces the pyrene rings into a more perpendicular arrangement (larger θ), thereby inhibiting excimer formation and promoting monomer emission.

This package provides the tools to test this hypothesis by quantifying the geometric properties of conformers generated through molecular dynamics or conformational searches.

## References

[1] Birks, J.B. (1970). *Photophysics of Aromatic Molecules*. Wiley-Interscience.

[2] Stevens, B. (1968). Excimer fluorescence XII. The pyrene crystal excimer interaction potential. *Proceedings of the Royal Society of London. Series A*, 305, 55-70. [https://doi.org/10.1098/rspa.1968.0087](https://doi.org/10.1098/rspa.1968.0087)

[3] Poisson, L., et al. (2017). The mechanism of excimer formation: an experimental and theoretical study on the pyrene dimer. *Physical Chemistry Chemical Physics*, 19(35), 23492-23506. [https://doi.org/10.1039/C7CP03990E](https://doi.org/10.1039/C7CP03990E)

[4] Ge, Y., et al. (2020). A key stacking factor for the effective formation of pyrene excimer in crystals: degree of π–π overlap. *Journal of Materials Chemistry C*, 8(30), 10223-10232. [https://doi.org/10.1039/D0TC02562C](https://doi.org/10.1039/D0TC02562C)

[5] Bains, G.K., et al. (2012). The extent of pyrene excimer fluorescence emission is a reflector of distance and flexibility. *Biochemistry*, 51(31), 6207-6219. [https://doi.org/10.1021/bi3005285](https://doi.org/10.1021/bi3005285)
