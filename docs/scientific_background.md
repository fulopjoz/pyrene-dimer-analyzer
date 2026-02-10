# Scientific Background

This document provides the scientific context for the `pyrene-dimer-analyzer` package, focusing on the photophysics of aromatic excimers and the geometric requirements for their formation across multiple chromophore systems.

## Aromatic Excimer Photophysics

### Monomers and Excimers

Polycyclic aromatic hydrocarbons (PAHs) such as pyrene, perylene, and anthracene are well-known for their unique fluorescent properties. When an aromatic molecule absorbs a photon, it is promoted to an excited singlet state (S1). From this state, it can return to the ground state (S0) by emitting a photon (fluorescence).

-   **Monomer Emission**: An isolated, excited aromatic molecule exhibits a characteristic, highly structured fluorescence spectrum. For pyrene, this appears in the violet region (approximately 375-400 nm).

-   **Excimer Emission**: If an excited aromatic molecule encounters a ground-state aromatic molecule in close proximity and with a suitable orientation, they can form a transient, excited-state dimer called an **excimer** (a portmanteau of "excited dimer"). This excimer is only stable in the excited state and dissociates upon returning to the ground state. The excimer has its own distinct, broad, and structureless fluorescence spectrum shifted to a longer wavelength. For pyrene, this appears in the blue/green region (approximately 460-500 nm).

The ratio of excimer to monomer emission intensity (Ie/Im) is highly sensitive to the local concentration and mobility of the aromatic molecules, making excimer-forming chromophores popular fluorescent probes.

## Geometric Requirements for Excimer Formation

The formation of a stable aromatic excimer is not arbitrary; it requires a very specific geometric arrangement between the two aromatic rings. The key parameters, which this package is designed to calculate, are the plane-plane angle (theta), the inter-plane distance (d), and the degree of pi-pi overlap.

### Classification Thresholds (Pyrene, v1.1.0)

| Property                      | Strong Excimer | Weak Excimer | Monomer Only |
| ----------------------------- | -------------- | ------------ | ------------ |
| **theta (Plane-plane angle)** | < 20 deg       | 20-60 deg    | > 60 deg     |
| **d (Inter-plane distance)**  | 3.3-3.7 A      | 3.7-4.5 A   | > 4.5 A      |
| **pi-overlap**                | > 50%          | 30-50%       | < 30%        |
| **Emission Wavelength**       | 460-500 nm     | Mixed        | 375-400 nm   |

> **Note (v1.1.0)**: The strong overlap threshold was corrected from 70% to 50% based on recent literature. See the "Literature Corrections" section below.

### 1. Plane-Plane Angle (theta)

The ideal geometry for excimer formation is a parallel, face-to-face arrangement of the two aromatic planes. As the angle between the planes increases, the efficiency of the electronic coupling required for excimer formation decreases dramatically.

> According to Birks (1970), an eclipsed geometry where the two aromatic molecules are perfectly aligned is optimal for excimer formation [1].

Recent work by Mazzeo et al. (2024) has shown that **twisted excimers** at 50-60 deg can actually be the dominant excimer species in certain systems, challenging the traditional view that only near-parallel arrangements form excimers [7]. These are captured by the "weak excimer" classification.

### 2. Inter-Plane Distance (d)

The distance between the two aromatic planes is critical. The optimal distance for pi-stacking interactions that lead to excimer formation is system-dependent:

-   **Pyrene**: 3.3-3.7 A (Birks & Kazzaz 1968: equilibrium 3.34 A [2]; Marazzi et al. 2024: 3.24 A by TD-DFT [8])
-   **Perylene**: 3.4-3.8 A (crystal stacking data)
-   **Anthracene**: 3.4-3.93 A (high-pressure excimer literature)

If the planes are too far apart, the necessary orbital overlap cannot be achieved. If they are too close, repulsive forces dominate.

#### High-Angle Distance Warning

**Important**: When theta > 60 deg, the inter-plane distance metric becomes unreliable. At high angles, the perpendicular distance from the centroid of one ring to the plane of the other ring measures the **edge-to-face distance** rather than the pi-stacking distance. This can produce physically misleading values (e.g., d = 0.52 A at theta = 86.3 deg).

As of v1.1.0, the tool emits a `UserWarning` when theta > 60 deg and records the warning in the `geometry_warnings` output field. Users should use `centroid_distance` instead for high-angle conformers.

### 3. pi-Overlap

Even with the correct angle and distance, significant overlap between the pi-orbitals of the two aromatic rings is necessary. The percentage of overlap is a direct measure of the electronic communication between the two moieties.

> Ge et al. (2020) showed that the degree of pi-overlap can be even more critical than the inter-plane distance for achieving effective excimer formation in crystalline systems, with effective excimers observed across a 40-80% overlap range [4].

## Literature Corrections in v1.1.0

### Overlap Threshold: 70% to 50%

The original v1.0 implementation used a strong overlap threshold of > 70%. This was corrected to > 50% based on:

1.  **Ge et al. (2020)** [4]: Studied pyrene excimer formation in crystals and found effective excimers across a broad 40-80% overlap range, with no sharp cutoff at 70%.

2.  **Basuroy et al. (2021)** [6]: Demonstrated in *J. Chem. Phys.* that a pyrene dimer geometry with only 42% pi-overlap still formed an excimer (confirmed by fluorescence spectroscopy). This directly contradicts a 70% threshold.

3.  **Mazzeo et al. (2024)** [7]: Showed that stacked-twisted excimers at 50-60 deg are the dominant excimer species in certain pyrene derivatives, further supporting that strict geometric criteria should be relaxed.

The corrected threshold of 50% captures the majority of experimentally confirmed excimer geometries while still distinguishing them from clearly non-excimeric arrangements.

### Twisted Excimers

Classical excimer theory (Birks 1970) assumed near-parallel stacking as a requirement. However, Mazzeo et al. (2024) demonstrated through computational and experimental evidence that twisted excimers (theta = 50-60 deg) can be energetically favorable and are sometimes the dominant excimer species [7]. These are correctly classified as "weak_excimer" by the tool.

## Multi-System Aromatic Excimers

The v1.1.0 release extends support beyond pyrene to multiple aromatic systems. Each system has its own classification thresholds based on available literature data:

| System       | Strong d (A) | Strong theta | Strong overlap | Weak d max | Min ring atoms | Key Reference                      |
| ------------ | ------------ | ------------ | -------------- | ---------- | -------------- | ---------------------------------- |
| Pyrene       | 3.3-3.7      | < 20 deg     | > 50%          | 4.5 A      | 10             | Birks 1968, Ge 2020, Basuroy 2021  |
| Perylene     | 3.4-3.8      | < 20 deg     | > 50%          | 4.5 A      | 14             | Crystal stacking data              |
| Anthracene   | 3.4-3.93     | < 20 deg     | > 50%          | 4.5 A      | 10             | High-pressure excimer lit.         |
| Naphthalene  | 3.3-3.7      | < 20 deg     | > 50%          | 4.5 A      | 8              | Approximated from pyrene           |
| Phenanthrene | 3.3-3.7      | < 20 deg     | > 50%          | 4.5 A      | 10             | Approximated from pyrene           |

**Note**: Naphthalene and phenanthrene thresholds are approximated from pyrene data due to limited crystal excimer structure data. Users working with these systems should validate against their own experimental data and can provide custom thresholds via the `custom_thresholds` parameter.

Custom aromatic systems can be registered at runtime using `register_system()` with user-defined SMARTS patterns and thresholds.

## Application to Covalently-Linked Dimers

In the context of this project, we are analyzing aromatic dimers that are covalently linked by a bridge. The chemical nature and steric bulk of the substituent groups (R-groups) on this bridge can impose conformational constraints that control the relative orientation of the two aromatic units.

**Hypothesis**: Bulkier R-groups (e.g., Ethyl < Isopropyl < Cyclohexyl < tert-Butyl) will introduce steric hindrance that forces the aromatic rings into a more perpendicular arrangement (larger theta), thereby inhibiting excimer formation and promoting monomer emission.

This package provides the tools to test this hypothesis by quantifying the geometric properties of conformers generated through molecular dynamics or conformational searches.

## References

[1] Birks, J.B. (1970). *Photophysics of Aromatic Molecules*. Wiley-Interscience.

[2] Stevens, B. (1968). Excimer fluorescence XII. The pyrene crystal excimer interaction potential. *Proceedings of the Royal Society of London. Series A*, 305, 55-70. [https://doi.org/10.1098/rspa.1968.0087](https://doi.org/10.1098/rspa.1968.0087)

[3] Poisson, L., et al. (2017). The mechanism of excimer formation: an experimental and theoretical study on the pyrene dimer. *Physical Chemistry Chemical Physics*, 19(35), 23492-23506. [https://doi.org/10.1039/C7CP03990E](https://doi.org/10.1039/C7CP03990E)

[4] Ge, Y., et al. (2020). A key stacking factor for the effective formation of pyrene excimer in crystals: degree of pi-pi overlap. *Journal of Materials Chemistry C*, 8(30), 10223-10232. [https://doi.org/10.1039/D0TC02562C](https://doi.org/10.1039/D0TC02562C)

[5] Bains, G.K., et al. (2012). The extent of pyrene excimer fluorescence emission is a reflector of distance and flexibility. *Biochemistry*, 51(31), 6207-6219. [https://doi.org/10.1021/bi3005285](https://doi.org/10.1021/bi3005285)

[6] Basuroy, K., et al. (2021). Excimer formation in pyrene crystals under pressure. *Journal of Chemical Physics*, 155, 234304. [https://doi.org/10.1063/5.0072498](https://doi.org/10.1063/5.0072498)

[7] Mazzeo, P., et al. (2024). Twisted excimers as dominant species in pyrene dimer photophysics. *ChemRxiv* (preprint).

[8] Marazzi, M., et al. (2024). TD-DFT study of pyrene excimer equilibrium geometry. *Journal of Physical Chemistry Letters*.
