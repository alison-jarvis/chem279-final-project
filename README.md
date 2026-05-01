# CNDO/2 Molecular Property Classification


## Dipole Moment Calculation

The permanent dipole can be calculated from the optimized molecular structure, as it is a property of the molecule with no external electric field applied. Therefore, a reasonable approximation for the permanent dipole can be calculated using the density matrices. For unrestricted CNDO/2 we get the that total density matrix is:

$$P_{\mu \nu} = P_{\mu \nu}^{\alpha} + P_{\mu \nu}^{\beta}$$

From the density matrix, for each atom $A$ we can calculate the net electron population on atom $A$ (which we call $N_A$) as:

$$N_A = \sum_{\mu \in A} P_{\mu \mu}$$

Then for every atom, we can calculate the net atomic charge, which is:

$$q_A = Z_A - N_A$$

Where $Z_A$ is the nuclear charge for that atom. We can then calculate the dipole moment, $\mu$, as:

$$\mathbf{\mu} = \sum_{A} q_A \mathbf{R_A}$$

Here, $\mathbf{R_A}$ is the center of each atom $A$. The dipole moment $\mathbf{\mu}$ will be a vector corresponding to the dipole in each cartesian direction. $\mu = |\mathbf{\mu}|$ represents the overall strength of the dipole. 


