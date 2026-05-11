# Theory

In order to develop our molecular classification pipeline, it was necessary to compute two properties for a given molecule, the dipole moment and the polarizability tensor. Conceptually, the dipole moment and the polarizability are the first and second derivatives, respectively, of the energy with respect to a change in external electric field. Using these quantities we can compute several derived properties, including dipole strength, isotropic polarizability, and anisotropic polarizability. These quantities are used in the molecular classification pipeline, as they provide tangible information about molecular properties. 

## Dipole Moment

There are several reasonable accurate approximations that can be used to calculate the dipole moment with CNDO/2 theory. Here, we explore the Mulliken approximation, as well as the full dipole integral method over atomic orbitals. 

### Dipole Approximation

The Mulliken approximation is atom-centered, and assumes the electron cloud for each atom is centered on the nucleus. We can use the optimized density matrix, $P_{\mu \nu}$, from the CNDO/2 SCF cycle, to calculate this. The total density matrix is the sum of the $\alpha$ and $\beta$ contributions from unrestricted Hartree Fock. 

$$P_{\mu \nu} = P^{\alpha}_{\mu \nu} + P^{\beta}_{\mu \nu}$$

Using this matrix, we can calculate the net electron population ($N_A$) on atom $A$ as the sum of the density matrix on that atom. 

$$N_A = \sum_{\mu \in A} P_{\mu \mu}$$

We can then calculate the net atomic charge on each atom ($q_A$) as the difference between the effective nuclear charge ($Z_A$) and the net electronic population. 

$$q_A = Z_A - N_A$$

The dipole moment is a function of charge and position. Therefore, our net dipole moment ($\mu$) is the following. 

$$\mathbf{\mu} = \sum_{A} q_A \mathbf{R_A}$$

### Dipole Integral

The dipole integral method is a more robust way to calculate the dipole moment under CNDO/2 theory. We start from first principles, using the concept that the dipole is the first derivative of the energy with respect to a change in external electric field, which we call $F$. We know that the formula for the CNDO/2 energy is the following. 

$$E_{CNDO/2} = \sum_{\mu \nu} P_{\mu \nu} (h_{\mu \nu} + f_{\mu \nu}) + E_{nuc}$$

The fock matrix, $f_{\mu \nu}$, is a function of $h_{\mu \nu}$, and additional terms which do not depend on the external electric field. Therefore, we consider the fock-specific contribution to be constant in the derivative. Using this concept and the fact that the density matrix itself is independent of an external electric field assuming converged SCF, we get the energy derivative as the following. Here, we are simplifying the derivation be for only the $x$ component, but the form will be identical for the other components. 

$$\frac{\partial E_{CNDO/2}}{\partial F_x} = \sum_{\mu \nu} P_{\mu \nu} \frac{\partial h_{\mu \nu}}{\partial F_x} + \frac{\partial E_{nuc}}{\partial F_x}$$

To determine $\frac{\partial h_{\mu \nu}}{\partial F_x}$, we need to understand how the hamiltonian depends on an external electric field. Below is the general notation for the hamiltonian as a function of $F$. 

$$\hat{H}(F) = \hat{H}_0 + F_x \hat{x} + F_y \hat{y} + F_z \hat{z}$$

Again, simplifying to just the $x$ component, if we take the derivative of this with respect to $F_x$, we get the following. 

$$\frac{\partial \hat{H}}{\partial F_x} = \hat{x}$$

Therefore, the derivative of the hamiltonian is only a function of the position operator. Putting this into the atomic orbital notation used for CNDO/2, we get the following, where we refer to the $\frac{\partial h_{\mu \nu}}{\partial F_x}$ as $D^x_{\mu \nu}$ for simplicity. 

$$D^x_{\mu \nu} = \langle \chi_{\mu} | x | \chi_{\nu} \rangle$$

Here, $\chi_{\mu}$ and $\chi_{\nu}$ are the atomic orbital basis functions and $x$ represents the position operator. Given that we are specifically using the Sto-3g basis set, which is comprised of contracted gaussians, the equation we use to compute $D^x_{\mu \nu}$ is as follows. 

$$D^x_{\mu \nu} = \sum_{k} \sum_{l} d_{k \mu} d_{l \nu} N_{k \mu} N_{l \nu} D_x^{kl}$$

Here, $k$ and $l$ index the primitive gaussians in the basis set. The $d_{k \mu}$ and the $d_{l \nu}$ terms are the contraction coefficients of primitive gaussians $k$ and $l$, respectively. Likewise, $N_{k \mu}$ and $N_{l \nu}$ are the normalization coefficients for these gaussians. Finally, $D_{kl}$ represents the dipole overlap integral for two primitive gaussians. For two arbitrary primitive gaussians $A$ and $B$, the dipole integral will be computed as follows. 

$$D_x^{AB} = \int g_A(r) \cdot x \cdot g_B(r) dr$$

Using the same techniques as for our gaussian overlap integral, we can separate this into multiplicative components, to get the following. 

$$D_x^{AB} = D_x^{AB,x} S_y^{AB} S_z^{AB}$$

The $S^{AB}$ terms are our overlap integrals for $y$ and $z$. The calculation of $D_x^{AB,x}$ works out to be the following. A full derivation of this is included in the appendix. 

$$D_x^{AB,x} = S_x^{{AB}} (l_{A} + 1, l_{B}) + X_{A} S_x^{{AB}}(l_{A}, l_{B})$$

Here, $l_A$ and $l_B$ refer to the angular momentum of primitive gaussians $A$ and $B$, respectively, and $X_A$ refers to the $x$ coordinate of the center of primitive gaussian $A$. $S_x^{AB}$ is the overlap integral term for $x$. We use these principles to calculate the $D_{\mu \nu}$ terms. 

Going back to our overall calculate of the dipole moment, the nuclear energy expressed as a function of $F$ is the following. 

$$E_{nuc} = \sum_{A} Z_A (F_x R_{A,x} + F_y R_{A,y} + F_z R_{A,z})$$

Here, $Z_A$ is the effective nuclear charge for a given atom. Therefore, the derivative of the nuclear energy with respect to a change in $F_x$ works out to be the following. 

$$\frac{\partial E_{nuc}}{\partial F_x} = \sum_{A} Z_A \mathbf{R_{A,x}}$$

We can substitute this equation into our expression for a single component of the dipole, $\mu_x$, to get the following. 

$$\mu_x = -\frac{\partial E}{\partial F_x} = -\sum_{\mu \nu} P_{\mu \nu} D^x_{\mu \nu} - \sum_{A} Z_A \mathbf{R_{A,x}}$$

We evaluate this for each component, to get the total dipole vector. 

$$\vec{\mu} = \left<\mu_x, \mu_y, \mu_z \right>$$

## Polarizability

The polarizability is conceptually the second derivative of the energy with respect to the change in external electric field. We label this quantity $\alpha$, where $\alpha = \frac{\partial^2 E}{\partial F^2}$. Given that it is significantly more difficult to derive a closed-form integral solution for this, we are using the central differences theorem to approximate the polarizability. This takes on the following theoretical formula, where $i$ and $j$ both represent cartesian directions. 

$$\alpha_{ij} = \left(\frac{\partial \mu_i}{\partial F_j}\right)_{\mathbf{F}=0}$$

Using the central differences approximation, we can evaluate this as the following. 

$$\alpha_{ij} \approx \frac{\mu_i(+F_j) - \mu_i(-F_j)}{2F}$$

In order to evaluate each $\alpha_{ij}$ term, we use our existing implementation of the CNDO/2 self-consistent field method. This requires being able to compute the dipole moment as a function of the ambient field. To accomplish this, we run our SCF cycle as normal, with a perturbation added to the hamiltonian, which is a function of the electric field and dipole integral terms. This takes on the following form. 

$$\Delta h_{\mu \nu} = F_x D^x_{\mu \nu} + F_y D^y_{\mu \nu} + F_z D^z_{\mu \nu}$$

This will yield a $3x3$ matrix for each combination of $i$ and $j$, to get the following polarizability tensor. 

$$\alpha = \begin{pmatrix}
\alpha_{xx} & \alpha_{xy} & \alpha_{xz} \\
\alpha_{yx} & \alpha_{yy} & \alpha_{yz} \\
\alpha_{zx} & \alpha_{zy} & \alpha_{zz} \\
\end{pmatrix}$$

## Derived Properties

We compute three scalar derived properties from the dipole ($\vec{\mu}$) and polarizability tensor ($\alpha$). These include the dipole strength, isotropic polarizability, and anisotropic polarizability. These are calculated as follows. 

### Dipole Strength ($\mu$)

The dipole strength ($\mu$) is computed as the magnitude of the dipole vector ($\vec{\mu}$). This represents the permanent polarity of the structurally optimized molecule. 

$$\mu = |\vec{\mu}|$$

### Isotropic Polarizability ($\bar{\alpha}$)

The isotropic polarizability is computed from the polarizability tensor. It represents the induced polarity of a given molecule when an external electric field is applied. 

$$\bar{\alpha} = \frac{\alpha_{xx} + \alpha_{yy} + \alpha_{zz}}{3}$$

### Anisotropic Polarizability ($\Delta \alpha$)

The anisotropic polarizability is also computed from the polarizability tensor. It quantifies how directional the induced dipole response to the external electric field is. 

$$\Delta \alpha = \left[ \frac{(\alpha_{xx} - \alpha_{yy})^2 + (\alpha_{yy} - \alpha_{zz})^2 + (\alpha_{zz} - \alpha_{xx})^2 + 6(\alpha_{xy}^2 + \alpha_{yz}^2 + \alpha_{xz}^2)}{2}  \right]^{\frac{1}{2}}$$

## Dipole Integral Derivation

The following represents the equations for two arbitrary primitive gaussians $A$ and $B$. 

$$g_A = (x - X_A)^{l_A} (y - Y_A)^{m_A} (z - Z_A)^{n_A} e^{\alpha |r - R_A|^2}$$

$$g_B = (x - X_B)^{l_B} (y - Y_B)^{m_B} (z - Z_B)^{n_B} e^{\beta |r - R_B|^2}$$

Therefore, the entire integral for the $D_x^{AB,x}$ term would be the following. 

$$D_x^{AB,x} = \int (x-X_A)^{l_A} x (x - X_B)^{l_B} e^{-\alpha (x-X_A)^2 - \beta (x-X_B)^2} dx$$