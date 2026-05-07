# CNDO/2 Molecular Property Classification

## Previous Dipole Moment Calculation

The permanent dipole is a static property of the optimized molecular structure. We can use an atom-centered approximation (for now) which essentially sums the electron cloud over each atom using the CNDO/2 SCF density matrices. For unrestricted CNDO/2 we get the that total density matrix is:

$$P_{\mu \nu} = P_{\mu \nu}^{\alpha} + P_{\mu \nu}^{\beta}$$

From the density matrix, for each atom $A$ we can calculate the net electron population on atom $A$ (which we call $N_A$) as:

$$N_A = \sum_{\mu \in A} P_{\mu \mu}$$

Then for every atom, we can calculate the net atomic charge, which is:

$$q_A = Z_A - N_A$$

Where $Z_A$ is the nuclear charge for that atom. We can then calculate the dipole moment, $\mu$, as:

$$\mathbf{\mu} = \sum_{A} q_A \mathbf{R_A}$$

Here, $\mathbf{R_A}$ is the vector center of each atom $A$. The dipole moment $\mathbf{\mu}$ will be a vector corresponding to the dipole in each cartesian direction. $\mu = |\mathbf{\mu}|$ represents the overall strength of the dipole. 

## Dipole Calculation - AO Integral Method

We want to modify the hamiltonian by applying an external electric field. 

Total molecular dipole:

$$\mu_{x} = - \frac{\partial E}{\partial F_x}$$

Which is the derivative of energy with respect to a change in external electric field in a given component, $F_x$. 

We can divide this into two components, the electronic and nuclear contribution. 

$$\mathbf{\mu}_x = \mathbf{\mu}_{x}^{nuc} + \mathbf{\mu}_{x}^{elec}$$

The nuclear contribution is sum of effective nuclear charges over all atomic positions:

$$\mathbf{\mu}_{x}^{nuc} = \sum_{A} Z_A \mathbf{R_{A,x}}$$

The electronic contribution is:

$$\mathbf{\mu}_x^{elec} = - \sum_{\mu \nu} P_{\mu \nu} D^x_{\mu \nu}$$

Where $D^x_{\mu \nu}$ is the change in the CNDO/2 hamiltonian with an external field perturbation $F_x$, i.e.:

$$\frac{\partial h_{\mu \nu}}{\partial F_x} = D^x_{\mu \nu}$$

How we get the $D^x_{\mu \nu}$ terms is:

$$D^x_{\mu \nu} = \sum_{kl} d_{k \mu} d_{l \nu} N_{k \mu} N_{l \nu} D_x^{kl}$$

Where:

$$D_x^{kl} = [S_x^{{kl}} (l_{\mu, x} + 1, l_{\nu, x}) + X_{\mu} S_x^{{kl}}(l_{\mu, x}, l_{\nu, x})] S_y^{kl} S_z^{kl}$$


## Polarizability Calculation

In order to calculate the polarizability tensor, we need to modify the core hamiltonian to account for a uniform external electric field. For a full ab initio method, the change in the hamiltonian due to a uniform external electric field $\mathbf{F}$ is:

$$\hat{H}^{field} = \mathbf{F} \cdot \mathbf{r}$$

In matrix and full operator form, this would be:

$$\hat{H}_{\mu \nu}^{field} = \langle \psi_{\mu} | \mathbf{F} \cdot \mathbf{r} | \psi_{\nu} \rangle$$

Because we are using neglect of differential overlap theory though, we assume that the orbitals on different atoms don't interact significantly. Therefore, for our CNDO/2 approximation, we can simplify to:

$$\Delta H_{\mu \nu}^{field} = 0, \mu \neq \nu$$

$$\Delta H_{\mu \mu}^{field} = \mathbf{F} \cdot \mathbf{R_A}$$

So $\Delta H_{\mu \mu}^{field}$ only contains diagonal terms in the following form, where $\mu \in A$:

$$\Delta H_{\mu \mu}^{field} = F_x x_A + F_y y_A + F_z z_A$$

We have to apply this change in the hamiltonian due to the external electric field to both the fock matrix and the core hamiltonian, so that the SCF solution will produce updated density matrices and energies. The modified fock and hamiltonian matrices as a function of $\mathbf{F}$ are:

$$F_{\mu \mu}(\mathbf{F}) = F_{\mu \mu}(0) + \Delta H_{\mu \mu}^{field}$$

$$H_{\mu \mu}(\mathbf{F}) = H_{\mu \mu}(0) + \Delta H_{\mu \mu}^{field}$$

So this allows us to run SCF with an arbitrary uniform external electric field applied, which will result in modified density matrices $P_{\mu \nu}^{\alpha}$ and $P_{\mu \nu}^{\beta}$. We can use these to calculate the dipole with an applied electric field. Using this principle, we can estimate the polarizability, which is conceptually:

$$\alpha_{ij} = \left(\frac{\partial \mu_i}{\partial F_j}\right)_{\mathbf{F}=0}$$

We don't have a closed form solution for this derivative, but we can use central finite differences to approximate it, in the form:

$$\alpha_{ij} \approx \frac{\mu_i(+F_j) - \mu_i(-F_j)}{2F}$$

The indices $i$ and $j$ refer to the cartesian direction components of the field/dipole vectors. Therefore, the final result will be a $3x3$ polarizability tensor in the format:

$$\alpha = \begin{pmatrix}
\alpha_{xx} & \alpha_{xy} & \alpha_{xz} \\
\alpha_{yx} & \alpha_{yy} & \alpha_{yz} \\
\alpha_{zx} & \alpha_{zy} & \alpha_{zz} \\
\end{pmatrix}$$

We can use this tensor to calculate additional derived properties. 


## Derived Properties Calculation

The first property is the dipole magnitude, which quantifies permanent charge asymmetry. As discussed, this is:

$$|\mu| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}$$

The second property is the isotropic polarizability, which quantifies how easily a molecule is polarized by an external electric field. This is calculated as:

$$\bar{\alpha} = \frac{\alpha_{xx} + \alpha_{yy} + \alpha_{zz}}{3}$$

The third property is the ansiotropic polarizability, which quantifies how directional the response to the external electric field is. This is calculated as:

$$\Delta \alpha = \left[ \frac{(\alpha_{xx} - \alpha_{yy})^2 + (\alpha_{yy} - \alpha_{zz})^2 + (\alpha_{zz} - \alpha_{xx})^2 + 6(\alpha_{xy}^2 + \alpha_{yz}^2 + \alpha_{xz}^2)}{2}  \right]^{\frac{1}{2}}$$

