# Lab Pack #2: Instructor Answer Key and Teaching Guide

**Lab Pack:** 2 - 3D Exploration, PES, and Orbitals
**Version:** 1.0
**Last Updated:** 2026-03-17
**Document Type:** Instructor Materials (CONFIDENTIAL -- Do Not Distribute to Students)

---

## Document Overview

This instructor key covers all 26 numbered questions (Q1.1-Q4.2) in the Lab Pack #2 student worksheet. It also provides:

- Expected outputs for every exercise (3D views, PES curves, orbital isosurfaces)
- Common student errors with remediation strategies for each section
- Performance task rubrics on 4-point scales (PES interpretation, orbital classification)
- FAQ for common 3D interaction issues
- Prior knowledge note for LO11b (graduate extension)
- Point allocation summary

For each question, the following information is provided:

- **Question Text:** Summary from the worksheet
- **Expected Answer:** The correct response with scientific values
- **Acceptable Range:** Tolerance for numerical answers
- **Common Misconceptions:** What students might get wrong
- **Pedagogical Intent:** What the question is testing (linked to LOs)
- **Grading Notes:** Criteria for full, partial, and no credit

### Scientific Reference Environment

All numerical values were computed using:
- PySCF 2.11.0
- SciPy 1.17.0
- NumPy 2.4.1
- Python 3.12
- Date verified: 2026-03-17

### Point Allocation Summary

| Section | Questions | Points | Notes |
|---------|-----------|--------|-------|
| Section 1: 3D Molecular Exploration | Q1.1-Q1.8 | 30 | LO9 |
| Section 2: PES Scanning | Q2.1-Q2.8 | 36 | LO7, LO10, LO11 |
| Section 3: Orbital Visualization | Q3.1-Q3.8 | 34 | LO8, LO12 |
| Section 4: Synthesis | Q4.1-Q4.2 | 10 | Integrative |
| **Total** | **26 questions** | **110** | |

**Note:** Point values may be scaled to 100 points by instructors using the formula: (raw score / 110) x 100.

---

## Section 1: 3D Molecular Exploration (Q1.1-Q1.8, 30 points)

### Q1.1: Predicting H2O Geometry and Bond Angle (4 points)

**Question Text:**
> Sketch the geometry of the water molecule (H2O). Label each atom. What bond angle do you predict between the two O-H bonds?

**Expected Answer:**

Sketch should show a bent molecule with oxygen at the apex and two hydrogens below, forming a V-shape. Predicted bond angle should be approximately 104.5 degrees.

**Acceptable Range:**
- Bond angle: 100-115 degrees for full credit
- 90-120 degrees for partial credit (recognizes bent geometry)

**Common Misconceptions:**
1. Students predict 180 degrees (linear molecule), forgetting lone pairs
2. Students predict exactly 109.5 degrees (tetrahedral) without considering lone pair compression
3. Students draw the molecule as linear H-O-H

**Pedagogical Intent:**
Activates prior knowledge from general chemistry (VSEPR theory). The prediction step establishes a mental model that will be compared to the 3D visualization. (LO9, POE-predict)

**Grading Notes:**
- 4 pts: Correct bent sketch with angle 100-115 degrees
- 3 pts: Correct bent sketch with angle 90-120 degrees
- 2 pts: Bent sketch but angle significantly wrong (e.g., 60 or 150 degrees)
- 1 pt: Linear sketch but recognizes uncertainty
- 0 pts: No sketch or no angle prediction

---

### Q1.2: Comparing 3D View to Prediction (3 points)

**Question Text:**
> Does the 3D geometry match your prediction from Q1.1? What is the approximate bond angle shown? If your prediction was different, what assumption led to the discrepancy?

**Expected Answer:**

The 3D view should confirm a bent geometry. The bond angle in the preset is approximately 104.5 degrees (depends on exact geometry in the IQCP preset; the preset uses atom coordinates O at origin, H at (0, +/-0.757, 0.587) in Angstrom, giving a bond angle of approximately 104.5 degrees).

Students who predicted 109.5 degrees should note that lone pair repulsion compresses the angle below the ideal tetrahedral value. Students who predicted linearity should recognize that two lone pairs on oxygen force the molecule into a bent geometry.

**Acceptable Range:**
- Bond angle observation: 100-110 degrees
- Accept qualitative "about 105 degrees" or "slightly less than tetrahedral"

**Common Misconceptions:**
1. Students may not know how to estimate angle from 3D view -- remind them to rotate the molecule
2. Students may report an angle far from 104.5 degrees due to perspective distortion

**Pedagogical Intent:**
POE-observe/explain stage. Confronts initial prediction with visual evidence. Students who were wrong learn more than those who were right. (LO9)

**Grading Notes:**
- 3 pts: Notes agreement/disagreement with prediction, gives reasonable angle, explains any discrepancy
- 2 pts: Notes agreement/disagreement but weak explanation
- 1 pt: Records angle but no reflection on prediction
- 0 pts: No answer

---

### Q1.3: Predicting Overlap Matrix Elements (4 points)

**Question Text:**
> Before examining the matrix: which pairs of atoms do you expect to have the largest overlap? Which pairs should have the smallest?

**Expected Answer:**

- **Largest overlap:** O-H bonded pairs (basis functions on oxygen and a bonded hydrogen), because bonded atoms are close together and their orbitals overlap significantly.
- **Smallest overlap:** H-H pair (the two hydrogens), because they are the farthest apart and not directly bonded. Also, the O 1s core orbital has very small overlap with hydrogen 1s because it is tightly contracted around the oxygen nucleus.

**PySCF Reference Values:**
- S(O 2s, H1 1s) = 0.4743 (largest off-diagonal element)
- S(O 2py, H1 1s) = 0.3107
- S(H1, H2) = 0.2517
- S(O 1s, H1 1s) = 0.0539 (smallest nonzero off-diagonal)

**Acceptable Range:**
- Must identify bonded O-H pairs as having largest overlap
- Must identify either H-H or O(1s)-H as smallest

**Common Misconceptions:**
1. Students may think all oxygen-hydrogen overlaps are equally large (they vary by shell type)
2. Students may expect H-H overlap to be zero (it is small but nonzero)

**Pedagogical Intent:**
Predict stage -- forces students to reason about overlap before seeing numbers. Connects spatial proximity to the mathematical overlap integral. (LO9)

**Grading Notes:**
- 4 pts: Correctly identifies O-H as largest AND H-H or O(1s)-H as smallest with spatial reasoning
- 3 pts: Correctly identifies one of largest/smallest with reasoning
- 2 pts: Reasonable prediction but incorrect assignment
- 1 pt: Attempt with minimal reasoning
- 0 pts: No prediction

---

### Q1.4: Reading Overlap Matrix Values (4 points)

**Question Text:**
> What is the approximate value of the overlap between: an O basis function and a bonded H (e.g., S_2,6)? The two hydrogens (S_6,7)? The diagonal elements?

**Expected Answer:**

| Matrix Element | Physical Meaning | PySCF Value |
|----------------|-----------------|-------------|
| S(2,6) = S(O 2s, H1 1s) | O valence -- bonded H | 0.4743 |
| S(6,7) = S(H1 1s, H2 1s) | H -- H nonbonded | 0.2517 |
| S(i,i) (any diagonal) | Self-overlap | 1.0000 |

**Acceptable Range:**
- S(2,6): 0.40-0.55
- S(6,7): 0.20-0.30
- Diagonal: Must be exactly 1.0

**Common Misconceptions:**
1. Students may report wrong matrix element (off by one) due to indexing confusion
2. Students may not realize basis function indexing: functions 1-5 are on O, 6 and 7 on H
3. Students may confuse overlap matrix with density matrix

**Pedagogical Intent:**
Observe stage -- reading numerical matrix data and connecting indices to atoms. Builds the symbolic-spatial bridge central to LO9.

**Grading Notes:**
- 4 pts: All three values reported within acceptable range
- 3 pts: Two of three correct
- 2 pts: One correct
- 1 pt: Attempt but values incorrect or confused with another matrix
- 0 pts: No values reported

---

### Q1.5: Why Diagonal Elements Equal 1.0 (4 points)

**Question Text:**
> Why are the diagonal elements of the overlap matrix all equal to 1.0? What would it mean physically if a diagonal element were not 1.0?

**Expected Answer:**

The diagonal element S(i,i) is the overlap of basis function i with itself:

S(i,i) = integral of |phi_i(r)|^2 dr

This equals 1.0 because basis functions are normalized -- the integral of the probability density over all space is 1. If a diagonal element were not 1.0, it would mean the basis function is not normalized (the "total probability" of finding an electron in that orbital would not sum to 1).

**Acceptable Range:**
- Must mention normalization
- Must connect to the integral of |phi|^2

**Common Misconceptions:**
1. Students may say "because it's a matrix and matrices have 1s on the diagonal" (confusing with identity matrix)
2. Students may not connect normalization to probability interpretation
3. Students may think overlap matrix is always the identity matrix

**Pedagogical Intent:**
Tests understanding of normalization as a physical constraint. Connects mathematical property (S_ii = 1) to the physical requirement that wavefunctions be normalizable. (LO9)

**Grading Notes:**
- 4 pts: Mentions normalization AND connects to integral of |phi|^2 or probability
- 3 pts: Mentions normalization but incomplete physical interpretation
- 2 pts: States "functions are normalized" without further explanation
- 1 pt: Incorrect reasoning but shows attempt
- 0 pts: No answer or "because the matrix is symmetric"

---

### Q1.6: Predicting H2 Off-Diagonal Overlap (3 points)

**Question Text:**
> The H2 overlap matrix is 2x2. Before looking: predict the approximate value of S_1,2. Will it be closer to 0.0, 0.5, or 1.0?

**Expected Answer:**

S_1,2 should be closer to 0.5 (actual value: 0.6593 at R = 1.4 bohr). The two 1s orbitals on the hydrogen atoms at the equilibrium bond length have significant but not complete overlap.

**PySCF Reference:** S_1,2 = 0.659318 at R = 1.4 bohr

**Acceptable Range:**
- Predicting "closer to 0.5" earns full credit
- Predicting "between 0.5 and 1.0" is also acceptable
- Predicting "closer to 0.0" is incorrect (atoms are bonded, substantial overlap)

**Common Misconceptions:**
1. Students predict 1.0 (confusing bonded atoms with identical orbitals)
2. Students predict 0.0 (thinking different atoms means no overlap)

**Pedagogical Intent:**
POE-predict for the simplest molecule. Establishes quantitative intuition about overlap magnitude for bonded atoms. (LO9)

**Grading Notes:**
- 3 pts: Predicts value near 0.5-0.7 with reasoning about bonded atoms
- 2 pts: Predicts a reasonable positive value with some reasoning
- 1 pt: Makes a prediction but reasoning is weak or absent
- 0 pts: No prediction

---

### Q1.7: Comparing H2 and H2O Hydrogen Overlaps (4 points)

**Question Text:**
> What is the actual value of S_1,2 for H2? How does this compare to S_6,7 in water? Explain the difference in terms of interatomic distances.

**Expected Answer:**

- H2: S_1,2 = 0.6593 (hydrogens are bonded, R approximately 1.4 bohr = 0.74 Angstrom)
- H2O: S_6,7 = 0.2517 (hydrogens are not bonded to each other, separated by approximately 1.51 Angstrom)

The H2 overlap is much larger because the hydrogen atoms in H2 are directly bonded and therefore much closer together. In H2O, the two hydrogens are separated by a greater distance (they are each bonded to oxygen, but not to each other), so their orbital overlap is smaller.

**Acceptable Range:**
- H2 S_1,2: 0.60-0.70
- Must correctly identify that H2 overlap > H2O H-H overlap
- Must connect the difference to interatomic distance

**Common Misconceptions:**
1. Students may think H-H overlap is the same in both molecules
2. Students may attribute the difference to "different molecules" without citing distance

**Pedagogical Intent:**
Central to LO9 -- demonstrates that overlap depends on distance. The comparison between bonded (H2) and non-bonded (H-H in H2O) hydrogen pairs makes the distance-overlap relationship concrete.

**Grading Notes:**
- 4 pts: Correct values for both, identifies H2 > H2O, explains via distance
- 3 pts: Correct comparison direction with distance explanation, approximate values
- 2 pts: Correct comparison but weak or missing distance explanation
- 1 pt: Values reported but no comparison or explanation
- 0 pts: No answer

---

### Q1.8: Predicting Effect of Bond Stretching on Overlap (4 points)

**Question Text:**
> If you were to increase the H-H bond length in H2, would S_1,2 increase or decrease? Explain using both the 3D view and the mathematical meaning of the overlap integral.

**Expected Answer:**

S_1,2 would **decrease** as R increases.

**3D reasoning:** As the atoms move apart, the electron clouds centered on each atom have less spatial region in common, so the overlap decreases.

**Mathematical reasoning:** The overlap integral S_1,2 = integral of phi_1(r) * phi_2(r) dr measures the product of two Gaussian-like functions centered at different points. As the centers move apart, the product of the two functions becomes smaller everywhere, reducing the integral.

**Limiting behavior:** As R approaches infinity, S_1,2 approaches 0 (no overlap). As R approaches 0, S_1,2 approaches 1 (functions become identical).

**Acceptable Range:**
- Must state "decrease"
- Must provide at least one line of reasoning (3D or mathematical)

**Common Misconceptions:**
1. Students may think overlap increases with distance (confusing with volume between atoms)
2. Students may think overlap is binary (either "overlapping" or "not") rather than continuous

**Pedagogical Intent:**
Tests predictive reasoning about the overlap integral. Synthesizes spatial and mathematical understanding. Prepares for PES section where changing geometry changes all matrix elements. (LO9)

**Grading Notes:**
- 4 pts: Correct prediction (decrease) with both 3D and mathematical reasoning
- 3 pts: Correct prediction with one type of reasoning
- 2 pts: Correct prediction with minimal reasoning
- 1 pt: Incorrect prediction but shows some reasoning
- 0 pts: No answer or "increase" with no reasoning

---

## Section 2: PES Scanning (Q2.1-Q2.8, 36 points)

### Q2.1: Predicting PES Curve Shape (4 points)

**Question Text:**
> Sketch what you think the energy vs. bond length curve looks like for H2.

**Expected Answer:**

The sketch should show:
- A **repulsive wall** at short R (energy rises steeply as R approaches 0)
- A **minimum** (well) at intermediate R (approximately 1.3-1.5 bohr)
- A **plateau or gradual approach** to a dissociation limit at large R
- The y-axis should decrease (become more negative) going from short R to the minimum

The overall shape should resemble a Morse potential or Lennard-Jones curve.

**PySCF Reference Points:**
- R = 0.5 bohr: E = -0.4033 Ha (repulsive wall)
- R = 1.35 bohr: E = -1.1175 Ha (near equilibrium)
- R = 5.0 bohr: E = -0.6864 Ha (dissociation region)

**Acceptable Range:**
- Sketch must show a minimum between 0.5 and 3.0 bohr
- Sketch must show energy rising at short R and leveling off at large R

**Common Misconceptions:**
1. Students may draw a symmetric parabola (harmonic approximation)
2. Students may draw the curve with the minimum at R = 0
3. Students may not label axes correctly

**Pedagogical Intent:**
POE-predict stage. Forces students to commit to a mental model of the PES before seeing the computed result. Wrong predictions are pedagogically valuable. (LO7)

**Grading Notes:**
- 4 pts: Correct qualitative shape with repulsive wall, minimum, and dissociation plateau
- 3 pts: Correct minimum with either repulsive wall or plateau
- 2 pts: Shows a minimum but overall shape incorrect
- 1 pt: Some attempt at a sketch with labeled axes
- 0 pts: No sketch

---

### Q2.2: Comparing Computed PES to Prediction (4 points)

**Question Text:**
> Compare the computed PES curve to your sketch. What features did you predict correctly? What surprised you?

**Expected Answer:**

Students should record numerical values from the scan:
- Equilibrium bond length: approximately 1.35-1.45 bohr (parabolic interpolation gives R_eq = 1.346 bohr)
- Equilibrium energy: approximately -1.117 Ha
- Energy at R = 0.5 bohr: approximately -0.40 Ha (repulsive)
- Energy at R = 5.0 bohr: approximately -0.69 Ha (dissociation region)

Students should note what they predicted correctly (usually the minimum) and what surprised them (commonly the asymmetry of the curve, the steepness of the repulsive wall, or the fact that energy at large R does not go to zero).

**Acceptable Range:**
- R_eq: 1.2-1.5 bohr
- E_eq: -1.05 to -1.15 Ha (depends on scan grid spacing)

**Common Misconceptions:**
1. Students may think the minimum should be at exactly 0.74 Angstrom (the experimental value) -- remind them that STO-3G gives a different optimum
2. Students may be surprised that energy is negative (expecting positive energies)

**Pedagogical Intent:**
POE-observe/explain. Self-assessment of prediction quality. Students who predicted wrong learn the most. (LO7)

**Grading Notes:**
- 4 pts: Records numerical values AND reflects on prediction accuracy
- 3 pts: Records values with brief reflection
- 2 pts: Records values but no reflection
- 1 pt: Qualitative comparison only, no values
- 0 pts: No answer

---

### Q2.3: Why Does the Energy Minimum Exist? (5 points)

**Question Text:**
> Why does the minimum exist? Consider two competing effects: nuclear-nuclear repulsion at short R and electron delocalization at intermediate R.

**Expected Answer:**

The energy minimum exists because of two competing effects:

1. **At very short R:** Nuclear-nuclear repulsion (V_nn = Z_A * Z_B / R) dominates. The positive nuclei repel each other, and this repulsion increases rapidly as R decreases (1/R dependence). Additionally, electron-electron repulsion increases as the electron clouds are forced together.

2. **At intermediate R (near equilibrium):** Electrons can delocalize over both nuclei, lowering the kinetic energy (electrons are confined to a larger effective region) and increasing the attractive electron-nuclear interaction (each electron "sees" both nuclei). This stabilization creates a net energy lowering compared to separated atoms.

3. **The balance:** The minimum occurs where the derivative of the total energy with respect to R is zero -- the point where the stabilizing effects of electron delocalization exactly balance the onset of nuclear repulsion.

**Acceptable Range:**
- Must mention nuclear repulsion at short R
- Must mention some form of electron stabilization (delocalization, orbital overlap, or lowered kinetic energy)
- Must convey that the minimum is a balance point

**Common Misconceptions:**
1. Students may only cite "attraction" without specifying what attracts what
2. Students may invoke "electron sharing" without explaining why sharing lowers energy
3. Students may think the minimum exists purely due to Coulomb attraction between electron and nucleus (missing the kinetic energy contribution)

**Pedagogical Intent:**
Tests conceptual understanding of the physical origin of chemical bonds. This is a core question for LO7. (LO7, Analyze)

**Grading Notes:**
- 5 pts: Both competing effects clearly explained, balance concept present
- 4 pts: Both effects mentioned, balance concept implied
- 3 pts: One effect clearly explained, other mentioned
- 2 pts: One effect explained
- 1 pt: Vague or mostly incorrect
- 0 pts: No answer

---

### Q2.4: Cause of the Repulsive Wall (4 points)

**Question Text:**
> The repulsive wall at short R rises steeply. Is this primarily due to (a) electron-electron repulsion, (b) nuclear-nuclear repulsion, or (c) both?

**Expected Answer:**

**(c) Both**, but nuclear-nuclear repulsion is the dominant contributor at very short R.

**Detailed explanation:** As R decreases:
- Nuclear-nuclear repulsion V_nn = Z_A * Z_B / R increases as 1/R and diverges as R approaches 0.
- Electron-electron repulsion also increases as the electron clouds are compressed into a smaller region.
- The kinetic energy also increases (via the uncertainty principle: confining electrons to a smaller region increases kinetic energy).

At the steepest part of the repulsive wall (very small R), V_nn dominates because it diverges as 1/R. However, the onset of the wall (near equilibrium) involves all three effects.

**Acceptable Range:**
- "(c) both" with some justification earns full credit
- "(b) nuclear" alone earns partial credit if justified
- "(a) electron" alone is insufficient

**Common Misconceptions:**
1. Students may attribute it entirely to electron repulsion (Pauli exclusion) -- this contributes but V_nn dominates at very short R
2. Students may forget that kinetic energy also increases

**Pedagogical Intent:**
Forces students to distinguish between different types of repulsion. Deepens understanding of what drives the repulsive wall. (LO7)

**Grading Notes:**
- 4 pts: Identifies both with nuclear dominance at short R, provides justification
- 3 pts: Identifies both but incomplete justification
- 2 pts: Identifies one correctly with justification
- 1 pt: Correct answer choice but no justification
- 0 pts: Incorrect or no answer

---

### Q2.5: Dissociation Limit Comparison (5 points)

**Question Text:**
> The energy of two isolated H atoms is 2 x E(H) = -0.9332 Ha. Look at the energy at R = 5.0 bohr. Is it close to -0.9332 Ha?

**Expected Answer:**

- Energy at R = 5.0 bohr (from PES scan): approximately **-0.6864 Ha**
- Expected for two isolated H atoms: **-0.9332 Ha**
- Difference: approximately **+0.247 Ha** (the RHF energy is too high by about 0.25 Ha)

The computed energy at large R is **not** close to the correct dissociation limit. It is substantially higher (less negative) than the energy of two isolated hydrogen atoms.

**PySCF Reference Values:**
- E(H atom, UHF/STO-3G) = -0.466582 Ha
- 2 * E(H) = -0.933164 Ha
- E(H2, RHF/STO-3G, R=5.0) = -0.686416 Ha
- Error: 0.247 Ha (too high)

**Acceptable Range:**
- Must record the energy at R = 5.0 bohr within 0.02 Ha of -0.69 Ha
- Must identify that it is NOT close to -0.9332 Ha
- Must note the discrepancy direction (computed energy too high / less negative)

**Common Misconceptions:**
1. Students may think "close enough" because both are negative
2. Students may not realize that 0.25 Ha is a large error (it is approximately 157 kcal/mol or 6.8 eV)
3. Students may attribute the error to basis set quality rather than the RHF method

**Pedagogical Intent:**
Sets up the RHF dissociation failure discussion. Students discover the problem empirically before learning the explanation. (LO11a)

**Grading Notes:**
- 5 pts: Correct energy recorded, discrepancy identified and quantified, direction noted
- 4 pts: Correct energy, discrepancy identified but not quantified
- 3 pts: Correct energy recorded, notes "not close" without quantification
- 2 pts: Energy recorded but no comparison
- 1 pt: Attempt but wrong energy or wrong comparison
- 0 pts: No answer

---

### Q2.6: RHF Dissociation Failure (5 points)

**Question Text:**
> The energy at large R does not approach the correct limit. What physical process does RHF fail to describe correctly at large R? Why might requiring alpha and beta electrons to share the same spatial orbitals cause problems when the bond is stretched?

**Expected Answer:**

RHF fails at large R because it uses **restricted orbitals** -- alpha and beta electrons are forced to occupy the same spatial orbital. At the equilibrium geometry, this is a good approximation because both electrons are shared between the nuclei. But as the bond stretches:

1. The correct physical picture at large R is two **neutral hydrogen atoms**, each with one electron. This requires different spatial distributions for the alpha and beta electrons -- one electron localized on atom A and the other on atom B.

2. RHF cannot describe this because it constrains both electrons to the same spatial orbital (sigma_g). The resulting wavefunction is a mixture of correct (H + H) and incorrect (H+ + H-) ionic configurations with equal weight. The ionic configuration (proton + hydride) is much higher in energy, raising the computed dissociation limit above the correct value.

3. In mathematical terms, the RHF wavefunction sigma_g(1)sigma_g(2) expands to [1s_A(1) + 1s_B(1)][1s_A(2) + 1s_B(2)] = [1s_A(1)1s_B(2) + 1s_B(1)1s_A(2)] + [1s_A(1)1s_A(2) + 1s_B(1)1s_B(2)]. The first bracket is the correct covalent component; the second bracket is the incorrect ionic component.

**Graduate Extension (LO11b):**
The correct wavefunction requires at least two determinants (or configurations): the sigma_g^2 and sigma_u*^2 configurations must be mixed. This is a multi-reference problem requiring methods like CASSCF, MCSCF, or UHF (which breaks spatial symmetry to capture the correct physics).

**Acceptable Range:**
- Must mention restriction of alpha/beta to same orbital
- Must connect this to incorrect behavior at large R
- Graduate extension: must mention multiple determinants or configurations

**Common Misconceptions:**
1. Students may blame basis set quality (STO-3G is small, but the dissociation error persists with any basis for RHF)
2. Students may confuse electron correlation with this issue (correlation exists at all R, but the dissociation failure is specifically about the single-determinant constraint)
3. Students may not understand why "same spatial orbital" is a problem

**Pedagogical Intent:**
Core question for LO11a (undergrad) and LO11b (graduate). Tests understanding of RHF's fundamental limitation. (LO11a/LO11b)

**Grading Notes:**
- 5 pts: Identifies restriction as the problem, explains covalent/ionic mixing or equivalent, connects to energy error
- 4 pts: Identifies restriction, partial explanation of the consequence
- 3 pts: Identifies restriction but minimal explanation of why it causes problems
- 2 pts: Mentions RHF limitation without connecting to orbital restriction
- 1 pt: Vague response about "method limitations"
- 0 pts: No answer or blames only basis set

**Graduate Extension Grading (bonus or separate):**
- Full credit: Explains single-determinant limitation and names multi-reference methods
- Partial credit: Mentions need for multiple configurations without naming methods

---

### Q2.7: Basis Set Comparison -- Variational Principle (5 points)

**Question Text:**
> Compare the equilibrium energy of H2 with STO-3G vs. 6-31G. Which gives lower energy? Why does a larger basis set always give a lower energy?

**Expected Answer:**

- **STO-3G:** E_eq approximately -1.1175 Ha (2 basis functions)
- **6-31G:** E_eq approximately -1.1268 Ha (4 basis functions)

The 6-31G energy is **lower** (more negative) because it has more basis functions, providing more variational freedom.

**Why:** The variational principle states that any trial wavefunction will yield an energy greater than or equal to the true ground-state energy:

E_trial >= E_exact

A larger basis set spans a larger subspace of the Hilbert space, allowing the optimization to find a lower-energy wavefunction. Since STO-3G is a subset of the functions that could be represented in 6-31G, the 6-31G result must be at least as low.

**PySCF Reference Values:**
- H2 STO-3G equilibrium energy: -1.117501 Ha (at R = 1.35 bohr)
- H2 6-31G equilibrium energy: -1.126828 Ha (at R = 1.38 bohr)
- Difference: 0.009 Ha (6-31G is lower by about 6 kcal/mol)

**Acceptable Range:**
- Must correctly identify 6-31G as giving lower energy
- Must invoke variational principle or equivalent reasoning

**Common Misconceptions:**
1. Students may think "lower energy" means STO-3G because the number is closer to zero
2. Students may attribute the difference to "better parameters" rather than more variational freedom
3. Students may think larger basis sets are always "more accurate" in all respects (they are for energy, but not necessarily for other properties)

**Pedagogical Intent:**
Tests understanding of the variational principle, which was introduced in Lab Pack #1. Extends the concept to basis set comparison. (LO10)

**Grading Notes:**
- 5 pts: Correct identification + variational principle explanation + approximate energy values
- 4 pts: Correct identification + variational principle, no numerical values
- 3 pts: Correct identification with partial explanation
- 2 pts: Correct identification but no explanation
- 1 pt: Incorrect identification or confused reasoning
- 0 pts: No answer

---

### Q2.8: Evaluating "Bigger Basis Always Better" (4 points)

**Question Text:**
> A student claims: "A bigger basis set always gives a better answer, so we should always use the biggest possible basis." Do you agree? What practical considerations might limit this approach?

**Expected Answer:**

**Partially agree, partially disagree.** A larger basis set does always give a lower energy (variational principle), and the energy is closer to the basis set limit. However, practical considerations include:

1. **Computational cost:** The number of two-electron integrals scales as O(N^4) where N is the number of basis functions. Doubling the basis set increases computation by roughly 16x. For large molecules, this becomes prohibitive.

2. **Diminishing returns:** The energy improvement from STO-3G to 6-31G is much larger than from 6-31G to cc-pVTZ. Each additional basis function contributes less to the total energy lowering.

3. **Other sources of error:** Using RHF with a very large basis set does not fix the method's inherent limitations (e.g., lack of electron correlation, RHF dissociation failure). The basis set limit of RHF is still above the true energy.

4. **"Better" depends on the property:** A larger basis set gives a lower energy but may not improve every molecular property equally. For some properties, a carefully chosen small basis can perform well.

**Acceptable Range:**
- Must disagree at least partially
- Must cite at least one practical consideration (cost or diminishing returns)

**Common Misconceptions:**
1. Students may fully agree without recognizing cost scaling
2. Students may think "better basis = correct answer" without understanding that method limitations persist

**Pedagogical Intent:**
Evaluate-level question (LO10). Requires weighing competing factors rather than applying a rule. (LO10, Evaluate)

**Grading Notes:**
- 4 pts: Partial disagreement with 2+ practical considerations, well-reasoned
- 3 pts: Partial disagreement with 1 practical consideration
- 2 pts: Agrees but mentions cost as a caveat
- 1 pt: Simple agree/disagree with no reasoning
- 0 pts: No answer

---

## Section 3: Orbital Visualization (Q3.1-Q3.8, 34 points)

### Q3.1: Predicting MO 1 Shape (3 points)

**Question Text:**
> H2O MO 1 is built primarily from the oxygen 1s atomic orbital. Predict: will MO 1 look like (a) a small sphere centered on oxygen, (b) a shape spread across the entire molecule, or (c) two lobes pointing along the O-H bonds?

**Expected Answer:**

**(a) A small sphere centered on oxygen.**

**Reasoning:** The oxygen 1s orbital is a core orbital. It is tightly bound (orbital energy approximately -20.24 Ha) and very compact. It does not significantly participate in bonding. The MO coefficient analysis confirms that MO 1 is 99.4% oxygen 1s character:

| Basis Function | Coefficient |
|----------------|-------------|
| O 1s | 0.9941 |
| O 2s | < 0.05 |
| H 1s | < 0.01 |

**Common Misconceptions:**
1. Students may predict (b) or (c) because they expect all molecular orbitals to involve bonding
2. Students may not distinguish core from valence orbitals

**Pedagogical Intent:**
POE-predict. Introduces the distinction between core and valence MOs. (LO8)

**Grading Notes:**
- 3 pts: Correct prediction (a) with reasoning about core orbital
- 2 pts: Correct prediction (a) with minimal reasoning
- 1 pt: Incorrect prediction with some reasoning
- 0 pts: No prediction

---

### Q3.2: Describing Core Orbital Shape (4 points)

**Question Text:**
> Describe the shape of MO 1. Is it centered on one atom or spread across the molecule? Does your observation match your prediction? Why is this orbital almost unaffected by bonding?

**Expected Answer:**

MO 1 is a **compact sphere centered on the oxygen atom**. It is not spread across the molecule. This matches prediction (a).

**Why unaffected by bonding:** The 1s core orbital on oxygen is at a much lower energy (-20.24 Ha) than the valence orbitals (-1.27 Ha for MO 2). The energy gap between the core and valence levels is enormous (approximately 19 Ha = 517 eV). For bonding to significantly mix two orbitals, they must be close in energy. The hydrogen 1s orbital (-0.5 Ha for a free H atom) is far too high in energy to interact meaningfully with the oxygen 1s.

**Acceptable Range:**
- Must describe compact/spherical shape centered on O
- Must explain why bonding does not affect it (energy gap or "core orbital" reasoning)

**Common Misconceptions:**
1. Students may think all orbitals in a molecule must be "molecular" (spread across atoms)
2. Students may not understand why energy matching matters for orbital mixing

**Pedagogical Intent:**
POE-observe/explain. Establishes that core orbitals are essentially atomic in character. (LO8)

**Grading Notes:**
- 4 pts: Correct shape description + comparison to prediction + energy-based explanation
- 3 pts: Correct shape + comparison, partial explanation
- 2 pts: Correct shape, no comparison to prediction
- 1 pt: Vague description
- 0 pts: No answer

---

### Q3.3: HOMO Shape and Classification (5 points)

**Question Text:**
> Describe the shape of the H2O HOMO. How is it different from MO 1? Classify the HOMO as bonding, nonbonding, or antibonding.

**Expected Answer:**

The H2O HOMO (MO 5, energy = -0.391 Ha) is the oxygen 2px lone pair orbital. It is a **p-orbital shape** -- two lobes on opposite sides of the oxygen atom, oriented perpendicular to the molecular plane.

**Key observations:**
- Unlike MO 1 (spherical, compact), the HOMO has two lobes with a nodal plane
- The HOMO is centered entirely on oxygen with virtually no hydrogen contribution (coefficient: O 2px = 1.0000)
- It extends further from the nucleus than MO 1

**Classification: Nonbonding.** The HOMO does not have significant electron density between the O and H atoms. It is essentially an atomic p-orbital on oxygen that does not participate in bonding. This is a lone pair.

**PySCF MO Coefficients for MO 5 (HOMO):**
| Basis Function | Coefficient |
|----------------|-------------|
| O 2px | 1.0000 |
| All others | < 0.01 |

**Acceptable Range:**
- Must describe two-lobed or p-orbital shape
- Must correctly classify as nonbonding (not bonding or antibonding)
- Must note difference from MO 1

**Common Misconceptions:**
1. Students may classify it as antibonding because it has two lobes (lobes do not automatically mean antibonding)
2. Students may confuse the HOMO with a bonding orbital because they expect the highest occupied orbital to be the "most bonding"
3. Students may call it "bonding" because it is occupied

**Pedagogical Intent:**
Core question for LO8. Tests ability to classify orbitals from shape. The HOMO being nonbonding is often surprising. (LO8, Analyze)

**Grading Notes:**
- 5 pts: Correct shape description + correct classification (nonbonding) + explanation referencing lone pair or O-only character
- 4 pts: Correct classification with partial explanation
- 3 pts: Correct shape description but incorrect classification
- 2 pts: Identifies difference from MO 1 but weak classification
- 1 pt: Vague response
- 0 pts: No answer

---

### Q3.4: Predicting Antibonding Orbital Shape (3 points)

**Question Text:**
> MO 1 is bonding (sigma_g) and MO 2 is antibonding (sigma_u*). Before viewing MO 2: predict how the antibonding orbital shape will differ from the bonding orbital. Where do you expect a node?

**Expected Answer:**

The antibonding orbital (sigma_u*) should have:
- **Two separate lobes**, one on each hydrogen atom (in contrast to the single merged lobe of sigma_g)
- A **nodal plane perpendicular to the bond axis**, located at the midpoint between the two atoms (where psi = 0)
- **Opposite signs** of the wavefunction on the two lobes (one positive, one negative)

The node is expected **between the two nuclei**, at the center of the bond.

**PySCF MO Coefficients for H2 MO 2 (sigma_u*):**
| Basis Function | Coefficient |
|----------------|-------------|
| H1 1s | +1.2115 |
| H2 1s | -1.2115 |

The equal-magnitude, opposite-sign coefficients confirm the antisymmetric (antibonding) character.

**Acceptable Range:**
- Must predict a node between the atoms
- Must predict two separate lobes or opposite signs

**Common Misconceptions:**
1. Students may expect the antibonding orbital to simply be "smaller" rather than having different topology
2. Students may place the node at one of the atoms rather than between them

**Pedagogical Intent:**
POE-predict. Tests knowledge of bonding vs. antibonding orbital topology before visualization. (LO8)

**Grading Notes:**
- 3 pts: Predicts node between atoms AND two lobes with opposite signs
- 2 pts: Predicts node OR two lobes, but not both
- 1 pt: Some prediction with partial correctness
- 0 pts: No prediction

---

### Q3.5: Sigma_g vs. Sigma_u* Comparison Table (5 points)

**Question Text:**
> Compare the shapes of MO 1 (sigma_g) and MO 2 (sigma_u*) for H2. Fill in the comparison table.

**Expected Answer:**

| Feature | MO 1 (sigma_g) | MO 2 (sigma_u*) |
|---------|----------------|------------------|
| Number of lobes | 1 (single merged lobe) | 2 (one per atom) |
| Electron density between nuclei | **High** | **Low** |
| Node between atoms? | **No** | **Yes** |
| Character | **Bonding** | **Antibonding** |

**Additional observations students may note:**
- The sigma_g orbital has both coefficients positive (+0.549, +0.549), meaning constructive interference between the atomic orbitals
- The sigma_u* orbital has opposite signs (+1.211, -1.211), meaning destructive interference
- The rendering shows positive lobes as solid and negative lobes as translucent

**Acceptable Range:**
- All four rows must be correctly filled
- Accept "merged"/"continuous" for 1 lobe in sigma_g

**Common Misconceptions:**
1. Students may say sigma_g has 2 lobes (confusing with a p-orbital)
2. Students may mark "Low" density between nuclei for both (not examining sigma_g carefully)
3. Students may confuse node between atoms with nodal planes of p-orbitals

**Pedagogical Intent:**
Structured comparison builds systematic classification skills. The table format ensures students examine each feature. (LO8)

**Grading Notes:**
- 5 pts: All 4 rows correct
- 4 pts: 3 of 4 rows correct
- 3 pts: 2 of 4 rows correct
- 2 pts: 1 row correct with attempt at others
- 1 pt: Attempt but mostly incorrect
- 0 pts: No attempt

---

### Q3.6: Predicting Isovalue Effect (3 points)

**Question Text:**
> As you decrease the isovalue, do you predict the isosurface will (a) shrink, (b) expand, or (c) stay the same size?

**Expected Answer:**

**(b) Expand.**

**Reasoning:** The isovalue is the threshold value of |psi(r)| at which the surface is drawn. Decreasing the isovalue means we are drawing the surface at points where |psi| is smaller. Since |psi| decreases with distance from the nucleus, a smaller isovalue threshold means the surface extends further out, encompassing a larger region of space.

**Analogy:** Think of topographic contour lines. A contour line at elevation 100m encloses a smaller area than a contour line at elevation 50m (which extends further from the peak).

**Common Misconceptions:**
1. Students predict (a) shrink, confusing "smaller value" with "smaller surface"
2. Students predict (c) stay the same, not understanding what the isovalue controls

**Pedagogical Intent:**
POE-predict for the isovalue concept. Targets misconception about orbital boundaries. (LO12)

**Grading Notes:**
- 3 pts: Correct prediction (expand) with reasoning about threshold meaning
- 2 pts: Correct prediction with minimal reasoning
- 1 pt: Incorrect prediction with some reasoning about isovalues
- 0 pts: No prediction or prediction with no reasoning

---

### Q3.7: Isovalue Observations and Orbital Boundaries (5 points)

**Question Text:**
> What happened when you decreased the isovalue to 0.01? When you increased it to 0.08? Does the orbital have a definite, physical edge?

**Expected Answer:**

**Observations:**
- **Isovalue = 0.01:** The isosurface **expanded** significantly, enclosing a much larger region of space. The surface extends further from the nuclei.
- **Isovalue = 0.08:** The isosurface **shrank** to a much smaller region close to the nuclei. It may almost disappear for diffuse parts of the orbital.
- **Isovalue = 0.03 (default):** Intermediate size, showing the "conventional" orbital shape.

**Key insight:** The orbital does **NOT** have a definite, physical edge. The boundary the student sees on screen is simply a chosen threshold -- a mathematical contour surface. The orbital wavefunction extends to infinity in all directions (it just gets exponentially smaller). Changing the isovalue changes where we "draw the line," but this line is our choice, not a physical boundary.

**Misconception correction:** "The orbital is the colored surface I see on screen. Electrons live inside that surface and cannot exist outside it." -- This is **wrong** because:
- The surface is a chosen threshold, not a physical barrier
- Electron probability density exists on both sides of the surface
- There is no quantum mechanical equivalent of a "wall" at the isosurface
- The electron has a nonzero probability of being found at any finite distance from the nucleus

**Acceptable Range:**
- Must describe expansion at low isovalue and contraction at high isovalue
- Must conclude that orbital has no definite edge
- Must explain the isovalue as a chosen threshold

**Common Misconceptions:**
1. "Orbitals have definite edges" -- the most targeted misconception of this section
2. "Electrons cannot exist outside the orbital" -- confusing boundary surface with physical barrier
3. "The conventional representation (isovalue = 0.03) is the 'true' orbital shape"

**Pedagogical Intent:**
Core question for LO12. Directly confronts the "orbital as solid object" misconception using hands-on manipulation. (LO12, Understand)

**Grading Notes:**
- 5 pts: Describes both observations + concludes no definite edge + explains isovalue as threshold
- 4 pts: Describes observations + concludes no definite edge, partial explanation
- 3 pts: Correct observations but incomplete conclusion about boundaries
- 2 pts: Partial observations, recognizes something about the threshold
- 1 pt: Describes one observation
- 0 pts: No answer or "the orbital has a definite edge"

---

### Q3.8: H2O Orbital Classification Table (6 points)

**Question Text:**
> Complete the orbital classification table for H2O (STO-3G) MOs 1-5.

**Expected Answer:**

| MO Index | Approximate Description | Bonding / Nonbonding / Antibonding | Key Visual Feature |
|----------|------------------------|-------------------------------------|-------------------|
| 1 | O 1s core orbital | Nonbonding (core) | Small sphere on O, no H contribution |
| 2 | O-H bonding (sigma, symmetric) | Bonding | Spread over O and both H atoms, density between O-H |
| 3 | O-H bonding (in-plane, antisymmetric) | Bonding | Electron density between O and H, lobed shape in molecular plane |
| 4 | O-H bonding + lone pair mix | Bonding (with lone pair character) | Directed along O-H bonds with some lone pair character on O |
| 5 (HOMO) | O 2p lone pair (out-of-plane) | Nonbonding | Two lobes on O perpendicular to molecular plane, no H contribution |

**PySCF MO Analysis:**

| MO | Energy (Ha) | Dominant AO Components |
|----|-------------|----------------------|
| 1 | -20.2420 | O 1s (0.994) |
| 2 | -1.2682 | O 2s (0.834), H1+H2 (0.159 each) |
| 3 | -0.6174 | O 2py (0.607), H1 (0.445), H2 (-0.445) |
| 4 | -0.4532 | O 2pz (0.776), O 2s (-0.537), H1+H2 (0.278 each) |
| 5 | -0.3913 | O 2px (1.000) -- pure lone pair |

**Acceptable Range:**
- MO 1: Must identify as core or nonbonding, centered on O
- MO 2: Must identify as bonding
- MO 3: Must identify as bonding
- MO 4: Accept "bonding" or "mixed bonding/nonbonding"
- MO 5: Must identify as nonbonding (lone pair)
- Visual descriptions should be qualitatively reasonable

**Common Misconceptions:**
1. Students classify MO 1 as bonding (all occupied MOs must be bonding)
2. Students confuse MO 4 and MO 5 classifications
3. Students may classify the HOMO as antibonding because of its two-lobed shape
4. Students may not recognize MO 3 as bonding because it has lobes of opposite sign

**Pedagogical Intent:**
Summative assessment of orbital classification skills. Requires systematic application of the bonding/nonbonding/antibonding framework to all five occupied MOs. (LO8, Analyze)

**Grading Notes:**
- 6 pts: All 5 MOs correctly classified with reasonable visual descriptions
- 5 pts: 4 of 5 correct
- 4 pts: 3 of 5 correct
- 3 pts: 2 of 5 correct with reasonable visual descriptions
- 2 pts: 1 of 5 correct with some attempt at others
- 1 pt: Some attempt, mostly incorrect
- 0 pts: No attempt

---

## Section 4: Synthesis (Q4.1-Q4.2, 10 points)

### Q4.1: Connecting Bonding Orbital Shape to PES Minimum (5 points)

**Question Text:**
> How does the shape of the bonding orbital explain why the energy minimum exists? What would happen to the orbital shape -- and to the stabilization energy -- if R were much larger?

**Expected Answer:**

**Connection:** The bonding orbital (sigma_g) places significant electron density between the two nuclei. This internuclear electron density creates an attractive electrostatic interaction: both nuclei are attracted toward the electron density concentrated between them, stabilizing the molecule. This stabilization is reflected in the PES minimum.

**At larger R:** As R increases:
1. The bonding orbital becomes more diffuse -- the electron density between the nuclei decreases as the atomic orbitals on each atom become more separated.
2. Eventually, at very large R, the electron density between the nuclei approaches zero, and the bonding orbital essentially "splits" into two separate atomic orbitals.
3. The stabilization energy (the depth of the PES well) diminishes because the electrons can no longer effectively "bridge" the two nuclei.
4. This is the physical reason the energy rises from the minimum toward the dissociation limit as R increases.

**Acceptable Range:**
- Must connect electron density between nuclei to energy lowering
- Must describe what happens to the orbital at large R

**Common Misconceptions:**
1. Students may say "the orbital gets bigger" without connecting to the PES
2. Students may not understand that delocalization is the key to bonding stabilization

**Pedagogical Intent:**
Synthesis question connecting Sections 2 and 3. Tests ability to integrate PES and orbital concepts into a unified picture. (LO7, LO8)

**Grading Notes:**
- 5 pts: Clear connection between orbital shape and energy + correct description of large-R behavior
- 4 pts: Good connection, partial large-R description
- 3 pts: Either connection or large-R description, not both
- 2 pts: Partial answer with some correct elements
- 1 pt: Vague attempt
- 0 pts: No answer

---

### Q4.2: Reflecting on Multiple Representations (5 points)

**Question Text:**
> How does combining all three representations (3D viewer, PES curve, orbital isosurface) give you a more complete picture of a molecule than any single representation alone?

**Expected Answer:**

Each representation shows a different aspect of molecular reality:

1. **3D molecular viewer** shows atomic positions and connectivity -- it answers "where are the atoms?" This is the most intuitive but says nothing about electronic structure or energetics.

2. **PES curve** shows how energy depends on geometry -- it answers "why is this geometry preferred?" and "how strong is the bond?" It captures the energetics but provides no spatial information about where electrons are.

3. **Orbital isosurfaces** show electron distribution -- they answer "where are the electrons?" and "how are electrons shared between atoms?" They explain the electronic mechanism behind bonding but do not directly show the energy consequences.

**Combining all three:** A molecule is simultaneously a collection of atoms in space (viewer), an energy minimum on a potential surface (PES), and a quantum mechanical electron distribution (orbitals). Understanding requires all three:
- The viewer shows the geometry that corresponds to the PES minimum
- The PES explains why that geometry is preferred (energy is lowest)
- The orbitals explain how electrons create that energy minimum (bonding orbitals stabilize, antibonding destabilize)

No single representation captures the full picture. A student who only sees the PES might think of atoms as balls on a spring; a student who only sees orbitals might not understand why a specific geometry is preferred.

**Acceptable Range:**
- Must identify what each representation uniquely contributes
- Must explain why combining them is superior to any single view
- Accept any reasonable synthesis

**Common Misconceptions:**
1. Students may list the three representations without explaining what each uniquely contributes
2. Students may focus on one as "best" rather than complementary

**Pedagogical Intent:**
Metacognitive reflection on representational competence (Kozma & Russell, 2005). This question assesses whether students can articulate why multiple representations matter. (LO7, LO8, LO9)

**Grading Notes:**
- 5 pts: Clearly articulates unique contribution of each representation AND explains why combination is superior
- 4 pts: Identifies contributions of 2-3 representations, good synthesis
- 3 pts: Lists representations with some synthesis
- 2 pts: Lists representations without meaningful synthesis
- 1 pt: Vague response
- 0 pts: No answer

---

## Common Student Errors and Remediation

### Section 1: 3D Molecular Exploration

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Predicting linear H2O geometry** | Common (20-30%) | Review VSEPR from general chemistry. Ask: "How many electron pairs are around oxygen? How does that affect geometry?" Point to the 3D view as direct evidence. |
| **Confusing basis function indices with atom indices** | Very common (40-50%) | Emphasize the mapping: functions 1-5 are on O, 6 is H1, 7 is H2. Write it on the board. Have students label the overlap matrix with atom identifiers. |
| **Thinking overlap = bond order** | Moderate (15-25%) | Clarify: overlap (S) measures how much two functions share space; it is not the same as bond order or bond strength. Two functions can overlap without contributing to bonding (e.g., nonbonding overlaps). |
| **Assuming S_ij = 0 for non-bonded atoms** | Moderate (15-25%) | Show that H-H overlap in H2O is 0.25, not zero. Explain that overlap decreases with distance but never reaches zero at finite distance. |

### Section 2: PES Scanning

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Drawing symmetric parabola for PES** | Common (30-40%) | After the scan, point out the asymmetry: the repulsive wall is steeper than the dissociation curve. Discuss why (1/R divergence vs. gradual orbital decoupling). |
| **Attributing dissociation error to basis set** | Very common (40-50%) | This is the most important error to address. Show that the dissociation error persists with larger basis sets. The issue is the RHF method, not the basis set. You may run a 6-31G PES scan to demonstrate. |
| **Confusing "more negative = higher energy"** | Moderate (20-30%) | Clarify sign convention: in atomic units, lower (more negative) energy = more stable. The energy minimum is the most negative point on the PES. |
| **Thinking equilibrium R should match experimental R_eq** | Moderate (15-25%) | Explain that computed equilibrium depends on basis set and method. STO-3G gives R_eq approximately 1.35 bohr; experimental is 1.40 bohr (0.74 Angstrom). |

### Section 3: Orbital Visualization

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Treating isosurface as physical boundary** | Very common (50-60%) | This is the primary target misconception. Have students adjust the isovalue slider and observe the surface change. Ask: "If the orbital has a real edge, why does the edge move when you change the slider?" |
| **Classifying all two-lobed orbitals as antibonding** | Common (30-40%) | Distinguish: two lobes can indicate a p-orbital (nonbonding), a sigma bond with a node (bonding with contributions from both atoms), or an antibonding orbital. The key diagnostic is whether there is a node between bonded atoms AND the orbital has contributions from both atoms. |
| **Confusing orbital shape with orbital energy** | Moderate (20-30%) | Remind students that shape (spatial distribution) and energy (eigenvalue) are different properties. A compact core orbital has the lowest energy; a diffuse valence orbital has higher energy. |
| **Thinking all occupied orbitals are bonding** | Common (25-35%) | Point to MO 1 (core, nonbonding) and MO 5 / HOMO (lone pair, nonbonding). Ask: "If MO 1 contributed to bonding, what would change about the molecule if we removed it from the bonding picture?" |

---

## Performance Task Rubrics

### PES Interpretation Rubric (4-point scale)

**Covers:** LO7 (geometry-energy connection), LO11a/LO11b (dissociation limits)

**Task description:** Given the computed H2 PES curve, students interpret the curve, explain the physical origin of the minimum, and analyze the dissociation limit behavior.

| Score | Level | Criteria |
|-------|-------|----------|
| **4** | **Exemplary** | **Equilibrium:** Correctly identifies the equilibrium bond length and energy. **Physical origin:** Explains the minimum as a balance between nuclear repulsion (short R) and electron delocalization stabilization (intermediate R). **Dissociation:** Recognizes that RHF gives an incorrect dissociation limit, identifies the energy discrepancy quantitatively, and explains that the single-determinant constraint (restricted orbital sharing) prevents correct dissociation to neutral atoms. May mention ionic contamination or multi-reference character. |
| **3** | **Proficient** | **Equilibrium:** Correctly identifies equilibrium with approximate values. **Physical origin:** Explains the minimum as a balance between competing effects, though explanation may be incomplete (e.g., mentions repulsion and attraction without specifying mechanisms). **Dissociation:** Recognizes RHF gives incorrect large-R behavior but explanation of why is partial (e.g., "RHF has limitations" without connecting to orbital restriction). |
| **2** | **Developing** | **Equilibrium:** Identifies the minimum from the curve but may not extract accurate numerical values. **Physical origin:** Provides a partial explanation (e.g., "atoms attract at medium distance") without distinguishing types of interactions. **Dissociation:** May not notice the dissociation limit problem, or notices it but cannot explain it. |
| **1** | **Beginning** | **Equilibrium:** Cannot accurately identify the equilibrium from the curve, or misinterprets the curve (e.g., identifies maximum as equilibrium). **Physical origin:** No meaningful explanation of why the minimum exists, or explanation is fundamentally incorrect. **Dissociation:** Does not address dissociation behavior. |

**Scoring notes:**
- Each sub-criterion (equilibrium, physical origin, dissociation) contributes roughly equally
- Award the score that best matches the overall response
- For undergraduate students, do not penalize for omitting LO11b content (graduate extension)
- Exemplary responses at the graduate level should address the sigma_g^2 / sigma_u*^2 configuration mixing

---

### Orbital Classification Rubric (4-point scale)

**Covers:** LO8 (orbital interpretation), LO12 (isovalue interpretation)

**Task description:** Given orbital isosurface visualizations, students classify orbitals as bonding/nonbonding/antibonding and demonstrate understanding of what the isovalue threshold means.

| Score | Level | Criteria |
|-------|-------|----------|
| **4** | **Exemplary** | **Classification:** Correctly classifies all examined orbitals (core, bonding, nonbonding/lone pair, antibonding) using shape-based evidence. Explains the diagnostic criteria (node between bonded atoms = antibonding, density between nuclei = bonding, localized on one atom = nonbonding). **Isovalue:** Correctly explains that the isovalue is a chosen threshold on |psi|, not a physical boundary. States that the orbital extends to infinity and has no definite edge. May connect to probability interpretation (e.g., "a given isovalue encloses a region containing some percentage of the electron density"). |
| **3** | **Proficient** | **Classification:** Correctly classifies most orbitals (4 of 5 for H2O, or all H2 orbitals) with shape-based reasoning. One classification may be incorrect or weakly justified. **Isovalue:** Understands that the isosurface changes with isovalue and that the boundary is not physical, but may not fully articulate why (e.g., "the orbital gets bigger when I decrease the isovalue" without connecting to the mathematical meaning of the threshold). |
| **2** | **Developing** | **Classification:** Correctly classifies some orbitals (2-3 of 5) but shows confusion between categories (e.g., calls lone pair "antibonding" because of two lobes, or calls core orbital "bonding" because it is occupied). **Isovalue:** Recognizes that the slider changes the surface but may believe the surface at some "correct" isovalue represents the "true" orbital boundary. May state that the orbital has a definite edge at the "right" isovalue. |
| **1** | **Beginning** | **Classification:** Cannot correctly classify orbitals from shapes, or classifies based on incorrect criteria (e.g., "occupied = bonding, virtual = antibonding"). Shows fundamental confusion between bonding/nonbonding/antibonding. **Isovalue:** Treats the isosurface as the physical boundary of the orbital. States or implies that electrons exist inside the surface and not outside it. May express the "planetary orbit" misconception. |

**Scoring notes:**
- Classification and isovalue understanding contribute roughly equally
- The misconception check answer in Q3.7 is a strong diagnostic for the isovalue score level
- Students at level 1 may benefit from one-on-one discussion after the lab

---

## FAQ for Common Student Issues with 3D Interaction

### "The 3D viewer is blank or shows 'WebGL required'"

**Cause:** The browser does not have WebGL enabled or hardware acceleration is disabled.

**Solutions:**
1. Try a different browser (Chrome or Firefox recommended; Safari may have WebGL disabled by default)
2. Enable hardware acceleration in browser settings:
   - Chrome: Settings > System > "Use hardware acceleration when available"
   - Firefox: Settings > Performance > uncheck "Use recommended performance settings," then check "Use hardware acceleration when available"
3. Update graphics drivers on the system
4. If using a Chromebook or thin client, WebGL may not be supported -- try a different machine

### "I can't rotate the molecule"

**Cause:** The user is not clicking and dragging within the 3D canvas area, or touch events are being intercepted.

**Solutions:**
1. Click directly on the 3D viewer area (not outside it) and drag while holding the mouse button
2. Scroll wheel to zoom in/out
3. Right-click drag for panning (translation)
4. On touchscreens: one finger to rotate, pinch to zoom, two fingers to pan
5. If using a trackpad, try clicking and dragging rather than two-finger gestures

### "The orbital surface doesn't appear after selecting an orbital"

**Cause:** The SCF calculation must converge before orbital data is available. Also, the orbital grid evaluation may take a moment.

**Solutions:**
1. Run the SCF calculation first and wait for convergence
2. After convergence, switch to the Orbitals tab
3. Select an orbital from the dropdown or list
4. Wait for the grid evaluation (a progress indicator may appear)
5. If the surface still does not appear, try adjusting the isovalue slightly (the default 0.03 may be too high for some diffuse orbitals)

### "The orbital colors look different from what the worksheet describes"

**Cause:** Browser rendering differences, color profile settings, or accessibility features.

**Solutions:**
1. IQCP uses solid rendering for positive lobes and translucent rendering for negative lobes -- the key distinction is opacity, not color
2. Check that your browser does not have a color filter or accessibility mode active
3. The exact appearance may vary between browsers; the shape and lobe count are the important features
4. If you cannot distinguish the two phases, try rotating the molecule -- the transparency difference is usually visible from multiple angles

### "The isovalue slider seems to have no effect"

**Cause:** Grid evaluation is running in the background, or the isovalue change is too small to see.

**Solutions:**
1. Make a large change (e.g., from 0.03 to 0.01) to see a visible effect
2. Wait a moment after moving the slider -- the marching cubes algorithm needs to recalculate the surface
3. If you set the isovalue too high (e.g., > 0.2), the surface may shrink to nothing -- decrease it back toward 0.03
4. If you set the isovalue too low (e.g., < 0.005), the surface may become very large and fill the viewer -- increase it back

### "The PES scan appears to hang or takes very long"

**Cause:** Each point on the PES scan requires a full SCF calculation. With 20 points, this may take 10-30 seconds total.

**Solutions:**
1. Be patient -- a progress indicator should show how many points have completed
2. If it appears stuck on one point, the SCF for that geometry may be having convergence difficulty (very common at short R)
3. Try reducing the number of scan points (e.g., from 20 to 10)
4. PES scanning is restricted to diatomic molecules in this version -- if you selected a polyatomic, the scan will not work
5. If the scan stalls, refresh the page and try again with fewer points

### "My deep links don't load the expected state"

**Cause:** Deep links may have been truncated when copying, or the IQCP version may have changed.

**Solutions:**
1. Copy the entire URL (make sure nothing is cut off)
2. Clear your browser cache and try again
3. If the deep link loads but shows different parameters, manually adjust the settings as described in the worksheet step
4. Report the issue to your instructor for the FAQ

---

## Prior Knowledge Note for LO11b (Graduate Extension)

### Prerequisites for Assigning the LO11b Extension

LO11b asks students to explain the RHF dissociation failure in terms of single-determinant limitations. This requires prior exposure to:

1. **Slater determinants:** Students must understand that the RHF wavefunction is a single antisymmetrized product of orbitals (a Slater determinant). This is typically covered in a graduate-level quantum chemistry course.

2. **Configuration interaction (CI) concept:** Students should have at least a conceptual understanding that the true wavefunction can be expressed as a sum of Slater determinants built from different orbital occupations (configurations).

3. **Electron correlation:** Students should understand that the difference between the Hartree-Fock energy and the exact energy (within a given basis set) is the correlation energy. At the dissociation limit, the correlation energy becomes very large for RHF.

4. **Covalent vs. ionic configurations:** Students should be able to distinguish the covalent configuration (one electron on each atom) from the ionic configuration (both electrons on one atom) and understand that RHF mixes them with equal weight.

### Guidance for Mixed-Level Classes

**Undergraduate-only sections:**
- Assign Q2.6 as the final dissociation question
- Skip the "Going Deeper" extension entirely
- Focus the discussion on the empirical observation: "RHF gives the wrong answer at large R"
- Learning target: LO11a only (identify the limitation, do not explain the mechanism)

**Graduate-only sections:**
- Assign both Q2.6 and the "Going Deeper" extension
- Expect students to discuss sigma_g^2 / sigma_u*^2 configuration mixing
- May extend discussion to UHF vs. RHF and symmetry breaking
- Learning targets: LO11a + LO11b

**Mixed undergraduate/graduate sections:**
- Present Q2.6 to all students
- Present the "Going Deeper" extension as optional extra credit for undergraduates
- Require the extension for graduate students
- During discussion, let graduate students explain the multi-reference concept to undergraduates -- peer teaching reinforces learning for both groups
- Grade undergraduates on LO11a only; grade graduates on LO11a + LO11b

**Instructor discussion prompts for LO11b:**
1. "If we allowed alpha and beta electrons to have different spatial orbitals (UHF), would the dissociation limit improve?"
2. "How many determinants would we need to correctly describe H2 dissociation?"
3. "Does this limitation of RHF matter near equilibrium? When does it become important?"

---

## Discussion Prompts by Section

### Section 1 Discussion Prompts

1. **After Q1.4:** "Looking at the overlap matrix, can you find a pair of basis functions that have zero overlap? Why is that overlap exactly zero?" (Answer: O 2px with O 2py or O 2pz -- orthogonal by symmetry.)

2. **After Q1.7:** "If you could design a basis set, would you want large or small overlap between functions on the same atom? What about between functions on different atoms?" (Explores the concept of basis set quality and linear dependence.)

3. **After Q1.8:** "How would the overlap matrix change if we used a bigger basis set like 6-31G instead of STO-3G?" (More basis functions means a larger matrix, but the physical overlap between atom pairs remains determined by geometry.)

### Section 2 Discussion Prompts

1. **After Q2.3:** "Is the chemical bond purely an electrostatic phenomenon, or does quantum mechanics play an essential role?" (The virial theorem and kinetic energy lowering are fundamentally quantum -- classical electrostatics alone cannot explain covalent bonding.)

2. **After Q2.6:** "Could we fix the RHF dissociation problem by using a very large basis set?" (No -- this is a method limitation, not a basis set limitation. Good transition to discussing electron correlation methods.)

3. **After Q2.8:** "In practice, how do computational chemists choose a basis set for their calculations?" (Balance of accuracy, cost, and the property being studied. Benchmark studies help. There is no universal "best" basis set.)

### Section 3 Discussion Prompts

1. **After Q3.3:** "Can you think of a situation where the HOMO being nonbonding matters for chemistry?" (Nucleophilic attack -- lone pairs are reactive. The HOMO determines where the molecule acts as an electron donor.)

2. **After Q3.7:** "Textbooks often show orbitals as solid shapes with definite boundaries. After using the isovalue slider, do you think this is a good representation? What would be better?" (Leads to discussion of probability density, contour plots, and the 90% probability surface convention.)

3. **After Q3.8:** "If we added electrons to the antibonding orbitals (MO 6 and MO 7 for H2O), what would happen to the molecule?" (Destabilization -- antibonding orbitals weaken or break bonds. Connects to excited states and photochemistry.)

---

## Expected Outputs for Exercises

### Section 1 Expected 3D Views

**H2O 3D View:**
- Bent molecule with O (red) at center/top and two H (white) atoms below
- Bond angle approximately 104.5 degrees
- CPK coloring: O = red, H = white
- With atom labels: "O", "H", "H" visible

**H2 3D View:**
- Linear molecule with two H (white) atoms
- Single bond between them
- Very simple -- just two white spheres connected by a bond

### Section 2 Expected PES Curves

**H2 PES Scan (STO-3G, 0.5-5.0 bohr, 20 points):**
- Characteristic Morse-like curve
- Steep repulsive wall on the left (R < 1.0 bohr)
- Minimum near R = 1.35 bohr, E approximately -1.117 Ha
- Gradual rise toward dissociation on the right
- Energy at R = 5.0 bohr approximately -0.686 Ha
- The equilibrium marker should appear at or near the minimum

**Key numerical values (PySCF reference):**

| R (bohr) | E (Ha) | Notes |
|-----------|--------|-------|
| 0.50 | -0.4033 | Repulsive wall |
| 0.97 | -1.0558 | Descending |
| 1.21 | -1.1114 | Approaching minimum |
| 1.35 | -1.1175 | Near equilibrium |
| 1.45 | -1.1149 | Past equilibrium |
| 2.63 | -0.9438 | Dissociation region |
| 5.00 | -0.6864 | Far from equilibrium |

**RHF dissociation limit (STO-3G):** Converges toward approximately -0.57 Ha at very large R (R = 20 bohr), compared to the correct limit of 2 x E(H) = -0.9332 Ha.

### Section 3 Expected Orbital Isosurfaces

**H2 MO 1 (sigma_g, bonding):**
- Single elongated lobe encompassing both hydrogen atoms
- Electron density concentrated between the nuclei
- Solid rendering (positive phase throughout)

**H2 MO 2 (sigma_u*, antibonding):**
- Two separate lobes, one on each hydrogen atom
- Nodal plane at the midpoint between nuclei
- One lobe solid (positive), one lobe translucent (negative)
- No electron density between nuclei

**H2O MO 1 (core):**
- Small compact sphere centered on oxygen
- No visible contribution from hydrogen atoms
- Very small at default isovalue (0.03) because the 1s orbital is compact

**H2O MO 5 / HOMO (O 2px lone pair):**
- Two lobes above and below the molecular plane (perpendicular to the plane defined by O-H-H)
- Entirely centered on oxygen
- Solid lobe on one side, translucent on the other (opposite phases)
- No electron density on hydrogen atoms

**H2O at various isovalues:**
- Isovalue 0.01: Large, diffuse surface extending far from nuclei
- Isovalue 0.03: "Standard" orbital shape, lobes clearly defined
- Isovalue 0.08: Small, compact surface close to nuclei; diffuse regions disappear

---

## Point Allocation Detail

### Per-Question Point Breakdown

| Question | Points | LO | Bloom's Level |
|----------|--------|-----|---------------|
| Q1.1 | 4 | LO9 | Apply |
| Q1.2 | 3 | LO9 | Analyze |
| Q1.3 | 4 | LO9 | Apply |
| Q1.4 | 4 | LO9 | Apply |
| Q1.5 | 4 | LO9 | Understand |
| Q1.6 | 3 | LO9 | Apply |
| Q1.7 | 4 | LO9 | Analyze |
| Q1.8 | 4 | LO9 | Analyze |
| **Section 1 Total** | **30** | | |
| Q2.1 | 4 | LO7 | Apply |
| Q2.2 | 4 | LO7 | Analyze |
| Q2.3 | 5 | LO7 | Analyze |
| Q2.4 | 4 | LO7 | Analyze |
| Q2.5 | 5 | LO11a | Analyze |
| Q2.6 | 5 | LO11a/LO11b | Analyze/Evaluate |
| Q2.7 | 5 | LO10 | Evaluate |
| Q2.8 | 4 | LO10 | Evaluate |
| **Section 2 Total** | **36** | | |
| Q3.1 | 3 | LO8 | Apply |
| Q3.2 | 4 | LO8 | Analyze |
| Q3.3 | 5 | LO8 | Analyze |
| Q3.4 | 3 | LO8 | Apply |
| Q3.5 | 5 | LO8 | Analyze |
| Q3.6 | 3 | LO12 | Understand |
| Q3.7 | 5 | LO12 | Understand |
| Q3.8 | 6 | LO8 | Analyze |
| **Section 3 Total** | **34** | | |
| Q4.1 | 5 | LO7, LO8 | Analyze |
| Q4.2 | 5 | LO7, LO8, LO9 | Evaluate |
| **Section 4 Total** | **10** | | |
| **Grand Total** | **110** | | |

### Cognitive Level Distribution

| Bloom's Level | Questions | Points | Percentage |
|---------------|-----------|--------|------------|
| Understand | Q1.5, Q3.6, Q3.7 | 12 | 10.9% |
| Apply | Q1.1, Q1.3, Q1.4, Q1.6, Q2.1, Q3.1, Q3.4 | 24 | 21.8% |
| Analyze | Q1.2, Q1.7, Q1.8, Q2.2, Q2.3, Q2.4, Q2.5, Q2.6, Q3.2, Q3.3, Q3.5, Q3.8, Q4.1 | 59 | 53.6% |
| Evaluate | Q2.7, Q2.8, Q4.2 | 14 | 12.7% |
| **Subtotals** | | | |
| Remember/Understand | | 12 | **10.9%** |
| Apply/Analyze | | 83 | **75.5%** |
| Evaluate | | 14 | **12.7%** |

**Note:** The cognitive level distribution is weighted toward Apply/Analyze (75.5%), reflecting the lab's emphasis on hands-on exploration and interpretation. The Evaluate level (12.7%) captures the basis set evaluation and synthesis questions. The distribution departs slightly from the 20/50/30 target in the development checklist because the POE structure naturally emphasizes Apply and Analyze activities.

### Converting to Course Grade

| Raw Score (out of 110) | Percentage | Suggested Grade |
|------------------------|------------|-----------------|
| 99-110 | 90-100% | A |
| 88-98 | 80-89% | B |
| 77-87 | 70-79% | C |
| 66-76 | 60-69% | D |
| 0-65 | < 60% | F |

---

## Timing Guidance

### Standard Pacing (60 minutes)

| Activity | Time | Notes |
|----------|------|-------|
| Section 1: 3D Exploration | 12-15 min | Q1.1-Q1.8 |
| Section 2: PES Scanning | 18-22 min | Q2.1-Q2.8, includes scan wait time |
| Section 3: Orbital Visualization | 20-25 min | Q3.1-Q3.8, includes rendering time |
| Section 4: Synthesis | 3-5 min | Q4.1-Q4.2 |
| **Total** | **53-67 min** | Target: 60 min |

### Adjustments by Class Level

**Introductory (general chemistry level):**
- Add 5-10 min buffer
- Skip "Going Deeper" extension in Q2.6
- Consider omitting Q2.8 (basis set evaluation) or making it a take-home question
- Provide more explicit guidance for 3D viewer interaction
- Target: 65-75 min

**Intermediate (physical chemistry level):**
- Standard pacing as above
- Assign "Going Deeper" as optional
- Target: 55-65 min

**Advanced (graduate level):**
- Require "Going Deeper" extension in Q2.6
- Expect more quantitative answers (e.g., energy values to more decimal places)
- Add discussion time after Section 2 for multi-reference methods
- Target: 55-65 min (faster on mechanics, slower on deep questions)

### Technology Setup Requirements

- Computer lab with modern web browsers (Chrome or Firefox preferred)
- WebGL-capable graphics (virtually all machines from 2015 onward)
- Internet access for iqcp.dev (or local deployment)
- Projector for instructor to demonstrate 3D interaction if needed
- Students should be instructed to clear browser cache before starting if they have used IQCP before

---

*Lab Pack #2 Instructor Key v1.0 | CONFIDENTIAL -- Instructor Use Only*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
