# Lab Pack #3: Instructor Answer Key and Teaching Guide

**Lab Pack:** 3 - Computational Layers of Quantum Chemistry
**Version:** 1.0
**Last Updated:** 2026-03-19
**Document Type:** Instructor Materials (CONFIDENTIAL -- Do Not Distribute to Students)

---

## Document Overview

This instructor key covers all 26 numbered questions (Q1.1-Q4.2) in the Lab Pack #3 student worksheet. It also provides:

- Expected outputs for every exercise (radial profiles, matrix heatmaps, density isosurfaces)
- Common student errors with remediation strategies for each section
- Performance task rubrics on 4-point scales (basis set analysis, integral interpretation)
- FAQ for Module A, Module B, and density-specific issues
- Point allocation summary
- Timing guidance with class-level adjustments

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
- NumPy 1.26.4
- Python 3.12
- Date verified: 2026-03-19

### Point Allocation Summary

| Section | Questions | Points | Primary LOs |
|---------|-----------|--------|-------------|
| Section 1: Basis Function Exploration | Q1.1-Q1.8 | 30 | LO13, LO14 |
| Section 2: Integral Inspection & Fock Tracing | Q2.1-Q2.8 | 34 | LO15, LO16 |
| Section 3: Electron Density & Difference Density | Q3.1-Q3.8 | 32 | LO17, LO18 |
| Section 4: Synthesis | Q4.1-Q4.2 | 10 | Integrative |
| **Total** | **26 questions** | **106** | |

**Note:** Point values may be scaled to 100 points by instructors using the formula: (raw score / 106) x 100.

---

## Section 1: Basis Function Anatomy and Comparison (Q1.1-Q1.8, 30 points)

### Q1.1: Predicting Contracted Gaussian Shape (3 points)

**Question Text:**
> A single Gaussian has the form g(r) = N exp(-alpha r^2). The STO-3G basis set for hydrogen uses three Gaussians combined. Predict: will the contracted function look like (a) a wider bell curve, (b) a sharper, more peaked curve, or (c) three separate bumps? Sketch the radial profile.

**Expected Answer:**

**(b) A sharper, more peaked curve** that approximates a Slater-type orbital (which has a cusp at r = 0 and an exponential tail).

The contracted function does not look like three separate bumps because the three primitives are added together with positive coefficients at every point in space. The result is a single smooth curve. It is sharper than any individual Gaussian because the tight (large-exponent) primitive contributes a sharp peak near the nucleus, while the diffuse (small-exponent) primitive extends the tail. The combination achieves a shape closer to the exponential decay of a Slater orbital than any single Gaussian could.

**STO-3G H 1s primitives (PySCF):**

| Primitive | Exponent (alpha) | Coefficient (d) |
|-----------|-----------------|-----------------|
| 1 (tight) | 3.42525091 | 0.15432897 |
| 2 (medium) | 0.62391373 | 0.53532814 |
| 3 (diffuse) | 0.16885540 | 0.44463454 |

**Acceptable Range:**
- Answer (b) with a sketch showing a single smooth curve (sharper than a bell curve near r = 0)
- Accept "something between (a) and (b)" if the student recognizes it is a single function

**Common Misconceptions:**
1. Students select (c), expecting three separate bumps -- they do not realize that contraction sums the primitives at every point
2. Students select (a), thinking more functions means a wider spread
3. Students draw a Gaussian rather than a Slater-like function

**Pedagogical Intent:**
POE-predict stage. Activates prior knowledge about Gaussian functions and sets up the observation of the actual radial profile. (LO13)

**Grading Notes:**
- 3 pts: Correct prediction (b) with a sketch showing a single smooth, Slater-like curve
- 2 pts: Correct prediction (b) with a sketch but shape is just a Gaussian bell curve
- 1 pt: Incorrect prediction but shows some reasoning about combining functions
- 0 pts: No prediction or no sketch

---

### Q1.2: Describing the STO-3G H 1s Radial Profile (4 points)

**Question Text:**
> Describe the radial profile you observe. How does the shape of the contracted function compare to the individual primitives? Which primitive contributes most near the nucleus? Which contributes most at larger distances?

**Expected Answer:**

The contracted function (solid line) is a smooth, single-peaked curve that rises sharply near r = 0 and decays approximately exponentially. It is sharper near the nucleus than a single Gaussian but broader in the tail.

- **Near the nucleus (r < 0.5 bohr):** Primitive 1 (alpha = 3.4253, the tightest) contributes the most. It has the largest exponent, so it is the most concentrated near r = 0. However, its coefficient (0.1543) is small, so it contributes a sharp but relatively small peak.
- **At intermediate distances (r ~ 0.5-1.5 bohr):** Primitive 2 (alpha = 0.6239) contributes the most. It has the largest coefficient (0.5353) and intermediate width.
- **At larger distances (r > 2 bohr):** Primitive 3 (alpha = 0.1689, the most diffuse) contributes the most. It has the smallest exponent and therefore decays most slowly.

**Acceptable Range:**
- Must correctly identify the tightest primitive as dominating near the nucleus
- Must correctly identify the most diffuse primitive as dominating at large r
- Description of overall shape as "smooth" and "single-peaked"

**Common Misconceptions:**
1. Students may confuse exponent magnitude with coefficient magnitude -- a large exponent means tight (concentrated near nucleus), not "bigger contribution everywhere"
2. Students may think the primitives add at some points and cancel at others (they do not cancel because all STO-3G coefficients for H 1s are positive)

**Pedagogical Intent:**
POE-observe/explain. Builds understanding of how contraction works -- different primitives serve different spatial regions. (LO13)

**Grading Notes:**
- 4 pts: Correctly identifies tight primitive at nucleus, diffuse at large r, and describes overall shape
- 3 pts: Two of three correct
- 2 pts: One correct with reasonable description
- 1 pt: Vague description without identifying primitive roles
- 0 pts: No answer

---

### Q1.3: Why Use 3 Primitives Instead of 1? (4 points)

**Question Text:**
> STO-3G uses 3 primitives to approximate a Slater-type orbital. Why not just use a single Gaussian? What does adding more primitives with different exponents achieve? What is the tradeoff?

**Expected Answer:**

**Why not a single Gaussian:** A single Gaussian has the wrong shape to represent an atomic orbital. Specifically:
- A Gaussian has zero slope at r = 0 (smooth maximum), but a Slater-type orbital has a cusp (discontinuous derivative) at the nucleus.
- A Gaussian decays as exp(-alpha r^2), which falls off too rapidly at large r compared to the exponential decay exp(-zeta r) of a Slater orbital.

**What multiple primitives achieve:** By combining three Gaussians with different exponents:
- The tight primitive captures the sharp peak near the nucleus
- The diffuse primitive extends the tail to larger distances
- The combination approximates the Slater orbital shape much better than any single Gaussian

**The tradeoff:** Each additional primitive increases the number of integrals that must be computed. For a basis with N contracted functions, each backed by K primitives, the number of primitive integral evaluations scales as K^4 for two-electron integrals. More primitives give a better fit to the Slater orbital but cost more computation.

**Acceptable Range:**
- Must mention shape mismatch (cusp or tail behavior)
- Must mention computational cost as the tradeoff
- Accept qualitative arguments about "better approximation"

**Common Misconceptions:**
1. Students may think a single Gaussian is "just as good" because it is also bell-shaped
2. Students may not recognize the cusp/tail distinction between Gaussians and Slater orbitals
3. Students may think the tradeoff is accuracy (less accurate with more primitives) rather than cost

**Pedagogical Intent:**
Tests understanding of why contraction exists -- the fundamental tension between mathematical convenience (Gaussians have analytic integrals) and physical accuracy (Slater orbitals have the right shape). (LO13)

**Grading Notes:**
- 4 pts: Mentions shape mismatch (cusp and/or tail) AND computational cost tradeoff
- 3 pts: Shape mismatch explained, cost mentioned briefly
- 2 pts: Either shape or cost, not both
- 1 pt: Vague answer about "better accuracy"
- 0 pts: No answer

---

### Q1.4: Misconception Check -- Basis Function vs. Atomic Orbital (4 points)

**Question Text:**
> A student says: "A basis function is the same thing as an atomic orbital." Do you agree? Consider: does the O 2p basis function in STO-3G have exactly the same shape as a hydrogen-like 2p orbital?

**Expected Answer:**

**Disagree.** A basis function and an atomic orbital are related but not the same.

**Key distinctions:**
1. **Atomic orbitals** are the exact eigenfunctions of the hydrogen-like atom (one-electron Schrodinger equation). They have specific shapes (1s, 2p, etc.) and exact radial forms involving Laguerre polynomials and exponential decay.
2. **Basis functions** are mathematical approximations chosen for computational convenience. In STO-3G, each basis function is a contracted Gaussian that approximates a Slater-type orbital, which itself is only an approximation to a hydrogen-like orbital.
3. The O 2p basis function in STO-3G does NOT have the same shape as a hydrogen-like 2p orbital. It has the correct angular part (proportional to x, y, or z) but the radial part is a sum of Gaussians, not a Laguerre polynomial times an exponential.
4. In a molecular context, the basis functions are centered on atoms, but the molecular orbitals that result from the SCF calculation are linear combinations of basis functions on ALL atoms.

**Acceptable Range:**
- Must disagree with the statement
- Must cite at least one distinction (shape, mathematical form, or role in calculation)

**Common Misconceptions:**
1. Students agree because basis functions are "named like atomic orbitals" (1s, 2p, etc.)
2. Students think basis functions ARE the orbitals of the atom in the molecule
3. Students confuse basis functions (inputs to the calculation) with molecular orbitals (outputs)

**Pedagogical Intent:**
Directly targets the "basis function = atomic orbital" misconception. This is one of the most common conceptual errors in introductory computational chemistry. (LO13)

**Grading Notes:**
- 4 pts: Correctly disagrees with at least two distinctions (shape, role, mathematical form)
- 3 pts: Correctly disagrees with one clear distinction
- 2 pts: Correctly disagrees but reasoning is vague
- 1 pt: Agrees but shows some uncertainty
- 0 pts: Fully agrees or no answer

---

### Q1.5: Predicting STO-3G vs. 6-31G for O 2s (3 points)

**Question Text:**
> The 6-31G basis set is "split-valence" -- it uses two sets of contractions for valence shells. Predict: will the O 2s radial profile in 6-31G extend further from the nucleus, closer to the nucleus, or be about the same as in STO-3G?

**Expected Answer:**

The O 2s radial profile in 6-31G will **extend further from the nucleus** compared to STO-3G. This is because 6-31G splits the valence into two contractions: a tighter 3-primitive contraction and a single diffuse (uncontracted) function with a small exponent (alpha = 0.2700). This diffuse function extends the reach of the basis set to larger distances, providing more flexibility to describe valence electron distribution.

**PySCF Reference:**
- STO-3G O 2s: 3 primitives (alpha = 5.033, 1.170, 0.380)
- 6-31G O 2s inner: 3 primitives (alpha = 15.54, 3.600, 1.014)
- 6-31G O 2s outer: 1 primitive (alpha = 0.270) -- this extends much further

**Acceptable Range:**
- "Further" or "extends to larger distances" for full credit
- "About the same but with more detail" for partial credit

**Common Misconceptions:**
1. Students predict "closer to the nucleus" because more basis functions = "more concentrated"
2. Students predict "same" because it is still oxygen

**Pedagogical Intent:**
POE-predict for the basis set comparison. Tests whether students understand what "split-valence" means before they see the comparison plot. (LO14)

**Grading Notes:**
- 3 pts: Predicts "further" with reasoning about the diffuse function or split-valence flexibility
- 2 pts: Predicts "further" without clear reasoning
- 1 pt: Incorrect prediction but some reasoning about basis set differences
- 0 pts: No prediction

---

### Q1.6: Describing STO-3G vs. 6-31G Differences for O 2s (4 points)

**Question Text:**
> Describe the differences between the O 2s radial profiles in STO-3G and 6-31G. How many contractions does each basis set use for the valence 2s? What additional flexibility does the split-valence approach provide?

**Expected Answer:**

**STO-3G O 2s:** Uses a single contraction of 3 primitives. The radial profile is a single smooth curve that decays relatively quickly.

**6-31G O 2s:** Uses two contractions for the valence -- a 3-primitive "inner" contraction and a 1-primitive "outer" (uncontracted) function. In the comparison plot:
- The 6-31G profile extends further from the nucleus due to the diffuse outer function
- The 6-31G effectively has two components that can be independently weighted by the SCF procedure

**Additional flexibility:** In STO-3G, the shape of the 2s radial function is completely fixed -- the SCF procedure can only adjust how much of this function to include in each MO. In 6-31G, the inner and outer parts are separate basis functions that can be independently weighted. This allows the SCF to find the optimal balance between "tight" and "diffuse" character, producing a better description of the valence region. This is especially important for describing bonding, where electron density shifts between atoms.

**PySCF Reference:**
- STO-3G O 2s: 3 primitives, 1 contraction -> 1 basis function
- 6-31G O 2s: (3 + 1) primitives, 2 contractions -> 2 basis functions

**Acceptable Range:**
- Must identify that 6-31G has 2 contractions vs. STO-3G's 1
- Must mention the flexibility of independent weighting or variational freedom

**Common Misconceptions:**
1. Students may think 6-31G simply has "more Gaussians" without understanding the split
2. Students may not grasp that the two contractions are independent basis functions

**Pedagogical Intent:**
POE-observe/explain. Tests comprehension of the split-valence concept through direct visual comparison. (LO14)

**Grading Notes:**
- 4 pts: Identifies contraction counts for both, describes visual differences, explains flexibility
- 3 pts: Identifies contraction counts and one of visual differences or flexibility
- 2 pts: Notes some visual difference without connecting to contraction structure
- 1 pt: Vague comparison
- 0 pts: No answer

---

### Q1.7: What Does "Split-Valence" Mean in 6-31G? (4 points)

**Question Text:**
> What does "split-valence" mean in 6-31G? The name encodes the structure: "6" primitives for the core, "31" means the valence is split into a 3-primitive contraction and a 1-primitive function. Why split the valence shells but not the core?

**Expected Answer:**

**Split-valence meaning:** The basis set name "6-31G" encodes its structure:
- **6:** The core shell (e.g., O 1s) is a single contraction of 6 primitive Gaussians. It is not split.
- **3-1:** Each valence shell is represented by two contractions: one with 3 primitives (inner/tight) and one with 1 primitive (outer/diffuse). The SCF procedure can independently weight these two components.

**Why split valence but not core:** Core electrons (e.g., O 1s) are:
1. **Tightly bound** to the nucleus with very low energy (-20.24 Ha for O 1s)
2. **Chemically inert** -- they do not significantly participate in bonding
3. **Essentially the same in every molecule** -- the O 1s orbital looks the same whether oxygen is in H2O, CO2, or O2

Since core electrons do not change significantly between different molecular environments, there is no need to give them extra variational flexibility. Splitting the core would add basis functions (increasing computational cost by N^4 scaling) with negligible improvement in energy or properties.

Valence electrons, by contrast, are the ones that form bonds, so they need the flexibility to redistribute between atoms. The split-valence approach gives this flexibility at minimal additional cost.

**Acceptable Range:**
- Must explain the naming convention (6 = core, 31 = split valence)
- Must explain why core is not split (chemically inert / tightly bound)

**Common Misconceptions:**
1. Students may think "31" means 31 total primitives
2. Students may think the core is not split because it is "too small" rather than chemically inert
3. Students may not connect splitting to variational freedom

**Pedagogical Intent:**
Tests understanding of basis set design philosophy -- cost vs. accuracy tradeoffs in different shells. (LO14)

**Grading Notes:**
- 4 pts: Correct naming convention explanation AND correct reasoning for not splitting core
- 3 pts: Correct naming convention with partial core reasoning
- 2 pts: Correct naming convention only
- 1 pt: Partial understanding of naming
- 0 pts: No answer or fundamentally wrong interpretation

---

### Q1.8: Misconception Check -- More Basis Functions = Always Better? (4 points)

**Question Text:**
> A student claims: "Adding more basis functions always makes the answer more accurate, so we should always use the biggest basis set available." Evaluate this claim. What practical consideration limits the benefit?

**Expected Answer:**

**Partially correct, but misleading.** The variational principle guarantees that a larger basis set will give an energy that is lower (or equal) compared to a smaller basis. In that narrow sense, "more accurate" is true for the total energy.

However, several practical considerations limit the benefit:

1. **Computational cost scaling:** The number of two-electron integrals scales as O(N^4) where N is the number of basis functions. Going from STO-3G (7 functions for H2O) to 6-31G (13 functions) increases the number of ERIs by a factor of approximately (13/7)^4 = 11.9. For larger molecules, this becomes prohibitive.

2. **Diminishing returns:** The energy improvement from STO-3G to 6-31G is large (about 1 Ha for H2O), but from 6-31G to cc-pVTZ the improvement is much smaller. Each additional function contributes less.

3. **Method limitations persist:** Even with a complete basis set (CBS limit), RHF still has the electron correlation error. A huge basis set with a limited method does not fix fundamental method limitations (as seen with RHF dissociation in Lab Pack #2).

4. **Linear dependence risk:** Very large basis sets can introduce near-linear dependencies in the overlap matrix, causing numerical instability.

**PySCF Reference:**
- H2O STO-3G (7 functions): -74.9630 Ha
- H2O 6-31G (13 functions): -75.9840 Ha
- Difference: -1.021 Ha (large improvement)

**Acceptable Range:**
- Must partially disagree
- Must cite at least one practical consideration (cost scaling or diminishing returns)

**Common Misconceptions:**
1. Students fully agree without reservation
2. Students think accuracy improvement is linear in the number of functions
3. Students confuse basis set completeness with method accuracy

**Pedagogical Intent:**
Evaluate-level question testing the ability to critically assess a claim about computational methodology. Reinforces the cost-accuracy tradeoff that is central to practical quantum chemistry. (LO14)

**Grading Notes:**
- 4 pts: Partial disagreement with 2+ practical considerations, well-reasoned
- 3 pts: Partial disagreement with 1 practical consideration
- 2 pts: Agrees but mentions cost as a caveat
- 1 pt: Simple agree/disagree with no reasoning
- 0 pts: No answer

---

## Section 2: Integral Inspection and Fock Tracing (Q2.1-Q2.8, 34 points)

### Q2.1: Predicting Dominant Off-Diagonal Matrix (4 points)

**Question Text:**
> For H2 (2 basis functions), which matrix do you predict will have the largest off-diagonal element (in absolute value): S, T, or V? Why?

**Expected Answer:**

**V (nuclear attraction)** has the largest off-diagonal element in absolute value.

**PySCF Reference Values for H2 STO-3G (R = 1.4 bohr):**

| Matrix | Off-diagonal element | Absolute value |
|--------|---------------------|---------------|
| S(1,2) | 0.6599 | 0.6599 |
| T(1,2) | 0.2370 | 0.2370 |
| V(1,2) | -1.1963 | **1.1963** |

V dominates because the nuclear attraction integral measures the Coulomb attraction of the electron density overlap (phi_1 * phi_2) with BOTH nuclei. Since both nuclei contribute attractive (negative) potential, the magnitude is large. The overlap integral S and kinetic energy integral T are significant but smaller.

**Acceptable Range:**
- Correct prediction of V with some reasoning about Coulomb attraction earns full credit
- Prediction of S is common and earns partial credit with good reasoning

**Common Misconceptions:**
1. Students predict S because "overlap" sounds largest
2. Students predict T because kinetic energy is a "big" quantity
3. Students do not realize V includes contributions from both nuclei

**Pedagogical Intent:**
POE-predict. Forces students to think about the physical meaning of each integral before seeing values. (LO15)

**Grading Notes:**
- 4 pts: Correct prediction (V) with reasoning about nuclear attraction and two-center contributions
- 3 pts: Correct prediction (V) with partial reasoning
- 2 pts: Incorrect prediction but well-reasoned (e.g., predicts S with good spatial argument)
- 1 pt: A prediction with minimal reasoning
- 0 pts: No prediction

---

### Q2.2: Identifying the Dominant Matrix and Why (4 points)

**Question Text:**
> Which matrix (S, T, or V) had the largest off-diagonal element for H2? Why does this matrix dominate?

**Expected Answer:**

**V (nuclear attraction)** has the largest off-diagonal element: V(1,2) = -1.1963 (absolute value 1.1963).

**Why V dominates:** The nuclear attraction integral V(1,2) = sum_A -Z_A integral phi_1(r) (1/r_A) phi_2(r) dr measures how much the overlap charge distribution (phi_1 * phi_2) is attracted to each nucleus. For H2:
- Both nuclei are hydrogen (Z = 1), contributing negative (attractive) potential
- The overlap charge distribution is significant because the atoms are bonded (close together)
- Each nucleus attracts the overlap density, so both contributions add (both are negative)

The kinetic energy off-diagonal element T(1,2) = 0.2370 is smaller because kinetic energy is a second-derivative operator, which is less sensitive to the spatial extent of the overlap. The overlap integral S(1,2) = 0.6599 is moderate because it simply measures the product of the two functions.

**Acceptable Range:**
- Must identify V as largest and give the approximate value
- Must provide a physical explanation

**Common Misconceptions:**
1. Students may be surprised that V is negative (attraction is a negative energy contribution)
2. Students may confuse the nuclear attraction integral with the nuclear repulsion energy

**Pedagogical Intent:**
POE-observe/explain. Builds physical intuition about why nuclear attraction integrals dominate the one-electron Hamiltonian. (LO15)

**Grading Notes:**
- 4 pts: Correct identification with physical explanation referencing Coulomb attraction
- 3 pts: Correct identification with partial explanation
- 2 pts: Correct identification but weak or absent explanation
- 1 pt: Wrong identification but reasonable value reporting
- 0 pts: No answer

---

### Q2.3: Predicting Effect of Distance on S(1,2) (4 points)

**Question Text:**
> If the H-H distance were doubled (from 1.4 to 2.8 bohr), predict how S(1,2) would change. Would it increase, decrease, or stay the same?

**Expected Answer:**

S(1,2) would **decrease** significantly.

**Reasoning:** The overlap integral S(1,2) = integral phi_1(r) phi_2(r) dr measures how much two basis functions "share space." As the atom centers move further apart, the product phi_1(r) * phi_2(r) becomes smaller at every point because each function is concentrated near its own nucleus. The overlap decreases approximately exponentially with increasing distance.

**PySCF can verify:** At R = 2.8 bohr, S(1,2) is approximately 0.32, compared to 0.66 at R = 1.4 bohr -- roughly halved.

**Limiting cases:**
- R -> 0: S(1,2) -> 1.0 (identical functions)
- R -> infinity: S(1,2) -> 0.0 (no overlap)

**Acceptable Range:**
- Must predict "decrease"
- Must provide reasoning based on spatial separation

**Common Misconceptions:**
1. Students predict "increase" because "more space between atoms means more room to overlap"
2. Students predict "stay the same" because the basis functions themselves do not change

**Pedagogical Intent:**
Tests understanding of the distance dependence of the overlap integral. Reinforces the spatial-symbolic bridge from Lab Pack #2. (LO15)

**Grading Notes:**
- 4 pts: Correct prediction (decrease) with reasoning about spatial separation and product integral
- 3 pts: Correct prediction with partial reasoning
- 2 pts: Correct prediction with minimal reasoning
- 1 pt: Incorrect prediction but shows spatial reasoning
- 0 pts: No prediction

---

### Q2.4: H2O Overlap Matrix Analysis (4 points)

**Question Text:**
> Identify the matrix elements corresponding to O-H bonded pairs and the H-H non-bonded pair. Record approximate values. Are these consistent with relative distances?

**Expected Answer:**

**Basis function mapping for H2O STO-3G:**

| Index | Basis Function | Atom |
|-------|---------------|------|
| 1 | O 1s | O |
| 2 | O 2s | O |
| 3 | O 2px | O |
| 4 | O 2py | O |
| 5 | O 2pz | O |
| 6 | H1 1s | H |
| 7 | H2 1s | H |

**Key overlap values (PySCF reference):**

| Atom Pair | Matrix Element | Value | Physical Interpretation |
|-----------|---------------|-------|------------------------|
| O(2s)-H1 (bonded) | S(2,6) | 0.4744 | Large: bonded atoms, short O-H distance |
| O(2py)-H1 (bonded) | S(4,6) | 0.3109 | Moderate: p-orbital overlap along bond |
| H1-H2 (non-bonded) | S(6,7) | 0.2515 | Smaller: non-bonded, larger distance |
| O(1s)-H1 (core-valence) | S(1,6) | 0.0539 | Very small: O 1s is very compact |

**Consistency with distances:** Yes. The O-H bond distance (~1.81 bohr) is shorter than the H-H distance (~2.86 bohr), so the O-H overlaps are larger. The O 1s overlap with H is tiny because the 1s orbital is extremely compact (tightly bound core).

**Acceptable Range:**
- Must identify at least one O-H bonded element (S(2,6), S(4,6), or S(5,6))
- Must identify S(6,7) as the H-H element
- Values within 0.05 of PySCF reference

**Common Misconceptions:**
1. Students confuse basis function indices with atom indices (very common, ~40-50%)
2. Students expect H-H overlap to be zero because the atoms are "not bonded"
3. Students may select diagonal elements as "largest" (those are always 1.0)

**Pedagogical Intent:**
Connects matrix indices to molecular geometry. Reinforces the spatial-symbolic bridge. (LO15)

**Grading Notes:**
- 4 pts: Both O-H and H-H correctly identified with approximate values and distance reasoning
- 3 pts: Both identified with values, minimal distance reasoning
- 2 pts: One correctly identified with value
- 1 pt: Attempts to read matrix but indices are confused
- 0 pts: No answer

---

### Q2.5: F = H^core + G(P) in Own Words (4 points)

**Question Text:**
> Write out the formula F = H^core + G(P) in your own words. What physical interactions does H^core capture? What does G(P) add?

**Expected Answer:**

**In words:** The Fock matrix is the sum of two contributions:

1. **H^core (core Hamiltonian)** captures the one-electron physics:
   - **Kinetic energy (T):** How much kinetic energy is associated with each pair of basis functions. This is an intrinsic property of the electron's motion.
   - **Nuclear attraction (V):** How strongly the nuclei attract electrons in each pair of basis functions. This depends on the positions of all nuclei and the spatial distribution of the basis functions.
   - H^core is **independent of the electron density** -- it depends only on the basis functions and nuclear positions.

2. **G(P) (two-electron matrix)** adds the electron-electron interactions, weighted by the density matrix P:
   - **Coulomb repulsion (J):** The classical repulsion between the electron density in one basis function pair and the density in all other pairs.
   - **Exchange interaction (K):** A quantum mechanical effect arising from the antisymmetry of the wavefunction. Exchange lowers the energy for electrons with parallel spin.
   - G(P) **depends on P** because the electron-electron interactions depend on where the electrons currently are (as described by the density matrix).

**Key insight:** Because G(P) depends on P, and P depends on the MOs, which depend on F, the equation is self-consistent -- hence "self-consistent field."

**Acceptable Range:**
- Must identify H^core as containing T and V (one-electron terms)
- Must identify G(P) as containing electron-electron terms
- Must note that G(P) depends on the density matrix

**Common Misconceptions:**
1. Students may describe H^core as "the Hamiltonian" without distinguishing one-electron from two-electron
2. Students may forget the exchange term in G(P)
3. Students may not understand why G depends on P

**Pedagogical Intent:**
Tests conceptual understanding of the Fock matrix decomposition. This is the central equation in HF theory. (LO16)

**Grading Notes:**
- 4 pts: Clear description of both H^core (T+V) and G(P) (Coulomb+exchange), notes P-dependence
- 3 pts: H^core and G(P) described, P-dependence noted, but exchange missing or unclear
- 2 pts: Both terms mentioned but description is superficial
- 1 pt: Only one term described
- 0 pts: No answer

---

### Q2.6: Why Rebuild F Every Iteration? (5 points)

**Question Text:**
> The G(P) matrix depends on P, which changes at each SCF iteration. Why must the Fock matrix be rebuilt every iteration? What would happen if you froze the density matrix?

**Expected Answer:**

**Why rebuild:** The Fock matrix F = H^core + G(P) depends on the density matrix P through the G(P) term. At each SCF iteration:
1. The density matrix P changes (because new MO coefficients are computed)
2. G(P) must be recomputed because it uses P to weight the electron repulsion integrals
3. A new Fock matrix F is formed
4. The new F is diagonalized to get new MO coefficients, which give a new P
5. The cycle repeats until P (and therefore F) stops changing -- self-consistency

H^core does NOT need to be recomputed because it depends only on the basis functions and nuclear positions, which are fixed during the SCF.

**If you froze the density matrix:** You would compute F = H^core + G(P_frozen) once and diagonalize it. The resulting MO coefficients would give a new P_new that is different from P_frozen. But since P_frozen is not updated, the F matrix would never change, and the energy would not improve beyond the first iteration. You would be computing the eigenvalues of an approximate Fock operator that is not self-consistent.

In effect, freezing P means solving a one-electron problem with fixed electron-electron interactions -- this is the "non-self-consistent" approximation, and it typically gives a poor energy.

**Acceptable Range:**
- Must explain that G(P) depends on P and P changes each iteration
- Must explain what happens if P is frozen (no self-consistency, poor energy)

**Common Misconceptions:**
1. Students may think the entire Fock matrix needs rebuilding (only G(P) changes; H^core is constant)
2. Students may think "convergence" means F is changing when it actually means F has STOPPED changing
3. Students may not connect freezing P to loss of self-consistency

**Pedagogical Intent:**
Core question for LO16. Tests understanding of the iterative nature of SCF and why self-consistency requires updating the Fock matrix. (LO16)

**Grading Notes:**
- 5 pts: Clear explanation of P-dependence + correct analysis of frozen-P scenario + notes H^core is constant
- 4 pts: P-dependence + frozen-P analysis, but does not note H^core is constant
- 3 pts: P-dependence explained, frozen-P analysis is partial
- 2 pts: Notes that F depends on P but does not analyze frozen-P scenario
- 1 pt: Vague statement about "the matrix changes"
- 0 pts: No answer

---

### Q2.7: Density Matrix Weighting of ERIs (5 points)

**Question Text:**
> In the Fock build trace, examine one ERI quartet that contributes to F(1,1). The contribution is weighted by P(lambda, sigma). How does the density matrix determine which ERI contributions matter most?

**Expected Answer:**

The two-electron contribution to the Fock matrix is:

G(mu,nu) = sum_{lambda,sigma} P(lambda,sigma) * [(mu nu | lambda sigma) - 0.5 * (mu lambda | nu sigma)]

The density matrix P acts as a **weighting function**: it determines how much each ERI contributes to the Fock matrix element.

**How P determines importance:**
- If P(lambda, sigma) is **large** (e.g., near 1.0), the corresponding ERI contributes significantly to G. These are pairs of basis functions where significant electron density resides.
- If P(lambda, sigma) is **near zero**, the corresponding ERI has almost no effect on G, regardless of how large the ERI value itself is.
- At the beginning of the SCF (when P comes from the initial guess), the weighting may be inaccurate. As P converges, the weighting becomes correct, and G accurately represents the electron-electron interactions.

**Example for H2 (PySCF):**
- P(1,1) = P(2,2) = P(1,2) = P(2,1) = 0.6025
- ERI (11|11) = 0.7746 -> Contribution: P(1,1) * 0.7746 = 0.467
- ERI (11|22) = 0.5700 -> Contribution: P(2,2) * 0.5700 = 0.343
- Both are significant because all density matrix elements are similar in magnitude for H2

**Acceptable Range:**
- Must explain that P weights the ERI contributions
- Must explain that small P means small contribution

**Common Misconceptions:**
1. Students think all ERIs contribute equally regardless of P
2. Students think the density matrix only appears in the energy expression, not in F
3. Students confuse the density matrix with the overlap matrix

**Pedagogical Intent:**
Tests understanding of how electron-electron interactions are self-consistently determined. The density matrix is the bridge between the current electronic state and the Fock operator. (LO16)

**Grading Notes:**
- 5 pts: Explains P as weighting function + explains large vs. small P effect + example or formula
- 4 pts: Explains P as weighting function + large vs. small P, no example
- 3 pts: Explains that P determines contributions but mechanism unclear
- 2 pts: Mentions P appears in the formula without explaining its role
- 1 pt: Vague mention of P
- 0 pts: No answer

---

### Q2.8: What Is Missing from "F = H^core + ERIs"? (4 points)

**Question Text:**
> A student says: "The Fock matrix is just H^core plus the electron repulsion integrals." What is missing from this statement?

**Expected Answer:**

Two things are missing:

1. **The density matrix weighting:** The ERIs are not added directly to H^core. They are weighted by the density matrix P. The correct statement is G(mu,nu) = sum_{lambda,sigma} P(lambda,sigma) [(mu nu | lambda sigma) - 0.5 (mu lambda | nu sigma)]. Without P, you would not know how much of each ERI to include.

2. **The exchange term:** The student's statement implies that only Coulomb repulsion is included. In the full Fock matrix, the two-electron contribution has two parts:
   - **Coulomb (J):** sum P(lambda,sigma) (mu nu | lambda sigma) -- classical electron-electron repulsion
   - **Exchange (K):** -0.5 * sum P(lambda,sigma) (mu lambda | nu sigma) -- quantum mechanical exchange interaction
   The exchange term lowers the energy for electrons with parallel spin and is essential for the correct description of electronic structure.

**Acceptable Range:**
- Must identify the density matrix weighting as missing
- Must identify the exchange term as missing (or at least note that ERIs enter as both Coulomb and exchange)

**Common Misconceptions:**
1. Students may think all that is needed is "just add the integrals" -- this ignores both P-weighting and exchange
2. Students may know about exchange but forget it enters with a -0.5 factor
3. Students may confuse exchange with electron correlation

**Pedagogical Intent:**
Tests precise understanding of Fock matrix construction. This misconception check ensures students do not oversimplify the G(P) term. (LO16)

**Grading Notes:**
- 4 pts: Identifies both missing elements (P-weighting AND exchange) with explanation
- 3 pts: Identifies both but explanation of one is weak
- 2 pts: Identifies one of the two missing elements clearly
- 1 pt: Vague answer about "something is missing" without specifying
- 0 pts: No answer or agrees with the student's statement

---

## Section 3: Electron Density and Difference Density (Q3.1-Q3.8, 32 points)

### Q3.1: Predicting H2O Density Distribution (3 points)

**Question Text:**
> Sketch where you expect the electron density to be highest for H2O. Consider: oxygen has 8 electrons and each hydrogen has 1. Where will the density concentrate?

**Expected Answer:**

The electron density should be **highest near the oxygen nucleus**. Oxygen has 8 electrons (compared to 1 on each hydrogen), and the core electrons (1s) are very tightly concentrated near the oxygen nucleus. The density peaks sharply at the nuclear positions, with the largest peak at oxygen.

**Expected sketch features:**
- Large density concentration at the oxygen position
- Smaller density concentrations at the hydrogen positions
- Some density in the O-H bonding regions (between O and each H)
- Very little density in the non-bonding region (opposite the O-H bonds)

**Common Misconceptions:**
1. Students predict even distribution because "the molecule has shared electrons"
2. Students forget about the core electrons on oxygen, which create a massive density peak
3. Students may think density is only in the bonding region

**Pedagogical Intent:**
POE-predict. Activates prior knowledge about electron count and electronegativity. (LO17)

**Grading Notes:**
- 3 pts: Predicts density highest near O with reasoning about electron count; sketch shows asymmetric distribution
- 2 pts: Predicts density near O but sketch is vague
- 1 pt: Predicts even distribution or density only on bonds
- 0 pts: No prediction

---

### Q3.2: Why Density Spikes Near Nuclei (4 points)

**Question Text:**
> In the 2D cross-section, there are sharp density peaks near the nuclear positions. Why does the electron density spike near the nuclei?

**Expected Answer:**

The electron density spikes near the nuclei because of the **strong Coulomb attraction** between the positively charged nuclei (Z = 8 for oxygen, Z = 1 for hydrogen) and the negatively charged electrons.

**Physical explanation:**
1. The Coulomb potential V = -Z/r diverges as r -> 0, creating an enormous attractive force near the nucleus.
2. Electrons in core orbitals (O 1s) are extremely tightly bound (orbital energy = -20.24 Ha) and concentrated in a very small volume near the oxygen nucleus.
3. Even valence orbitals have significant density near the nucleus because their wavefunctions (and thus density contributions) are largest where the potential is most attractive.
4. Mathematically, the density rho(r) = sum P_ij phi_i(r) phi_j(r) includes contributions from all basis functions, and the tight (large-exponent) primitives contribute enormous values near r = 0.

**The oxygen peak is much larger than the hydrogen peaks** because:
- Oxygen has 8 electrons vs. 1 per hydrogen
- The O 1s core orbital alone contributes approximately 2 electrons to a tiny volume
- The nuclear charge Z = 8 creates a much stronger Coulomb attraction than Z = 1

**Acceptable Range:**
- Must mention Coulomb attraction between electrons and nuclei
- Must note that core electrons are especially concentrated

**Common Misconceptions:**
1. Students may attribute the peaks to "nuclear charge" without connecting to electron attraction
2. Students may think the peaks are artifacts of the basis set rather than physical features
3. Students may not realize that the density at the nuclear position is actually infinite for a point nucleus (cusp condition)

**Pedagogical Intent:**
Tests understanding of why density accumulates near nuclei -- the fundamental Coulomb attraction that drives atomic and molecular structure. (LO17)

**Grading Notes:**
- 4 pts: Coulomb attraction explanation + core electron concentration + O peak > H peak reasoning
- 3 pts: Coulomb attraction + one of the other two points
- 2 pts: Mentions attraction but explanation is incomplete
- 1 pt: Notes peaks exist without physical explanation
- 0 pts: No answer

---

### Q3.3: Bonding vs. Non-Bonding Density Regions (4 points)

**Question Text:**
> Compare the density in the O-H bonding region to the density in the H-H non-bonding region. Is there more density in the bonding region? What does this tell you about bonding?

**Expected Answer:**

**Yes, there is more density in the O-H bonding region** than in the non-bonding region (the region opposite the O-H bonds, away from the hydrogen atoms).

**Observations from the cross-section:**
- Between oxygen and each hydrogen, there is a "bridge" of electron density connecting the nuclear peaks
- The H-H non-bonding region has much lower density
- The density between O and H is not as high as at the nuclear positions, but it is clearly elevated compared to other regions

**What this tells us about bonding:** Chemical bonds involve a concentration of electron density between the bonded atoms. This density:
- Experiences attractive Coulomb interaction with BOTH nuclei (it is in the potential well between them)
- Stabilizes the molecule by lowering the total energy
- Is the physical manifestation of "electron sharing" -- the electrons are not on one atom or the other but distributed between them

This is consistent with the bonding molecular orbitals seen in Lab Pack #2, which showed electron density concentrated between nuclei.

**Acceptable Range:**
- Must state that bonding region has more density
- Must connect density between nuclei to bonding/stabilization

**Common Misconceptions:**
1. Students may expect density to be uniform because "electrons are everywhere"
2. Students may confuse the density map with the orbital visualization from LP2
3. Students may not realize that the density between nuclei is what "glues" the molecule together

**Pedagogical Intent:**
Connects the density visualization to the concept of covalent bonding. Reinforces the geometry-energy connection from LO7 (LP2). (LO17)

**Grading Notes:**
- 4 pts: Correct observation (more in bonding region) + connection to Coulomb stabilization or bonding
- 3 pts: Correct observation + partial connection to bonding
- 2 pts: Correct observation without bonding connection
- 1 pt: Incorrect or vague comparison
- 0 pts: No answer

---

### Q3.4: Misconception Check -- Isovalue and Electron Count (4 points)

**Question Text:**
> When you decreased the isovalue, the isosurface expanded. Does this mean there are more electrons? Explain what the isovalue threshold actually represents.

**Expected Answer:**

**No, there are NOT more electrons.** The number of electrons is fixed (10 for H2O) and does not change when you adjust the isovalue.

**What actually happens:** The isovalue is a threshold on the density value rho(r). The isosurface connects all points where rho(r) = isovalue.

- **Higher isovalue (e.g., 0.10):** The surface encloses only the high-density regions close to the nuclei. The enclosed volume is small but the density within it is high.
- **Lower isovalue (e.g., 0.01):** The surface extends further out to include regions of lower density. The enclosed volume is larger, but the density at the surface itself is lower.
- **The electron density extends beyond any finite isovalue surface.** The surface is a chosen threshold, not a physical boundary.

**Analogy:** Think of a topographic map. Lowering the elevation contour line does not add more mountain -- it just shows more of the terrain that was always there. Similarly, lowering the isovalue does not add electrons -- it reveals more of the density distribution that was always present.

**Acceptable Range:**
- Must state "no, same number of electrons"
- Must explain isovalue as a density threshold
- Must note that density extends beyond the surface

**Common Misconceptions:**
1. "More surface = more electrons" -- confusing enclosed volume with electron count
2. "The isovalue controls the number of electrons included" -- partially true in the sense of enclosed probability, but misleading
3. "Electrons can only exist inside the isosurface" -- the orbital boundary misconception from LP2

**Pedagogical Intent:**
Directly targets the misconception that isosurface size correlates with electron count. Reinforces the isovalue concept from LP2 (LO12) in the new context of density visualization. (LO17)

**Grading Notes:**
- 4 pts: Correctly states "no" + explains isovalue as threshold + notes density extends beyond surface
- 3 pts: Correctly states "no" + threshold explanation, omits density extension
- 2 pts: States "no" but explanation is vague
- 1 pt: Says "yes, more electrons" but shows some understanding of isovalue
- 0 pts: Says "yes" with no qualification or no answer

---

### Q3.5: Predicting H2 Difference Density Pattern (4 points)

**Question Text:**
> When two hydrogen atoms form H2, electrons rearrange. Predict: where will electrons accumulate (compared to the promolecule)? Where will they deplete? Sketch your prediction.

**Expected Answer:**

**Accumulation (Delta-rho > 0):** Electrons should accumulate in the **bonding region between the two nuclei**. When atoms form a covalent bond, electron density shifts from around each individual atom toward the internuclear region. This is the physical basis of covalent bonding.

**Depletion (Delta-rho < 0):** Electrons should deplete from the **outer regions** of each atom -- the regions behind each nucleus (away from the bond). The density that moves into the bonding region must come from somewhere, and it comes from the non-bonding side of each atom.

**Expected sketch:** An elongated accumulation region (solid) between the two H nuclei, with two depletion regions (translucent) on the outside of each atom, forming a "dumbbell" pattern with accumulation in the middle.

**Acceptable Range:**
- Must predict accumulation between the atoms
- Must predict depletion from the outer regions
- Sketch should show the correct qualitative pattern

**Common Misconceptions:**
1. Students predict accumulation at the nuclei (that is where total density is highest, but the DIFFERENCE density is zero or slightly depleted there because the promolecule already has nuclear density)
2. Students predict no change because "it is still H atoms"
3. Students predict depletion between the atoms (opposite of the correct answer)

**Pedagogical Intent:**
POE-predict for difference density. Forces students to think about how bonding rearranges electrons before seeing the visualization. (LO18)

**Grading Notes:**
- 4 pts: Correctly predicts both accumulation (between) and depletion (outside) with reasoning about bonding
- 3 pts: Correctly predicts accumulation between atoms, depletion is incorrect or absent
- 2 pts: Predicts some rearrangement but location is wrong
- 1 pt: Predicts no change or vague answer
- 0 pts: No prediction

---

### Q3.6: Why Use the Promolecule as Reference? (4 points)

**Question Text:**
> What is a promolecule, and why is it the right reference for computing the difference density?

**Expected Answer:**

**What is a promolecule:** The promolecule density is the sum of spherically-averaged free-atom densities placed at the molecular geometry positions, WITHOUT allowing them to interact. It represents what the electron density would look like if the atoms were at their molecular positions but not forming bonds.

**Why it is the right reference:**
1. **Isolates bonding effects:** By subtracting the promolecule (non-interacting atoms at molecular positions), the difference density Delta-rho shows ONLY the changes due to chemical bonding. If we subtracted isolated atoms at infinite separation, we would also see changes due to simple overlap (bringing atoms close without bonding), which would obscure the bonding signal.

2. **Physical meaning:** Delta-rho > 0 means electrons have accumulated in that region DUE TO bonding interactions. Delta-rho < 0 means electrons have been depleted from that region and redistributed elsewhere DUE TO bonding.

3. **Zero reference:** At large interatomic distances where bonding is negligible, the molecular density should approach the promolecule density, giving Delta-rho approximately zero -- this is physically correct.

**Why NOT isolated atoms at infinite separation:**
Using infinitely separated atoms as the reference would include the trivial density change from "atoms are now nearby" (overlap of atomic densities) mixed with the chemically interesting change from "atoms are now bonding." The promolecule cleanly separates these effects.

**Acceptable Range:**
- Must define promolecule as non-interacting atoms at molecular positions
- Must explain why this reference isolates bonding effects

**Common Misconceptions:**
1. Students think the promolecule is "the molecule before bonds form" -- it is a theoretical construct, not a physical intermediate
2. Students may not distinguish between "atoms nearby" and "atoms bonding"
3. Students may think any reference would give the same result

**Pedagogical Intent:**
Tests understanding of the reference frame for interpreting difference density maps. Essential for correct interpretation of bonding-induced charge redistribution. (LO18)

**Grading Notes:**
- 4 pts: Correct promolecule definition + explains why it isolates bonding + contrasts with infinite-separation reference
- 3 pts: Correct definition + isolation reasoning, no contrast
- 2 pts: Correct definition but purpose is unclear
- 1 pt: Vague or incorrect definition
- 0 pts: No answer

---

### Q3.7: Comparing H2 and H2O Difference Density (5 points)

**Question Text:**
> Compare the H2 and H2O difference density maps. Where does charge accumulate in H2O? Do you observe asymmetry in the accumulation along the O-H bonds?

**Expected Answer:**

**H2 difference density:** Symmetric accumulation between the two H nuclei (in the bonding region), with symmetric depletion on the outer sides of each atom. The pattern is perfectly symmetric because the two atoms are identical.

**H2O difference density:**
- **Accumulation:** Charge accumulates in the O-H bonding regions (between O and each H), similar to H2 but with important differences.
- **Asymmetry:** The accumulation is NOT symmetric along the O-H bond. It is shifted toward the oxygen atom. This reflects the higher electronegativity of oxygen -- oxygen attracts electron density more strongly than hydrogen.
- **Oxygen lone pair regions:** There may be additional accumulation in the regions corresponding to oxygen's lone pairs (perpendicular to the molecular plane and behind the oxygen relative to the H atoms).
- **Hydrogen depletion:** The regions behind each hydrogen (away from the bond) show depletion, similar to H2 but more pronounced because oxygen pulls electron density more strongly.

**Comparison:**
- H2: symmetric accumulation because identical atoms
- H2O: asymmetric accumulation because O is more electronegative than H
- Both show accumulation in bonding regions and depletion in non-bonding regions

**Acceptable Range:**
- Must describe accumulation in O-H bonding regions for H2O
- Must note the asymmetry (shifted toward O) and connect to electronegativity

**Common Misconceptions:**
1. Students expect symmetric accumulation in O-H bonds (ignoring electronegativity)
2. Students think the difference density shows total electron count rather than redistribution
3. Students confuse the lone pair accumulation with bonding accumulation

**Pedagogical Intent:**
Applies difference density interpretation to a polar molecule. The asymmetry demonstrates electronegativity effects, connecting computational output to general chemistry concepts. (LO18)

**Grading Notes:**
- 5 pts: Describes both molecules, notes asymmetry in H2O, connects to electronegativity, systematic comparison
- 4 pts: H2O accumulation described with asymmetry noted, but comparison to H2 is superficial
- 3 pts: Accumulation in bonding regions described without noting asymmetry
- 2 pts: Partial description without clear comparison
- 1 pt: Vague observations
- 0 pts: No answer

---

### Q3.8: Misconception Check -- Difference Density Regions (4 points)

**Question Text:**
> A student says: "The solid regions in the difference density map contain all the bonding electrons; the translucent regions contain none." Correct this statement.

**Expected Answer:**

**This statement is incorrect** on both counts.

**Correction:**
1. The solid regions (Delta-rho > 0) do NOT "contain all the bonding electrons." They show where the electron density has **increased** relative to the promolecule. There were already electrons in these regions in the promolecule -- the solid region shows the ADDITIONAL density that appeared due to bonding. Quantitatively, these regions typically account for a fraction of an electron (e.g., 0.1-0.5 electrons for a single bond), not all the bonding electrons.

2. The translucent regions (Delta-rho < 0) are NOT empty. They show where the electron density has **decreased** relative to the promolecule. There are still electrons in these regions -- just fewer than in the promolecule. The density has been redistributed, not removed entirely.

**Key insight:** The difference density shows CHANGES in density, not absolute density values. Both solid and translucent regions still contain electrons. The difference density integrates to zero over all space (total electron count is conserved -- electrons are redistributed, not created or destroyed).

**Acceptable Range:**
- Must correct both parts of the statement
- Must explain that difference density shows changes, not absolute values
- Must note electron conservation (density redistributed, not created/destroyed)

**Common Misconceptions:**
1. Students treat solid regions as "where electrons are" and translucent as "where they are not"
2. Students do not realize the difference density integrates to zero
3. Students confuse difference density with total density

**Pedagogical Intent:**
Directly targets the most common misconception about difference density maps. Ensures students understand that Delta-rho represents redistribution, not absolute location. (LO18)

**Grading Notes:**
- 4 pts: Corrects both parts + explains "changes not absolute" + mentions conservation/redistribution
- 3 pts: Corrects both parts + explains "changes not absolute"
- 2 pts: Corrects one part clearly
- 1 pt: Recognizes something is wrong but correction is vague
- 0 pts: Agrees with the student or no answer

---

## Section 4: Synthesis (Q4.1-Q4.2, 10 points)

### Q4.1: Tracing the Computational Path (5 points)

**Question Text:**
> Suppose you changed the basis set from STO-3G to 6-31G for H2O. At each layer, describe one specific thing that would change.

**Expected Answer:**

**Basis functions:** The number of basis functions increases from 7 to 13. Each valence shell (O 2s, O 2p, H 1s) is now represented by TWO contractions instead of one (split-valence), and the O core shell uses 6 primitives instead of 3. The radial profiles extend further from the nuclei due to the additional diffuse functions.

**Integrals:** The integral matrices grow from 7x7 to 13x13. The number of unique two-electron integrals increases from approximately (7^4)/8 = 300 to (13^4)/8 = 3,570 (approximately 12x more). The overlap matrix will have more nonzero elements because there are more basis function pairs. Individual integral values will also change because the basis functions have different exponents and coefficients.

**Fock matrix / SCF result:** The Fock matrix grows from 7x7 to 13x13. The SCF energy decreases (becomes more negative) from -74.963 Ha to -75.984 Ha -- an improvement of about 1.02 Ha. This is guaranteed by the variational principle. The MO coefficients and orbital energies will change. The SCF may converge in a different number of iterations.

**Electron density:** The density will be more accurate -- the additional variational freedom allows the SCF to find a better representation of the electron distribution. Differences will be most noticeable in the valence/bonding regions (where the split-valence functions provide extra flexibility). The core density near the oxygen nucleus will be largely unchanged because core electrons are insensitive to bonding.

**PySCF Reference:**
- H2O STO-3G (7 functions): -74.9630 Ha
- H2O 6-31G (13 functions): -75.9840 Ha

**Acceptable Range:**
- Must provide a specific change at each of the four layers
- Numerical values welcome but not required
- Must convey the cascade: basis set choice -> integral computation -> SCF result -> density

**Common Misconceptions:**
1. Students may think only the energy changes, not the density or integrals
2. Students may not realize the integral count scales as N^4
3. Students may think the core density changes significantly (it does not)

**Pedagogical Intent:**
Synthesis question connecting all three sections. Tests understanding of the computational pipeline. (LO13, LO14, LO15, LO16, LO17)

**Grading Notes:**
- 5 pts: Specific, correct change at all 4 layers with cascade logic
- 4 pts: Correct changes at 3 of 4 layers
- 3 pts: Correct changes at 2 of 4 layers
- 2 pts: Correct change at 1 layer
- 1 pt: Vague changes without specifics
- 0 pts: No answer

---

### Q4.2: Which Representation Best Supported Understanding? (5 points)

**Question Text:**
> You explored mathematical, graphical, and spatial representations. Which was most useful for understanding bonding, and why? Was there a concept that only became clear in a particular representation?

**Expected Answer:**

This is an open-ended reflection question. There is no single correct answer, but strong responses will:

1. **Identify a specific representation** and explain why it was most useful for a specific concept:
   - Mathematical (exponent tables, matrix elements): Best for understanding the precision and quantitative structure of the computation
   - Graphical (radial profiles, heatmaps, cross-sections): Best for seeing patterns and comparisons
   - Spatial (3D isosurfaces): Best for connecting to molecular geometry and spatial intuition

2. **Give a specific example** of a concept that became clear through a particular representation:
   - "I understood contraction when I saw the radial profile with the three primitives" (graphical)
   - "The difference density isosurface was the first time I really 'saw' electron redistribution in bonding" (spatial)
   - "Looking at the overlap matrix values helped me understand why some bonds are stronger" (mathematical)
   - "The Fock build trace with step-by-step matrices made the self-consistency loop click" (graphical/mathematical)

3. **Articulate the value of multiple representations:**
   - No single representation captures everything
   - Mathematical gives precision; graphical gives patterns; spatial gives intuition
   - Together they provide complementary perspectives

**Acceptable Range:**
- Must name at least one representation and explain its value for a specific concept
- Must give a concrete example (not just "the graphical one was good")

**Common Misconceptions:**
1. Students list all three without articulating what each uniquely contributes
2. Students focus on aesthetics ("the 3D view looked cool") rather than understanding

**Pedagogical Intent:**
Metacognitive reflection on representational competence (Kozma & Russell, 2005). Assesses whether students can articulate how different representations support understanding of the same underlying physics. (Integrative)

**Grading Notes:**
- 5 pts: Names representation + explains value for specific concept + concrete example + articulates complementarity
- 4 pts: Names representation + explains value + concrete example
- 3 pts: Names representation + partial explanation
- 2 pts: Lists representations without clear explanation
- 1 pt: Vague response
- 0 pts: No answer

---

## Common Student Errors and Remediation

### Section 1: Basis Function Exploration

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Expecting three separate bumps from contraction** | Common (25-35%) | Show the radial profile with primitive decomposition. Ask: "Are all the coefficients positive? If so, can the sum of positive functions ever have separate bumps?" The answer is no -- the primitives add constructively everywhere. |
| **Confusing exponent with coefficient** | Common (30-40%) | Clarify: exponent controls WIDTH (large exponent = tight/narrow); coefficient controls how much of each primitive is included. Use the exponent slider to demonstrate. |
| **Thinking basis function = atomic orbital** | Very common (40-50%) | This is the target misconception for Q1.4. Emphasize: basis functions are computational inputs (mathematical tools); atomic orbitals are physical solutions to the Schrodinger equation. The former approximate the latter. |
| **Thinking "more basis functions = always better" without qualification** | Common (30-40%) | Discuss N^4 scaling. Show the PySCF reference: going from 7 to 13 basis functions increases ERIs by 12x. Ask: "Would you rather wait 12 seconds or 1 second?" |

### Section 2: Integral Inspection and Fock Tracing

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Confusing basis function indices with atom indices** | Very common (40-50%) | Write the mapping on the board: functions 1-5 = O, 6 = H1, 7 = H2. Have students label the matrix axes before reading values. |
| **Thinking S, T, V are equally important** | Moderate (20-30%) | Point to the PySCF values. V(1,2) = -1.20 is twice as large as S(1,2) = 0.66. Ask: "Which physical interaction is strongest for electrons near two nuclei?" |
| **Thinking F is computed once (not iteratively)** | Common (25-35%) | Step through the Fock build trace twice: once with the initial P, once with the converged P. Show how G(P) changes. Ask: "If P changes, does F stay the same?" |
| **Omitting exchange from the Fock matrix** | Common (30-40%) | Write the formula: G = sum P [(mu nu \| lambda sigma) - 0.5 (mu lambda \| nu sigma)]. Point to the -0.5 term. Explain: "This is the quantum mechanical contribution that distinguishes HF from classical electrostatics." |
| **Confusing density matrix P with overlap matrix S** | Moderate (15-25%) | Side-by-side comparison: S is symmetric with diagonal = 1.0 and depends only on basis functions (geometry). P is symmetric with diagonal != 1.0 in general and depends on the electronic state (MO coefficients). |

### Section 3: Electron Density and Difference Density

| Error | Frequency | Remediation |
|-------|-----------|-------------|
| **Thinking isovalue change creates/destroys electrons** | Common (30-40%) | Ask: "Did you add or remove electrons? No -- you just changed the threshold at which the surface is drawn. The density is there whether or not you draw a surface around it." |
| **Confusing difference density with total density** | Very common (40-50%) | Emphasize: Delta-rho = rho_mol - rho_pro. The solid/translucent surfaces show CHANGES, not where electrons ARE. Ask: "Is there zero density in the translucent regions?" (No -- there is just LESS density than the promolecule.) |
| **Treating solid regions as "where all bonding electrons are"** | Common (25-35%) | The solid regions typically contain a fraction of an electron (0.1-0.5 e), not all the bonding electrons. The difference density integrates to zero. Show: solid = redistribution IN, translucent = redistribution OUT. |
| **Not recognizing asymmetry in polar bond difference density** | Moderate (20-30%) | Compare H2 (symmetric) to H2O (asymmetric). Ask: "Oxygen is more electronegative. Where should the accumulation shift?" |

---

## Performance Task Rubrics

### Basis Set Analysis Rubric (4-point scale)

**Covers:** LO13 (contracted Gaussian anatomy), LO14 (basis set comparison)

**Task description:** Given radial profile data for a contracted Gaussian basis function and a comparison between two basis sets, students analyze the structure of the contraction, explain the role of different primitives, and evaluate the tradeoff between basis set size and computational cost.

| Score | Level | Criteria |
|-------|-------|----------|
| **4** | **Exemplary** | **Anatomy:** Correctly explains contraction as a fixed linear combination of primitives; identifies the role of tight (near-nucleus), medium, and diffuse (tail) primitives. Distinguishes basis functions from atomic orbitals. **Comparison:** Correctly explains what "split-valence" means -- two independent contractions for the valence. Explains why splitting provides additional variational freedom (the SCF can independently weight each contraction). Evaluates the cost-accuracy tradeoff quantitatively (e.g., references N^4 scaling or specific energy differences). |
| **3** | **Proficient** | **Anatomy:** Correctly explains contraction and identifies at least two primitive roles. Recognizes basis functions as approximations. **Comparison:** Identifies that split-valence provides more functions with some explanation of flexibility. Mentions cost-accuracy tradeoff but without quantitative detail. |
| **2** | **Developing** | **Anatomy:** Understands that contraction combines multiple Gaussians but cannot explain the roles of different primitives. May confuse basis functions with atomic orbitals. **Comparison:** Notes that 6-31G has "more functions" but cannot explain what the split accomplishes. Cost-accuracy reasoning is absent or incorrect. |
| **1** | **Beginning** | **Anatomy:** Cannot explain contraction or confuses it with other concepts (e.g., thinks each primitive is a separate basis function). Cannot distinguish basis functions from atomic orbitals. **Comparison:** No meaningful comparison between basis sets, or incorrect conclusions (e.g., "STO-3G is better because it is simpler"). |

**Scoring notes:**
- Anatomy and comparison contribute roughly equally
- The misconception check answer in Q1.4 is a strong diagnostic for the anatomy score level
- Students who can explain the exponent-width relationship earn at least "Developing"

---

### Integral Interpretation Rubric (4-point scale)

**Covers:** LO15 (S, T, V physical meaning), LO16 (Fock matrix construction)

**Task description:** Given integral matrix values and a Fock build trace, students interpret the physical meaning of one-electron integrals, explain how the Fock matrix is assembled from H^core and G(P), and describe the role of the density matrix in weighting electron repulsion contributions.

| Score | Level | Criteria |
|-------|-------|----------|
| **4** | **Exemplary** | **Integral meaning:** Correctly explains S (spatial overlap), T (kinetic energy), and V (nuclear attraction) with physical reasoning for relative magnitudes. Predicts how integrals change with geometry. **Fock construction:** Correctly traces F = H^core + G(P); explains that H^core is one-electron (fixed) and G(P) is two-electron (density-dependent). Identifies both Coulomb and exchange contributions. Explains why F must be rebuilt each iteration (self-consistency). Describes the density matrix as a weighting function for ERIs. |
| **3** | **Proficient** | **Integral meaning:** Correctly identifies the physical meaning of at least two of S, T, V. Predicts distance dependence of overlap correctly. **Fock construction:** Traces F = H^core + G(P) with correct identification of each component's physical content. Mentions density-dependence but may not fully explain Coulomb/exchange distinction or iterative rebuilding. |
| **2** | **Developing** | **Integral meaning:** Can read matrix values but struggles to connect them to physical meaning. Distance dependence prediction may be incorrect. **Fock construction:** Knows F is built from components but cannot clearly distinguish H^core from G(P) or explain why iteration is needed. May omit exchange. |
| **1** | **Beginning** | **Integral meaning:** Cannot identify the physical meaning of the integral matrices. May confuse S, T, V with each other. **Fock construction:** Cannot explain how F is assembled. Does not understand the iterative nature of SCF or the role of the density matrix. |

**Scoring notes:**
- Integral meaning and Fock construction contribute roughly equally
- Q2.8 (missing exchange/density weighting) is a strong diagnostic for the Fock construction score
- Students who correctly trace the Fock build steps earn at least "Developing"

---

## Expected Outputs for Exercises

### Section 1 Expected Radial Profiles

**H 1s STO-3G Radial Profile:**
- Three dashed lines showing the three primitives (tight, medium, diffuse)
- Tight primitive (alpha = 3.425): narrow peak near r = 0, small coefficient
- Medium primitive (alpha = 0.624): moderate width, largest coefficient
- Diffuse primitive (alpha = 0.169): wide, extends to r > 3 bohr
- Solid line: contracted function -- single smooth curve, sharper than any individual Gaussian near r = 0, with extended tail

**O 2s Comparison (STO-3G vs. 6-31G):**
- STO-3G: single contraction with moderate extent
- 6-31G: two components visible -- inner (3-primitive) contraction concentrated near nucleus, outer (1-primitive, alpha = 0.270) function extending much further
- The 6-31G profile should visibly extend further from the nucleus

### Section 2 Expected Integral Matrices

**H2 STO-3G Matrices:**

| Matrix | (1,1) | (1,2) | (2,2) |
|--------|-------|-------|-------|
| S | 1.0000 | 0.6599 | 1.0000 |
| T | 0.7600 | 0.2370 | 0.7600 |
| V | -1.8810 | -1.1963 | -1.8810 |
| H^core | -1.1210 | -0.9594 | -1.1210 |
| G(P) | 0.7549 | 0.3651 | 0.7549 |
| F | -0.3660 | -0.5943 | -0.3660 |

**H2O STO-3G Overlap Matrix:** 7x7 symmetric matrix with block structure. Key off-diagonal elements: S(2,6) = 0.4744, S(4,6) = 0.3109, S(6,7) = 0.2515.

### Section 3 Expected Density Visualizations

**H2O Total Density Isosurface:**
- Elongated surface with largest extent around oxygen
- Smaller "bumps" at hydrogen positions
- Overall shape reflects the bent molecular geometry

**H2O Density Cross-Section (XZ plane):**
- Sharp peaks at nuclear positions (O peak much larger than H peaks)
- Density bridges between O and each H (bonding regions)
- Lower density in the non-bonding region opposite the H atoms

**H2 Difference Density:**
- Solid (accumulation): elongated region between the two nuclei
- Translucent (depletion): two regions on the outer sides of each atom
- Symmetric because both atoms are identical

---

## Point Allocation Detail

### Per-Question Point Breakdown

| Question | Points | LO | Bloom's Level |
|----------|--------|-----|---------------|
| Q1.1 | 3 | LO13 | Apply |
| Q1.2 | 4 | LO13 | Analyze |
| Q1.3 | 4 | LO13 | Analyze |
| Q1.4 | 4 | LO13 | Analyze |
| Q1.5 | 3 | LO14 | Apply |
| Q1.6 | 4 | LO14 | Analyze |
| Q1.7 | 4 | LO14 | Analyze |
| Q1.8 | 4 | LO14 | Evaluate |
| **Section 1 Total** | **30** | | |
| Q2.1 | 4 | LO15 | Apply |
| Q2.2 | 4 | LO15 | Analyze |
| Q2.3 | 4 | LO15 | Apply |
| Q2.4 | 4 | LO15 | Apply |
| Q2.5 | 4 | LO16 | Understand |
| Q2.6 | 5 | LO16 | Analyze |
| Q2.7 | 5 | LO16 | Analyze |
| Q2.8 | 4 | LO16 | Analyze |
| **Section 2 Total** | **34** | | |
| Q3.1 | 3 | LO17 | Apply |
| Q3.2 | 4 | LO17 | Analyze |
| Q3.3 | 4 | LO17 | Analyze |
| Q3.4 | 4 | LO17 | Analyze |
| Q3.5 | 4 | LO18 | Apply |
| Q3.6 | 4 | LO18 | Analyze |
| Q3.7 | 5 | LO18 | Analyze |
| Q3.8 | 4 | LO18 | Analyze |
| **Section 3 Total** | **32** | | |
| Q4.1 | 5 | Integrative | Analyze |
| Q4.2 | 5 | Integrative | Evaluate |
| **Section 4 Total** | **10** | | |
| **Grand Total** | **106** | | |

### Cognitive Level Distribution (Worksheet)

| Bloom's Level | Questions | Points | Percentage |
|---------------|-----------|--------|------------|
| Understand | Q2.5 | 4 | 3.8% |
| Apply | Q1.1, Q1.5, Q2.1, Q2.3, Q2.4, Q3.1, Q3.5 | 25 | 23.6% |
| Analyze | Q1.2, Q1.3, Q1.4, Q1.6, Q1.7, Q2.2, Q2.6, Q2.7, Q2.8, Q3.2, Q3.3, Q3.4, Q3.6, Q3.7, Q3.8, Q4.1 | 67 | 63.2% |
| Evaluate | Q1.8, Q4.2 | 9 | 8.5% |
| **Subtotals** | | | |
| Remember/Understand | | 4 | **3.8%** |
| Apply/Analyze | | 92 | **86.8%** |
| Evaluate | | 9 | **8.5%** |

**Note:** The cognitive level distribution heavily favors Apply/Analyze (86.8%) because the lab emphasizes hands-on exploration and interpretation through the POE framework. The Evaluate level (8.5%) captures the critical evaluation questions (Q1.8 and Q4.2). This distribution is appropriate for a lab activity where students are building understanding through guided inquiry rather than recalling facts.

### Converting to Course Grade

| Raw Score (out of 106) | Percentage | Suggested Grade |
|------------------------|------------|-----------------|
| 96-106 | 90-100% | A |
| 85-95 | 80-89% | B |
| 75-84 | 70-79% | C |
| 64-74 | 60-69% | D |
| 0-63 | < 60% | F |

---

## Timing Guidance

### Standard Pacing (60 minutes)

| Activity | Time | Notes |
|----------|------|-------|
| Section 1: Basis Function Exploration | 12-15 min | Q1.1-Q1.8, includes radial profile inspection |
| Section 2: Integral Inspection & Fock Tracing | 18-22 min | Q2.1-Q2.8, includes Fock build stepping |
| Section 3: Electron Density & Difference Density | 22-25 min | Q3.1-Q3.8, includes density rendering time |
| Section 4: Synthesis | 3-5 min | Q4.1-Q4.2 |
| **Total** | **55-67 min** | Target: 60 min |

### Adjustments by Class Level

**Introductory (general chemistry level):**
- Add 5-10 min buffer
- Simplify Q1.7 (split-valence) to a brief explanation
- Consider making Q2.7 (density matrix weighting) and Q2.8 (misconception check on exchange) take-home
- Provide more explicit basis function index mapping
- Target: 65-75 min

**Intermediate (physical chemistry level):**
- Standard pacing as above
- Expect quantitative answers for integral values
- Target: 55-65 min

**Advanced (graduate level):**
- Expect deeper analysis (e.g., discuss linear dependence in Q1.8)
- May add discussion of exchange integral interpretation
- Ask students to compute difference density integrals mentally
- Target: 50-60 min

### Technology Setup Requirements

- Computer lab with modern web browsers (Chrome or Firefox preferred)
- Internet access for iqcp.dev (or local deployment)
- IQCP Modules A, B, and E (density tab) must be functional
- Projector for instructor to demonstrate Module A/B controls if needed

---

## FAQ for Common Issues

### Module A (Basis Explorer) Issues

**"The radial profile is empty or does not show primitives"**
- Verify the correct element and basis set are selected
- Click on a shell row in the shell table to select it -- the profile appears for the selected shell
- The primitive decomposition (dashed lines) appears alongside the contracted function (solid line)

**"The comparison mode does not show two profiles"**
- Toggle comparison mode ON in the controls
- Select a second basis set from the comparison dropdown
- Both basis sets must support the selected element (e.g., both must have oxygen)

**"The exponent slider does not seem to do anything"**
- Move the slider by a larger amount to see a visible change
- The radial profile updates in real-time -- look for the dashed line corresponding to the modified primitive
- Some primitives (especially the diffuse ones) may show subtle changes

### Module B (Integral Inspector) Issues

**"The integral matrix shows all zeros"**
- The integrals require an SCF calculation to be complete. Select a pre-computed system from the dropdown, or click "Run SCF" first
- Check the browser console for errors if the calculation does not start

**"I cannot find the correct matrix element in the heatmap"**
- Remember the basis function indexing: for H2O STO-3G, functions 1-5 are on oxygen, 6-7 on hydrogen
- Click directly on a cell to see its value in the detail panel
- Use the axis labels to identify rows and columns

**"The Fock build trace does not show steps"**
- Use the numbered step buttons or step slider in the Fock Build panel
- Step 1: H^core = T + V; Step 2: G(P); Step 3: F = H^core + G(P)
- If the steps are not visible, try expanding the panel or scrolling down

### Density Visualization Issues

**"The density isosurface is invisible or very small"**
- Try decreasing the isovalue from the default (e.g., from 0.05 to 0.02)
- Ensure the SCF calculation has converged (green indicator)
- The density tab may need a moment to compute the grid

**"The difference density shows nothing"**
- Difference density values are much smaller than total density. Use a lower isovalue (e.g., 0.005 or 0.002)
- Ensure "Difference density" mode is selected (not "Total density")
- Try rotating the molecule to see the accumulation/depletion regions from different angles

**"I cannot tell which regions are accumulation vs. depletion"**
- Solid surfaces = accumulation (Delta-rho > 0)
- Translucent surfaces = depletion (Delta-rho < 0)
- If colors are hard to distinguish, try rotating -- the transparency difference is usually visible from multiple angles

---

## Discussion Prompts by Section

### Section 1 Discussion Prompts

1. **After Q1.3:** "Gaussians have no cusp at r = 0 but Slater orbitals do. Why does this matter for the energy? Which integral is most affected by the cusp?" (Answer: The kinetic energy integral and the nuclear attraction integral are most affected. The cusp produces a finite kinetic energy contribution that Gaussians cannot reproduce exactly.)

2. **After Q1.4:** "If basis functions are not atomic orbitals, then what ARE the molecular orbitals? How do they relate to both?" (MOs are linear combinations of basis functions, optimized by the SCF. Basis functions approximate AOs; MOs are built from basis functions but span the whole molecule.)

3. **After Q1.8:** "A researcher wants to compute energies to 'chemical accuracy' (1 kcal/mol = 0.0016 Ha). Is the STO-3G to 6-31G improvement (1 Ha) enough? What about 6-31G to cc-pVTZ?" (Motivates the concept of basis set convergence.)

### Section 2 Discussion Prompts

1. **After Q2.2:** "The nuclear attraction integral V(1,2) is negative. Is that good or bad for bonding? What would happen if V(1,2) were zero?" (Negative V means electrons in the overlap region are attracted to nuclei -- this stabilizes the molecule. Zero V would mean no nuclear attraction of overlap density, weakening the bond.)

2. **After Q2.6:** "DIIS (from Lab Pack #1) extrapolates the Fock matrix from previous iterations. How does this connect to what you just learned about F depending on P?" (DIIS accelerates convergence by predicting the converged F without fully rebuilding G(P) each time.)

3. **After Q2.8:** "In density functional theory (DFT), what replaces the exchange term?" (The exchange-correlation functional -- a preview for students going on to study DFT.)

### Section 3 Discussion Prompts

1. **After Q3.2:** "The electron density is highest at the nuclei. But the probability of finding an electron at exactly one point is zero. How do you reconcile these statements?" (rho(r) is a probability density, not a probability. The probability of finding an electron in a volume dV is rho(r) dV. The density can be large but the probability in an infinitesimal volume is still infinitesimal.)

2. **After Q3.5:** "The difference density for H2 shows accumulation between the atoms. Is this the cause of bonding, or a consequence of it?" (This is a deep question. Modern understanding: the kinetic energy lowering from delocalization is the primary driving force; the charge accumulation is a consequence of the orbital optimization.)

3. **After Q3.8:** "Could you determine bond strength from the difference density alone?" (Qualitatively yes -- more accumulation typically correlates with stronger bonds. Quantitatively, you would need to integrate the density and compute energy contributions.)

---

*Lab Pack #3 Instructor Key v1.0 | CONFIDENTIAL -- Instructor Use Only*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
