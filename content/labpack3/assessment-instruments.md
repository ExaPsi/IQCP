# Lab Pack #3: Assessment Instruments

**Lab Pack:** 3 - Computational Layers of Quantum Chemistry
**Version:** 1.0
**Last Updated:** 2026-03-19
**Document Type:** Assessment Portfolio Document
**Target Publication:** J. Chem. Educ. Technology Report

---

## Executive Overview

This document consolidates all assessment instruments for Lab Pack #3 "Computational Layers of Quantum Chemistry" of the Interactive Quantum Chemistry Playground (IQCP). The assessment portfolio measures student learning across six learning outcomes (LO13--LO18) spanning three Phase 3 modules: basis function exploration (Module A), integral inspection (Module B), and electron density visualization (Module E density tab).

### Assessment Philosophy

The assessment framework follows the principles established in Lab Packs #1 and #2, extended to the computational infrastructure underlying SCF calculations:

1. **Constructive alignment:** All assessment items trace directly to stated learning outcomes LO13--LO18 (Wiggins & McTighe, 2005). Item stems reference authentic computational scenarios that students encounter in the IQCP worksheet activities.
2. **Multiple measures:** Conceptual understanding (concept checks), procedural and interpretive skills (worksheet), and integrative reasoning (performance tasks) are assessed through complementary instruments.
3. **Formative and summative:** Pre/post concept checks enable learning gain measurement; the worksheet provides formative guidance; performance tasks provide summative evaluation of applied competency.
4. **Misconception targeting:** Distractors in multiple-choice items are drawn from documented student misconceptions in the chemistry education literature and common errors observed in computational chemistry instruction, ensuring diagnostic value beyond simple right/wrong scoring.

### Total Point Allocation

| Instrument | Points | Percentage | Purpose |
|------------|--------|------------|---------|
| Concept Check (Pre) | 6 | 3.7% | Baseline knowledge assessment |
| Concept Check (Post) | 8 | 4.9% | Learning gain measurement |
| Worksheet | 106 | 65.0% | Guided POE exploration assessment |
| Performance Tasks | 43 | 26.4% | Applied competency assessment |
| **Total Portfolio** | **163** | **100%** | |

### Learning Outcomes Assessed

| ID | Learning Outcome | Bloom's Level | Primary Instruments |
|----|------------------|---------------|---------------------|
| **LO13** | Describe components of a contracted Gaussian basis function | Understand | CC (P1, Q1), WS (Q1.1-Q1.4) |
| **LO14** | Compare basis sets by examining radial profiles | Analyze | CC (P2, Q2, Q3), WS (Q1.5-Q1.8), PT-Basis |
| **LO15** | Identify physical meaning of S, T, V integrals | Analyze | CC (P3, Q4), WS (Q2.1-Q2.4) |
| **LO16** | Trace Fock matrix construction from H^core and G(P) | Analyze | CC (P4, Q5, Q6), WS (Q2.5-Q2.8), PT-Integral |
| **LO17** | Interpret density isosurfaces and cross-sections | Analyze | CC (P5, Q7), WS (Q3.1-Q3.4) |
| **LO18** | Analyze difference density maps for bonding | Analyze | CC (P6, Q8), WS (Q3.5-Q3.8) |

---

## Pre-Activity Concept Check

**Administration:** Before students access IQCP. No IQCP use permitted during the check.
**Time:** 5--7 minutes
**Total Points:** 6 (1 point per item)
**Format:** 5 MC items + 1 SA item

The pre-check establishes baseline understanding of basis functions, integrals, density, and Fock matrix concepts that students bring from prior coursework (Lab Packs #1 and #2, physical chemistry). Items target lower Bloom's levels (Remember, Understand) appropriate for a pre-assessment.

---

### Item P1: What Is a Contracted Gaussian?

**Primary LO:** LO13 | **Secondary LO:** -- | **Bloom's Level:** Remember | **Type:** MC

In computational chemistry, a "contracted Gaussian basis function" is:

**(a)** A single Gaussian function with an exponent that contracts (decreases) during the SCF calculation.

**(b)** A fixed linear combination of multiple Gaussian primitive functions, each with a different exponent and coefficient.

**(c)** A Gaussian function that has been compressed to fit inside the molecular boundary surface.

**(d)** A Gaussian function that describes only the core electrons of an atom, not the valence electrons.

---

### Item P2: Basis Set Size and Energy

**Primary LO:** LO14 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

A student computes the SCF energy for H2O using STO-3G (7 basis functions) and 6-31G (13 basis functions). Which statement is correct?

**(a)** STO-3G gives a lower energy because it uses fewer basis functions and therefore has less computational noise.

**(b)** Both basis sets give the same energy because the molecule is the same.

**(c)** 6-31G gives a lower energy because more basis functions provide greater variational freedom to minimize the energy.

**(d)** 6-31G gives a higher energy because more basis functions introduce more electron repulsion terms.

---

### Item P3: Overlap Integral Meaning

**Primary LO:** LO15 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

The overlap integral S_ij between two basis functions phi_i and phi_j is defined as S_ij = integral phi_i(r) phi_j(r) dr. What does this integral measure?

**(a)** The energy of the interaction between the two basis functions.

**(b)** The spatial overlap between the two functions -- how much they share the same region of space.

**(c)** The probability that an electron is located exactly between the two atoms.

**(d)** The chemical bond strength between the atoms carrying the two basis functions.

---

### Item P4: Fock Matrix Components

**Primary LO:** LO16 | **Secondary LO:** -- | **Bloom's Level:** Remember | **Type:** MC

The Fock matrix in Hartree-Fock theory is given by F = H^core + G(P). Which of the following correctly identifies what H^core and G(P) represent?

**(a)** H^core = electron-electron repulsion; G(P) = kinetic energy + nuclear attraction.

**(b)** H^core = kinetic energy + nuclear attraction (one-electron terms); G(P) = electron-electron interactions (Coulomb and exchange, weighted by the density matrix).

**(c)** H^core = the Hamiltonian for all electrons; G(P) = a correction for electron correlation.

**(d)** H^core = overlap matrix; G(P) = kinetic energy matrix.

---

### Item P5: What Does Electron Density Represent?

**Primary LO:** LO17 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

The electron density rho(r) at a point r in a molecule represents:

**(a)** The number of electrons located exactly at point r.

**(b)** The probability density for finding an electron at point r -- a quantity proportional to the likelihood of observing an electron in a small volume around r.

**(c)** The electric charge stored at point r, measured in coulombs per cubic meter.

**(d)** The energy of an electron at point r, which determines whether the electron is bound or free.

---

### Item P6: What Is a Difference Density?

**Primary LO:** LO18 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** SA

When forming a molecule from atoms, the electron density changes. The "difference density" (or deformation density) is computed by subtracting a reference density from the molecular density. In 1--2 sentences, describe what the difference density shows and what the reference (promolecule) density represents.

*Expected response length: 1--2 sentences*

---

## Post-Activity Concept Check

**Administration:** Immediately after worksheet completion. No IQCP access during the check.
**Time:** 7--10 minutes
**Total Points:** 8 (1 point per item)
**Format:** 5 MC items + 3 SA items

The post-check assesses understanding at higher cognitive levels (Apply, Analyze) than the pre-check, reflecting the learning that should have occurred during the IQCP activities. Items target the same LOs as the pre-check but require deeper reasoning.

---

### Item Q1: Explaining Why Contraction Uses Multiple Exponents

**Primary LO:** LO13 | **Secondary LO:** -- | **Bloom's Level:** Apply | **Type:** SA

In the IQCP Basis Explorer, you saw that the STO-3G hydrogen 1s basis function uses three Gaussian primitives with exponents 3.43, 0.62, and 0.17. In 2--3 sentences, explain why three different exponents are needed rather than a single Gaussian. What aspect of the atomic orbital does a single Gaussian fail to reproduce?

*Expected response length: 2--3 sentences*

---

### Item Q2: Predicting Basis Set Behavior for a Novel Atom

**Primary LO:** LO14 | **Secondary LO:** -- | **Bloom's Level:** Apply | **Type:** MC

You learned that 6-31G splits the valence shell into two contractions (inner + outer). If you compare the STO-3G and 6-31G radial profiles for the nitrogen 2p shell, you predict that:

**(a)** The 6-31G profile will be identical to STO-3G because nitrogen and oxygen have the same number of p orbitals.

**(b)** The 6-31G profile will extend further from the nucleus than STO-3G because the outer diffuse function adds coverage at larger distances.

**(c)** The 6-31G profile will be more compact than STO-3G because 6-31G uses more primitives, concentrating the function near the nucleus.

**(d)** The 6-31G profile will show three separate peaks corresponding to the three 2p orbitals (2px, 2py, 2pz).

---

### Item Q3: Evaluating a Basis Set Choice

**Primary LO:** LO14 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** SA

A researcher needs to compute the binding energy of a protein-ligand complex (over 1000 atoms). They consider using cc-pVQZ (a very large basis set). A colleague suggests using 6-31G instead. In 2--3 sentences, advise the researcher. Your response should reference computational cost scaling and the concept of diminishing returns.

*Expected response length: 2--3 sentences*

---

### Item Q4: Predicting How an Integral Changes with Geometry

**Primary LO:** LO15 | **Secondary LO:** -- | **Bloom's Level:** Apply | **Type:** MC

In the IQCP Integral Inspector, you examined overlap integrals for H2 at R = 1.4 bohr and found S(1,2) = 0.66. If the bond length is decreased to R = 0.7 bohr (half the original distance), what do you predict will happen to S(1,2)?

**(a)** S(1,2) will decrease because the basis functions are now too close and start to cancel each other.

**(b)** S(1,2) will increase toward 1.0 because the basis functions are centered closer together and overlap more.

**(c)** S(1,2) will remain 0.66 because the overlap integral depends on the basis function type, not on the distance.

**(d)** S(1,2) will become negative because at very short distances the basis functions interfere destructively.

---

### Item Q5: Explaining Why the Fock Matrix Is Rebuilt

**Primary LO:** LO16 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** SA

In the IQCP Fock build trace, you saw that the Fock matrix F = H^core + G(P) is rebuilt at every SCF iteration. In 2--3 sentences, explain why the Fock matrix changes from one iteration to the next. Which component (H^core or G(P)) changes, and why?

*Expected response length: 2--3 sentences*

---

### Item Q6: Identifying Missing Components in the Fock Matrix

**Primary LO:** LO16 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** MC

A student writes the Fock matrix element as: F(mu,nu) = H^core(mu,nu) + sum_{lambda,sigma} (mu nu | lambda sigma). What is missing from this formula?

**(a)** Nothing is missing -- this is the correct Fock matrix formula.

**(b)** The density matrix P(lambda,sigma) as a weighting factor, and the exchange integral contribution -0.5 sum P(lambda,sigma) (mu lambda | nu sigma).

**(c)** The overlap matrix S(mu,nu) must be subtracted from the Fock matrix.

**(d)** The nuclear repulsion energy E_nuc must be added to each Fock matrix element.

---

### Item Q7: Interpreting a Density Cross-Section

**Primary LO:** LO17 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** MC

A student examines a 2D electron density cross-section for H2O in the molecular plane. They observe sharp peaks at the nuclear positions and elevated density between each O and H. However, the oxygen peak is much taller than the hydrogen peaks. The best explanation for this difference is:

**(a)** The oxygen atom has more electrons (8 vs. 1 per hydrogen) and the core electrons create a very high density concentration near the oxygen nucleus.

**(b)** The hydrogen atoms are farther from the center of the cross-section, so perspective distortion makes them appear smaller.

**(c)** The basis functions on hydrogen are less accurate than those on oxygen, producing a lower density.

**(d)** The SCF calculation converged incorrectly, putting too many electrons on oxygen.

---

### Item Q8: Interpreting Difference Density for a Polar Bond

**Primary LO:** LO18 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** MC

When examining the difference density (Delta-rho) for H2O, you observe that the accumulation region (Delta-rho > 0) along each O-H bond is shifted toward the oxygen atom rather than being centered between O and H. The best interpretation of this observation is:

**(a)** The calculation has a systematic error that places too much density on oxygen.

**(b)** Oxygen is more electronegative than hydrogen, so it attracts more electron density toward itself during bond formation.

**(c)** The STO-3G basis set has more functions on oxygen, artificially shifting the density toward it.

**(d)** The promolecule reference density was computed incorrectly, causing an artificial shift.

---

## Alignment Matrix

### Item-to-LO Mapping

| Item | Content Focus | Primary LO | Secondary LO | Bloom's Level | Item Type | Points | Targeted Misconception |
|------|---------------|------------|--------------|---------------|-----------|--------|----------------------|
| P1 | Contracted Gaussian definition | LO13 | -- | Remember | MC | 1 | "Contraction = compression" |
| P2 | Basis set size and energy | LO14 | -- | Understand | MC | 1 | "Fewer functions = less noise" |
| P3 | Overlap integral meaning | LO15 | -- | Understand | MC | 1 | "Overlap = bond strength" |
| P4 | Fock matrix components | LO16 | -- | Remember | MC | 1 | "H^core = all-electron Hamiltonian" |
| P5 | Electron density meaning | LO17 | -- | Understand | MC | 1 | "Density = electron count at a point" |
| P6 | Difference density concept | LO18 | -- | Understand | SA | 1 | -- |
| Q1 | Why multiple exponents needed | LO13 | -- | Apply | SA | 1 | -- |
| Q2 | Split-valence for novel atom | LO14 | -- | Apply | MC | 1 | "More primitives = more compact" |
| Q3 | Basis set choice evaluation | LO14 | -- | Analyze | SA | 1 | "Always use biggest basis" |
| Q4 | Overlap change with geometry | LO15 | -- | Apply | MC | 1 | "Overlap is geometry-independent" |
| Q5 | Why Fock matrix is rebuilt | LO16 | -- | Analyze | SA | 1 | -- |
| Q6 | Missing Fock matrix components | LO16 | -- | Analyze | MC | 1 | "F = H^core + ERIs (no P, no exchange)" |
| Q7 | Density cross-section peaks | LO17 | -- | Analyze | MC | 1 | "Peak difference = calculation error" |
| Q8 | Polar bond difference density | LO18 | -- | Analyze | MC | 1 | "Asymmetry = calculation error" |

### LO Coverage Verification

| LO | Pre-Check Items | Post-Check Items | Total Items | Coverage Status |
|----|-----------------|------------------|-------------|-----------------|
| LO13 | P1 | Q1 | 2 | Adequate |
| LO14 | P2 | Q2, Q3 | 3 | Strong |
| LO15 | P3 | Q4 | 2 | Adequate |
| LO16 | P4 | Q5, Q6 | 3 | Strong |
| LO17 | P5 | Q7 | 2 | Adequate |
| LO18 | P6 | Q8 | 2 | Adequate |
| **Total** | **6** | **8** | **14** | **All LOs >= 2 items** |

---

## Cognitive Level Distribution Analysis

### Distribution Across Concept Check Items (Pre + Post)

| Bloom's Level | Items | Count | Percentage |
|---------------|-------|-------|------------|
| Remember | P1, P4 | 2 | 14.3% |
| Understand | P2, P3, P5, P6 | 4 | 28.6% |
| Apply | Q1, Q2, Q4 | 3 | 21.4% |
| Analyze | Q3, Q5, Q6, Q7, Q8 | 5 | 35.7% |

**Aggregate summary (using 2-level grouping):**

| Grouping | Target % | Actual % | Items |
|----------|----------|----------|-------|
| Remember/Understand | ~25% | 42.9% (6/14) | P1, P2, P3, P4, P5, P6 |
| Apply/Analyze | ~75% | 57.1% (8/14) | Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8 |

**Design note:** The aggregate distribution across concept checks alone shows more Remember/Understand items than the 25% target because the 6-item pre-check intentionally targets lower cognitive levels (baseline assessment). When the full portfolio is considered (concept checks + worksheet + performance tasks), the distribution shifts substantially:

### Distribution Across Full Portfolio

| Grouping | Concept Checks (14 items) | Worksheet (26 items) | Performance Tasks (4 parts) | Combined (44 items) | Combined % |
|----------|--------------------------|---------------------|---------------------------|---------------------|------------|
| Remember/Understand | 6 | 1 (Q2.5) | 0 | 7 | 15.9% |
| Apply/Analyze | 8 | 23 | 4 | 35 | 79.5% |
| Evaluate | 0 | 2 (Q1.8, Q4.2) | 0 | 2 | 4.5% |

The full portfolio achieves approximately **16% Remember/Understand and 84% Apply/Analyze/Evaluate**, exceeding the 75% target for higher-order thinking. The elevation in Apply/Analyze reflects the POE framework's emphasis on prediction, observation, and explanation -- activities that inherently require application and analysis.

---

## Performance Task Rubrics

The following performance tasks may be administered as part of the worksheet debrief, as a separate assessment session, or as take-home assignments. Each task is scored on a 4-point analytic rubric.

### Performance Task 1: Basis Set Analysis (PT-Basis)

**Target LOs:** LO13 (primary), LO14 (secondary)
**Total Points:** 21 (3 dimensions x 7-point scale)

#### Task Description

A student is given the following data for the carbon 2p basis function in two basis sets:

**STO-3G Carbon 2p (1 contraction, 3 primitives):**

| Primitive | Exponent | Coefficient |
|-----------|----------|-------------|
| 1 | 2.9412 | 0.1559 |
| 2 | 0.6835 | 0.6077 |
| 3 | 0.2222 | 0.3920 |

**6-31G Carbon 2p (2 contractions: 3 primitives + 1 primitive):**

Inner contraction:

| Primitive | Exponent | Coefficient |
|-----------|----------|-------------|
| 1 | 9.4398 | 0.0381 |
| 2 | 2.0024 | 0.2095 |
| 3 | 0.5460 | 0.5085 |

Outer function:

| Primitive | Exponent | Coefficient |
|-----------|----------|-------------|
| 1 | 0.1517 | 1.0000 |

**Part A (7 points):** Explain the structure of each basis set. How many independent basis functions does each provide for the 2p shell? Why does 6-31G have an additional function?

**Part B (7 points):** Predict the qualitative difference in the radial profiles. Which basis set will extend further from the nucleus? Identify the specific primitive responsible.

**Part C (7 points):** A researcher computes the SCF energy for methane (CH4) with both basis sets and finds STO-3G gives -39.727 Ha and 6-31G gives -40.195 Ha. Explain why 6-31G gives a lower energy. Would continuing to add more basis functions eventually give the exact energy? Why or why not?

#### PT-Basis Rubric

**Part A -- Basis Set Structure (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly identifies STO-3G as 1 contraction = 1 basis function per p component; 6-31G as 2 contractions = 2 basis functions per p component (inner 3-primitive + outer 1-primitive). Explains that the additional function provides variational freedom -- the SCF can independently weight the tight and diffuse components to optimize the description of bonding. |
| 3 (Proficient) | Correctly identifies contraction counts and basis function counts. Notes the extra function provides flexibility but explanation is partial. |
| 2 (Developing) | Identifies that 6-31G has more functions but cannot clearly explain the contraction structure or why the extra function matters. |
| 1 (Beginning) | Cannot identify the structure of either basis set or confuses contractions with primitives. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part B -- Radial Profile Prediction (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly predicts 6-31G extends further. Identifies the outer function (exponent = 0.1517) as the specific primitive responsible because it has the smallest exponent and therefore the widest spatial extent. May note that this exponent is smaller than any STO-3G primitive (0.2222 is the smallest), confirming the prediction. |
| 3 (Proficient) | Correctly predicts 6-31G extends further and mentions the diffuse function. Does not identify the specific exponent or compare quantitatively. |
| 2 (Developing) | Predicts 6-31G extends further but cannot identify the responsible primitive or gives incorrect reasoning. |
| 1 (Beginning) | Incorrect prediction or no prediction offered. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part C -- Energy and Variational Principle (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly explains that 6-31G gives lower energy due to the variational principle -- more basis functions span a larger subspace of Hilbert space. Correctly states that adding more functions would continue to lower the energy toward the complete basis set (CBS) limit, but this limit is NOT the exact energy because RHF still misses electron correlation. Distinguishes basis set error from method error. |
| 3 (Proficient) | Invokes the variational principle correctly. Notes that more functions help but may not clearly distinguish CBS limit from exact energy, or does not mention correlation. |
| 2 (Developing) | States that more functions give lower energy but cannot explain why (does not invoke variational principle) or incorrectly states that infinite functions give the exact answer. |
| 1 (Beginning) | Cannot explain the energy difference or states that STO-3G should give a lower energy. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

---

### Performance Task 2: Integral Interpretation (PT-Integral)

**Target LOs:** LO15 (primary), LO16 (secondary)
**Total Points:** 22 (3 dimensions, 7-7-8 point scale)

#### Task Description

A student examines the following one-electron integral data for LiH (STO-3G, R = 3.0 bohr), which has 6 basis functions (Li: 1s, 2s; H: 1s, plus the 2p shells on Li are omitted for simplicity in this task):

**Simplified 3x3 one-electron matrices (Li 1s, Li 2s, H 1s):**

| | Li 1s | Li 2s | H 1s |
|---|-------|-------|------|
| **S** | 1.000 | 0.237 | 0.009 |
| | 0.237 | 1.000 | 0.507 |
| | 0.009 | 0.507 | 1.000 |

| | Li 1s | Li 2s | H 1s |
|---|-------|-------|------|
| **T** | 1.548 | -0.164 | -0.001 |
| | -0.164 | 0.155 | 0.087 |
| | -0.001 | 0.087 | 0.760 |

| | Li 1s | Li 2s | H 1s |
|---|-------|-------|------|
| **V** | -4.792 | -0.826 | -0.339 |
| | -0.826 | -1.449 | -0.805 |
| | -0.339 | -0.805 | -1.252 |

**Part A (7 points):** The overlap between Li 1s and H 1s is 0.009 -- nearly zero. Explain why this overlap is so small despite the atoms being bonded at 3.0 bohr. Contrast this with the Li 2s - H 1s overlap of 0.507.

**Part B (7 points):** Examine the kinetic energy matrix T. Why is T(Li 1s, Li 1s) = 1.548 much larger than T(Li 2s, Li 2s) = 0.155? Relate your answer to the spatial extent of these basis functions.

**Part C (8 points):** Using the V matrix, compute H^core(Li 2s, H 1s) = T(Li 2s, H 1s) + V(Li 2s, H 1s). The Fock matrix element F(Li 2s, H 1s) at convergence is approximately -0.30. Given that F = H^core + G(P), what is the approximate value of G(Li 2s, H 1s)? Is G positive or negative? What does the sign of G tell you about whether electron-electron interactions strengthen or weaken the Li-H bond in this matrix element?

#### PT-Integral Rubric

**Part A -- Overlap Analysis (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly explains that Li 1s is a very tight core orbital (large exponents, small spatial extent), so it has negligible overlap with H 1s despite the atoms being bonded. Contrasts with Li 2s, which is a diffuse valence orbital extending much further from the Li nucleus, producing significant overlap with H 1s at 3.0 bohr. Notes that overlap depends on spatial extent, not just distance. |
| 3 (Proficient) | Identifies Li 1s as compact/core and Li 2s as diffuse/valence. Correct contrast but explanation of why spatial extent matters is partial. |
| 2 (Developing) | Notes the overlap difference but attributes it only to distance or "different types of orbitals" without connecting to spatial extent. |
| 1 (Beginning) | Cannot explain the overlap difference or attributes it to a calculation error. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part B -- Kinetic Energy Analysis (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly explains that the kinetic energy operator involves the second derivative (Laplacian). Tight (compact) functions like Li 1s have sharper curvature, producing larger kinetic energy. Diffuse functions like Li 2s have gentler curvature, producing smaller kinetic energy. May connect to the uncertainty principle: more spatially confined electrons have higher kinetic energy. |
| 3 (Proficient) | Identifies that tighter functions have larger kinetic energy and connects to curvature or confinement. Does not invoke the Laplacian or uncertainty principle explicitly. |
| 2 (Developing) | Notes T(Li 1s) > T(Li 2s) but cannot clearly explain why. May incorrectly attribute it to "more energy" without connecting to spatial extent. |
| 1 (Beginning) | Cannot explain the kinetic energy difference or confuses T with V or S. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part C -- Fock Matrix Decomposition (8 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplary) | Correctly computes H^core(Li 2s, H 1s) = 0.087 + (-0.805) = -0.718. Correctly computes G = F - H^core = -0.30 - (-0.718) = +0.418. Identifies G as positive. Correctly interprets: G > 0 means electron-electron interactions make this off-diagonal element less negative, weakening the effective one-electron attraction that drives bonding. In other words, electron repulsion partially offsets the nuclear attraction that creates the bond. |
| 3 (Proficient) | Correct computation of H^core and G. Identifies G as positive. Interpretation is partial -- notes that G opposes H^core but does not clearly connect to bond strengthening/weakening. |
| 2 (Developing) | Computation of H^core is correct or nearly correct. G computation has an error or is absent. Sign interpretation missing. |
| 1 (Beginning) | Cannot perform the computation or confuses the matrix elements. No interpretation. |

**Scoring conversion:** 4 -> 8, 3 -> 6, 2 -> 3, 1 -> 1

---

## Concept Check Answer Key

### Pre-Check Answers

#### P1: What Is a Contracted Gaussian?

**Correct answer: (b)** A fixed linear combination of multiple Gaussian primitive functions, each with a different exponent and coefficient.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses "contracted" with "contracting during calculation"; basis set is fixed, not optimized during SCF | **"Contraction = dynamic compression"** |
| **(b)** | **Correct.** Identifies the static, predetermined nature of the combination | -- |
| (c) | Confuses contraction with spatial compression; relates to orbital boundary misconception | **"Compressed to fit inside boundary"** |
| (d) | Confuses contracted basis functions with core orbitals; both core and valence can be contracted | Core/valence confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P2: Basis Set Size and Energy

**Correct answer: (c)** 6-31G gives a lower energy because more basis functions provide greater variational freedom.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses "fewer" with "cleaner" -- misunderstands variational principle | **"Fewer functions = less noise"** |
| (b) | Fails to recognize that basis set quality affects computed energy | "Same molecule = same energy" |
| **(c)** | **Correct.** Applies variational principle | -- |
| (d) | Confuses number of basis functions with number of physical interactions | "More functions = more repulsion" |

**Scoring:** 1 point for (c), 0 for all others.

---

#### P3: Overlap Integral Meaning

**Correct answer: (b)** The spatial overlap between the two functions -- how much they share the same region of space.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses the overlap integral (a geometric measure) with an energy integral | Overlap = energy confusion |
| **(b)** | **Correct.** Identifies the spatial/geometric nature of the overlap | -- |
| (c) | Confuses overlap with probability; the overlap is not a probability at a point | **"Overlap = probability between atoms"** |
| (d) | Equates overlap with bond strength; overlap is a necessary but not sufficient condition for bonding | **"Overlap = bond strength"** |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P4: Fock Matrix Components

**Correct answer: (b)** H^core = kinetic energy + nuclear attraction (one-electron terms); G(P) = electron-electron interactions weighted by the density matrix.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Reverses the roles of H^core and G(P) | Component role reversal |
| **(b)** | **Correct.** Identifies one-electron vs. two-electron decomposition | -- |
| (c) | Confuses H^core with the full Hamiltonian; G(P) is not a correlation correction (HF has no explicit correlation) | **"H^core = full Hamiltonian"** |
| (d) | Confuses matrix types entirely; S is not part of F, and T alone is not G | Matrix identity confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P5: What Does Electron Density Represent?

**Correct answer: (b)** The probability density for finding an electron at point r.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses probability density (continuous function) with a discrete count; rho(r) can exceed 1 | **"Density = electron count at a point"** |
| **(b)** | **Correct.** Identifies the probabilistic interpretation | -- |
| (c) | Confuses electron density (quantum chemistry) with charge density (electrostatics); units are different | Unit/concept confusion |
| (d) | Confuses density with potential energy; rho is a spatial distribution, not an energy | Density = energy confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P6: What Is a Difference Density?

**Correct answer (exemplar):** The difference density shows how the electron distribution changes when atoms form a molecule. It is computed by subtracting the promolecule density -- the sum of non-interacting atomic densities placed at the molecular positions -- from the actual molecular density. Positive regions indicate electron accumulation due to bonding; negative regions indicate depletion.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies the difference density as showing changes in density upon molecule formation AND describes the promolecule as non-interacting atoms at molecular positions. Partial credit for identifying it as "molecular minus atomic" without specifying that the reference atoms are at molecular positions. |
| 0 | Response is absent, incoherent, or describes the difference density as "total density" or "the density of the bond." |

---

### Post-Check Answers

#### Q1: Why Multiple Exponents?

**Correct answer (exemplar):** Three exponents are needed because a single Gaussian cannot reproduce the correct shape of an atomic orbital. A single Gaussian has a smooth, rounded peak at the nucleus (zero slope), but a real atomic orbital has a cusp (discontinuous derivative) at the nucleus. Additionally, a single Gaussian decays as exp(-alpha r^2), which falls off too rapidly at large distances compared to the correct exponential decay exp(-zeta r). Multiple Gaussians with different exponents -- tight, medium, and diffuse -- combine to approximate both the cusp behavior and the correct long-range decay.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies at least one failure of a single Gaussian (cusp, tail, or overall shape mismatch) AND explains that multiple exponents with different widths address different spatial regions. |
| 0 | Response does not identify any specific limitation of a single Gaussian, or states that multiple Gaussians are used "for accuracy" without specifying what aspect is improved. |

---

#### Q2: Split-Valence for Novel Atom

**Correct answer: (b)** The 6-31G profile will extend further from the nucleus because the outer diffuse function adds coverage at larger distances.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Assumes same number of p orbitals means same profile; confuses orbital count with basis function shape | "Same orbital type = same profile" |
| **(b)** | **Correct.** Applies split-valence understanding to a novel element | -- |
| (c) | Reverses the relationship; more primitives in 6-31G include a diffuse function, not just tight ones | **"More primitives = more compact"** |
| (d) | Confuses 2px/2py/2pz components with radial profile features; the radial profile is the same for all three p components | Angular/radial confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### Q3: Basis Set Choice Evaluation

**Correct answer (exemplar):** For a system with over 1000 atoms, cc-pVQZ would be computationally prohibitive. The number of two-electron integrals scales as N^4, and cc-pVQZ has approximately 5x as many basis functions per atom as 6-31G, making the calculation roughly 625 times more expensive per integral batch. I would recommend 6-31G (or possibly 6-31G*) as a pragmatic choice: it provides the essential split-valence flexibility at manageable cost. The energy improvement from 6-31G to cc-pVQZ follows diminishing returns -- most of the basis set error is captured by the first step beyond minimal basis.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response recommends a basis set with a justified cost-accuracy argument. Must reference N^4 scaling or computational cost as a limiting factor. Must mention diminishing returns or basis set convergence. Accept any reasonable recommendation (6-31G, 6-31G*, or even a mixed strategy) as long as the reasoning is sound. |
| 0 | Recommends cc-pVQZ without mentioning cost, states "always use the biggest," or provides no justification. |

---

#### Q4: Overlap Change with Geometry

**Correct answer: (b)** S(1,2) will increase toward 1.0 because the basis functions are centered closer together and overlap more.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Incorrect; s-type functions on the same atom do not cancel | "Too close = cancellation" |
| **(b)** | **Correct.** Closer centers produce more spatial overlap | -- |
| (c) | Assumes overlap is geometry-independent; it strongly depends on distance | **"Overlap depends only on orbital type"** |
| (d) | s-type overlap between same-type functions is always non-negative | Sign confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### Q5: Why the Fock Matrix Is Rebuilt

**Correct answer (exemplar):** The Fock matrix F = H^core + G(P) changes from iteration to iteration because the G(P) component depends on the density matrix P, which is updated at each iteration based on the new MO coefficients. H^core (kinetic energy + nuclear attraction) does not change because it depends only on the basis functions and nuclear positions, which are fixed. Since G(P) encodes the electron-electron repulsion weighted by the current electron distribution, it must be recomputed whenever P changes to maintain self-consistency.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies G(P) as the component that changes AND correctly states that it depends on the density matrix P, which updates each iteration. Must note that H^core is constant. |
| 0 | Does not identify which component changes, or states that both components change, or does not mention the density matrix. |

---

#### Q6: Missing Fock Matrix Components

**Correct answer: (b)** The density matrix P(lambda,sigma) as a weighting factor, and the exchange integral contribution.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Accepts the incomplete formula as correct | **"F = H^core + ERIs (no P, no exchange)"** |
| **(b)** | **Correct.** Identifies both missing elements | -- |
| (c) | Confuses the Roothaan-Hall equation (FC = SCe) with the Fock matrix definition | S confusion |
| (d) | Confuses a scalar (E_nuc) with a matrix element; nuclear repulsion does not enter the Fock matrix | Scalar/matrix confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### Q7: Density Cross-Section Peaks

**Correct answer: (a)** Oxygen has more electrons and core electrons create a high concentration near the oxygen nucleus.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| **(a)** | **Correct.** Identifies electron count and core concentration | -- |
| (b) | Attributes a physical feature to a visualization artifact | **"Peak difference = perspective distortion"** |
| (c) | Blames basis set quality; basis functions are comparable in quality per atom | "Bad basis = wrong density" |
| (d) | Attributes a physically correct feature to a calculation error | **"Asymmetric density = calculation error"** |

**Scoring:** 1 point for (a), 0 for all others.

---

#### Q8: Polar Bond Difference Density

**Correct answer: (b)** Oxygen is more electronegative, attracting more density toward itself during bond formation.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Attributes a physically correct feature to a systematic error | **"Asymmetry = calculation error"** |
| **(b)** | **Correct.** Connects electronegativity to density shift | -- |
| (c) | Blames basis set asymmetry; this is a method-level effect, not a basis set artifact | "More functions on O = artificial shift" |
| (d) | Blames the reference calculation; the promolecule is computed correctly | "Wrong reference = wrong difference" |

**Scoring:** 1 point for (b), 0 for all others.

---

## Point Allocation

### Summary by Instrument

| Instrument | Points | Percentage |
|------------|--------|------------|
| Pre-Activity Concept Check (P1--P6) | 6 | 3.7% |
| Post-Activity Concept Check (Q1--Q8) | 8 | 4.9% |
| Student Worksheet (Q1.1--Q4.2) | 106 | 65.0% |
| Performance Task: Basis Set Analysis (PT-Basis) | 21 | 12.9% |
| Performance Task: Integral Interpretation (PT-Integral) | 22 | 13.5% |
| **Total** | **163** | **100%** |

### Worksheet Points by Section

| Section | Questions | Points | Primary LOs |
|---------|-----------|--------|-------------|
| Section 1: Basis Function Exploration | Q1.1--Q1.8 | 30 | LO13, LO14 |
| Section 2: Integral Inspection & Fock Tracing | Q2.1--Q2.8 | 34 | LO15, LO16 |
| Section 3: Electron Density & Difference Density | Q3.1--Q3.8 | 32 | LO17, LO18 |
| Section 4: Synthesis | Q4.1--Q4.2 | 10 | Integrative |
| **Worksheet Total** | **26 items** | **106** | |

### Performance Task Points by Dimension

| Task | Part | Points | Primary LO |
|------|------|--------|------------|
| PT-Basis | A: Basis set structure | 7 | LO13 |
| PT-Basis | B: Radial profile prediction | 7 | LO14 |
| PT-Basis | C: Energy and variational principle | 7 | LO14 |
| PT-Integral | A: Overlap analysis | 7 | LO15 |
| PT-Integral | B: Kinetic energy analysis | 7 | LO15 |
| PT-Integral | C: Fock matrix decomposition | 8 | LO16 |
| **Performance Total** | | **43** | |

### Recommended Grade Weighting

For courses using Lab Pack #3 as a graded assignment:

```
Lab Grade = (0.65 x Worksheet%) + (0.09 x ConceptCheck%) + (0.26 x PerformanceTask%)
```

This weighting reflects the centrality of the guided exploration (worksheet) while ensuring that higher-order reasoning (performance tasks) contributes meaningfully to the grade.

---

## Sample Responses

### SA Item Q1: Exemplar, Adequate, and Inadequate

**Exemplar (1 point):**
"A single Gaussian has a smooth, rounded peak at r = 0 with zero slope, but a real atomic orbital has a cusp there. Also, a single Gaussian decays as exp(-alpha r^2), which drops off too fast at large r. Using three Gaussians with different exponents (tight, medium, diffuse) approximates both the near-nucleus cusp and the correct long-range tail."

**Adequate (1 point, borderline):**
"A single Gaussian is too wide or too narrow -- it can't get both the peak and the tail right. Three Gaussians with different widths can approximate the Slater-type shape better."

**Inadequate (0 points):**
"Three Gaussians are used because they are more accurate than one. More is always better."

---

### SA Item Q5: Exemplar, Adequate, and Inadequate

**Exemplar (1 point):**
"G(P) changes because it depends on the density matrix P, which is updated every iteration when new MO coefficients are computed. H^core stays the same because it only depends on basis functions and nuclear positions, which don't change during SCF."

**Adequate (1 point, borderline):**
"The part with the electron repulsion (G) has to be recalculated because the density changes. The H^core part doesn't change."

**Inadequate (0 points):**
"The Fock matrix is rebuilt because the calculation hasn't converged yet."

---

### Performance Task PT-Integral Part C: Exemplar, Adequate, and Inadequate

**Exemplar (8 points):**
"H^core(Li 2s, H 1s) = T(Li 2s, H 1s) + V(Li 2s, H 1s) = 0.087 + (-0.805) = -0.718. Since F = H^core + G(P), we get G = F - H^core = -0.30 - (-0.718) = +0.418. G is positive, meaning electron-electron interactions make this Fock matrix element less negative than H^core alone. Since a more negative off-diagonal Fock element promotes bonding (it contributes to lower orbital energies), the positive G partially counteracts the bonding effect of H^core. In other words, electron repulsion weakens the Li-H bonding interaction in this matrix element."

**Adequate (6 points):**
"H^core = 0.087 - 0.805 = -0.718. G = -0.30 - (-0.718) = +0.418. G is positive, which means it opposes the negative H^core. Electron repulsion makes bonding weaker."

**Inadequate (3 points):**
"H^core = T + V = some negative number. G must be the rest of F. I think G is positive."

---

## Validity Argument (Kane, 2006)

This section presents the validity argument for the Lab Pack #3 assessment portfolio, structured according to Kane's (2006) argument-based approach to validation. The argument addresses three inferential links: scoring, generalization, and extrapolation.

### 1. Scoring Inference

**Claim:** Observed scores accurately reflect the quality of student responses.

**Evidence and warrants:**

*Multiple-choice items:*
- MC items (P1--P5, Q2, Q4, Q6, Q7, Q8) are scored dichotomously (1 or 0). Each item has one unambiguously correct answer verified against established quantum chemistry knowledge and cross-checked with PySCF 2.11.0 reference calculations.
- Distractors are constructed from documented misconceptions and common student errors in computational chemistry instruction, ensuring that incorrect responses carry diagnostic meaning.

*Short-answer items:*
- SA items (P6, Q1, Q3, Q5) include explicit scoring rubrics with criteria for full credit (1 point) and no credit (0 points). Each rubric specifies essential content elements.
- Sample responses at exemplar, adequate, and inadequate levels are provided for scorer calibration.

*Performance tasks:*
- Both PT-Basis and PT-Integral use 4-point analytic rubrics with dimension-specific criteria. Each rubric level has concrete behavioral descriptors.
- Scoring conversion formulas (4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1) differentiate performance levels.

*Inter-rater reliability protocol:*
- For SA items and performance tasks, double-blind scoring of a random 20% sample with adjudication of disagreements exceeding 1 rubric level. Target inter-rater reliability: Cohen's kappa >= 0.70.

**Potential threats:**
- SA items scored as 0/1 may lack sensitivity for partial understanding. Mitigated by the worksheet items (multi-point rubrics) assessing the same LOs.
- Scorer drift over large batches. Mitigated by clear rubric anchors and calibration examples.

### 2. Generalization Inference

**Claim:** Scores generalize to the broader content domain defined by LO13--LO18.

**Evidence and warrants:**

*Content coverage:*
- The alignment matrix demonstrates every LO (LO13--LO18) is assessed by at least 2 concept check items (pre + post combined).
- Across the full portfolio (concept checks + worksheet + performance tasks), each LO is assessed by 4--10 items spanning multiple cognitive levels and formats.

*Cognitive level coverage:*
- Items span Remember through Analyze/Evaluate. The full portfolio distribution (16% R/U, 84% Ap/An/Ev) ensures generalization across cognitive levels.

*Item sampling:*
- The concept check contains 14 items sampling from 6 LOs. While individual LO sub-scale scores (2--3 items) are unreliable, the aggregate score provides a reliable measure of overall learning.

*Reliability targets:*
- Target Cronbach's alpha >= 0.70 for the combined concept check (14 items).

**Potential threats:**
- Content sampling may not fully represent each LO's breadth with only 2--3 items per LO. Mitigated by the 26-item worksheet providing dense coverage.
- Item difficulty may cluster. Post-pilot item analysis should examine difficulty indices.

### 3. Extrapolation Inference

**Claim:** Performance on these instruments indicates transferable understanding of computational chemistry infrastructure.

**Evidence and warrants:**

*Authentic tasks:*
- Concept check items reference scenarios that parallel authentic computational chemistry reasoning (choosing basis sets, interpreting matrix elements, analyzing density maps). These test conceptual understanding applicable to any quantum chemistry software.
- Performance tasks present novel data (carbon 2p, LiH) not encountered in the IQCP worksheet, requiring transfer.

*Misconception targeting:*
- Items target misconceptions documented in QC education: "basis function = atomic orbital," "more functions = always better," "overlap = bond strength," "Fock = H^core + ERIs," "density asymmetry = error." Reduction in misconception endorsement indicates conceptual change extending beyond IQCP.

*Transfer evidence (planned):*
- A future classroom pilot (N >= 30) should include transfer items presenting novel molecular systems and basis sets not used in the IQCP activities.

**Potential threats:**
- Students may learn to interpret IQCP-specific interfaces without developing transferable understanding. Mitigated by the POE framework and SA items requiring verbal explanation.
- The controlled lab setting may not predict open-ended performance.

### 4. Content Validity: SME Panel Review

**Status:** Planned but not yet completed.

**Protocol:** The assessment instruments should be reviewed by a panel of at least 3 subject matter experts (SMEs) with expertise in:
1. Quantum chemistry (computational methods, basis sets, integral evaluation)
2. Chemistry education research (assessment design, misconception identification)
3. Undergraduate instruction (content appropriateness, timing, difficulty)

**Review criteria:**
- Content accuracy: All correct answers are scientifically accurate
- Distractor quality: Distractors represent plausible misconceptions, not "obviously wrong" answers
- Cognitive level: Bloom's level assignments are appropriate
- Coverage: LO mapping is complete and balanced
- Clarity: Item stems and response options are unambiguous
- Accessibility: Language is appropriate for the target population (upper-division undergraduates)

**Review format:** Each SME independently rates each item on a 4-point scale (1 = not relevant, 2 = somewhat relevant, 3 = quite relevant, 4 = highly relevant) for content validity. The Content Validity Index (CVI) should be >= 0.80 across all items.

### Limitations of the Validity Argument

1. **No pilot data yet.** This validity argument is based on content analysis and expert review framework, not empirical data. Item statistics will be available after classroom administration.
2. **Single expert development.** Items were developed by one domain expert. The planned 3-expert panel review will strengthen content validity.
3. **Transfer evidence.** Direct evidence of extrapolation is planned but not yet available.
4. **Small number of items per LO.** With 2--3 concept check items per LO, individual LO sub-scale scores are unreliable. The portfolio-level gain target (normalized gain >= 0.3) aggregates across all items.

---

## Administration Guidelines

### Pre-Check Administration

1. Distribute the pre-check (P1--P6) at the start of the session, **before** students open IQCP.
2. Allow 5--7 minutes. Students should work individually without notes or discussion.
3. Collect all pre-check forms before distributing the worksheet or granting IQCP access.
4. **Important:** Do not review pre-check answers before the activity.

### Worksheet Administration

1. Distribute the worksheet and direct students to https://iqcp.dev.
2. Students work through Sections 1--4 at their own pace.
3. Target time: 60 minutes. Most students complete in 55--65 minutes.
4. Circulate to answer procedural questions (IQCP controls) but avoid giving conceptual answers.
5. Collect worksheets before distributing the post-check.

### Post-Check Administration

1. Distribute the post-check (Q1--Q8) immediately after worksheet collection.
2. Allow 7--10 minutes. Students should work individually without IQCP access, notes, or worksheets.
3. Collect all post-check forms.
4. **Optional:** Brief concept review after collection.

### Performance Task Administration

| Format | Timing | Recommended For |
|--------|--------|-----------------|
| Same session | After post-check (adds 20--25 min) | Extended lab periods (85+ min) |
| Separate session | Within 1 week of the lab | Standard lab periods |
| Take-home | Due within 1 week | Large enrollment courses |

### Data Collection for Publication

When collecting data for J. Chem. Educ. publication:

1. **IRB approval** required before data collection.
2. **Pre/post pairing:** Use anonymous identifiers to pair responses.
3. **Record administration conditions:** Date, time, section size, irregularities.
4. **Double-score** SA and PT items for 20% of the sample.

### Recommended Sample Sizes

| Analysis Type | Minimum N | Recommended N |
|---------------|-----------|---------------|
| Descriptive statistics | 15 | 30+ |
| Pre/post paired t-test | 20 | 40+ |
| Normalized gain analysis | 25 | 50+ |
| Item-level analysis | 30 | 75+ |
| Cronbach's alpha estimation | 30 | 100+ |

---

## Appendix: Summary Statistics Template

Use this template to record and analyze assessment data after classroom administration.

### A. Concept Check Statistics

#### Pre-Check Item Analysis

| Item | N | Mean | SD | Difficulty | Discrimination |
|------|---|------|----|----|----------------|
| P1 | | | | | |
| P2 | | | | | |
| P3 | | | | | |
| P4 | | | | | |
| P5 | | | | | |
| P6 | | | | | |
| **Total** | | | | | |

#### Post-Check Item Analysis

| Item | N | Mean | SD | Difficulty | Discrimination |
|------|---|------|----|----|----------------|
| Q1 | | | | | |
| Q2 | | | | | |
| Q3 | | | | | |
| Q4 | | | | | |
| Q5 | | | | | |
| Q6 | | | | | |
| Q7 | | | | | |
| Q8 | | | | | |
| **Total** | | | | | |

#### Learning Gain Analysis

| Metric | Value |
|--------|-------|
| Pre-check mean (out of 6) | |
| Post-check mean (out of 8) | |
| Normalized pre (out of 1.0) | |
| Normalized post (out of 1.0) | |
| Normalized gain g | |
| Effect size (Cohen's d) | |
| Paired t-test p-value | |

**Normalized Gain Formula (for unequal pre/post scales):**
```
g = (Post_normalized - Pre_normalized) / (1.0 - Pre_normalized)
where Pre_normalized = Pre_mean / 6, Post_normalized = Post_mean / 8
```

**Interpretation Guide:**
- g < 0.3: Low gain
- 0.3 <= g < 0.7: Medium gain
- g >= 0.7: High gain

### B. Worksheet Statistics

| Section | N | Mean | SD | Min | Max |
|---------|---|------|----|----|-----|
| Section 1: Basis Functions (30 pts) | | | | | |
| Section 2: Integrals & Fock (34 pts) | | | | | |
| Section 3: Density (32 pts) | | | | | |
| Section 4: Synthesis (10 pts) | | | | | |
| **Total (106 pts)** | | | | | |

### C. Performance Task Statistics

| Task | N | Mean | SD | Min | Max |
|------|---|------|----|----|-----|
| PT-Basis (21 pts) | | | | | |
| PT-Integral (22 pts) | | | | | |
| **Total (43 pts)** | | | | | |

### D. Reliability Analysis

| Instrument | Cronbach's alpha | KR-20 | Notes |
|------------|------------------|-------|-------|
| Pre-check (6 items) | | | Target >= 0.60 |
| Post-check (8 items) | | | Target >= 0.60 |
| Combined CC (14 items) | | | Target >= 0.70 |
| Worksheet (26 items) | | N/A | Target >= 0.80 |

### E. Inter-Rater Reliability

| Item | Cohen's kappa | % Exact Agreement | Notes |
|------|---------------|-------------------|-------|
| P6 | | | Target >= 0.70 |
| Q1 | | | Target >= 0.70 |
| Q3 | | | Target >= 0.70 |
| Q5 | | | Target >= 0.70 |
| PT-Basis (all parts) | | | Target >= 0.70 |
| PT-Integral (all parts) | | | Target >= 0.70 |

---

## References

Hehre, W.J., Stewart, R.F., and Pople, J.A. (1969). Self-consistent molecular-orbital methods. I. Use of Gaussian expansions of Slater-type atomic orbitals. *Journal of Chemical Physics*, 51(6), 2657--2664.

Kane, M.T. (2006). Validation. In R.L. Brennan (Ed.), *Educational Measurement* (4th ed., pp. 17--64). American Council on Education.

Nakhleh, M.B. (1992). Why some students don't learn chemistry: Chemical misconceptions. *Journal of Chemical Education*, 69(3), 191--196.

Tsaparlis, G., and Papaphotis, G. (2009). High-school students' conceptual difficulties and attempts at conceptual change: The case of basic quantum chemical concepts. *International Journal of Science Education*, 31(7), 895--930.

Wiggins, G., and McTighe, J. (2005). *Understanding by Design* (2nd ed.). Association for Supervision and Curriculum Development.

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-03-19 | IQCP Team | Initial release |

---

*Lab Pack #3 Assessment Instruments v1.0*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
