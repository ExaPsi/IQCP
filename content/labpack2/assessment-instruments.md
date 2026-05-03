# Lab Pack #2: Assessment Instruments

**Lab Pack:** 2 - 3D Exploration, PES, and Orbitals
**Version:** 1.0
**Last Updated:** 2026-03-17
**Document Type:** Assessment Portfolio Document
**Target Publication:** J. Chem. Educ. Technology Report

---

## Executive Overview

This document consolidates all assessment instruments for Lab Pack #2 "3D Exploration, PES, and Orbitals" of the Interactive Quantum Chemistry Playground (IQCP). The assessment portfolio measures student learning across six learning outcomes (LO7--LO12) spanning three Phase 2 modules: 3D molecular visualization, potential energy surface scanning, and orbital isosurface rendering.

### Assessment Philosophy

The assessment framework follows the same principles established in Lab Pack #1, extended to the spatial and visual domains of Phase 2:

1. **Constructive alignment:** All assessment items trace directly to stated learning outcomes LO7--LO12 (Wiggins & McTighe, 2005). Item stems reference authentic computational scenarios that students encounter in the IQCP worksheet activities.
2. **Multiple measures:** Conceptual understanding (concept checks), procedural and interpretive skills (worksheet), and integrative reasoning (performance tasks) are assessed through complementary instruments.
3. **Formative and summative:** Pre/post concept checks enable learning gain measurement; the worksheet provides formative guidance; performance tasks provide summative evaluation of applied competency.
4. **Misconception targeting:** Distractors in multiple-choice items are drawn from documented student misconceptions in the chemistry education literature (Tsaparlis & Papaphotis, 2009; Harle & Towns, 2011), ensuring diagnostic value beyond simple right/wrong scoring.

### Total Point Allocation

| Instrument | Points | Percentage | Purpose |
|------------|--------|------------|---------|
| Concept Check (Pre) | 5 | 3.2% | Baseline knowledge assessment |
| Concept Check (Post) | 5 | 3.2% | Learning gain measurement |
| Worksheet | 110 | 71.0% | Guided POE exploration assessment |
| Performance Tasks | 35 | 22.6% | Applied competency assessment |
| **Total Portfolio** | **155** | **100%** | |

### Learning Outcomes Assessed

| ID | Learning Outcome | Bloom's Level | Primary Instruments |
|----|------------------|---------------|---------------------|
| **LO7** | Geometry-energy connection: predict energy vs. distance; explain PES minimum | Analyze | CC (P1, Q1), WS (Q2.1-Q2.4), PT-PES |
| **LO8** | Orbital interpretation: distinguish bonding/nonbonding/antibonding MOs | Analyze | CC (P2, Q2), WS (Q3.1-Q3.5, Q3.8), PT-Orbital |
| **LO9** | Spatial-symbolic bridging: connect overlap matrix to 3D geometry | Apply | CC (P3, Q3), WS (Q1.3-Q1.8) |
| **LO10** | Basis set evaluation: cost-accuracy tradeoff, variational principle | Evaluate | CC (P4, Q4), WS (Q2.7-Q2.8) |
| **LO11a** | Dissociation limit awareness: identify RHF failure at large R | Analyze | CC (Q5), WS (Q2.5-Q2.6) |
| **LO11b** | Dissociation explanation: single-determinant limitations (graduate) | Evaluate | WS (Q2.6 extension) |
| **LO12** | Isovalue interpretation: explain isovalue meaning, predict shape changes | Understand | CC (P5, Q2 secondary), WS (Q3.6-Q3.7) |

---

## Pre-Activity Concept Check

**Administration:** Before students access IQCP. No IQCP use permitted during the check.
**Time:** 5--7 minutes
**Total Points:** 5 (1 point per item)
**Format:** 4 MC items + 1 SA item

The pre-check establishes baseline understanding of spatial, energetic, and orbital concepts that students bring from prior coursework (general chemistry, physical chemistry). Items target lower Bloom's levels (Remember, Understand) appropriate for a pre-assessment.

---

### Item P1: Energy vs. Bond Distance

**Primary LO:** LO7 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

Consider a diatomic molecule like H2. As you bring two hydrogen atoms closer together from a large separation, what happens to the total energy of the system?

**(a)** The energy decreases continuously -- the closer the atoms, the lower the energy.

**(b)** The energy first decreases to a minimum and then increases sharply at very short distances.

**(c)** The energy increases continuously -- bringing atoms together always costs energy.

**(d)** The energy remains constant regardless of the distance between the atoms.

---

### Item P2: What Is a Molecular Orbital?

**Primary LO:** LO8 | **Secondary LO:** LO12 | **Bloom's Level:** Remember | **Type:** MC

Which of the following best describes a molecular orbital?

**(a)** A fixed path that an electron follows as it orbits the nuclei, similar to a planet orbiting the sun.

**(b)** A mathematical function that describes the probability amplitude for finding an electron at each point in space within a molecule.

**(c)** A physical boundary surface that electrons cannot cross, enclosing the region where electrons exist.

**(d)** The region of space between two bonded atoms where electrons are shared equally.

---

### Item P3: Overlap Matrix Meaning

**Primary LO:** LO9 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

In an SCF calculation, the overlap matrix element S_ij measures the spatial overlap between basis functions i and j. If two atoms are bonded and close together, their basis function overlap S_ij is expected to be:

**(a)** Exactly zero, because different atoms have different basis functions.

**(b)** Exactly 1.0, because bonded atoms share electrons completely.

**(c)** A positive value between 0 and 1, with larger values for closer atoms.

**(d)** A negative value, because electrons repel each other.

---

### Item P4: Basis Set and Energy

**Primary LO:** LO10 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** MC

A student runs an SCF calculation on H2 using two different basis sets: STO-3G (2 basis functions) and 6-31G (4 basis functions). Which statement about the computed energies is correct?

**(a)** STO-3G will give a lower energy because fewer basis functions means less computational noise.

**(b)** 6-31G will give a lower energy because more basis functions provide greater variational freedom.

**(c)** Both basis sets will give exactly the same energy because they describe the same molecule.

**(d)** The larger basis set always gives a higher energy because it includes more electron repulsion terms.

---

### Item P5: What Does an Isovalue Represent?

**Primary LO:** LO12 | **Secondary LO:** -- | **Bloom's Level:** Understand | **Type:** SA

When visualizing a molecular orbital as a 3D isosurface, you must choose an "isovalue." In 1--2 sentences, explain what the isovalue represents and what the displayed surface shows.

*Expected response length: 1--2 sentences*

---

## Post-Activity Concept Check

**Administration:** Immediately after worksheet completion. No IQCP access during the check.
**Time:** 5--7 minutes
**Total Points:** 5 (1 point per item)
**Format:** 3 MC items + 2 SA items

The post-check assesses understanding at higher cognitive levels (Apply, Analyze, Evaluate) than the pre-check, reflecting the learning that should have occurred during the IQCP activities. Items target the same LOs as the pre-check but require deeper reasoning.

---

### Item Q1: Explaining the PES Minimum

**Primary LO:** LO7 | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** SA

You ran a PES scan for H2 and observed a minimum in the energy vs. bond length curve. In 2--3 sentences, explain why the energy minimum exists. Your explanation should reference at least two competing physical effects.

*Expected response length: 2--3 sentences*

---

### Item Q2: Classifying an Orbital from Its Description

**Primary LO:** LO8 | **Secondary LO:** LO12 | **Bloom's Level:** Apply | **Type:** MC

A student examines a molecular orbital of H2O and observes the following: the isosurface consists of two lobes on opposite sides of the oxygen atom, oriented perpendicular to the molecular plane. There is virtually no electron density on the hydrogen atoms. This orbital is best classified as:

**(a)** Bonding, because it is an occupied orbital and all occupied orbitals contribute to bonding.

**(b)** Antibonding, because it has two lobes with a nodal plane, which is the hallmark of antibonding character.

**(c)** Nonbonding, because electron density is localized on oxygen with negligible hydrogen contribution, indicating a lone pair.

**(d)** Core, because it is a small, tightly bound orbital that does not extend beyond the boundary of the oxygen atom.

---

### Item Q3: Predicting Overlap Change with Distance

**Primary LO:** LO9 | **Secondary LO:** -- | **Bloom's Level:** Apply | **Type:** MC

You are viewing a 3D model of H2 alongside its 2x2 overlap matrix. The off-diagonal element S_1,2 = 0.66 at the equilibrium bond length. If you increase the bond length by 50%, what do you predict will happen to S_1,2?

**(a)** S_1,2 will increase toward 1.0, because the basis functions will have more space to overlap.

**(b)** S_1,2 will decrease toward 0.0, because the basis functions are centered further apart and their spatial overlap diminishes.

**(c)** S_1,2 will remain approximately 0.66, because the overlap integral depends only on the type of orbital, not on the distance.

**(d)** S_1,2 will become negative, because at large distances the orbitals interfere destructively.

---

### Item Q4: Evaluating Basis Set Choice for a Research Scenario

**Primary LO:** LO10 | **Secondary LO:** -- | **Bloom's Level:** Evaluate | **Type:** SA

A researcher needs to compute the SCF energy for 500 different molecular geometries of a small molecule to map its potential energy surface. They must choose between STO-3G (fast, less accurate) and cc-pVTZ (slow, more accurate). In 2--3 sentences, advise the researcher on which basis set to use and justify your recommendation. Consider both accuracy and computational cost.

*Expected response length: 2--3 sentences*

---

### Item Q5: Identifying RHF Limitation from PES Behavior

**Primary LO:** LO11a | **Secondary LO:** -- | **Bloom's Level:** Analyze | **Type:** MC

A student computes the PES for H2 using RHF/STO-3G and observes that at very large bond distances (R > 4 bohr), the energy is -0.69 Ha rather than the expected -0.93 Ha for two isolated hydrogen atoms. What is the most likely explanation for this discrepancy?

**(a)** The STO-3G basis set is too small to describe hydrogen atoms at large distances. Using a larger basis set would eliminate the error.

**(b)** RHF forces alpha and beta electrons to occupy the same spatial orbital, which produces an incorrect mixture of covalent and ionic configurations at large R.

**(c)** The SCF calculation did not converge properly at large R, so the energy is simply a numerical artifact.

**(d)** RHF always gives the correct dissociation limit for any molecule; the expected value of -0.93 Ha must be wrong.

---

## Alignment Matrix

### Item-to-LO Mapping

| Item | Content Focus | Primary LO | Secondary LO | Bloom's Level | Item Type | Points |
|------|---------------|------------|--------------|---------------|-----------|--------|
| P1 | Energy vs. bond distance concept | LO7 | -- | Understand | MC | 1 |
| P2 | Molecular orbital definition | LO8 | LO12 | Remember | MC | 1 |
| P3 | Overlap matrix interpretation | LO9 | -- | Understand | MC | 1 |
| P4 | Basis set effect on energy | LO10 | -- | Understand | MC | 1 |
| P5 | Isovalue meaning | LO12 | -- | Understand | SA | 1 |
| Q1 | Physical origin of PES minimum | LO7 | -- | Analyze | SA | 1 |
| Q2 | Orbital classification from shape | LO8 | LO12 | Apply | MC | 1 |
| Q3 | Overlap change prediction | LO9 | -- | Apply | MC | 1 |
| Q4 | Basis set cost-accuracy tradeoff | LO10 | -- | Evaluate | SA | 1 |
| Q5 | RHF dissociation failure diagnosis | LO11a | -- | Analyze | MC | 1 |

### LO Coverage Verification

| LO | Pre-Check Items | Post-Check Items | Total Items | Coverage Status |
|----|-----------------|------------------|-------------|-----------------|
| LO7 | P1 | Q1 | 2 | Adequate |
| LO8 | P2 | Q2 | 2 | Adequate |
| LO9 | P3 | Q3 | 2 | Adequate |
| LO10 | P4 | Q4 | 2 | Adequate |
| LO11a | -- | Q5 | 1 | Adequate (post-only by design) |
| LO11b | -- | -- | 0 | Assessed via worksheet Q2.6 extension only |
| LO12 | P5 | Q2 (secondary) | 1 + 1 secondary | Adequate |

**Design note on LO11a:** The pre-check intentionally omits an RHF dissociation item because students are not expected to know about this limitation before the activity. Including a pre-item would be invalid (testing knowledge not yet taught) and could prime students to look for the answer rather than discover it through the PES scan. The single post-check item (Q5) combined with worksheet items Q2.5--Q2.6 provides adequate assessment coverage.

**Design note on LO11b:** This is a graduate extension objective. It is assessed through the worksheet (Q2.6 extended response) and discussed in the instructor key but is not included in the concept check, which is designed for the undergraduate population.

---

## Detailed Item-to-LO Mapping

### LO7 -- Geometry-Energy Connection

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | P1 (Understand, MC) | Tests baseline understanding of how energy varies with distance |
| Post-check | Q1 (Analyze, SA) | Requires explanation of competing physical effects producing the minimum |
| Worksheet | Q2.1--Q2.4 (4 items) | Full POE cycle: predict PES shape, observe scan, explain minimum and repulsive wall |
| Performance | PT-PES (4-point rubric) | Interpret a novel PES curve and identify key features |

**Cognitive progression:** Understand (P1) -> Analyze (Q1, worksheet) -> Analyze/Evaluate (PT-PES)

### LO8 -- Orbital Interpretation

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | P2 (Remember, MC) | Tests recall of molecular orbital definition; misconception distractor: "orbital = orbit" |
| Post-check | Q2 (Apply, MC) | Requires classifying an orbital from descriptive features; misconception distractor: "two lobes = antibonding" |
| Worksheet | Q3.1--Q3.5, Q3.8 (6 items) | Core, bonding, nonbonding, and antibonding orbital examination and classification |
| Performance | PT-Orbital (4-point rubric) | Classify orbitals in a novel molecular system |

**Cognitive progression:** Remember (P2) -> Apply (Q2) -> Analyze (worksheet, PT-Orbital)

### LO9 -- Spatial-Symbolic Bridging

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | P3 (Understand, MC) | Tests understanding of what overlap matrix elements represent |
| Post-check | Q3 (Apply, MC) | Requires predicting how overlap changes with distance |
| Worksheet | Q1.3--Q1.8 (6 items) | Full POE cycle connecting matrix elements to 3D molecular structure |

**Cognitive progression:** Understand (P3) -> Apply (Q3) -> Apply/Analyze (worksheet)

### LO10 -- Basis Set Evaluation

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | P4 (Understand, MC) | Tests understanding of variational principle in basis set context; misconception distractor: "bigger basis always better" |
| Post-check | Q4 (Evaluate, SA) | Requires weighing cost vs. accuracy for a practical scenario |
| Worksheet | Q2.7--Q2.8 (2 items) | Compare basis sets, evaluate "always use biggest basis" claim |

**Cognitive progression:** Understand (P4) -> Evaluate (Q4, worksheet)

### LO11a -- Dissociation Limit Awareness

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | -- | Not assessed pre-activity (students discover the limitation during the lab) |
| Post-check | Q5 (Analyze, MC) | Tests ability to diagnose RHF dissociation failure from PES data |
| Worksheet | Q2.5--Q2.6 (2 items) | Students observe and explain the dissociation limit problem |

**Cognitive progression:** Discovery (worksheet) -> Analyze (Q5)

### LO12 -- Isovalue Interpretation

| Assessment | Items | Rationale |
|------------|-------|-----------|
| Pre-check | P5 (Understand, SA) | Tests baseline understanding of isovalue concept |
| Post-check | Q2 (secondary) | Distractor in Q2 targets "orbital boundary" misconception related to LO12 |
| Worksheet | Q3.6--Q3.7 (2 items) | Predict and observe isovalue changes; confront "orbital = solid object" misconception |

**Cognitive progression:** Understand (P5) -> Understand/Apply (worksheet)

---

## Concept Check Answer Key

### Pre-Check Answers

#### P1: Energy vs. Bond Distance

**Correct answer: (b)** The energy first decreases to a minimum and then increases sharply at very short distances.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses "closer = more bonding = more stable" without considering nuclear repulsion | Incomplete understanding of bonding energetics |
| **(b)** | **Correct.** Recognizes the balance between attractive and repulsive interactions | -- |
| (c) | May recall electron repulsion but not electron delocalization stabilization | Overgeneralization of repulsion |
| (d) | Fails to recognize any interaction between atoms | Fundamental misunderstanding |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P2: What Is a Molecular Orbital?

**Correct answer: (b)** A mathematical function describing the probability amplitude for finding an electron at each point in space.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Conflates "orbital" with "orbit"; common in general chemistry students | **"Orbitals are like planetary orbits"** |
| **(b)** | **Correct.** Identifies the probabilistic, mathematical nature of MOs | -- |
| (c) | Treats the orbital isosurface as a physical boundary that electrons cannot cross | **"Orbitals have sharp boundaries"** |
| (d) | Restricts MOs to the bonding region; ignores antibonding, nonbonding, and lone pair MOs | Incomplete orbital concept |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P3: Overlap Matrix Meaning

**Correct answer: (c)** A positive value between 0 and 1, with larger values for closer atoms.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses "different atoms" with "no overlap"; fails to understand that Gaussians extend in space | Orthogonality confusion |
| (b) | Confuses normalization (diagonal = 1.0) with off-diagonal overlap | Diagonal vs. off-diagonal confusion |
| **(c)** | **Correct.** Recognizes distance dependence and bounded range | -- |
| (d) | Confuses overlap (always non-negative for s-type functions at reasonable distances) with repulsion | Sign confusion |

**Scoring:** 1 point for (c), 0 for all others.

---

#### P4: Basis Set and Energy

**Correct answer: (b)** 6-31G gives a lower energy because more basis functions provide greater variational freedom.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Assumes fewer basis functions means "cleaner" calculation | Misunderstanding of variational principle |
| **(b)** | **Correct.** Applies the variational principle correctly | -- |
| (c) | Fails to distinguish molecular identity from basis set representation | **"Bigger basis always gives the same answer"** |
| (d) | Confuses more basis functions with more repulsion terms | Conflation of basis size with physics content |

**Scoring:** 1 point for (b), 0 for all others.

---

#### P5: What Does an Isovalue Represent?

**Correct answer (exemplar):** The isovalue is a threshold value of the orbital wavefunction magnitude |psi|. The displayed surface connects all points in 3D space where the wavefunction equals that threshold value, showing where the orbital amplitude reaches a specified level.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies the isovalue as a threshold or cutoff value AND connects the surface to points where |psi| (or electron density) equals that value. Partial credit for identifying it as a threshold without specifying what quantity is thresholded. |
| 0 | Response is absent, incoherent, or describes the isovalue as a physical boundary, orbital size, or energy level. |

**Common incorrect responses:**
- "The isovalue is the energy of the orbital." (Confuses isovalue with orbital energy.)
- "The isovalue controls how big the orbital is." (Conflates the visualization parameter with a physical property of the orbital.)
- "The isovalue is the probability of finding an electron inside the surface." (Closer, but confuses the wavefunction threshold with integrated probability.)

---

### Post-Check Answers

#### Q1: Explaining the PES Minimum

**Correct answer (exemplar):** The PES minimum exists because of a balance between two competing effects. At intermediate distances, electron delocalization over both nuclei stabilizes the system by lowering kinetic energy and increasing electron-nuclear attraction. At very short distances, nuclear-nuclear repulsion (which scales as 1/R) dominates and drives the energy up sharply. The equilibrium bond length is where these effects balance.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies at least two competing physical effects (e.g., nuclear repulsion at short R and electron stabilization at intermediate R) AND conveys that the minimum results from their balance. Accept: "attraction vs. repulsion," "delocalization vs. nuclear repulsion," "lowered kinetic energy vs. Coulomb repulsion." |
| 0 | Response names only one effect, attributes the minimum to a single cause (e.g., "atoms attract each other"), or is absent/incoherent. |

---

#### Q2: Classifying an Orbital from Its Description

**Correct answer: (c)** Nonbonding, because electron density is localized on oxygen with negligible hydrogen contribution, indicating a lone pair.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Assumes all occupied orbitals are bonding | "Occupied = bonding" misconception |
| (b) | Confuses the p-orbital two-lobe shape with antibonding character | **"Two lobes = antibonding"** (also related to **"orbitals have sharp boundaries"** -- students may focus on lobe count rather than electron density distribution) |
| **(c)** | **Correct.** Identifies lone pair character from localization on O and absence of H contribution | -- |
| (d) | Confuses a valence lone pair orbital with a core orbital; the description explicitly states two lobes, which is not characteristic of a core 1s orbital | Core vs. valence confusion |

**Scoring:** 1 point for (c), 0 for all others.

---

#### Q3: Predicting Overlap Change with Distance

**Correct answer: (b)** S_1,2 will decrease toward 0.0, because the basis functions are centered further apart and their spatial overlap diminishes.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Confuses "more space" with "more overlap"; does not understand that overlap measures the product integral of two functions | Spatial reasoning error |
| **(b)** | **Correct.** Correctly predicts that increasing distance decreases the overlap integral | -- |
| (c) | Believes overlap is an intrinsic property of the orbital type, independent of geometry | Distance-independence misconception |
| (d) | Confuses overlap (always non-negative for same-type s orbitals) with interference | Sign confusion |

**Scoring:** 1 point for (b), 0 for all others.

---

#### Q4: Evaluating Basis Set Choice

**Correct answer (exemplar):** For 500 single-point calculations mapping a PES, I would recommend starting with STO-3G. The qualitative shape of the PES (location of minima, saddle points) is usually captured even with a minimal basis, and the 16-fold cost reduction (from N^4 scaling with twice as many basis functions) is significant when multiplied over 500 geometries. If quantitative accuracy is needed for specific points (e.g., the equilibrium energy), the researcher could re-compute those select geometries with cc-pVTZ.

**Scoring rubric:**

| Score | Criteria |
|-------|----------|
| 1 | Response identifies the cost-accuracy tradeoff AND makes a justified recommendation. Must mention computational cost scaling or time as a factor. Accept either recommendation (STO-3G for efficiency, cc-pVTZ for accuracy, or a two-stage strategy) as long as the justification is internally consistent. |
| 0 | Response recommends a basis set without mentioning cost, states "always use the biggest basis," or is absent/incoherent. |

---

#### Q5: Identifying RHF Limitation from PES Behavior

**Correct answer: (b)** RHF forces alpha and beta electrons to occupy the same spatial orbital, which produces an incorrect mixture of covalent and ionic configurations at large R.

| Choice | Why students select it | Misconception targeted |
|--------|----------------------|----------------------|
| (a) | Blames the basis set rather than the method; this error persists with any basis in RHF | **"Bigger basis always fixes the problem"** (variant of basis set misconception) |
| **(b)** | **Correct.** Identifies the single-determinant restriction as the root cause | -- |
| (c) | Attributes a systematic error to a numerical artifact; the energy at large R is converged, not artifactual | Confusing method limitation with numerical error |
| (d) | Denies the existence of RHF limitations entirely | **"RHF always gives correct dissociation"** |

**Scoring:** 1 point for (b), 0 for all others.

---

## Performance Task Rubrics

The following performance tasks may be administered as part of the worksheet debrief, as a separate assessment session, or as take-home assignments. Each task is scored on a 4-point analytic rubric.

### Performance Task 1: PES Interpretation (PT-PES)

**Target LOs:** LO7 (primary), LO11a (secondary)
**Total Points:** 21 (3 dimensions x 7-point scale, described below as a 4-point rubric per dimension)

#### Task Description

A student is given the following PES data for LiH (STO-3G) computed with IQCP:

| R (bohr) | Energy (Ha) | Converged? |
|-----------|-------------|------------|
| 1.0 | -7.578 | Yes |
| 1.5 | -7.826 | Yes |
| 2.0 | -7.862 | Yes |
| 2.5 | -7.863 | Yes |
| 3.0 | -7.855 | Yes |
| 4.0 | -7.837 | Yes |
| 6.0 | -7.814 | Yes |
| 8.0 | -7.802 | Yes |

The correct dissociation limit for Li(2S) + H(2S) is approximately -7.797 Ha.

**Part A (7 points):** Identify the equilibrium bond length and energy from this data. Explain why the energy minimum occurs at that distance.

**Part B (7 points):** Describe the behavior of the PES at large R. Does the RHF energy approach the correct dissociation limit? If not, explain why.

**Part C (7 points):** A colleague suggests that switching to a larger basis set (e.g., 6-31G*) would fix the incorrect dissociation behavior. Evaluate this claim.

#### PT-PES Rubric

**Part A -- Equilibrium Identification and Explanation (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplar) | Correctly identifies R_eq near 2.0--2.5 bohr from the data. Explains that the minimum arises from the balance between nuclear repulsion at short R and loss of electron stabilization at large R. Mentions electron delocalization or electron-nuclear attraction as the stabilizing effect. |
| 3 (Proficient) | Correctly identifies R_eq. Provides a partial explanation referencing either repulsion at short R or stabilization at intermediate R, but not both in a balanced way. |
| 2 (Developing) | Identifies a minimum in the data but places it at the wrong R value, or provides only a vague explanation ("atoms like to be at a certain distance"). |
| 1 (Beginning) | Attempts to identify a feature of the PES but does not identify the minimum correctly or provides no physical explanation. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part B -- Dissociation Limit Analysis (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplar) | Notes that the RHF energy at large R (-7.802 Ha) is close to but slightly above the correct limit (-7.797 Ha). Correctly identifies RHF's single-determinant constraint as the source of error. Recognizes that for LiH, the dissociation error is smaller than for H2 because the ionic Li+ + H- configuration is lower in energy than H+ + H-. |
| 3 (Proficient) | Notes the energy at large R and compares to the correct limit. Identifies RHF limitation but does not explain the mechanism or note the molecule-dependent magnitude. |
| 2 (Developing) | Notes that the curve levels off at large R but does not compare to the correct limit or attributes the behavior to basis set quality alone. |
| 1 (Beginning) | Does not analyze the large-R behavior or incorrectly states that RHF gives the correct limit. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part C -- Basis Set Claim Evaluation (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplar) | Correctly evaluates the claim as incorrect. Explains that the RHF dissociation error is a method limitation (single determinant), not a basis set limitation. Notes that a larger basis set would lower the energy at all R (variational principle) but would not correct the qualitative dissociation behavior. May suggest UHF or multi-reference methods as alternatives. |
| 3 (Proficient) | Correctly disagrees with the claim and identifies it as a method limitation rather than a basis set issue. Does not elaborate on what would fix it. |
| 2 (Developing) | Partially agrees with the claim (e.g., "it would help but not fully fix it") or correctly disagrees but with weak justification. |
| 1 (Beginning) | Agrees with the claim or provides no evaluation. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

---

### Performance Task 2: Orbital Classification (PT-Orbital)

**Target LOs:** LO8 (primary), LO12 (secondary)
**Total Points:** 14 (2 dimensions x 7-point scale)

#### Task Description

A student runs an SCF calculation on H2 (STO-3G) using IQCP and examines the two molecular orbitals.

**Part A (7 points):** The student observes MO 1 (sigma_g) as a single, continuous isosurface encompassing both hydrogen nuclei. MO 2 (sigma_u*) shows two separate lobes, one on each atom, rendered with different opacity (positive lobe solid, negative lobe translucent). Using these observations:

1. Classify each MO as bonding or antibonding.
2. Explain how the visual features (number of lobes, electron density between nuclei, node presence) support your classification.
3. Predict: if you decrease the isovalue from 0.05 to 0.01, how would the appearance of each orbital change?

**Part B (7 points):** The student then adjusts the isovalue slider for MO 1.

1. At isovalue 0.01, the isosurface is very large, extending far from the nuclei. At isovalue 0.08, the isosurface is small and concentrated near the nuclei. Does MO 1 have a physical edge? Justify your answer.
2. A classmate says: "The orbital is the colored shape on the screen. Electrons can only exist inside that shape." Evaluate this claim.

#### PT-Orbital Rubric

**Part A -- Orbital Classification and Isovalue Prediction (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplar) | Correctly classifies MO 1 as bonding and MO 2 as antibonding. Cites all three visual features (lobe count, internuclear density, node) to support classification. Correctly predicts that decreasing isovalue expands both isosurfaces, and explains that the threshold is lowered so more of the wavefunction is enclosed. |
| 3 (Proficient) | Correct classification of both MOs with at least two supporting visual features. Correct isovalue prediction but with incomplete explanation. |
| 2 (Developing) | Correct classification of one MO, or both classified correctly but with only one supporting feature. Isovalue prediction incorrect or absent. |
| 1 (Beginning) | Attempts classification but both incorrect, or correct classification with no visual evidence. No isovalue prediction. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

**Part B -- Isovalue Interpretation and Misconception Evaluation (7 points)**

| Score | Criteria |
|-------|----------|
| 4 (Exemplar) | Correctly states MO 1 has no physical edge. Explains that the isosurface is a chosen threshold, not a physical boundary, and that the wavefunction extends to infinity (decaying exponentially). Correctly evaluates the classmate's claim as wrong, explaining that electrons have nonzero probability of being found at any finite distance from the nuclei. |
| 3 (Proficient) | States no physical edge and identifies the isosurface as a threshold. Evaluates the classmate's claim as wrong but with incomplete justification. |
| 2 (Developing) | Recognizes that the isosurface changes with isovalue (implying no fixed edge) but does not explicitly state "no physical edge." Evaluation of the classmate's claim is vague. |
| 1 (Beginning) | States the orbital has a physical edge, agrees with the classmate, or provides no analysis. |

**Scoring conversion:** 4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1

---

## Cognitive Level Distribution Analysis

### Distribution Across Concept Check Items (Pre + Post)

| Bloom's Level | Target % | Items | Count | Actual % |
|---------------|----------|-------|-------|----------|
| Remember | 10% | P2 | 1 | 10% |
| Understand | 40% | P1, P3, P4, P5 | 4 | 40% |
| Apply | 20% | Q2, Q3 | 2 | 20% |
| Analyze | 20% | Q1, Q5 | 2 | 20% |
| Evaluate | 10% | Q4 | 1 | 10% |

**Aggregate summary (using 3-level grouping):**

| Grouping | Target % | Actual % | Items |
|----------|----------|----------|-------|
| Remember/Understand | ~20% | 50% (5/10) | P1, P2, P3, P4, P5 |
| Apply/Analyze | ~50% | 40% (4/10) | Q1, Q2, Q3, Q5 |
| Evaluate | ~30% | 10% (1/10) | Q4 |

**Design note:** The aggregate distribution across concept checks alone skews toward Remember/Understand because the pre-check intentionally targets lower cognitive levels (baseline assessment). When the full portfolio is considered (concept checks + worksheet + performance tasks), the distribution shifts substantially toward Apply/Analyze and Evaluate:

### Distribution Across Full Portfolio

| Grouping | Target % | Worksheet (26 items) | Performance Tasks (5 parts) | Combined (41 items) | Combined % |
|----------|----------|---------------------|---------------------------|---------------------|------------|
| Remember/Understand | ~20% | 6 items (Q1.1, Q1.2, Q1.4, Q3.1, Q3.6, P5 equivalent) | 0 | 11 | 27% |
| Apply/Analyze | ~50% | 14 items (Q1.3, Q1.5-Q1.8, Q2.1-Q2.6, Q3.2-Q3.5) | 3 parts (PT-PES A&B, PT-Orbital A) | 21 | 51% |
| Evaluate | ~30% | 6 items (Q2.7, Q2.8, Q3.7, Q3.8, Q4.1, Q4.2) | 2 parts (PT-PES C, PT-Orbital B) | 9 | 22% |

The full portfolio achieves approximately 27% Remember/Understand, 51% Apply/Analyze, and 22% Evaluate -- reasonably close to the 20/50/30 target. The slight elevation in Remember/Understand reflects the intentional use of lower-level items in the pre-check and prediction steps of the POE framework.

---

## Point Allocation

### Summary by Instrument

| Instrument | Points | Percentage |
|------------|--------|------------|
| Pre-Activity Concept Check (P1--P5) | 5 | 3.2% |
| Post-Activity Concept Check (Q1--Q5) | 5 | 3.2% |
| Student Worksheet (Q1.1--Q4.2) | 110 | 71.0% |
| Performance Task: PES Interpretation (PT-PES) | 21 | 13.5% |
| Performance Task: Orbital Classification (PT-Orbital) | 14 | 9.0% |
| **Total** | **155** | **100%** |

### Worksheet Points by Section

| Section | Questions | Points | Primary LOs |
|---------|-----------|--------|-------------|
| Section 1: 3D Molecular Exploration | Q1.1--Q1.8 | 30 | LO9 |
| Section 2: PES Scanning | Q2.1--Q2.8 | 36 | LO7, LO10, LO11 |
| Section 3: Orbital Visualization | Q3.1--Q3.8 | 34 | LO8, LO12 |
| Section 4: Synthesis | Q4.1--Q4.2 | 10 | Integrative |
| **Worksheet Total** | **26 items** | **110** | |

### Performance Task Points by Dimension

| Task | Part | Points | Primary LO |
|------|------|--------|------------|
| PT-PES | A: Equilibrium identification | 7 | LO7 |
| PT-PES | B: Dissociation limit analysis | 7 | LO11a |
| PT-PES | C: Basis set claim evaluation | 7 | LO10 |
| PT-Orbital | A: Classification + isovalue prediction | 7 | LO8 |
| PT-Orbital | B: Isovalue interpretation + misconception | 7 | LO12 |
| **Performance Total** | | **35** | |

### Recommended Grade Weighting

For courses using Lab Pack #2 as a graded assignment:

```
Lab Grade = (0.70 x Worksheet%) + (0.07 x ConceptCheck%) + (0.23 x PerformanceTask%)
```

This weighting reflects the centrality of the guided exploration (worksheet) while ensuring that higher-order reasoning (performance tasks) contributes meaningfully to the grade.

---

## Validity Argument (Kane, 2006)

This section presents the validity argument for the Lab Pack #2 assessment portfolio, structured according to Kane's (2006) argument-based approach to validation. The argument addresses three inferential links in the interpretive chain: scoring, generalization, and extrapolation.

### 1. Scoring Inference

**Claim:** Observed scores accurately reflect the quality of student responses.

**Evidence and warrants:**

*Multiple-choice items:*
- MC items (P1--P4, Q2, Q3, Q5) are scored dichotomously (1 or 0). Each item has one unambiguously correct answer verified against established quantum chemistry knowledge and cross-checked with PySCF reference calculations.
- Distractors were constructed from documented misconceptions in the chemistry education research literature (Tsaparlis & Papaphotis, 2009; Nakhleh, 1992), ensuring that incorrect responses carry diagnostic meaning rather than reflecting random guessing.

*Short-answer items:*
- SA items (P5, Q1, Q4) include explicit scoring rubrics with criteria for full credit (1 point) and no credit (0 points). Each rubric specifies the essential content elements required.
- Sample responses at each score level are provided for scorer calibration (see Concept Check Answer Key section above).

*Performance tasks:*
- Both PT-PES and PT-Orbital use 4-point analytic rubrics with dimension-specific criteria. Each rubric level includes concrete behavioral descriptors that minimize subjective judgment.
- Scoring conversion formulas (4 -> 7, 3 -> 5, 2 -> 3, 1 -> 1) produce point values with deliberate spacing to differentiate performance levels.

*Inter-rater reliability protocol:*
- For SA items and performance tasks, the recommended protocol is double-blind scoring of a random 20% sample, with adjudication of disagreements exceeding 1 rubric level. Target inter-rater reliability: Cohen's kappa >= 0.70.

**Potential threats:**
- SA items scored as 0/1 may lack sensitivity for partial understanding. This is mitigated by the complementary worksheet items (which use multi-point rubrics) assessing the same LOs.
- Scorer drift over large batches. Mitigated by clear rubric anchors and periodic calibration checks.

### 2. Generalization Inference

**Claim:** Scores obtained from these particular items generalize to the broader content domain defined by LO7--LO12.

**Evidence and warrants:**

*Content coverage:*
- The alignment matrix (above) demonstrates that every LO (LO7--LO12) is assessed by at least one concept check item, with the exception of LO11b, which is a graduate extension assessed only through the worksheet.
- Across the full portfolio (concept checks + worksheet + performance tasks), each LO is assessed by 2--8 items spanning multiple cognitive levels and item formats.
- The content domain is well-defined by the learning objectives, which were derived from the Phase 2 PRD and aligned to specific IQCP features (3D viewer, PES scanner, orbital visualizer).

*Cognitive level coverage:*
- Items span Remember through Evaluate, with the full portfolio distribution approximately matching the 20/50/30 target (actual: 27/51/22). This ensures that generalization extends across cognitive levels, not just content topics.

*Item sampling:*
- The concept check contains 10 items sampling from 6 LOs. While individual LO sub-scale scores (1--2 items each) are unreliable, the aggregate score across all items provides a reliable measure of overall learning. This is consistent with the portfolio-level gain target (normalized gain >= 0.3 across all LO7--LO12 items combined).

*Reliability targets:*
- Target Cronbach's alpha >= 0.70 for the combined concept check (10 items). Individual sub-scale reliability is not targeted due to the small number of items per LO.

**Potential threats:**
- With only 10 concept check items, the content sampling may not fully represent the breadth of each LO. This is mitigated by the 26-item worksheet, which provides dense coverage of each LO with multiple items per topic.
- Item difficulty may cluster (too easy or too hard), reducing discrimination. Post-pilot item analysis should examine difficulty indices and flag items outside the 0.30--0.90 range.

### 3. Extrapolation Inference

**Claim:** Performance on these instruments indicates understanding that transfers to broader quantum chemistry contexts beyond IQCP.

**Evidence and warrants:**

*Authentic tasks:*
- Concept check items use scenarios that parallel authentic computational chemistry reasoning (interpreting PES curves, classifying orbitals, evaluating basis set choices). These are not "IQCP navigation" questions -- they test conceptual understanding that applies to any quantum chemistry software.
- Performance tasks (PT-PES, PT-Orbital) present data in a format consistent with computational chemistry practice (energy tables, isosurface descriptions) and require reasoning that transfers to research contexts.

*Misconception targeting:*
- Items specifically target misconceptions documented in the broader QC education literature ("orbitals as orbits," "orbital boundaries are physical," "bigger basis always better," "RHF always correct"). Reduction in misconception endorsement is evidence of conceptual change that extends beyond the specific IQCP context.

*Connection to established assessments:*
- The assessment design follows best practices from validated chemistry concept inventories (e.g., ACS standardized exams, Quantum Chemistry Concept Inventory -- see Tsaparlis, 2005). While the specific items are novel, their structure and cognitive demands are consistent with established instruments in the field.

*Transfer evidence (planned):*
- A future properly powered classroom pilot (N >= 30) should include at least 2--3 transfer items that present novel molecular systems not encountered in the IQCP activities (e.g., "Given the PES for N2, identify the equilibrium and predict whether RHF dissociation is correct"). This would provide direct evidence of extrapolation.

**Potential threats:**
- Students may learn to interpret IQCP-specific representations without developing transferable understanding. Mitigated by the POE framework (which requires prediction before observation, engaging prior knowledge) and by including SA items that require verbal explanation rather than recognition.
- The controlled setting of a lab session may not predict performance in more open-ended contexts. This is an inherent limitation of timed assessments and is documented as a scope constraint.

### Limitations of the Validity Argument

1. **No pilot data yet.** This validity argument is based on content analysis and expert review, not empirical data. Item statistics (difficulty, discrimination, reliability) will be available after the first classroom administration.
2. **Single expert development.** Items were developed by one domain expert. Ideally, items should be reviewed by 2--3 additional experts and pilot-tested with a think-aloud protocol before publication.
3. **LO11b coverage.** The graduate extension objective (LO11b) is assessed only through the worksheet, not through the concept check or performance tasks. This limits the generalizability of conclusions about LO11b learning.
4. **Transfer evidence.** Direct evidence of extrapolation (transfer items, delayed post-tests) is planned for a future pilot but is not yet available.

---

## Administration Guidelines

### Pre-Check Administration

1. Distribute the pre-check (P1--P5) at the start of the session, **before** students open IQCP or any web browser.
2. Allow 5--7 minutes. Students should work individually without notes or discussion.
3. Collect all pre-check forms before distributing the worksheet or granting IQCP access.
4. **Important:** Do not review the pre-check answers at this point. Students should discover the concepts through the IQCP activities.

### Worksheet Administration

1. Distribute the worksheet and direct students to https://iqcp.dev.
2. Students work through Sections 1--4 at their own pace, following the POE prompts.
3. Target time: 60 minutes. Most students complete in 50--60 minutes.
4. Circulate to answer procedural questions (how to use IQCP controls) but avoid giving conceptual answers.
5. Collect worksheets (or confirm digital submission) before distributing the post-check.

### Post-Check Administration

1. Distribute the post-check (Q1--Q5) immediately after worksheet collection.
2. Allow 5--7 minutes. Students should work individually without IQCP access, notes, or worksheets.
3. Collect all post-check forms.
4. **Optional:** After collection, briefly review key concepts and invite questions. This is the teachable moment.

### Performance Task Administration

Performance tasks can be administered in three ways:

| Format | Timing | Recommended For |
|--------|--------|-----------------|
| Same session | After post-check (adds 15--20 min) | Extended lab periods (80+ min) |
| Separate session | Within 1 week of the lab | Standard lab periods |
| Take-home | Due within 1 week | Large enrollment courses |

For the same-session format, allow IQCP access during performance tasks (students may use IQCP to verify their reasoning). For take-home format, students should have IQCP access.

### Data Collection for Publication

When collecting data for J. Chem. Educ. publication:

1. **IRB approval** required before data collection. Use consent forms that allow publication of aggregate data with anonymized identifiers.
2. **Pre/post pairing:** Use anonymous identifiers (e.g., last 4 digits of student ID) to pair pre- and post-check responses for gain analysis.
3. **Record administration conditions:** Date, time, section size, any irregularities (e.g., technology failures, time extensions).
4. **Score all items independently** before entering into analysis. Double-score SA items and PTs for 20% of the sample.

### Recommended Sample Sizes

| Analysis Type | Minimum N | Recommended N |
|---------------|-----------|---------------|
| Descriptive statistics | 15 | 30+ |
| Pre/post paired t-test | 20 | 40+ |
| Normalized gain analysis | 25 | 50+ |
| Item-level analysis | 30 | 75+ |
| Cronbach's alpha estimation | 30 | 100+ |

### Data Quality Checks

Before analysis, verify:

- [ ] All students completed both pre and post checks
- [ ] No more than 10% missing data per item
- [ ] Pre-check administered before any IQCP exposure
- [ ] Post-check administered same day as the worksheet activity
- [ ] Inter-rater reliability >= 0.70 for short-answer and performance task items

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
| **Total** | | | | | |

#### Post-Check Item Analysis

| Item | N | Mean | SD | Difficulty | Discrimination |
|------|---|------|----|----|----------------|
| Q1 | | | | | |
| Q2 | | | | | |
| Q3 | | | | | |
| Q4 | | | | | |
| Q5 | | | | | |
| **Total** | | | | | |

#### Learning Gain Analysis

| Metric | Value |
|--------|-------|
| Pre-check mean (out of 5) | |
| Post-check mean (out of 5) | |
| Raw gain (Post - Pre) | |
| Normalized gain g | |
| Effect size (Cohen's d) | |
| Paired t-test p-value | |

**Normalized Gain Formula:**
```
g = (Post_mean - Pre_mean) / (5 - Pre_mean)
```

**Interpretation Guide:**
- g < 0.3: Low gain
- 0.3 <= g < 0.7: Medium gain
- g >= 0.7: High gain

### B. Worksheet Statistics

| Section | N | Mean | SD | Min | Max |
|---------|---|------|----|----|-----|
| Section 1: 3D Exploration (30 pts) | | | | | |
| Section 2: PES Scanning (36 pts) | | | | | |
| Section 3: Orbital Visualization (34 pts) | | | | | |
| Section 4: Synthesis (10 pts) | | | | | |
| **Total (110 pts)** | | | | | |

### C. Performance Task Statistics

| Task | N | Mean | SD | Min | Max |
|------|---|------|----|----|-----|
| PT-PES (21 pts) | | | | | |
| PT-Orbital (14 pts) | | | | | |
| **Total (35 pts)** | | | | | |

### D. Reliability Analysis

| Instrument | Cronbach's alpha | KR-20 | Notes |
|------------|------------------|-------|-------|
| Pre-check (5 items) | | | Target >= 0.60 |
| Post-check (5 items) | | | Target >= 0.60 |
| Combined CC (10 items) | | | Target >= 0.70 |
| Worksheet (26 items) | | N/A | Target >= 0.80 |

### E. Inter-Rater Reliability

| Item | Cohen's kappa | % Exact Agreement | Notes |
|------|---------------|-------------------|-------|
| P5 | | | Target >= 0.70 |
| Q1 | | | Target >= 0.70 |
| Q4 | | | Target >= 0.70 |
| PT-PES (all parts) | | | Target >= 0.70 |
| PT-Orbital (all parts) | | | Target >= 0.70 |

---

## References

Harle, M., and Towns, M.H. (2011). A review of spatial ability literature, its connection to chemistry, and implications for instruction. *Journal of Chemical Education*, 88(3), 351--360.

Kane, M.T. (2006). Validation. In R.L. Brennan (Ed.), *Educational Measurement* (4th ed., pp. 17--64). American Council on Education.

Nakhleh, M.B. (1992). Why some students don't learn chemistry: Chemical misconceptions. *Journal of Chemical Education*, 69(3), 191--196.

Tsaparlis, G. (2005). Non-algorithmic quantitative problem solving in university physical chemistry: A correlation study of the role of selective cognitive factors. *Research in Science and Technological Education*, 23(2), 125--148.

Tsaparlis, G., and Papaphotis, G. (2009). High-school students' conceptual difficulties and attempts at conceptual change: The case of basic quantum chemical concepts. *International Journal of Science Education*, 31(7), 895--930.

Wiggins, G., and McTighe, J. (2005). *Understanding by Design* (2nd ed.). Association for Supervision and Curriculum Development.

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-03-17 | IQCP Team | Initial release |

---

*Lab Pack #2 Assessment Instruments v1.0*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
