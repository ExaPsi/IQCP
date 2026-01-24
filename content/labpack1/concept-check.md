# Lab Pack #1: Pre/Post Concept Check Assessment

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18

---

## Administration Guidelines

### Purpose

This concept check measures student understanding before and after completing Lab Pack #1. The parallel structure of pre- and post-items enables direct comparison of learning gains.

### Administration Protocol

| Phase | Timing | Duration | Notes |
|-------|--------|----------|-------|
| **Pre-Check** | Before starting the lab | 10-12 minutes | Assesses prior knowledge |
| **Post-Check** | Immediately after lab | 10-12 minutes | Assesses learning gains |

### Instructions for Students

- Answer all questions to the best of your ability
- For multiple-choice items, select the BEST answer
- For short-answer items, write 1-3 sentences
- You may NOT use external resources or the IQCP tool during this assessment
- Time limit: approximately 10-12 minutes

### Scoring Summary

| Item Type | Count | Points Each | Total |
|-----------|-------|-------------|-------|
| Multiple Choice (MC) | 4 pre + 5 post = 9 | 1 point | 9 points |
| Short Answer (SA) | 2 pre + 2 post = 4 | 2 points | 8 points |
| **Pre-Check Total** | 6 items | - | **8 points** |
| **Post-Check Total** | 7 items | - | **9 points** |
| **Assessment Total** | 13 items | - | **17 points** |

---

## Pre-Activity Concept Check

*Administer BEFORE students begin the lab activity.*

### Instructions

Answer the following questions based on your current understanding. These questions assess background knowledge before you begin the interactive exploration.

---

### Item P1 (Multiple Choice - Quadrature Concepts)

**When evaluating a definite integral numerically using quadrature, increasing the number of quadrature points generally:**

A) Increases accuracy but increases computational cost

B) Decreases accuracy but increases computational cost

C) Increases accuracy and decreases computational cost

D) Has no effect on accuracy

---

### Item P2 (Multiple Choice - Iterative Methods)

**In an iterative numerical method, "convergence" means:**

A) The method runs faster with each iteration

B) Successive approximations approach a fixed answer

C) The method uses less memory over time

D) The initial guess was correct

---

### Item P3 (Multiple Choice - Gaussian Integrals)

**The integral of exp(-x^2) from 0 to infinity is:**

A) Undefined (the integral diverges)

B) A finite positive value

C) Zero

D) Negative infinity

---

### Item P4 (Multiple Choice - Hartree-Fock Theory)

**In quantum chemistry, the Hartree-Fock method calculates:**

A) The exact molecular energy

B) An approximate molecular energy (variational upper bound)

C) Only the kinetic energy of electrons

D) Only the nuclear repulsion energy

---

### Item P5 (Short Answer - Numerical Regimes)

**In one or two sentences, explain why numerical algorithms might use different computational methods (e.g., series expansion vs. asymptotic formula) for different parameter ranges.**

*Your answer:*

_______________________________________________

_______________________________________________

_______________________________________________

---

### Item P6 (Short Answer - Convergence Tolerance)

**What does it mean for an iterative calculation to "converge to a tolerance of 10^-6"? Explain in one or two sentences.**

*Your answer:*

_______________________________________________

_______________________________________________

_______________________________________________

---

## Post-Activity Concept Check

*Administer IMMEDIATELY AFTER students complete the lab activity.*

### Instructions

Answer the following questions based on what you learned during the interactive exploration. These questions assess your understanding after completing the lab.

---

### Item Q1 (Multiple Choice - Rys Quadrature)

**In Rys quadrature for molecular integrals, increasing the quadrature order from n=5 to n=8:**

A) Increases accuracy by allowing exact integration of higher-order polynomial moments

B) Decreases accuracy because more points introduce additional round-off error

C) Has no effect because the integration domain [0,1] is fixed

D) Increases accuracy only for small T values

---

### Item Q2 (Multiple Choice - DIIS Acceleration)

**DIIS (Direct Inversion in the Iterative Subspace) improves SCF convergence by:**

A) Making each individual iteration compute faster

B) Extrapolating an optimal Fock matrix from previous iterations

C) Reducing the number of basis functions used

D) Changing the molecular geometry during optimization

---

### Item Q3 (Multiple Choice - Boys Function Asymptotics)

**For large values of T, the Boys function F_m(T) approaches:**

A) Infinity

B) One

C) 1/(2m+1)

D) Zero

---

### Item Q4 (Multiple Choice - Method Selection)

**IQCP uses a Taylor series for F_m(T) when T is small, but switches to a recurrence relation (erf + upward recurrence) when T exceeds a threshold. This is primarily because:**

A) The Taylor series converges too slowly (requires too many terms) for large T

B) The recurrence relation is always more accurate than the series

C) The Taylor series is computationally faster

D) Large T values cause integer overflow errors

---

### Item Q5 (Multiple Choice - Root Count Rule)

**In Rys quadrature for electron repulsion integrals, the number of quadrature roots needed depends on the shell quartet. For a (pp|pp) shell quartet (L=4), how many roots are required according to the formula n_r = floor(L/2) + 1?**

A) 2

B) 3

C) 4

D) 5

---

### Item Q6 (Short Answer - Quadrature Order Selection)

**Based on your observations in the Rys module, explain how the optimal quadrature order depends on both the parameter T and the target accuracy. (2-3 sentences)**

*Your answer:*

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

### Item Q7 (Short Answer - DIIS Effect on Convergence)

**Describe in 2-3 sentences how DIIS changes the SCF convergence behavior compared to standard (non-accelerated) iteration. What specific differences did you observe?**

*Your answer:*

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

## Scoring Key

### Pre-Check Answers

| Item | Type | Correct Answer | Justification |
|------|------|----------------|---------------|
| **P1** | MC | **A** | More quadrature points improve accuracy but require more function evaluations, increasing cost. This is the fundamental trade-off in numerical integration. |
| **P2** | MC | **B** | Convergence means successive iterations produce values that stabilize toward a fixed result. The other options describe computational efficiency, not mathematical convergence. |
| **P3** | MC | **B** | The Gaussian integral evaluates to sqrt(pi)/2, a well-known finite positive value. This is fundamental to quantum chemistry because Gaussian basis functions have similar forms. |
| **P4** | MC | **B** | Hartree-Fock provides an approximate energy that is guaranteed to be at or above the true ground state energy (variational principle). It neglects electron correlation. |

**P5 Scoring Rubric (2 points):**

| Points | Criteria |
|--------|----------|
| 2 | Correctly explains that different methods work better in different ranges (e.g., series converges quickly for small values but slowly for large; asymptotic formulas are accurate for large values but diverge for small). May mention numerical stability, convergence rate, or computational efficiency. |
| 1 | Partially correct explanation that mentions different methods for different ranges but lacks clear reasoning about WHY (convergence, stability, or efficiency). |
| 0 | Incorrect or no meaningful response. |

**P6 Scoring Rubric (2 points):**

| Points | Criteria |
|--------|----------|
| 2 | Correctly explains that consecutive iterations differ by less than 10^-6 (or that the quantity of interest is within 10^-6 of its final value). May mention energy difference, density matrix change, or residual norm. |
| 1 | Mentions the tolerance threshold but does not clearly explain what quantity is being compared or how convergence is determined. |
| 0 | Incorrect or no meaningful response. |

---

### Post-Check Answers

| Item | Type | Correct Answer | Justification |
|------|------|----------------|---------------|
| **Q1** | MC | **A** | Higher quadrature order allows exact integration of polynomials up to degree 2n-1. This directly improves accuracy for the polynomial moments needed in integral evaluation. Round-off error (B) is negligible at these orders. |
| **Q2** | MC | **B** | DIIS combines Fock matrices from previous iterations using coefficients that minimize the error vector, producing an extrapolated Fock matrix closer to the converged solution. It does not speed up individual iterations (A) or change the basis (C) or geometry (D). |
| **Q3** | MC | **D** | As T increases, the exponential weight exp(-T*t^2) decays rapidly, suppressing the integrand and causing F_m(T) to approach zero. Students should observe this directly in the Boys module. |
| **Q4** | MC | **A** | Taylor series require many terms to converge when T is large because the function value is small but the series starts from large intermediate terms. The recurrence relation starting from erf(sqrt(T)) works efficiently for moderate-to-large T. |
| **Q5** | MC | **B** | For L=4, n_r = floor(4/2) + 1 = 2 + 1 = 3. This is the minimum number of quadrature roots needed to exactly integrate the polynomial terms in the (pp|pp) electron repulsion integral. |

**Q6 Scoring Rubric (2 points):**

| Points | Criteria |
|--------|----------|
| 2 | Correctly explains BOTH relationships: (1) Higher accuracy targets require higher quadrature order, AND (2) Larger T values generally require higher order to maintain the same accuracy (or notes that the optimal order depends on both factors). May include specific observations from the Rys module. |
| 1 | Correctly explains ONE relationship (either accuracy-order OR T-order) but not both, or explains both incompletely. |
| 0 | Incorrect or no meaningful response. |

**Q7 Scoring Rubric (2 points):**

| Points | Criteria |
|--------|----------|
| 2 | Correctly describes DIIS effects including: (1) Reduces number of iterations, AND (2) Produces smoother/faster convergence (avoids oscillations or slow linear convergence). May include specific iteration counts observed in the lab. |
| 1 | Mentions one aspect of DIIS improvement (fewer iterations OR smoother convergence) but not both, or description is vague. |
| 0 | Incorrect or no meaningful response. |

---

## Distractor Analysis

### Pre-Check Distractors

| Item | Distractor | Why Students Might Choose It |
|------|------------|------------------------------|
| P1-B | "Decreases accuracy but increases cost" | Confusion about numerical precision vs. mathematical accuracy |
| P1-C | "Increases accuracy and decreases cost" | Wishful thinking; not understanding trade-offs |
| P1-D | "No effect" | Lack of understanding of quadrature |
| P2-A | "Runs faster" | Confusing convergence with computational efficiency |
| P2-C | "Uses less memory" | Confusing convergence with resource usage |
| P2-D | "Initial guess correct" | Misunderstanding what convergence means |
| P3-A | "Diverges" | Incorrect assumption that exp(-x^2) does not decay fast enough |
| P3-C | "Zero" | Confusing integral of exp(-x^2) with integral of exp(-x) or other functions |
| P4-A | "Exact energy" | Not understanding that HF is an approximation |
| P4-C | "Kinetic only" | Confusing Hartree-Fock with a component of the Hamiltonian |
| P4-D | "Nuclear repulsion only" | Confusing HF with classical electrostatics |

### Post-Check Distractors

| Item | Distractor | Why Students Might Choose It |
|------|------------|------------------------------|
| Q1-B | "More points = more error" | Over-generalizing floating-point concerns |
| Q1-C | "No effect - domain fixed" | Not understanding how order affects polynomial accuracy |
| Q1-D | "Only for small T" | Incomplete understanding of T-dependence |
| Q2-A | "Faster iterations" | Confusing acceleration of convergence with speed per iteration |
| Q2-C | "Fewer basis functions" | Confusing DIIS with basis set truncation |
| Q5-A | "2 roots" | May compute L/2 without the +1 |
| Q5-C | "4 roots" | May think n_r = L instead of floor(L/2)+1 |
| Q5-D | "5 roots" | May add L to something incorrectly |
| Q2-D | "Changes geometry" | Confusing DIIS with geometry optimization |
| Q3-A | "Infinity" | Incorrect extrapolation of integral behavior |
| Q3-B | "One" | Confusing with F_m(0) limit |
| Q3-C | "1/(2m+1)" | This is F_m(0), not F_m(large T) |
| Q4-B | "Recurrence always better" | Not understanding when each method is appropriate |
| Q4-C | "Series is faster" | Series is NOT faster for large T |
| Q4-D | "Overflow errors" | Technically possible but not the primary reason |

---

## Learning Outcome Alignment

| Learning Outcome | Pre-Check Items | Post-Check Items |
|------------------|-----------------|------------------|
| **LO1:** Boys function behavior and numerical methods | P3 (Gaussian integral), P5 (numerical methods) | Q3 (asymptotics), Q4 (method selection) |
| **LO2:** Rys quadrature order-accuracy relationship | P1 (quadrature concepts) | Q1 (Rys order), Q5 (root count rule), Q6 (order selection) |
| **LO3:** SCF convergence and DIIS acceleration | P2 (iterative methods), P6 (convergence tolerance) | Q2 (DIIS mechanism), Q7 (DIIS effect) |
| **LO4:** Parameter-outcome connections | P5 (method rationale), P6 (tolerance meaning) | Q4 (method selection), Q6 (order-T-accuracy) |

---

## Administration Notes for Instructors

### Timing Recommendations

- Allow 10-12 minutes per phase (pre and post)
- If students finish early, have them review their answers
- If time is limited, consider reducing to 4 items per phase (drop P3/Q3 and P4/Q4)

### Common Student Difficulties

| Item | Common Issue | Intervention |
|------|--------------|--------------|
| P5/Q6 | Vague or incomplete responses | Prompt: "Be specific about WHAT changes and WHY" |
| P6/Q7 | Confusing tolerance with precision | Clarify: "Tolerance is about HOW CLOSE successive values are" |
| Q4 | Not understanding series convergence | In post-discussion, explain partial sum behavior |
| Q5 | Forgetting the +1 in the formula | Review n_r = floor(L/2) + 1; emphasize that even L=0 needs 1 root |

### Interpreting Results

**Learning Gain Calculation:**

```
Normalized Gain = (Post - Pre) / (Max - Pre)
```

| Gain Range | Interpretation |
|------------|----------------|
| g > 0.7 | High gain (excellent learning) |
| 0.3 < g < 0.7 | Medium gain (typical for interactive activities) |
| g < 0.3 | Low gain (may indicate issues with lab or assessment) |

**Item-Level Analysis:**

- Compare pre/post performance on parallel items (P1 vs Q1, etc.)
- Items with <50% post-check accuracy may indicate concepts needing additional instruction
- Items with high pre-check accuracy (>80%) may be too easy or assessing prior knowledge only

### Adaptations

**For Shorter Sessions:**
- Use 4-item version: P1, P2, P5, P6 (pre) and Q1, Q2, Q5, Q6 (post)
- Focus on quadrature and SCF concepts, skip Boys-specific items

**For Advanced Students:**
- Add follow-up prompts to short-answer items
- Require quantitative justification (e.g., "cite specific iteration counts")

**For Assessment Research:**
- Administer pre-check at least one class period before the lab
- Consider delayed post-check (1-2 weeks) to assess retention

---

## Printable Student Versions

### Pre-Check Student Form

```
LAB PACK #1: PRE-ACTIVITY CONCEPT CHECK
Name: _________________________ Date: _____________

Instructions: Answer all questions. MC: select best answer. SA: 1-3 sentences.
Time: ~10 minutes. No external resources.

P1. When evaluating a definite integral numerically using quadrature,
    increasing the number of quadrature points generally:

    A) Increases accuracy but increases computational cost
    B) Decreases accuracy but increases computational cost
    C) Increases accuracy and decreases computational cost
    D) Has no effect on accuracy

    Answer: [ ]

P2. In an iterative numerical method, "convergence" means:

    A) The method runs faster with each iteration
    B) Successive approximations approach a fixed answer
    C) The method uses less memory over time
    D) The initial guess was correct

    Answer: [ ]

P3. The integral of exp(-x^2) from 0 to infinity is:

    A) Undefined (the integral diverges)
    B) A finite positive value
    C) Zero
    D) Negative infinity

    Answer: [ ]

P4. In quantum chemistry, the Hartree-Fock method calculates:

    A) The exact molecular energy
    B) An approximate molecular energy (variational upper bound)
    C) Only the kinetic energy of electrons
    D) Only the nuclear repulsion energy

    Answer: [ ]

P5. In one or two sentences, explain why numerical algorithms might use
    different computational methods for different parameter ranges.

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________

P6. What does it mean for an iterative calculation to "converge to a
    tolerance of 10^-6"? Explain in one or two sentences.

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________
```

### Post-Check Student Form

```
LAB PACK #1: POST-ACTIVITY CONCEPT CHECK
Name: _________________________ Date: _____________

Instructions: Answer all questions. MC: select best answer. SA: 1-3 sentences.
Time: ~10 minutes. No external resources.

Q1. In Rys quadrature for molecular integrals, increasing the quadrature
    order from n=5 to n=8:

    A) Increases accuracy by allowing exact integration of higher-order
       polynomial moments
    B) Decreases accuracy because more points introduce additional
       round-off error
    C) Has no effect because the integration domain [0,1] is fixed
    D) Increases accuracy only for small T values

    Answer: [ ]

Q2. DIIS (Direct Inversion in the Iterative Subspace) improves SCF
    convergence by:

    A) Making each individual iteration compute faster
    B) Extrapolating an optimal Fock matrix from previous iterations
    C) Reducing the number of basis functions used
    D) Changing the molecular geometry during optimization

    Answer: [ ]

Q3. For large values of T, the Boys function F_m(T) approaches:

    A) Infinity
    B) One
    C) 1/(2m+1)
    D) Zero

    Answer: [ ]

Q4. IQCP uses a Taylor series for F_m(T) when T is small, but switches
    to a recurrence relation (erf + upward recurrence) when T exceeds
    a threshold. This is primarily because:

    A) The Taylor series converges too slowly (requires too many terms)
       for large T
    B) The recurrence relation is always more accurate than the series
    C) The Taylor series is computationally faster
    D) Large T values cause integer overflow errors

    Answer: [ ]

Q5. Based on your observations in the Rys module, explain how the optimal
    quadrature order depends on both the parameter T and the target
    accuracy. (2-3 sentences)

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________

Q6. Describe in 2-3 sentences how DIIS changes the SCF convergence
    behavior compared to standard (non-accelerated) iteration. What
    specific differences did you observe?

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________

    _________________________________________________________________
```

---

*IQCP Lab Pack #1 v1.0 | Concept Check Assessment | https://iqcp.dev*
