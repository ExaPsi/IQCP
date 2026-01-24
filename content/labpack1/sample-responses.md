# Lab Pack #1: Sample Student Responses and Scoring Guide

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18

---

## Purpose

This document provides sample student responses for the four short-answer items (P5, P6, Q5, Q6) in the Lab Pack #1 concept check assessment. Each item includes three response levels with detailed annotations explaining the scoring rationale. Instructors can use this guide to:

1. Calibrate their scoring across multiple graders
2. Understand the range of acceptable responses
3. Identify common misconceptions
4. Provide meaningful feedback to students

---

## Scoring Framework

### Point Distribution

| Score | Level | Description |
|-------|-------|-------------|
| **2 points** | Exemplary | Demonstrates complete understanding with accurate, well-articulated reasoning |
| **1 point** | Adequate | Shows partial understanding; may be incomplete, vague, or contain minor errors |
| **0 points** | Inadequate | Incorrect, irrelevant, or demonstrates fundamental misunderstanding |

### Scoring Principles

1. **Focus on conceptual understanding**, not perfect terminology
2. **Award partial credit** for responses that show some correct reasoning
3. **Do not penalize** minor grammatical or spelling errors
4. **Consider the context** of what students observed in the IQCP modules
5. **Look for key concepts** listed in the scoring annotations

---

## Pre-Check Item P5: Numerical Regimes

### Question

**In one or two sentences, explain why numerical algorithms might use different computational methods (e.g., series expansion vs. asymptotic formula) for different parameter ranges.**

### Key Concepts Expected

- Convergence rate varies with parameter value
- Numerical stability differs across ranges
- Efficiency/accuracy trade-offs
- Different approximations are accurate in different domains

---

### Exemplary Response (2 points) - Sample A

> "Different computational methods are optimal in different parameter ranges because the rate of convergence changes. A Taylor series converges rapidly when the parameter is small but needs many terms for large parameters, whereas an asymptotic expansion works best when the parameter is large and would diverge for small values."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | YES - explains convergence rate changes |
| Mentions accuracy/stability | YES - notes divergence for inappropriate ranges |
| Specific example | YES - references series vs. asymptotic appropriately |
| **Score** | **2 points** |

**Why this earns full credit:** The response correctly explains that convergence behavior depends on parameter magnitude, and provides accurate reasoning about when each method is appropriate. The student demonstrates understanding that methods have domains of validity.

---

### Exemplary Response (2 points) - Sample B

> "Numerical algorithms use different methods because each method has a range where it is most accurate and stable. For instance, series expansions work well near a starting point but become unstable or require too many terms far from that point, while asymptotic formulas give good approximations in limiting cases."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | YES - accuracy and stability vary |
| Mentions accuracy/stability | YES - explicitly states both |
| Specific example | YES - describes behavior of each method type |
| **Score** | **2 points** |

**Why this earns full credit:** The response identifies both accuracy and numerical stability as reasons, and correctly characterizes the domains of validity for different methods.

---

### Adequate Response (1 point) - Sample C

> "Different methods work better in different ranges. A series expansion is used when the parameter is small, and asymptotic formulas are used when it is large."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | PARTIAL - says "work better" but does not explain why |
| Mentions accuracy/stability | NO - lacks explanation |
| Specific example | YES - correctly identifies when each is used |
| **Score** | **1 point** |

**Why this earns partial credit:** The response correctly identifies THAT different methods are used in different ranges, but does not explain WHY one method is better in a given range. The reasoning is missing.

---

### Adequate Response (1 point) - Sample D

> "Using different methods improves efficiency and accuracy. Some methods converge faster in certain parameter ranges, so switching methods saves computation time."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | YES - mentions convergence speed |
| Mentions accuracy/stability | PARTIAL - mentions accuracy but not stability |
| Specific example | NO - does not distinguish series from asymptotic |
| **Score** | **1 point** |

**Why this earns partial credit:** The response correctly identifies efficiency and convergence as factors, but lacks specificity about what makes each method appropriate for its range.

---

### Inadequate Response (0 points) - Sample E

> "Computers use different methods because some algorithms are faster than others."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | NO - speed alone does not explain regime-dependent switching |
| Mentions accuracy/stability | NO |
| Specific example | NO |
| **Score** | **0 points** |

**Why this earns no credit:** The response misses the key concept that method choice depends on parameter value. It suggests methods are categorically faster, rather than conditionally appropriate.

---

### Inadequate Response (0 points) - Sample F

> "Different methods give different answers, so you pick the one that gives the answer you want."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Identifies WHY | NO - suggests arbitrary choice or answer manipulation |
| Mentions accuracy/stability | NO |
| Specific example | NO |
| **Score** | **0 points** |

**Why this earns no credit:** This response suggests that method choice is arbitrary or allows selecting preferred answers, which is a fundamental misconception. Correct numerical methods should converge to the same answer; the choice is about efficiency and numerical behavior, not outcome.

---

### Common Misconceptions for P5

| Misconception | Correct Understanding |
|--------------|----------------------|
| "One method is always better" | Each method has a domain where it excels |
| "Methods give different answers" | Correct methods converge to the same value |
| "It's about programming convenience" | It's about mathematical and numerical properties |
| "Faster computer eliminates need for regimes" | Convergence issues exist regardless of speed |

---

## Pre-Check Item P6: Convergence Tolerance

### Question

**What does it mean for an iterative calculation to "converge to a tolerance of 10^-6"? Explain in one or two sentences.**

### Key Concepts Expected

- Successive iterations produce values that differ by less than the tolerance
- OR: The computed quantity is within the tolerance of its final/true value
- Understanding that convergence is about stability/agreement, not absolute correctness

---

### Exemplary Response (2 points) - Sample A

> "It means that the difference between the results of consecutive iterations is less than 10^-6. When this threshold is reached, we consider the calculation finished because additional iterations would not significantly change the answer."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | YES - consecutive iteration difference |
| Explains implication | YES - stopping criterion, no significant change |
| Clear and complete | YES |
| **Score** | **2 points** |

**Why this earns full credit:** The response accurately defines what convergence to a tolerance means (difference between successive iterations) and explains the practical implication (stopping criterion).

---

### Exemplary Response (2 points) - Sample B

> "Converging to a tolerance of 10^-6 means that the iterative process continues until the change in the calculated quantity (such as energy) from one iteration to the next is smaller than 0.000001. This indicates the solution has stabilized."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | YES - change between iterations |
| Explains implication | YES - solution stabilization |
| Clear and complete | YES - includes example (energy) |
| **Score** | **2 points** |

**Why this earns full credit:** Provides accurate definition with a concrete example and correctly interprets the meaning (stabilization).

---

### Adequate Response (1 point) - Sample C

> "The calculation stops when the answer is accurate to 6 decimal places."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | PARTIAL - conflates tolerance with decimal precision |
| Explains implication | NO - does not explain convergence concept |
| Clear and complete | NO - missing key concept |
| **Score** | **1 point** |

**Why this earns partial credit:** The response demonstrates some understanding that 10^-6 relates to precision, but incorrectly equates tolerance with decimal place accuracy and misses the concept of comparing successive iterations.

---

### Adequate Response (1 point) - Sample D

> "It means the error in the calculation is less than 10^-6."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | PARTIAL - "error" is ambiguous |
| Explains implication | NO - does not clarify what error means |
| Clear and complete | NO - too vague |
| **Score** | **1 point** |

**Why this earns partial credit:** The response uses the word "error" which could mean the difference between iterations (correct) or the difference from the true value (problematic, since the true value is unknown). The ambiguity prevents full credit.

---

### Inadequate Response (0 points) - Sample E

> "The calculation runs for 10^6 iterations."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | NO - confuses tolerance with iteration count |
| Explains implication | NO |
| Clear and complete | NO |
| **Score** | **0 points** |

**Why this earns no credit:** Fundamental misunderstanding - tolerance is about the magnitude of change, not the number of iterations.

---

### Inadequate Response (0 points) - Sample F

> "The result is exactly correct to 6 significant figures."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Defines tolerance correctly | NO - conflates convergence with exact correctness |
| Explains implication | NO |
| Clear and complete | NO |
| **Score** | **0 points** |

**Why this earns no credit:** Convergence to a tolerance does NOT guarantee the answer is correct to that precision relative to the true value; it only means the iteration has stabilized. The response demonstrates a misconception about what convergence guarantees.

---

### Common Misconceptions for P6

| Misconception | Correct Understanding |
|--------------|----------------------|
| "Tolerance = number of iterations" | Tolerance is about change magnitude, not count |
| "Convergence means exact answer" | Convergence means stabilization, not correctness |
| "10^-6 means 6 decimal places" | 10^-6 is a threshold for change, not decimal precision |
| "Smaller tolerance is always better" | Smaller tolerance takes longer and may exceed numerical precision |

---

## Post-Check Item Q5: Quadrature Order Selection

### Question

**Based on your observations in the Rys module, explain how the optimal quadrature order depends on both the parameter T and the target accuracy. (2-3 sentences)**

### Key Concepts Expected

1. Higher accuracy requirements need higher quadrature order
2. The relationship between T and optimal order (larger T may need higher order)
3. References to actual observations in the Rys module

---

### Exemplary Response (2 points) - Sample A

> "From the Rys module, I observed that higher target accuracy requires higher quadrature order because more points are needed to exactly integrate the polynomial moments. Additionally, for larger values of T, higher order is generally needed to maintain the same accuracy because the integrand behavior changes. For example, at T=5 with order 3, I saw errors around 10^-8, but at T=20 with the same order, errors were larger."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | YES - clearly explained |
| T-order relationship | YES - correctly describes that larger T needs higher order |
| Module observations | YES - cites specific example with values |
| **Score** | **2 points** |

**Why this earns full credit:** The response addresses both required relationships (accuracy-order and T-order) and includes specific observations from the module, demonstrating genuine engagement with the material.

---

### Exemplary Response (2 points) - Sample B

> "The Rys module showed that quadrature order and accuracy are directly related - doubling the order roughly squares the accuracy for polynomial integrands. I also noticed that when T increases, the same quadrature order produces slightly larger errors, so higher T values benefit from increasing the order. The error curves in the module made this tradeoff visible."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | YES - quantitative observation |
| T-order relationship | YES - correctly identifies T-error relationship |
| Module observations | YES - references error curves |
| **Score** | **2 points** |

**Why this earns full credit:** Demonstrates understanding of both relationships with quantitative insight and references specific module features.

---

### Adequate Response (1 point) - Sample C

> "Higher accuracy needs more quadrature points. I saw in the module that increasing the order from 3 to 5 reduced the error significantly."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | YES - correctly stated |
| T-order relationship | NO - not addressed |
| Module observations | YES - includes specific observation |
| **Score** | **1 point** |

**Why this earns partial credit:** Correctly explains the accuracy-order relationship with a module observation, but does not address how T affects the optimal order.

---

### Adequate Response (1 point) - Sample D

> "The parameter T and the quadrature order both affect accuracy. When T is large, you need higher order, and when you want more accuracy, you also need higher order."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | YES - stated but not explained |
| T-order relationship | YES - stated but not explained |
| Module observations | NO - no specific observations cited |
| **Score** | **1 point** |

**Why this earns partial credit:** Identifies both relationships but lacks explanation of WHY and provides no evidence of module engagement.

---

### Inadequate Response (0 points) - Sample E

> "Quadrature order should always be as high as possible to get the best accuracy."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | INCORRECT - ignores efficiency tradeoff |
| T-order relationship | NO - not addressed |
| Module observations | NO |
| **Score** | **0 points** |

**Why this earns no credit:** Misses the concept of "optimal" order - there is a tradeoff between accuracy and computational cost. The response also ignores T entirely.

---

### Inadequate Response (0 points) - Sample F

> "T and quadrature order are independent. You choose order based on how fast you want the calculation to run."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Accuracy-order relationship | NO - does not connect to accuracy |
| T-order relationship | INCORRECT - claims they are independent |
| Module observations | NO |
| **Score** | **0 points** |

**Why this earns no credit:** Incorrect claim that T and optimal order are independent. The module clearly shows their relationship.

---

### Common Misconceptions for Q5

| Misconception | Correct Understanding |
|--------------|----------------------|
| "Higher order is always better" | Higher order has diminishing returns and higher cost |
| "T does not affect order choice" | Larger T generally requires higher order for same accuracy |
| "Order only affects speed" | Order primarily affects accuracy |
| "All T values need the same order" | Optimal order varies with T |

---

## Post-Check Item Q6: DIIS Effect on Convergence

### Question

**Describe in 2-3 sentences how DIIS changes the SCF convergence behavior compared to standard (non-accelerated) iteration. What specific differences did you observe?**

### Key Concepts Expected

1. DIIS reduces the number of iterations needed to converge
2. DIIS produces smoother convergence (avoids oscillations)
3. References to specific observations (iteration counts, energy curves)

---

### Exemplary Response (2 points) - Sample A

> "DIIS dramatically reduces the number of iterations required for SCF convergence. In the H2O preset, standard iteration took about 18 iterations to converge, while DIIS converged in only 8 iterations. The energy vs. iteration plot also showed that DIIS produces a smooth, monotonic decrease to the final energy, whereas standard iteration showed some oscillation before settling down."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | YES - specific numbers cited |
| Smoother convergence | YES - mentions monotonic decrease vs. oscillation |
| Specific observations | YES - cites H2O preset with iteration counts |
| **Score** | **2 points** |

**Why this earns full credit:** Addresses both key aspects (fewer iterations and smoother convergence) with specific numerical observations from the module.

---

### Exemplary Response (2 points) - Sample B

> "With DIIS enabled, convergence was much faster - the SCF calculation finished in roughly half the iterations compared to the non-accelerated version. I also noticed that without DIIS, the energy plot showed a 'zigzag' pattern in the first few iterations as if it was overshooting, but DIIS eliminated this oscillation and converged directly."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | YES - "roughly half" |
| Smoother convergence | YES - describes zigzag elimination |
| Specific observations | YES - describes energy plot behavior |
| **Score** | **2 points** |

**Why this earns full credit:** Correctly describes both effects with vivid description of the oscillation behavior observed.

---

### Adequate Response (1 point) - Sample C

> "DIIS makes the SCF calculation converge faster. It took fewer iterations with DIIS turned on."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | YES - stated correctly |
| Smoother convergence | NO - not mentioned |
| Specific observations | PARTIAL - no specific numbers |
| **Score** | **1 point** |

**Why this earns partial credit:** Correctly identifies that DIIS reduces iterations but does not address the qualitative change in convergence behavior (smoothness, oscillation reduction).

---

### Adequate Response (1 point) - Sample D

> "DIIS extrapolates from previous Fock matrices to speed up convergence. The convergence was smoother and there were no oscillations in the energy."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | NO - not explicitly stated |
| Smoother convergence | YES - mentions smoother, no oscillations |
| Specific observations | PARTIAL - describes pattern but no numbers |
| **Score** | **1 point** |

**Why this earns partial credit:** Correctly describes the smoothing effect but does not quantify the iteration reduction. The mention of "Fock matrices" is accurate but the question asks for observed effects, not mechanism.

---

### Inadequate Response (0 points) - Sample E

> "DIIS is a type of acceleration that makes the computer run faster."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | NO - confuses fewer iterations with faster computation |
| Smoother convergence | NO |
| Specific observations | NO |
| **Score** | **0 points** |

**Why this earns no credit:** Fundamental misconception - DIIS does not make the computer faster; it reduces the number of iterations needed. Each iteration takes the same time.

---

### Inadequate Response (0 points) - Sample F

> "I didn't notice much difference between DIIS on and off. The final energy was the same."

**Scoring Annotation:**

| Criterion | Assessment |
|-----------|------------|
| Fewer iterations | NO - claims no difference observed |
| Smoother convergence | NO |
| Specific observations | INCORRECT - should have observed clear differences |
| **Score** | **0 points** |

**Why this earns no credit:** The response indicates the student did not engage with the comparison activity or failed to observe the obvious differences. While correct that the final energy is the same, this misses the point of the question about convergence BEHAVIOR.

---

### Common Misconceptions for Q6

| Misconception | Correct Understanding |
|--------------|----------------------|
| "DIIS makes the computer faster" | DIIS reduces iteration count, not per-iteration speed |
| "DIIS changes the final answer" | DIIS reaches the same answer, just faster |
| "DIIS always helps" | DIIS can occasionally fail for difficult cases |
| "No visible difference" | Clear differences in iteration count and convergence pattern |

---

## Grading Tips for Inter-Rater Consistency

### Before Grading

1. **Read all sample responses** in this document first
2. **Identify the key concepts** for each item before scoring
3. **Agree on borderline cases** with co-graders if applicable

### During Grading

1. **Score holistically** - look for conceptual understanding, not keyword matching
2. **Read the entire response** before assigning a score
3. **Do not penalize** unusual phrasing if the concept is correct
4. **Award the higher score** when genuinely uncertain between two levels

### Borderline Decision Rules

| Situation | Recommendation |
|-----------|----------------|
| Correct concept but poorly explained | Award **1 point** |
| One concept correct, one incorrect | Award **1 point** |
| Both concepts present but very brief | Award **1-2 points** based on clarity |
| Correct but no module observations (Q5/Q6) | Award **1 point** maximum |
| Minor technical error but right idea | Award **1-2 points** based on severity |

### Common Grading Errors to Avoid

| Error | Correction |
|-------|------------|
| Requiring exact terminology | Accept synonyms and descriptions |
| Penalizing length | Judge content, not word count |
| Over-rewarding correct buzzwords | Ensure understanding, not just keywords |
| Inconsistent standards across students | Re-read samples periodically |

### Calibration Exercise

Before grading student work, have all graders independently score these sample responses, then compare. Discuss any discrepancies until agreement is reached.

---

## Summary Table: Sample Responses by Item and Score

| Item | Exemplary (2) | Adequate (1) | Inadequate (0) |
|------|---------------|--------------|----------------|
| **P5** | Explains convergence/stability varies with parameter range | States methods differ by range, lacks WHY | Claims methods are arbitrarily chosen or always-better |
| **P6** | Defines as successive iteration difference < threshold | Mentions tolerance but vague on comparison | Confuses with iteration count or exact correctness |
| **Q5** | Addresses both accuracy-order AND T-order relationships with observations | Addresses one relationship or lacks specificity | Claims order is independent of T or always maximal |
| **Q6** | Notes fewer iterations AND smoother convergence with specific observations | Notes one aspect (count OR smoothness) | Confuses computational speed or claims no difference |

---

## Appendix: Quick Reference Scoring Checklist

### P5 - Numerical Regimes

- [ ] Explains WHY different methods are used (convergence, stability, efficiency)
- [ ] Correctly characterizes when each method type is appropriate
- [ ] Does NOT suggest methods give different answers

### P6 - Convergence Tolerance

- [ ] Defines tolerance as difference between successive iterations
- [ ] Distinguishes from iteration count
- [ ] Does NOT claim tolerance guarantees exact correctness

### Q5 - Quadrature Order Selection

- [ ] Connects higher accuracy to higher order
- [ ] Addresses how T affects optimal order
- [ ] Includes specific observation from Rys module

### Q6 - DIIS Effect

- [ ] Notes reduction in iteration count
- [ ] Describes smoother/non-oscillatory convergence
- [ ] Includes specific observation from SCF module

---

*IQCP Lab Pack #1 v1.0 | Sample Responses Scoring Guide | https://iqcp.dev*
