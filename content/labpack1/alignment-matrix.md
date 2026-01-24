# Lab Pack #1: Assessment Alignment Matrix

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18

---

## Learning Outcomes Summary

The following six learning outcomes (LOs) are targeted by Lab Pack #1:

| ID | Learning Outcome |
|----|------------------|
| **LO1** | Describe the qualitative behavior of Boys functions F_m(T) as T varies, and explain why different computational regimes are necessary for numerical stability. |
| **LO2** | Explain the relationship between Rys quadrature order and integration accuracy, and select an appropriate order to meet a specified error target. |
| **LO3** | Interpret SCF energy convergence plots and explain how DIIS (Direct Inversion in the Iterative Subspace) accelerates convergence. |
| **LO4** | Connect numerical parameters (quadrature order, convergence thresholds, DIIS settings) to computational outcomes and cost. |
| **LO5** | Export run artifacts that document computational explorations for reproducibility. |
| **LO6** | Apply physical and mathematical reasoning to predict computational behavior before running calculations. |

---

## Concept Check Alignment Matrix

### Pre-Activity Items (P1-P6)

| Item | Content Focus | Primary LO | Secondary LO | Bloom's Level | Item Type |
|------|---------------|------------|--------------|---------------|-----------|
| **P1** | Quadrature points vs. accuracy/cost tradeoff | LO2 | LO4 | Understand | MC |
| **P2** | Convergence definition in iterative methods | LO3 | - | Remember | MC |
| **P3** | Gaussian integral value (exp(-x^2)) | LO1 | - | Remember | MC |
| **P4** | Hartree-Fock approximation nature | LO3 | - | Understand | MC |
| **P5** | Why different computational methods for different parameter ranges | LO1 | LO4 | Understand | SA |
| **P6** | What convergence to tolerance 10^-6 means | LO3 | LO4 | Understand | SA |

### Post-Activity Items (Q1-Q6)

| Item | Content Focus | Primary LO | Secondary LO | Bloom's Level | Item Type |
|------|---------------|------------|--------------|---------------|-----------|
| **Q1** | Rys quadrature order effect on accuracy | LO2 | - | Understand | MC |
| **Q2** | DIIS mechanism for SCF acceleration | LO3 | - | Understand | MC |
| **Q3** | Boys function large-T limit (asymptotic behavior) | LO1 | - | Understand | MC |
| **Q4** | Regime selection rationale (series vs. asymptotic) | LO1 | LO4 | Analyze | MC |
| **Q5** | Quadrature order depends on T and target accuracy | LO2 | LO4 | Apply | SA |
| **Q6** | DIIS effect on convergence behavior | LO3 | LO4 | Apply | SA |

---

## Detailed Item-to-LO Mapping

### Learning Outcome 1 (Boys Function Behavior and Regimes)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| P3 | Pre | Assesses foundational knowledge of Gaussian integrals, which underlie Boys function definition. |
| P5 | Pre | Directly probes understanding of why numerical algorithms switch methods at regime boundaries. |
| Q3 | Post | Tests knowledge of Boys function asymptotic behavior learned through module exploration. |
| Q4 | Post | Requires analysis of regime selection criteria, demonstrating understanding of convergence properties. |

**Coverage Assessment:** Strong coverage with 4 items (2 pre, 2 post) spanning Remember, Understand, and Analyze levels.

### Learning Outcome 2 (Rys Quadrature Order-Accuracy)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| P1 | Pre | Assesses general understanding of quadrature point-accuracy tradeoffs as prerequisite knowledge. |
| Q1 | Post | Tests specific understanding of Rys quadrature order effects gained from module exploration. |
| Q5 | Post | Requires application of order-accuracy relationship to justify quadrature selection. |

**Coverage Assessment:** Adequate coverage with 3 items (1 pre, 2 post) spanning Understand and Apply levels.

### Learning Outcome 3 (SCF Convergence and DIIS)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| P2 | Pre | Assesses foundational understanding of iterative convergence concept. |
| P4 | Pre | Assesses background knowledge of Hartree-Fock as an approximation method. |
| P6 | Pre | Probes understanding of convergence tolerance interpretation. |
| Q2 | Post | Tests knowledge of DIIS mechanism learned through exploration. |
| Q6 | Post | Requires application of DIIS understanding to describe observed convergence behavior. |

**Coverage Assessment:** Strong coverage with 5 items (3 pre, 2 post) spanning Remember, Understand, and Apply levels.

### Learning Outcome 4 (Parameter-Outcome Connections)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| P1 | Pre | Connects quadrature points to computational cost (secondary alignment). |
| P5 | Pre | Connects parameter ranges to method selection (secondary alignment). |
| P6 | Pre | Connects tolerance parameter to convergence outcome (secondary alignment). |
| Q4 | Post | Connects T value to regime selection (secondary alignment). |
| Q5 | Post | Connects both T and accuracy target to quadrature order choice (secondary alignment). |
| Q6 | Post | Connects DIIS settings to convergence behavior (secondary alignment). |

**Coverage Assessment:** LO4 is assessed as a secondary outcome across 6 items, reflecting its integrative nature.

### Learning Outcome 5 (Run Artifact Export)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| - | - | Not directly assessed by concept check items. |

**Coverage Assessment:** LO5 is a procedural skill assessed through the worksheet artifact export checkpoints, not the concept check. This is appropriate for a practical competency.

### Learning Outcome 6 (Predictive Reasoning)

| Item | Pre/Post | Alignment Rationale |
|------|----------|---------------------|
| - | - | Not directly assessed by concept check items. |

**Coverage Assessment:** LO6 is assessed through the synthesis questions (Q5.1-Q5.3) in the worksheet, which require prediction and reasoning. The concept check focuses on content knowledge that supports this reasoning skill.

---

## Cognitive Level Distribution

### Bloom's Taxonomy Classification

| Level | Definition | Pre Items | Post Items | Total |
|-------|------------|-----------|------------|-------|
| **Remember** | Recall facts and basic concepts | P2, P3 | - | 2 |
| **Understand** | Explain ideas or concepts | P1, P4, P5, P6 | Q1, Q2, Q3 | 7 |
| **Apply** | Use information in new situations | - | Q5, Q6 | 2 |
| **Analyze** | Draw connections among ideas | - | Q4 | 1 |
| **Evaluate** | Justify a decision or position | - | - | 0 |
| **Create** | Produce new or original work | - | - | 0 |

### Cognitive Progression Analysis

The pre-check items cluster at Remember (17%) and Understand (66%) levels, appropriate for assessing prior knowledge before the activity. The post-check items shift toward Understand (50%), Apply (33%), and Analyze (17%) levels, reflecting the deeper engagement expected after hands-on exploration.

```
Pre-Check:  [Remember ██] [Understand ████████████] [Apply] [Analyze]
Post-Check: [Remember] [Understand ██████] [Apply ████] [Analyze ██]
```

This progression is pedagogically appropriate: students begin with foundational recall and comprehension, then demonstrate application and analysis after experiential learning.

---

## Learning Outcome Coverage Summary

| LO | Pre Items | Post Items | Total Items | Cognitive Range |
|----|-----------|------------|-------------|-----------------|
| **LO1** | P3, P5 | Q3, Q4 | 4 | Remember - Analyze |
| **LO2** | P1 | Q1, Q5 | 3 | Understand - Apply |
| **LO3** | P2, P4, P6 | Q2, Q6 | 5 | Remember - Apply |
| **LO4** | (P1, P5, P6)* | (Q4, Q5, Q6)* | 6* | Understand - Apply |
| **LO5** | - | - | 0** | Procedural (worksheet) |
| **LO6** | - | - | 0** | Higher-order (synthesis) |

*LO4 is assessed as a secondary outcome
**LO5 and LO6 are assessed through worksheet activities rather than concept check

---

## Pre/Post Parallel Structure

The concept check uses parallel item pairs to enable direct measurement of learning gains:

| Pre Item | Post Item | Parallel Content Domain |
|----------|-----------|------------------------|
| P1 | Q1 | Quadrature points and accuracy |
| P2 | Q2 | Iterative method convergence mechanism |
| P3 | Q3 | Integral/function limiting behavior |
| P4 | Q4 | Method rationale and selection |
| P5 | Q5 | Parameter-method relationships (open response) |
| P6 | Q6 | Convergence behavior description (open response) |

This parallel structure supports normalized gain calculations:

```
g = (Post_score - Pre_score) / (Max_score - Pre_score)
```

Item-level parallel analysis enables identification of specific conceptual gains.

---

## Content Validity Argument

### Evidence of Content Representativeness

The concept check items comprehensively sample the content domain defined by the six learning outcomes. Each of the three computational topics (Boys functions, Rys quadrature, SCF/DIIS) receives dedicated coverage in both pre- and post-phases, ensuring that learning gains can be measured across all primary content areas. The item distribution (4 items for LO1, 3 for LO2, 5 for LO3) appropriately weights topics according to their emphasis in the worksheet activities, where Section 2 (Boys) receives 15-20 minutes, Section 3 (Rys) receives 15-20 minutes, and Section 4 (SCF) receives 20-30 minutes.

The cognitive levels employed span from Remember through Analyze, with the majority at the Understand and Apply levels. This distribution aligns with the educational goals of Lab Pack #1, which aims to build conceptual understanding through guided exploration rather than to develop expert-level evaluation or creation skills. The progression from lower-order pre-check items to higher-order post-check items mirrors the expected learning trajectory and provides evidence that gains reflect genuine understanding rather than memorization.

### Evidence of Item Quality

Each item was constructed following assessment best practices: multiple-choice items include plausible distractors based on documented student misconceptions, short-answer items specify expected response length, and all items can be answered within the 10-12 minute time allocation. The distractor analysis (provided in the scoring key) documents the reasoning behind incorrect options, demonstrating that alternatives were designed to diagnose specific misunderstandings rather than to trick students. Short-answer rubrics specify distinct criteria for full, partial, and no credit, enabling consistent scoring across graders.

### Limitations and Scope

Two learning outcomes (LO5: artifact export; LO6: predictive reasoning) are not directly assessed by concept check items. This is a deliberate design choice: LO5 represents a procedural skill best assessed through direct observation of artifact submission, while LO6 represents a higher-order integrative skill assessed through the worksheet synthesis questions. The concept check serves as a focused instrument for measuring content knowledge gains, complemented by other assessment components that address procedural and higher-order outcomes.

---

## Assessment Component Summary

| Component | Items | Purpose | LOs Assessed |
|-----------|-------|---------|--------------|
| **Pre-Check** | P1-P6 | Measure prior knowledge | LO1, LO2, LO3, LO4 |
| **Post-Check** | Q1-Q6 | Measure learning gains | LO1, LO2, LO3, LO4 |
| **Worksheet Questions** | Q2.1-Q5.3 | Guide exploration, assess understanding | LO1-LO6 |
| **Artifact Checkpoints** | 3 exports | Document reproducibility skills | LO5 |
| **Synthesis Questions** | Q5.1-Q5.3 | Assess predictive reasoning | LO6 |

---

*IQCP Lab Pack #1 v1.0 | Assessment Alignment Matrix | https://iqcp.dev*
