# Lab Pack #1: Assessment Instruments Master Document

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18
**Document Type:** Assessment Portfolio Master Document
**Target Publication:** J. Chem. Educ. Technology Report

---

## Executive Overview

This document consolidates all assessment instruments for Lab Pack #1 "From Boys to Orbitals" of the Interactive Quantum Chemistry Playground (IQCP). The assessment portfolio is designed to measure student learning across six learning outcomes spanning three computational chemistry modules: Boys functions, Rys quadrature, and SCF convergence with DIIS acceleration.

### Assessment Philosophy

The assessment framework follows established principles from science education research:

1. **Constructive alignment:** All assessment items trace directly to stated learning outcomes
2. **Multiple measures:** Conceptual understanding, procedural skills, and integrative reasoning are assessed through complementary instruments
3. **Formative and summative:** Pre/post assessments enable learning gain measurement while worksheet and performance tasks provide summative evaluation
4. **Authentic tasks:** Performance tasks require students to use IQCP to solve novel problems, mirroring authentic computational chemistry practice

### Total Point Allocation

| Instrument | Points | Percentage | Purpose |
|------------|--------|------------|---------|
| Concept Check (Pre) | 8 | 5.1% | Baseline knowledge assessment |
| Concept Check (Post) | 8 | 5.1% | Learning gain measurement |
| Worksheet | 100 | 63.3% | Guided exploration assessment |
| Performance Tasks | 42 | 26.5% | Applied competency assessment |
| **Total Portfolio** | **158** | **100%** | |

### Learning Outcomes Assessed

| ID | Learning Outcome | Primary Instruments |
|----|------------------|---------------------|
| **LO1** | Describe Boys function behavior and explain computational regime necessity | CC (P3, P5, Q3, Q4), WS (Q2.1-Q2.5), PT-Boys |
| **LO2** | Explain Rys quadrature order-accuracy relationship | CC (P1, Q1, Q5), WS (Q3.1-Q3.5), PT-Rys |
| **LO3** | Interpret SCF convergence and DIIS acceleration | CC (P2, P4, P6, Q2, Q6), WS (Q4.1-Q4.7), PT-SCF |
| **LO4** | Connect numerical parameters to computational outcomes | CC (P1, P5, P6, Q4-Q6), WS (all sections), PT (all) |
| **LO5** | Export reproducible run artifacts | WS (checkpoints), PT-Rys, PT-SCF |
| **LO6** | Apply predictive reasoning to computational behavior | WS (Q5.1-Q5.3), PT-Boys, PT-SCF |

---

## Table of Contents

### Component Documents

| Document | File | Purpose |
|----------|------|---------|
| 1. Alignment Matrix | [alignment-matrix.md](./alignment-matrix.md) | LO-item mapping and cognitive levels |
| 2. Concept Check | [concept-check.md](./concept-check.md) | 12-item pre/post assessment |
| 3. Concept Check Key | [concept-check-key.md](./concept-check-key.md) | Detailed scoring key with distractor analysis |
| 4. Student Worksheet | [worksheet-student.md](./worksheet-student.md) | 60-90 minute guided exploration |
| 5. Grading Rubric | [grading-rubric.md](./grading-rubric.md) | 100-point worksheet rubric |
| 6. Answer Key | [answer-key.md](./answer-key.md) | Expected worksheet responses |
| 7. Performance Rubrics | [performance-rubrics.md](./performance-rubrics.md) | PT-Boys, PT-Rys, PT-SCF rubrics |
| 8. Sample Responses | [sample-responses.md](./sample-responses.md) | Exemplar responses for calibration |

### Sections in This Document

1. [Executive Overview](#executive-overview)
2. [Learning Outcome Coverage Summary](#learning-outcome-coverage-summary)
3. [Assessment Sequence and Timeline](#assessment-sequence-and-timeline)
4. [Data Collection Protocols](#data-collection-protocols)
5. [Validity Argument](#validity-argument)
6. [Administration Timeline](#administration-timeline)
7. [Appendix: Summary Statistics Template](#appendix-summary-statistics-template)

---

## Learning Outcome Coverage Summary

### Coverage Matrix

The following matrix shows how each learning outcome is assessed across instruments, with cognitive level (Bloom's taxonomy) indicated.

| LO | Concept Check Pre | Concept Check Post | Worksheet | Performance Tasks |
|----|-------------------|-------------------|-----------|-------------------|
| **LO1** | P3 (R), P5 (U) | Q3 (U), Q4 (An) | Q2.1-Q2.5 (U-An) | PT-Boys (Ap-An) |
| **LO2** | P1 (U) | Q1 (U), Q5 (Ap) | Q3.1-Q3.5 (U-Ap) | PT-Rys (Ap-An) |
| **LO3** | P2 (R), P4 (U), P6 (U) | Q2 (U), Q6 (Ap) | Q4.1-Q4.7 (U-Ev) | PT-SCF (Ap-An) |
| **LO4** | P1, P5, P6 (secondary) | Q4, Q5, Q6 (secondary) | All sections | All tasks |
| **LO5** | - | - | Checkpoints (Ap) | PT-Rys, PT-SCF (Ap) |
| **LO6** | - | - | Q5.1-Q5.3 (An-Ev) | PT-Boys, PT-SCF (An) |

**Key:** R = Remember, U = Understand, Ap = Apply, An = Analyze, Ev = Evaluate

### Cognitive Level Distribution

| Level | Pre-Check | Post-Check | Worksheet | Performance Tasks | Total |
|-------|-----------|------------|-----------|-------------------|-------|
| Remember | 2 (33%) | 0 (0%) | 3 (15%) | 0 (0%) | 5 (13%) |
| Understand | 4 (67%) | 3 (50%) | 8 (40%) | 0 (0%) | 15 (39%) |
| Apply | 0 (0%) | 2 (33%) | 5 (25%) | 6 (67%) | 13 (34%) |
| Analyze | 0 (0%) | 1 (17%) | 3 (15%) | 3 (33%) | 7 (18%) |
| Evaluate | 0 (0%) | 0 (0%) | 1 (5%) | 0 (0%) | 1 (3%) |

The cognitive progression from pre-check (lower-order) through post-check and worksheet (middle-order) to performance tasks (higher-order) reflects the expected learning trajectory.

---

## Assessment Sequence and Timeline

### Recommended Implementation Order

```
PRE-CHECK (10-12 min)
    |
    v
WORKSHEET (60-90 min)
    |-- Section 1: Warm-up (~5 min)
    |-- Section 2: Boys Functions (~15-20 min)
    |-- Section 3: Rys Quadrature (~15-20 min)
    |-- Section 4: SCF Convergence (~20-30 min)
    |-- Section 5: Synthesis (~5-10 min)
    |
    v
POST-CHECK (10-12 min)
    |
    v
PERFORMANCE TASKS (30-45 min) [Optional/Extension]
    |-- PT-Boys (~10-15 min)
    |-- PT-Rys (~10-15 min)
    |-- PT-SCF (~10-15 min)
```

### When to Use Each Instrument

| Instrument | When to Administer | Purpose |
|------------|-------------------|---------|
| **Pre-Check** | Start of session, before any IQCP use | Establish baseline knowledge |
| **Worksheet** | Main activity period | Guide exploration, formative feedback |
| **Post-Check** | Immediately after worksheet completion | Measure learning gains |
| **Performance Tasks** | After post-check OR separate session | Assess transfer and application |

### Timing Recommendations

| Session Format | Total Time | Allocation |
|----------------|------------|------------|
| **Standard (90 min)** | 90 min | Pre (10) + Worksheet (60) + Post (10) + Buffer (10) |
| **Extended (120 min)** | 120 min | Pre (10) + Worksheet (70) + Post (10) + PT (30) |
| **Two-Session** | 90 + 45 min | Session 1: Pre + Worksheet + Post; Session 2: Performance Tasks |
| **Compressed (75 min)** | 75 min | Pre (8) + Worksheet (55) + Post (8) + Buffer (4) |

---

## Data Collection Protocols

### For Course Assessment

Standard grading uses the following components and weights:

| Component | Raw Points | Suggested Weight |
|-----------|------------|------------------|
| Worksheet | 100 | 60-70% |
| Concept Checks | 16 | 10-15% |
| Performance Tasks | 42 | 20-30% |

**Recommended grade calculation:**
```
Course Grade = (0.65 × Worksheet%) + (0.10 × CC%) + (0.25 × PT%)
```

### For Publication Data Collection

When collecting data for J. Chem. Educ. publication, follow these protocols:

#### Required Documentation

1. **IRB/Ethics Approval**
   - Obtain institutional review board approval before data collection
   - Use consent forms that allow publication of aggregate data
   - Ensure student identifiers are removed or pseudonymized

2. **Pre/Post Administration**
   - Administer pre-check at least 5 minutes before students access IQCP
   - Administer post-check immediately after worksheet completion
   - Do not allow IQCP access during concept checks
   - Record administration date, time, and any irregularities

3. **Data Recording**
   - Use standardized scoring sheets
   - Score short-answer items using provided rubrics
   - Record raw scores (not percentages) for statistical analysis
   - Note any items with inter-rater disagreement

#### Recommended Sample Sizes

| Analysis Type | Minimum N | Recommended N |
|---------------|-----------|---------------|
| Descriptive statistics | 15 | 30+ |
| Pre/post paired t-test | 20 | 40+ |
| Normalized gain analysis | 25 | 50+ |
| Item-level analysis | 30 | 75+ |

#### Data Quality Checks

Before analysis, verify:

- [ ] All students completed both pre and post checks
- [ ] No more than 10% missing data per item
- [ ] Pre-check administered before IQCP exposure
- [ ] Post-check administered same day as activity
- [ ] Inter-rater reliability > 0.7 for short-answer items

---

## Validity Argument

This section presents the validity argument for the Lab Pack #1 assessment portfolio, structured according to Kane's (2006) argument-based approach to validation.

### Content Validity Evidence

#### 1. Domain Representativeness

The assessment portfolio comprehensively samples the content domain defined by the six learning outcomes:

**Coverage Analysis:**
- LO1 (Boys functions): 4 concept check items + 5 worksheet questions + 1 performance task
- LO2 (Rys quadrature): 3 concept check items + 5 worksheet questions + 1 performance task
- LO3 (SCF/DIIS): 5 concept check items + 7 worksheet questions + 1 performance task
- LO4 (Parameter-outcome): Secondary coverage across 6+ items per instrument
- LO5 (Artifact export): 3 worksheet checkpoints + 2 performance task components
- LO6 (Predictive reasoning): 3 synthesis questions + 2 performance task components

**Content Expert Review:**
All items were developed by a computational chemistry domain expert and reviewed for accuracy and pedagogical appropriateness. Each item maps to specific theoretical content from established references (Shavitt, 1963; Dupuis et al., 1976; Pulay, 1980, 1982).

#### 2. Cognitive Level Alignment

Items span Bloom's taxonomy levels appropriate to learning objectives:
- Pre-check items emphasize Remember (17%) and Understand (67%) as appropriate for baseline assessment
- Post-check items shift toward Apply (33%) and Analyze (17%)
- Performance tasks require primarily Apply and Analyze skills
- This progression matches the expected learning trajectory from novice to developing competence

#### 3. Item Quality

Multiple-choice items include:
- One unambiguously correct answer
- Plausible distractors based on documented misconceptions
- Clear, unambiguous stems
- Appropriate reading level for upper-undergraduate students

Short-answer items include:
- Clear prompts specifying expected response length
- Detailed rubrics with criteria for full, partial, and no credit
- Sample responses at each score level for calibration

### Construct Validity Evidence

#### 1. Theoretical Foundation

The assessment constructs align with established theories of computational chemistry understanding:
- Boys function regimes reflect numerical analysis principles (convergence, stability)
- Rys quadrature items assess understanding of Gaussian quadrature theory
- SCF items probe understanding of iterative methods and variational principles
- DIIS items assess understanding of acceleration techniques

#### 2. Internal Structure

Pre/post parallel structure enables valid gain calculation:
- Each pre-check item has a corresponding post-check item in the same content domain
- Parallel items assess the same LO but may differ in specific context
- This design supports normalized gain calculation: g = (post - pre) / (max - pre)

#### 3. Convergent Evidence

Multiple assessment methods converge on the same constructs:
- Concept check MC items, short-answer items, and worksheet questions address overlapping content
- Performance task scores should correlate with corresponding worksheet section scores
- Artifact export tasks provide procedural verification of declarative knowledge

### Consequential Validity Considerations

#### 1. Intended Consequences

- Students gain conceptual understanding of computational chemistry foundations
- Instructors receive actionable feedback about student learning
- Assessment data support evidence-based improvement of IQCP
- Publication contributes to chemistry education literature

#### 2. Unintended Consequence Mitigation

- Items assess understanding, not IQCP navigation skill
- Multiple assessment modes accommodate different learning styles
- Partial credit options reward developing competence
- Clear rubrics ensure consistent, fair scoring

### Limitations and Scope

1. **LO5 and LO6 Coverage:** These higher-order outcomes receive less direct assessment in the concept check. This is by design: LO5 (artifact export) is procedural and assessed through checkpoints; LO6 (predictive reasoning) is assessed through synthesis questions and performance tasks.

2. **Sample Specificity:** Performance tasks use specific IQCP parameter combinations. Results generalize to understanding of the underlying concepts, not to specific parameter values.

3. **Timing Constraints:** The compressed format (75 min) may not allow full completion, potentially affecting validity of incomplete assessments.

---

## Administration Timeline

### Session Preparation (1-2 Days Before)

| Task | Time | Notes |
|------|------|-------|
| Print pre/post concept checks | 15 min | One per student, separate sheets |
| Print worksheets | 30 min | Can be digital if preferred |
| Test IQCP deep links | 15 min | Verify all links load correctly |
| Review answer key | 20 min | Familiarize with expected responses |
| Prepare submission method | 10 min | LMS, email, or paper collection |

### Session Day

| Time | Activity | Materials |
|------|----------|-----------|
| T+0 | Welcome, distribute pre-check | Pre-check forms |
| T+2 | Pre-check administration | Timer visible |
| T+12 | Collect pre-checks, distribute worksheets | Worksheets |
| T+15 | Students begin IQCP exploration | Computers/laptops |
| T+30 | Circulate, answer procedural questions | - |
| T+45 | Mid-point check-in (optional) | - |
| T+60 | 15-minute warning | - |
| T+70 | Begin collecting worksheets | Post-check forms |
| T+75 | Distribute post-check | Timer visible |
| T+85 | Collect post-checks | - |
| T+90 | Session complete | - |

### Post-Session (Within 1 Week)

| Task | Time | Notes |
|------|------|-------|
| Score pre/post concept checks | 1-2 hr | Use scoring key |
| Score worksheets | 3-4 hr | Use grading rubric |
| Enter scores to gradebook | 30 min | - |
| Calculate class statistics | 30 min | Use template in Appendix |
| Identify common difficulties | 30 min | For follow-up instruction |

### Performance Tasks (Optional Second Session)

| Time | Activity |
|------|----------|
| T+0 | Distribute PT prompts |
| T+5 | Students begin PT-Boys |
| T+20 | Transition to PT-Rys |
| T+35 | Transition to PT-SCF |
| T+50 | Collect all materials |

---

## Appendix: Summary Statistics Template

Use this template to record and analyze assessment data for publication or course improvement.

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
| **Total** | | | | | |

#### Learning Gain Analysis

| Metric | Value |
|--------|-------|
| Pre-check mean (out of 8) | |
| Post-check mean (out of 8) | |
| Raw gain (Post - Pre) | |
| Normalized gain g | |
| Effect size (Cohen's d) | |
| Paired t-test p-value | |

**Normalized Gain Formula:**
```
g = (Post_mean - Pre_mean) / (8 - Pre_mean)
```

**Interpretation Guide:**
- g < 0.3: Low gain
- 0.3 <= g < 0.7: Medium gain
- g >= 0.7: High gain

### B. Worksheet Statistics

| Section | N | Mean | SD | Min | Max |
|---------|---|------|----|----|-----|
| Section 2: Boys (25 pts) | | | | | |
| Section 3: Rys (25 pts) | | | | | |
| Section 4: SCF (35 pts) | | | | | |
| Section 5: Synthesis (15 pts) | | | | | |
| **Total (100 pts)** | | | | | |

#### Item-Level Analysis (Record for items with < 60% or > 90% success)

| Question | Mean | Issue Identified | Action Taken |
|----------|------|------------------|--------------|
| | | | |
| | | | |
| | | | |

### C. Performance Task Statistics

| Task | N | Mean | SD | Min | Max |
|------|---|------|----|----|-----|
| PT-Boys (14 pts) | | | | | |
| PT-Rys (14 pts) | | | | | |
| PT-SCF (14 pts) | | | | | |
| **Total (42 pts)** | | | | | |

#### Dimension-Level Analysis

| Dimension | Mean | Notes |
|-----------|------|-------|
| Prediction/Selection | | |
| Evidence Quality | | |
| Reasoning/Explanation | | |
| Artifact Quality | | |

### D. Correlation Analysis

| Comparison | r | p | Interpretation |
|------------|---|---|----------------|
| Pre vs. Post | | | Expected positive |
| Worksheet vs. PT | | | Convergent validity |
| CC Post vs. Worksheet | | | Convergent validity |
| Boys section vs. PT-Boys | | | Within-topic consistency |
| Rys section vs. PT-Rys | | | Within-topic consistency |
| SCF section vs. PT-SCF | | | Within-topic consistency |

### E. Reliability Analysis

| Instrument | Cronbach's alpha | KR-20 | Notes |
|------------|------------------|-------|-------|
| Pre-check MC (4 items) | | | Target > 0.6 |
| Post-check MC (4 items) | | | Target > 0.6 |
| Combined MC (8 items) | | | Target > 0.7 |
| Worksheet (20 items) | | N/A | Target > 0.8 |

### F. Inter-Rater Reliability (Short-Answer Items)

| Item | Cohen's kappa | % Exact Agreement | Notes |
|------|---------------|-------------------|-------|
| P5 | | | Target > 0.7 |
| P6 | | | Target > 0.7 |
| Q5 | | | Target > 0.7 |
| Q6 | | | Target > 0.7 |
| Worksheet SA items | | | Average across items |

### G. Student Feedback Summary (Optional)

| Question | Strongly Disagree | Disagree | Neutral | Agree | Strongly Agree |
|----------|-------------------|----------|---------|-------|----------------|
| The lab helped me understand Boys functions | | | | | |
| The lab helped me understand Rys quadrature | | | | | |
| The lab helped me understand SCF convergence | | | | | |
| IQCP was easy to use | | | | | |
| The worksheet instructions were clear | | | | | |
| I would recommend this lab to other students | | | | | |

### H. Instructor Notes and Observations

**What worked well:**

1.
2.
3.

**What could be improved:**

1.
2.
3.

**Specific items needing revision:**

| Item | Issue | Proposed Change |
|------|-------|-----------------|
| | | |
| | | |

**Technical issues encountered:**

| Issue | Frequency | Resolution |
|-------|-----------|------------|
| | | |
| | | |

---

## References

Dupuis, M., Rys, J., and King, H.F. (1976). Evaluation of molecular integrals over Gaussian basis functions. *Journal of Chemical Physics*, 65(1), 111-116.

Hake, R.R. (1998). Interactive-engagement versus traditional methods: A six-thousand-student survey of mechanics test data for introductory physics courses. *American Journal of Physics*, 66(1), 64-74.

Kane, M.T. (2006). Validation. In R.L. Brennan (Ed.), *Educational Measurement* (4th ed., pp. 17-64). American Council on Education.

Pulay, P. (1980). Convergence acceleration of iterative sequences: The case of SCF iteration. *Chemical Physics Letters*, 73(2), 393-398.

Pulay, P. (1982). Improved SCF convergence acceleration. *Journal of Computational Chemistry*, 3(4), 556-560.

Shavitt, I. (1963). The Gaussian Function in Calculations of Statistical Mechanics and Quantum Mechanics. *Methods in Computational Physics*, Vol. 2, pp. 1-45.

Wiggins, G., and McTighe, J. (2005). *Understanding by Design* (2nd ed.). Association for Supervision and Curriculum Development.

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-18 | IQCP Team | Initial release |

---

*Lab Pack #1 Assessment Instruments Master Document v1.0*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
