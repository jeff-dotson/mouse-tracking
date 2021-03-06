
---
title: "Dynamically Assessing Respondent Quality in Conjoint Studies"
output: github_document
---
Dynamically Assessing Respondent Quality in Conjoint Studies
================

## Abstract

Although the advent of online consumer panels has reduced the cost of
data collection for survey-based marketing research projects, it has
introduced a whole host of additional problems that must be addressed in
order to make valid inference. Central to these concerns are issues
related to respondent quality. In this paper we introduce a new
technology, mouse tracking, that allows for passive surveillance of
mouse movements. Prior research has shown that these movements can be
linked to a variety of psychological states including engagement,
conflict, and affect.We show through a series of experiments that mouse
tracking can be used to assess respondent engagement.

### Keywords
mouse tracking, conjoint analysis, error scale modeling, Hierarchical Multinomial Logit 

### Disciplines
Behavioral Economics, Business, Marketing, Statistics 

## Introduction

Literature suggests that simple tracking and physiological measurements can be informative in modeling different behaviors. Some of these methods include survey task response times, heart rate tracking, galvanic skin response, and eye tracking. Something about how these so far are all easy to measure and collect but are decidedly noisy as a signal, and so we have settled on using mouse tracking.

The potential applications of this in conjoint include being able to determine engagement of respondents, how are choices being made, are respondents incorporating available data into their decisions, and can this be informative of their preferences.



## Empirical Methods and Models

In the Generalized Hierarchical Multinomial Logit (MNL) framework, the utility of respondent *i* for alternative *j* whithin choice c is the same as the standard Hierarchical Bayes model: 
$U_{ijc} = {X^\prime_{ijc}\beta_{i} + \epsilon_{ijc}}$ 

with $\epsilon_{ijc}$ distributed independent and identically distributed random variables (IID) in a Gumbel distribution.  We know that $\beta_{i}$ cannot be separately identified from the scale parameter of the errors ($\lambda$), so this is typically fixed to 1.

If we allow a Heterogenous Error scale, the generalized MNL model specifies a joint distribution for $\beta_{i}$ and $\lambda_{i}$, where $\beta_{i} \sim \mathcal{N}(\Delta^\prime{z_{i}}, V_{\beta})$ and $\log(\lambda_{i}) \sim \mathcal{N}(\delta^\prime{z_{i}}, \sigma^2)$.

**Not sure here what this is,find out what it is and introduce\expound on it**
Identification is acheived centering the covariates $(z_{i})$ across respondents. 

This model can be extended to incorporate task-level tracking data into the error term. This is done by a deterministic deviation in the error scale in task level differences,  
$\log(\lambda_{ic}) = \omega^{\prime}tr_{ic} + \delta^{\prime}z_{i} + \xi_{i}$. In this instance, identification is achieved by centering the task-level tracking variables, $tr_{ic}$ across tasks for each individual.


## Interpretation of results

## Conclusions and future research

## References


For general details on GitHub usage, project organization, and project
workflow, see [Research Assistant
Training](https://github.com/marcdotson/ra-training).

