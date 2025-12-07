# trait_responsibility

Data, code, and materials for: Xiang, Y., Hsu, S., & Gershman, S. J. (2025). _Intuitions about trait responsibilities._

## Project Overview

Many popular depictions of mental life imagine personality traits as if they were agents in a small society—each contributing to a person's behavior, much like individuals contribute to a group's actions. Do people actually think about traits in this way? If so, then they should assign responsibility to these traits similar to how they would assign responsibility to individual people within a group.

To investigate this, we conducted two experiments. In **Experiment 1**, participants judged how trait scores predicted a person's actions. In **Experiment 2**, participants judged how responsible each trait was for those same actions. Using the predictive judgments from Experiment 1, we constructed and evaluated computational models of trait responsibility attribution, adapted from models of responsibility in groups.

Our results showed that responsibility judgments were best explained by a **production-style model** that assigns responsibility based on a trait's predictive contribution to behavior. Unlike studies of personal responsibility in groups, **counterfactual-style models** did not account for the data as well. These findings reveal both parallels and differences between how people attribute responsibility to traits within a person and individuals within a group.

## Links

Preprint: https://osf.io/preprints/psyarxiv/meyrx_v1  
Experiment 1: https://gershmanlab.com/experiments/yang/pc/experiment1.html  
Experiment 2: https://gershmanlab.com/experiments/yang/pc/experiment2_resp.html  
Preregistration: https://aspredicted.org/hh3sx5.pdf

---

## Code

Located in the `code` folder:

- `1_compute_trait_weights.R` — Analysis script for Experiment 1. Computes and visualizes trait weights.
- `2_responsibility_analysis.R` — Analysis script for Experiment 2. Simulates model predictions, compares models to each other, and evaluates them against behavioral data.
- `helper.R` — Helper functions for model simulation, visualization, and analysis.
- `output/` — Contains intermediate output files (e.g., saved data frames and model results) generated during the analyses. These files allow the scripts to skip time-consuming computations on subsequent runs.

---

## Data

Located in the `data` folder:

- `exp1.csv` — Data from Experiment 1.
  - `scenario`: Scenario index (see Supplement for the list of 20 scenarios).
  - `person`: Person index from 1–5 (each scenario shown five times with trait profiles randomly sampled).
  - `action_probability`: Participant-reported probability of taking the action.
  - `total_failed_attention_checks`: 0 or 1 (participants who failed 2 checks were asked to leave the study and their data were not saved).
  - Remaining five columns: Big Five trait profile for each person.

- `exp2.csv` — Data from Experiment 2.  
  Contains all columns from `exp1.csv` except `action_probability`, plus five `responsibility_*` columns representing participants' trait responsibility judgments.

- `population_data.csv` — Population Big Five data retrieved from:  
  https://github.com/automoto/big-five-data

---

## Experiment Materials

Located in the `experiments` folder:

- `exp1.html` and `exp2_resp.html` — Experiment 1 and Experiment 2 task scripts.
- `consent_exp1.html` and `consent_exp2.html` — Consent forms.
- `save_data.php` — Script for writing data to server.
- Experiments were built using jsPsych v8.2.2, available at:  
  https://github.com/jspsych/jsPsych/releases/tag/jspsych%408.2.2  
- Note: To run the experiments locally, the consent form needs to be commented out (i.e., comment out `timeline.push(consent);`) and jsPsych library files need to be added to this folder.
