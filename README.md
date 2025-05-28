# Functional Connectivity Analysis in EEG using Coherence (MATLAB)

This project analyzes functional connectivity in EEG signals collected under three distinct states: active mental calculation with eyes closed, resting with eyes closed, and resting with eyes open.

## Objective

To evaluate differences in EEG channel connectivity through coherence analysis, with a focus on specific frequency bands (theta, alpha, beta), and to represent the results using coherence matrices to identify patterns indicative of the three functional states.

## Methodology

- Preprocessing with offset removal and filtering (notch, high-pass, low-pass)
- Band-pass filtering (theta: 4–8 Hz, alpha: 8–13 Hz, etc.)
- Coherence calculation between all pairs of EEG channels using `mscohere`
- Generation of heatmaps to visualize spatial connectivity patterns

## Technologies Used

- MATLAB
- EEGLAB (for `topoplot` and electrode location loading)
- Real EEG signals (.bdf format)

## Contents

- `main.m`: main script that executes the entire processing pipeline
- `apply_filters.m`: function to apply multiple filters to EEG signals
- `calculate_coherence.m`: function that computes the coherence matrix between channels
- `projeto_final_v3.pdf`: final project report (written in Portuguese)

## Sample Results

Includes:
- Coherence matrices for each condition
- Topographic maps of average coherence per channel
- Coherence differences between conditions (e.g., active vs. resting)

## Author

João Mota — Biomedical Engineering student at NOVA FCT (Lisbon, Portugal)

## Context

This project was developed as part of the **Electrophysiology** course at NOVA School of Science and Technology (FCT).
