# Aperiodic Brain Activity, Interoceptive Self-Regulation, Psychological Distress, and Perceived Stress

Code and data to reproduce the analyses reported in the paper: estimation of the aperiodic (1/f) exponent of resting-state EEG with FOOOF, and a structural equation mediation model linking interoceptive self-regulation (MAIA), perceived stress (PSS), psychological distress (GHQ-12), and the aperiodic exponent.

## Repository structure

```
Data/
  Behavioral_data.csv          Behavioral scores (record_id, PSS_DIRt, GHQ_DIRt, MAIA_Autorregulacion_DIRd)
  Prepro_data/Chile/CN_02/     Preprocessed resting-state EEG (EEGLAB .set, 32 channels, 128 Hz), one folder per subject
fooof/
  processing.py                PSD (Welch) + FOOOF per channel -> global exponent per subject -> merge with behavioral data
Analisis/
  df.csv                       Analysis dataset (output of processing.py)
  Fixed_model.R                Mediation model (lavaan, 5000 BCa bootstrap), results tables, correlations, and descriptives
  power_MonteCarlo.R           A priori Monte Carlo power analysis (simulated data only; does not use the real data)
```

## Reproducing the analyses

Run everything from the repository root.

1. Aperiodic exponent (Python 3.10+):

   ```bash
   pip install -r requirements.txt
   python fooof/processing.py
   ```

   This writes `Analisis/df.csv` and, as an intermediate output, `Results_Aperiodic_2026/Chile/CN_02/CN_02_Exponent_Global.csv`.

2. Mediation model (R, working directory = repository root):

   ```r
   source("Analisis/Fixed_model.R")
   ```

   Required packages: `lavaan`, `tidyverse`, `flextable`, `officer`, `pagedown` (needs Chrome/Chromium for PDF export), `Hmisc`, `corrplot`, `psych`. Tables are saved to `Analisis/Tablas/`.

3. Power analysis (optional; about 1.5 h with `NSIM = 5000`):

   ```r
   source("Analisis/power_MonteCarlo.R")
   ```

   Required packages: `lavaan`, `MASS`, `dplyr`, `ggplot2`, `tidyr`. Results are saved to `Analisis/MonteCarlo/`.

## Variables

| Variable | Description |
|---|---|
| `record_id` | Pseudonymized participant ID |
| `Exponent__Global` | Aperiodic exponent (FOOOF, fixed mode, 2-40 Hz), averaged across the 32 channels |
| `PSS_DIRt` | Perceived Stress Scale, total score |
| `GHQ_DIRt` | General Health Questionnaire (GHQ-12), total score |
| `MAIA_Autorregulacion_DIRd` | MAIA Self-Regulation subscale score |

## Data notes

- 94 EEG recordings and 103 participants with behavioral data; the 92 participants with both form the analytic sample.
- Data are pseudonymized: the `Subject_IDxxx` identifiers cannot be used to identify participants. Recording dates and times, original file paths, and processing timestamps were removed from the `.set` files. The EEG signal itself was not modified.
