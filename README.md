# Early Albumin Administration and Mortality in ICU Patients With Heart Failure and Hypoalbuminemia: A Target Trial Emulation

## About This Project

This repository contains the SQL and R code for the published open-access article **"Early albumin administration and mortality in ICU patients with heart failure and hypoalbuminemia: a target trial emulation"**. The study used **MIMIC-IV version 3.1**.

## Published Article

- **Journal:** *BMC Cardiovascular Disorders*
- **Published:** 26 August 2026
- **Article page:** [Springer Nature Link](https://link.springer.com/article/10.1186/s12872-026-06514-0)
- **DOI:** [10.1186/s12872-026-06514-0](https://doi.org/10.1186/s12872-026-06514-0)
- **Open access:** [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)

**Full citation:** Yingfang She, Xu Zheng, Hongwu Guo, Wendi Xiang, Liang Luo, and Yide Li. "Early albumin administration and mortality in ICU patients with heart failure and hypoalbuminemia: a target trial emulation." *BMC Cardiovascular Disorders* (2026). https://doi.org/10.1186/s12872-026-06514-0

## Abstract

This study emulated a parallel-group target trial using MIMIC-IV version 3.1 to evaluate albumin initiation from 6 to less than 48 hours after ICU admission in adults with heart failure and hypoalbuminemia. Among 1,308 eligible patients, 80 initiated albumin and 1,228 did not. Weighted 28-day cumulative mortality was 32.9% in the albumin group and 30.7% in the no-albumin group (risk difference, 2.2%; 95% CI, -7.6% to 12.0%; hazard ratio, 1.037; 95% CI, 0.890 to 1.207). Albumin initiation during this window did not demonstrate lower 28-day mortality, although the small treated group and wide confidence intervals do not exclude clinically meaningful benefit or harm.

## Data Availability

The SQL scripts used to construct the study cohort and extract the analysis variables are available in the [`SQL`](./SQL) directory. Access to the underlying MIMIC-IV version 3.1 data is governed by the PhysioNet credentialing and data use requirements.

## Running the Analysis

1. Run the scripts in the [`SQL`](./SQL) directory in numerical order, from `00_first_hf_lab.sql` through `11_albumin_events_48h.sql`, to construct the study cohort and extract the required analysis variables from MIMIC-IV version 3.1.
2. After completing the SQL data preparation, run the analysis scripts in the [`R`](./R) directory. Start with `landmark_analysis_v2.R`; the secondary and figure-generation scripts use data or saved results produced by the main analysis.

The SQL data-preparation stage must be completed before running the R analysis.

## Authors

- Yingfang She
- Xu Zheng
- Hongwu Guo
- Wendi Xiang
- Liang Luo
- Yide Li

Yingfang She, Xu Zheng, and Hongwu Guo contributed equally and share first authorship. Liang Luo and Yide Li are the corresponding authors.

## License

This project is for academic research purposes only.

---

*This repository accompanies the published article "Early albumin administration and mortality in ICU patients with heart failure and hypoalbuminemia: a target trial emulation," based on MIMIC-IV version 3.1.*
