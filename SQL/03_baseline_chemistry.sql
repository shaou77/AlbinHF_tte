-- v2 baseline chemistry extraction
-- Keeps the original -12h to +6h window and adds a pre-ICU/time-zero-only
-- window (-12h to 0h) for sensitivity analyses addressing overadjustment.

WITH target_window AS (
    SELECT
        c.stay_id,
        c.subject_id,
        c.icu_intime,
        c.icu_intime - INTERVAL '12 HOUR' AS baseline_start,
        c.icu_intime + INTERVAL '6 HOUR' AS baseline_end,
        c.icu_intime - INTERVAL '12 HOUR' AS preicu_start,
        c.icu_intime AS preicu_end
    FROM cohort_heart_failure_albumin_v2 c
),
chem_data AS (
    SELECT
        tw.stay_id,
        CASE
            WHEN le.charttime >= tw.baseline_start AND le.charttime <= tw.baseline_end THEN 1
            ELSE 0
        END AS in_baseline_window,
        CASE
            WHEN le.charttime >= tw.preicu_start AND le.charttime <= tw.preicu_end THEN 1
            ELSE 0
        END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (le.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - le.charttime)) AS time_before_icu,
        le.albumin,
        le.globulin,
        le.total_protein,
        le.aniongap,
        le.bicarbonate,
        le.bun,
        le.calcium,
        le.chloride,
        le.creatinine,
        le.glucose,
        le.sodium,
        le.potassium
    FROM target_window tw
    LEFT JOIN mimiciv_derived.chemistry le
      ON le.subject_id = tw.subject_id
     AND le.charttime >= tw.baseline_start
     AND le.charttime <= tw.baseline_end
)
SELECT
    stay_id,

    MIN(albumin) FILTER (WHERE in_baseline_window = 1) AS albumin_min,
    (ARRAY_AGG(albumin ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND albumin IS NOT NULL))[1] AS albumin_closest,
    MAX(creatinine) FILTER (WHERE in_baseline_window = 1) AS creatinine_max,
    (ARRAY_AGG(bun ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND bun IS NOT NULL))[1] AS bun_baseline,
    (ARRAY_AGG(glucose ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND glucose IS NOT NULL))[1] AS glucose_baseline,
    (ARRAY_AGG(sodium ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND sodium IS NOT NULL))[1] AS sodium_baseline,
    (ARRAY_AGG(potassium ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND potassium IS NOT NULL))[1] AS potassium_baseline,
    (ARRAY_AGG(chloride ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND chloride IS NOT NULL))[1] AS chloride_baseline,
    (ARRAY_AGG(bicarbonate ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND bicarbonate IS NOT NULL))[1] AS bicarbonate_baseline,
    (ARRAY_AGG(aniongap ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND aniongap IS NOT NULL))[1] AS aniongap_baseline,
    (ARRAY_AGG(calcium ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND calcium IS NOT NULL))[1] AS calcium_baseline,
    (ARRAY_AGG(total_protein ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND total_protein IS NOT NULL))[1] AS total_protein_baseline,
    (ARRAY_AGG(globulin ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND globulin IS NOT NULL))[1] AS globulin_baseline,

    MIN(albumin) FILTER (WHERE in_preicu_window = 1) AS albumin_min_preicu,
    (ARRAY_AGG(albumin ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND albumin IS NOT NULL))[1] AS albumin_closest_preicu,
    MAX(creatinine) FILTER (WHERE in_preicu_window = 1) AS creatinine_max_preicu,
    (ARRAY_AGG(bun ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND bun IS NOT NULL))[1] AS bun_preicu,
    (ARRAY_AGG(glucose ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND glucose IS NOT NULL))[1] AS glucose_preicu,
    (ARRAY_AGG(sodium ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND sodium IS NOT NULL))[1] AS sodium_preicu,
    (ARRAY_AGG(potassium ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND potassium IS NOT NULL))[1] AS potassium_preicu,
    (ARRAY_AGG(chloride ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND chloride IS NOT NULL))[1] AS chloride_preicu,
    (ARRAY_AGG(bicarbonate ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND bicarbonate IS NOT NULL))[1] AS bicarbonate_preicu,
    (ARRAY_AGG(aniongap ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND aniongap IS NOT NULL))[1] AS aniongap_preicu,
    (ARRAY_AGG(calcium ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND calcium IS NOT NULL))[1] AS calcium_preicu,
    (ARRAY_AGG(total_protein ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND total_protein IS NOT NULL))[1] AS total_protein_preicu,
    (ARRAY_AGG(globulin ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND globulin IS NOT NULL))[1] AS globulin_preicu
FROM chem_data
GROUP BY stay_id
ORDER BY stay_id;
