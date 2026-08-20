-- v2 baseline arterial blood gas extraction
-- Keeps original -12h to +6h fields and adds -12h to 0h pre-ICU fields.

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
bg_data AS (
    SELECT
        tw.stay_id,
        CASE
            WHEN bg.charttime >= tw.baseline_start AND bg.charttime <= tw.baseline_end THEN 1
            ELSE 0
        END AS in_baseline_window,
        CASE
            WHEN bg.charttime >= tw.preicu_start AND bg.charttime <= tw.preicu_end THEN 1
            ELSE 0
        END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (bg.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - bg.charttime)) AS time_before_icu,
        bg.ph,
        bg.po2,
        bg.pco2,
        bg.lactate,
        bg.baseexcess,
        bg.bicarbonate,
        bg.totalco2,
        bg.pao2fio2ratio
    FROM target_window tw
    LEFT JOIN mimiciv_derived.bg bg
      ON bg.subject_id = tw.subject_id
     AND bg.charttime >= tw.baseline_start
     AND bg.charttime <= tw.baseline_end
     AND bg.specimen = 'ART.'
)
SELECT
    stay_id,
    (ARRAY_AGG(ph ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND ph IS NOT NULL))[1] AS ph_baseline,
    (ARRAY_AGG(po2 ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND po2 IS NOT NULL))[1] AS po2_baseline,
    (ARRAY_AGG(pco2 ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND pco2 IS NOT NULL))[1] AS pco2_baseline,
    (ARRAY_AGG(lactate ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND lactate IS NOT NULL))[1] AS lactate_bg_baseline,
    (ARRAY_AGG(baseexcess ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND baseexcess IS NOT NULL))[1] AS baseexcess_baseline,
    (ARRAY_AGG(bicarbonate ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND bicarbonate IS NOT NULL))[1] AS bicarbonate_bg_baseline,
    (ARRAY_AGG(totalco2 ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND totalco2 IS NOT NULL))[1] AS totalco2_baseline,
    (ARRAY_AGG(pao2fio2ratio ORDER BY time_diff_abs ASC) FILTER (WHERE in_baseline_window = 1 AND pao2fio2ratio IS NOT NULL))[1] AS pao2fio2_baseline,

    (ARRAY_AGG(ph ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND ph IS NOT NULL))[1] AS ph_preicu,
    (ARRAY_AGG(po2 ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND po2 IS NOT NULL))[1] AS po2_preicu,
    (ARRAY_AGG(pco2 ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND pco2 IS NOT NULL))[1] AS pco2_preicu,
    (ARRAY_AGG(lactate ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND lactate IS NOT NULL))[1] AS lactate_bg_preicu,
    (ARRAY_AGG(baseexcess ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND baseexcess IS NOT NULL))[1] AS baseexcess_preicu,
    (ARRAY_AGG(bicarbonate ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND bicarbonate IS NOT NULL))[1] AS bicarbonate_bg_preicu,
    (ARRAY_AGG(totalco2 ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND totalco2 IS NOT NULL))[1] AS totalco2_preicu,
    (ARRAY_AGG(pao2fio2ratio ORDER BY time_before_icu ASC) FILTER (WHERE in_preicu_window = 1 AND pao2fio2ratio IS NOT NULL))[1] AS pao2fio2_preicu
FROM bg_data
GROUP BY stay_id
ORDER BY stay_id;
