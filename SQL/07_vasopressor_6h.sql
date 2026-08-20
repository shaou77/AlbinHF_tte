-- v2 vasopressor extraction
-- Keeps old-compatible fields and adds cumulative / any-use fields in 0-6h.

WITH target_window AS (
    SELECT
        c.stay_id,
        c.icu_intime,
        c.icu_intime + INTERVAL '6 HOUR' AS window_end
    FROM cohort_heart_failure_albumin_v2 c
),
vaso_clipped AS (
    SELECT
        tw.stay_id,
        neq.norepinephrine_equivalent_dose AS rate_mcg_kg_min,
        GREATEST(neq.starttime, tw.icu_intime) AS effective_start,
        LEAST(neq.endtime, tw.window_end) AS effective_end
    FROM target_window tw
    JOIN mimiciv_derived.norepinephrine_equivalent_dose neq
      ON tw.stay_id = neq.stay_id
    WHERE neq.starttime <= tw.window_end
      AND neq.endtime > tw.icu_intime
      AND neq.norepinephrine_equivalent_dose > 0
),
vaso_agg AS (
    SELECT
        stay_id,
        MAX(rate_mcg_kg_min) AS max_neq_rate_6h,
        SUM(
            rate_mcg_kg_min * (EXTRACT(EPOCH FROM (effective_end - effective_start)) / 60.0)
        ) AS total_neq_mcg_kg_6h,
        MIN(effective_start) AS first_vaso_time_6h
    FROM vaso_clipped
    GROUP BY stay_id
),
vaso_at_admission AS (
    SELECT
        tw.stay_id,
        CASE
            WHEN MIN(neq.starttime) < tw.icu_intime + INTERVAL '1 HOUR' THEN 1
            ELSE 0
        END AS vaso_use_at_admission
    FROM target_window tw
    JOIN mimiciv_derived.norepinephrine_equivalent_dose neq
      ON tw.stay_id = neq.stay_id
    WHERE neq.starttime <= tw.window_end
      AND neq.endtime > tw.icu_intime
      AND neq.norepinephrine_equivalent_dose > 0
    GROUP BY tw.stay_id, tw.icu_intime
)
SELECT
    c.stay_id,
    ROUND(COALESCE(v.max_neq_rate_6h, 0)::numeric, 4) AS vaso_ne_dose_max_6h,
    COALESCE(va.vaso_use_at_admission, 0) AS vaso_use_at_admission,
    ROUND(COALESCE(v.max_neq_rate_6h, 0)::numeric, 4) AS vaso_ne_rate_max_6h,
    ROUND(COALESCE(v.total_neq_mcg_kg_6h, 0)::numeric, 1) AS vaso_ne_cumulative_6h,
    CASE WHEN COALESCE(v.max_neq_rate_6h, 0) > 0 THEN 1 ELSE 0 END AS vaso_active_6h_flag,
    v.first_vaso_time_6h
FROM cohort_heart_failure_albumin_v2 c
LEFT JOIN vaso_agg v
  ON c.stay_id = v.stay_id
LEFT JOIN vaso_at_admission va
  ON c.stay_id = va.stay_id
ORDER BY c.stay_id;
