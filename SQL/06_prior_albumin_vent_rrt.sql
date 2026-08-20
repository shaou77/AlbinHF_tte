-- v2 prior albumin, mechanical ventilation, and RRT extraction
-- Keeps original -12h to +6h baseline flags and adds pre-ICU-only flags.

WITH target_window AS (
    SELECT
        c.stay_id,
        c.subject_id,
        c.hadm_id,
        c.icu_intime,
        c.icu_intime - INTERVAL '12 HOUR' AS baseline_start,
        c.icu_intime + INTERVAL '6 HOUR' AS baseline_end,
        c.icu_intime - INTERVAL '12 HOUR' AS preicu_start,
        c.icu_intime AS preicu_end
    FROM cohort_heart_failure_albumin_v2 c
),
prior_alb AS (
    SELECT
        tw.stay_id,
        1 AS has_prior_albumin
    FROM target_window tw
    JOIN mimiciv_hosp.prescriptions p
      ON tw.hadm_id = p.hadm_id
    WHERE (p.drug ILIKE '%albumin%' OR p.drug_type ILIKE '%albumin%')
      AND p.starttime < tw.icu_intime
    GROUP BY tw.stay_id
),
baseline_vent AS (
    SELECT
        tw.stay_id,
        1 AS is_vent_baseline
    FROM target_window tw
    JOIN mimiciv_derived.ventilation v
      ON v.stay_id = tw.stay_id
    WHERE v.ventilation_status IN ('InvasiveVent', 'Tracheostomy')
      AND v.starttime <= tw.baseline_end
      AND v.endtime >= tw.baseline_start
    GROUP BY tw.stay_id
),
preicu_vent AS (
    SELECT
        tw.stay_id,
        1 AS is_vent_preicu
    FROM target_window tw
    JOIN mimiciv_derived.ventilation v
      ON v.stay_id = tw.stay_id
    WHERE v.ventilation_status IN ('InvasiveVent', 'Tracheostomy')
      AND v.starttime <= tw.preicu_end
      AND v.endtime >= tw.preicu_start
    GROUP BY tw.stay_id
),
baseline_rrt AS (
    SELECT
        tw.stay_id,
        1 AS is_rrt_baseline
    FROM target_window tw
    JOIN mimiciv_derived.rrt r
      ON r.stay_id = tw.stay_id
    WHERE r.charttime >= tw.baseline_start
      AND r.charttime <= tw.baseline_end
    GROUP BY tw.stay_id
),
preicu_rrt AS (
    SELECT
        tw.stay_id,
        1 AS is_rrt_preicu
    FROM target_window tw
    JOIN mimiciv_derived.rrt r
      ON r.stay_id = tw.stay_id
    WHERE r.charttime >= tw.preicu_start
      AND r.charttime <= tw.preicu_end
    GROUP BY tw.stay_id
)
SELECT
    c.stay_id,
    COALESCE(pa.has_prior_albumin, 0) AS prior_albumin_use,
    COALESCE(v.is_vent_baseline, 0) AS mech_vent_baseline,
    COALESCE(r.is_rrt_baseline, 0) AS rrt_baseline,
    COALESCE(vp.is_vent_preicu, 0) AS mech_vent_preicu,
    COALESCE(rp.is_rrt_preicu, 0) AS rrt_preicu
FROM cohort_heart_failure_albumin_v2 c
LEFT JOIN prior_alb pa
  ON c.stay_id = pa.stay_id
LEFT JOIN baseline_vent v
  ON c.stay_id = v.stay_id
LEFT JOIN baseline_rrt r
  ON c.stay_id = r.stay_id
LEFT JOIN preicu_vent vp
  ON c.stay_id = vp.stay_id
LEFT JOIN preicu_rrt rp
  ON c.stay_id = rp.stay_id
ORDER BY c.stay_id;
