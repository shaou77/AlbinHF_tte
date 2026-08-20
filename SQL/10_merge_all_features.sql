-- v2 core cohort/demographics extraction.
-- Feature SQL files 02-09 should be exported and joined by stay_id in R/Python.
-- This file deliberately avoids duplicating the full feature logic so each feature
-- extraction remains auditable and easy to rerun.

WITH cohort AS (
    SELECT
        subject_id,
        hadm_id,
        stay_id,
        icu_intime,
        icu_outtime,
        albumin_min,
        albumin_lt_3_flag,
        albumin_lt_35_flag
    FROM cohort_heart_failure_albumin_v2
)
SELECT
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    i.gender,
    i.los_hospital,
    i.admission_age,
    i.race,
    i.hospital_expire_flag,
    i.los_icu,
    i.dod,
    c.icu_intime,
    c.icu_outtime,
    c.albumin_min,
    c.albumin_lt_3_flag,
    c.albumin_lt_35_flag,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 28
        THEN CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0)::integer
        ELSE 28
    END AS time_28,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 28
        THEN 1 ELSE 0
    END AS status_28,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 90
        THEN CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0)::integer
        ELSE 90
    END AS time_90,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 90
        THEN 1 ELSE 0
    END AS status_90,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 365
        THEN CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0)::integer
        ELSE 365
    END AS time_365,
    CASE
        WHEN i.dod IS NOT NULL
         AND CEIL(EXTRACT(EPOCH FROM (i.dod::timestamp - c.icu_intime)) / 86400.0) BETWEEN 0 AND 365
        THEN 1 ELSE 0
    END AS status_365,
    gc.gcs_min AS gcs,
    fh.height,
    fw.weight,
    ch.charlson_comorbidity_index AS charlson,
    ii.first_careunit AS icu_type
FROM cohort c
LEFT JOIN mimiciv_derived.icustay_detail i
  ON c.stay_id = i.stay_id
LEFT JOIN mimiciv_derived.first_day_gcs gc
  ON c.stay_id = gc.stay_id
LEFT JOIN mimiciv_derived.first_day_height fh
  ON c.stay_id = fh.stay_id
LEFT JOIN mimiciv_derived.first_day_weight fw
  ON c.stay_id = fw.stay_id
LEFT JOIN mimiciv_derived.charlson ch
  ON c.hadm_id = ch.hadm_id
LEFT JOIN mimiciv_icu.icustays ii
  ON c.stay_id = ii.stay_id
ORDER BY c.stay_id;
