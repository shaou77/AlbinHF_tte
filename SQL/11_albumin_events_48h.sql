-- v2 event-level albumin administrations in the 0-48h exposure window.
-- This file is for audit/descriptive tables and concentration-specific checks.

WITH albumin_conc AS (
    SELECT 220861 AS itemid, '20%' AS concentration, 0.20::numeric AS conc_factor_g_per_ml, 1 AS hyperoncotic_flag UNION ALL
    SELECT 220862 AS itemid, '25%' AS concentration, 0.25::numeric AS conc_factor_g_per_ml, 1 AS hyperoncotic_flag UNION ALL
    SELECT 220863 AS itemid, '4%'  AS concentration, 0.04::numeric AS conc_factor_g_per_ml, 0 AS hyperoncotic_flag UNION ALL
    SELECT 220864 AS itemid, '5%'  AS concentration, 0.05::numeric AS conc_factor_g_per_ml, 0 AS hyperoncotic_flag
)
SELECT
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.icu_intime,
    ie.starttime,
    ie.endtime,
    EXTRACT(EPOCH FROM (ie.starttime - c.icu_intime)) / 3600.0 AS hours_from_icu_intime,
    ie.itemid,
    ac.concentration,
    ac.hyperoncotic_flag,
    ie.amount AS amount_ml,
    ie.amountuom,
    ie.rate,
    ie.rateuom,
    ie.ordercategoryname,
    ie.patientweight,
    (ie.amount * ac.conc_factor_g_per_ml) AS dose_g
FROM cohort_heart_failure_albumin_v2 c
JOIN mimiciv_icu.inputevents ie
  ON ie.stay_id = c.stay_id
JOIN albumin_conc ac
  ON ie.itemid = ac.itemid
WHERE ie.starttime >= c.icu_intime
  AND ie.starttime < c.icu_intime + INTERVAL '48 HOUR'
  AND ie.amount IS NOT NULL
  AND ie.amount > 0
ORDER BY c.stay_id, ie.starttime, ie.itemid;
