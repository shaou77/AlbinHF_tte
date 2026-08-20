-- v2 albumin exposure extraction
-- Adds concentration-specific exposure details needed for reviewer comments.

WITH albumin_conc AS (
    SELECT 220861 AS itemid, '20%' AS concentration, 0.20::numeric AS conc_factor_g_per_ml, 1 AS hyperoncotic_flag UNION ALL
    SELECT 220862 AS itemid, '25%' AS concentration, 0.25::numeric AS conc_factor_g_per_ml, 1 AS hyperoncotic_flag UNION ALL
    SELECT 220863 AS itemid, '4%'  AS concentration, 0.04::numeric AS conc_factor_g_per_ml, 0 AS hyperoncotic_flag UNION ALL
    SELECT 220864 AS itemid, '5%'  AS concentration, 0.05::numeric AS conc_factor_g_per_ml, 0 AS hyperoncotic_flag
),
albumin_events_48h AS (
    SELECT
        c.stay_id,
        ie.starttime,
        ie.endtime,
        ie.itemid,
        ac.concentration,
        ac.hyperoncotic_flag,
        ie.amount AS amount_ml,
        ie.amountuom,
        (ie.amount * ac.conc_factor_g_per_ml) AS dose_g
    FROM mimiciv_icu.inputevents ie
    JOIN cohort_heart_failure_albumin_v2 c
      ON ie.stay_id = c.stay_id
    JOIN albumin_conc ac
      ON ie.itemid = ac.itemid
    WHERE ie.starttime >= c.icu_intime
      AND ie.starttime < c.icu_intime + INTERVAL '48 HOUR'
      AND ie.amount IS NOT NULL
      AND ie.amount > 0
),
albumin_48h_usage AS (
    SELECT
        stay_id,
        SUM(dose_g) AS total_albumin_g_exact,
        COUNT(*) AS albumin_doses_48h,
        MIN(starttime) AS first_albumin_time,
        (ARRAY_AGG(itemid ORDER BY starttime, itemid))[1] AS first_albumin_itemid,
        (ARRAY_AGG(concentration ORDER BY starttime, itemid))[1] AS first_albumin_concentration,
        (ARRAY_AGG(amount_ml ORDER BY starttime, itemid))[1] AS first_albumin_amount_ml,
        (ARRAY_AGG(dose_g ORDER BY starttime, itemid))[1] AS first_albumin_dose_g,
        SUM(dose_g) FILTER (WHERE itemid = 220861) AS albumin_g_20pct_48h,
        SUM(dose_g) FILTER (WHERE itemid = 220862) AS albumin_g_25pct_48h,
        SUM(dose_g) FILTER (WHERE itemid = 220863) AS albumin_g_4pct_48h,
        SUM(dose_g) FILTER (WHERE itemid = 220864) AS albumin_g_5pct_48h,
        SUM(dose_g) FILTER (WHERE hyperoncotic_flag = 1) AS hyperoncotic_albumin_g_48h,
        SUM(dose_g) FILTER (WHERE hyperoncotic_flag = 0) AS iso_low_albumin_g_48h,
        COUNT(*) FILTER (WHERE hyperoncotic_flag = 1) AS hyperoncotic_albumin_doses_48h,
        COUNT(*) FILTER (WHERE hyperoncotic_flag = 0) AS iso_low_albumin_doses_48h
    FROM albumin_events_48h
    GROUP BY stay_id
)
SELECT
    c.stay_id,
    ROUND(COALESCE(u.total_albumin_g_exact, 0)) AS albumin_g_48h,
    COALESCE(u.total_albumin_g_exact, 0) AS albumin_g_48h_exact,
    CASE WHEN COALESCE(u.total_albumin_g_exact, 0) > 0 THEN 1 ELSE 0 END AS albumin_use_48h_flag,
    u.first_albumin_time,
    u.albumin_doses_48h,
    u.first_albumin_itemid,
    u.first_albumin_concentration,
    u.first_albumin_amount_ml,
    u.first_albumin_dose_g,
    COALESCE(u.albumin_g_20pct_48h, 0) AS albumin_g_20pct_48h,
    COALESCE(u.albumin_g_25pct_48h, 0) AS albumin_g_25pct_48h,
    COALESCE(u.albumin_g_4pct_48h, 0) AS albumin_g_4pct_48h,
    COALESCE(u.albumin_g_5pct_48h, 0) AS albumin_g_5pct_48h,
    COALESCE(u.hyperoncotic_albumin_g_48h, 0) AS hyperoncotic_albumin_g_48h,
    COALESCE(u.iso_low_albumin_g_48h, 0) AS iso_low_albumin_g_48h,
    COALESCE(u.hyperoncotic_albumin_doses_48h, 0) AS hyperoncotic_albumin_doses_48h,
    COALESCE(u.iso_low_albumin_doses_48h, 0) AS iso_low_albumin_doses_48h,
    CASE WHEN COALESCE(u.hyperoncotic_albumin_g_48h, 0) > 0 THEN 1 ELSE 0 END AS hyperoncotic_albumin_48h_flag,
    CASE WHEN COALESCE(u.iso_low_albumin_g_48h, 0) > 0 THEN 1 ELSE 0 END AS iso_low_albumin_48h_flag
FROM cohort_heart_failure_albumin_v2 c
LEFT JOIN albumin_48h_usage u
  ON c.stay_id = u.stay_id
ORDER BY c.stay_id;
