-- v2 baseline hematology/coagulation/enzyme extraction
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
cbc_data AS (
    SELECT
        tw.stay_id,
        CASE WHEN le.charttime >= tw.baseline_start AND le.charttime <= tw.baseline_end THEN 1 ELSE 0 END AS in_baseline_window,
        CASE WHEN le.charttime >= tw.preicu_start AND le.charttime <= tw.preicu_end THEN 1 ELSE 0 END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (le.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - le.charttime)) AS time_before_icu,
        le.wbc,
        le.hemoglobin,
        le.hematocrit,
        le.platelet
    FROM target_window tw
    LEFT JOIN mimiciv_derived.complete_blood_count le
      ON le.subject_id = tw.subject_id
     AND le.charttime >= tw.baseline_start
     AND le.charttime <= tw.baseline_end
),
diff_data AS (
    SELECT
        tw.stay_id,
        CASE WHEN le.charttime >= tw.baseline_start AND le.charttime <= tw.baseline_end THEN 1 ELSE 0 END AS in_baseline_window,
        CASE WHEN le.charttime >= tw.preicu_start AND le.charttime <= tw.preicu_end THEN 1 ELSE 0 END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (le.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - le.charttime)) AS time_before_icu,
        le.neutrophils_abs,
        le.lymphocytes_abs
    FROM target_window tw
    LEFT JOIN mimiciv_derived.blood_differential le
      ON le.subject_id = tw.subject_id
     AND le.charttime >= tw.baseline_start
     AND le.charttime <= tw.baseline_end
),
coag_data AS (
    SELECT
        tw.stay_id,
        CASE WHEN le.charttime >= tw.baseline_start AND le.charttime <= tw.baseline_end THEN 1 ELSE 0 END AS in_baseline_window,
        CASE WHEN le.charttime >= tw.preicu_start AND le.charttime <= tw.preicu_end THEN 1 ELSE 0 END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (le.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - le.charttime)) AS time_before_icu,
        le.inr,
        le.pt,
        le.ptt,
        le.fibrinogen,
        le.d_dimer
    FROM target_window tw
    LEFT JOIN mimiciv_derived.coagulation le
      ON le.subject_id = tw.subject_id
     AND le.charttime >= tw.baseline_start
     AND le.charttime <= tw.baseline_end
),
enz_data AS (
    SELECT
        tw.stay_id,
        CASE WHEN le.charttime >= tw.baseline_start AND le.charttime <= tw.baseline_end THEN 1 ELSE 0 END AS in_baseline_window,
        CASE WHEN le.charttime >= tw.preicu_start AND le.charttime <= tw.preicu_end THEN 1 ELSE 0 END AS in_preicu_window,
        ABS(EXTRACT(EPOCH FROM (le.charttime - tw.icu_intime))) AS time_diff_abs,
        EXTRACT(EPOCH FROM (tw.icu_intime - le.charttime)) AS time_before_icu,
        le.alt,
        le.ast,
        le.alp,
        le.ggt,
        le.ld_ldh,
        le.bilirubin_total,
        le.bilirubin_direct,
        le.ck_cpk,
        le.ck_mb
    FROM target_window tw
    LEFT JOIN mimiciv_derived.enzyme le
      ON le.subject_id = tw.subject_id
     AND le.charttime >= tw.baseline_start
     AND le.charttime <= tw.baseline_end
)
SELECT
    tw.stay_id,

    (ARRAY_AGG(cbc.wbc ORDER BY cbc.time_diff_abs ASC) FILTER (WHERE cbc.in_baseline_window = 1 AND cbc.wbc IS NOT NULL))[1] AS wbc_baseline,
    (ARRAY_AGG(cbc.hemoglobin ORDER BY cbc.time_diff_abs ASC) FILTER (WHERE cbc.in_baseline_window = 1 AND cbc.hemoglobin IS NOT NULL))[1] AS hemoglobin_baseline,
    (ARRAY_AGG(cbc.hematocrit ORDER BY cbc.time_diff_abs ASC) FILTER (WHERE cbc.in_baseline_window = 1 AND cbc.hematocrit IS NOT NULL))[1] AS hematocrit_baseline,
    (ARRAY_AGG(cbc.platelet ORDER BY cbc.time_diff_abs ASC) FILTER (WHERE cbc.in_baseline_window = 1 AND cbc.platelet IS NOT NULL))[1] AS platelet_baseline,
    (ARRAY_AGG(diff.neutrophils_abs ORDER BY diff.time_diff_abs ASC) FILTER (WHERE diff.in_baseline_window = 1 AND diff.neutrophils_abs IS NOT NULL))[1] AS neutrophils_abs_baseline,
    (ARRAY_AGG(diff.lymphocytes_abs ORDER BY diff.time_diff_abs ASC) FILTER (WHERE diff.in_baseline_window = 1 AND diff.lymphocytes_abs IS NOT NULL))[1] AS lymphocytes_abs_baseline,
    (ARRAY_AGG(coag.inr ORDER BY coag.time_diff_abs ASC) FILTER (WHERE coag.in_baseline_window = 1 AND coag.inr IS NOT NULL))[1] AS inr_baseline,
    (ARRAY_AGG(coag.pt ORDER BY coag.time_diff_abs ASC) FILTER (WHERE coag.in_baseline_window = 1 AND coag.pt IS NOT NULL))[1] AS pt_baseline,
    (ARRAY_AGG(coag.ptt ORDER BY coag.time_diff_abs ASC) FILTER (WHERE coag.in_baseline_window = 1 AND coag.ptt IS NOT NULL))[1] AS aptt_baseline,
    (ARRAY_AGG(coag.fibrinogen ORDER BY coag.time_diff_abs ASC) FILTER (WHERE coag.in_baseline_window = 1 AND coag.fibrinogen IS NOT NULL))[1] AS fibrinogen_baseline,
    (ARRAY_AGG(coag.d_dimer ORDER BY coag.time_diff_abs ASC) FILTER (WHERE coag.in_baseline_window = 1 AND coag.d_dimer IS NOT NULL))[1] AS d_dimer_baseline,
    (ARRAY_AGG(enz.alt ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.alt IS NOT NULL))[1] AS alt_baseline,
    (ARRAY_AGG(enz.ast ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.ast IS NOT NULL))[1] AS ast_baseline,
    (ARRAY_AGG(enz.alp ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.alp IS NOT NULL))[1] AS alp_baseline,
    (ARRAY_AGG(enz.ggt ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.ggt IS NOT NULL))[1] AS ggt_baseline,
    (ARRAY_AGG(enz.ld_ldh ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.ld_ldh IS NOT NULL))[1] AS ldh_baseline,
    (ARRAY_AGG(enz.bilirubin_total ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.bilirubin_total IS NOT NULL))[1] AS bilirubin_total_enz_baseline,
    (ARRAY_AGG(enz.bilirubin_direct ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.bilirubin_direct IS NOT NULL))[1] AS bilirubin_direct_baseline,
    (ARRAY_AGG(enz.ck_cpk ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.ck_cpk IS NOT NULL))[1] AS ck_cpk_baseline,
    (ARRAY_AGG(enz.ck_mb ORDER BY enz.time_diff_abs ASC) FILTER (WHERE enz.in_baseline_window = 1 AND enz.ck_mb IS NOT NULL))[1] AS ck_mb_baseline,

    (ARRAY_AGG(cbc.wbc ORDER BY cbc.time_before_icu ASC) FILTER (WHERE cbc.in_preicu_window = 1 AND cbc.wbc IS NOT NULL))[1] AS wbc_preicu,
    (ARRAY_AGG(cbc.hemoglobin ORDER BY cbc.time_before_icu ASC) FILTER (WHERE cbc.in_preicu_window = 1 AND cbc.hemoglobin IS NOT NULL))[1] AS hemoglobin_preicu,
    (ARRAY_AGG(cbc.hematocrit ORDER BY cbc.time_before_icu ASC) FILTER (WHERE cbc.in_preicu_window = 1 AND cbc.hematocrit IS NOT NULL))[1] AS hematocrit_preicu,
    (ARRAY_AGG(cbc.platelet ORDER BY cbc.time_before_icu ASC) FILTER (WHERE cbc.in_preicu_window = 1 AND cbc.platelet IS NOT NULL))[1] AS platelet_preicu,
    (ARRAY_AGG(diff.neutrophils_abs ORDER BY diff.time_before_icu ASC) FILTER (WHERE diff.in_preicu_window = 1 AND diff.neutrophils_abs IS NOT NULL))[1] AS neutrophils_abs_preicu,
    (ARRAY_AGG(diff.lymphocytes_abs ORDER BY diff.time_before_icu ASC) FILTER (WHERE diff.in_preicu_window = 1 AND diff.lymphocytes_abs IS NOT NULL))[1] AS lymphocytes_abs_preicu,
    (ARRAY_AGG(coag.inr ORDER BY coag.time_before_icu ASC) FILTER (WHERE coag.in_preicu_window = 1 AND coag.inr IS NOT NULL))[1] AS inr_preicu,
    (ARRAY_AGG(coag.pt ORDER BY coag.time_before_icu ASC) FILTER (WHERE coag.in_preicu_window = 1 AND coag.pt IS NOT NULL))[1] AS pt_preicu,
    (ARRAY_AGG(coag.ptt ORDER BY coag.time_before_icu ASC) FILTER (WHERE coag.in_preicu_window = 1 AND coag.ptt IS NOT NULL))[1] AS aptt_preicu,
    (ARRAY_AGG(coag.fibrinogen ORDER BY coag.time_before_icu ASC) FILTER (WHERE coag.in_preicu_window = 1 AND coag.fibrinogen IS NOT NULL))[1] AS fibrinogen_preicu,
    (ARRAY_AGG(coag.d_dimer ORDER BY coag.time_before_icu ASC) FILTER (WHERE coag.in_preicu_window = 1 AND coag.d_dimer IS NOT NULL))[1] AS d_dimer_preicu,
    (ARRAY_AGG(enz.alt ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.alt IS NOT NULL))[1] AS alt_preicu,
    (ARRAY_AGG(enz.ast ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.ast IS NOT NULL))[1] AS ast_preicu,
    (ARRAY_AGG(enz.alp ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.alp IS NOT NULL))[1] AS alp_preicu,
    (ARRAY_AGG(enz.ggt ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.ggt IS NOT NULL))[1] AS ggt_preicu,
    (ARRAY_AGG(enz.ld_ldh ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.ld_ldh IS NOT NULL))[1] AS ldh_preicu,
    (ARRAY_AGG(enz.bilirubin_total ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.bilirubin_total IS NOT NULL))[1] AS bilirubin_total_enz_preicu,
    (ARRAY_AGG(enz.bilirubin_direct ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.bilirubin_direct IS NOT NULL))[1] AS bilirubin_direct_preicu,
    (ARRAY_AGG(enz.ck_cpk ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.ck_cpk IS NOT NULL))[1] AS ck_cpk_preicu,
    (ARRAY_AGG(enz.ck_mb ORDER BY enz.time_before_icu ASC) FILTER (WHERE enz.in_preicu_window = 1 AND enz.ck_mb IS NOT NULL))[1] AS ck_mb_preicu
FROM target_window tw
LEFT JOIN cbc_data cbc
  ON tw.stay_id = cbc.stay_id
LEFT JOIN diff_data diff
  ON tw.stay_id = diff.stay_id
LEFT JOIN coag_data coag
  ON tw.stay_id = coag.stay_id
LEFT JOIN enz_data enz
  ON tw.stay_id = enz.stay_id
GROUP BY tw.stay_id
ORDER BY tw.stay_id;
