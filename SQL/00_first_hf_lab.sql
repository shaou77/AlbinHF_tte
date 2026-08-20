-- Reproduce the locally generated mimiciv_derived.first_hf_lab table.
--
-- This query adapts the official MIMIC Code first_day_lab concept by using
-- the study baseline window: 12 hours before through 6 hours after ICU
-- admission. It returns one row per ICU stay.
--
-- Prerequisites: the official MIMIC-IV v3.1 raw tables and derived concepts
-- complete_blood_count, chemistry, blood_differential, coagulation, and
-- enzyme must already exist.
--
-- This statement intentionally does not overwrite an existing table. Drop or
-- rename mimiciv_derived.first_hf_lab explicitly before rebuilding it.

CREATE TABLE mimiciv_derived.first_hf_lab AS
WITH cbc AS (
    SELECT
        ie.stay_id,
        MIN(le.hematocrit) AS hematocrit_min,
        MAX(le.hematocrit) AS hematocrit_max,
        MIN(le.hemoglobin) AS hemoglobin_min,
        MAX(le.hemoglobin) AS hemoglobin_max,
        MIN(le.platelet) AS platelets_min,
        MAX(le.platelet) AS platelets_max,
        MIN(le.wbc) AS wbc_min,
        MAX(le.wbc) AS wbc_max
    FROM mimiciv_icu.icustays ie
    LEFT JOIN mimiciv_derived.complete_blood_count le
      ON le.subject_id = ie.subject_id
     AND le.charttime >= ie.intime - INTERVAL '12 HOUR'
     AND le.charttime <= ie.intime + INTERVAL '6 HOUR'
    GROUP BY ie.stay_id
),
chem AS (
    SELECT
        ie.stay_id,
        MIN(le.albumin) AS albumin_min,
        MAX(le.albumin) AS albumin_max,
        MIN(le.globulin) AS globulin_min,
        MAX(le.globulin) AS globulin_max,
        MIN(le.total_protein) AS total_protein_min,
        MAX(le.total_protein) AS total_protein_max,
        MIN(le.aniongap) AS aniongap_min,
        MAX(le.aniongap) AS aniongap_max,
        MIN(le.bicarbonate) AS bicarbonate_min,
        MAX(le.bicarbonate) AS bicarbonate_max,
        MIN(le.bun) AS bun_min,
        MAX(le.bun) AS bun_max,
        MIN(le.calcium) AS calcium_min,
        MAX(le.calcium) AS calcium_max,
        MIN(le.chloride) AS chloride_min,
        MAX(le.chloride) AS chloride_max,
        MIN(le.creatinine) AS creatinine_min,
        MAX(le.creatinine) AS creatinine_max,
        MIN(le.glucose) AS glucose_min,
        MAX(le.glucose) AS glucose_max,
        MIN(le.sodium) AS sodium_min,
        MAX(le.sodium) AS sodium_max,
        MIN(le.potassium) AS potassium_min,
        MAX(le.potassium) AS potassium_max
    FROM mimiciv_icu.icustays ie
    LEFT JOIN mimiciv_derived.chemistry le
      ON le.subject_id = ie.subject_id
     AND le.charttime >= ie.intime - INTERVAL '12 HOUR'
     AND le.charttime <= ie.intime + INTERVAL '6 HOUR'
    GROUP BY ie.stay_id
),
diff AS (
    SELECT
        ie.stay_id,
        MIN(le.basophils_abs) AS abs_basophils_min,
        MAX(le.basophils_abs) AS abs_basophils_max,
        MIN(le.eosinophils_abs) AS abs_eosinophils_min,
        MAX(le.eosinophils_abs) AS abs_eosinophils_max,
        MIN(le.lymphocytes_abs) AS abs_lymphocytes_min,
        MAX(le.lymphocytes_abs) AS abs_lymphocytes_max,
        MIN(le.monocytes_abs) AS abs_monocytes_min,
        MAX(le.monocytes_abs) AS abs_monocytes_max,
        MIN(le.neutrophils_abs) AS abs_neutrophils_min,
        MAX(le.neutrophils_abs) AS abs_neutrophils_max,
        MIN(le.atypical_lymphocytes) AS atyps_min,
        MAX(le.atypical_lymphocytes) AS atyps_max,
        MIN(le.bands) AS bands_min,
        MAX(le.bands) AS bands_max,
        MIN(le.immature_granulocytes) AS imm_granulocytes_min,
        MAX(le.immature_granulocytes) AS imm_granulocytes_max,
        MIN(le.metamyelocytes) AS metas_min,
        MAX(le.metamyelocytes) AS metas_max,
        MIN(le.nrbc) AS nrbc_min,
        MAX(le.nrbc) AS nrbc_max
    FROM mimiciv_icu.icustays ie
    LEFT JOIN mimiciv_derived.blood_differential le
      ON le.subject_id = ie.subject_id
     AND le.charttime >= ie.intime - INTERVAL '12 HOUR'
     AND le.charttime <= ie.intime + INTERVAL '6 HOUR'
    GROUP BY ie.stay_id
),
coag AS (
    SELECT
        ie.stay_id,
        MIN(le.d_dimer) AS d_dimer_min,
        MAX(le.d_dimer) AS d_dimer_max,
        MIN(le.fibrinogen) AS fibrinogen_min,
        MAX(le.fibrinogen) AS fibrinogen_max,
        MIN(le.thrombin) AS thrombin_min,
        MAX(le.thrombin) AS thrombin_max,
        MIN(le.inr) AS inr_min,
        MAX(le.inr) AS inr_max,
        MIN(le.pt) AS pt_min,
        MAX(le.pt) AS pt_max,
        MIN(le.ptt) AS ptt_min,
        MAX(le.ptt) AS ptt_max
    FROM mimiciv_icu.icustays ie
    LEFT JOIN mimiciv_derived.coagulation le
      ON le.subject_id = ie.subject_id
     AND le.charttime >= ie.intime - INTERVAL '12 HOUR'
     AND le.charttime <= ie.intime + INTERVAL '6 HOUR'
    GROUP BY ie.stay_id
),
enz AS (
    SELECT
        ie.stay_id,
        MIN(le.alt) AS alt_min,
        MAX(le.alt) AS alt_max,
        MIN(le.alp) AS alp_min,
        MAX(le.alp) AS alp_max,
        MIN(le.ast) AS ast_min,
        MAX(le.ast) AS ast_max,
        MIN(le.amylase) AS amylase_min,
        MAX(le.amylase) AS amylase_max,
        MIN(le.bilirubin_total) AS bilirubin_total_min,
        MAX(le.bilirubin_total) AS bilirubin_total_max,
        MIN(le.bilirubin_direct) AS bilirubin_direct_min,
        MAX(le.bilirubin_direct) AS bilirubin_direct_max,
        MIN(le.bilirubin_indirect) AS bilirubin_indirect_min,
        MAX(le.bilirubin_indirect) AS bilirubin_indirect_max,
        MIN(le.ck_cpk) AS ck_cpk_min,
        MAX(le.ck_cpk) AS ck_cpk_max,
        MIN(le.ck_mb) AS ck_mb_min,
        MAX(le.ck_mb) AS ck_mb_max,
        MIN(le.ggt) AS ggt_min,
        MAX(le.ggt) AS ggt_max,
        MIN(le.ld_ldh) AS ld_ldh_min,
        MAX(le.ld_ldh) AS ld_ldh_max
    FROM mimiciv_icu.icustays ie
    LEFT JOIN mimiciv_derived.enzyme le
      ON le.subject_id = ie.subject_id
     AND le.charttime >= ie.intime - INTERVAL '12 HOUR'
     AND le.charttime <= ie.intime + INTERVAL '6 HOUR'
    GROUP BY ie.stay_id
)
SELECT
    ie.subject_id,
    ie.stay_id,
    cbc.hematocrit_min,
    cbc.hematocrit_max,
    cbc.hemoglobin_min,
    cbc.hemoglobin_max,
    cbc.platelets_min,
    cbc.platelets_max,
    cbc.wbc_min,
    cbc.wbc_max,
    chem.albumin_min,
    chem.albumin_max,
    chem.globulin_min,
    chem.globulin_max,
    chem.total_protein_min,
    chem.total_protein_max,
    chem.aniongap_min,
    chem.aniongap_max,
    chem.bicarbonate_min,
    chem.bicarbonate_max,
    chem.bun_min,
    chem.bun_max,
    chem.calcium_min,
    chem.calcium_max,
    chem.chloride_min,
    chem.chloride_max,
    chem.creatinine_min,
    chem.creatinine_max,
    chem.glucose_min,
    chem.glucose_max,
    chem.sodium_min,
    chem.sodium_max,
    chem.potassium_min,
    chem.potassium_max,
    diff.abs_basophils_min,
    diff.abs_basophils_max,
    diff.abs_eosinophils_min,
    diff.abs_eosinophils_max,
    diff.abs_lymphocytes_min,
    diff.abs_lymphocytes_max,
    diff.abs_monocytes_min,
    diff.abs_monocytes_max,
    diff.abs_neutrophils_min,
    diff.abs_neutrophils_max,
    diff.atyps_min,
    diff.atyps_max,
    diff.bands_min,
    diff.bands_max,
    diff.imm_granulocytes_min,
    diff.imm_granulocytes_max,
    diff.metas_min,
    diff.metas_max,
    diff.nrbc_min,
    diff.nrbc_max,
    coag.d_dimer_min,
    coag.d_dimer_max,
    coag.fibrinogen_min,
    coag.fibrinogen_max,
    coag.thrombin_min,
    coag.thrombin_max,
    coag.inr_min,
    coag.inr_max,
    coag.pt_min,
    coag.pt_max,
    coag.ptt_min,
    coag.ptt_max,
    enz.alt_min,
    enz.alt_max,
    enz.alp_min,
    enz.alp_max,
    enz.ast_min,
    enz.ast_max,
    enz.amylase_min,
    enz.amylase_max,
    enz.bilirubin_total_min,
    enz.bilirubin_total_max,
    enz.bilirubin_direct_min,
    enz.bilirubin_direct_max,
    enz.bilirubin_indirect_min,
    enz.bilirubin_indirect_max,
    enz.ck_cpk_min,
    enz.ck_cpk_max,
    enz.ck_mb_min,
    enz.ck_mb_max,
    enz.ggt_min,
    enz.ggt_max,
    enz.ld_ldh_min,
    enz.ld_ldh_max
FROM mimiciv_icu.icustays ie
LEFT JOIN cbc
  ON ie.stay_id = cbc.stay_id
LEFT JOIN chem
  ON ie.stay_id = chem.stay_id
LEFT JOIN diff
  ON ie.stay_id = diff.stay_id
LEFT JOIN coag
  ON ie.stay_id = coag.stay_id
LEFT JOIN enz
  ON ie.stay_id = enz.stay_id;
