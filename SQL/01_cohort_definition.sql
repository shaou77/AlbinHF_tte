-- v2 cohort definition
-- Purpose:
--   1) support the manuscript definition of hypoalbuminemia as albumin < 3.5 g/dL
--   2) preserve albumin_min in the cohort view so R can run
--      sensitivity analyses using either <3.0 or <3.5 g/dL without re-extracting.

DROP MATERIALIZED VIEW IF EXISTS cohort_heart_failure_albumin_v2;

CREATE MATERIALIZED VIEW cohort_heart_failure_albumin_v2 AS
WITH hf_codes (icd_code, icd_version) AS (
    VALUES
    ('39891', 9), ('40201', 9), ('40211', 9), ('40291', 9), ('40401', 9), ('40403', 9),
    ('40411', 9), ('40413', 9), ('40491', 9), ('40493', 9), ('4280', 9), ('4281', 9),
    ('42820', 9), ('42821', 9), ('42822', 9), ('42823', 9), ('42830', 9), ('42831', 9),
    ('42832', 9), ('42833', 9), ('42840', 9), ('42841', 9), ('42842', 9), ('42843', 9),
    ('4289', 9), ('I0981', 10), ('I110', 10), ('I130', 10), ('I132', 10), ('I50', 10),
    ('I502', 10), ('I5020', 10), ('I5021', 10), ('I5022', 10), ('I5023', 10),
    ('I503', 10), ('I5030', 10), ('I5031', 10), ('I5032', 10), ('I5033', 10),
    ('I504', 10), ('I5040', 10), ('I5041', 10), ('I5042', 10), ('I5043', 10),
    ('I508', 10), ('I5081', 10), ('I50810', 10), ('I50811', 10), ('I50812', 10),
    ('I50813', 10), ('I50814', 10), ('I5082', 10), ('I5083', 10), ('I5084', 10),
    ('I5089', 10), ('I509', 10), ('I9713', 10), ('I97130', 10), ('I97131', 10),
    ('T8622', 10), ('T8632', 10)
),
cirrhosis_codes (icd_code, icd_version) AS (
    VALUES
    ('5712', 9), ('5715', 9), ('5716', 9), ('K703', 10), ('K7030', 10), ('K7031', 10),
    ('K717', 10), ('K74', 10), ('K743', 10), ('K744', 10), ('K745', 10), ('K746', 10),
    ('K7460', 10), ('K7469', 10), ('P7881', 10)
),
shock_codes (icd_code, icd_version) AS (
    VALUES
    ('78552', 9), ('9584', 9), ('99802', 9), ('R6521', 10), ('T7501', 10),
    ('T7501XA', 10), ('T7501XD', 10), ('T7501XS', 10), ('T794', 10), ('T794XXA', 10),
    ('T794XXD', 10), ('T794XXS', 10)
),
target_admissions AS (
    SELECT DISTINCT d.hadm_id
    FROM mimiciv_hosp.diagnoses_icd d
    JOIN hf_codes hf
      ON d.icd_code = hf.icd_code
     AND d.icd_version = hf.icd_version
    WHERE NOT EXISTS (
        SELECT 1
        FROM mimiciv_hosp.diagnoses_icd e
        JOIN cirrhosis_codes c
          ON e.icd_code = c.icd_code
         AND e.icd_version = c.icd_version
        WHERE e.hadm_id = d.hadm_id
    )
      AND NOT EXISTS (
        SELECT 1
        FROM mimiciv_hosp.diagnoses_icd e
        JOIN shock_codes s
          ON e.icd_code = s.icd_code
         AND e.icd_version = s.icd_version
        WHERE e.hadm_id = d.hadm_id
    )
),
ranked_stays AS (
    SELECT
        id.subject_id,
        id.hadm_id,
        id.stay_id,
        id.icu_intime,
        id.icu_outtime,
        fdl.albumin_min,
        ROW_NUMBER() OVER (PARTITION BY id.subject_id ORDER BY id.icu_intime ASC) AS rn
    FROM mimiciv_derived.icustay_detail id
    JOIN target_admissions ta
      ON id.hadm_id = ta.hadm_id
    JOIN mimiciv_derived.first_hf_lab fdl
      ON id.stay_id = fdl.stay_id
    WHERE fdl.albumin_min < 3.5
)
SELECT
    subject_id,
    hadm_id,
    stay_id,
    icu_intime,
    icu_outtime,
    albumin_min,
    CASE WHEN albumin_min < 3.0 THEN 1 ELSE 0 END AS albumin_lt_3_flag,
    CASE WHEN albumin_min < 3.5 THEN 1 ELSE 0 END AS albumin_lt_35_flag
FROM ranked_stays
WHERE rn = 1;

CREATE INDEX idx_cohort_hf_albumin_v2_stay_id ON cohort_heart_failure_albumin_v2(stay_id);
CREATE INDEX idx_cohort_hf_albumin_v2_hadm_id ON cohort_heart_failure_albumin_v2(hadm_id);
