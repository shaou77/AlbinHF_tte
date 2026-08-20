WITH target_window AS (
    SELECT
        c.stay_id,
        c.icu_intime,
        -- 严格限定为入科后 6 小时
        c.icu_intime + INTERVAL '6 HOUR' AS window_end
    FROM cohort_heart_failure_albumin_v2 c
),

-- 2. 提取并标准化利尿剂记录 (单位换算)
raw_diuretics AS (
    SELECT
        ie.stay_id,
        ie.starttime,
        ie.endtime,
        ie.amount,

        -- 计算总时长 (秒)
        EXTRACT(EPOCH FROM (ie.endtime - ie.starttime)) AS duration_sec,

        -- 转换系数: 将所有药物转换为 Furosemide Equivalents (mg)
        CASE
            -- Furosemide (呋塞米): 系数 1.0
            WHEN ie.itemid IN (221794, 228340) THEN 1.0
            -- Bumetanide (布美他尼): 1mg Bumex = 40mg Lasix -> 系数 40.0
            WHEN ie.itemid IN (221829, 229639) THEN 40.0
            -- Torsemide (托拉塞米): 20mg Torsemide = 40mg Lasix -> 系数 2.0 (虽然少见，建议加上)
            WHEN ie.itemid = 222315 THEN 2.0
            ELSE 0
        END AS conversion_factor
    FROM mimiciv_icu.inputevents ie
    JOIN target_window tw ON ie.stay_id = tw.stay_id
    WHERE
        ie.itemid IN (
            221794, 228340, -- Furosemide
            221829, 229639 -- Bumetanide
        )
        AND ie.amountuom = 'mg'
        -- 只要跟 0-6h 窗口有交集的记录都拿进来
        AND ie.starttime < tw.window_end
        AND ie.endtime > tw.icu_intime
),

-- 3. 计算 6小时内的有效剂量 (处理推注 vs 泵入)
diuretics_calc AS (
    SELECT
        rd.stay_id,

        CASE
            -- === 情况A: 瞬时推注 (Bolus) 或 极短时间给药 (< 2分钟) ===
            -- 逻辑：只要开始时间落在 0-6h 窗口内，就全部计入
            WHEN rd.duration_sec < 120 THEN
                CASE
                    WHEN rd.starttime >= tw.icu_intime AND rd.starttime < tw.window_end
                    THEN rd.amount * rd.conversion_factor
                    ELSE 0
                END

            -- === 情况B: 持续泵入 (Infusion) ===
            -- 逻辑：精确计算时间重叠比例
            ELSE
                (
                    GREATEST(0, EXTRACT(EPOCH FROM (
                        LEAST(rd.endtime, tw.window_end) - GREATEST(rd.starttime, tw.icu_intime)
                    )))
                    / NULLIF(rd.duration_sec, 0)
                ) * (rd.amount * rd.conversion_factor)
        END AS dose_6h_mg

    FROM raw_diuretics rd
    JOIN target_window tw ON rd.stay_id = tw.stay_id
),

-- 4. 聚合总剂量
diuretics_agg AS (
    SELECT
        stay_id,
        SUM(dose_6h_mg) AS total_furosemide_mg_6h
    FROM diuretics_calc
    GROUP BY stay_id
)

-- 5. 主查询
SELECT
    c.stay_id,

    -- 0-6h 总等效剂量 (mg)，保留1位小数
    ROUND(COALESCE(da.total_furosemide_mg_6h, 0)::numeric, 1) AS furosemide_eq_dose_6h,

    -- 是否使用Flag (作为二分类协变量也很重要)
    CASE
        WHEN da.total_furosemide_mg_6h > 0 THEN 1
        ELSE 0
    END AS diuretics_use_baseline

FROM cohort_heart_failure_albumin_v2 c
LEFT JOIN diuretics_agg da
    ON c.stay_id = da.stay_id
ORDER BY c.stay_id;
