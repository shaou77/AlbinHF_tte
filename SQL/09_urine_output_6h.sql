WITH target_window AS (
    SELECT
        c.stay_id,
        c.icu_intime,
        -- 窗口：入科 0 到 +6 小时
        c.icu_intime + INTERVAL '6 HOUR' AS win_end
    FROM cohort_heart_failure_albumin_v2 c
),

-- 应用官方逻辑提取窗口内的原始记录
uo_raw AS (
    SELECT
        tw.stay_id,
        oe.charttime,
        -- 【官方逻辑】处理膀胱冲洗 (GU Irrigant)
        -- 如果是冲洗液输入(227488)，记为负值；其他记为正值
        CASE
            WHEN oe.itemid = 227488 AND oe.value > 0 THEN -1 * oe.value
            ELSE oe.value
        END AS urineoutput
    FROM target_window tw
    JOIN mimiciv_icu.outputevents oe
        ON tw.stay_id = oe.stay_id
    WHERE
        oe.itemid IN (
            226559, -- Foley (留置导尿)
            226560, -- Void (自解)
            226561, -- Condom Cath (避孕套导尿)
            226584, -- Ileoconduit (回肠导管)
            226563, -- Suprapubic (耻骨上造瘘)
            226564, -- R Nephrostomy (右肾造瘘)
            226565, -- L Nephrostomy (左肾造瘘)
            226567, -- Straight Cath (直导)
            226557, -- R Ureteral Stent (右输尿管支架)
            226558, -- L Ureteral Stent (左输尿管支架)
            227488, -- GU Irrigant Volume In (冲洗入量，需扣除)
            227489  -- GU Irrigant/Urine Volume Out (冲洗出量)
        )
        -- 时间限制：入科 0-6h
        AND oe.charttime >= tw.icu_intime
        AND oe.charttime <= tw.win_end
)

SELECT
    tw.stay_id,

    -- 1. 0-6h 总尿量 (mL)
    -- 使用 COALESCE 处理没有记录的情况，默认为 0
    ROUND(COALESCE(SUM(u.urineoutput), 0)::NUMERIC, 1) AS urineoutput_6h_total,

    -- 2. 0-6h 平均尿量流速 (mL/hr)
    -- 既然是固定 6h 窗口，直接除以 6 即可
    ROUND((COALESCE(SUM(u.urineoutput), 0) / 6.0)::NUMERIC, 1) AS urineoutput_6h_rate

FROM target_window tw
LEFT JOIN uo_raw u ON tw.stay_id = u.stay_id
GROUP BY tw.stay_id
ORDER BY tw.stay_id;
