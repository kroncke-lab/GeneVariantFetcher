-- DuckDB queries used to reproduce the bounded analytical datasets rendered in
-- the Codex high-carrier 48-paper report artifact. Run from the repository root.

CREATE OR REPLACE TEMP VIEW codex_paper_metrics AS
SELECT *
FROM read_csv_auto(
  'benchmarks/codex_paper_eval/runs/codex_20260723_highcarrier48_routed/paper_metrics.csv',
  header = true
);

CREATE OR REPLACE TEMP VIEW codex_model_comparison AS
SELECT *
FROM read_csv_auto(
  'benchmarks/codex_paper_eval/runs/codex_20260723_highcarrier48_routed/model_comparison.csv',
  header = true
);

-- Overview metrics.
WITH totals AS (
  SELECT
    sum(tp)::DOUBLE AS tp,
    sum(fp)::DOUBLE AS fp,
    sum(fn)::DOUBLE AS fn,
    sum(total_tokens)::DOUBLE AS total_tokens,
    sum(elapsed_seconds)::DOUBLE AS summed_seconds,
    sum(carriers_predicted)::DOUBLE AS carrier_values,
    sum(carriers_gold_asserted)::DOUBLE AS carrier_assertions,
    sum(affected_predicted)::DOUBLE AS affected_values,
    sum(affected_gold_asserted)::DOUBLE AS affected_assertions,
    sum(unaffected_predicted)::DOUBLE AS unaffected_values,
    sum(unaffected_gold_asserted)::DOUBLE AS unaffected_assertions
  FROM codex_paper_metrics
)
SELECT
  tp / (tp + fp) AS precision,
  tp / (tp + fn) AS recall,
  2 * tp / (2 * tp + fp + fn) AS f1,
  tp::INTEGER AS tp,
  fp::INTEGER AS fp,
  fn::INTEGER AS fn,
  (tp + fn)::INTEGER AS gold_variants,
  (tp + fp)::INTEGER AS predicted_variants,
  carrier_values / carrier_assertions AS carriers_count_recall,
  affected_values / affected_assertions AS affected_count_recall,
  unaffected_values / unaffected_assertions AS unaffected_count_recall,
  total_tokens::BIGINT AS total_tokens,
  total_tokens / 48 AS tokens_per_paper,
  summed_seconds / 60 AS summed_paper_minutes
FROM totals;

-- Per-gene metrics and count coverage.
WITH gene_totals AS (
  SELECT
    gene,
    sum(tp)::DOUBLE AS tp,
    sum(fp)::DOUBLE AS fp,
    sum(fn)::DOUBLE AS fn,
    sum(carriers_predicted)::DOUBLE AS carrier_values,
    sum(carriers_gold_asserted)::DOUBLE AS carrier_assertions,
    sum(affected_predicted)::DOUBLE AS affected_values,
    sum(affected_gold_asserted)::DOUBLE AS affected_assertions,
    sum(unaffected_predicted)::DOUBLE AS unaffected_values,
    sum(unaffected_gold_asserted)::DOUBLE AS unaffected_assertions
  FROM codex_paper_metrics
  GROUP BY gene
)
SELECT
  gene,
  tp::INTEGER AS tp,
  fp::INTEGER AS fp,
  fn::INTEGER AS fn,
  tp / (tp + fp) AS precision,
  tp / (tp + fn) AS recall,
  2 * tp / (2 * tp + fp + fn) AS f1,
  carrier_values / carrier_assertions AS carriers_count_recall,
  affected_values / affected_assertions AS affected_count_recall,
  unaffected_values / unaffected_assertions AS unaffected_count_recall
FROM gene_totals
ORDER BY gene;

-- Tidy gene precision/recall/F1 rows for the grouped bar chart.
WITH gene_totals AS (
  SELECT gene, sum(tp)::DOUBLE AS tp, sum(fp)::DOUBLE AS fp, sum(fn)::DOUBLE AS fn
  FROM codex_paper_metrics
  GROUP BY gene
),
gene_metrics AS (
  SELECT
    gene,
    tp,
    fp,
    fn,
    tp + fn AS gold_variants,
    tp / (tp + fp) AS precision,
    tp / (tp + fn) AS recall,
    2 * tp / (2 * tp + fp + fn) AS f1
  FROM gene_totals
)
SELECT gene, 'Precision' AS metric, precision AS value, tp, fp, fn, gold_variants
FROM gene_metrics
UNION ALL
SELECT gene, 'Recall', recall, tp, fp, fn, gold_variants
FROM gene_metrics
UNION ALL
SELECT gene, 'F1', f1, tp, fp, fn, gold_variants
FROM gene_metrics
ORDER BY gene, metric;

-- Overall count coverage and conditional MAE/RMSE.
WITH count_rows AS (
  SELECT
    'Carriers' AS field,
    carriers_predicted AS supplied,
    carriers_gold_asserted AS gold_asserted,
    carriers_mae AS mae,
    carriers_rmse AS rmse
  FROM codex_paper_metrics
  UNION ALL
  SELECT
    'Affected',
    affected_predicted,
    affected_gold_asserted,
    affected_mae,
    affected_rmse
  FROM codex_paper_metrics
  UNION ALL
  SELECT
    'Unaffected',
    unaffected_predicted,
    unaffected_gold_asserted,
    unaffected_mae,
    unaffected_rmse
  FROM codex_paper_metrics
)
SELECT
  field,
  sum(supplied)::INTEGER AS supplied,
  sum(gold_asserted)::INTEGER AS gold_asserted,
  sum(supplied)::DOUBLE / sum(gold_asserted) AS count_recall,
  sum(coalesce(mae, 0) * supplied) / nullif(sum(supplied), 0) AS mae,
  sqrt(
    sum(coalesce(rmse, 0) * coalesce(rmse, 0) * supplied)
    / nullif(sum(supplied), 0)
  ) AS rmse
FROM count_rows
GROUP BY field
ORDER BY count_recall DESC;

-- Complete per-paper table.
SELECT
  gene,
  pmid::VARCHAR AS pmid,
  tool,
  tp,
  fp,
  fn,
  precision,
  recall,
  f1,
  carriers_recall AS carriers_count_recall,
  affected_recall AS affected_count_recall,
  unaffected_recall AS unaffected_count_recall,
  elapsed_seconds AS seconds,
  total_tokens AS tokens
FROM codex_paper_metrics
ORDER BY fn DESC, gene, pmid;

-- Codex and Claude aggregate comparison.
SELECT
  system,
  score_view,
  precision,
  recall,
  f1,
  extracted_variants AS extracted,
  total_tokens,
  tokens_per_paper,
  elapsed_minutes AS wall_minutes,
  comparability_note
FROM codex_model_comparison
ORDER BY recall DESC;
