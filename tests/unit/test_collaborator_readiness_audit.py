import json
import sqlite3

from scripts.collaborator_readiness_audit import audit_gene


def _build_run(tmp_path, *, gene="BRCA1"):
    run_dir = tmp_path / gene / "run"
    run_dir.mkdir(parents=True)
    manifest = tmp_path / f"{gene}.txt"
    manifest.write_text("111\n222\n")
    (run_dir / f"{gene}_pmids.txt").write_text("111\n222\n")
    extractions = run_dir / "extractions"
    extractions.mkdir()
    for pmid in ("111", "222"):
        (extractions / f"{gene}_PMID_{pmid}.json").write_text(
            json.dumps(
                {
                    "paper_metadata": {"pmid": pmid, "title": f"Paper {pmid}"},
                    "variants": [] if pmid == "222" else [{"variant": "c.1A>G"}],
                    "extraction_metadata": {"total_variants_found": int(pmid == "111")},
                }
            )
        )
    (run_dir / "RUN_STATUS.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "exit_code": 0,
                "stage_failures": [],
                "stage_warnings": [],
                "active_db": f"{gene}.db",
                "llm_trace": {
                    "integrity_level": "write_time_verified",
                    "integrity_errors": [],
                    "llm_call_count": 2,
                    "decision_event_count": 2,
                },
            }
        )
    )
    (run_dir / "source_completeness.json").write_text(
        json.dumps(
            {
                "total_pmids_discovered": 2,
                "papers_with_fulltext": 1,
                "papers_abstract_only": 1,
                "abstract_only_pmids": ["222"],
            }
        )
    )
    (run_dir / f"{gene}_workflow.log").write_text(
        "vf-enrich: {'enriched': True, 'quarantined': True}\n"
    )

    con = sqlite3.connect(run_dir / f"{gene}.db")
    con.executescript(
        """
        CREATE TABLE papers (pmid TEXT PRIMARY KEY, title TEXT);
        CREATE TABLE variants (
          variant_id INTEGER PRIMARY KEY, gene_symbol TEXT, cdna_notation TEXT,
          protein_notation TEXT, genomic_position TEXT, structural_description TEXT
        );
        CREATE TABLE variant_papers (variant_id INTEGER, pmid TEXT);
        CREATE TABLE penetrance_data (
          penetrance_id INTEGER PRIMARY KEY, variant_id INTEGER, pmid TEXT,
          total_carriers_observed INTEGER, affected_count INTEGER,
          unaffected_count INTEGER, uncertain_count INTEGER,
          penetrance_percentage REAL, trust_tier TEXT, field_trust TEXT
        );
        CREATE TABLE age_dependent_penetrance (penetrance_percentage REAL);
        CREATE TABLE vf_enrichment (
          variant_id INTEGER, matched INTEGER, fp_class TEXT
        );
        CREATE TABLE quarantined_variants (fp_class TEXT);
        INSERT INTO papers VALUES ('111', 'Human BRCA1 cohort');
        INSERT INTO papers VALUES ('222', 'Human BRCA1 negative cohort');
        INSERT INTO variants VALUES (1, 'BRCA1', 'c.1A>G', 'p.Met1Val', NULL, NULL);
        INSERT INTO variant_papers VALUES (1, '111');
        INSERT INTO penetrance_data VALUES (
          1, 1, '111', 1, 1, 1, NULL, NULL, 'quarantine',
          '{"total_carriers":"quarantine","affected":"quarantine","unaffected":"quarantine"}'
        );
        INSERT INTO vf_enrichment VALUES (1, 1, NULL);
        INSERT INTO quarantined_variants VALUES ('wrong_gene_residue_mismatch');
        """
    )
    con.commit()
    con.close()
    return run_dir, manifest


def test_audit_passes_when_raw_contradiction_is_quarantined(tmp_path):
    run_dir, manifest = _build_run(tmp_path)

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is True
    assert result["cohort"]["exact_order"] is True
    assert result["counts"]["impossible_partitions_raw"] == 1
    assert result["counts"]["impossible_partitions_unquarantined"] == 0
    assert result["variant_features"]["quarantined_variants"] == {
        "wrong_gene_residue_mismatch": 1
    }
    assert result["variant_features"]["uncovered_live_variants"] == 0
    assert result["variant_features"]["unmatched_live_variants"] == 0
    assert result["variant_features"]["blocked_live_variants"] == 0
    assert result["vf_gate_passed"] is True


def test_audit_uses_primary_trace_after_empty_validation_resume(tmp_path):
    run_dir, manifest = _build_run(tmp_path)
    (run_dir / "RUN_STATUS.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "exit_code": 0,
                "stage_failures": [],
                "stage_warnings": [],
                "active_db": "BRCA1.db",
                "llm_trace": {
                    "run_id": "validation-resume",
                    "integrity_level": "generated_now",
                    "integrity_errors": [],
                    "llm_call_count": 0,
                    "decision_event_count": 0,
                },
            }
        )
    )
    trace_dir = run_dir / "llm_traces"
    trace_dir.mkdir()
    trace_manifest = trace_dir / "trace_manifest.json"
    trace_manifest.write_text(
        json.dumps(
            {
                "run_id": "original-extraction",
                "trace_root": str(trace_dir),
                "llm_call_count": 351,
                "decision_event_count": 378,
                "verification": {
                    "level": "write_time_verified",
                    "errors": [],
                },
            }
        )
    )

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is True
    assert result["trace"] == {
        "run_id": "original-extraction",
        "evidence_source": str(trace_manifest),
        "integrity_level": "write_time_verified",
        "integrity_errors": [],
        "llm_calls": 351,
        "decision_events": 378,
        "missing_decision_links": [],
    }


def test_audit_fails_manifest_drift_and_exposed_bad_rows(tmp_path):
    run_dir, manifest = _build_run(tmp_path)
    (run_dir / "BRCA1_pmids.txt").write_text("222\n111\n")
    (run_dir / "extractions" / "BRCA1_PMID_222.json").unlink()
    con = sqlite3.connect(run_dir / "BRCA1.db")
    con.execute("INSERT INTO variants VALUES (2, 'BRCA1', NULL, NULL, NULL, NULL)")
    con.execute(
        """INSERT INTO penetrance_data VALUES
           (2, 1, '111', 1, 1, 1, NULL, NULL, 'trusted', '{}')"""
    )
    con.commit()
    con.close()

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is False
    failures = " ".join(result["structural_gate_failures"])
    assert "membership or order" in failures
    assert "extraction JSON" in failures
    assert "enrichment/quarantine" in failures
    assert "nameless" in failures
    assert "unquarantined" in failures


def test_audit_fails_when_species_scoped_gene_paper_retains_links(tmp_path):
    run_dir, manifest = _build_run(tmp_path)
    abstract_dir = run_dir / "abstract_json"
    abstract_dir.mkdir()
    (abstract_dir / "222.json").write_text(
        json.dumps(
            {
                "metadata": {
                    "title": "Single nucleotide variation in canine BRCA1 tissue"
                }
            }
        )
    )
    con = sqlite3.connect(run_dir / "BRCA1.db")
    con.execute("UPDATE papers SET title='' WHERE pmid='222'")
    con.execute("INSERT INTO variants VALUES (2, 'BRCA1', 'c.2A>G', NULL, NULL, NULL)")
    con.execute("INSERT INTO variant_papers VALUES (2, '222')")
    con.commit()
    con.close()

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is False
    assert result["paper_scope"]["nonhuman_target_title_links"] == [
        {
            "pmid": "222",
            "title": "Single nucleotide variation in canine BRCA1 tissue",
            "variant_links": 1,
        }
    ]
    assert "non-human target-gene" in " ".join(result["structural_gate_failures"])


def test_audit_fails_empty_or_unbound_extraction_claims(tmp_path):
    run_dir, manifest = _build_run(tmp_path)
    (run_dir / "extractions" / "BRCA1_PMID_222.json").write_text("{}")
    con = sqlite3.connect(run_dir / "BRCA1.db")
    con.execute("DELETE FROM papers WHERE pmid='222'")
    con.commit()
    con.close()

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is False
    failures = " ".join(result["structural_gate_failures"])
    assert "empty extraction payload" in failures
    assert "papers table PMID set" in failures


def test_audit_fails_unverified_empty_trace_and_unmatched_live_variant(tmp_path):
    run_dir, manifest = _build_run(tmp_path)
    status = json.loads((run_dir / "RUN_STATUS.json").read_text())
    status["llm_trace"].update(
        {"integrity_level": "generated_now", "llm_call_count": 0}
    )
    (run_dir / "RUN_STATUS.json").write_text(json.dumps(status))
    con = sqlite3.connect(run_dir / "BRCA1.db")
    con.execute("INSERT INTO variants VALUES (2, 'BRCA1', 'c.2A>G', NULL, NULL, NULL)")
    con.execute("INSERT INTO variant_papers VALUES (2, '111')")
    con.execute("INSERT INTO vf_enrichment VALUES (2, 0, 'residue_mismatch')")
    con.commit()
    con.close()

    result = audit_gene("BRCA1", run_dir, manifest)

    assert result["structural_gate_passed"] is False
    assert result["vf_gate_passed"] is False
    assert result["variant_features"]["unmatched_live_variants"] == 1
    assert result["variant_features"]["blocked_live_variants"] == 1
    failures = " ".join(result["structural_gate_failures"])
    assert "not write-time verified" in failures
    assert "no model calls" in failures
