from scripts.recall_audit.fn_root_cause import render


def test_render_labels_paper_derived_primary_lane() -> None:
    rendered = render(
        {"run_id": "paper-run", "primary_score_lane": "paper_derived"}, []
    )

    assert "0 paper-derived primary false negatives" in rendered
    assert "| paper lane | linkage lane |" in rendered


def test_render_labels_legacy_trusted_projection() -> None:
    rendered = render({"run_id": "legacy-run"}, [])

    assert "0 legacy trusted-projection false negatives" in rendered
    assert "| scored lane | linkage lane |" in rendered
