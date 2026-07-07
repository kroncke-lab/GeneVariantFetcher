from scripts.recall_audit.run_claim_debate_pilot import select_records


def test_select_records_filters_to_escalation_queue_keys():
    records = [
        {
            "model": "azure_ai/gpt-5.4",
            "gene": "KCNH2",
            "pmid": "1",
            "variant": "p.Arg1Trp",
            "failure_class": "count_semantics",
        },
        {
            "model": "azure_ai/gpt-5.4",
            "gene": "KCNH2",
            "pmid": "2",
            "variant": "p.Arg2Trp",
            "failure_class": "count_semantics",
        },
    ]

    selected = select_records(
        records,
        baseline_models=None,
        pmids=None,
        failure_classes=None,
        queue_keys={("KCNH2", "2", "p.Arg2Trp")},
        limit_cards=0,
    )

    assert selected == [records[1]]
