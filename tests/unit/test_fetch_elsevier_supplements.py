from pathlib import Path

import pytest

from harvesting.elsevier_api import ElsevierAPIClient
from harvesting.supplement_fold import FOLD_BEGIN
from scripts.fetch_elsevier_supplements import (
    MANIFEST_NAME,
    PaperTarget,
    _cached_complete_refs,
    _doi_for,
    _load_input,
    augment_paper,
)


class FakeElsevierClient:
    def __init__(self, xml: str):
        self.xml = xml
        self.download_calls: list[list[str]] = []
        self.fulltext_calls = 0

    def get_fulltext_by_doi(self, doi: str):
        self.fulltext_calls += 1
        return self.xml, None

    extract_supplement_refs = staticmethod(ElsevierAPIClient.extract_supplement_refs)
    supplement_local_name = staticmethod(ElsevierAPIClient.supplement_local_name)

    def download_supplements(self, xml_content: str, dest_dir: Path):
        refs = self.extract_supplement_refs(xml_content)
        downloaded: list[str] = []
        dest_dir.mkdir(parents=True, exist_ok=True)
        for ref in refs:
            path = dest_dir / self.supplement_local_name(ref)
            if path.exists() and path.stat().st_size > 1000:
                continue
            path.write_text(
                "variant,carriers\nc.2A>G,4\n" * 60,
                encoding="utf-8",
            )
            downloaded.append(path.name)
        self.download_calls.append(downloaded)
        return [dest_dir / name for name in downloaded]


def test_load_input_rejects_missing_file(tmp_path):
    missing = tmp_path / "missing.csv"
    with pytest.raises(SystemExit, match=f"Input file does not exist: {missing}"):
        _load_input(missing)


def test_doi_lookup_ignores_non_object_and_invalid_utf8_json(tmp_path):
    pmid = "12345678"
    (tmp_path / f"{pmid}_artifacts.json").write_text('["not", "an object"]')
    (tmp_path / "result.json").write_bytes(b"\xff\xfe invalid json")
    (tmp_path / f"{pmid}_FULL_CONTEXT.md").write_text(
        "# Article\n\nDOI: 10.1016/j.hrthm.2018.01.014\n"
    )

    assert _doi_for(tmp_path, pmid) == "10.1016/j.hrthm.2018.01.014"


def test_manifest_lookup_ignores_non_object_and_invalid_utf8_json(tmp_path):
    client = FakeElsevierClient("")
    manifest = tmp_path / MANIFEST_NAME
    manifest.write_text('["not", "an object"]')
    assert _cached_complete_refs(client, tmp_path, "10.1016/example") == []

    manifest.write_bytes(b"\xff\xfe invalid json")
    assert _cached_complete_refs(client, tmp_path, "10.1016/example") == []


def test_partial_elsevier_set_fetches_only_missing_and_self_folds(tmp_path):
    pmid = "29325976"
    original = "# MAIN TEXT\n\n## Results\n\nUsable article body.\n"
    (tmp_path / f"{pmid}_FULL_CONTEXT.md").write_text(original, encoding="utf-8")
    supp_dir = tmp_path / f"{pmid}_supplements"
    supp_dir.mkdir()
    (supp_dir / "mmc1.txt").write_text("existing supplement\n" * 80)
    xml = "1-s2.0-S1547527118300146-mmc1.txt 1-s2.0-S1547527118300146-mmc2.txt"
    client = FakeElsevierClient(xml)
    target = PaperTarget(pmid, tmp_path, "10.1016/j.hrthm.2018.01.014")

    result = augment_paper(client, target)

    assert result["new_files"] == ["mmc2.txt"]
    assert result["complete"] is True
    assert result["folded"] is True
    assert client.download_calls == [["mmc2.txt"]]
    folded = (tmp_path / f"{pmid}_FULL_CONTEXT.md").read_text(encoding="utf-8")
    assert folded.split(FOLD_BEGIN, 1)[0].rstrip() == original.rstrip()
    assert "mmc1.txt" in folded and "mmc2.txt" in folded

    result_again = augment_paper(client, target)

    assert result_again["new_files"] == []
    assert result_again["complete"] is True
    assert result_again["folded"] is False
    assert client.download_calls == [["mmc2.txt"]]
    assert client.fulltext_calls == 1
    assert (tmp_path / f"{pmid}_FULL_CONTEXT.md").read_text() == folded


def test_reuses_sibling_gene_supplement_before_network(tmp_path):
    pmid = "20541041"
    target_dir = tmp_path / "SCN5A" / pmid
    sibling_supp = tmp_path / "KCNQ1" / pmid / f"{pmid}_supplements"
    target_dir.mkdir(parents=True)
    sibling_supp.mkdir(parents=True)
    (target_dir / f"{pmid}_FULL_CONTEXT.md").write_text("# MAIN\n\nbody\n")
    (sibling_supp / "mmc1.txt").write_text("shared supplement\n" * 80)
    client = FakeElsevierClient("1-s2.0-S1547527118300146-mmc1.txt")
    target = PaperTarget(
        pmid,
        target_dir,
        "10.1016/j.hrthm.2010.06.013",
        (sibling_supp,),
    )

    result = augment_paper(client, target)

    assert result["complete"] is True
    assert result["new_files"] == ["mmc1.txt"]
    assert client.download_calls == []


def test_xml_failure_still_folds_local_supplements(tmp_path):
    pmid = "20541041"
    (tmp_path / f"{pmid}_FULL_CONTEXT.md").write_text("# MAIN\n\nbody\n")
    supp_dir = tmp_path / f"{pmid}_supplements"
    supp_dir.mkdir()
    (supp_dir / "mmc1.txt").write_text("variant,carriers\nc.2A>G,4\n" * 60)
    client = FakeElsevierClient("")
    target = PaperTarget(pmid, tmp_path, "10.1016/j.hrthm.2010.06.013")

    result = augment_paper(client, target)

    assert result["complete"] is False
    assert result["folded"] is True
    assert "c.2A>G,4" in (tmp_path / f"{pmid}_FULL_CONTEXT.md").read_text()
