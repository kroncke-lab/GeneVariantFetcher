from types import SimpleNamespace
from unittest.mock import MagicMock

from harvesting.browser_html.strategies.generic import GenericStrategy


class _FakePage:
    url = "https://karger.com/crd/article-abstract/133/2/73/97064/example"

    def goto(self, *args, **kwargs):
        return None

    def wait_for_timeout(self, *args, **kwargs):
        return None

    def content(self):
        return (
            "<html><body><article>"
            + ("KCNH2 full text " * 100)
            + "</article></body></html>"
        )


def test_generic_strategy_routes_karger_to_karger_supplement_scraper(monkeypatch):
    strategy = GenericStrategy()
    monkeypatch.setattr(strategy, "accept_cookies", lambda page: False)
    monkeypatch.setattr(
        strategy, "extract_via_scraper", lambda *args, **kwargs: "body " * 300
    )
    monkeypatch.setattr(strategy, "download_figures", lambda *args, **kwargs: [])

    scraper = MagicMock()
    scraper.scrape_karger_supplements.return_value = [
        {"url": "https://ndownloader.figshare.com/files/1", "name": "table.doc"}
    ]
    ctx = SimpleNamespace(
        doi="10.1159/000440608",
        timeout_s=1,
        scraper=scraper,
        output_dir=None,
        pmid="26496715",
    )

    result = strategy.fetch(_FakePage(), ctx)

    scraper.scrape_karger_supplements.assert_called_once()
    scraper.scrape_generic_supplements.assert_not_called()
    assert result.supp_files[0]["name"] == "table.doc"
