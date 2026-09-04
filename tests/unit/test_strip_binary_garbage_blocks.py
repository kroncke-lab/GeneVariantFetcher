"""Unit tests for the base-block degarbage script (harvest-time leaks)."""

from harvesting.supplement_fold import FOLD_BEGIN, FOLD_END
from scripts.strip_binary_garbage_blocks import strip_garbage_supplement_blocks


def _soup(n: int) -> str:
    unit = "\x00\x01g\x00\x00\x0b|ˇˇ\x02\x03\x04 –œ\x11‡°±\x1a·"
    return (unit * (n // len(unit) + 1))[:n]


PAPER = "# MAIN TEXT\n\nRYR2 proband cohort description.\n\n"
GOOD_BLOCK = "# SUPPLEMENTAL FILE 1: tableS1.csv\n\nvariant,carriers\nc.1A>G,3\n" + (
    "real rows of table text\n" * 60
)


def test_strips_garbage_block_keeps_paper_and_good_block():
    garbage_block = "# SUPPLEMENTAL FILE 2: leaked.doc\n\n" + _soup(50_000)
    text = PAPER + GOOD_BLOCK + "\n\n" + garbage_block
    cleaned, removed = strip_garbage_supplement_blocks(text)
    assert removed == ["leaked.doc"]
    assert "# MAIN TEXT" in cleaned
    assert "c.1A>G,3" in cleaned
    assert "\x00" not in cleaned


def test_no_change_when_all_blocks_are_clean():
    text = PAPER + GOOD_BLOCK
    cleaned, removed = strip_garbage_supplement_blocks(text)
    assert removed == []
    assert cleaned == text


def test_sentinel_fold_block_is_preserved_verbatim():
    garbage_block = "# SUPPLEMENTAL FILE 1: leaked.doc\n\n" + _soup(50_000)
    sentinel = (
        f"{FOLD_BEGIN}\n\n# FOLDED SUPPLEMENTS (re-extraction aid)\n"
        "# SUPPLEMENTAL FILE 1: clean.csv\n\nvariant,carriers\nc.2T>C,4\n"
        f"\n\n{FOLD_END}\n"
    )
    text = PAPER + garbage_block + "\n\n" + sentinel
    cleaned, removed = strip_garbage_supplement_blocks(text)
    assert removed == ["leaked.doc"]
    assert sentinel.strip() in cleaned
    assert "\x00" not in cleaned


def test_nested_file_granularity_keeps_clean_archive_members():
    archive_block = (
        "# SUPPLEMENTAL FILE 1: bundle.zip\n\n"
        "#### Nested file: good_member.csv\n\n"
        + ("variant,carriers\nc.1A>G,3\n" * 80)
        + "\n#### Nested file: leaked_member.doc\n\n"
        + _soup(60_000)
    )
    text = PAPER + archive_block
    cleaned, removed = strip_garbage_supplement_blocks(text)
    assert removed == ["bundle.zip :: leaked_member.doc"]
    assert "good_member.csv" in cleaned
    assert "c.1A>G,3" in cleaned
    assert "\x00" not in cleaned
