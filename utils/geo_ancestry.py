"""Deterministic carrier-ethnicity and cohort-origin extraction.

Three cheap, no-LLM signals, tried in order of trust:

1. Explicit ethnicity/ancestry stated in the text (``find_ethnicity``).
2. Explicit country/region of origin stated in the text (``find_origin``).
3. **Last resort** — the corresponding author's country, parsed from a PubMed
   affiliation string (``country_from_affiliation``). Clearly weaker than a
   stated cohort origin, so callers mark it as inferred.

Used to backfill ``individual_records.ethnicity`` / ``geographic_origin`` and a
paper-level cohort origin at migrate time, so existing extractions gain this
without a re-run. The gazetteers are deliberately curated (common study
populations), not exhaustive — extend as new cohorts appear.
"""

from __future__ import annotations

import re
from typing import List, Optional, Tuple

# Ethnicity / ancestry terms. Longer, more specific phrases first so
# "Ashkenazi Jewish" wins over "Jewish" and "African American" over "African".
ETHNICITY_TERMS: Tuple[str, ...] = (
    "Ashkenazi Jewish",
    "Ashkenazi",
    "Sephardic Jewish",
    "Sephardic",
    "African American",
    "African-American",
    "Han Chinese",
    "East Asian",
    "South Asian",
    "Southeast Asian",
    "Native American",
    "Pacific Islander",
    "Middle Eastern",
    "Hispanic",
    "Latino",
    "Latina",
    "Latinx",
    "Caucasian",
    "European",
    "African",
    "Chinese",
    "Japanese",
    "Korean",
    "Vietnamese",
    "Thai",
    "Filipino",
    "Malay",
    "Indian",
    "Pakistani",
    "Iranian",
    "Turkish",
    "Arab",
    "Bedouin",
    "Druze",
    "Roma",
    "Finnish",
    "Amish",
    "Mennonite",
    "Māori",
    "Maori",
    "Aboriginal",
    "Indigenous",
    "Jewish",
    "Slavic",
    "White",
    "Black",
    "Asian",
)

# Country / region gazetteer for stated origins and affiliation parsing.
# Value is the canonical display name; keys include common variants/demonyms.
_COUNTRY_CANON = {
    "usa": "United States",
    "u.s.a": "United States",
    "u.s.a.": "United States",
    "u.s.": "United States",
    "us": "United States",
    "united states": "United States",
    "united states of america": "United States",
    "america": "United States",
    "uk": "United Kingdom",
    "u.k.": "United Kingdom",
    "united kingdom": "United Kingdom",
    "great britain": "United Kingdom",
    "england": "United Kingdom",
    "scotland": "United Kingdom",
    "wales": "United Kingdom",
    "northern ireland": "United Kingdom",
    "the netherlands": "Netherlands",
    "holland": "Netherlands",
    "netherlands": "Netherlands",
    "p.r. china": "China",
    "pr china": "China",
    "p.r.china": "China",
    "people's republic of china": "China",
    "prc": "China",
    "china": "China",
    "republic of korea": "South Korea",
    "south korea": "South Korea",
    "korea": "South Korea",
    "republic of ireland": "Ireland",
    "ireland": "Ireland",
    "russian federation": "Russia",
    "russia": "Russia",
    "czech republic": "Czech Republic",
    "czechia": "Czech Republic",
}
# Straightforward single-token countries (canonical == title-cased key).
_COUNTRIES = (
    "Japan",
    "Germany",
    "France",
    "Italy",
    "Spain",
    "Portugal",
    "Belgium",
    "Switzerland",
    "Austria",
    "Sweden",
    "Norway",
    "Denmark",
    "Finland",
    "Iceland",
    "Poland",
    "Hungary",
    "Greece",
    "Turkey",
    "Israel",
    "Egypt",
    "Iran",
    "Iraq",
    "India",
    "Pakistan",
    "Bangladesh",
    "Thailand",
    "Vietnam",
    "Malaysia",
    "Singapore",
    "Indonesia",
    "Philippines",
    "Taiwan",
    "Canada",
    "Mexico",
    "Brazil",
    "Argentina",
    "Chile",
    "Colombia",
    "Australia",
    "Zealand",
    "Africa",
    "Nigeria",
    "Kenya",
    "Morocco",
    "Tunisia",
    "Algeria",
    "Croatia",
    "Serbia",
    "Slovenia",
    "Slovakia",
    "Romania",
    "Bulgaria",
    "Ukraine",
    "Lithuania",
    "Latvia",
    "Estonia",
    "Lebanon",
    "Jordan",
    "Saudi",
    "Qatar",
    "Emirates",
    "Kuwait",
    "Oman",
)
for _c in _COUNTRIES:
    _COUNTRY_CANON.setdefault(_c.lower(), _c)
# A couple of multi-word canonicalizations.
_COUNTRY_CANON.setdefault("new zealand", "New Zealand")
_COUNTRY_CANON.setdefault("south africa", "South Africa")
_COUNTRY_CANON.setdefault("saudi arabia", "Saudi Arabia")
_COUNTRY_CANON.setdefault("united arab emirates", "United Arab Emirates")

_MULTIWORD_ORIGINS = (
    "United States",
    "United Kingdom",
    "New Zealand",
    "South Korea",
    "South Africa",
    "Saudi Arabia",
    "United Arab Emirates",
    "Czech Republic",
    "Northern Italy",
    "Southern Italy",
    "Eastern Europe",
    "Western Europe",
    "Northern Europe",
    "Southern Europe",
    "North America",
    "South America",
    "Central Europe",
    "Mainland China",
)


def _compile_alt(terms) -> re.Pattern:
    alts = sorted(terms, key=len, reverse=True)
    body = "|".join(re.escape(t) for t in alts)
    return re.compile(rf"(?<![A-Za-z])({body})(?![A-Za-z])", re.IGNORECASE)


_ETHNICITY_RE = _compile_alt(ETHNICITY_TERMS)
_ORIGIN_RE = _compile_alt(
    list(_MULTIWORD_ORIGINS) + [c.title() for c in _COUNTRY_CANON] + list(_COUNTRIES)
)


def _canonical_origin(raw: str) -> str:
    key = re.sub(r"\.\s*$", "", raw.strip()).lower()
    if key in _COUNTRY_CANON:
        return _COUNTRY_CANON[key]
    # Multi-word regions keep their title-cased form.
    return " ".join(w.capitalize() if w.islower() else w for w in raw.split())


def find_ethnicity(*texts: Optional[str]) -> Optional[str]:
    """First stated race/ethnicity/ancestry across the given text fields, or None."""
    for text in texts:
        if not text:
            continue
        match = _ETHNICITY_RE.search(str(text))
        if match:
            return match.group(1)
    return None


def find_origin(*texts: Optional[str]) -> Optional[str]:
    """First stated country/region of origin across the given text fields."""
    for text in texts:
        if not text:
            continue
        match = _ORIGIN_RE.search(str(text))
        if match:
            return _canonical_origin(match.group(1))
    return None


def country_from_affiliation(affiliation: Optional[str]) -> Optional[str]:
    """Best-effort country from a PubMed affiliation string.

    Affiliations conventionally end with the country ("..., City, Country."), so
    scan the trailing comma-separated tokens for a known country. Also matches a
    country appearing anywhere as a fallback.
    """
    if not affiliation:
        return None
    text = str(affiliation).strip().rstrip(".")
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    for token in reversed(tokens[-3:]):  # country is almost always at the tail
        canon = _COUNTRY_CANON.get(token.lower())
        if canon:
            return canon
    match = _ORIGIN_RE.search(text)
    return _canonical_origin(match.group(1)) if match else None


def enrich_ethnicity_origin(
    *,
    stated_ethnicity: Optional[str] = None,
    stated_origin: Optional[str] = None,
    text_fields: Optional[List[Optional[str]]] = None,
    author_affiliation: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """Resolve (ethnicity, origin) using stated values, then text, then author country.

    A value derived only from the author's affiliation is tagged
    "(inferred from author affiliation)" so downstream consumers can flag its lower
    confidence. Returns (ethnicity, origin); either may be None.
    """
    text_fields = text_fields or []
    ethnicity = stated_ethnicity or find_ethnicity(*text_fields)
    origin = stated_origin or find_origin(*text_fields)
    if not origin:
        author_country = country_from_affiliation(author_affiliation)
        if author_country:
            origin = f"{author_country} (inferred from author affiliation)"
    return ethnicity, origin
