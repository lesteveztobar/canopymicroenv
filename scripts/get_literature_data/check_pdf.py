"""
check_pdf.py
------------
Scans a literature PDF (species checklist / epiphyte ecology paper) and
stages candidate observation rows in the data/csv/combinedv3.csv schema,
for a human to review before merging.

WHY A STAGING TOOL AND NOT A FULLY AUTOMATIC EXTRACTOR
This is a *checker*, not an oracle. The papers we need data from format
their tables wildly differently (Johansson-zone abundance matrices,
host x epiphyte matrices, AGH-range columns, per-record Darwin-Core-style
"Observations." citations, indicator-value tables with no height data at
all, ...). Regex/table heuristics cannot reliably disambiguate all of
that. So this script surfaces what it *can* find, tags every row with
where each value came from and what's missing, and leaves judgment calls
to a human. Always eyeball the staging CSV against the PDF before
merging rows into data/csv/combinedv3.csv.

WHAT COUNTS AS "COMPLETE"
A row is complete if it has all four of:
    1. lat/lon               (site-level coordinates count; per-record beats them)
    2. a date                (day-precision beats month beats year-range)
    3. a species ID          (Genus + species epithet)
    4. a height reference    (meters above ground, OR a Johansson-zone
                              abundance breakdown recorded in `note`)
Rows missing any of these are still written (so you can see what the PDF
*does* say about a species) but are flagged in `_missing_required`.

USAGE
    python check_pdf.py paper.pdf
    python check_pdf.py paper.pdf --source "Komada et al. 2025" --site "Jardin Botanique A, Ankarafantsika NP"
    python check_pdf.py *.pdf --outdir data/literature_staging
    python check_pdf.py paper.pdf --dump-candidates     # raw regex hits, for manual cross-checking

OUTPUT
    One staging CSV per PDF: <outdir>/<pdf_stem>_staging.csv
    Columns = combinedv3.csv columns, plus QC columns prefixed with "_"
    (_page, _confidence, _missing_required, _lat_lon_basis, _date_basis,
    _height_basis, _raw_evidence). Drop the "_"-prefixed columns before
    appending reviewed rows into data/csv/combinedv3.csv.

DEPENDENCIES
    pip install pdfplumber pandas
"""

import argparse
import os
import re
import sys
from datetime import date

try:
    import pdfplumber
except ImportError:
    pdfplumber = None

import pandas as pd


# ── OUTPUT SCHEMA ────────────────────────────────────────────────────────────

COMBINEDV3_COLUMNS = [
    'Source', 'Area_or_Site', 'lat', 'lon', 'Elevation_m', 'FieldID',
    'Abundance', 'Height_m', 'CanopyHeight_m', 'note', 'Genus', 'species',
    'FinalID', 'datetime',
]
QC_COLUMNS = [
    '_page', '_confidence', '_missing_required',
    '_lat_lon_basis', '_date_basis', '_height_basis', '_raw_evidence',
]
STAGING_COLUMNS = COMBINEDV3_COLUMNS + QC_COLUMNS


# ── STOPLISTS (false-positive filters for the binomial-name regex) ─────────

GENUS_STOPWORDS = {
    "the", "this", "these", "those", "table", "figure", "fig", "section",
    "results", "methods", "discussion", "introduction", "conclusion",
    "conclusions", "acknowledgements", "acknowledgments", "references",
    "appendix", "supplementary", "additional", "data", "material",
    "materials", "observation", "observations", "identification",
    "remarks", "distribution", "host", "hosts", "note", "notes",
    "abstract", "keywords", "key", "january", "february", "march", "april",
    "may", "june", "july", "august", "september", "october", "november",
    "december", "national", "botanical", "check", "list", "journal",
    "ecological", "biotropica", "plant", "plants", "ecology", "science",
    "sciences", "annals", "we", "in", "on", "at", "as", "for", "with",
    "from", "were", "was", "are", "is", "and", "or", "but", "study",
    "area", "species", "site", "sites", "forest", "forests", "each",
    "all", "some", "many", "most", "epiphyte", "epiphytes", "epiphytic",
    "vascular", "several", "these", "andean", "montane", "tropical",
    "field", "guide", "conservation", "biology", "biological", "review",
    "diversity", "richness", "canopy", "zone", "zones", "figure1",
    "author", "authors", "corresponding", "university", "institute",
    "department", "faculty", "received", "accepted", "published", "doi",
    # Common Spanish/French/German function words and title-case nouns
    # that keep showing up as fake "genus" hits from reference-list
    # entries, institution names, and PDF download watermarks.
    "del", "de", "la", "las", "los", "el", "un", "una", "en", "y", "und",
    "des", "der", "die", "das",
    "informe", "secretaria", "secretaría", "centro", "centenario",
    "sustentabilidad", "protección", "proteccion", "especies", "nativas",
    "drought", "tolerance", "neotropical", "rican", "orchid", "life",
    "comparisons", "among", "indicate", "four", "type", "form", "final",
    "statement", "ethical", "archiv", "asuntos", "universitats",
    "landesbibliothek", "download", "downloaded", "terms", "conditions",
    "wiley", "library", "online", "elsevier", "springer", "copyright",
    "editor", "editors", "editorial",
}
SPECIES_EPITHET_STOPWORDS = {
    "et", "al", "sp", "spp", "cf", "aff", "nov", "var", "subsp", "ex",
    "auct", "sensu", "lato", "stricto",
    "del", "de", "la", "las", "los", "el", "un", "una", "en", "y", "und",
    "des", "der", "die", "das",
    "informe", "secretaria", "secretaría", "centro", "centenario",
    "sustentabilidad", "protección", "proteccion", "especies", "nativas",
    "drought", "tolerance", "neotropical", "rican", "orchid", "life",
    "comparisons", "among", "indicate", "four", "fruiting", "bark",
    "type", "form", "final", "statement", "ethical", "naturwissenschaftlicher",
    "download", "downloaded", "terms", "conditions", "ambiental",
}
# Real genus/species epithets essentially never run this long as a single
# word; anything past it is almost certainly a rotated-page spacing
# reconstruction failure (multiple words merged with no space at all).
MAX_GENUS_LEN = 20
MAX_EPITHET_LEN = 25
AUTHOR_TAIL_RE = re.compile(
    r'^\s*(?:\([A-ZÀ-Ý][A-Za-zà-ÿ\.]*\)\s*)?'
    r'[A-ZÀ-Ý][A-Za-zà-ÿ]*\.?(?:\s*[&,]\s*[A-ZÀ-Ý][A-Za-zà-ÿ\.]*)*'
)
# A parenthetical family name right after a binomial -- "Baudouinia
# fluggeiformis (Fabaceae)" -- looks like AUTHOR_TAIL_RE's optional
# "(X.)" author-abbreviation group, but it isn't one: it's the family
# annotation used throughout "Host taxa." lists. Whatever capitalized
# word happens to follow (often just the *next* list item's genus) must
# not be mistaken for a continuing author citation.
FAMILY_PAREN_RE = re.compile(r'^\s*\([A-ZÀ-Ý][a-zà-ÿ]+(?:aceae|idae)\)')


def _has_author_tail(tail):
    if FAMILY_PAREN_RE.match(tail):
        return False
    return bool(AUTHOR_TAIL_RE.match(tail))

# ── DATE / MONTH TABLES ──────────────────────────────────────────────────────

MONTH_EN = {m.lower(): i + 1 for i, m in enumerate([
    'January', 'February', 'March', 'April', 'May', 'June', 'July',
    'August', 'September', 'October', 'November', 'December'])}
MONTH_ES = {
    'enero': 1, 'febrero': 2, 'marzo': 3, 'abril': 4, 'mayo': 5,
    'junio': 6, 'julio': 7, 'agosto': 8, 'septiembre': 9, 'setiembre': 9,
    'octubre': 10, 'noviembre': 11, 'diciembre': 12,
}
ROMAN_MONTHS = {
    'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6,
    'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12,
}

_MONTH_EN_ALT = '|'.join(MONTH_EN)
_MONTH_ES_ALT = '|'.join(MONTH_ES)

DATE_EXCLUDE_RE = re.compile(
    r'\b(received|accepted|revised|published|available\s+online|'
    r'recibido|aceptado|revisado|copyright|doi|issn|citar\s+como|'
    r'cite\s+as|editor)\b', re.I)
# Botanical/zoological nomenclature citations end "... Journal Vol: Page
# (Year)" -- e.g. "Ann. Fac. Sc. Marseille21: 215 (1912)". That year is a
# publication date for the *name*, not an observation date. Exclude bare
# years immediately preceded by this "N: M (" citation shape.
NOMENCLATURE_CITATION_RE = re.compile(r'\d+[a-z]?\s*:\s*\d+[a-z]?\s*\($')

DATE_DMY_ROMAN_RE = re.compile(r'\b(\d{1,2})\.([IVXLCDM]{1,4})\.(\d{4})\b')
DATE_DMY_EN_RE = re.compile(
    rf'\b(\d{{1,2}})\s+({_MONTH_EN_ALT})\.?,?\s+(\d{{4}})\b', re.I)
DATE_DMY_ES_RE = re.compile(
    rf'\b(\d{{1,2}})\s+de\s+({_MONTH_ES_ALT})\s+de\s+(\d{{4}})\b', re.I)
DATE_MY_EN_RE = re.compile(rf'\b({_MONTH_EN_ALT})\s+(\d{{4}})\b', re.I)
DATE_MY_ES_RE = re.compile(rf'\b({_MONTH_ES_ALT})\s+(?:de\s+)?(\d{{4}})\b', re.I)
YEAR_RANGE_RE = re.compile(r'\b((?:19|20)\d{2})\s*[-–—/]\s*((?:19|20)?\d{2})\b')
YEAR_RE = re.compile(r'\b((?:19|20)\d{2})\b')

# Precision ranking, higher = better.
PRECISION_RANK = {'day': 4, 'month': 3, 'year_range': 2, 'year': 1}

# Keywords that mark a sentence as describing fieldwork (used to pick the
# best *document-level* fallback date when no per-species date is found).
FIELDWORK_KEYWORDS_RE = re.compile(
    r'\b(survey|surveys|sampled|sampling|collect(?:ed|ion)?|fieldwork|'
    r'field\s+work|census|conducted|field\s+survey|study\s+period|'
    r'muestre[oó]|recolect|trabajo\s+de\s+campo)\b', re.I)


# ── COORDINATE PARSING ────────────────────────────────────────────────────────

_NUM = r'\d{1,3}(?:\.\d+)?'
# Prime / double-prime / accent glyphs vary wildly across PDF text
# extractors (', ′, ´, ́ , curly quotes...). Match generically: any short
# run of non-alnum, non-whitespace punctuation -- except dash variants,
# which must NOT be swallowed here: a dash between two minute/second
# values means "this is a range" (see _COORD_RANGE_TOKEN_RE below), and
# treating it as a decorative mark would misparse "27'-32' S" as if 32
# were seconds of a single point.
_MARK = r'[^\dA-Za-z\s\.,;\-–—‐]{0,3}'
# The degree glyph itself sometimes fails to round-trip through a PDF's
# font encoding and comes out as a literal "(cid:N)" placeholder (pdfminer's
# fallback when a glyph has no ToUnicode mapping) instead of °/º.
_DEG = r'(?:[°º]|\(cid:\d+\))'

_COORD_TOKEN_RE = re.compile(
    rf'({_NUM})\s*{_DEG}\s*(?:({_NUM})\s*{_MARK}\s*(?:({_NUM})\s*{_MARK}\s*)?)?'
    rf'\s*([NSEWnsew])(?![a-zA-Z])'
)
# Sites are sometimes given as a small extent rather than a point, e.g.
# "15°27'-32' S" (latitude spans 15°27'S to 15°32'S). Same idea as above,
# but the minutes are a range; we take the midpoint as the representative
# value.
_COORD_RANGE_TOKEN_RE = re.compile(
    rf'({_NUM})\s*{_DEG}\s*({_NUM}){_MARK}\s*[-–—]\s*({_NUM}){_MARK}'
    rf'\s*([NSEWnsew])(?![a-zA-Z])'
)
# Fallback: plain decimal degree pairs, e.g. "16.319, -46.807" or "0.121497,-78.629857"
_DECIMAL_PAIR_RE = re.compile(
    r'(-?\d{1,2}\.\d{3,})\s*[°]?\s*[NSns]?\s*[,;]\s*(-?\d{1,3}\.\d{3,})\s*[°]?\s*[EWew]?'
)


def _dms_to_decimal(deg, minute, second, hemi):
    val = float(deg) + (float(minute) if minute else 0) / 60 \
        + (float(second) if second else 0) / 3600
    if hemi.upper() in ('S', 'W'):
        val = -val
    return round(val, 6)


def find_coordinates(text):
    """
    Finds DMS coordinate tokens and pairs adjacent lat (N/S) + lon (E/W)
    tokens that occur within ~200 characters of each other (i.e. clearly
    part of the same locality mention, e.g. "16 19 S, 046 48 E").
    Returns a list of dicts: lat, lon, raw, start, end.
    """
    tokens = []
    for m in _COORD_TOKEN_RE.finditer(text):
        deg, minute, second, hemi = m.groups()
        tokens.append({
            'value': _dms_to_decimal(deg, minute, second, hemi),
            'axis': 'lat' if hemi.upper() in 'NS' else 'lon',
            'raw': m.group(0), 'start': m.start(), 'end': m.end(),
        })
    claimed = [(t['start'], t['end']) for t in tokens]
    for m in _COORD_RANGE_TOKEN_RE.finditer(text):
        if any(s <= m.start() < e for s, e in claimed):
            continue
        deg, min1, min2, hemi = m.groups()
        mid_minute = (float(min1) + float(min2)) / 2
        tokens.append({
            'value': _dms_to_decimal(deg, mid_minute, None, hemi),
            'axis': 'lat' if hemi.upper() in 'NS' else 'lon',
            'raw': m.group(0), 'start': m.start(), 'end': m.end(),
        })
    tokens.sort(key=lambda t: t['start'])

    pairs = []
    pending = {}
    for tok in tokens:
        pending[tok['axis']] = tok
        if 'lat' in pending and 'lon' in pending:
            lat_tok, lon_tok = pending['lat'], pending['lon']
            if abs(lat_tok['start'] - lon_tok['start']) < 200:
                pairs.append({
                    'lat': lat_tok['value'], 'lon': lon_tok['value'],
                    'raw': text[min(lat_tok['start'], lon_tok['start']):
                                max(lat_tok['end'], lon_tok['end'])],
                    'start': min(lat_tok['start'], lon_tok['start']),
                    'end': max(lat_tok['end'], lon_tok['end']),
                })
            pending = {}

    if not pairs:
        for m in _DECIMAL_PAIR_RE.finditer(text):
            lat, lon = float(m.group(1)), float(m.group(2))
            if -90 <= lat <= 90 and -180 <= lon <= 180:
                pairs.append({
                    'lat': lat, 'lon': lon, 'raw': m.group(0),
                    'start': m.start(), 'end': m.end(),
                })
    return pairs


# ── ELEVATION PARSING ────────────────────────────────────────────────────────

ELEV_RANGE_RE = re.compile(
    r'(?:elevation|altitud|altitude|alt\.)s?\s*(?:of|range)?\s*[:\-]?\s*'
    r'(\d{2,5})\s*(?:[-–—]|to)\s*(\d{2,5})\s*m\b', re.I)
ELEV_POINT_RE = re.compile(r'(\d{2,5})\s*m\s*alt\.', re.I)
ELEV_SINGLE_RE = re.compile(
    r'(?:elevation|altitud|altitude|alt\.)\s*(?:of)?\s*[:\-]?\s*(\d{2,5})\s*m\b', re.I)


def find_elevation(text):
    """Best elevation guess in `text`: (value_m, raw) or (None, None)."""
    m = ELEV_RANGE_RE.search(text)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return round((lo + hi) / 2, 1), m.group(0)
    m = ELEV_POINT_RE.search(text)
    if m:
        return float(m.group(1)), m.group(0)
    m = ELEV_SINGLE_RE.search(text)
    if m:
        return float(m.group(1)), m.group(0)
    return None, None


# ── DATE PARSING ─────────────────────────────────────────────────────────────

def _excluded(text, start, end, window=40):
    ctx = text[max(0, start - window):end + window]
    if DATE_EXCLUDE_RE.search(ctx):
        return True
    before = text[max(0, start - 25):start]
    return bool(NOMENCLATURE_CITATION_RE.search(before))


def find_dates(text):
    """
    Returns a list of date candidates, each a dict:
        iso        : best-effort "YYYY-MM-DD" string
        precision  : 'day' | 'month' | 'year_range' | 'year'
        raw        : matched substring
        start, end : character offsets in `text`
        range_note : human-readable range text, set when precision < 'day'
    Filters out citation/publication-metadata dates (Received/Accepted/...).
    """
    out = []

    for m in DATE_DMY_ROMAN_RE.finditer(text):
        day, roman, year = m.groups()
        month = ROMAN_MONTHS.get(roman.upper())
        if month and not _excluded(text, m.start(), m.end()):
            try:
                iso = date(int(year), month, int(day)).isoformat()
            except ValueError:
                continue
            out.append({'iso': iso, 'precision': 'day', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(), 'range_note': None})

    for rx in (DATE_DMY_EN_RE,):
        for m in rx.finditer(text):
            day, mon, year = m.groups()
            month = MONTH_EN.get(mon.lower())
            if month and not _excluded(text, m.start(), m.end()):
                try:
                    iso = date(int(year), month, int(day)).isoformat()
                except ValueError:
                    continue
                out.append({'iso': iso, 'precision': 'day', 'raw': m.group(0),
                            'start': m.start(), 'end': m.end(), 'range_note': None})

    for m in DATE_DMY_ES_RE.finditer(text):
        day, mon, year = m.groups()
        month = MONTH_ES.get(mon.lower())
        if month and not _excluded(text, m.start(), m.end()):
            try:
                iso = date(int(year), month, int(day)).isoformat()
            except ValueError:
                continue
            out.append({'iso': iso, 'precision': 'day', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(), 'range_note': None})

    for m in DATE_MY_EN_RE.finditer(text):
        mon, year = m.groups()
        month = MONTH_EN.get(mon.lower())
        if month and not _excluded(text, m.start(), m.end()):
            iso = date(int(year), month, 1).isoformat()
            out.append({'iso': iso, 'precision': 'month', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(),
                        'range_note': m.group(0)})

    for m in DATE_MY_ES_RE.finditer(text):
        mon, year = m.groups()
        month = MONTH_ES.get(mon.lower())
        if month and not _excluded(text, m.start(), m.end()):
            iso = date(int(year), month, 1).isoformat()
            out.append({'iso': iso, 'precision': 'month', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(),
                        'range_note': m.group(0)})

    for m in YEAR_RANGE_RE.finditer(text):
        y1, y2 = m.group(1), m.group(2)
        if len(y2) == 2:
            y2 = y1[:2] + y2
        if not _excluded(text, m.start(), m.end()):
            iso = date(int(y1), 1, 1).isoformat()
            out.append({'iso': iso, 'precision': 'year_range', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(),
                        'range_note': f"{y1}-{y2}"})

    claimed = [(o['start'], o['end']) for o in out]
    for m in YEAR_RE.finditer(text):
        if any(s <= m.start() < e for s, e in claimed):
            continue
        if not _excluded(text, m.start(), m.end()):
            iso = date(int(m.group(1)), 1, 1).isoformat()
            out.append({'iso': iso, 'precision': 'year', 'raw': m.group(0),
                        'start': m.start(), 'end': m.end(),
                        'range_note': m.group(0)})

    return out


def best_date(candidates):
    """Highest-precision candidate (ties broken by first occurrence)."""
    if not candidates:
        return None
    return sorted(candidates, key=lambda d: (-PRECISION_RANK[d['precision']], d['start']))[0]


# ── SPECIES NAME PARSING ─────────────────────────────────────────────────────

BINOMIAL_RE = re.compile(r'\b([A-Z][a-zà-ÿ]{2,})\s+([a-zà-ÿ][a-zà-ÿ\-]{2,})\b')


def _is_plausible_binomial(genus, epithet):
    """Stopword + length sanity check shared by both species scanners."""
    g, e = genus.lower(), epithet.lower().rstrip('-')
    if g in GENUS_STOPWORDS:
        return False
    if e in GENUS_STOPWORDS or e in SPECIES_EPITHET_STOPWORDS:
        return False
    if len(genus) > MAX_GENUS_LEN or len(epithet) > MAX_EPITHET_LEN:
        return False
    return True


def find_species_candidates(text):
    """
    Returns a list of dicts: genus, species, raw, start, end, has_author.
    `has_author` is True when the binomial is immediately followed by an
    author-citation-looking tail (e.g. "Rchb.f.", "(L.) Rusby"), which is
    a strong signal this is a real taxonomic name and not running prose.
    """
    out = []
    for m in BINOMIAL_RE.finditer(text):
        genus, epithet = m.groups()
        if not _is_plausible_binomial(genus, epithet):
            continue
        tail = text[m.end():m.end() + 30]
        has_author = _has_author_tail(tail)
        out.append({
            'genus': genus, 'species': epithet, 'raw': m.group(0),
            'start': m.start(), 'end': m.end(), 'has_author': has_author,
        })
    return out


def species_heading_candidates(text):
    """
    Finds binomials that look like a taxonomic-treatment section heading:
    "Genus species Author(s), Citation" starting a line, as used by
    Check List / PhytoKeys / Novon-style papers. Highest-confidence
    species signal available, and anchors per-record coordinate/date
    extraction (see `find_observation_records`).
    """
    out = []
    for line_match in re.finditer(r'^[ \t]*(.+)$', text, re.M):
        line = line_match.group(1).strip()
        if not (3 < len(line) < 140):
            continue
        m = BINOMIAL_RE.match(line)
        if not m:
            continue
        genus, epithet = m.groups()
        if not _is_plausible_binomial(genus, epithet):
            continue
        tail = line[m.end():m.end() + 40]
        if _has_author_tail(tail):
            out.append({
                'genus': genus, 'species': epithet,
                'start': line_match.start(1), 'end': line_match.end(1),
                'raw': line,
            })
    return out


# ── "Observations." / "Material examined." per-record citation blocks ──────
# Common in Darwin-Core-flavoured taxonomic treatments (e.g. Check List
# journal): "Observations. COUNTRY — Region • Locality; lat, lon; elev m
# alt.; date; collector obs." Gives a genuine PER-RECORD lat/lon + date.

OBS_RECORD_RE = re.compile(
    r'(?:Observations?|Material\s+examined)\.\s*'
    r'(?P<locality>[A-Za-zÀ-ÿ0-9°′″\'\".,;:•\-–—()\s]{0,300}?)'
    r'(?P<lat>\d{1,3}\s*[°º][^,;]{0,15}?[NSns])\s*,?\s*'
    r'(?P<lon>\d{1,3}\s*[°º][^,;]{0,15}?[EWew])\s*;\s*'
    r'(?P<elev>\d{2,5})\s*m\s*alt\.?\s*;\s*'
    r'(?P<date>\d{1,2}\.[IVXLCDM]{1,4}\.\d{4}|\d{1,2}\s+\w+\.?\s+\d{4})',
    re.S,
)


def find_observation_records(text):
    """
    Returns a list of dicts: locality, lat, lon, elev_m, date (raw), start,
    end. Falls back gracefully (returns []) if the paper doesn't use this
    citation style -- most don't, that's expected.
    """
    out = []
    for m in OBS_RECORD_RE.finditer(text):
        coord_txt = f"{m.group('lat')}, {m.group('lon')}"
        pairs = find_coordinates(coord_txt)
        if not pairs:
            continue
        out.append({
            'locality': re.sub(r'\s+', ' ', m.group('locality')).strip(' •-'),
            'lat': pairs[0]['lat'], 'lon': pairs[0]['lon'],
            'elev_m': float(m.group('elev')),
            'date_raw': m.group('date'),
            'start': m.start(), 'end': m.end(),
        })
    return out


# ── HEIGHT / JOHANSSON-ZONE / AGH-RANGE PARSING (free text) ────────────────

HEIGHT_RANGE_RE = re.compile(
    r'(\d{1,3}(?:\.\d+)?)\s*[-–—]\s*(\d{1,3}(?:\.\d+)?)\s*m\b'
    r'(?!\s*(?:alt|elevation|asl))', re.I)
HEIGHT_POINT_RE = re.compile(
    r'(?:at|to|of)\s+(?:ca\.?\s*)?(\d{1,3}(?:\.\d+)?)\s*m\s*(?:high|tall|'
    r'above[\s-]ground)?', re.I)
AGH_LABEL_RE = re.compile(r'\bAGH\b', re.I)
JZ_LABEL_RE = re.compile(r'\bJZ\s*[1-5]\b|\bZ[IVX]{1,3}\b|Johansson\s+zone', re.I)
# Some papers give vertical position only as a qualitative classification
# in running prose rather than meters or a table (e.g. "trunk-restricted",
# "shade epiphyte", "crown-centered"). No Height_m follows from this, but
# it does satisfy the height/position requirement per our schema.
POSITION_CATEGORY_RE = re.compile(
    r'\b(?:trunk|crown|canopy|base)[\s-]restricted\b|'
    r'\bcrown[\s-]cent(?:er|re)d\b|'
    r'\b(?:trunk|canopy|crown)\s+epiphyte(?:s)?\b|'
    r'\b(?:shade|sun)\s+epiphyte(?:s)?\b|'
    r'\bgeneralist(?:\s+epiphyte(?:s)?)?\b|'
    r'\bwhole[\s-]profile\b|'
    r'\brestricted\s+to\s+the\s+trunk\s+base\b|'
    r'\bevenly\s+distributed\b', re.I)


# A height mention near these phrases describes the *forest/tree* canopy
# at the site level (like CanopyHeight_m), not this species' own
# above-ground position -- e.g. "mean canopy height ... ranged from 9 to
# 16 m" is a site statistic that would otherwise get misattributed as a
# co-occurring epiphyte's Height_m purely because they share a block.
CANOPY_CONTEXT_RE = re.compile(
    r'canopy\s+height|mean\s+canopy|canopy\s+cover|tree\s+height|'
    r'forest\s+canopy|phorophyte\s+height', re.I)


def find_height_evidence(text):
    """
    Best-effort free-text height evidence in `text`:
    returns dict(height_m, height_min, height_max, basis, raw) or None.
    """
    for m in HEIGHT_RANGE_RE.finditer(text):
        if CANOPY_CONTEXT_RE.search(text[max(0, m.start() - 150):m.start()]):
            continue
        lo, hi = float(m.group(1)), float(m.group(2))
        if lo < 100 and hi < 100:  # guard against catching elevation ranges
            basis = 'AGH_range' if AGH_LABEL_RE.search(text) else 'height_range'
            return {'height_m': round((lo + hi) / 2, 2), 'height_min': lo,
                    'height_max': hi, 'basis': basis, 'raw': m.group(0)}
    for m in HEIGHT_POINT_RE.finditer(text):
        if CANOPY_CONTEXT_RE.search(text[max(0, m.start() - 150):m.start()]):
            continue
        val = float(m.group(1))
        if val < 100:
            return {'height_m': val, 'height_min': None, 'height_max': None,
                    'basis': 'point_estimate', 'raw': m.group(0)}
    return None


# ── TABLE PARSING (Johansson-zone / AGH / abundance columns) ────────────────

ZONE_HEADER_RE = re.compile(r'\bJZ\s*[1-5]\b|\bZ[IVX]{1,3}\b', re.I)
ROMAN_ZONE_SET = {"I", "II", "III", "IV", "V"}
HEIGHT_HEADER_RE = re.compile(r'\bAGH\b|height\s*\(m\)|altura\s*\(m\)', re.I)
ABUND_HEADER_RE = re.compile(
    r'abundance|individuals|no\.?\s*of\s*ind|n[uú]mero\s+de\s+individuos|'
    r'\brec\.?\b', re.I)
# A coarser (non-metric) vertical-position system some papers use instead of
# meters or Johansson zones: presence/absence on trunk vs. canopy, or a
# categorical classification (shade/sun/generalist epiphyte). Doesn't give
# a Height_m value, but it *is* height/position evidence per our schema.
POSITION_HEADER_RE = re.compile(
    r'\btrunk\b|\bcanopy\b|present\s*(?:on|in)|absent\s*(?:on|in)', re.I)


def _clean_cell(cell):
    return (cell or '').replace('\n', ' ').strip()


# A numeric range without a required unit suffix, for use *only* on a cell
# already known (via its column header) to be a height/AGH column.
NUMERIC_RANGE_RE = re.compile(r'(\d{1,3}(?:\.\d+)?)\s*[-–—]\s*(\d{1,3}(?:\.\d+)?)')
NUMERIC_RE = re.compile(r'\d{1,3}(?:\.\d+)?')


def _header_score(row):
    """How strongly a row looks like a table header for our purposes."""
    score = 0
    for cell in row:
        c = _clean_cell(cell)
        if not c:
            continue
        if ZONE_HEADER_RE.search(c) or c in ROMAN_ZONE_SET:
            score += 1
        if HEIGHT_HEADER_RE.search(c):
            score += 1
        if ABUND_HEADER_RE.search(c):
            score += 1
        if POSITION_HEADER_RE.search(c):
            score += 1
    return score


def parse_zone_tables(tables):
    """
    tables: list of (page_num, table) where table is pdfplumber's
    `.extract_tables()` row-list-of-cells output.
    Returns dict[(genus_lower, species_lower)] -> evidence dict with
    zone_counts (dict label->int), abundance, height_min/max, raw, page.
    Best-effort: pdfplumber's table detection often mangles the ragged,
    two-column layouts these papers use, so an empty result here is
    common and expected -- free-text scanning is the primary fallback.
    """
    evidence = {}
    for page_num, table in tables:
        if not table or len(table) < 2:
            continue

        header_idx, header_score = None, 0
        for i, row in enumerate(table[:3]):
            score = _header_score(row)
            if score > header_score:
                header_idx, header_score = i, score
        if header_idx is None:
            continue  # no zone/height/abundance header found anywhere

        col_labels, height_col, abund_col = {}, None, None
        has_position_header = False
        for j, cell in enumerate(table[header_idx]):
            c = _clean_cell(cell)
            if not c:
                continue
            if c in ROMAN_ZONE_SET:
                col_labels[j] = f"Z{c}"
            elif ZONE_HEADER_RE.search(c):
                col_labels[j] = c.upper().replace(' ', '')
            elif HEIGHT_HEADER_RE.search(c):
                height_col = j
            elif ABUND_HEADER_RE.search(c):
                abund_col = j
            if POSITION_HEADER_RE.search(c):
                has_position_header = True

        for row in table[header_idx + 1:]:
            row_text = ' '.join(_clean_cell(c) for c in row)
            sp_hits = find_species_candidates(row_text)
            if not sp_hits:
                continue
            genus, species = sp_hits[0]['genus'], sp_hits[0]['species']
            key = (genus.lower(), species.lower())
            ent = evidence.setdefault(key, {
                'genus': genus, 'species': species, 'zone_counts': {},
                'height_min': None, 'height_max': None, 'abundance': None,
                'raw': [], 'page': page_num, 'position_evidence': False,
            })
            ent['raw'].append(row_text)
            if has_position_header:
                ent['position_evidence'] = True

            for j, label in col_labels.items():
                if j < len(row):
                    val = _clean_cell(row[j])
                    if re.fullmatch(r'\d+', val):
                        ent['zone_counts'][label] = int(val)

            if height_col is not None and height_col < len(row):
                cell = _clean_cell(row[height_col])
                hrange = NUMERIC_RANGE_RE.search(cell)
                if hrange:
                    lo, hi = float(hrange.group(1)), float(hrange.group(2))
                    ent['height_min'], ent['height_max'] = lo, hi
                else:
                    hnum = NUMERIC_RE.search(cell)
                    if hnum:
                        ent['height_min'] = ent['height_max'] = float(hnum.group(0))

            if abund_col is not None and abund_col < len(row):
                cell = _clean_cell(row[abund_col])
                anum = NUMERIC_RE.search(cell)
                if anum:
                    ent['abundance'] = int(float(anum.group(0)))
    return evidence


# ── ROTATED-TABLE RECOVERY ───────────────────────────────────────────────────
#
# Wide tables are often typeset sideways (rotated ~90°) to fit a portrait
# page -- a common layout for a "Appendix: per-species counts by zone/site"
# table. pdfplumber's default extract_text()/extract_tables() assume
# upright text: they sort characters by ascending (top, then x0), which for
# a 90°-rotated run produces every line *character-reversed* ("Appendix"
# comes out as "xidneppA") and extract_tables() finds no table structure at
# all (0 rows). Both a naive keyword search and the regex scanners above
# are therefore completely blind to this content, even though it's fully
# present as real (searchable, non-OCR) text in the PDF.
#
# The fix works directly on pdfplumber's per-character data (`page.chars`),
# which carries each character's true Unicode identity untouched by
# rotation (only its drawing position/orientation is rotated) plus an
# `upright` flag pdfminer sets per character:
#   1. Detect a rotated page: mostly-non-upright characters.
#   2. For a 90°-rotated run, all characters of one original line/row share
#      (almost) the same x0 -- rotation turned "left-to-right along a row"
#      into "top-to-bottom along a column". So cluster characters into
#      x0-bands; each band is one original table row, and bands sorted by
#      ascending x0 recover the original top-to-bottom row order (verified
#      against known content: family names come out in alphabetical order).
#   3. Within a band, sort characters by *descending* top to get correct
#      reading order (this, and only this, needs reversing relative to
#      pdfplumber's default ascending-top sort).
#   4. Reconstruct whitespace from the gap between characters. The naive
#      gap (consecutive 'top' values) is contaminated by each character's
#      own width -- for rotated text, a character's rendered "top..bottom"
#      span *is* its original (horizontal) advance width, so e.g. a wide
#      'm' looks like a gap even with zero actual whitespace. The fix:
#      gap = previous_char['top'] - next_char['bottom'], which cancels out
#      each character's own width and isolates true inter-glyph whitespace
#      (verified: this gap is ~0 within a word, a few points for a real
#      space, and 15+ points at a table-cell boundary).
#
# This recovers the table with high fidelity (validated against two real
# papers with this exact layout) but is still heuristic, so reconstructed
# pages are also handed to the normal free-text scanners (species/height/
# date) as a safety net, and the raw reconstructed row is always kept in
# `_raw_evidence` for a human to double check.

ROTATED_FRACTION_THRESHOLD = 0.5   # share of non-upright chars to call a page "rotated"
ROTATED_MIN_CHARS = 20
BAND_X_TOLERANCE = 1.5             # points; characters within this share a row-band
CELL_DELIM = '\t'                  # field separator used in synthetic reconstructed tables


def _page_is_rotated(page):
    chars = page.chars
    if len(chars) < ROTATED_MIN_CHARS:
        return False
    n_non_upright = sum(1 for c in chars if not c.get('upright', True))
    return (n_non_upright / len(chars)) >= ROTATED_FRACTION_THRESHOLD


def _cluster_row_bands(chars, x_tol=BAND_X_TOLERANCE):
    """Groups rotated characters into x0-bands (each = one original row)."""
    chars_sorted = sorted(chars, key=lambda c: c['x0'])
    bands = []
    for c in chars_sorted:
        if bands and abs(c['x0'] - bands[-1]['x0']) <= x_tol:
            bands[-1]['chars'].append(c)
        else:
            bands.append({'x0': c['x0'], 'chars': [c]})
    return bands


def _reconstruct_row(band_chars):
    """
    Returns (row_text, row_fields) for one rotated row-band, with
    whitespace and field breaks reinserted (see module-level comment for
    the gap math). row_fields is None if no clear field boundary was found
    (i.e. this band is just a caption/title line, not a tabular row).
    """
    cs = sorted(band_chars, key=lambda c: -c['top'])
    if len(cs) == 1:
        return cs[0]['text'], None
    avg_size = sum(c['bottom'] - c['top'] for c in cs) / len(cs)
    word_cutoff = max(0.3 * avg_size, 1.0)
    cell_cutoff = max(4.5 * avg_size, 15.0)

    parts = [cs[0]['text']]
    fields = [[cs[0]['text']]]
    saw_cell_break = False
    for i in range(1, len(cs)):
        gap = cs[i - 1]['top'] - cs[i]['bottom']
        if gap >= cell_cutoff:
            parts.append(CELL_DELIM)
            fields.append([])
            saw_cell_break = True
        elif gap >= word_cutoff:
            parts.append(' ')
            fields[-1].append(' ')
        parts.append(cs[i]['text'])
        fields[-1].append(cs[i]['text'])
    row_text = ''.join(parts)
    row_fields = [''.join(f) for f in fields] if saw_cell_break else None
    return row_text, row_fields


def _reconstruct_rotated_page(page, page_num):
    """
    Returns (text, table_or_None) for a rotated page: `text` is the
    de-rotated page content (one reconstructed row per line, in original
    top-to-bottom order) suitable for the normal free-text scanners;
    `table` is a synthetic row-list (like pdfplumber's extract_tables()
    output) built only from bands that showed a clear cell boundary, for
    parse_zone_tables() to consume -- or None if nothing looked tabular.
    """
    chars = [c for c in page.chars if not c.get('upright', True)]
    bands = _cluster_row_bands(chars)
    lines, table_rows = [], []
    for b in bands:
        row_text, row_fields = _reconstruct_row(b['chars'])
        lines.append(row_text)
        if row_fields:
            table_rows.append(row_fields)
    text = '\n'.join(lines)
    table = table_rows if len(table_rows) >= 2 else None
    return text, table


# ── PDF LOADING ──────────────────────────────────────────────────────────────

def _strip_running_headers_footers(pages, min_page_frac=0.3, min_pages=3):
    """
    Repository download watermarks ("Downloaded from ... by Universitats
    und Landesbibliothek ...", running page headers repeating the paper
    title, etc.) show up on nearly every page and are indistinguishable
    from real content to the regex scanners -- they've produced fake
    "species" rows in practice. Detect them generically: any line (after
    normalizing digits, which absorb page numbers/dates) that recurs on a
    large fraction of pages is boilerplate, not paper content, and is
    stripped from every page before further processing.
    """
    n_pages = len(pages)
    if n_pages < min_pages:
        return pages
    norm_counts = {}
    per_page_lines = []
    for _, text in pages:
        lines = text.split('\n')
        per_page_lines.append(lines)
        seen = set()
        for line in lines:
            stripped = line.strip()
            if len(stripped) < 15:
                continue
            norm = re.sub(r'\d+', '#', stripped)
            if norm not in seen:
                norm_counts[norm] = norm_counts.get(norm, 0) + 1
                seen.add(norm)
    threshold = max(min_pages, int(n_pages * min_page_frac))
    boilerplate = {norm for norm, cnt in norm_counts.items() if cnt >= threshold}
    if not boilerplate:
        return pages
    cleaned = []
    for (pnum, _text), lines in zip(pages, per_page_lines):
        kept = [l for l in lines if re.sub(r'\d+', '#', l.strip()) not in boilerplate]
        cleaned.append((pnum, '\n'.join(kept)))
    return cleaned


def load_pdf(pdf_path):
    """Returns (pages, tables): pages=[(num,text)], tables=[(num,table)]."""
    if pdfplumber is None:
        sys.exit(
            "pdfplumber is required: pip install pdfplumber pandas"
        )
    pages, tables, rotated_pages = [], [], []
    with pdfplumber.open(pdf_path) as pdf:
        for i, page in enumerate(pdf.pages, start=1):
            if _page_is_rotated(page):
                text, table = _reconstruct_rotated_page(page, i)
                rotated_pages.append(i)
                if table:
                    tables.append((i, table))
            else:
                text = page.extract_text() or ''
                for table in (page.extract_tables() or []):
                    tables.append((i, table))
            pages.append((i, text))
    if rotated_pages:
        print(f"  [i] recovered {len(rotated_pages)} rotated/sideways page(s): "
              f"{rotated_pages} (de-rotated via character-position reconstruction)")
    pages = _strip_running_headers_footers(pages)
    return pages, tables


# ── DOCUMENT-LEVEL METADATA GUESSES ─────────────────────────────────────────

CITE_AS_RE = re.compile(
    r'(?:Citar\s+como|Cite\s+as|How\s+to\s+cite)\s*[:.]?\s*(.+)', re.I)
FIRST_AUTHOR_YEAR_RE = re.compile(
    r'\b([A-ZÀ-Ý][a-zà-ÿ\-]+(?:[\s,]+[A-Z]\.){1,3})[,.]?\s*(?:et al\.?)?\s*'
    r'\(?((?:19|20)\d{2})\)?')


def guess_source(first_page_text):
    m = CITE_AS_RE.search(first_page_text)
    if m:
        snippet = re.sub(r'\s+', ' ', m.group(1)).strip()
        ym = re.search(r'(19|20)\d{2}', snippet)
        m2 = re.match(r'([A-ZÀ-Ý][a-zà-ÿ\-]+)', snippet)
        if m2 and ym:
            return f"{m2.group(1)} et al. {ym.group(0)}"
    m = FIRST_AUTHOR_YEAR_RE.search(first_page_text[:2000])
    if m:
        surname = m.group(1).split(',')[0].split()[0]
        return f"{surname} et al. {m.group(2)}"
    return None


def guess_fallback_date(all_text):
    """
    Best document-level date to use when a species has no local date
    evidence: prefer a date near fieldwork-keywords, else the
    highest-precision date found anywhere (excluding publication dates,
    already filtered out by find_dates).
    """
    candidates = find_dates(all_text)
    if not candidates:
        return None
    near_fieldwork = [
        c for c in candidates
        if FIELDWORK_KEYWORDS_RE.search(
            all_text[max(0, c['start'] - 150):c['start']])
    ]
    pool = near_fieldwork or candidates
    return best_date(pool)


# ── ROW ASSEMBLY ─────────────────────────────────────────────────────────────

# "Host taxa. Genus species (Family), ..." lists the phorophyte(s) of the
# epiphyte a treatment section is actually about -- not new observations.
HOST_BLOCK_RE = re.compile(r'^\s*Host\s+tax(?:a|on)\b', re.I)
# Darwin-Core-style treatment sections often glue several labeled
# sub-parts together with no blank line between them (e.g. a species
# heading immediately followed by "Remarks. ... Host taxa. ..." as one
# block). Catching only whole-block-starts-with-"Host taxa" (above)
# misses these; this finds each label's position *within* a block so we
# can tell which labeled sub-part a given match actually falls under.
SECTION_LABEL_RE = re.compile(
    r'(?:^|\n)[ \t]*(Identification|Remarks|Host\s+tax(?:a|on)|Distribution|'
    r'Material\s+examined|Observations?)\s*\.', re.I)


def _preceding_section_label(block, pos):
    """Lowercased label of the sub-section `pos` falls in, or None."""
    label = None
    for m in SECTION_LABEL_RE.finditer(block):
        if m.start() > pos:
            break
        label = m.group(1).lower()
    return label


def _split_blocks(text):
    """Splits page text into paragraph-ish blocks on blank lines."""
    blocks, cur, start = [], [], 0
    pos = 0
    for line in text.split('\n'):
        if line.strip() == '':
            if cur:
                blocks.append((start, '\n'.join(cur)))
                cur = []
            pos += len(line) + 1
            start = pos
            continue
        cur.append(line)
        pos += len(line) + 1
    if cur:
        blocks.append((start, '\n'.join(cur)))
    return blocks


def _collect_host_species(pages):
    """
    (genus_lower, species_lower) pairs named in a "Host taxa."/"Host
    taxon." labeled sub-section anywhere in the document. These same
    host trees are often *also* named elsewhere with a fully valid
    author citation of their own (in a host-species summary table, or in
    running prose praising a "diverse host tree") -- genuinely
    well-formed taxonomic names, just not the epiphyte a treatment is
    about. Once we've seen a name explicitly labeled as a host, treat
    every other mention of that same name in this document as the same
    host too, not a coincidental new epiphyte.
    """
    hosts = set()
    for _, text in pages:
        for _start, block in _split_blocks(text):
            if not SECTION_LABEL_RE.search(block):
                continue
            for c in find_species_candidates(block):
                label = _preceding_section_label(block, c['start'])
                if label and label.startswith('host'):
                    hosts.add((c['genus'].lower(), c['species'].lower()))
    return hosts


def build_rows(pages, tables, source, site, verbose=False):
    all_text = '\n'.join(t for _, t in pages)
    fallback_coords = find_coordinates(all_text)
    fallback_coord = fallback_coords[0] if fallback_coords else None
    fallback_elev, _ = find_elevation(all_text)
    fallback_date = guess_fallback_date(all_text)
    known_hosts = _collect_host_species(pages)

    rows = {}  # (genus_lower, species_lower) -> row dict

    def get_row(genus, species):
        key = (genus.lower(), species.lower())
        if key in known_hosts:
            return None
        if key not in rows:
            rows[key] = {
                'Source': source, 'Area_or_Site': site,
                'lat': None, 'lon': None, 'Elevation_m': fallback_elev,
                'FieldID': 'NA', 'Abundance': None, 'Height_m': None,
                'CanopyHeight_m': None, 'note': [], 'Genus': genus,
                'species': species, 'FinalID': f"{genus} {species}",
                'datetime': None,
                '_page': set(), '_lat_lon_basis': 'site-level (fallback)',
                '_date_basis': 'document-level (fallback)',
                '_height_basis': None, '_raw_evidence': [],
            }
            if fallback_coord:
                rows[key]['lat'] = fallback_coord['lat']
                rows[key]['lon'] = fallback_coord['lon']
            if fallback_date:
                rows[key]['datetime'] = fallback_date['iso']
                if fallback_date['range_note']:
                    rows[key]['note'].append(
                        f"date approx (paper reports {fallback_date['range_note']})")
        return rows[key]

    # 1) Table-driven evidence (Johansson zones / AGH ranges / abundance).
    zone_evidence = parse_zone_tables(tables)
    for (g_lo, s_lo), ent in zone_evidence.items():
        row = get_row(ent['genus'], ent['species'])
        if row is None:
            continue
        row['_page'].add(ent['page'])
        if ent['zone_counts']:
            zc = ', '.join(f"{k}={v}" for k, v in sorted(ent['zone_counts'].items()))
            row['note'].append(f"Johansson zone counts: {zc}")
            row['_height_basis'] = 'johansson_zone_counts (see note)'
            row['Abundance'] = row['Abundance'] or sum(ent['zone_counts'].values())
        if ent['height_min'] is not None:
            row['Height_m'] = round((ent['height_min'] + ent['height_max']) / 2, 2)
            row['note'].append(f"height range {ent['height_min']}-{ent['height_max']} m")
            row['_height_basis'] = 'height_range (table)'
        if ent['abundance'] is not None:
            row['Abundance'] = ent['abundance']
        if (ent.get('position_evidence') and not ent['zone_counts']
                and ent['height_min'] is None and row['_height_basis'] is None):
            row['note'].append(f"vertical-position row (unparsed, verify manually): "
                                f"{ent['raw'][0][:200]}")
            row['_height_basis'] = 'trunk_canopy_position (raw row, see note)'
        row['_raw_evidence'].extend(ent['raw'][:2])

    # 2) Darwin-Core-style "Observations." per-record citations, anchored
    #    to the nearest preceding species-heading line.
    for page_num, text in pages:
        headings = species_heading_candidates(text)
        records = find_observation_records(text)
        for rec in records:
            heading = None
            for h in headings:
                if h['start'] <= rec['start']:
                    heading = h
                else:
                    break
            if heading is None:
                continue
            row = get_row(heading['genus'], heading['species'])
            if row is None:
                continue
            row['_page'].add(page_num)
            row['lat'], row['lon'] = rec['lat'], rec['lon']
            row['_lat_lon_basis'] = 'per-observation record'
            row['Elevation_m'] = rec['elev_m']
            row['Abundance'] = (row['Abundance'] or 0) + 1
            dcands = find_dates(rec['date_raw'])
            d = best_date(dcands)
            if d:
                row['datetime'] = d['iso']
                row['_date_basis'] = 'per-observation record'
            row['_raw_evidence'].append(
                f"Observations.: {rec['locality']} {rec['lat']},{rec['lon']}; "
                f"{rec['elev_m']} m alt.; {rec['date_raw']}")

    # 3) Free-text, block-scoped evidence. Two ways a species earns a row
    #    here: it heads a taxonomic-treatment section ("Genus species
    #    Author, Citation" at the start of a line), or it's named
    #    mid-sentence immediately followed by an author citation (e.g.
    #    "... Anisotes venosus T.F.Daniel, Letsara & Martin-Bravo was
    #    found ..."). Either way we require the author-citation anchor --
    #    a binomial-shaped word pair with no citation attached is too
    #    likely a false positive, or just a bare second mention of a
    #    species already captured via its heading/table/Observations
    #    record elsewhere. Once anchored, nearby coordinates / dates /
    #    height mentions within the same paragraph-ish block are pulled
    #    in as that species' local evidence.
    for page_num, text in pages:
        for start, block in _split_blocks(text):
            if HOST_BLOCK_RE.match(block):
                # "Host taxa. Baudouinia fluggeiformis (Fabaceae), ..." --
                # these are the phorophyte(s) of the epiphyte this section
                # is actually about, not new epiphyte observations.
                continue
            candidates = {}
            for h in species_heading_candidates(block):
                candidates[(h['genus'], h['species'])] = h
            for c in find_species_candidates(block):
                if c['has_author']:
                    candidates.setdefault((c['genus'], c['species']), c)
            if not candidates:
                continue

            for genus, species in candidates:
                cand = candidates[(genus, species)]
                label = _preceding_section_label(block, cand['start'])
                if label and label.startswith('host'):
                    # Named only inside this block's "Host taxa."
                    # sub-part -- it's the phorophyte, not a new
                    # epiphyte observation.
                    continue
                row = get_row(genus, species)
                if row is None:
                    continue
                row['_page'].add(page_num)

                coords = find_coordinates(block)
                if coords and row['_lat_lon_basis'] != 'per-observation record':
                    row['lat'], row['lon'] = coords[0]['lat'], coords[0]['lon']
                    row['_lat_lon_basis'] = 'per-observation (block text)'

                dcands = find_dates(block)
                d = best_date(dcands)
                if d and row['_date_basis'] != 'per-observation record':
                    row['datetime'] = d['iso']
                    row['_date_basis'] = 'per-observation (block text)'
                    if d['range_note']:
                        row['note'].append(f"date approx ({d['range_note']})")

                elev, elev_raw = find_elevation(block)
                if elev is not None:
                    row['Elevation_m'] = elev

                if row['_height_basis'] is None:
                    hev = find_height_evidence(block)
                    pos_match = POSITION_CATEGORY_RE.search(block)
                    if hev:
                        row['Height_m'] = hev['height_m']
                        row['_height_basis'] = hev['basis']
                        row['note'].append(f"height evidence: {hev['raw']}")
                    elif pos_match:
                        row['_height_basis'] = 'categorical_position (see note)'
                        row['note'].append(f"vertical-position category: {pos_match.group(0)}")
                    elif JZ_LABEL_RE.search(block):
                        row['_height_basis'] = 'johansson_zone_mentioned (no counts parsed)'

                row['_raw_evidence'].append(block[:200].replace('\n', ' '))

    return list(rows.values())



# _height_basis values that mean "we know something is nearby but didn't
# actually capture a height/zone value" -- these should NOT count as
# satisfying the height requirement.
HEIGHT_BASIS_UNRESOLVED = {None, 'johansson_zone_mentioned (no counts parsed)'}


def finalize_rows(rows):
    """Computes _confidence / _missing_required and flattens list fields."""
    out = []
    for row in rows:
        missing = []
        if row['lat'] is None or row['lon'] is None:
            missing.append('lat_lon')
        if row['datetime'] is None:
            missing.append('date')
        if row['Height_m'] is None and row['_height_basis'] in HEIGHT_BASIS_UNRESOLVED:
            missing.append('height')
        # species is always present (it's the row key)

        row = dict(row)
        row['_missing_required'] = ', '.join(missing) if missing else ''
        row['_confidence'] = f"{4 - len(missing)}/4"
        row['_page'] = ', '.join(str(p) for p in sorted(row['_page']))
        row['note'] = '; '.join(dict.fromkeys(row['note']))  # dedupe, keep order
        row['_raw_evidence'] = ' || '.join(row['_raw_evidence'][:5])
        if row['_height_basis'] is None:
            row['_height_basis'] = ''
        out.append({c: row.get(c) for c in STAGING_COLUMNS})
    return out


# ── CANDIDATE DUMP (manual cross-checking mode) ─────────────────────────────

def dump_candidates(pages):
    all_text = '\n'.join(t for _, t in pages)
    print("\n=== Coordinate candidates ===")
    for c in find_coordinates(all_text):
        print(f"  lat={c['lat']:.6f} lon={c['lon']:.6f}  raw={c['raw']!r}")
    print("\n=== Date candidates (publication metadata excluded) ===")
    for d in sorted(find_dates(all_text), key=lambda d: d['start']):
        print(f"  {d['iso']}  ({d['precision']})  raw={d['raw']!r}")
    print("\n=== Species-heading candidates (treatment-style papers) ===")
    seen = set()
    for page_num, text in pages:
        for h in species_heading_candidates(text):
            key = (h['genus'], h['species'])
            if key in seen:
                continue
            seen.add(key)
            print(f"  p.{page_num}  {h['genus']} {h['species']}")
    print("\n=== 'Observations.' per-record citations ===")
    for page_num, text in pages:
        for r in find_observation_records(text):
            print(f"  p.{page_num}  {r['lat']:.5f},{r['lon']:.5f}  "
                  f"{r['elev_m']} m alt.  {r['date_raw']}  ({r['locality'][:60]})")


# ── CLI ──────────────────────────────────────────────────────────────────────

def process_pdf(pdf_path, outdir, source, site, dump=False, verbose=False):
    print(f"\n--- {pdf_path} ---")
    pages, tables = load_pdf(pdf_path)
    print(f"  {len(pages)} pages, {len(tables)} tables extracted")

    if dump:
        dump_candidates(pages)
        return None

    stem = os.path.splitext(os.path.basename(pdf_path))[0]
    src = source or guess_source(pages[0][1]) or stem
    if not source:
        print(f"  [!] --source not given, guessed: {src!r} -- verify this")
    st = site or stem
    if not site:
        print(f"  [!] --site not given, using filename stem: {st!r} -- set --site")

    rows = build_rows(pages, tables, src, st, verbose=verbose)
    rows = finalize_rows(rows)
    df = pd.DataFrame(rows, columns=STAGING_COLUMNS)
    df = df.sort_values(['_confidence', 'Genus', 'species'],
                         ascending=[False, True, True])

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{stem}_staging.csv")
    df.to_csv(out_path, index=False, na_rep='')

    n = len(df)
    complete = (df['_missing_required'] == '').sum()
    print(f"  {n} candidate species rows -> {out_path}")
    print(f"  Complete (all 4 required fields): {complete}/{n}")
    if n:
        for field in ('lat_lon', 'date', 'height'):
            n_missing = df['_missing_required'].str.contains(field, na=False).sum()
            print(f"    missing {field:8s}: {n_missing}")
    print("  Review the staged CSV before merging into data/csv/combinedv3.csv.")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Stage candidate epiphyte observations from a literature "
                    "PDF into the data/csv/combinedv3.csv schema, for review.")
    parser.add_argument("pdfs", nargs='+', help="One or more PDF paths.")
    parser.add_argument("--outdir", default="data/literature_staging",
                        help="Directory for staging CSVs "
                             "(default: data/literature_staging)")
    parser.add_argument("--source", default=None,
                        help="Value for the Source column (e.g. 'Komada et al. 2025'). "
                             "If omitted, guessed from the PDF and flagged for review.")
    parser.add_argument("--site", default=None,
                        help="Value for the Area_or_Site column. "
                             "If omitted, the PDF filename is used and flagged for review.")
    parser.add_argument("--dump-candidates", action="store_true",
                        help="Print raw regex hits (coordinates/dates/species) "
                             "instead of writing a staging CSV. Useful for "
                             "sanity-checking extraction on a new paper layout.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    failures = []
    for pdf_path in args.pdfs:
        if not os.path.isfile(pdf_path):
            print(f"  [!] not found, skipping: {pdf_path}")
            failures.append((pdf_path, "not found"))
            continue
        # A batch run over a real literature corpus will include corrupted,
        # encrypted, or scanned-image-only PDFs -- one bad file shouldn't
        # abort the rest of the run.
        try:
            process_pdf(
                pdf_path, args.outdir, args.source, args.site,
                dump=args.dump_candidates, verbose=args.verbose,
            )
        except Exception as exc:
            print(f"  [!] failed to process {pdf_path}: {exc}")
            failures.append((pdf_path, str(exc)))

    if len(args.pdfs) > 1 and failures:
        print(f"\n{len(failures)}/{len(args.pdfs)} file(s) failed:")
        for path, reason in failures:
            print(f"  - {path}: {reason}")


if __name__ == "__main__":
    main()
