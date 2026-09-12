"""Read a molecule's force constants out of Svehla, NASA TR R-132, Table I(a).

    from tests.constants_integrity.svehla import force_constants
    force_constants(page_text, 'CH4') -> {'sigma': '3.758', 'eps_k': '148.6'}

Svehla (NASA TR R-132, 1962, a US Government work) tabulates the Lennard-Jones force
constants sigma (angstrom) and epsilon/k (K) for several hundred molecules across PDF
pp.22-26. It is the classic source for these constants, and GAS_MOLECULAR_PARAMS is
sourced to it.

The PDF is a 1962 SCAN, and its text layer is damaged in three ways, each measured on the
on-disk file rather than assumed:

  SPACED DIGITS  a number is printed with spaces inside it - "3. 711", "10. 22", "2. 900".
                 A run is joined back up; the value is unchanged.

  TWO COLUMN ORDERS  pp.22, 24 and 25 print  molecule sigma METHOD eps/k METHOD page,
                 while p.23 prints  molecule sigma eps/k METHOD page  - the method column
                 for sigma is hoisted into a separate block further down the page. Both
                 orders are tried, longest-first, and a molecule is read only if one of
                 them fits.

  MANGLED TOKENS  the scan renders several formulae wrongly and consistently: Cl2 as
                 "C_", N2 as "N_", SF6 as "SF_", O2 as "02", C3H8 as "C3Hs", n-C4H10 as
                 "n-C4Hlo". The caller passes the token AS THE SCAN PRINTS IT; this module
                 does not guess what a molecule "should" look like.

Two things this module will NOT do, both learned by getting them wrong:

  It will not read a molecule POSITIONALLY. p.26's concluded block separates the molecule
  names from the numbers entirely, and indexing the Nth number tuple does not track the
  Nth molecule, because OCR-broken entries ("& 977 18 390.2 8 107") fail to match a clean
  tuple and silently shift the alignment. Indexing p.26's 20th tuple for Xenon returns
  2.608 / 10.22 / page 114, which is HELIUM's second determination. Xenon has no token and
  is therefore not readable here at all.

  It will not accept a token that is not UNIQUE on its page. The first parser took the
  first match; a token matching twice would then be read from whichever came first with no
  indication that a choice had been made.

A value is returned as the STRING the page prints, so the caller compares digits, and a
page that fits neither layout raises ValueError rather than yielding a guess.

Note on OCR-corrupt values: this module reports what it parses, it does not sanitise.
Ammonia's eps/k is printed "55& 3" (it is 558.3) and parses as "55". Nothing here special-
cases that - a citation claiming 558.3 will simply fail against the 55 this returns, which
is the correct outcome: the reader refuses to confirm what the scan does not legibly say.
"""
from __future__ import annotations

import re

#: a number, possibly with the scan's spaces inside it
NUM = r'\d+\.\s?\d+|\d+'

#: the pages of Table I(a) in the on-disk PDF
PAGES = (22, 23, 24, 25, 26)


def normalize(text):
    """Collapse the text layer's whitespace to single spaces."""
    return ' '.join(text.split())


def _join(tok):
    """'3. 711' -> '3.711'. The scan spaces digits; the value is unchanged."""
    return tok.replace(' ', '')


def occurrences(page_text, token):
    """How many times `token` stands as its own word on this page."""
    return len(re.findall(_anchor(token), normalize(page_text)))


def _anchor(token):
    # not preceded by a word character, '(' or '-', so 'He' does not match inside
    # 'Helium' and 'C_' does not match inside a longer chlorine formula
    return r'(?<![A-Za-z0-9(\-])' + re.escape(token) + r'\s+'


def force_constants(page_text, token):
    """{'sigma': str, 'eps_k': str} - the digits Table I(a) prints for one molecule.

    `token` is the molecule AS THE SCAN PRINTS IT ('C_' for Cl2, '02' for O2).
    Raises ValueError if the token is absent, is not unique on the page, or is
    followed by no sigma/eps-k pair in either of the table's two column orders.
    """
    s = normalize(page_text)
    n = len(re.findall(_anchor(token), s))
    if n == 0:
        raise ValueError(f'token {token!r} is not on this page')
    if n > 1:
        raise ValueError(f'token {token!r} occurs {n} times on this page; not unique, '
                         f'and a force constant is never taken from the first of several')
    anchor = _anchor(token)
    for order, pat in (
            ('sigma METHOD eps/k', anchor + r'(' + NUM + r')\s+\d+\s+(' + NUM + r')(?![\d.])'),
            ('sigma eps/k', anchor + r'(' + NUM + r')\s+(' + NUM + r')(?![\d.])')):
        m = re.search(pat, s)
        if not m:
            continue
        sigma, eps_k = _join(m.group(1)), _join(m.group(2))
        try:
            float(sigma)
            float(eps_k)
        except ValueError:
            continue
        return {'sigma': sigma, 'eps_k': eps_k, 'order': order}
    raise ValueError(f'no sigma/eps-k pair follows {token!r} in either column order')
