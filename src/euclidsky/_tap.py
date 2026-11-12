"""A small synchronous TAP client for the ESA Euclid Science Archive.

Anonymous access is allowed. A query that times out on the server or meets a
network error is retried with backoff; a query the server REJECTS (a parse
error, an unknown column) is raised at once, since retrying it cannot help.
"""

from __future__ import annotations

import csv
import io
import math
import time
import urllib.error
import urllib.parse
import urllib.request

TAP_SYNC = "https://eas.esac.esa.int/tap-server/tap/sync"
USER_AGENT = "fwh-euclid/0.1 (catalogue survey; anonymous)"


class QueryError(RuntimeError):
    """The archive rejected the query (do not retry)."""


class TransientError(RuntimeError):
    """The query failed in a way a retry may fix (timeout, 5xx, network)."""


def _once(adql: str, timeout: float) -> str:
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}
    ).encode()
    req = urllib.request.Request(TAP_SYNC, data=body, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            text = r.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        text = e.read().decode("utf-8", "replace")
        if e.code >= 500 or e.code in (408, 429):
            raise TransientError(f"HTTP {e.code}: {text[:300]}") from e
        if "QUERY_STATUS" in text:
            raise QueryError(_votable_error(text)) from e
        raise TransientError(f"HTTP {e.code}: {text[:300]}") from e
    except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
        raise TransientError(str(e)) from e
    if text.lstrip().startswith("<"):
        msg = _votable_error(text)
        low = msg.lower()
        # 'Unknown table' for a table that exists: the server's schema cache
        # is briefly empty (seen on release day, intermittently, for every
        # table). A retry a few seconds later succeeds.
        if "timeout" in low or "aborted" in low or "unknown table" in low:
            raise TransientError(msg)
        raise QueryError(msg)
    return text


def _votable_error(text: str) -> str:
    import re

    m = re.search(r'name="QUERY_STATUS" value="ERROR">(.*?)</INFO>', text, re.DOTALL)
    return (m.group(1) if m else text[:500]).strip()


def query(adql: str, *, retries: int = 5, timeout: float = 600.0) -> list[dict[str, str]]:
    """Run ``adql`` and return the rows as dicts of strings (CSV cells)."""
    delay = 5.0
    for attempt in range(retries + 1):
        try:
            text = _once(adql, timeout)
            return list(csv.DictReader(io.StringIO(text)))
        except TransientError:
            if attempt == retries:
                raise
            time.sleep(delay)
            delay = min(delay * 2, 300.0)
    raise AssertionError("unreachable")


def num(s: str | None) -> float | None:
    """A CSV cell as a float; empty and NaN become None."""
    if s is None or s.strip() in ("", "--") or s.lower() == "nan":
        return None
    v = float(s)
    return None if math.isnan(v) else v
