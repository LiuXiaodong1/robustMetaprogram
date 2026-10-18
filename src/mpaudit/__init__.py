"""mpaudit: pilot audit (P1) of NMF-based tumour metaprograms."""
__version__ = "0.1.0"  # before the imports: code_version() below hashes it

from . import decoy as _decoy  # noqa: E402,F401  (its decoy generator is hashed below)
from .data import Keys, read_gmt, write_gmt  # noqa: E402
from .gavish import GavishParams, GavishResult, NMFCache, code_version as _code_version, run_gavish  # noqa: E402

__all__ = ["Keys", "GavishParams", "GavishResult", "NMFCache", "run_gavish", "read_gmt", "write_gmt"]

# The code hashes in the NMF cache keys, taken now from the sources this process imported, not at first use from
# files that may have been edited since (mpaudit.gavish.code_version)
_code_version(decoy_data=True)
