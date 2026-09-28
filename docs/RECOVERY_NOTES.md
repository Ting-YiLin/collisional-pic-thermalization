# Source recovery notes

The original C++ source was printed in Appendix B of the 2010 thesis PDF. The repository copy was recovered from that appendix.

Only non-scientific recovery/portability edits were made:

1. PDF page-wrap artifacts inside C++ statements were repaired.
2. A fallback definition for `M_PI` was added for compilers that do not expose it by default.
3. The Windows-only `system("PAUSE")` line was disabled.
4. No collision model, particle update, field solve, weighting formula, statistical routine, or research parameter logic was intentionally changed.

The thesis PDF remains the archival source of truth for the historical listing.
