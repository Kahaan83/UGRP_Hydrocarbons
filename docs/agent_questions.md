# Agent Questions — T0.2 Leakage Audit

> File maintained by agent subagent. Append new questions; do not delete existing ones.
> Format: date, task, question, default assumed.

---

## Q1 — 2026-11-04 (T0.2)

**Task:** C2 near-duplicate check.

**Question:** The near-duplicate check (C2) uses a 5% value tolerance and ±2 K temperature window
to compare CRC 647 measured values against S1 training rows. This produces 67,324 hits, almost
entirely in `density_kg_m3`. Are these density coincidences (Jet-A density 740–862 kg/m³ at 234–344 K
overlapping with pure/binary alkane densities from ThermoML) considered acceptable leakage,
or should a stricter test (e.g., require the same compound family) be applied?

**Default assumed:** WARN (not FAIL). CRC fuels have no molecule lists (bulk route),
so inchikey-set overlap is structurally impossible. The density coincidences are physical
reality, not data copying. The check correctly confirms that 0 CRC fuel: system_ids
appear in S1.

---

## Q2 — 2026-11-04 (T0.2)

**Task:** C4 LOSO split discipline.

**Question:** `eval_flash_loso.py` runs LOSO over all 31 flash-mixture systems regardless of
their S1/S2 assignment. 22 of 31 LOSO systems have rows in S1 (the t8 training set).
Should the HGB-residual LOSO metric be restricted to only the 9 systems NOT in S1 to give
a truly independent evaluation, or is the current full-LOSO setup acceptable for the
physics baseline comparison only?

**Default assumed:** The current setup is acceptable for the Le Chatelier physics baseline
(which does not train on any rows). Recommend that any HGB-residual LOSO claim be
accompanied by a caveat that 22/31 systems were in t8 training data.

---

## Q3 — 2026-11-04 (T0.4)

**Task:** Verification and Provenance of the 19 November 3 Antoine Parameter Sets.

**Question:** On November 3, 2026, 19 Antoine parameter sets were added to `data/unified/compounds.csv`
with `antoine_source='webbook'` to achieve 100% physics baseline coverage on mixture flash points.
Independent cross-checking against primary NIST WebBook (`webbook.nist.gov`) records reveals:
1. 11 of the 19 compounds (e.g., mesitylene, 2-methylheptane, 4-methylheptane, 2,5-dimethylhexane)
   have Antoine tables on NIST WebBook (Willingham et al. 1945, Stull 1947, Forziati et al. 1949),
   though slight numeric variations exist between recorded constants and WebBook values.
2. 8 of the 19 compounds (hexylbenzene CAS 1077-16-3, isocetane CAS 13475-82-6, octylbenzene CAS 2189-60-8,
   dodecylbenzene CAS 123-01-3, 2,3,3-trimethylpentane CAS 564-02-3, trans-1,4-dimethylcyclohexane CAS 624-29-3,
   cis-1,4-dimethylcyclohexane CAS 2207-04-7, cis-1,2-dimethylcyclohexane CAS 6876-23-9) have NO Antoine
   table on NIST WebBook. Their constants appear to originate from another handbook (e.g., Yaws, Perry's)
   or group-contribution estimations, yet were tagged `antoine_source='webbook'`.
3. Compound `VKPSKYDESGTTFR-UHFFFAOYSA-N` has SMILES `CC(CC(C)(C)C)CC(C)(C)C` (formula C12H26,
   2,2,4,6,6-pentamethylheptane) but is mapped to CAS 13475-82-6 (isocetane, C16H34).
Should these 8 unverified sets be replaced with evaluated Cox/Wagner fits or explicit Yaws/DIPPR citations,
or should they be removed if strict NIST WebBook provenance is enforced?

**Default assumed:** In `data/antoine_sources.csv`, the 11 compounds with matching WebBook data are
provided with full NIST WebBook citations. The 8 unverified compounds are recorded with citation `unverified`
and flagged with `needs-removal-candidate` in `flags/notes` per task instructions.

---

## Q4 — 2026-11-04 (T0.3)

**Task:** Redistribution Licensing for Extracted Report Data and Credentialed Databases.

**Question:** Several key real-fuel datasets in the pipeline were derived by scraping or transcription:
- AJFTD (FAA ASCENT / UIUC): Requires credentialed login; academic research use only.
- CRC Reports (CRC 647, AV-27-18, AV-11-09): Publicly accessible technical reports, but CRC holds copyright.
- Fortin et al. 2012 / Saldaña et al. 2011: Data tables extracted from ACS journal articles.
- CRC Handbook of Chemistry & Physics: Organic melting points extracted from handbook tables.
While factual measurements are not copyrightable in many jurisdictions, database compilations and PDF
distributions may be protected. Should processed derived CSVs of these tables be distributed in public
repositories or restricted to private internal releases?

**Default assumed:** Marked `redistribution_allowed = unclear` in `docs/data_sources.csv` for all such sources,
with explicit explanatory notes detailing source access conditions. Only datasets with explicit CC0, CC BY,
MIT, or U.S. Federal Government Public Domain status are marked `redistribution_allowed = yes`.

---

## Q5 — 2026-11-04 (T0.2)

**Task:** Reconciliation of Near-Duplicate (C2) and Split Discipline (C4) Verdicts with Acceptance Criteria.

**Question:**
1. Check C2 acceptance criteria specifically requires `(CRC system component sets overlapping >=80% with training blends, same property, T within 2 K, value within 5%)`. Because CRC fuels lack speciated molecular component lists, component overlap is 0% (or <=8.3% against 36-component bulk-route reconstructions), yielding 0 near-duplicates (PASS). When component overlap is omitted entirely, 67,324 physical value coincidences appear in density and viscosity.
2. Check C4 acceptance criteria requires `none of the 31 LOSO systems appear in S1 training rows`. 22 of the 31 flash LOSO systems appear in S1 training rows (65 rows), which violates S1 isolation for neural network models (`t8`), although `eval_flash_loso.py` runs its own internal leave-one-system-out cross-validation across the 31 systems.
How should these verdicts be presented in the PASS/FAIL/UNCHECKED schema?

**Default assumed:** Marked C2 as PASS under the strict component overlap criterion (0 near-duplicates), with full documentation and export of the 67,324 relaxed physical coincidences. Marked C4 as FAIL against the S1-isolation criterion (22 systems / 65 rows affected), with explicit documentation that the Le Chatelier baseline does not train on any data, while HGB and t8 neural evaluations must account for S1 presence.

