# v8 G5 reading: design re-run not performed

* G1 failed on dev for density, viscosity, flash and freezing (`results/v8/substitution/READING.md`), so no ML-filled constant is admitted.
* G2 failed for freezing on real fuels (`results/v8/represented_margins/READING.md`: no skill over the constant). A design row certifies only if every constraint does, so freezing blocks every P1 to P3 row.
* Tasks 4 to 6 (engine switch, pool P4, design re-run) were therefore skipped: with no admitted ML constant and no held v8 margin they cannot change any certification.
* Certified count stays at v5: 3 of 16. Pool P4 was not formed (no ML constant kind passed G1; no Antoine refit was attempted).
* G3 holds trivially: no v5 code, config or result was changed (`reproduce_base.py` is run in the final check).
