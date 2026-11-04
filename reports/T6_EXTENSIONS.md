T6 extension decision - implement or exclude per property
model_table rows=42662 splits S1=28604 S2=8002 S3=826

## data inventory per candidate
net_heat_content                   pure=    0 mix=    0 excluded_by_R1=0
gross_heat_of_combustion           pure=    0 mix=    0 excluded_by_R1=0
gcv                                pure=    0 mix=    0 excluded_by_R1=0
mon                                pure=  267 mix=  106 excluded_by_R1=0
dcn                                pure=   80 mix=  346 excluded_by_R1=0
freezing_point_C                   pure=   25 mix=    0 excluded_by_R1=2278
pour_point_C                       pure=    0 mix=    0 excluded_by_R1=60
cloud_point_C                      pure=    0 mix=    0 excluded_by_R1=3

## decisions
- GCV: EXCLUDE - zero pure-compound heats of combustion in open HC data (T0); the physics chain cannot start.
- freezing_point_C: EXCLUDE for now - 25 pure rows < 200-row training floor (same rule raised it in WP4).
- pour_point_C / cloud_point_C: EXCLUDE - exists only as excluded L0 real-fuel rows, no molecular anchor.
- mon / dcn: IMPLEMENT - pure and mixture rows both exist (Leenhouts), so provide null + pure-anchored baselines.

== mon: S1 n=306 S2 n=31 (S2 systems=25), range [-2.0, 120.0]
   null baseline  (S1 median = 81.0): MAE=16.40 RMSE=25.53
   pure anchors: 220 molecules (S1 pure median per inchikey)
   ideal-mix: no scoreable S2 mixture rows (skipped=8)

== dcn: S1 n=342 S2 n=67 (S2 systems=20), range [2.8, 104.4]
   null baseline  (S1 median = 27.8): MAE=17.33 RMSE=23.00
   pure anchors: 62 molecules (S1 pure median per inchikey)
   ideal-mix: no scoreable S2 mixture rows (skipped=56)
