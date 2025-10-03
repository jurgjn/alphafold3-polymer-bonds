remove solvent
align 1df6_baseline_model, 1DF6-assembly1
align 1df6_polybonds_model, 1DF6-assembly1

set_name 1DF6-assembly1, 1df6_baseline_ref
copy 1df6_polybonds_ref, 1df6_baseline_ref

set grid_mode,1
set grid_slot, 1, 1df6_polybonds_model
set grid_slot, 1, 1df6_polybonds_ref
set grid_slot, 2, 1df6_baseline_model
set grid_slot, 2, 1df6_baseline_ref

color gray
color green, 1df6_polybonds_model
color blue, 1df6_baseline_model