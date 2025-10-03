remove solvent
align 6oq1_baseline_model, 6OQ1-assembly1
align 6oq1_polybonds_model, 6OQ1-assembly1

set_name 6OQ1-assembly1, 6oq1_baseline_ref
copy 6oq1_polybonds_ref, 6oq1_baseline_ref

set grid_mode,1
set grid_slot, 1, 6oq1_polybonds_model
set grid_slot, 1, 6oq1_polybonds_ref
set grid_slot, 2, 6oq1_baseline_model
set grid_slot, 2, 6oq1_baseline_ref

color gray
color green, 6oq1_polybonds_model
color blue, 6oq1_baseline_model