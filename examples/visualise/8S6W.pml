
remove solvent

align 8S6W_bondsonly_model, 8S6W-assembly1
align 8S6W_encoded_model, 8S6W-assembly1

set_name 8S6W-assembly1, 8S6W_bondsonly_ref
copy 8S6W_encoded_ref, 8S6W_bondsonly_ref

set grid_mode,1
set grid_slot, 1, 8S6W_encoded_model
set grid_slot, 1, 8S6W_encoded_ref
set grid_slot, 2, 8S6W_bondsonly_model
set grid_slot, 2, 8S6W_bondsonly_ref

color gray
color green, 8S6W_encoded_model
color blue, 8S6W_bondsonly_model

### cut below here and paste into script ###
set_view (\
     0.999066532,    0.020660367,   -0.037749916,\
    -0.019877352,    0.999576986,    0.020975810,\
     0.038164970,   -0.020206964,    0.999061823,\
     0.000000000,    0.000000000, -315.108367920,\
   139.736404419,  139.794433594,  136.205902100,\
   254.742858887,  375.473846436,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0
png 8S6W.png, width=4cm, height=2cm, dpi=600, ray=1
