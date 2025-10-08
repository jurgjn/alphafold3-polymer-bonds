
remove solvent

align 6OQ1_bondsonly_model, 6OQ1-assembly1
align 6OQ1_encoded_model, 6OQ1-assembly1

set_name 6OQ1-assembly1, 6OQ1_bondsonly_ref
copy 6OQ1_encoded_ref, 6OQ1_bondsonly_ref

set grid_mode,1
set grid_slot, 1, 6OQ1_encoded_model
set grid_slot, 1, 6OQ1_encoded_ref
set grid_slot, 2, 6OQ1_bondsonly_model
set grid_slot, 2, 6OQ1_bondsonly_ref

color gray
color green, 6OQ1_encoded_model
color blue, 6OQ1_bondsonly_model

### cut below here and paste into script ###
set_view (\
     0.213880911,   -0.784625769,    0.581910551,\
    -0.973274529,   -0.120165512,    0.195698038,\
    -0.083622634,   -0.608217061,   -0.789356589,\
     0.000000000,    0.000000000, -164.296997070,\
    18.210874557,    7.732210159,    6.063646317,\
   132.511642456,  196.082351685,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0
png visualise_6OQ1.png, width=4cm, height=2cm, dpi=600, ray=1
