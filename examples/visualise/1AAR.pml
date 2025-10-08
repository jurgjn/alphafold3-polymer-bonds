
remove solvent

align 1AAR_bondsonly_model, 1AAR-assembly1
align 1AAR_encoded_model, 1AAR-assembly1

set_name 1AAR-assembly1, 1AAR_bondsonly_ref
copy 1AAR_encoded_ref, 1AAR_bondsonly_ref

set grid_mode,1
set grid_slot, 1, 1AAR_encoded_model
set grid_slot, 1, 1AAR_encoded_ref
set grid_slot, 2, 1AAR_bondsonly_model
set grid_slot, 2, 1AAR_bondsonly_ref

color gray
color green, 1AAR_encoded_model
color blue, 1AAR_bondsonly_model

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
png visualise_1AAR.png, width=4cm, height=2cm, dpi=600, ray=1
