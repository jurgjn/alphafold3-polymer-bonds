
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

select (model 6OQ1_bondsonly_model | model 6OQ1_encoded_model) & chain F & (resid 11 | resid 48)
show sticks, (sele)
color orange, (sele)

select (chain AA | chain CA)
show sticks, (sele)
color red, (sele)

select (model 6OQ1_bondsonly_model) & (chain A | chain C) & (resid 76)
show sticks, (sele)
color red, (sele)

### cut below here and paste into script ###
set_view (\
     0.578283727,   -0.257416517,    0.774161041,\
     0.803969860,    0.341055453,   -0.487146914,\
    -0.138631150,    0.904111981,    0.404181808,\
     0.000000000,    0.000000000, -165.200439453,\
    68.374526978,    4.010967255,   14.196598053,\
   132.428390503,  197.972488403,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0
png 6OQ1.png, width=4cm, height=2cm, dpi=600, ray=1
