
remove solvent

align 1DF6_bondsonly_model, 1DF6-assembly1
align 1DF6_encoded_model, 1DF6-assembly1

set_name 1DF6-assembly1, 1DF6_bondsonly_ref
copy 1DF6_encoded_ref, 1DF6_bondsonly_ref

set grid_mode,1
set grid_slot, 1, 1DF6_encoded_model
set grid_slot, 1, 1DF6_encoded_ref
set grid_slot, 2, 1DF6_bondsonly_model
set grid_slot, 2, 1DF6_bondsonly_ref

color gray
color green, 1DF6_encoded_model
color blue, 1DF6_bondsonly_model

select (resn cys)
show sticks, (sele)

### cut below here and paste into script ###
set_view (\
     0.813714802,    0.579730511,    0.042154476,\
     0.543022990,   -0.732314229,   -0.410901546,\
    -0.207340419,    0.357247800,   -0.910705507,\
     0.000000447,   -0.000003036,  -67.083122253,\
     0.335362315,   -1.565861225,   -0.539414167,\
    47.674209595,   86.492042542,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0
png 1DF6.png, width=4cm, height=2cm, dpi=600, ray=1
