# Using gmt potential module to calculate gravity of a cube
cube_len_x=5
cube_len_y=5
cube_len_z=3
z_top_cube=-3
x0=0
y0=0
z0=`echo ${z_top_cube} ${cube_len_z} | awk '{print $1-$2}'`
range=-15/15/-15/15
dxdy=0.1
gmt gmtgravmag3d -Fsites.xy -C1000 -M+sprism,$cube_len_x/$cube_len_y/$cube_len_z/$z0/$x0/$y0 -L-3 
gmt gmtgravmag3d -R$range -I$dxdy -C1000 -M+sprism,$cube_len_x/$cube_len_y/$cube_len_z/$z0/$x0/$y0 -L-0 -Gcube_grav.grd
gmt grd2cpt  cube_grav.grd -E20 -D > m.cpt
gmt begin cube_grav pdf
    gmt grdimage cube_grav.grd -Cm.cpt -JX12c -Ba 
    gmt colorbar -DJCB+o0/1c -Cm.cpt -Bxaf+l"Gravity of cube" -By+l"mGal"
gmt end 

rm m.cpt gmt.history 