siteFile_thermal=MAR_Atlantis_sites_103097.xyz
vtu_rho=vp_clip.vtu
outPath_thermal=.
thermalmodel=vp
../install/vtk2grav -p $siteFile_thermal -F density -D 3300 -i $vtu_rho -o ${outPath_thermal}/grav_${thermalmodel}.txt