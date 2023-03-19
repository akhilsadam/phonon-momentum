echo "Merge '$cwd/STL/slice$a.$i.stl';\nSurface Loop(1) = {1};\nVolume(1) = {1};" >$cwd/GEO/slice$a.$i.geo
        # convert .stl to .msh
        gmsh/gmsh $cwd/GEO/slice$a.$i.geo -3 -o $cwd/slice$a.$i.msh -format "msh22"