#!/bin/bash
python -c "import gmsh; gmsh.initialize(); gmsh.open('fivecoils.geo'); gmsh.write('fivecoils.msh')"