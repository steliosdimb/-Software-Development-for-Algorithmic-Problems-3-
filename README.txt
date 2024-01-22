cmake -DCGAL_DIR=$CMAKE_INSTALLED_PREFIX/lib/CGAL-DCMAKE_BUILD_TYPE=Release .

make

./PROJECT_3 -i data/uniform/ -o PROJECT_output.txt
