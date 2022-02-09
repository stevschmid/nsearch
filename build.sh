rm -rf build
mkdir build
cd build
cmake -Wno-dev -D CMAKE_BUILD_TYPE=Release ..
make
# cmake -Wno-dev ..
# cmake --build . --config Release

# Additionally, either copy the nsearch file from build/nsearch to /usr/local/bin or add the directory to path
# sudo cp nsearch/nsearch /usr/local/bin/
# something like echo $PATH:$PWD/build/nsearch added to bashrc or something
