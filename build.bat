rmdir /s /q build_win
mkdir build_win
cd build_win
cmake -Wno-dev ..
cmake --build . --config Release
cd ..