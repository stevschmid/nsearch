rmdir /s /q build_win
mkdir build_win
cd build_win
cmake -Wno-dev ..
cmake --build . --config Release
cd ..

## Add the folder with nsearch.exe to system path
## Temporarily:
## SET PATH=%PATH%;%cd%\build_win\nsearch\Release
## Permanently:
## Either add the path above permanently using setx or the windows registry editor. Read about risks involved beforehand.
## Or copy the nsearch.exe executable to a folder already in PATH