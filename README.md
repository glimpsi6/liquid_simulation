for Release:
```bash
mkdir build && cd build
cmake ..
make -j
```

for ASan:
```bash
mkdir build-asan && cd build-asan
cmake .. -DENABLE_ASAN=ON -DCMAKE_BUILD_TYPE=Debug
make -j
```

for Native:
```bash
mkdir build && cd build
cmake .. -DENABLE_NATIVE=ON -DCMAKE_BUILD_TYPE=Release
make -j
```