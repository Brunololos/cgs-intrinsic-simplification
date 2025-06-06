Running on Windows:
1. navigate to this directory & mkdir build
2. cd build
3. cmake ..
4. open .sln in MS Visual Studio & compile as Release
5. run via:
    ./bin/Release/iSIMP_gui.exe (and pass args like ../meshes/spot.obj)
    ./bin/Release/iSIMP_gui.exe ../meshes/spot.obj tex=2000 snap=200
    ./bin/Release/iSIMP_gui.exe

   args: (if no prefixes are provided, this same argument order is expected.)
    mesh=[.obj mesh file path]
    tex=[integer texture resolution]
    snap=[integer number of vertex removals]
    verbose=[true/false]
