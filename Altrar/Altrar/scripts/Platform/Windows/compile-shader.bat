rmdir /s /q bin\shaders
mkdir bin
mkdir bin\shaders
..\vendor\bin\Windows\glslc.exe shaders/shader.vert -o bin/shaders/vert.spv
..\vendor\bin\Windows\glslc.exe shaders/shader.frag -o bin/shaders/frag.spv
