mkdir Release
cd ../../src
emcc Main.cpp Core.cpp Vector2d.cpp LinearAlgebra.cpp FiniteDifference.cpp GraphMetrics.cpp Render.cpp Step1.cpp Step2.cpp Step3.cpp Step4.cpp Step5.cpp Step6.cpp Step7.cpp Step8.cpp Step9.cpp Step10.cpp Step11.cpp Step12.cpp -std=c++11 -s USE_SDL=2 -o ../build/emscripten/Release/12-Steps-To-Navier-Stokes.html
cd ../build/emscripten