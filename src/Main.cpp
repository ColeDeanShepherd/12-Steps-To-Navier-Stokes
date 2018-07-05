#include <math.h>
#include <sstream>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#include "App.h"

int main(int argc, char* argv[])
{
    App app;
    app.run();

    return 0;
}