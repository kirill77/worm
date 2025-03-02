#include "simulation/World.h"

int main()
{
    World world;
    for (uint32_t u = 0; u < 1000; ++u)
    {
        world.simulateStep(0.1);
    }
}
