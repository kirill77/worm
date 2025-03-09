#include "simulation/World.h"
#include "simulation/Worm.h"

int main()
{
    std::shared_ptr<Worm> pWorm = std::make_shared<Worm>();
    World world(pWorm);
    for (uint32_t u = 0; u < 1000; ++u)
    {
        world.simulateStep(0.1);
    }
}
