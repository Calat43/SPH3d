#include "utils.h"

int max(int a, int b) {
    return a > b ? a : b;
}

int min(int a, int b) {
    return a < b ? a : b;
}

/*
//TODO IT MIGHT BE INEFFICIENT!!!
for(int i = 0; i < next_grid.x_size; ++i)
{
for(int j = 0; j < next_grid.y_size; ++j)
{
for(int k = 0; k < next_grid.z_size; ++k)
{
Cell & cell = next_grid.cells[i][j][k];
std::vector<Particle> remaining_particles;
for(auto particle: cell.gas_particles)
{
if(particle.x >= params.left && particle.x <= params.right &&
        particle.y >= params.right && particle.z <= params.left &&
        particle.z >= params.right && particle.z <= params.left)
{
continue;
}
else
{
remaining_particles.push_back(particle);
}
}
cell.gas_particles.erase(cell.gas_particles.begin(), cell.gas_particles.end());
cell.gas_particles.insert(cell.gas_particles.begin(), remaining_particles.begin(), remaining_particles.end());

}
}
}
*/
