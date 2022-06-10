

#include "mesh.hpp"

int main()
{
	Raster2DMesh<10, 5> mesh1;
	Spiral2DMesh<10, 5> mesh2;

	mesh1.generate();
	mesh2.generate();

	mesh1.export("mesh1.dat");
	mesh2.export("mesh2.dat");
}

class 