//  __  __           _
// |  \/  | ___  ___| |__
// | |\/| |/ _ \/ __| '_ \
// | |  | |  __/\__ \ | | |
// |_|  |_|\___||___/_| |_|

/** @brief Describes objects for generation of meshes. */


#pragma once

#include <random>
#include <vector>
#include <cmath>
#include <cstdint>

#include "vforvector.hpp"



//Requirements from V
/*
	1. V::from_str()
	2. V.Z
	3. V.str() [Add seperator functionalities. ]
*/


/** @brief Generic `Mesh` object that contains common functions. This class contains pure 
 * virtual functions and it cannot be instantiated without being inherited and explicit
 * implementation of the present virtual functions. */
class Mesh
{
public:

	
	std::vector<V> mesh;	// Mesh container.
	std::uint64_t ptr = 0;	// Pointer for next() function.

	/** @brief Trivial and blank constructor. */
	Mesh()
	{}


	/** @brief Export the mesh to the given file in simple ASCII format.
	 * @param filename Name of the file to export. 
	 * @param (optional) precision Precison of floating point numbers to output to file. */
	V export(const std::string filename, unsigned int precision = 3) const
	{
		std::ostringstream buffer;
		buffer << std::setprecision(precision);

		unsigned int sizex = size();
		for(unsigned int i = 0; i < sizex; i++)
		{
			buffer << mesh[i].str() << '\n';
		}

		std::ofstream file(filename, std::ios::out);
		file << buffer.str();
		file.close();
	}

	/** @brief Generating function wrapper for mesh object. Calling of this function
	 *  sequentially returns different/ sequential mesh positions. */
	V next()
	{
		ptr_t ptr_now = ptr;
		ptr = (ptr+1 >= size()) * 0 + (ptr+1 < size()) * (ptr+1); // Check if the logic is correct.
		return mesh[ptr_now];
	}

	/** @brief Alias of `next()` member function. */
	V inline operator()()
	{
		return next();
	}

	/** @brief Returns the **current** size of the mesh. */
	unsigned int inline size() const
	{
		return mesh.size();
	}

	/** @brief Returns `true` if the mesh is empty. `false` otherwis. 
	 * \Note Use `generate()` to fill the mesh. */
	bool is_empty() const
	{
		return mesh.size() == 0;
	}

	/** @brief Randomise the mesh.
	 * @param PRNG prng A pseudo random number genertor object that is STL compatible. */
	template <class PRNG>
	void randomise(PRNG &prng)
	{
		std::shuffle(mesh.begin(), mesh.end(), prng);
	}


	/** @brief Returns a vector of the generated mesh. */
	const std::vector<V> get_mesh()
	{
		return mesh;
	}

// Virtual functions ↓

	virtual void generate() = 0;

};



class BlankMesh : public Mesh
{
public:

	BlankMesh(std::string load_file = ""): Mesh()
	{
		if(load_file != "")
		{
			this->import(load_file);
		}
	}


	/** @brief Import a file to mesh. */
	void import(const std::string load_file)
	{
		std::ifstream file(load_file, std::ios::in);
		
		if (file.is_open()) 
		{
		    std::string line;
		    while (std::getline(file, line)) 
		    {
		        this->mesh.emplace_back(V::from_str(line));
		    }
		}

		file.close();
	}

	/** @brief This function does nothing. Use `import()` to correctly use this class
	 *  implementation. */
	void generate() override
	{
		return;
	}

};


/** @brief Describes a rectangular 2D Mesh object. Sequential generaion of positions will
 *  reset x-axis coordinate after each scan while increasing y.
 * @param [Template] XLEN Mesh length along x-axis.
 * @param [Template] YLEN Mesh length along y-axis.
 * */
template <std::uint64_t XLEN, std::uint64_t YLEN>
class Raster2DMesh
{
	
	const std::uint64_t xlen = XLEN;
	const std::uint64_t ylen = YLEN;
public:

	/** @brief Trivial / Blank consructor.*/
	Raster2DMesh() : Mesh()
	{}

	/** @brief `generate()` implementation. */
	void generate() override
	{
		std::uint64_t x = 0;
		std::uint64_t y = 0;

		const std::uint64_t capacity = XLEN * YLEN;
		mesh.reserve(capacity);

		for(std::uint64_t i = 0; i < capacity; i++)
		{
			if(x >= xlen-1)
			{
				x = 0;
				y = (y+1) * (y+1 < ylen) + 0 * !(y+1 < ylen);
			}

			else
			{
				x++;
			}
					
			mesh.emplace_back(V(x, y, 0));
		}
	}

};


/** @brief Describes a rectangular 2D Mesh object. Sequential generaion of positions
 *  happens in a spiral fashion.
 * @param [Template] XLEN Mesh length along x-axis.
 * @param [Template] YLEN Mesh length along y-axis.
 * */
template <std::uint64_t XLEN, std::uint64_t YLEN>
class Spiral2DMesh
{
	
	const std::uint64_t xlen = XLEN;
	const std::uint64_t ylen = YLEN;

public:

	/** @brief Trivial / Blank consructor.*/
	Spiral2DMesh() : Mesh()
	{}


	/** @brief `generate()` implementation. */
	void generate() override
	{
		static int delx = 1;
		std::uint64_t x = 0;
		std::uint64_t y = 0;

		const std::uint64_t capacity = XLEN * YLEN;
		mesh.reserve(capacity);

		for(std::uint64_t i = 0; i < capacity; i++)
		{
			if(x >= xlen || x < 0)
			{
				delx = delx * -1;
				y = (y+1) * (y+1 < ylen) + 0 * !(y+1 < ylen);
				x+=delx;

				// Grid scanning has finished - repeat
				if(y == 0 && ylen>0)
				{
					x = 0;
					delx = +1;
				}

				// For 1D case, avoid taking two consequitive images.
				if(y==0 && ylen==0)
				{
					x += delx;
				}
			}

			mesh.emplace_back(V(x, y, 0));
			x+=delx;
		}
	}

};


/** @brief Cuboidal mesh object that generates 3D meshes.
 * @param [Template] ICNT Mesh length along x-axis.
 * @param [Template] JCNT Mesh length along y-axis.
 * @param [Template] KCNT Mesh length along z-axis.
 * @param [Optional & Template] Mesh2DType Type of 2D mesh used to generate the planes
 *  of the 3D mesh. Default typ is `Raster2DMesh`. */
template <std::uint64_t ICNT, std::uint64_t JCNT, std::uint64_t KCNT, class Mesh2DType=Raster2DMesh> 
class CuboidalMesh : public Mesh
{
	Mesh2DType<ICNT, JCNT> mesh2D;	//2D Mesh object.


public:

	/** @brief Trivial / Blank consructor.*/
	CuboidalMesh(): Mesh()
	{}

	/** @brief `generate()` implementation. */
	void generate() override
	{
		const std::uint64_t capacity = ICNT * JCNT * KCNT;
		mesh.reserve(capacity);

		mesh2D.generate();
		std::vector<V> meshII = mesh2D.get_mesh(); 
		V turtle;

		for(unsigned int z = 0; z < KCNT; z++)
		{
			for(unsigned int xy = 0; xy < mesh2D.size(); xy++)
			{
				turtle = meshII[xy];
				turtle.Z = z;
				mesh.emplace_back(turtle);
			}
			
		}


	}
};


template <double Radius, unsigned int RadiusCnt, unsigned int ThetaCnt> 
class SphericalMesh
{
	SphericalMesh(): Mesh()
	{}

	void generate() override
	{
		//Calculate parameters

		std::uint64_t d_radius = Radius / RadiusCnt;
		std::uint64_t d_angle = 2 * Constants::PI / ThetaCnt;
		double radius = d_radius;


		mesh.emplace_back(V()); //Add (0,0,0)
		for(std::uint64_t r = 0; r < RadiusCnt; r++)
		{
			radius += d_radius;
			V axis = V::unit_x();
			
			for(std::uint64_t th = 0; th < ThetaCnt; th++) // Iterate from 0 to 2PI
			{
				V local_it = V::unit_z(); //Make a copy of the z-iterator.
				axis.rotate(V::unit_z(), d_angle); //Increment Theta.

				for(std::uint64_t phi = 0; phi < ThetaCnt/2; phi++) //Iterate from +z to -z
				{
					mesh.emplace_back(local_it * radius); //Scaled by radius
					local_it.rotate(axis, d_angle); //Angle is rotated
				}
			}
		} 
	}
};





