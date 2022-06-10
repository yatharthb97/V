#include <iostream>
#include <iomanip>
#include "V.hpp"


double gen1() { return 3.14; }
double gen2() { static double i = 0; i=i+1; return i; }

double op(double a, double b) { return a+b; }


V vec1;
V vec2(1.23456,2.34567,3.45678);
V vec3(vec2);
V vec4(2, 2, 2);




void test2()
{
	std::cout << "vec4 | size: " << vec4.size() << std::endl;
	std::cout << "vec4 | size^2: " << vec4.size_sq() << std::endl;
	std::cout << "vec4 | reduce: " << vec4.reduce() << std::endl;

	std::cout << "vec4 | sq: " << vec4.square() << std::endl;
	std::cout << "vec4 | sq_reduce: " << vec4.sq_reduce() << std::endl;

	std::cout << "vec4 | norm: " << vec4.norm() << "\t size: " << vec4.norm().size() <<  std::endl;
	std::cout << "vec4 | neg: " << vec4.neg() << std::endl;
	std::cout << "vec4 | abs: " << vec4.neg().abs() << std::endl;



}



void test1()
{
	std::cout << std::boolalpha;
	std::cout << "vec1 | " << vec1.rep() << std::endl;
	std::cout << "vec2 | " << vec2.rep() << std::endl;
	std::cout << "vec2 | " << vec3.rep() << std::endl;
	std::cout << "vec2 | " << vec2.format(" ::: ", "[{}]") << std::endl;
	std::cout << "vec3 | " << vec2.str() << std::endl;


	std::cout << "vec2 | x: " << (vec2[0] == vec2.X) << std::endl;
	std::cout << "vec2 | y: " << (vec2['y'] == vec2.Y) << std::endl;
	std::cout << "vec2 | z: " << (vec2[zz] == vec2.Z)<< std::endl;

	std::cout << "vec2 | x: " << (vec2[0] == vec2.X) << std::endl;
	std::cout << "vec2 | y: " << (vec2[1] == vec2.Y) << std::endl;
	std::cout << "vec2 | z: " << (vec2[2] == vec2.Z) << std::endl;


	std::cout << "unit-x | " << V::unit_x().rep() << std::endl;
	std::cout << "unit-y | " << V::unit_y().rep() << std::endl;
	std::cout << "unit-z | " << V::unit_z().rep() << std::endl;

	//Stream tests
	std::cout << "stream-operator | " << vec2 << std::endl;

	//Template tests

	V vec; vec.xy(1.23456789);
	std::cout << "xy  | "<< vec << std::endl;

	vec.make_null(); vec.yz(1.23456789);
	std::cout << "yz  | "<< vec << std::endl;


	vec.make_null(); vec.zx(1.23456789);
	std::cout << "zx  | "<< vec << std::endl;


	vec.make_null(); vec.xyz(1.23456789);
	std::cout << "xyz | "<< vec << std::endl;


	// Generator tests

	vec.make_null(); vec.generate_add(gen1, 2);
	std::cout << "gen-add | "<< vec << std::endl;

	vec.make_null(); vec.generate(gen1, 2);
	std::cout << "gen | "<< vec << std::endl;

	vec.make_null(); vec.generator(gen1, 2, op);
	std::cout << "gener | "<< vec << std::endl;

	vec.make_null(); vec.generator(gen1, gen1, op);
	std::cout << "gener2 | "<< vec << std::endl;

}

int main()
{
	test2();
}