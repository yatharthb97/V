# V for Vector

```V``` for Vector is a simple implementation of a 3-vector in C++. It is restricted to the type ```double``` as it is primarily intended for scientific computing use. The aim of this library is to provide a set of diverse operations that are quite trivial in its nature but their error free implementation is essential for many problems in scientific computing. 

The library is published as a stand-alone entity as it finds ubiquitous use in many projects (except for use of STL headers). At the same time, it needs to be reliably maintained and expanded for needs of newer projects.

The library is also a **personal project** and primarily caters to my own use. However, for the same reason — comprehensibility, efficiency, and accuracy are a primary concern. That said, the library is under continuous development and might be intrinsically broken. 

The`dev` directory contains the new developments.



### Development Stages

1. `vforvector` is compiling without errors.
2. Add "New Ideas" & compile.
3. Complete `Mesh` modue.
4. Complete `Lattice` module.



### Members

```c++
X 	//-> x-component
Y 	//-> y-component
Z 	//-> z-component
static V::tolerance //-> Tolerance used in floating point comparisions (allowed error)
static V::floatfield //-> Output formatting for string representation
static V::precision //-> Precision used for formatting string output
```

### Categories of functions

1. Constructors

2. Info functions

3. Component Accessors

4. Const Operators

5. Size functions

6. Norm functions

7. Negative functions

8. Check functions

9. Product functions

10. Scalar functions

11. Templates and Generators

12. Non-const operators (mutating)

13. Line (linear geometry) functions

19. Component-wise mutating operations (same operation performed for each component)

    ---

    **<u>New Ideas</u>** ↓

    ---

20. Random Orthogonal functions

16. Random Unit Generators

17. Rotation functions

18. Comparator functions

19. Projection functions

20. Linear algebra functions

21. V-Set Generators

22. Generic Summations and Multiplication functions

23. Inverse functions

23. From string construction

## Planned Features:

1. Single header file include library. [DONE]

2. **Explicit Vectorization** (this is a perfect use case for trivial vectorization, pun intended).

3. Has a consistent library grammar (clear syntactic separation between accessors and modifiers).

4. Has unit tests.

5. **Experiments and implements all accessor forms ↓**

   ```c++
   V vec(1.0, 2.0, 3.0);
   std::cout << vec.x();  //working
   std::cout << vec.X;    //working
   std::cout << vec[0];   //working
   std::cout << vec[xx];  //working
   std::cout << vec['x']; //working
   
   //All statements output ==> 1.0
   ```

   

## Library Grammar

1. The library has 3 members x, y, z for the three vector components. It also has some other static members like ```V::tolerance```.

2. All operations on the object usually returns a new ```V``` object (`const` operations).

3. Methods that start with "comp_xxx()," example : ```V::comp_square()``` — modify the current object. The suffix *"comp"* stands for *component* and performs the same operation on all components.

   ```c++
   V eg_square(0.5, 0, 4);
   eg_square.comp_square();
   std::cout << eg_square.info();
   ==> Output : (0.25, 0, 16)
   ```

4. Any Boolean operations that check a particular condition either start with the suffix `is_` or `has_`.  Example: ```V::is_null()``` checks if the given object is a null vector.

4. Any function with the suffix `f` as in `V::fdot()` , whose counterpart is `V::dot()`, represents **fast** and implies that it is the optimized version of its counterpart and lacks certain functionality (e.g. tolerance based decisions). 

### Todo

DONE. Change `info()` to `str()` .

DONE. Compile and test struct:

   ```c++
   typedef struct
   {
       double X, Y, Z = 0.0;
   }__attribute__((packed, aligned)) vec;
   
   class V : vec
   {
       V(double x, double y, double z)
       {
           X = x; Y = y; Z = z;
       }
       
       // Other definations
   }
   
   // In another function
   V eg(1.1, 2.2, 3.3);
   double x_ = eg[0];
   double y_ = eg[1];
   double z_ = eg[2];
   ```

   

DONE. Introduce Secondary classes:

1. `Q` for Quat (Quaternions)
2. `C` or `Ci` for Complex Numbers

3. Generic Summations and Multiplication functions

   ```c++
   V V::Sum(const Container<V> &c)
   {
       //Returns the component-wise sum of the container of vectors
   }
   
   V V::Multiply(const Container<V> &c)
   {
       //Returns the component-wide product of the container of vectors
   }
   ```
4. Float locale memeber that defines string representation → std::fixed, std::scientific, etc...


## Grammar

1. The templates and generator functions modify the current member function.
2. All minor operations return a modified vector.



### Printing Objects

1. `str()` : 1-space separated values.
2. `rep()` :  Representation - comma separated values in appropriate brackets.
3. `expand()` : Components added with literals. C(1.0, 2.0) → 1.0 + i2.0 .
4. `format()` : Returns a custom formatted output (only available for `V` class.)



## Mesh

Mesh module provides a set of objects that can be used to generate specific meshes. `Mesh` object defines a virtual interface for all mesh implementations. `BlankMesh` object can be used to import coordinates from a file (`BlankMesh::import()`). Mesh objects follow a **Generating Function** like interface. 

```c++
Spiral2DMesh<1000, 1500> mesh;
mesh.generate(); //Compute mesh values and store.
std::vector<V> coordinates = mesh.get_mesh(); //Read the coordiantes all at once -> same as "mesh.mesh".
V coord = mesh.next(); //Read coordinates successively.
```

### 1D Meshes

2D rectangular meshes can be used to generate 1D meshes by setting the `Y` component/count to zero.

### 2D Meshes

1. `Raster2DMesh`
2. `Spiral2DMesh`
3. `Circular2DMesh`

### 3D Meshes

1. `CuboidalMesh` (Can be generated using a `Raster2DMesh` or a `Spiral2DMesh` object)
2. `SphericalMesh`



## Lattice

`Lattice` module provides a set of objects to generate standard Bravais  lattice fields.

** Planned for stage 3**