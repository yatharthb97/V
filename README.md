# V for Vector

```V``` for Vector is a simple implementation of a 3-vector in C++. It is restricted to the type ```double``` as it is primarily intended for scientific computing use. The aim of this library is to provide a set of diverse operations that are quite trivial in its nature but their error free implementation is essential for many problems in scientific computing. 

The library is published as a stand-alone entity as it finds ubiquitous use in many projects (except for use of STL headers). At the same time, it needs to be reliably maintained and expanded for needs of newer projects.

The library is also a **personal project** and primarily caters to my own use. However, for the same reason — comprehensibility, efficiency, and accuracy are a primary concern. That said, the library is under continuous development and might be intrinsically broken. 

The`dev` directory contains the new developments.

### Members

```c++
X 	//-> x-component
Y 	//-> y-component
Z 	//-> z-component
static V::tolerance //-> Tolerance used in floating point comparisions (allowed error)
static V::floatfield //-> Output formatting for string representation
static V::precision //-> Precision used for formatting string output
```

### Categories

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

14. Random Orthogonal functions

15. Random Unit Generators

16. Rotation functions

17. Comparator functions

18. Projection functions

19. Component-wise mutating operations (same operation performed for each component)

    ---

    **<u>New Ideas</u>** ↓

    ---

20. Linear algebra functions

21. V-Set Generators

22. Generic Summations and Multiplication functions

23. Inverse functions

## Planned Features:

1. Single header file include library.

2. **Explicit Vectorization** (this is a perfect use case for trivial vectorization, pun intended).

3. Has a consistent library grammar (clear syntactic separation between accessors and modifiers).

4. Has unit tests.

5. **Experiments and implements all accessor forms ↓**

   ```c++
   V vec(1.0, 2.0, 3.0);
   std::cout << vec.x();
   std::cout << vec.X;
   std::cout << vec[0];
   std::cout << vec[xx];
   std::cout << vec['x'];
   
   //All statements output ==> 1.0
   ```

   

## Library Grammar

1. The library has 3 members x, y, z for the three vector components. It also has some other static members like ```V::tolerance```.

2. All operations on the object usually returns a new ```V``` object.

3. Methods that start with "comp_xxx()," example : ```V::comp_square()``` — modify the current object. The suffix *"comp"* stands for *component* and performs the same operation on all components.

   ```c++
   V eg_square(0.5, 0, 4);
   eg_square.comp_square();
   std::cout << eg_square.info();
   ==> Output : (0.25, 0, 16)
   ```

4. Any Boolean operations that check a particular condition either start with the suffix `*"is_" or "has_"*`.  Example: ```V::is_null()``` checks if the given object is a null vector.

4. Any function with the suffix `f` as in `V::fdot()` , whose counterpart is `V::dot()`, represents **fast** and implies that it is the optimized version of its counterpart and lacks certain functionality. 

### Todo

1. Change `info()` to `str()` .

2. Compile and test struct:

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
   x_ = eg[0];
   y_ = eg[1];
   z_ = eg[2];
   ```
   
   

3. Introduce Secondary classes:

   1. `Q` for Quat (Quaternions)
   2. `C` or `Ci` for Complex Numbers

4. Generic Summations and Multiplication functions

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

   

