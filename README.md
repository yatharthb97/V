# V for Vector

```V``` for Vector is a simple implementation of a 3-vector in C++. It is restricted to the type ```double``` as it is primarily intended for scientific computing use. The aim of this library is to provide a set of diverse operations that are quite trivial in its nature but are essential for many problems in scientific computing. 

The library is published as a stand-alone entity as it finds ubiquitous use in many projects. At the same time, it needs to be reliably maintained and expanded for needs of newer projects.

The library is also a **personal project** and primarily caters to my own use. However, for the same reason comprehensibility, efficiency, and accuracy are a primary concern. That said, the library is under continuous development and will not make a version release anytime soon. 

## Planned Features:

1. Single header file include library.

2. Explicit Vectorization.

3. Has a consistent library grammar (clear syntactic separation between accessors and modifiers).

4. Has unit tests.

5. **Experiment with all accessor forms ↓**

   ```c++
   V vec(1.0, 2.0, 3.0);
   std::cout << vec.x();
   std::cout << vec.x;
   std::cout << vec[0];
   std::cout << vec[xx];
   std::cout << vec['x'];
   
   //All statements output ==> 1.0
   ```

   

## Library Grammar

1. The library has 3 members x, y, z for the three vector components. It also has some other static members like ```V::tolerance```.

2. All operations on the object usually returns a new ```V``` object.

3. Methods that start with "comp_xxx()," example : ```V::comp_square()``` — modify the current object. The suffix *"comp"* stands for *component* and hence modifies the components of the current vector.

   ```c++
   V eg_square(0.5, 0, 0);
   eg_square.comp_square();
   std::cout << eg_square.info();
   ==> Output : (0.25, 0, 0)
   ```

4. Any Boolean operations that check a particular condition either start with the suffix *"is" or "has_"*.  Example: ```V::is_null()``` checks if the given object is a null vector.

