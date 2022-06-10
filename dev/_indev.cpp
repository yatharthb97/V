



/** @brief Casts a 2D vector to a 3D vector. z-value is assigned zero.
 * @param x and y values. */
inline V cast_2D(double x, double y)
{
    return V(x, y, 0);
}

/** @brief Casts a 2D vector to a 3D vector. z-value is assigned zero.
 * @param Array of size 2 with x and y values. */
inline V cast_2D(double* array)
{
    return V(array[0], array[1], 0);
}







//Mutators

    //17
    //! Does uniform scaling of all the components.
    /*
    \param double - scaling constant
    */
    inline void scale(double scale) 
    {
        x = x * scale; y = y * scale, z = z * scale;
    } //End of scale()



    
/*    //36
    //! Overloaded Array Subindex Operator. Accepts three valid character: 'x', 'y', and 'z', and generates a compiler error if any other character is passed. Correct character param returns the corresponding component of the vector. x → 120, y → 121, z → 122; X → 88, Y → 89, Z → 90
    inline double operator[] (const char* index) const
    {
        if(*index == 'x')
            return this->x;
        else if(*index == 'y')
            return this->y;
        else if(*index == 'z')
            return this->z;
        else
            exit(1);
            //static_assert(false, "Vector V index out of bounds.");
    } //End of operator[]
   // error: invalid conversion from ‘const char*’ to ‘int’ [-fpermissive]
    
    //37
    //! Overloaded Array Subindex Operator for int parameter. Accepts three valid indices: 0, 1, and 2, and generates a compiler error if any other index is passed. Correct index param returns the corresponding component of the vector.
    inline double operator[] (int index) const
    {
        static_assert((index <= 2 && index >= 0), "Vector V index out of bounds.");

        if(index == 0)
            return this->x;
        else if(index == 1)
            return this->y;
        else if(index == 2)
            return this->z;
    } //End of operator[] overloaded*/



//////=========================== END OF OPERATORS ========================================



        //16
    //! Converts the given vector into its reflection : Negative. Alias "reflect" also defined for this function.
    inline void comp_neg()
    {
        this->x = -x;
        this->y = -y;
        this->z = -z;
    } //End of make_neg()

        //Alia "reflect"
    void __attribute__((always_inline)) comp_reflect()
    {
        comp_neg();
    }



    


    //22
    //! Multiply corresponding components of a vector.
    inline V comp_mul(V &other) const
    {
        return V(x*other.x, y*other.y, z*other.z);
    } //End of compmul()


    //23
    //! Convert the Vector into a null Vector.
    inline void comp_null()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    } //End of null()