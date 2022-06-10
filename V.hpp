#pragma once
//Header file class V for Vector → <V.hpp>

//><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
//><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

 /*$    /$$   ||  --- V for Vector ---
| $$   | $$   ||  (A scientific 3-vector library) 
| $$   | $$   ||  Developed by: Yatharth Bhasin       || Licence: *******************
|  $$ / $$/   ||  Computational Physicist             || Compile : NOk
 \  $$ $$/    ||  IIT Indore  | TIFR Hyderabd         || Tested : NOK
  \  $$$/     ||  (yatharth1997+git@gmail.com)        || Documentation : NOK
   \  $/      ||  (github: yatharthb97)               ||      
    \*/  
                      
//><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

//Preprocessor Commands


#include <cmath>    //! For basic math operations.
#include <sstream>  //! For info functions for printing.
#include <iomanip>  //! For 


#ifndef DISABLE_VFORVECTOR_MACROS
    //!!!!! --> No spaces  ○  Used with operator[]. 
    #define xx 0    // Used as - vect[xx] vs vect['x']
    #define yy 1    // Used as - vect[yy] vs vect['y']
    #define zz 2    // Used as - vect[zz] vs vect['z']
#endif

//**********************************************************************************

//! Represents 3D Vectors. A class of variables that defines a three component Vector.
class V
{
public:

//—————————————————————————————— MEMBERS ————————————————————————————————————————
	#pragma pack(push, 8) //Structure alignment
        double X; //! x-component of V.
    	double Y; //! y-component of V.
    	double Z; //! z-component of V.
    #pragma pack(pop)

    inline static double tolerance = 1e-10; //!< Tolerance limit for the class.
    inline const static unsigned int cardinality = 3; //!< Number of elements in each vector.
    inline static unsigned int precision = 5;//!<  Precision used for formatting string output.

    //Friend Declarations
    friend std::ostream& operator<< (std::ostream &stream, const V &vec);

//—————————————————————————————— CONSTRUCTORS ————————————————————————————————————————

    //0
    /**
     * @brief Default class constructor. Initialises each component to 0.0 .*/
    V():X(0.0), Y(0.0), Z(0.0)
    {}


    //1
    /*
    @brief class constructor. Copy a vector(or construt from a vector).
    @param V &other - const passed by reference that gets copied. */
    V(V const &other): X(other.X), Y(other.Y), Z(other.Z)
    {}


    //2
    /**
     * @brief Overloaded class constructor. Initialises each component to passed parameters.
     * @param x,y,z components. */
    V(double X, double Y, double Z): X(X), Y(Y), Z(Z)
    {} 

//—————————————————————————————— COMPONENT ACCESSOR  ———————————————————————————————————————
    //3
    //! Accessor for component -> x.
    double __attribute__((always_inline)) x() const
    {
        return X;
    }

    //4
    //! Accessor for component -> y.
    double __attribute__((always_inline)) y() const
    {
        return Y;
    }

    //5
    //! Accessor for component -> z.
    double __attribute__((always_inline)) z() const
    {
        return Z;
    }

    //6
    //! Accessor for component -> x squared.
    double __attribute__((always_inline)) x_sq() const
    {
        return X * X;
    }

    //7
    //! Accessor for component -> y squared.
    double __attribute__((always_inline)) y_sq() const
    {
        return Y * Y;
    }

    //8
    //! Accessor for component -> z squared.
    double __attribute__((always_inline)) z_sq() const
    {
        return Z * Z;
    }

    //9
    /**
     * \brief Operator overload that supports unsigned integer based array-like indexing.
     * \example V vec; x = vec[0].
     * \warning This definaion does not protect against memory leak.
     * */
    double operator[](int idx)
    {
        double* ptr = &(this->X);
        return ptr[idx];
    }

    //10
    /** @brief Operator overload that supports character based indices.
     * \warning This definaion does not protect against memory leak.
     * \example V vec; x = vec['x']. */
    double operator[](char idx)
    {
        // If safety is enabled - enforce lowercase.
        #if VFORVECTOR_SAFETY == 1
            idx = std::tolower(idx);
        #endif
        double* ptr = &(this->X);
        return ptr[static_cast<unsigned int>(idx) - static_cast<unsigned int>('x')];
    }


//————————————————————————————————— TEMPLATES & GENERATORS ———————————————————————————————————————————————

    //11
    /**
     * \brief Initalize x, y, and z with the given arguement.
     * \param double initalization value. */
    inline void xyz(double val)
    {
        X = val; Y = val; Z = val; 
    }

    //12
    /**
     * \brief Initalize x and y with the given arguement, and z is set to zero.
     * \param double initalization value. */
    inline void xy(double val)
    {
        X = val; Y = val; Z = 0.0;
    }

    //13
    /**
     * \brief Initalize y and z with the given arguement, and z is set to zero.
     * \param double initalization value. */
    inline void yz(double val)
    {
        X = 0.0; Y = val; Z = val;
    }

    //14
    /**
     * \brief  Initalize z and x with the given arguement, and z is set to zero.
     * \param double initalization value. */
    inline void zx(double val)
    {
        X = val; Y = 0.0; Z = val;   
    }

    //15
    /**
     * \brief Generates a vector by adding a given scalar and generating function.
     * \param Scalar common to all the three components. Default is 0.
     * \param A generating function pointer that returns a double.
     * \return A generated 3 vector V.*/
    V generate_add(double(*generating_fn)(), double scalar = 0)
    {
        //V tmp;
        this->X = scalar + generating_fn();
        this->Y = scalar + generating_fn();
        this->Z = scalar + generating_fn();
        //return tmp;
    }

    //15
    //! Alias of `generate()`.
    void generate(double(*generating_fn)(), double scalar = 0)
    {
        generate_add(generating_fn, scalar);
    }

    //16
    /**
     * /brief Generates a vector by multiplying a given scalar and generating function.
     * /param Scalar common to all the three components. Default is 1.
     * /param A generating function pointer that returns a double. */
    void generate_mul(double (*generating_fn)(), double scalar=1)
    {
        //V temp;
        this->X = scalar * generating_fn();
        this->Y = scalar * generating_fn();
        this->Z = scalar * generating_fn();
        //return tmp;
    }

    //17
    /**
     * /brief Generates a vector by operating the scalar and the generating
     *  function with the passed operation (3rd parameter).
     * /param Scalar common to all the three components.
     * /param A generating function pointer that returns a double.
     * /param The operation to perform between the scalar and the vector. */
    void generator(double (*generating_fn)(), double scalar, double (*op)(double, double))
    {
        this->X = op(scalar, generating_fn());
        this->Y = op(scalar, generating_fn());
        this->Z = op(scalar, generating_fn());
    }

    //18
    /**
     * /brief Generates a vector by operating the two generating
     *  functions with the passed operation (3rd parameter).
     * /param First generating function pointer that returns a double.
     * /param Second generating function pointer that returns a double.
     * /param The operation to perform between the two generating functions. */
    V generator(double(*generating_fn1)(), double (*generating_fn2)(), double (*op)(double, double))
    {
        this->X = op(generating_fn1(), generating_fn2());
        this->Y = op(generating_fn1(), generating_fn2());
        this->Z = op(generating_fn1(), generating_fn2());
    }



//————————————————————————————————— SCALAR FUNCTIONS ————————————————————————————————————————————————

    //19
    /**
     * \brief Converts the object to a scalar mask. 
     * The scalar is assigned to the X component.
     * The rest of the components are set to zero.
     * \param Scalar to set. */
    void inline make_scalar(double scalar)
    {
        this->X = scalar;
        this->Y = 0.0;
        this->Z = 0.0;
    }

    //20
    /**
     * @brief Retrives the scalar value. It is just a wrapper around V::x() accessor. */
    double inline get_scalar() const
    {
        return this->x();
    }

//—————————————————————————————— INFO FUNCTIONS ————————————————————————————————————————

    //21
    /**
    * @brief Print Vector Info function. Prints the vector with comma seperated components and enclosed in brackets. */
    std::string rep() const 
    {
        std::ostringstream stream;
        stream << "(" << X << ", " << Y << ", " << Z << ")";
        return stream.str();
    }

    //20
    /**
    * @brief Print Vector Info function seperated by 'sep'. Prints the vector with custom seperator seperated components and enclosed in brackets. */
    std::string str(const char* sep=" ") const 
    {
        std::ostringstream stream;
        stream << std::setprecision(V::precision);
        stream << X << sep << Y << sep << Z;
        return stream.str();
    }

     
    //21
    /**
    * \brief Format the vector with a given seperator and given bounds.
    * \param sep Seperator used between the components.
    * \param bounds The boundaries used for wrapping the vector string. 
    * The passed string is divided by two and appended on each side. */
    std::string format(std::string sep, std::string bounds) const 
    {
        std::ostringstream stream;
        stream << std::setprecision(V::precision);
        stream << bounds.substr(0, bounds.size()/2) 
               << X << sep << Y << sep << Z 
               << bounds.substr(bounds.size()/2, bounds.size()-1);
        return stream.str();
    }


    /** 
     * /todo i, j, k hat forms.*/
    std::string expand() const 
     {
        std::ostringstream stream;
        stream << std::setprecision(V::precision);
        stream << x() << "i + " << y() << "j + "<< z() << "k";
         return stream.str();
     }


//————————————————————————————————— COMMON UTILITY ————————————————————————————————————————————————
    
    //22
    /** \brief Function returns size of the vector.*/
    double inline size() const 
    {
        return std::sqrt(this->x_sq() + this->y_sq() + this->z_sq());
    }

    //23
    /** \brief Function returns the size squared of the vector. */
    double inline size_sq() const
    {
        return (this->x_sq() + this->y_sq() + this->z_sq());
    }

    //24
    /** \brief Returns a single sum of all components. */ 
    double inline reduce() const
    {
        return x() + y() + z();
    }

    //25
    /** \brief Returns a vector with each component squared. */
    V inline square() const
    {
        return V(x()*x(), y()*y(), z()*z());
    }

    //25
    /** \brief Returns a vector with each component squared. */
    double inline sq_reduce() const
    {
        V sq = this->square();
        return sq.reduce();
    }

    //27.1
    /**
     * @brief Return a normalised vector that has unit length.
     * Only normalises if the vector has non-zero length(greater than
     * tolerance limit), else all components are zero.
     * \warning If the size of the vctor is very small, then  */    
    V inline normalise() const
    {
        double tot = 1/this->size(); //Always positive
        //tot = (tot > 0) * (1.0 / tot) + 0; //Is this meaningless?
        return V(this->X*tot, this->Y*tot, this->Z*tot);

    } //End of normalise()

    //27.2
     /**
     * \brief Alias of V::normalise(). */ 
    V inline norm() const
    {
        return normalise();
    }

    //28
    //! Returns the negative of the given vector
    V inline neg() const 
    {
        return V(-X, -Y, -Z);
    }

    /** \brief Returns the absolute value of each */
    V inline abs() const
    {
        return V(std::fabs(x()), std::fabs(y()), std::fabs(z()));
    }


//————————————————————————————————— PRODUCT FUNCTIONS ————————————————————————————————————————————————

    //32
    /**
     * @brief Calculates the dot product between two vectors. Returns 0 if the answer is below tolerance.
     * @return double dot product. */
    double inline dot(const V& other) const
    {
        double temp = (x() * other.x()) + (y() * other.y()) + (z() * other.z());
        
        return (std::fabs(temp) > V::tolerance) * temp + 0;
    }

    
    //33
    /**
     * @brief Fast Dot - Calculates the dot product between two vectors. Does not check for tolerance.
     * @return return double dot product. */
    double inline fdot(const V& other) const
    {
        return (x() * other.x() + y() * other.y() + z() * other.z());
    }


    //34
    /**
     * @Cross product of two vectors.
     * @param B V other vector.
     * @return cross product - V. */
    V inline cross(const V& other) const 
    {
        return V(y() * other.z() - z() * other.y(),
                -x() * other.z() + z() * other.x(),
                 x() * other.y() - y() * other.x());
    }



    //35
    /**
     * @brief Returns the scalar triple product of three vectors.
     * @param Three vectors - A, B, C; such that - Scalar triple product: A • (B x C) . */
    double inline scalar_tri_prod(const V &A, const V &B, const V &C)
    {
        V tmp = B.cross(C);
        return A.dot(tmp);
    }

    //36
    /**
     * @brief Returns the vector triple product of three vectors.
     * @param Three vectors - A, B, C; such that - Vector triple product: A x (B x C) = B * A•C - C * A•B . */
    V inline vector_tri_prod(const V &A, const V &B, const V &C)
    {
        return (B * (A.dot(C))) - (C * A.dot(B));
    }


    //37.1 - Tolerance
    //! Calculates the angle between two vectors using arccos(cos inverse).
    double inline arccos(const V& other) const
    {
        double temp = this->dot();
        return std::acos(temp);
    } //End of arccos()


   //37.2
   //! Alias of V::arccos().  Calculates the angle between two vectors using arccos(cos inverse).
    double inline ang_in_rad(const V& other) const
   {
        return arccos(other);
   }


//————————————————————————————————— CONST OPERATORS  ————————————————————————————————————————————————

    
    //=, +=, -=, *=, /= operate on the existing Vectors.


    //13
    /**
     * @brief Returns a component wise subtraction result.
     * @param V other vector.
     * @return subracted vector. */
    inline V operator- (const V& other) const 
    {
        return V(x() - other.x(), y() - other.y(), z() - other.z());
    }


    inline V operator- (const double scalar) const 
    {
        return V(x() - scalar, y() - scalar, z() - scalar);
    }


    //14
    /**
     * @Returns the component wise addition result.
     * @param V other vector.
     * @return added vector. */
    inline V operator+ (const V &other) const
    {
         return V(x() + other.x(), y() + other.y(), z() + other.z());
    }

    inline V operator+ (const double scalar) const
    {
         return V(x() + scalar, y() + scalar, z() + scalar);
    }


    //15
    /**
     * \brief  Multiplication Operator. Scales the vector components with a scalar.
     * \param scale scalar -> template.
     * \return V scaled vector. */
    template <typename NumericType>
    V inline operator* (NumericType scalar) const
    {
      return V(X * scalar, Y * scalar, Z * scalar);     
    }


    //16
    /**
     * \brief  Division Operator. Scales the vector components with a scalar.
     * \param scale scalar -> template.
     * \return V scaled vector. */
    template <typename NumericType>
    inline V operator/ (NumericType scalar) const
    {
      return V(x()/scalar, y()/scalar, z()/scalar);     
    } //End of Operator/



    //17
    /**
     * @brief  Equality Operator with default tolerance.
     * @param other V other vector.
     * @return bool comparision result. */
    inline bool operator== (const V &other)  const
    {
        return ((x() - other.x()) < V::tolerance && 
                (y() - other.y()) < V::tolerance && 
                (z() - other.z()) < V::tolerance);
    } //End of operator==



    //18
    /**
     * @brief Unequality Operator with set tolerance.
     * @param o V other vector.
     * @return bool comparision result. */
    inline bool operator!= (V &other) const
    {
        return !(operator==(other));
    } //End of operator!=


//—————————————————————————————— Non-const Operators (mutators) ————————————————————————————————————————

    //49
    /**
     * @brief Assignment operator
     * @param V other - assigner. */
    void operator= (const V& other) inline
    {
        X = other.X;
        Y = other.Y;
        Z = other.Z;
    }

    //50
    /**
     * @brief  Iterative Subtraction Operator.
     * @param V other vector. */
    void operator-= (const V& other) inline 
    {
        X -= other.X;
        Y -= other.Y;
        Z -= other.Z;
    }


    //51
    /**
     * @brief  Iterative Addition Operator.
     * @param V other vector. */
    inline void operator+= (const V& other) 
    {
            X += other.X;
            Y += other.Y;
            Z += other.Z;
    }


    //52
    /**
     * @param Iterative Scaling Operator - double.
     * @param scale double -> template. */
    template <typename T>
    inline void operator*= (T scale) //MUTATOR
    {
        X *= scale; Y *= scale; Z *= scale;
    }


    //53
    /**
     * @brief Iterative Scaling Division Operator - double.
     * @param scale double -> template. */
    template <typename T>
    inline void operator/= (T scale)
    {
        X /= scale; Y /= scale; Z /= scale;
    }


















//————————————————————————————————— CHECK FUNCTIONS ————————————————————————————————————————————————
    
    //29
    /**
     * \brief Checks if the vector is a unit vector.
     * \return True if the vector is a unit vector. */
    bool inline is_unit() const
     {
        return (std::fabs(this->size()) - 1.0) <= V::tolerance;
    }

    //30
    /**
     * \brief Fast check if the vector is a unit vector. Needs to be used
     *  with causion as the operation is ususlly meaningless in most cases
     *  due to floating point rounding error.
     * \return bool answer.*/
    bool inline is_unit_f() const
     {
        return (this->size() == 1.0);
    }

    //31
    /**
     * \brief Checks if one or more components is NaN. 
     * Returns true if any of the component is NaN.
     * \return True if any component is NaN. */    
    bool inline is_nan() const
    {
        return (std::isnan(x()) ||
                std::isnan(y()) || 
                std::isnan(z()));
    }

    //32 - Tolerance
    /**
     * \brief Checks if the vector is a null vector with set tolerance
     * (All component must be below the tolerance).
     * \param (optional) tolerance. 
     * \return bool -> true is null. */
    inline bool is_null() const
    {
        V temp = this->abs();
        return (temp.x() <= V::tolerance && 
                temp.y() <= V::tolerance &&
                temp.z() <= V::tolerance);
    }


    //33
    /**
     * \brief Returns true if the vectors are orthogonal to each other.
     * The dot product must be zero.
     * */
    bool inline is_ortogonal(V &other) const
    {
        double dot = this->fdot(other);
        return dot <=V::tolerance;
    }


    //34
    /**
     * \brief Returns if the vector is a scalar. Scalar mask can be used to treat
     * a V object as a scalar. In that case, the X component assumes the scalar 
     * and the rest two components must be below the tolerance. */
    bool inline is_scalar() const
    {
        return std::fabs(this->X) >= V::tolerance && 
               std::fabs(this->Y) <  V::tolerance && 
               std::fabs(this->Z) <  V::tolerance;
    }


    //35
    /**
     * \brief Compares the two vectors without any error tolerance.
     * \param Other vector for comparision. */
    inline bool is_exactly_eq(const V &other) const
    {
        return (x() == other.x()) && (y() == other.y()) && (z() == other.z());
    }

    //36
    /**
     * \brief Returns whether the given vector is parallel to the other
     * vector with some tolerance. Equivalent to cross product is zero. */
    bool is_parallel(const V &other) const
    {
        return (std::fabs(this->cross(other)) - 1.0) <= V::tolerance;
    }


//—————————————————————————————— STANDARD REFS ————————————————————————————————————————

    void inline make_null()
    {
        this->X = 0; this->Y = 0; this->Z = 0;
    }

    //19
    //! Returns a unit vector along x-axis.
    static inline V unit_x() { return V(1, 0, 0); }

    //20
    //! Returns a unit vector along y-axis.
    static inline V unit_y() { return V(0, 1, 0); }

    //21
    //! Returns a unit vector along z-axis.
    static inline V unit_z() { return V(0, 0, 1); }

    //
    //! Returns a null vector.
    static inline V null() { return V(0, 0, 0); }




    //! Maximum precision for output is set to *15*.
    constexpr static unsigned int MAX_PRECISION = 15;

    /** \brief Sets the V::precision value.
     *  \param precision number of significant digits which is printed. 
     *  Default is `V::MAX_PRECISION`. Use the function without any param
     *  to use the maximum allowed precision. */
    static void set_precision(int precision = V::MAX_PRECISION)
    {
        V::precision = precision;
    }

private:

    bool within_tol(double comp) const
    {
        return std::fabs(comp) <= V::tolerance;
    }

};

//61 
//! Stream operator overload that puts V::info() to the ostream.
std::ostream& operator<< (std::ostream &stream, const V &vect)
{
    stream << vect.rep();
    return stream;
} //End of friend overload operator<<