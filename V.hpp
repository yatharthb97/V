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

//Macro Definaions
#define __VFORVECTOR_DEFAULT_TOLERANCE__ 1E-10 //! Defines the default tolerance value for the class.

#ifndef __DISABLE_VFORVECTOR_MACROS__
    //!!!!! --> No spaces   Useful with operator[]. 
    #define xx 0    // Used as - vect[xx] vs vect['x']
    #define yy 1    // Used as - vect[yy] vs vect['y']
    #define zz 2    // Used as - vect[zz] vs vect['z']
#endif

//**********************************************************************************

/**
 *  @brief - Represents 3D Vectors. A class of variables that defines a three component Vector.*/
class V
{
public:

//—————————————————————————————— MEMBERS ————————————————————————————————————————
	double X; //! x-component of V.
	double Y; //! y-component of V.
	double Z; //! z-component of V.
    
    constexpr static double tolerance = __VFORVECTOR_DEFAULT_TOLERANCE__; //! Tolerance limit for the class.


//—————————————————————————————— CONSTRUCTORS ————————————————————————————————————————

    //1
    /**
     * @brief Default class constructor. Initialises each component to 0.0 .*/
    V():X(0.0), Y(0.0), Z(0.0)
    {}


    //2
    /*
    @brief class constructor. Copy a vector(or construt from a vector).
    @param V &other - const passed by reference that gets copied. */
    V(V const &other): X(other.X), Y(other.Y), Z(other.Z)
    {}


    //3
    /**
     * @brief Overloaded class constructor. Initialises each component to passed parameters.
     * @param x,y,z components. */
    V(double X, double Y, double Z): X(X), Y(Y), Z(Z)
    {} 


//—————————————————————————————— INFO FUNCTIONS ————————————————————————————————————————

    //4
    /**
     * @brief Print Vector Info function. Prints the vector with comma seperated components and enclosed in brackets. */
   std::string info() inline const 
    {
        std::ostringstream stream;
        stream << "(" << X << ", " << Y << ", " << Z << ")";
        return stream.str();
    } //End of info()

    //5
    /**
     * @brief Print Vector Info function seperated by 'sep'. Prints the vector with custom seperator seperated components and enclosed in brackets. */
    std::string info(char sep) inline const 
    {
        std::ostringstream stream;
        stream << X << sep << Y << sep << Z;
        return stream.str();
    } //End of info()

    
    //6
    /**
     * @brief Print vector raw info function. Prints the vector with space seperated components and without brackets. */
    std::string info_raw() inline const 
    {
         std::ostringstream stream;
         stream << X << " " << Y << " " << Z;
         return stream.str();
    } //End of info_raw()



//—————————————————————————————— COMPONENT ACCESSOR  ————————————————————————————————————————
    //7
    //! Accessor for component -> x.
    double __attribute__((always_inline)) x() const
    {
        return X;
    }

    //8
    //! Accessor for component -> y.
    double __attribute__((always_inline)) y() const
    {
        return Y;
    }

    //9
    //! Accessor for component -> z.
    double __attribute__((always_inline)) z() const
    {
        return Z;
    }

    //10
    //! Accessor for component -> x squared.
    double __attribute__((always_inline)) x_sq() const
    {
        return X * X;
    }

    //11
    //! Accessor for component -> y squared.
    double __attribute__((always_inline)) y_sq() const
    {
        return Y * Y;
    }

    //12
    //! Accessor for component -> z squared.
    double __attribute__((always_inline)) z_sq() const
    {
        return Z * Z;
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
        return V(X - other.X, Y - other.Y, Z - other.Z);
    }



    //14
    /**
     * @Returns the component wise addition result.
     * @param V other vector.
     * @return added vector. */
    inline V operator+ (const V &other) const
    {
         return V( X + other.X, Y + other.Y, Z + other.Z );
    } //End of Operator+


    //15
    /**
     * @brief  Multiplication Operator. Scales the vector components with a scalar.
     * @param scale scalar -> template.
     * @return V scaled vector. */
    template <typename T>
    inline V operator* (T scale) const
    {
      return V(X * scale, Y * scale, Z * scale);     
    } //End of Operator*


    //16
    /**
     * @brief  Division Operator. Scales the vector components with a scalar.
     * \param scale scalar -> template.
     * @return V scaled vector. */
    template <typename U>
    inline V operator/ (U divscale) const
    {
      return V(X/divscale, Y/divscale, Z/divscale);     
    } //End of Operator/



    //17
    /**
     * @brief  Equality Operator with default tolerance.
     * @param other V other vector.
     * @return bool comparision result. */
    inline bool operator== (const V &other)  const
    {
        return ((X - other.X) < V::tolerance && 
                (Y - other.Y) < V::tolerance && 
                (Z - other.Z) < V::tolerance);
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


//————————————————————————————————— SIZE FUNCTIONS ————————————————————————————————————————————————
    //19
    /**
     * @brief Function returns size of the vector.*/
    double size() inline const 
    {
        return std::sqrt(this->x_sq() + this->y_sq() + this->z_sq());
    }

    //20
    /**
     * @brief Function returns the size squared of the vector. */
    double size_sq() inline const
    {
        return (this->x_sq() + this->y_sq() + this->z_sq());
    }


//————————————————————————————————— NORM FUNCTIONS ————————————————————————————————————————————————
    //21
    /**
     * @brief Return a normalised vector that has unit length. Only normalises if the vector has non-zero length(greater than tolerance limit), else all components are zero. */    
    V normalise() inline const
    {
        double tot = this->size(); //Always positive
        //tot = (tot > V::tolerance) * (1.0 / tot) + 0; //Is this meaningless?
        return V(this->X*tot, this->Y*tot, this->Z*tot);

    } //End of normalise()

    //22
     /**
     * @brief Alias of V::normalise(). Return a normalised vector that has unit length. Only normalises if the vector has non-zero length(greater than tolerance limit), else all components are zero. */ 
    V __attribute__((always_inline)) norm() const
    {
        return normalise();
    }



//————————————————————————————————— NEGATIVE FUNCTIONS ————————————————————————————————————————————————

    //23
    //! Returns the negative of the given vector
    inline V neg() inline const 
    {
        return V(-X, -Y, -Z);
    }

//————————————————————————————————— CHECK FUNCTIONS ————————————————————————————————————————————————
    
    //24
    /**
     * @brief Checks if the vector is a unit vector.
     * \return True if the vector is a unit vector. */
    bool is_unit(const double tolerance = V::tolerance) inline const
     {
        return (std::fabs(this->size() - 1.0) <= tolerance);
    }

    //25
    /**
     * @brief Fast check if the vector is a unit vector. Needs to be used with causion as the operation is ususlly meaningless in most cases due to floating point rounding error.
     * \return bool answer.*/
    bool fis_unit() inline const
     {
        return (this->size() == 1.0);
    }

    //26
    /**
     * @brief Checks if one or more components is NaN. Returns true if any of the component is NaN.
     * @return True if any component is NaN. */    
    inline bool is_nan() inline const
    {
        return (std::isnan(X) ||
                std::isnan(Y) || 
                std::isnan(Z));
    }

    //27
    /**
     * @brief Checks if the vector is a null vector with set tolerance (All component must be below the tolerance).
     * @param (optional) tolerance. 
     * @return bool -> true is null. */
    inline bool is_null(double tolerance = V::tolerance) inline const
    {
        V temp = this->abs();
        return (temp.X <= tolerance && 
                temp.Y <= tolerance &&
                temp.Z <=tolerance);
    }


    //28
    /**
     * @brief Returns true if the vectors are orthogonal to each other. The dot product must be zero.
     * */
    bool inline is_ortogonal(V &other, double tolerance = V::tolerance) inline const
    {
        double dot = this->fdot(other);
        return dot <=tolerance;
    }


    //29
    /**
     * @brief Returns if the vector is a scalar. Scalar mask can be used to treat a V object as a scalar. In that case, the X component assumes the scalar and the rest two components must be below the tolerance. */
    inline bool is_scalar(double tolerance = V::tolerance) inline const
    {
        return std::fabs(this->X) >= V::tolerance && 
               std::fabs(this->Y) < V::tolerance && 
               std::fabs(this->Z) < V::tolerance;
    }


    //30
    /**
     * @brief Compares the two vectors without any error tolerance.
     * @param Other vector for comparision. */
    inline bool is_exact_equal(const V &other) inline const
    {
        return (X == other.X) && (Y == other.Y) && (Z == other.Z);
    }

//————————————————————————————————— PRODUCT FUNCTIONS ————————————————————————————————————————————————

    //31
    /**
     * @brief Calculates the dot product between two vectors. Returns 0 if the answer is below tolerance.
     * @return double dot product. */
    inline double dot(const V& other, double tolerance = V::tolerance) const
    {
        double temp = X * other.X + Y * other.Y + Z * other.Z;
        
        return (std::fabs(temp) > tolerance) * temp + 0;
    }

    
    //32
    /**
     * @brief Fast Dot - Calculates the dot product between two vectors. Does not check for tolerance.
     * @return return double dot product. */
    inline double fdot(const V& other) const
    {
        return (X * other.X + Y * other.Y + Z * other.Z);
    }


    //33
    /**
     * @Cross product of two vectors.
     * @param B V other vector.
     * @return cross product - V. */
    inline V cross(const V& other) const 
    {
        return V(Y * other.Z - Z * other.Y,
                -X * other.Z + Z * other.X,
                 X * other.Y - Y * other.X);
    }



    //34
    /**
     * @brief Returns the scalar triple product of three vectors.
     * @param Three vectors - A, B, C; such that - Scalar triple product: A • (B x C) . */
    inline double scalar_tri_prod(const V &A, const V &B, const V &C)
    {
        V tmp = B.cross(C);
        return A.dot(tmp);
    }

    //35
    /**
     * @brief Returns the vector triple product of three vectors.
     * @param Three vectors - A, B, C; such that - Vector triple product: A x (B x C) = B * A•C - C * A•B . */
    inline V vector_tri_prod(const V &A, const V &B, const V &C)
    {
        return (B * (A.dot(C))) - (C * A.dot(B));
    }


    //36
    //! Calculates the angle between two vectors using arccos(cos inverse).
    inline double arccos(const V& other) const
    {
        double temp = this->fdot();
        return std::acos(temp);
    } //End of arccos()


   //37
   /**
    * @brief Alias of V::arccos().  Calculates the angle between two vectors using arccos(cos inverse). */
    double __attribute__((always_inline)) radian_angle(const V& other) const
   {
        return arccos(other);
   }

//————————————————————————————————— SCALAR FUNCTIONS ————————————————————————————————————————————————

    //38
    /**
     * @brief Converts the object to a scalar mask. The scalar is assigned to the X component. The rest of the components are set to zero.
     * @param Scalar to set. */
    void make_scalar(double scalar) inline
    {
        this->x = scalar;
        this->y = 0.0;
        this->z = 0.0;
    }

    //39
    /**
     * @brief Retrives the scalar value. It is just a wrapper around V::x() accessor. */
    double get_scalar() inline const
    {
        return this->x();
    }

//————————————————————————————————— TEMPLATES & GENERATORS ————————————————————————————————————————————————

    //40
    /**
     * @brief Initalize x, y, and z with the given arguement.
     * @param double initalization value. */
    void set_xyz(double val)  inline const
    {
        X = val; Y = val; Z = val; 
    }

    //41
    /**
     * @brief Initalize x and y with the given arguement, and z is set to zero.
     * @param double initalization value. */
    void set_xy(double val) inline const
    {
        X = val; Y = val; Z = 0.0;
    }

    //42
    /**
     * @brief Initalize y and z with the given arguement, and z is set to zero.
     * @param double initalization value. */
    void set_yz(double val) inline const
    {
        X = 0.0; Y = val; Z = val;
    }

    //43
    /**
     * @brief  Initalize z and x with the given arguement, and z is set to zero.
     * @param double initalization value. */
    void set_zx(double val) inline const
    {
        X = val; Y = 0.0; Z = val;   
    }

    //44
    /**
     * @brief Generates a vector by adding a given scalar and generating function.
     * @param Scalar common to all the three components.
     * @param A generating function pointer that returns a double.
     * @return A generated 3 vector V.*/
    V generate_add(double scalar = 0, double (*generating_fn)())
    {
        return V(scalar + generating_fn(), scalar + generating_fn(), scalar + generating_fn());
    }

    //45
    /**
     * @brief Generates a vector by multiplying a given scalar and generating function.
     * @param Scalar common to all the three components.
     * @param A generating function pointer that returns a double.
     * @return A generated 3 vector V.*/
    V generate_mul(double scalar = 0, double (*generating_fn)())
    {
        return V(scalar * generating_fn(), scalar * generating_fn(), scalar * generating_fn());
    }

    //46
    /**
     * @brief Generates a vector by operating the scalar and the generating function with the passed operation (3rd parameter).
     * @param Scalar common to all the three components.
     * @param A generating function pointer that returns a double.
     * @param The operation to perform between the scalar and the vector.
     * @return A generated 3 vector V.*/
    V generate(double scalar = 0, double (*generating_fn)(), double (*op)(double, double))
    {
        return V( op(scalar, generating_fn()),
                  op(scalar, generating_fn()),
                  op(scalar, generating_fn()) );
    }



    //47
    /**
     * @brief Generates a vector by operating the two generating functions with the passed operation (3rd parameter).
     * @param First generating function pointer that returns a double.
     * @param Second generating function pointer that returns a double.
     * @param The operation to perform between the two generating functions.
     * @return A generated 3 vector V.*/
    V generate( double(*generating_fn1)(), double (*generating_fn2)(), double (*op)(double, double))
    {
        return V( op(generating_fn1(), generating_fn2()),
                  op(generating_fn1(), generating_fn2()),
                  op(generating_fn1(), generating_fn2()) );
    }



//—————————————————————————————— Non-const Operators (mutators) ————————————————————————————————————————

    //48
    /**
     * @brief Assignment operator
     * @param V other - assigner. */
    void operator= (const V& other) inline
    {
        X = other.X;
        Y = other.Y;
        Z = other.Z;
    }

    //49
    /**
     * @brief  Iterative Subtraction Operator.
     * @param V other vector. */
    void operator-= (const V& other) inline 
    {
        X -= other.X;
        Y -= other.Y;
        Z -= other.Z;
    }


    //50
    /**
     * @brief  Iterative Addition Operator.
     * @param V other vector. */
    inline void operator+= (const V& other) 
    {
            X += other.X;
            Y += other.Y;
            Z += other.Z;
    }


    //51
    /**
     * @param Iterative Scaling Operator - double.
     * @param scale double -> template. */
    template <typename T>
    inline void operator*= (T scale) //MUTATOR
    {
        X *= scale; Y *= scale; Z *= scale;
    }


    //52
    /**
     * @brief Iterative Scaling Division Operator - double.
     * @param scale double -> template. */
    template <typename T>
    inline void operator/= (T scale)
    {
        X /= scale; Y /= scale; Z /= scale;
    }



//——————————————————————————————  ————————————————————————————————————————


   /*
    Include code for int functors, double functors and, mt_19937 generator

   */
    //41 - RND
    //! Orthagonalises the given vector to the other vector randomly. Uses random number generation.
    /*
    \param B V other vector
    \param usig_RND_db -> Function pointer that generates unsigned double values
    \return vector A which is orthogonalized
    */
    inline static V rnd_orthogonal_of(const V& B, double(*usig_RND_db)()) 
    {
        //extern int Rndm(int, int);
        V temp;
        temp.x = (0.5 - usig_RND_db())*2;
        temp.y = (0.5 - usig_RND_db())*2;
        temp.z = (-1*temp.x * B.x - temp.y * B.y) / B.z;
        temp.normalise();
        return temp;
    } //End of rnd_orthogonal_of()

    //42 - RND
    //! Orthagonalises the given vector to the other vector randomly. Uses random number generation.
    /*
    \param B V other vector
    \param ar,br -> Two random unsigned double values
    \return vector A which is orthogonalized
    */
    inline static V rnd_orthogonal_of(const V& B, double ar, double br) 
    {
        //extern int Rndm(int, int);
        V temp;
        temp.x = (0.5 - ar)*2;
        temp.y = (0.5 - br)*2;
        temp.z = (-1*temp.x * B.x - temp.y * B.y) / B.z;
        temp.normalise();
        return temp;
    } //End of rnd_orthogonal_of()


    //43 - RND
    //! Orthagonalises the given vector to the other vector randomly. Uses random number generation.
    /*
    \param B V other vector
    \param rndengine -> std::mt19937 object by reference
    \return vector A which is orthogonalized
    */
    inline static V rnd_orthogonal_of(const V& B, std::mt19937 &rndengine) 
    {
        //extern int Rndm(int, int);
        
        V temp;
        temp.x = (0.5 - double(rndengine())/double(rndengine.max()))*2;
        temp.y = (0.5 - double(rndengine())/double(rndengine.max()))*2;
        temp.z = (-1*temp.x * B.x - temp.y * B.y) / B.z;
        temp.normalise();
        return temp;
    } //End of rnd_orthogonal_of()
//**********************************************************************************

    //44 - RND
    //! Converts the given vector to a random unit vector. Uses a function.
    /*
    \param std::mt19937 engine by reference
    */    
    inline void rndUnit(std::mt19937 &rndengine)
    {
        this->x = (0.5 - double(rndengine())/double(rndengine.max()))*2;
        this->y = (0.5 - double(rndengine())/double(rndengine.max()))*2;
        this->z = (0.5 - double(rndengine())/double(rndengine.max()))*2;
        this->normalise();
    } //End of rndUnit()


    //45
    //! Converts the given vector to a random unit vector. Uses a function.
    /*
    \param function pointer usig_RND_db that gives unsigned double values
    */
    inline void rndUnit(double(*usig_RND_db)())
    {
        this->x = (0.5 - usig_RND_db())*2;
        this->y = (0.5 - usig_RND_db())*2;
        this->z = (0.5 - usig_RND_db())*2;
        this->normalise();
    } //End of rndUnit()


    //46
    //! Converts the given vector to a random unit vector. Uses a function.
    /*
    \param a,b,c double values that are unsigned and randomly generated
    */
    inline void rndUnit(double a, double b, double c)
    {
        this->x = (0.5 - a)*2;
        this->y = (0.5 - b)*2;
        this->z = (0.5 - c)*2;
        this->normalise();
    } //End of rndUnit()
//**********************************************************************************
   

///><><><><><><><><><><><>< ROTATION FUNCTIONS ><><><><><><><><><><><><><><><><><><><

    //47
    //! Axis-Angle rotation of vector.
    /*
    \param &Axis V unit axis.
    */
    inline void rotate(const V& Axis, double angle) 
    {
        V axis = Axis; //If Axis is not a unit vector safety.
        if(!axis.isUnit())
        {
            axis.normalise();
        }

        double c,s,c_1;
        c = cos(angle); s = sin(angle); c_1 = 1-c;

        double t1 =  axis.x * axis.x * c_1 + c;
        double t2 =  axis.x * axis.y * c_1 - axis.z * s;
        double t3 =  axis.x * axis.z * c_1 + axis.y * s;
        double t4 =  axis.y * axis.x * c_1 + axis.z * s;
        double t5 =  axis.y * axis.y * c_1 + c;
        double t6 =  axis.y * axis.z * c_1 - axis.x * s;
        double t7 =  axis.z * axis.x * c_1 - axis.y * s;
        double t8 =  axis.z * axis.y * c_1 + axis.x * s;
        double t9 =  axis.z * axis.z * c_1 + c;
            
        double newx = t1*this->x + t2*this->y + t3*this->z;
        double newy = t4*this->x + t5*this->y + t6*this->z;
        double newz = t7*this->x + t8*this->y + t9*this->z;

        this->x = newx;
        this->y = newy;
        this->z = newz;


    } //End of rotate()
//**********************************************************************************
    
    //48
    //! Quarternion rotation method. Overloaded rotate().
    /*
    \param &q Quarternion Q provided for rotation.
    */
    inline void rotate(Q &q) 
    {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;

            //    t1 = quat.w * quat.w;
            t2 =  q.a * q.b;
            t3 =  q.a * q.c;
            t4 =  q.a * q.d;
            t5 = -q.b * q.b;
            t6 =  q.b * q.c;
            t7 =  q.b * q.d;
            t8 = -q.c * q.c;
            t9 =  q.c * q.d;
            t10 = -q.d * q.d;

            newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
            newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
            newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

            x = newx;
            y = newy;
            z = newz;
    } //End of rotate() Overloaded
//**********************************************************************************

//END
///><><><><><><><><><><><>< ROTATION FUNCTIONS ><><><><><><><><><><><><><><><><><><><





    //51
    //! Returns a unit vector that points from the initial to the terminal point. Point vector = (terminal point) - (initial point) . {{ initial point --> terminal point }}
    /*
    \param other V other vector - terminal point.
    \return point unit vector
    */
    inline V points_to(V const &other) const //say "points from 'this->' to 'other'."
    {
            V temp;
            temp.x = other.x - x;
            temp.y = other.y - y;
            temp.z = other.z - z;
            temp.normalise();
            return temp;
    } //End of point()
//**********************************************************************************
    //52
    //! Reurns the largest component of the three. 0:x, 1:y, 2:z. If there is a tie, the first element is returned.
    /*
    \param &vec V vector
    \return unsigned int component
    */
/*    inline unsigned int max_component(V const &vec)
    {
        int max = 0;
            
        if(vec.y > vec.x)
             max = 1;
        if(vec.z > vec[max])
            max = 2;

        return max;
    } //End of max_component()
//**********************************************************************************

    //53
    //! Reurns the largest component of the three. 0:x, 1:y, 2:z. If there is a tie, the first element is returned.
    /*
    \param &vec V vector
    \return unsigned int component
    */
/*    inline unsigned int min_component(V const &vec)
    {
        int min = 0;
            
        if(vec.y < vec.x)
             min = 1;
        if(vec.z < vec[min])
            min = 2;

        return min;
    } //End of max_component()*/
//**********************************************************************************

    //54
    //! Returns the projection of vec on the plane with norm planenorm.
    /*
    \param &vec V vector
    \param &planenorm V normal to the given plane
    */
    V static plane_projection(const V &vec, const V &planenorm) //Vec must be terminating on the plane.
    {
        V tempplanenorm = planenorm;
        V temp = tempplanenorm*(planenorm.dot(vec)); //Find projection of vec with the normal
        return V(vec - temp); //return the other component

    } //End of plane_projection()

    //55
    //! Returns the projection of vec on the plane defined by the other two vectors.
    /*
    \param &vec V vector
    \param &planevec1, &planevec2 are vectors that reside in the said plane
    */
    V static plane_projection(const V &vec, const V &planevec1, const V &planevec2)
    {
        V planenorm = planevec1.cross(planevec2); //Find the normal using the two vectors
        V temp = planenorm*planenorm.dot(vec); //Find projection of vec with the normal
        return V(vec - temp); //return the other component

    } //End of plane_projection() Overloaded


    //56
    //! Returns the projection of vec parallel to the given plane.
    /*
    \param &vec V vector
    \param &planenorm V normal to the given plane
    */
    V static plane_parallelprojection(const V &vec, const V &planenorm) //Vec must be terminating on the plane.
    {
        V tempplanenorm = planenorm;
        V temp = tempplanenorm*(planenorm.dot(vec)); //Find projection of vec with the normal
        return temp; //return the projection

    } //End of plane_parallelprojection()



    //57
    //! Returns the size of the corresponding linesegment
    /*
    \param other V& other
    \return line segment length as double
    */
    inline double segment_len(const V &other) const
    {
        V temp(x - other.x, y - other.y, z - other.z);
        return temp.size();
    }






//------------------------------------        Scalar functions    ----------------------------------



//--------------------------------------- Component Mutator Functions ------------------------------   
    inline void comp_divide(const V &other)
    {
        this->x /= other.x;
        this->y /= other.y;
        this->z /= other.z;
    }


    inline void comp_square()
    {
        this->x = this->x*this->x;
        this->y = this->y*this->y;
        this->z = this->z*this->z;
    }

    //How to represent Infinitesimal Vector?

    inline double accumulate() const
    {
        return this->x + this->y + this->z;
    }

    inline V square() const
    {
        return V(this->x*this->x, this->y*this->y, this->z*this->z)
    }

    //Multi Vector Operations


    //Friend Declarations
    friend std::ostream& operator<< (std::ostream &stream, const V &vec); //61

}; //end of class V
//**********************************************************************************
 
//61 
//! Stream operator overload that puts V::info() to the ostream.
std::ostream& operator<< (std::ostream &stream, const V &vect)
{
    stream << vect.info();
    return stream;
} //End of friend overload operator<<

////////////////////////////////END OF TRANSLATION UNIT/////////////////////////


/*   N   O   T   E   S
    
    • Typically, we won’t be able to use a member overload if the left operand is either not a class (e.g. int), or it is a class that we can’t modify (e.g. std::ostream).

    • The normal or friend function version has the added benefit of “symmetry”, as all operands become explicit parameters (instead of the left operand becoming *this and the right operand becoming an explicit parameter).

*/


//////|||||||||||||||||||||||||| NEW UPCOMING FEATURES|||||||||||||

  /*  bool linearly_dependent(const V &A, const V &B, const V &C)
    {
        double x_indices[3];
        double y_indices[3];
        double z_indices[3];

        double coeff[3] = {0.0};        
    }*/

    /*Orthonormal_Set()
    {
        //Inner product is dirac_delta()


    }


    Make Variadic
    Gen_Matrix_3XN<N>()
    {

    }

    Gen_Matrix_NX3<N>()
    {
        
    }

    Direct_Sum(const V &A, const V &B)

    Direct_Product(const V &A, const V &B)


    //Derivative and Integration si wrt to all the elements of matrix individually.

    //***
    //Returns the inverse of vector such that {{ this->V * returned V = I (Identity) }}
    V const Inverse() //==> 
    {
        V temp = ()
        
    }*/

//////|||||||||||||||||||||||||| NEW UPCOMING FEATURES|||||||||||||