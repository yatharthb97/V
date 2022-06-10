

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
    inline bool is_exactly_equal(const V &other) inline const
    {
        return (X == other.X) && (Y == other.Y) && (Z == other.Z);
    }

    //31
    /**
     * @brief Returns whether the given vector is parallel to the other vector with some tolerance. Equivalent to cross product is zero. */
    bool is_parallel(const V &other, double tolerance = V::tolerance) const
    {
        return std::fabs(this->cross(other) - 1.0) <= tolerance;
    }













//—————————————————————————————— LINE FUNCTIONS ————————————————————————————————————————


    //54
    /**
     * @brief  Returns a vector that points from "this" vector to the "other" vector. Vector = (terminal point) - (initial point) . {{ this --> other }}.
     * @param other V other vector - terminal point.
     * @return point vector. */
    V points_to(const V &other) inline const
    {
            V temp;
            temp.X = other.X - X;
            temp.Y = other.Y - Y;
            temp.Z = other.Z - Z;
            temp.normalise();
            return temp;
    }


    //55
    /**
     * @brief  Returnes vector originates at the "other vector" and points to "this" vector. {{ Other --> This}}. The behaviour is the same as V::points_to(&V) but with the negative of its result.
     * @param other V& other vector - origin point.
     * @return point unit vector. */
    V __attribute__((always-inline)) point_from(const V &other) const
    {
        return this->points_to().neg();
    }


    //56
    /**
     * @brief  Unit vector version of V::points_to(). Returns a unit vector that points from "this" vector to the "other" vector. Vector = (terminal point) - (origin point) . {{ this --> other }}.
     * @param other V other vector - terminal point.
     * @return point unit vector. */
    V __attribute__((always-inline)) unit_points_to(const V &other) const
    {
            return this->points_to().normalise();
    }


    //57
    /**
     * @brief Unit vector version of V::points_from(). Returns unit vector originates at the "other vector" and points to "this" vector. {{ Other --> This}}. The behaviour is the same as V::points_to(&V) but with the negative of its result.
     * @param other V& other vector - origin point.
     * @return point unit vector. */
    V __attribute__((always-inline)) unit_point_from(const V &other) const
    {
        return this->unit_points_to().neg();
    }


    //58
    /**
     * @brief Returns the size of the corresponding linesegment.
     * @param other V& other vector.
     * @return line segment length as double. */
    inline double segment_len(const V &other) const
    {
        V temp(x - other.x, y - other.y, z - other.z);
        return temp.size();
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
    inline void rotate(const V& Axis, const double angle) 
    {
        V axis = Axis; //If Axis is not a unit vector safety.
        
        #ifdef __VFORVECTOR_SAFETY__
        if(!axis.is_unit())
        {
            axis.normalise();
        }
        #endif

        double c,s,c_1;
        c = std::cos(angle); s = std::sin(angle); c_1 = 1-c;

        double t1 =  axis.x() * axis.x() * c_1 + c;
        double t2 =  axis.x() * axis.y() * c_1 - axis.z() * s;
        double t3 =  axis.x() * axis.z() * c_1 + axis.y() * s;
        double t4 =  axis.y() * axis.x() * c_1 + axis.z() * s;
        double t5 =  axis.y() * axis.y() * c_1 + c;
        double t6 =  axis.y() * axis.z() * c_1 - axis.x() * s;
        double t7 =  axis.z() * axis.x() * c_1 - axis.y() * s;
        double t8 =  axis.z() * axis.y() * c_1 + axis.x() * s;
        double t9 =  axis.z() * axis.z() * c_1 + c;
            
        double newx = t1*x() + t2*y() + t3*z();
        double newy = t4*x() + t5*y() + t6*z();
        double newz = t7*x() + t8*y() + t9*z();

        this->X = newx;
        this->Y = newy;
        this->Z = newz;


    } //End of rotate()
//**********************************************************************************
    
    //48
    //! Quarternion rotation method. Overloaded rotate().
    /*
    \param &q Quarternion Q provided for rotation.
    */
    inline void rotate(const Q &q) 
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

            newx = 2.0 * ( (t8+t10) * X + (t6-t4)  * Y + (t3+t7) * z ) + X;
            newy = 2.0 * ( (t4+t6)  * X + (t5+t10) * Y + (t9-t2) * z ) + Y;
            newz = 2.0 * ( (t7-t3)  * X + (t2+t9)  * Y + (t5+t8) * z ) + Z;

            X = newx;
            Y = newy;
            Z = newz;
    } //End of rotate() Overloaded
//**********************************************************************************



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

    //
    /**
     * @brief Returns the projection of this vector on the other vector.
     * @param V other vector on which the projection is computed. */
    V projection_on(const V &other) inline const
    {
        V temp = other * (this->fdot(other));
        return temp;
    }

    //
    /**
     * @brief  Returns the projection of vec on a plane described by its normal "planenorm". The function assumes that this vector originates on the given plane.
     * @param V normal describing the plane. */
    V plane_projection(const V &planenorm) inline const
    {
        V temp = this->projection_on(planenorm); //Find projection of vec with the normal
        return V(this - temp); //return the other component

    }

    //
    /**
     * @brief Returns the projection of this vector on the plane defined by the other two vectors.
     * @param Two vectors that are coplaner and specifies the plane. The two vectors must not be parallel(no safety included). */
    V  plane_projection(const V &planevec1, const V &planevec2) inline const
    {
        V planenorm = planevec1.cross(planevec2); //Find the normal using the two vectors
        V temp = this->projection_on(planenorm); //Find projection of vec with the normal
        return V(this - temp); //return the other component

    } //End of plane_projection() Overloaded



    std::tuple<V, V> plane_proj_pair(const V &planenorm) inline const
    {
        V normal_proj = this->projection_on(planenorm); //Find projection of vec with the normal
        V plane_proj = V(this - temp);
        
        return std::make_tuple(plane_proj, normal_proj); //return the other component
    }


    std::tuple<V, V> plane_proj_pair(const V &planevec1, const V &planevec2) inline const
    {
        V planenorm = planevec1.cross(planevec2);
        return plane_proj_pair(planenorm);
    }



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



    //Multi Vector Operations

}; //end of class V
//**********************************************************************************
 

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

//Macro Definaions
#define __VFORVECTOR_DEFAULT_TOLERANCE__ 1E-10 //! Defines the default tolerance value for the class.
