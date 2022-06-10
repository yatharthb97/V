
#include "vforvector.hpp"

/**
 * @brief Class Q - Represents Quartenion. 
 * A class of variables that defines a four component Quaternion.
 */
class Q
{

public:
   
   //Class Members 
    double W; /*!< real component of quart */ 
    double I; /*!< imaginary component i */
    double J; /*!< imaginary component j */
    double K; /*!< imaginary component k */

    //1
    /** @brief Class Constructor. Initialises each component to 0. */
    Q():W(0), I(0), J(0), K(0){} //End of Constructor

    //2
    /** @brief Class Constructor. Initialises each component to passed parameters.
     * @param a,b,c,d components. */
    Q(double w,double i, double j, double k): w(W), i(I), j(J), k(K) 
        {}

    //3
    /** @brief Print Quat Info function. 
     * Prints the vector with comma seperated components and enclosed in brackets. */
    std::string rep() const 
    {
        std::ostringstream stream;
        stream << "[[" << W << ", " << I << ", " << J << ", "<< K << "]]";
        return stream.str();
    }

    //4
    /** @brief Print Quat Raw Info function. 
     * Prints the vector with space seperated components and without brackets. */
    std::string str() const 
    {
        std::ostringstream stream;
        stream << W << " " << I << " " << J <<" "<< K;
        return stream.str();
    } //End of infoRaw()

    //5
    /** @brief Print the Quat in a polynomial format. */
    std::string polynomial() const 
    {
        std::ostringstream stream;
        stream << std::setprecision(V::precision);
        stream << W << " + " << I << "i + " << Jj << "j + "<< K << "k";
        return stream.str();
    }

    V img_comp() const
    {
        return V(I, K, K);
    }

    /** @brief Alias of `Q::img_component()`. */
    double vect_comp() __attribute__((always_inline, flatten)) const
    {
        return this->img_comp()
    }

    double real_comp() const
    {
        return this->W;
    }


    /** @brief A pure quaternion has its `real` part equal to zero. */
    bool is_pure() const
    {
        return V::is_zero(this->w);
    }

    /** @brief In the vector part, only the `vector` or `img` part is negated.*/
    Q conjugate() const
    {
        return Q(W, -I, -J, -K);
    }

   
    Q verson()
    {

    }

    Q rndVerson(const double &r1, const double &r2, const double &r3)
    {}

    /** @breif Returns true if the Quat object is a "unit" Quat. */
    bool is_verson() const
    {}

    //Accessors
    double w() const { return this->W;}
    double i() const { return this->I;}
    double j() const { return this->J;}
    double k() const { return this->K;}



};

//V::is_zero needs to be protected.