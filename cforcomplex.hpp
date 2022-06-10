


/** @brief A very light class that represents a complex number. */


class C
{
public:


	double X;
	double I;

	C(): X(0), I(0)
	{}

	C(double x, double i): x(X), i(I)
	{}

	std::string rep() const
	{
		std::ostringstream stream;
		stream << std::setprecision(V::precision);
		stream << "(" << X << ", " << I << ")";
	}
	
	std::string str() const 
	{
		std::ostringstream stream;
		stream << std::setprecision(V::precision);
		stream << X << " " << I;
	}

	std::string expand() const
	{
		std::ostringstream stream;
		stream << std::setprecision(V::precision);
		stream << << X << "+ i" << I;
	}


	double x() const { return this->X; }
	double i() const { return this->I; }

	double img() const { return i(); }
	double real() const { return x(); }

};