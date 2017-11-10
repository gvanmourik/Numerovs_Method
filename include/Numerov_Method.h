#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>

class Numerov_Method {
	
private:
	bool parabolic, wedge;
	double a;								// Halfwidth of the well (fm)
	double V_o;								// Potential in the well (MeV)
	double c;								// 2*m/h_bar^2
	double E;								// Constant energy (MeV)
	double h;								// Step size
	double range;							// User_defined range (fm)
	double epsilon;
	double initial_psi_0;
	double initial_psi_1;
	double numberOfPoints;					// User-defined number of points
	std::vector<double> x_left;				// x-points left side of the well
	std::vector<double> x_right;			// x-points right side of the well
	std::vector<double> psi_left;			// psi for foward integration
	std::vector<double> psi_right;			// psi for backwards integration


public:

	// Constructors and destructor
	Numerov_Method(): a(2.0), V_o(-83.0), c(0.04819), range(4.0), numberOfPoints(500), parabolic(true), wedge(false)
	{	
		initial_psi_0 = 1.0e-10;
		initial_psi_1 = 1.0e-9;
		populate_x();
		populate_psi();
	}
	Numerov_Method(double user_range, double user_numberOfPoints, std::string type): a(2.0), V_o(-83.0), c(0.04819) 
	{
		initial_psi_0 = 1.0e-10;
		initial_psi_1 = 1.0e-9;
		range = user_range;
		numberOfPoints = user_numberOfPoints;

		set_potential_type(type);
		populate_x();
		populate_psi();
	}
	~Numerov_Method(){}


	/**
	 * @brief      Sets the potential type.
	 *
	 * @param[in]  type  The type for the potential
	 */
	void set_potential_type(std::string type)
	{
		if (type == "parabolic")
		{
			parabolic = true;
			wedge = false;
		}
		else if (type == "wedge")
		{
			parabolic = false;
			wedge = true;
		}
	}

	/**
	 * @brief      Function of energy
	 *
	 * @param[in]  energy  The energy
	 * @param[in]  even    Denotes if the returned value corresponds to even solns
	 *
	 * @return     The difference between the left and right gradients
	 */
	double f_E_(double energy, bool even)
	{
		if ( !parabolic && !wedge )
		{
			printf("Potential type has NOT been set...try again\n");
			return 0.0;
		}

		E = energy;
		left_integration();
		right_integration(even);
		return gradient_difference();
	}

	/**
	 * @brief      Integrate with Numerov's fron the left
	 */
	void left_integration()
	{
		int max = x_left.size();

		// Explicit step calculation. Should already be set.
		find_h();

		// Add two initial values based on psi described in report
		psi_left[0] = initial_psi_0;
		psi_left[1] = initial_psi_1;

		// printf("in left_integration()...\n");
		for (int i = 2; i < max; i++)
		{
			psi_left[i] = psi_next(i, h, "left");
		}
	}

	/**
	 * @brief      Integrate with Numerov's fron the right
	 *
	 * @param[in]  isEven  Denotes if the returned value corresponds to even solns
	 */
	void right_integration(bool isEven)
	{
		int max = x_right.size();

		// Explicit step calculation. Should already be set.
		find_h();

		// Add two initial values based on psi described in report
		if ( isEven )
		{
			psi_right[0] = initial_psi_0;
			psi_right[1] = initial_psi_1;
		}
		// Negative initial values needed for odd solution
		else
		{
			psi_right[0] = -initial_psi_0;
			psi_right[1] = -initial_psi_1;
		}

		for (int i = 2; i < max; i++)
		{
			psi_right[i] = psi_next(i, h, "right");
		}
	}

	/**
	 * @brief      Numerov's method to find the next psi value
	 *
	 * @param[in]  n     The index of next
	 * @param[in]  step  The step
	 * @param[in]  type  Specifies if left or right integration
	 *
	 * @return     The next value for psi
	 */
	double psi_next(int n, double step, std::string type)
	{
		double v1, v2, v3;
		double q = 1.0/12.0;

		if ( type == "left")
		{
			v1 = 2.0 * psi_left[n-1] * (1.0 - (5 * q) * pow(step,2) * k_x(x_left[n-1]) );
			v2 = psi_left[n-2] * (1.0 + q * pow(step,2) * k_x(x_left[n-2]) );
			v3 = 1.0 + ( q * pow(step,2) * k_x(x_left[n]) );
		}
		if ( type == "right" )
		{
			v1 = 2.0 * psi_right[n-1] * (1.0 - (5 * q) * pow(step,2) * k_x(x_right[n-1]) );
			v2 = psi_right[n-2] * (1.0 + q * pow(step,2) * k_x(x_right[n-2]) );
			v3 = 1.0 + ( q * pow(step,2) * k_x(x_right[n]) );
		}

		return (v1 - v2) / v3;
	}

	/**
	 * @brief      k(x)
	 *
	 * @param[in]  x     The position
	 *
	 * @return     A value based on the position
	 */
	double k_x(double x)
	{
		return c * ( E - V_x(x) );
	}

	/**
	 * @brief      Potential function V(x)
	 *
	 * @param[in]  x     The position
	 *
	 * @return     The value of the potential based on the position
	 */
	double V_x(double x)
	{
		if (parabolic)
		{
			if ( std::abs(x) > a )
			{
				return 0.0;
			}
			else
			{
				return V_o * ( pow(x, 2) - pow( a,2) ) / pow(a, 2);
			}
		}

		if (wedge)
		{
			if ( -a <= x && x <= 0.0 )
			{
				return -V_o * (x + a) / a;
			}
			if ( 0.0 <= x && x <= a )
			{
				return V_o * (x - a) / a;
			}
			return 0.0;
		}

		return 0.0;
	}

	/**
	 * @brief      Evenly populate the x-point vector
	 */
	void populate_x()
	{
		find_h();
		int length = numberOfPoints/2;

		for (int i = 0; i < length+1; i++)
		{
			x_left.push_back( (h * i) - range );
			x_right.push_back( range - (h * i) );
		}
	}

	/**
	 * @brief      Populate left and right psi to work with indicies
	 */
	void populate_psi()
	{
		int length = numberOfPoints/2;
		double initial_value = 0.0;
		for (int i = 0; i < length+1; ++i)
		{
			psi_left.push_back(initial_value);
			psi_right.push_back(initial_value);
		}
	}

	/**
	 * @brief      Find the step based on the range and the number of points
	 */
	void find_h()
	{
		h = (range * 2) / (numberOfPoints - 1);
	}

	/**
	 * @brief      Find the difference in the left and right gradients
	 *
	 * @return     The difference between the gradients
	 */
	double gradient_difference()
	{
		int zero = psi_left.size() - 1;

		double left_gradient = (psi_left[zero-1] - psi_left[zero]) / h;
		double right_gradient = (psi_right[zero-1] - psi_right[zero]) / -h;

		return left_gradient - right_gradient + psi_left[zero] - psi_right[zero];
	}

	void plot_wavefunction(double energy, bool even, std::string label)
	{
		f_E_(energy, even);
		export_vector_to_file(x_left, psi_left, "psi_left_" + label + ".data");
		export_vector_to_file(x_right, psi_right, "psi_right_" + label + ".data");
	}

	/**
 	 * @brief      Export two vectors to a data file
  	 *
 	 * @param[in]  x      A vector of the x-points
 	 * @param[in]  y      A vector of the y-points
 	 * @param[in]  fname  The file name
 	 */
	void export_vector_to_file(std::vector<double> x, std::vector<double> y, std::string fname)
	{
		std::ofstream file("data/" + fname);
		// std::copy(x.begin(), x.end(), std::ostream_iterator<double>(file));
		// std::copy(y.begin(), y.end(), std::ostream_iterator<double>(file, "\n"));	

		if ( file.is_open() )
	    {
   		 	int length = y.size();
       		for (int i = 0; i < length; i++)
       		{
         		file << x[i] << " " << y[i] << "\n";
        	}
    	}
    	else
    	{
        	perror("ERROR Unable to export to CSV");
    	}
	}

};