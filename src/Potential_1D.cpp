#include <Numerov_Method.h>

// Declarations
std::vector<double> generate_energies(double E_range, double stepSize);
std::vector<double> generate_f_E(bool evenSolutions, std::vector<double> energies, std::string type);
void export_vector_to_file(std::vector<double> x, std::vector<double> y, std::string fname);
double bisection_method(double a, double b, double epsilon, bool even, std::string type);
bool acceptable_error(double x_1, double x_2, double epsilon);
int sign(double x_1);


int main( int argc, char* argv[] )
{
	// Variables
	std::vector<double> energies, f_E_even, f_E_odd;
	Numerov_Method parabolic_potential, wedge_potential;
	double root1_p, root2_p, root3_p, root4_p, root1_w, root2_w, root3_w, root4_w;
	double epsilon = 1.0e-5;
	bool even = true;
	bool odd = false;
	std::string p = "parabolic";
	std::string w = "wedge";

	// Set potential types
	parabolic_potential.set_potential_type(p);
	wedge_potential.set_potential_type(w);

	// Generate f as a function of E
	/* UNCOMMENT BELOW SECTION TO FIND APPX ROOTS */
	// Parabolic energies
	energies = generate_energies(85.0, 0.1);
	f_E_even = generate_f_E(even, energies, p);
	f_E_odd = generate_f_E(odd, energies, p);
	parabolic_potential.export_vector_to_file(energies, f_E_even, "f_E_even_p.data");
	parabolic_potential.export_vector_to_file(energies, f_E_odd, "f_E_odd_p.data");

	// Wedge energies
	f_E_even = generate_f_E(even, energies, w);
	f_E_odd = generate_f_E(odd, energies, w);
	wedge_potential.export_vector_to_file(energies, f_E_even, "f_E_even_w.data");
	wedge_potential.export_vector_to_file(energies, f_E_odd, "f_E_odd_w.data");

	// From observation of f(E) use the following energies
	root1_p = bisection_method(23.7, 22.7, epsilon, even, p);
	root2_p = bisection_method(72.7, 71.7, epsilon, even, p);
	root3_p = bisection_method(23.7, 22.7, epsilon, odd, p);
	root4_p = bisection_method(81.7, 79.7, epsilon, odd, p);

	root1_w = bisection_method(20.5, 19.5, epsilon, even, w);
	root2_w = bisection_method(63.3, 60.8, epsilon, even, w);
	root3_w = bisection_method(20.8, 20.5, epsilon, odd, w);
	root4_w = bisection_method(72.8, 69.6, epsilon, odd, w);

	// Print out root values for user
	printf("Parabolic root1 (even) = %f MeV\n", root1_p);
	printf("Parabolic root2 (even) = %f MeV\n", root2_p);
	printf("Parabolic root3 (odd)  = %f MeV\n", root3_p);
	printf("Parabolic root4 (odd)  = %f MeV\n\n", root4_p);
	printf("Wedge root1 (even) = %f MeV\n", root1_w);
	printf("Wedge root2 (even) = %f MeV\n", root2_w);
	printf("Wedge root3 (odd)  = %f MeV\n", root3_w);
	printf("Wedge root4 (odd)  = %f MeV\n", root4_w);

	// printf("root3 = %f MeV\n", root3);

	// Show progress towards convergence for parabolic potential
	parabolic_potential.plot_wavefunction(35.0, even, "35");
	parabolic_potential.plot_wavefunction(29.0, even, "29");
	parabolic_potential.plot_wavefunction(23.0, even, "23");

	// // Show converged solutions
	parabolic_potential.plot_wavefunction(root1_p, even, "root1_p");
	parabolic_potential.plot_wavefunction(root2_p, even, "root2_p");
	parabolic_potential.plot_wavefunction(root3_p, odd, "root3_p");
	parabolic_potential.plot_wavefunction(root4_p, odd, "root4_p");
	wedge_potential.plot_wavefunction(root1_w, even, "root1_w");
	wedge_potential.plot_wavefunction(root2_w, even, "root2_w");
	wedge_potential.plot_wavefunction(root3_w, odd, "root3_w");
	wedge_potential.plot_wavefunction(root4_w, odd, "root4_w");


	return 0;
}


/**
 * @brief      Bisection Method - Find a root within the range a-->b
 *
 * @param[in]  a     Min value of the range 
 * @param[in]  b     Max value of the range
 *
 * @return     Energy in MeV corresponding to the root found
 */
double bisection_method(double a, double b, double epsilon, bool even, std::string type)
{
	Numerov_Method potential;
	potential.set_potential_type(type);
	double c, result_a, result_c;

	while ( 1 )
	{
		c = (a+b)/2;
		result_c = potential.f_E_(c, even);
		result_a = potential.f_E_(a, even);

		if ( result_c == 0 || acceptable_error(a, b, epsilon) )
		{
			return c;
		}
		if ( sign(result_c) == sign(result_a) )
		{
			a = c;
		}
		else
		{
			b = c;			
		}
	}
}

/**
 * @brief      Determines if the error between two values is acceptable
 *
 * @param[in]  x_1   The first value
 * @param[in]  x_2   The second value
 *
 * @return     Bool answering if the error is acceptable
 */
bool acceptable_error(double x_1, double x_2, double epsilon)
{
	if ( std::abs((x_2 - x_1) / x_1) < epsilon )
	{
		return true;
	}
	return false;
}

/**
 * @brief      Determines the sign of a value
 *
 * @param[in]  x_1   The value
 *
 * @return     -1 if the value is negative, 1 otherwise
 */
int sign(double x_1)
{
	if ( x_1 < 0 )
	{
		return -1;
	}
	return 1;
}

/**
 * @brief      Generate energies to test
 *
 * @param[in]  E_range   The energy range
 * @param[in]  stepSize  The step size
 *
 * @return     A vector containing the energies
 */
std::vector<double> generate_energies(double E_range, double stepSize)
{
	double numberOfSteps = E_range / stepSize + 1;
	double energy;
	std::vector<double> energies;

	for (int i = 0; i < numberOfSteps; i++)
	{
		energies.push_back( i * stepSize - E_range ); // -E_range to zero
	}

	return energies;
}

/**
 * @brief      Generate a collection of the differences in the gradients for each energy
 *
 * @param[in]  evenSolutions  Denotes looking for the gradient differences of even solns
 * @param[in]  energies       The energies
 *
 * @return     A vector with the collection of gradient differences 
 */
std::vector<double> generate_f_E(bool evenSolutions, std::vector<double> energies, std::string type)
{
	Numerov_Method potential;
	potential.set_potential_type(type);
	std::vector<double> f_E;

	for (int i = 0; i < energies.size(); i++)
	{
		f_E.push_back( potential.f_E_(-energies[i], evenSolutions) );
	}

	return f_E;
}