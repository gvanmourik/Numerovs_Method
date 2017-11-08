#include <Numerov_Method.h>

// Declarations
std::vector<double> generate_energies(double E_range, double stepSize);
std::vector<double> generate_f_E(bool evenSolutions, std::vector<double> energies);
void export_vector_to_file(std::vector<double> x, std::vector<double> y, std::string fname);
double bisection_method(double a, double b, double epsilon, bool even);
bool acceptable_error(double x_1, double x_2, double epsilon);
int sign(double x_1);


int main( int argc, char* argv[] )
{
	// Variables
	std::vector<double> energies, f_E_even, f_E_odd;
	Numerov_Method wedge_potential;
	double root1, root2, root3;
	double epsilon = 1.0e-5;
	bool even = true;
	bool odd = false;

	// Generate f as a function of E
	/* UNCOMMENT BELOW SECTION TO FIND APPX ROOTS */
	energies = generate_energies(183.0, 0.1);
	f_E_even = generate_f_E(even, energies);
	f_E_odd = generate_f_E(odd, energies);
	wedge_potential.export_vector_to_file(energies, f_E_even, "f_E_even.data");
	wedge_potential.export_vector_to_file(energies, f_E_odd, "f_E_odd.data");

	// From observation of f(E) use the following energies
	// root1 = bisection_method(-75.4, -74.4, epsilon, even);
	// root2 = bisection_method(-52.8, -50.8, epsilon, odd);
	// root3 = bisection_method(-16.4, -15.8, epsilon, even);

	// printf("root1 = %f MeV\n", root1);
	// printf("root2 = %f MeV\n", root2);
	// printf("root3 = %f MeV\n", root3);

	// // Show progress towards convergence
	// finiteSqWell.plot_wavefunction(-30.0, even, "30");
	// finiteSqWell.plot_wavefunction(-20.0, even, "20");
	// finiteSqWell.plot_wavefunction(-16.0, even, "16");

	// // Show converged solutions
	// finiteSqWell.plot_wavefunction(root1, even, "root1");
	// finiteSqWell.plot_wavefunction(root2, odd, "root2");
	// finiteSqWell.plot_wavefunction(root3, even, "root3");


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
double bisection_method(double a, double b, double epsilon, bool even)
{
	Numerov_Method finiteSqWell;
	double c, result_a, result_c;

	while ( 1 )
	{
		c = (a+b)/2;
		result_c = finiteSqWell.f_E_(c, even);
		result_a = finiteSqWell.f_E_(a, even);

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
std::vector<double> generate_f_E(bool evenSolutions, std::vector<double> energies)
{
	Numerov_Method finiteSqWell;
	std::vector<double> f_E;

	for (int i = 0; i < energies.size(); i++)
	{
		f_E.push_back( finiteSqWell.f_E_(energies[i], evenSolutions) );
	}

	return f_E;
}