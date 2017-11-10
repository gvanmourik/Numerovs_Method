#include <Numerov_Method.h>

// Declarations
std::vector<double> generate_energies(double E_range, double stepSize);


int main( int argc, char* argv[] )
{
	// Variables
	std::vector<double> energies, f_E_even, f_E_odd;
	Numerov_Method parabolic_potential, wedge_potential;
	double root1_p, root2_p, root3_p, root4_p, root1_w, root2_w, root3_w, root4_w;
	bool even = true;
	bool odd = false;
	std::string p = "parabolic";
	std::string w = "wedge";

	// Set potential types
	parabolic_potential.set_potential_type(p);
	wedge_potential.set_potential_type(w);

	// Generate a range of energies in order to observe the roots
	energies = generate_energies(75.0, 0.1);

	// Parabolic energies
	f_E_even = parabolic_potential.generate_f_E_values(even, energies);
	f_E_odd = parabolic_potential.generate_f_E_values(odd, energies);
	parabolic_potential.export_xy_to_file(energies, f_E_even, "f_E_even_p.data");
	parabolic_potential.export_xy_to_file(energies, f_E_odd, "f_E_odd_p.data");
	// Wedge energies
	f_E_even = wedge_potential.generate_f_E_values(even, energies);
	f_E_odd = wedge_potential.generate_f_E_values(odd, energies);
	wedge_potential.export_xy_to_file(energies, f_E_even, "f_E_even_w.data");
	wedge_potential.export_xy_to_file(energies, f_E_odd, "f_E_odd_w.data");

	// Export k values
	// parabolic_potential.export_k_("k_parabolic.data");
	// wedge_potential.export_k_("k_wedge.data");

	// From observation of f(E) use the following energies
	root1_p = parabolic_potential.bisection_method(-23.0, -22.3, odd);
	root2_p = parabolic_potential.bisection_method(-62.7, -61.7, even);
	root1_w = wedge_potential.bisection_method(-10.5, -7.0, odd);
	root2_w = wedge_potential.bisection_method(-49.0, -50.0, even);

	// Print out root values for user
	printf("Parabolic root1 (odd) = %f MeV\n", root1_p);
	printf("Parabolic root2 (even) = %f MeV\n", root2_p);
	printf("Wedge root1 (odd) = %f MeV\n", root1_w);
	printf("Wedge root2 (even) = %f MeV\n", root2_w);

	// Show progress towards convergence for parabolic potential
	parabolic_potential.export_wavefunction(-35.0, odd, "35");
	parabolic_potential.export_wavefunction(-29.0, odd, "29");
	parabolic_potential.export_wavefunction(-23.0, odd, "23");

	// // Show converged solutions
	parabolic_potential.export_wavefunction(root1_p, odd, "root1_p");
	parabolic_potential.export_wavefunction(root2_p, even, "root2_p");
	wedge_potential.export_wavefunction(root1_w, odd, "root1_w");
	wedge_potential.export_wavefunction(root2_w, even, "root2_w");


	return 0;
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
	std::vector<double> energies;

	for (int i = 0; i < numberOfSteps; i++)
	{
		energies.push_back( i * stepSize - E_range ); // -E_range to zero
	}

	return energies;
}