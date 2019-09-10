
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// hierher kill headers?
// To generate (normally distributed) random numbers
//#include <iostream>
//#include <chrono>
//#include <random>

//========================================================================
/// Helper namespace 
//========================================================================
namespace PowerIterationHelperNamespace
{

 /* /// \short Construct a trivial random generator engine from a time-based */
 /* /// seed. Used to construct an initial guess vector (which isn't orthogonal */
 /* /// to the dominant eigenvector) */
 /* unsigned Seed=140; // hierher std::chrono::system_clock::now().time_since_epoch().count(); */
 
 /// \short The maximum number of iterations to use inside the power (and
 /// inverse power) iteration
 unsigned Max_iteration=100000;
 
 // The tolerance for convergence
 double Tolerance=1.0e-12;

 /// \short Boolean to indicate whether we want to compute the singular
 /// values (for the condition number) or the eigenvalues
 bool Compute_singular_values=true;
 
 /// Function that returns the sign of a double (returns zero if x=0)
 int sgn(const double& x)
 {
  // If x is less than zero then we get false-true=-1 and if x is
  // greater than zero we get true-false=1. The final case is if
  // x is in fact zero. In this case we get false-false=0.
  return (x>0.0)-(x<0.0);
 } // End of sgn
 

 /// \short Function that populates the entries of the DoubleVector 
 /// random_vector with (normally-distributed) random numbers
 void randn(DoubleVector& random_vector)
 {
  // Create a trivial random generator engine from a time-based seed
  //hierher std::default_random_engine generator(Seed);

  // Set the distribution to have zero mean and unit variance 
  // hierher std::normal_distribution<double> distribution(0.0,1.0);

  // How many rows?
  unsigned n_row=random_vector.nrow();
 
  // initialize random seed: 
  srand (time(NULL));

  // Loop over 
  for (unsigned i=0;i<n_row;i++)
   {
    // Generate a random number and store it
    random_vector[i]=double(rand())/double(RAND_MAX);
    //oomph_info << "random vector: " << random_vector[i] << std::endl;
  }

 } // End of randn
} // End of PowerIterationHelperNamespace



//========================================================================
/// \short Power method: used to compute the *largest* eigenvalue (and
/// corresponding eigenvector) of the matrix pointed to by matrix_pt
/// Based on implementation described here:
///             http://www.cs.huji.ac.il/~csip/tirgul2.pdf
//========================================================================
double power_method(CRDoubleMatrix* const matrix_pt, DoubleVector& eigenvector)
{
 // Grab the LinearAlgebraDistribution associated with the input matrix
 LinearAlgebraDistribution* dist_pt=matrix_pt->distribution_pt();
 
 // If the eigenvector hasn't been built by the user
 if (eigenvector.built()==0)
 {
  // Build it!
  eigenvector.build(dist_pt,false);
 } 

 // Boolean to indicate convergence
 bool has_converged=false;
 
 // Allocate space for the initial guess
 DoubleVector q_old(dist_pt,0.0);
 
 // Allocate space for the updated guess
 DoubleVector q_new(dist_pt,0.0);
 
 // Allocate space for the matrix-vector product A*q_old
 DoubleVector z_new(dist_pt,0.0);
 
 // Temporary vector for intermediate calculations (only needed if we're
 // computing singular values instead of the eigenvalues)
 DoubleVector sing_calc_vec(dist_pt,0.0);
 
 // Populate the entries of q_old
 PowerIterationHelperNamespace::randn(q_old);
 
 // Normalise it
 q_old/=q_old.norm();

 // Storage for the eigenvalue
 double lambda=0.0;
 
 // The difference between the previous and current guess (for the eigenvector)
 double q_diff_norm=0.0;

 // How many rows in the matrix?
 unsigned n_row=matrix_pt->nrow();
 
 // How many iterations did it take to converge?
 unsigned n_iteration=0;
 
 // Iterate until convergence
 for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)
 {  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute q_new=A*q_old 
   matrix_pt->multiply(q_old,q_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{T}*A*q_old so also multiply by A^{T} 
  else
  {
   // Compute sing_calc_vec=A*q_old 
   matrix_pt->multiply(q_old,sing_calc_vec);
   
   // Compute q_new=A^{T}*A*q_old=A^{T}*sing_calc_vec
   matrix_pt->multiply_transpose(sing_calc_vec,q_new);
  }

  // Normalise
  q_new/=q_new.norm();
  
  // If we're computing eigenvalues compute q_new=A*q_new
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute z_new=A*q_new 
   matrix_pt->multiply(q_new,z_new);
   
   // ...and now compute the eigenvalue
   lambda=q_new.dot(z_new);
  }
  // If we're computing singular values we actually need to compute
  // z_new=A^{T}*A*q_new so also multiply by A^{T} 
  else
  {
   // Compute sing_calc_vec=A*q_new 
   matrix_pt->multiply(q_new,sing_calc_vec);
   
   // Compute z_new=A^{T}*A*q_new=A^{T}*sing_calc_vec
   matrix_pt->multiply_transpose(sing_calc_vec,z_new);
   
   // ...and now compute the singular value
   lambda=std::sqrt(q_new.dot(z_new));
  }

  // Increment the iteration counter
  n_iteration++;

  // Scale q_old and q_new so they have the same sign in the first entry as
  // the eigenvectors can be invariant up to a sign change. Scale q_old first
  q_old*=PowerIterationHelperNamespace::sgn(q_old[0]);
  
  // Now scale q_new in the same way
  q_new*=PowerIterationHelperNamespace::sgn(q_new[0]);
  
  // The norm of the difference between the previous and current guess
  q_diff_norm=0.0;

  // Loop over the entries of the vector
  for (unsigned j=0;j<n_row;j++)
  {
   // Update the value of q_diff_norm
   q_diff_norm+=std::pow(q_new[j]-q_old[j],2);
  }

  // Now square root it
  q_diff_norm=std::sqrt(q_diff_norm);
  
  // Check if the convergence tolerance has been met
  if (q_diff_norm<PowerIterationHelperNamespace::Tolerance)
  {
   // We've converged
   has_converged=true;
   
   // We're done so jump out now
   break;   
  }

  // If we're not done yet, copy the updated vector into the old one
  q_old=q_new;
 } // for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)

 // If we've actually converged
 if (has_converged)
 {  
  // Document the convergence
  oomph_info << "\nPower method converged within " << n_iteration
	     << " iterations to accuracy: " << q_diff_norm
	     << std::endl;
 }
 // Otherwise we've run out of iterations
 else
 {
  // Document the (non-)convergence
  oomph_info << "\nPower method has not converged within " << n_iteration
	     << " iterations (with\ntolerance: "
	     << PowerIterationHelperNamespace::Tolerance
	     << "). Difference between previous and current\nguess: "
	     << q_diff_norm << std::endl;
  abort();
 }

 // Copy the approximation to the dominant eigenvector over
 eigenvector=q_new;
 
 // Return the largest eigenvalue
 return lambda; 
} // End of power_method


//========================================================================
/// \short Inverse power method: used to compute the *smallest* eigenvalue
/// (and corresponding eigenvector) of the matrix pointed to by matrix_pt
/// Based on implementation described here:
///             http://www.cs.huji.ac.il/~csip/tirgul2.pdf
//========================================================================
double inverse_power_method(CRDoubleMatrix* const matrix_pt,
			    DoubleVector& eigenvector)
{
 // Grab the LinearAlgebraDistribution associated with the input matrix
 LinearAlgebraDistribution* dist_pt=matrix_pt->distribution_pt();
 
 // If the eigenvector hasn't been built by the user
 if (eigenvector.built()==0)
 {
  // Build it!
  eigenvector.build(dist_pt,false);
 } 

 // Handle to the LinearSolver associated with the input matrix
 LinearSolver* solver_pt=matrix_pt->linear_solver_pt();
   
 // Enable resolves for later
 solver_pt->enable_resolve();

 // ...and be quiet!
 solver_pt->disable_doc_time();
 
 // Allocate space for the transpose of the matrix (only used if we're
 // computing singular values instead of the eigenvalues)
 CRDoubleMatrix* matrix_transpose_pt=new CRDoubleMatrix;

 // Handle to the LinearSolver associated with the transpose of the input matrix
 LinearSolver* transpose_solver_pt=0;
 
 // If we're computing singular values we need the matrix transpose
 if (PowerIterationHelperNamespace::Compute_singular_values)
 {
  // Compute the transpose of the input matrix and store it
  matrix_pt->get_matrix_transpose(matrix_transpose_pt);
  
  // Get the LinearSolver associated with the transpose of the input matrix
  transpose_solver_pt=matrix_transpose_pt->linear_solver_pt();
 
  // Enable resolves for later (using the transposed matrix solver)
  transpose_solver_pt->enable_resolve();

  // ...and, again, be quiet!
  transpose_solver_pt->disable_doc_time();
 }

 // Make sure everything has been set up to do resolves
 {
  // Dummy RHS vector 
  DoubleVector rhs_temp(dist_pt,0.0);
 
  // Dummy LHS vector
  DoubleVector lhs_temp(dist_pt,0.0);
  
  // Solve the system once to set everything necessary for a resolve up
  matrix_pt->solve(rhs_temp,lhs_temp);
  
  // If we're computing singular values we need the matrix transpose
  if (PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Solve the system once to set everything necessary for a resolve up
   matrix_transpose_pt->solve(rhs_temp,lhs_temp);
  }
 } // Set up everything necessary for a resolve
 
 // Boolean to indicate convergence
 bool has_converged=false;
 
 // Allocate space for the initial guess
 DoubleVector q_old(dist_pt,0.0);
 
 // Allocate space for the updated guess
 DoubleVector q_new(dist_pt,0.0);
 
 // Allocate space for the matrix-vector product A^{-1}*q_old
 DoubleVector z_new(dist_pt,0.0);
 
 // Temporary vector for intermediate calculations (only needed if we're
 // computing singular values instead of the eigenvalues)
 DoubleVector sing_calc_vec(dist_pt,0.0);
 
 // Populate the entries of q_old
 PowerIterationHelperNamespace::randn(q_old);
 
 // Normalise it
 q_old/=q_old.norm();

 // Storage for the eigenvalue
 double lambda=0.0;
 
 // The difference between the previous and current guess (for the eigenvector)
 double q_diff_norm=0.0;

 // How many rows in the matrix?
 unsigned n_row=matrix_pt->nrow();

 // How many iterations did it take to converge?
 unsigned n_iteration=0;
 
 // Iterate until convergence
 for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)
 {  
  // Reinitialise
  q_new.initialise(0.0);
  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute q_new=A^{-1}*q_old 
   solver_pt->resolve(q_old,q_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{T}*A*q_old so also multiply by A^{T} 
  else
  {
   // Compute q_new=A^{-1}*q_old 
   solver_pt->resolve(q_old,sing_calc_vec);
   
   // Compute q_new=A^{-T}*A^{-1}*q_old=A^{-T}*sing_calc_vec
   transpose_solver_pt->resolve(sing_calc_vec,q_new);
  }
  
  // Normalise
  q_new/=q_new.norm();
  
  // Reinitialise
  z_new.initialise(0.0);
  
  // If we're computing eigenvalues compute q_new=A*q_old
  if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Compute z_new=A^{-1}*q_new
   solver_pt->resolve(q_new,z_new);
   
   // ...and now compute the singular value (recalling that we're actually
   // calculating 1/lambda with the inverse power method)
   lambda=1.0/q_new.dot(z_new);
  }
  // If we're computing singular values we actually need to compute
  // q_new=A^{-T}*A^{-1}*q_old so also multiply by A^{-T} 
  else
  {
   // Compute z_new=A^{-1}*q_new
   solver_pt->resolve(q_new,sing_calc_vec);
   
   // Compute z_new=A^{-T}*A^{-1}*q_old=A^{-T}*sing_calc_vec
   transpose_solver_pt->resolve(sing_calc_vec,z_new);
  
   // ...and now compute the singular value (recalling that we're actually
   // calculating 1/lambda with the inverse power method)
   lambda=std::sqrt(1.0/q_new.dot(z_new));
  }
  
  // Increment the iteration counter
  n_iteration++;

  // Scale q_old and q_new so they have the same sign in the first entry as
  // the eigenvectors can be invariant up to a sign change. Scale q_old first
  q_old*=PowerIterationHelperNamespace::sgn(q_old[0]);
  
  // Now scale q_new in the same way
  q_new*=PowerIterationHelperNamespace::sgn(q_new[0]);
  
  // The norm of the difference between the previous and current guess
  q_diff_norm=0.0;

  // Loop over the entries of the vector
  for (unsigned j=0;j<n_row;j++)
  {
   // Update the value of q_diff_norm
   q_diff_norm+=std::pow(q_new[j]-q_old[j],2);
  }

  // Now square root it
  q_diff_norm=std::sqrt(q_diff_norm);

  // Check if the convergence tolerance has been met
  if (q_diff_norm<PowerIterationHelperNamespace::Tolerance)
  {
   // We've converged
   has_converged=true;
   
   // We're done so jump out now
   break;   
  }

  // If we're not done yet, copy the updated vector into the old one
  q_old=q_new;
 } // for (unsigned i=0;i<PowerIterationHelperNamespace::Max_iteration;i++)

 // If we've actually converged
 if (has_converged)
 {
  // Document the convergence
  oomph_info << "\nInverse power method converged within " << n_iteration
	     << " iterations to accuracy: " << q_diff_norm
	     << std::endl;
 }
 // Otherwise we've run out of iterations
 else
 {
  // Document the (non-)convergence
  oomph_info << "\nInverse power method has not converged within "
	     << n_iteration << " iterations (with\ntolerance: "
	     << PowerIterationHelperNamespace::Tolerance
	     << "). Difference between previous and current\nguess: "
	     << q_diff_norm << std::endl;
  abort();
 }

 // Copy the approximation to the dominant eigenvector over
 eigenvector=q_new;
 
 // Return the smallest eigenvalue
 return lambda; 
} // End of inverse_power_method

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
