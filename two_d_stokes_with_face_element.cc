//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//Driver for Navier-Stokes in backward step domain -- meshed with triangle

//Generic includes
#include "generic.h"

// QUEHACERES delete
//#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

// QUEHACERES delete once we've got the stokes element working
// #include "new_poisson_sing_face_element.h"

#include "navier_stokes_sing_face_element.h"

#include "power_method.h"

// wrapper for locate zeta to allow us to visualise along a given line
#include "generic/line_visualiser.h"

using namespace std;
using namespace oomph;

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
  // directory for output 
  string dir = "RESLT";
  
  // dimensionality of the problem
  const unsigned Dim = 2;
  
  // Constant term in singular fct
  double Constant_on_singular_boundary = 0.0;

  /// Number that identifies which
  /// case we're doing when checking the condition number
  /// 0: real problem; 1: FE-only; 2: fake r_c equation
  unsigned Problem_type_for_check_condition_number = 0;

  // Default uniform element area
  double Uniform_element_area = 0.1;

  // element area in region around corner
  double High_res_element_area = 1e-3;
  
  // If true: Traction bc on outflow boundary
  bool Do_traction_problem = true;

  // Default rescaling factor for domain
  double Scaling_factor_for_domain = 1.0;

  // fake singular amplitude for debug
  double singular_amplitude_for_debug = 1.0;
  
  // Dimensionless domain values
  // ---------------------------

  /// Dimless width of inflow channel
  double H_up = 1.0; 

  /// Dimless length of channel upstream of step
  double L_up = 2.0; 

  /// \short Dimless length of channel downstream of step
  double L_down = 2.0; 

  /// Dimless width of outflow channel
  double H_down = 2.0; 

  /// Radius of internal boundary (surrounding the singularity)
  double Radius_of_internal_boundary = 0.2; // 0.5;

  // half the convex angle of the step
  double alpha = 3.0 * MathematicalConstants::Pi / 4.0;
  
  // IDs for the boundaries
  enum
  {
    Inflow_boundary_id                   = 0,
    Top_boundary_id                      = 1,
    Outflow_boundary_id                  = 2,
    Bottom_boundary_id                   = 3,
    Vertical_step_boundary_id            = 4,
    Corner_boundary_id                   = 5,
    Horizontal_step_boundary_id          = 6,
    Enriched_region_lower_boundary_id    = 7,
    Enriched_region_upper_boundary_id    = 8,
    Enriched_region_dividing_boundary_id = 9,
  };

  enum
  {
    Enriched_region_upper_id = 1,
    Enriched_region_lower_id = 2
  };
  
  // hierher kill
  // /// Test for assigning C directly
  // double Imposed_amplitude = 0.0;

  // Bit of a hack but it facilitates reuse...
#include "unstructured_backward_step_mesh.h"

  // Kronecker delta function
  int delta(unsigned i, unsigned j)
  {
    return (i == j);
  }

  // sign function
  template <typename T>
  int sgn(T val)
  {
    return (T(0) < val) - (val < T(0));
  }
  
  /// \short Newtonian stress tensor
  DenseMatrix<double> get_stress(const DenseMatrix<double>& strain_rate,
				 const double& p)
  {
    // \tau_{ij}
    DenseMatrix<double> stress(2,2);

    for(unsigned i=0; i<2; i++)
    {
      for(unsigned j=0; j<2; j++)
      {
	// Newtonian constitutive relation
	stress(i,j) = -p*delta(i,j) + 2.0*strain_rate(i,j);
      }
    }

    return stress;
  }
  
  /// Non-singular part of the solution for testing
  void u_non_singular(const Vector<double>& x, Vector<double>& u);

  /// \short "Singular" function and gradient
  void singular_fct_and_gradient(const Vector<double>& x,
				 Vector<double>& u, DenseMatrix<double>& du_dx)
  {    
    // big number to numerically represent infinity; value of 1/r at 1/10th of
    // the element lengthscale
    const double infinity = 10.0 / sqrt(2.0 * Uniform_element_area);

    // Radius & polar angle, accounting for the fact the singularity is not at the
    // origin but at (L_up, 0)
    double r = sqrt( (x[0]-L_up) * (x[0]-L_up) + x[1]*x[1] );

    bool at_origin = (r == 0);

    // QUEHACERES delete
    // // catch the singularity point
    // if (at_origin)
    // {
    //   // shift the radius to 1/10th of the element lengthscale
    //   r = 0.1*sqrt(2.0 * Uniform_element_area);      
    // }   
    
    // make sure we've got enough storage
    u.resize(3);
            
    // Little hack to make sure atan2 doesn't overshoot
    double y = x[1];
    double tol_y = -1.0e-12;
    if ((y<0.0) && (y>tol_y)) y = 0.0;

    // angle w.r.t. x-axis from the singularity at the corner, coords (L_up, 0)
    double theta = atan2(y, (x[0]-L_up));

    // angle w.r.t. to axis of symmetry between the two walls of the step,
    // the angle the analytic solution is given in terms of
    double phi = theta + alpha - MathematicalConstants::Pi;
    
    // streamfunction exponent, from Mathematica
    double lambda = 1.544483736782746;

    // polar components in r and phi directions
    double ur = 0, v = 0;

    // ------------------------------------------------------------------
    // components in polar coordinates
    // ------------------------------------------------------------------
    ur = pow(r, lambda-1.0)*((lambda-2.0)*cos(alpha*lambda)*sin((lambda-2.0)*phi)
			     -lambda*cos(alpha*(lambda-2.0))*sin(lambda*phi));
  
    v  = pow(r,lambda-1.0)*(-lambda*cos(alpha*(lambda-2.0))*cos(lambda*phi) +
			    lambda*cos((lambda-2.0)*phi)*cos(alpha*lambda) );

    // ------------------------------------------------------------------
    // derivatives of polar components w.r.t. polar components
    // ------------------------------------------------------------------
    double dudr =  pow(r,-2 + lambda)*(-1 + lambda) *
      ((-2 + lambda)*cos(alpha*lambda)*sin(phi*(-2 + lambda)) - 
       lambda*cos(alpha*(-2 + lambda))*sin(phi*lambda));

    double dudphi = pow(r,-1 + lambda) *
      (-(pow(lambda,2) * cos(alpha*(-2 + lambda))*cos(phi*lambda)) + 
       pow(-2 + lambda,2)*cos(phi*(-2 + lambda))*cos(alpha*lambda));

    double dvdr = pow(r,-2 + lambda) *
      (-1 + lambda)*(-(lambda*cos(alpha*(-2 + lambda))*cos(phi*lambda)) + 
		     lambda*cos(phi*(-2 + lambda))*cos(alpha*lambda));

    double dvdphi = pow(r,-1 + lambda)*lambda *
      (-((-2 + lambda)*cos(alpha*lambda)*sin(phi*(-2 + lambda))) + 
       lambda*cos(alpha*(-2 + lambda))*sin(phi*lambda));
    
    // Cartesian components (in polar coords, i.e. ux(r,theta) and uy(r,theta) )
    // note we're now using the angle w.r.t. the x-axis
    u[0] = ur * cos(theta) - v * sin(theta);
    u[1] = ur * sin(theta) + v * cos(theta);

    // singular pressure \hat p
    u[2] = 4.0*pow(r,lambda-2.0)*(lambda-1.0)*cos(alpha*lambda)*sin((lambda-2.0)*phi);
    
    // ------------------------------------------------------------------
    // derivatives of Cartesian components w.r.t. Cartesian coordinates
    // ------------------------------------------------------------------

    // make sure we've got a 2x2
    du_dx.resize(2);
        
    // QUEHACERES not sure these are right, probably delete
    // // S is tensor \bm{\grad u} in polar coordinates
    // DenseMatrix<double> S(2,2);
    // S(0,0) = dudr;               // (\bm{\grad u})_{rr}
    // S(0,1) = (1.0/r)*dudphi;     // (\bm{\grad u})_{r\theta}
    // S(1,0) = dvdr;               // (\bm{\grad u})_{\theta r}
    // S(1,1) = (1.0/r)*dvdphi;     // (\bm{\grad u})_{\theta\theta}

    // // dux_dx
    // du_dx(0,0) = cos(theta) * ( S(0,0)*cos(theta) - S(0,1)*sin(theta) ) -
    //   sin(theta) * ( S(1,0)*cos(theta) - S(1,1)*sin(theta) );

    // // dux_dy
    // du_dx(0,1) = cos(theta) * ( S(0,0)*sin(theta) + S(0,1)*cos(theta) ) -
    //   sin(theta) * ( S(1,0)*sin(theta) + S(1,1)*cos(theta) );

    // // duy_dx
    // du_dx(1,0) = sin(theta) * ( S(0,0)*cos(theta) - S(0,1)*sin(theta) ) +
    //   cos(theta) * ( S(1,0)*cos(theta) - S(1,1)*sin(theta) );

    // // duy_dy
    // du_dx(1,1) = sin(theta) * ( S(0,0)*sin(theta) + S(0,1)*cos(theta) ) +
    //   cos(theta) * ( S(1,0)*sin(theta) + S(1,1)*cos(theta) );
    
    // dux_dx
    du_dx(0,0) = cos(theta) * (dudr*cos(theta) - dvdr*sin(theta))
      -(1.0/r)*sin(theta) * (dudphi*cos(theta) - ur*sin(theta) -
    			     dvdphi*sin(theta) - v*cos(theta));
  
    // dux_dy
    du_dx(0,1) = sin(theta) * (dudr*cos(theta) - dvdr*sin(theta))
      +(1.0/r)*cos(theta) * (dudphi*cos(theta) - ur*sin(theta) -
    			     dvdphi*sin(theta) - v*cos(theta));

    // duy_dx
    du_dx(1,0) = cos(theta) * (dudr*sin(theta) + dvdr*cos(theta))
      -(1.0/r)*sin(theta) * (dudphi*sin(theta) + ur*cos(theta) +
    			     dvdphi*cos(theta) - v*sin(theta));

    // duy_dy
    du_dx(1,1) = sin(theta) * (dudr*sin(theta) + dvdr*cos(theta))
      +(1.0/r)*cos(theta) * (dudphi*sin(theta) + ur*cos(theta) +
    			     dvdphi*cos(theta) - v*sin(theta));
    
    // catch the point exactly at the origin
    if(at_origin)
    {
      // from BCs
      u[0] = 0.0;
      u[1] = 0.0; 
      u[2] = sgn(sin((2.0 - lambda)*phi)) * infinity;

      du_dx(0,0) =  infinity;
      du_dx(0,1) =  infinity;
      du_dx(1,0) = -infinity;
      du_dx(1,1) = -infinity;
    }
    
    // QUEHACERES what dis?
    // Now add constant
    // u+=Constant_on_singular_boundary;
  }

  /// \short "Singular" function
  Vector<double> singular_fct(const Vector<double>& x)
  {
    Vector<double> u(3);
    DenseMatrix<double> du_dx(2,2);
    
    singular_fct_and_gradient(x, u, du_dx);
    
    return u;
  }

  /// \short Gradient of "Singular" function
  DenseMatrix<double> gradient_of_singular_fct(const Vector<double>& x)
  {
    Vector<double> u(3, 0.0);
    DenseMatrix<double> du_dx(2,2);    
    
    singular_fct_and_gradient(x, u, du_dx);

    return du_dx;
  }

  /// Exact solution
  void u_exact(const Vector<double>& x, Vector<double>& u)
  {
    // make sure we have the right amount of storage
    u.resize(3);
    
    if (!CommandLineArgs::command_line_flag_has_been_set
	("--suppress_sing_in_exact_soln"))
    {
      Vector<double> u_sing(3);
      u_sing = singular_fct(x);
      
      for(unsigned i=0; i<3; i++)
      {
	u[i] += u_sing[i];
      }
    }

    // Add non-singular part of the solution
    Vector<double> u_reg(3);
    u_non_singular(x, u_reg);
      
    for(unsigned i=0; i<2; i++)
    {
      u[i] += u_reg[i];
    }
  }

  /// Non-singular part of the solution for testing
  void u_and_gradient_non_singular(const Vector<double>& x,
				   Vector<double>& u,
				   DenseMatrix<double>& dudx)
  {
    // QUEHACERES needs updating
    
  }


  /// Non-singular part of the solution for testing
  void u_non_singular(const Vector<double>& x, Vector<double>& u)
  {
    u.resize(3);
    DenseMatrix<double> dudx(2,2);
    
    u_and_gradient_non_singular(x, u, dudx);
  }

  void prescribed_traction(const Vector<double>& x, Vector<double>& traction)
  {
    double tol=1.0e-6;
    Vector<double> normal(2);

    // rightmost boundary
    if(fabs(x[0]-2.0-2.0*Scaling_factor_for_domain) < Scaling_factor_for_domain*tol)
    {
      normal[0] = 1.0;
      normal[1] = 0.0;
    }

    // Vertical boundary below the step
    if(fabs(x[0]-2.0)<Scaling_factor_for_domain*tol)
    {
      normal[0] = -1.0;
      normal[1] =  0.0;
    }


    // 
    if(fabs(x[0]-2.0*(1.0-Scaling_factor_for_domain)) < Scaling_factor_for_domain*tol)
    {
      normal[0] = -1.0;
      normal[1] =  0.0;
    }

    if(fabs(x[1]-1.0*Scaling_factor_for_domain) < Scaling_factor_for_domain*tol)
    {
      normal[0] = 0.0;
      normal[1] = 1.0;
    }

    if(fabs(x[1]) < Scaling_factor_for_domain*tol)
    {
      normal[0] =  0.0;
      normal[1] = -1.0;
    }

    if(fabs(x[1]+1.0*Scaling_factor_for_domain) < Scaling_factor_for_domain*tol)
    {
      normal[0] =  0.0;
      normal[1] = -1.0;
    }

    // make sure we've got enough storage
    traction.resize(Dim);
    
    // enforce parallel outflow with parabolic profile
    traction[0] = 0;
    traction[1] = x[1] * normal[0];
    
    // regular and singular parts of the solution        
    Vector<double> u_reg;
    u_non_singular(x, u_reg);
    
    Vector<double> u_sing = singular_fct(x);

    // QUEHACERES come back to the reg/sing shizzle above
  }

  /// Function to specify boundary conditions
  Vector<double> u_BC(const Vector<double>& x, const unsigned& boundary_id)
  {
    // no-slip as default
    Vector<double> u(Dim, 0.0);

    if(boundary_id == Inflow_boundary_id)
    {
      // apply parabolic velocity profile to the x-component,
      // i.e. standard flow profile in a pipe
      u[0] = x[1]*(H_up - x[1]);      
    }
    
    return u;
  }

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Navier-Stokes in backward-facing step. Dirichlet boundary conditions along
/// all boundaries.
//====================================================================
template<class ELEMENT>
class StepProblem : public Problem
{

public:


  /// Constructor
  StepProblem();

  /// Destructor 
  ~StepProblem()
    {
      // hierher: at some point delete things properly and do memory leak
      // check
      delete Bulk_mesh_pt->spatial_error_estimator_pt();
      delete Bulk_mesh_pt;
    }
 
  /// Update the after solve (empty)
  void actions_after_newton_solve(){}
 
  /// \short Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}
 
  // Perform actions after mesh adaptation
  void actions_after_adapt()
    {
      // Recreate face elements
      create_face_elements();
   
      // Complete problem setup
      complete_problem_setup();

      // Rebuild global mesh
      rebuild_global_mesh();
    }
 
  /// Perform actions after mesh adaptation (empty)
  void actions_before_adapt()
    {
      // Kill face elements
      delete_face_elements();

      // Rebuild global mesh
      rebuild_global_mesh();
    }
 
  /// Access function for the specific mesh
  RefineableTriangleMesh<ELEMENT>* mesh_pt()
    {
      return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
    }
 
  /// Doc the solution
  void doc_solution();

  /// Validation
  void impose_amplitude_runs();

  /// hierher more validation
  void check_residual_for_exact_non_singular_fe_soln();

  /// Check condition number
  void check_condition_number();

  /// Assign nodal values to be the exact singular solution
  void set_values_to_singular_solution();

  // function which validates the singular stress 
  void validate_stress();
  
  DocInfo* doc_info()
    {
      return &Doc_info;
    }
  
private:

  /// Do what it says
  void complete_problem_setup();
 
  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();


  /// Helper function to duplicate the corner node, to allow for a
  //dipole-like pressure singularity
  void duplicate_corner_node();
  
  /// hierher Delete face elements and flush meshes
  void delete_face_elements()
  {
    // Loop over the flux elements
    unsigned n_element = Traction_boundary_condition_mesh_pt->nelement();
    for(unsigned e=0;e<n_element;e++)
    {
      // Kill
      delete Traction_boundary_condition_mesh_pt->element_pt(e);
    }
   
    // Wipe the mesh
    Traction_boundary_condition_mesh_pt->flush_element_and_node_storage();

    if (CommandLineArgs::command_line_flag_has_been_set
	("--dont_subtract_singularity"))
    {
      return;
    }

    // hierher
    // // Loop over the flux jump elements
    // unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
    // for(unsigned e=0;e<n_element;e++)
    //  {
    //   delete Face_mesh_for_flux_jump_pt->element_pt(e);
    //  }
   
    // // hierher: actually kill nodes too because they've been duplicated

    // // Wipe the mesh
    // Face_mesh_for_flux_jump_pt->flush_element_and_node_storage();

    // Loop over the bc elements
    n_element = Face_mesh_for_bc_pt->nelement();
    for(unsigned e=0;e<n_element;e++)
    {
      // Kill
      delete Face_mesh_for_bc_pt->element_pt(e);
    }
   
    // Wipe the mesh
    Face_mesh_for_bc_pt->flush_element_and_node_storage();

    // Loop over the integral face elements
    n_element = Face_mesh_for_singularity_integral_pt->nelement();
    for(unsigned e=0;e<n_element;e++)
    {
      delete Face_mesh_for_singularity_integral_pt->element_pt(e);
    }
    Face_mesh_for_singularity_integral_pt->flush_element_and_node_storage();
  }

  /// Create face elements
  void create_face_elements()
  {   
    // Flux boundary conditions?
    if (Global_Physical_Variables::Do_traction_problem)
    {
      // Only doing the outflow boundary has traction conditions
      unsigned i_bound = Global_Physical_Variables::Outflow_boundary_id;
	
      // Traction boundaries
      // for(unsigned i_bound = 0; i_bound<4; i_bound++)
      {
	unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
	for(unsigned e=0; e<n_element; e++)
	{
	  //Create Pointer to bulk element adjacent to the boundary
	  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
	    (Bulk_mesh_pt->boundary_element_pt(i_bound, e));
         
	  //Get Face index of boundary in the bulk element
	  int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound,e);
         
	  //Create corresponding face element
	  NavierStokesWithSingularityTractionElement<ELEMENT>* traction_element_pt =
	    new NavierStokesWithSingularityTractionElement<ELEMENT>(
	      bulk_elem_pt, face_index);
         
	  // Set the pointer to the prescribed traction function
	  traction_element_pt->traction_fct_pt() = 
	    &Global_Physical_Variables::prescribed_traction;
         
	  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
	  {
	    // We pass the pointer of singular function element to the 
	    // face element (Set function because it also declares 
	    // the amplitude to be external data for that element).
	    traction_element_pt->set_navier_stokes_sing_el_pt(
	      dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
		Singular_fct_element_mesh_pt->element_pt(0)));
	  }
	    
	  //Attach it to the mesh
	  Traction_boundary_condition_mesh_pt->add_element_pt(traction_element_pt);
	}
      }
    }
    // Dirichlet boundary conditions
    else
    {
      if (CommandLineArgs::command_line_flag_has_been_set
	  ("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
      {
	// BC elements live on the whole outer boundary
	for(unsigned i_bound=0;
	    i_bound<Global_Physical_Variables::Enriched_region_lower_boundary_id; i_bound++)
	{	    
	  unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
          
	  // Loop over the bulk elements adjacent to boundary b
	  for(unsigned e=0; e<n_element; e++)
	  {
	    // Get pointer to the bulk element that is adjacent to boundary b
	    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
	      Bulk_mesh_pt->boundary_element_pt(i_bound, e));
	      
	    //Find the index of the face of element e along boundary b 
	    int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound, e);
            
	    // Build the corresponding bc element
	    NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_element_pt =
	      new NavierStokesWithSingularityBCFaceElement<ELEMENT>
	      (bulk_elem_pt, face_index, BC_el_id);
            
	    // Tell the element about the singular fct
	    if (!CommandLineArgs::command_line_flag_has_been_set
		("--dont_subtract_singularity"))
	    {
	      bc_element_pt->set_navier_stokes_sing_el_pt(
		dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
		  Singular_fct_element_mesh_pt->element_pt(0)));
	    }
            
	    //Add the bc element to the surface mesh
	    Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
	  }
	}
      }
    }

    if (CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
    {
      return;
    }
      
    // QUEHACERES need to update if we re-enable
    // hierher renable


    // // Map keeps a running count of duplicate nodes already created;
    // // existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
    // std::map<Node*,Node*> existing_duplicate_node_pt;


    // Flux jump elements on boundary 7 of region 1:
    //----------------------------------------------
    // NOTE: Since these duplicate nodes, these elements must be
    //----------------------------------------------------------
    //       constructed first!
    //       ------------------
    // {
    //  // Where are we? Inside region 1 on boundary 7
    //  unsigned b=7;
    //  unsigned region_id=1;
    //  unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
    //  for (unsigned e=0;e<nel;e++)
    //   {
    //    FiniteElement* el_pt=
    //     Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
      
    //    // What is the index of the face of the bulk element at the boundary
    //    int face_index = Bulk_mesh_pt->
    //     face_index_at_boundary_in_region(b,region_id,e);
      
    //    // Build the corresponding flux jump element
    //    PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* flux_jump_element_pt 
    //     = new PoissonWithSingularityFluxJumpFaceElement<ELEMENT>
    //     (el_pt,face_index,existing_duplicate_node_pt,Flux_jump_el_id);

    //    //Add the flux jump element to the mesh
    //    Face_mesh_for_flux_jump_pt->add_element_pt(flux_jump_element_pt);
    //   }
    // }
   
    // Now add all new (duplicated) nodes to mesh
    // for (std::map<Node*,Node*>::iterator it=
    //       existing_duplicate_node_pt.begin();
    //      it!=existing_duplicate_node_pt.end();it++)
    //  {
    //   Face_mesh_for_flux_jump_pt->add_node_pt((*it).second);
    //  }


    // Create the face elements needed to compute the amplitude of
    // the singular function,
   
    // Only outer boundaries
    for(unsigned i_bound=0;
	i_bound<Global_Physical_Variables::Enriched_region_lower_boundary_id; i_bound++)
    {
      unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
      for(unsigned e=0; e<n_element; e++)
      {
	//Create Pointer to bulk element adjacent to the boundary
	ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
	  (Bulk_mesh_pt->boundary_element_pt(i_bound, e));
       
	//Get Face index of boundary in the bulk element
	int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound, e);
       
	//Create corresponding face element
	NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>* 
	  boundary_integral_face_element_pt =
	  new NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>(
	    bulk_elem_pt, face_index);

	//We pass the pointer of singular function element to the face element
	boundary_integral_face_element_pt->navier_stokes_sing_el_pt() =
	  dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	    Singular_fct_element_mesh_pt->element_pt(0));
       
	//Attach it to the mesh
	Face_mesh_for_singularity_integral_pt->add_element_pt(boundary_integral_face_element_pt);
      }
    }
   
    // Update the pointer to the face elements (currently needed so
    // this GeneralisedElement can assemble the contributions to the
    // r_C residual from the face elements!
    dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
      Singular_fct_element_mesh_pt->element_pt(0))->
      set_mesh_of_face_elements(Face_mesh_for_singularity_integral_pt);


    // // Now loop over bulk elements 
    // //----------------------------
    // // and swap over any of their nodes have been replaced
    // //----------------------------------------------------
    // unsigned n_el=Bulk_mesh_pt->nelement();
    // for (unsigned e=0;e<n_el;e++)
    //  {
    //   ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
    //    Bulk_mesh_pt->element_pt(e));
   
    //   // Loop over all nodes and check if they're amongst the replaced
    //   // ones
    //   unsigned nnod=bulk_el_pt->nnode();
    //   for (unsigned j=0;j<nnod;j++)
    //    {
    //     Node* nod_pt=bulk_el_pt->node_pt(j);
     
    //     // Find this original node in the map; if we find it
    //     // it's already been duplicated
    //     std::map<Node*,Node*>::iterator it=
    //      existing_duplicate_node_pt.find(nod_pt);
    //     if (it!=existing_duplicate_node_pt.end())
    //      {
    //       // Use the existing duplicate node
    //       bulk_el_pt->node_pt(j)=(*it).second;
    //      }
    //    }   
    //  }
  }
 
  /// Pointer to the bulk mesh
  RefineableTriangleMesh<ELEMENT> *Bulk_mesh_pt;
 
  /// Face element mesh for jump in interior of domain
  // hierher Mesh* Face_mesh_for_flux_jump_pt;

  /// Face element mesh for BC (Lagrange multiplier!) 
  Mesh* Face_mesh_for_bc_pt;

  /// \short Face elements used to compute the amplitude of the singular
  /// function
  Mesh* Face_mesh_for_singularity_integral_pt;
 
  /// Mesh for (single) element containing singular fct
  Mesh* Singular_fct_element_mesh_pt;
 
  /// Mesh of face elements applying traction bc 
  Mesh* Traction_boundary_condition_mesh_pt;

  /// \short Enumeration for IDs of FaceElements (used to figure out
  /// who's added what additional nodal data...)
  enum
  {
    Flux_jump_el_id,
    BC_el_id
  };

  // document some tingz
  DocInfo Doc_info;

  // output stream 
  ofstream c_integral_debug_ofstream;
  
}; // end_of_problem_class

template<class ELEMENT>
void StepProblem<ELEMENT> :: duplicate_corner_node()
{
  
  // Want elements in the upper region on the plate
  unsigned b         = Global_Physical_Variables::Enriched_region_dividing_boundary_id;
  unsigned region_id = Global_Physical_Variables::Enriched_region_upper_id;

  // number of elements in the upper enriched region
  unsigned nel = Bulk_mesh_pt->nboundary_element_in_region(b, region_id);

  // =========================================
  
  oomph_info << "Total number of nodes in mesh before duplication: "
	     << Bulk_mesh_pt->nnode() << "\n\n" << std::flush; 
  
  // map to keep track of nodes we've already duplicated;
  // map is original node -> new node
  std::map<Node*,Node*> existing_duplicate_node_pt;
  
  // loop over the elements above the plate
  for (unsigned e=0; e<nel; e++)
  {
    // get a pointer to this element
    FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_in_region_pt(b, region_id, e);

    // get the number of nodes in this element
    unsigned nnode = el_pt->nnode();

    // loop over nodes
    for(unsigned j=0; j<nnode; j++)
    {
      // get a pointer to the node
      Node* nod_pt = el_pt->node_pt(j);

      // check if we've found the corner node, i.e. the one which is on both the
      // region-dividing boundary and the corner boundary
      bool is_corner_node = nod_pt->is_on_boundary(Global_Physical_Variables::Corner_boundary_id);

      // QUEHACERES debug
      Vector<double> x(2);
      x[0] = nod_pt->x(0);
      x[1] = nod_pt->x(1);

      if(is_corner_node)
      {	    
	// Find this original node in the map; if we find it
	// it's already been duplicated earlier; in that case
	// use that duplicate rather than creating another one.
	std::map<Node*,Node*>::iterator it =
	  existing_duplicate_node_pt.find(nod_pt);
	    
	unsigned n_dim           = nod_pt->ndim();
	unsigned n_position_type = nod_pt->nposition_type();
	unsigned n_value         = nod_pt->nvalue();

	// check if we've already duplicated this node
	if (it == existing_duplicate_node_pt.end())
	{	    
	  // create a new node
	  el_pt->node_pt(j) = new BoundaryNode<Node>(n_dim, n_position_type, n_value);
         
	  // It has the same coordinate; hierher add history values too
	  // when implementing as Navier Stokes
	  for (unsigned i=0; i<n_dim; i++)
	  {
	    el_pt->node_pt(j)->x(i) = nod_pt->x(i);
	  }
         
	  // ...and the same values
	  for (unsigned i=0; i<n_value; i++)
	  {
	    el_pt->node_pt(j)->set_value(i, nod_pt->value(i));
	  }

	  // =================
	  
	  // add new node to the new boundary
	  // el_pt->node_pt(j)->add_to_boundary(new_boundary_id);
	  // calling this both calls the node pointers function to tell it it's on the boundary,
	  // and also updates the mesh's list of boundary nodes
	  Bulk_mesh_pt->add_boundary_node( Global_Physical_Variables::Corner_boundary_id,
					   el_pt->node_pt(j) );
	  Bulk_mesh_pt->add_boundary_node( Global_Physical_Variables::Enriched_region_dividing_boundary_id,
					   el_pt->node_pt(j) );
	  
	  // Get/set boundary coordinates
	  if ( nod_pt->boundary_coordinates_have_been_set_up() )
	  {
	    // get number of coordinates on the original plate boundary
	    unsigned ncoords = nod_pt->ncoordinates_on_boundary(b);
		
	    Vector<double> boundary_zeta(ncoords);

	    // get 'em
	    nod_pt->get_coordinates_on_boundary(
	      Global_Physical_Variables::Corner_boundary_id, boundary_zeta);
	    
	    // set 'em
	    el_pt->node_pt(j)->set_coordinates_on_boundary(
	      Global_Physical_Variables::Corner_boundary_id, boundary_zeta);

	    // get 'em
	    nod_pt->get_coordinates_on_boundary(
	      Global_Physical_Variables::Enriched_region_dividing_boundary_id, boundary_zeta);
	    
	    // set 'em
	    el_pt->node_pt(j)->set_coordinates_on_boundary(
	      Global_Physical_Variables::Enriched_region_dividing_boundary_id, boundary_zeta);

	  }
	  else
	  {
	    // hierher throw? (Doesn't happen at the moment, i.e. 
	    // when this diagnostic was finally commented out)
             
	    oomph_info << "No boundary coordinates have been set up"
		       << " for new local node " << j
		       << " at : "
		       << nod_pt->x(0) << " "
		       << nod_pt->x(1)
		       << std::endl;
	  }	    
	      
	  // add  new node to the bulk mesh
	  Bulk_mesh_pt->add_node_pt(el_pt->node_pt(j));
	      
	  // Keep track, add entry mapping from the old node to the new node
	  existing_duplicate_node_pt[nod_pt] = el_pt->node_pt(j);
	}
      }
    }
  }
  
  // QUEHACERES probably need this since we've fiddled the nodes on the plate boundary
  Bulk_mesh_pt->setup_boundary_element_info();

  // ================================================================================================
  // Now we've fixed the boundary elements, i.e. the ones who have a edge along the boundary,
  // we need to fix the adjacent elements who only have a single vertex on the boundary (which
  // are not considered boundary elements so don't get returned by boundary_element_in_region_pt(...) )

  // Now want element in the upper region on the corner boundary
  b = Global_Physical_Variables::Corner_boundary_id;
  
  // number of elements in the upper enriched region
  nel = Bulk_mesh_pt->nboundary_element_in_region(b, region_id);
  
  // loop over the elements above the plate
  for (unsigned e=0; e<nel; e++)
  {
    // get a pointer to this element
    FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_in_region_pt(b, region_id, e);

    // get the number of nodes in this element
    unsigned nnode = el_pt->nnode();

    // loop over nodes
    for(unsigned j=0; j<nnode; j++)
    {
      // get a pointer to the node
      Node* node_pt = el_pt->node_pt(j);

      // check if we've found the corner node, i.e. the one which is on both the
      // region-dividing boundary and the corner boundary
      bool is_corner_node = node_pt->is_on_boundary(
	Global_Physical_Variables::Enriched_region_dividing_boundary_id);
  
      // get coordinates of this node
      unsigned x_coord = node_pt->x(0);
      unsigned y_coord = node_pt->x(1);

      if(is_corner_node)
      {
	// check our map to find out the correct new node	      
	std::map<Node*,Node*>::iterator it =
	  existing_duplicate_node_pt.find(node_pt);

	// did we find it?
	if (it != existing_duplicate_node_pt.end())
	{
	  // update the node pointer
	  el_pt->node_pt(j) = existing_duplicate_node_pt[node_pt];

	  // should be all we need to do, the boundary info, coordinates, values etc. have already
	  // been sorted out when we created the new node
	}	
	else
	{
	  // if we don't find it something weird has happened, because we have a node
	  // on the original plate boundary that didn't get duplicated
	  oomph_info << "\nNode: (" << x_coord << ", " << y_coord << "); "
		     << "Weird tingz have happened here; I've found a node which \n"
		     <<" is apparently on the original corner boundary, but \n"
		     << "which didn't get duplicated\n\n";
	}
      }
    }    
  }

  // ====================================================

  // QUEHACERES debug, final check
  oomph_info << "Total number of nodes in mesh after duplication: "
	     << Bulk_mesh_pt->nnode() << "\n\n"; 
      
} // end duplicate_plate_nodes_and_add_boundary()


//==start_of_constructor==================================================
/// Constructor for StepProblem problem
//========================================================================
template<class ELEMENT>
StepProblem<ELEMENT>::StepProblem()
{  
  // Build the mesh
  Bulk_mesh_pt = Global_Physical_Variables::build_the_mesh<ELEMENT>
    (Global_Physical_Variables::Uniform_element_area);

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
  
  // Set element size limits
  Bulk_mesh_pt->max_element_size()    = 0.1;
  Bulk_mesh_pt->min_element_size()    = 1e-30;
  Bulk_mesh_pt->max_permitted_error() = 0.005;
  Bulk_mesh_pt->min_permitted_error() = 0.0;

  // Rescale to shrink domain
  if(Global_Physical_Variables::Scaling_factor_for_domain != 1.0)
  {
    unsigned nnode = Bulk_mesh_pt->nnode();
    for(unsigned i=0; i<nnode; i++)
    {
      Node* nod_pt = Bulk_mesh_pt->node_pt(i);
      nod_pt->x(0) = Global_Physical_Variables::Scaling_factor_for_domain*
	(nod_pt->x(0)-2.0) + 2.0;
      nod_pt->x(1) = Global_Physical_Variables::Scaling_factor_for_domain*
	(nod_pt->x(1));
    }
  }
  
  // Let's have a look at the boundary enumeration
  Bulk_mesh_pt->output_boundaries("boundaries.dat");

  // duplicate the corner node to allow for a dipole-like pressure discontinuity
  duplicate_corner_node();
  
  // Add sub-mesh
  add_sub_mesh(Bulk_mesh_pt);

  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {
    
    // Create element that stores the singular fct and its amplitude
    //---------------------------------------------------------------
    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
      new ScalableSingularityForNavierStokesElement<ELEMENT>;
    
    // Pass fct pointers:
    el_pt->unscaled_singular_fct_pt() = &Global_Physical_Variables::singular_fct;
    el_pt->gradient_of_unscaled_singular_fct_pt() =
      &Global_Physical_Variables::gradient_of_singular_fct;
    
    char filename[100];
    
    // Output contributions to the reciprocal integral which determines the amplitude C
    sprintf(filename,"%s/c_integral_contributions.dat",Doc_info.directory().c_str() );
    c_integral_debug_ofstream.open(filename);
    el_pt->c_boundary_integral_ofstream_pt() = &c_integral_debug_ofstream; 
    
    // Add to mesh
    Singular_fct_element_mesh_pt = new Mesh;
    Singular_fct_element_mesh_pt->add_element_pt(el_pt);
    add_sub_mesh(Singular_fct_element_mesh_pt);


    // Create face elements that compute contribution
    //-----------------------------------------------
    // to r_c
    //-------
    Face_mesh_for_singularity_integral_pt = new Mesh; 

    
    // Create face elements for flux jump
    //-----------------------------------
    // hierher Face_mesh_for_flux_jump_pt=new Mesh;
  }
  
  // Create face elements for imposition of BC
  Face_mesh_for_bc_pt = new Mesh;
  
  // Traction boundary condition 
  Traction_boundary_condition_mesh_pt = new Mesh;

  // Build the face elements
  create_face_elements();
  
  // Add to mesh
  add_sub_mesh(Face_mesh_for_bc_pt);
  add_sub_mesh(Traction_boundary_condition_mesh_pt);

  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {
    // hierher add_sub_mesh(Face_mesh_for_flux_jump_pt);
    
    // hierher currently not needed because the contributions from these elements are accumulated in
    // Singular_fct_element_mesh_pt but will need to add this back in 
    // when we optimise the code
    // add_sub_mesh(Face_mesh_for_singularity_integral_pt); 
  }

  // Build global mesh
  build_global_mesh();
  
  // Complete problem setup
  complete_problem_setup();
    
  // // hierher kill
  // {
  //  oomph_info << "imposing amplitude; removve thiis!\n";
   
  //  ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt=
  //   dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>
  //   (Singular_fct_element_mesh_pt->element_pt(0));
   
   
  //  // Change r_C so that C is assigned directly
  //  double imposed_amplitude=1.0;
  //  el_pt->impose_singular_fct_amplitude(imposed_amplitude);
  // }


  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " 
             << this->assign_eqn_numbers() 
             << std::endl;
  
} // end_of_constructor


//==start_of_complete======================================================
/// Set boundary condition, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::complete_problem_setup()
{ 
  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {   
    // Loop over the elements to set up element-specific
    // things that cannot be handled by constructor
    unsigned n_el=Bulk_mesh_pt->nelement();
    for (unsigned e=0; e<n_el; e++)
    {
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
	Bulk_mesh_pt->element_pt(e));

      // Tell bulk element which function computes the stress
      bulk_el_pt->stress_fct_pt() = &Global_Physical_Variables::get_stress;
      
      // Tell the bulk element about the singular fct
      bulk_el_pt->navier_stokes_sing_el_pt() =
	dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
	  Singular_fct_element_mesh_pt->element_pt(0));
      
      bulk_el_pt->exact_non_singular_fct_pt() =
	&Global_Physical_Variables::u_and_gradient_non_singular;
    }
   
    // hierher move this to create_face_elements as for the others
    // Flux jump elements
    // unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
    // for(unsigned e=0;e<n_element;e++)
    //  {
    //   // Upcast from GeneralisedElement to the present element
    //   PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* el_pt = dynamic_cast<
    //    PoissonWithSingularityFluxJumpFaceElement<ELEMENT>*>(
    //     Face_mesh_for_flux_jump_pt->element_pt(e));
     
    //   // Tell the element about the singular fct
    //   el_pt->set_poisson_sing_el_pt(
    //    dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
    //     Singular_fct_element_mesh_pt->element_pt(0)));
    //  }

    //BC elements
    unsigned n_element = Face_mesh_for_bc_pt->nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      NavierStokesWithSingularityBCFaceElement<ELEMENT>* el_pt = dynamic_cast<
	NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
	  Face_mesh_for_bc_pt->element_pt(e));
     
      // Tell the element about the singular fct
      el_pt->set_navier_stokes_sing_el_pt(
	dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	  Singular_fct_element_mesh_pt->element_pt(0)));
    }
  }
 
  // Apply bcs
  apply_boundary_conditions(); 
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::apply_boundary_conditions()
{
  // shorthand
  const unsigned Dim = Global_Physical_Variables::Dim;
  
  //---------------------------------------------------------------
  // DEFAULT SETUP FOR STANDARD NAVIER-STOKES PROBLEM WITHOUT LAGRANGE
  // MULTIPLIERS (CHANGED BELOW IF THE SINGULARITY IS INCLUDED)
  //---------------------------------------------------------------

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned num_bound = Bulk_mesh_pt->nboundary();
  for(unsigned ibound=0; ibound<num_bound; ibound++)
  {   
    bool pin_it = false;
   
    // Traction problem
    if (Global_Physical_Variables::Do_traction_problem)
    {
      if ( (ibound != Global_Physical_Variables::Outflow_boundary_id) &&
	   (ibound != Global_Physical_Variables::Enriched_region_lower_boundary_id) &&
	   (ibound != Global_Physical_Variables::Enriched_region_upper_boundary_id) )
      {
	pin_it = true;
      }
    }    
    // Otherwise pin everywhere but leave internal boundary alone
    else
    {
      if (ibound != Global_Physical_Variables::Enriched_region_lower_boundary_id &&
	  ibound != Global_Physical_Variables::Enriched_region_upper_boundary_id) 
      {
	pin_it = true;
      }
      // except if this is the outflow boundary and we're applying a "do nothing"
      // condition
      if(ibound == Global_Physical_Variables::Outflow_boundary_id
	 && CommandLineArgs::command_line_flag_has_been_set("--do_nothing_outflow_bc"))
      {
	pin_it = false;
      }

      //... but over-rule the lot if we use Lagrange multipliers
      if (CommandLineArgs::command_line_flag_has_been_set
	  ("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
      {
	pin_it = false;
      }
    }
    
    // Pin it?
    if (pin_it)
    { 
      unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
      for (unsigned inod=0; inod<num_nod; inod++)
      {
	// pin both components of the velocity
	for(unsigned i=0; i<Dim; i++)
	{
	  Bulk_mesh_pt->boundary_node_pt(ibound, inod)->pin(i);
	}
      }
    }
  } // end loop over boundaries
  
  // Now set boundary values
  for (unsigned ibound=0; ibound<num_bound; ibound++)
  {    
    // Leave internal boundary alone
    if (ibound != Global_Physical_Variables::Enriched_region_lower_boundary_id &&
	ibound != Global_Physical_Variables::Enriched_region_upper_boundary_id )
    {
      // loop over the nodes on this boundary
      unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
      for (unsigned inod=0; inod<num_nod; inod++)
      {
	// grab a pointer to this node
	Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound, inod);

	// position
	Vector<double> x(2);
	x[0] = nod_pt->x(0);
	x[1] = nod_pt->x(1);

	// get the boundary conditions for this location and boundary
	Vector<double> u(Dim);
	u = Global_Physical_Variables::u_BC(x, ibound);

	for(unsigned i=0; i<Dim; i++)
	{
	  nod_pt->set_value(i, u[i]);
	}
      }
    }
  }

  // Free values on boundary if Lagrange multiplier is used
  // to enforce bc for sum of FE soln and singuar fct
  // hierher
  // if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  if (CommandLineArgs::command_line_flag_has_been_set("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
  {
    // Now unpin nodal values where the bc conditions are enforced
    // by Lagrange multipliers to ensure that the sum of FE and singular
    // solution is correct
    unsigned nel = Face_mesh_for_bc_pt->nelement();
    for (unsigned e=0; e<nel; e++)
    {
      // Get element
      NavierStokesWithSingularityBCFaceElement<ELEMENT>* el_pt =
	dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
	  Face_mesh_for_bc_pt->element_pt(e));
     
      // Specify desired nodal values for compound solution
      unsigned nnod = el_pt->nnode();

      // matrix to store velocities at each boundary node
      DenseMatrix<double> nodal_boundary_value(nnod, Dim);

      // Unpin the FE part of the solution
      for (unsigned j=0; j<nnod; j++)
      {
	for(unsigned i=0; i<Dim; i++)
	{
	  el_pt->unpin_u_fe_at_specified_local_node(j, i);
	}
	// QUEHACERES this needs to change - if we have the enriched region,
	// we probably have 3 dofs for u,v,p, and 2 Lagrange multipliers for each of
	// u and v (for continuity across the jump and for the Dirichlet BCs), so 3+2x2 = 7 dofs.
	// leaving it for now as it won't be 3 anyways
	// -----
	// Now here's another subtle one! If this node stores three values
	// they come from: original Poisson dof; one Lagrange multiplier
	// to ensure continuity of solution across "jump"; one Lagrange
	// multiplier from enforcing Dirichlet boundary conditions. 
	// The latter two enforce the same thing, so we can only apply
	// the constraint once --> Pin the Lagrange multiplier that tries to 
	// enforce the Dirichlet BC
	unsigned nval = el_pt->node_pt(j)->nvalue();
	if (nval == 3)
        {
	  oomph_info << "hierher this shouldn't currently be used \n"
		     << " and ought to be double checked when \n"
		     << "flux jump comes back in...\n";
	  abort();

	  for(unsigned i=0; i<Dim; i++)
	  {
	    el_pt->pin_lagrange_multiplier_at_specified_local_node(j, i);
	  }
        }

	// if we're using Lagrange multipliers to enforce the BCs but we also want a
	// do nothing condition on the outflow, need to pin the lagrange multipliers on
	// the outflow boundary
	if(el_pt->node_pt(j)->is_on_boundary(Global_Physical_Variables::Outflow_boundary_id)
	  && !el_pt->node_pt(j)->is_on_boundary(Global_Physical_Variables::Top_boundary_id)
	    && !el_pt->node_pt(j)->is_on_boundary(Global_Physical_Variables::Bottom_boundary_id) )
	{
	  // pin the lagrange multplier for the x-component of the velocity on the outflow,
	  // u_x can be whatever it wants.
	  el_pt->pin_lagrange_multiplier_at_specified_local_node(j, 0);

	  // if we're doing parallel outflow (default) then we'll leave the lagrange multiplier for
	  // u_y unpinned, but if we're doing a do-nothing outflow condition then need to pin it too.
	  if(CommandLineArgs::command_line_flag_has_been_set("--do_nothing_outflow_bc") )
	  {
	    el_pt->pin_lagrange_multiplier_at_specified_local_node(j, 1);
	  }
	}
	
	Node* nod_pt = el_pt->node_pt(j);

	Vector<double> x(2);
	x[0] = nod_pt->x(0);
	x[1] = nod_pt->x(1);

	// find the boundaries this node is on
	std::set<unsigned>* boundaries_set;
	el_pt->node_pt(j)->get_boundaries_pt(boundaries_set);

	// assuming BCs are continuous so doesn't matter if node is on multiple
	// boundaries, just take the first
	std::set<unsigned>::iterator boundaries_iter = (*boundaries_set).begin();
	unsigned first_boundary = *boundaries_iter;
	
	// get velocity from boundary conditions at this point
	Vector<double> u(Dim);
	u = Global_Physical_Variables::u_BC(x, first_boundary);
	
	// assign to the matrix of nodal values
	for(unsigned i=0; i<Dim; i++)
	{
	  nodal_boundary_value(j,i) = u[i];
	}
      }
     
      // Tell the element about these nodal boundary values
      el_pt->set_nodal_boundary_values(nodal_boundary_value);
    }
  }

} // end set bc

//==start of set_values_to_singular_solution =============================
/// Function to assign the singular solution to all nodes of the mesh
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::set_values_to_singular_solution()
{
  // get the number of nodes in the mesh
  unsigned nnode = Bulk_mesh_pt->nnode();
  
  for(unsigned i=0; i<nnode; i++)
  {
    // get a pointer to this node
    Node* node_pt = Bulk_mesh_pt->node_pt(i);

    // get the position of this node
    Vector<double> x(2, 0.0);
    x[0] = node_pt->x(0);
    x[1] = node_pt->x(1);

    // get the singular solution at this point
    Vector<double> u(3, 0.0);

    u = Global_Physical_Variables::singular_fct(x);
    
    // assign the velocities
    node_pt->set_value(0, u[0]);
    node_pt->set_value(1, u[1]);

    // this is a bit naughty, if Lagrange multpliers are used it won't work!
    if(node_pt->nvalue() == 3)
    {
      node_pt->set_value(2, u[2]);
    }
  }
}

//== start of validate_stress ============================================
/// Function to validate the singular stress function by assigning the singular velocity field
// to the nodes of the mesh, then computing the "FE" stress via the navier-stokes
// helper functions and comparing the two
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::validate_stress()
{
  // assign \hat u_i to the nodal values
  set_values_to_singular_solution();
  
  // loop over all the elements in the mesh to compute the error in the stress
  const unsigned nel = Bulk_mesh_pt->nelement();

  char filename[100];

  sprintf(filename, "%s/error_in_singular_stress.dat", Doc_info.directory().c_str());
  
  // open the output file to record the error (tecplot format)
  ofstream stress_error_output(filename);
  
  // open the output file to record the error (plain format)
  sprintf(filename, "%s/error_in_singular_stress_plain.dat", Doc_info.directory().c_str());
  ofstream stress_error_output_plain(filename);
  
  // number of plot points per side
  unsigned nplot = 5;
  
  for(unsigned e=0; e<nel; e++)
  {
    // get a pointer to this element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // dimension of this element
    const unsigned dim = elem_pt->dim();
  
    // write the tecplot header for this element
    stress_error_output << elem_pt->tecplot_zone_string(nplot);
    
    //Set the Vector to hold local coordinates
    Vector<double> s(dim);
 
    // Loop over plot points    
    unsigned num_plot_points = elem_pt->nplot_points(nplot);
    for (unsigned iplot=0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      elem_pt->get_s_plot(iplot, nplot, s);

      // global coordinates
      Vector<double> x(dim, 0.0);
      
      // get interpolated global coordinates
      for(unsigned i=0; i<dim; i++)
      { 
	x[i] = elem_pt->interpolated_x(s,i);
	stress_error_output << x[i] << " ";
	stress_error_output_plain << x[i] << " ";
      }

      // -----------------------------------------
      // singular stuff
      // -----------------------------------------
      
      // singular solution at this knot (don't care about velocity, just need the pressure for the stress)
      Vector<double> u_sing(3, 0.0);
      u_sing = Global_Physical_Variables::singular_fct(x);

      // extract the singular pressure
      double p_sing = u_sing[2];
	
      // get the singular velocity gradient
      DenseMatrix<double> du_dx_sing(dim, dim, 0.0);    
      du_dx_sing = Global_Physical_Variables::gradient_of_singular_fct(x);

      // compute the singular strain rate
      DenseMatrix<double> strain_rate_sing(dim, dim, 0.0);

      for(unsigned i=0; i<dim; i++)
      {
	for(unsigned j=0; j<dim; j++)
	{
	  strain_rate_sing(i,j) = 0.5*(du_dx_sing(i,j) + du_dx_sing(j,i));
	}
      }

      // get the singular stress
      DenseMatrix<double> stress_sing(dim, dim, 0.0);
      stress_sing = Global_Physical_Variables::get_stress(strain_rate_sing, p_sing);

      // -----------------------------------------
      // "FE" stuff
      // -----------------------------------------

      // FE pressure
      double p_fe = elem_pt->interpolated_p_nst(s);
      
      // compute the "FE" strain-rate
      DenseMatrix<double> strain_rate_fe(dim, dim, 0.0);

      elem_pt->strain_rate(s, strain_rate_fe);

      // compute the "FE" stress
      DenseMatrix<double> stress_fe(dim, dim, 0.0);
      stress_fe = Global_Physical_Variables::get_stress(strain_rate_fe, p_fe);

      // -----------------------------------------
      // Error
      // -----------------------------------------
	
      // compute the error
      for(unsigned i=0; i<dim; i++)
      {
	for(unsigned j=0; j<dim; j++)
	{
	  // compute the error between the interpolated "FE" stress and the exact singular stress
	  double error = stress_fe(i,j) - stress_sing(i,j);

	  // output it
	  stress_error_output       << error << " ";
	  stress_error_output_plain << error << " ";
	}
      }

      stress_error_output       << std::endl;
      stress_error_output_plain << std::endl;
      
    } // end loop over plot point

    stress_error_output       << std::endl;    
    
    // Write tecplot footer (e.g. FE connectivity lists)
    elem_pt->write_tecplot_zone_footer(stress_error_output, nplot);
  } // end loop over elements
  
  // done, close the output file
  stress_error_output.close();
  stress_error_output_plain.close();
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::doc_solution()
{

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts=10;

  // Output solution
  sprintf(filename, "%s/soln%i.dat", Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts); // output_paraview
  some_file.close();


  // Output solution just using vertices so we can see the mesh
  sprintf(filename, "%s/coarse_soln%i.dat", Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  npts=2;
  Bulk_mesh_pt->output(some_file, npts);
  some_file.close();
  

  // Plot "extended solution" showing contributions; also work out
  // average element size
  double av_el_size=0.0;
  sprintf(filename,"%s/extended_soln%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  
  unsigned nel=Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<nel; e++)
  {
    npts=10;
    ELEMENT* el_pt= dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    
    el_pt->output_with_various_contributions(some_file, npts);
    av_el_size += el_pt->size();
  }
  some_file.close();
  av_el_size/=double(nel);

  // Get error
  double error,norm; 
  sprintf(filename, "%s/error%i.dat", Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename); 
  Bulk_mesh_pt->compute_error(some_file,
                              Global_Physical_Variables::u_exact,
                              error,norm);
  some_file.close();

  // Doc error norm:
  oomph_info << "\n av el size, av h, Ndof, # bulk els, Norm of error    : "   
             << av_el_size << " " 
             << sqrt(av_el_size) << " " 
             << ndof() << " " 
             << Bulk_mesh_pt->nelement() << " " 
             << sqrt(error) << std::endl;
  oomph_info << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
  oomph_info << std::endl;
  
  // Exact solution
  sprintf(filename, "%s/exact_soln%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  unsigned nplot=5;
  Bulk_mesh_pt->output_fct(some_file,nplot,
                           Global_Physical_Variables::u_exact);
  some_file.close();

  // Non-singular part of the solution
  sprintf(filename, "%s/non_singular_part_of_soln%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  nplot=5;
  Bulk_mesh_pt->output_fct(some_file,nplot,
                           Global_Physical_Variables::u_non_singular);
  some_file.close();

  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_subtract_singularity"))
  {
    // hierher put into one file and flush trace file for value of C
    sprintf(filename,"%s/suivi_C.dat",Doc_info.directory().c_str());
    
    some_file.open(filename,std::ios::out | std::ios::app);
    
    some_file << sqrt(av_el_size) <<" " 
              << dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
		Singular_fct_element_mesh_pt->element_pt(0)
		)->amplitude_of_singular_fct()<<endl;
    some_file.close();
  }

  // DoubleVector r;
  // CRDoubleMatrix jac;

  // sprintf(filename,"%s/most_recent_jacobian%i.dat",Doc_info.directory().c_str(),Doc_info.number());
  // some_file.open(filename);
  // oomph_info << "SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
  // get_jacobian(r,jac);
  // oomph_info << "DONE SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
  // jac.sparse_indexed_output(some_file);

  // sprintf(filename,"%s/suivi_derivee.dat",Doc_info.directory().c_str());
  // some_file.open(filename,std::ios::out | std::ios::app);
  // // We're looking for a bulk element that has a node on the
  // // the corner of the backward step; all candidate elements live on
  // // boundary 4
  // unsigned b = 4;
  // unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
  // // Loop over the bulk elements adjacent to boundary b
  // for(unsigned e=0;e<n_element;e++)
  // {
  //  // Get pointer to the bulk element that is adjacent to boundary b
  //  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
  //   Bulk_mesh_pt->boundary_element_pt(b,e));
  //  unsigned sortie = 0;
     
  //  unsigned nnod=bulk_elem_pt->nnode();
  //  for (unsigned j=0;j<nnod;j++)
  //  {
  //   Node* nod_pt=bulk_elem_pt->node_pt(j);
      
  //   // Are we at the corner?
  //   double distance=sqrt((nod_pt->x(0)-Global_Physical_Variables::L_up)*
  //                         (nod_pt->x(0)-Global_Physical_Variables::L_up)+
  //                         (nod_pt->x(1)-0.0)*(nod_pt->x(1)-0.0));
  //   double tol=1.0e-10;
  //   if (distance<tol)
  //   {
  //    Vector<double> outward_vector(2);
  //    outward_vector[0] = 1.0/sqrt(2.0);
  //    outward_vector[1] = 1.0/sqrt(2.0);
  //    Vector<double> Flux(1);
  //    bulk_elem_pt->get_flux(outward_vector,Flux);
  //    some_file << sqrt(av_el_size) << " " << Flux[0] << endl;
  //    sortie = 1;
  //    break;
  //   }
  //  }
  //  if (sortie==1) break;
  // }
  // some_file.close();


  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_subtract_singularity"))
  {
    oomph_info 
      << "Amplitude of singular function: "
      << dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
	Singular_fct_element_mesh_pt->element_pt(0))->
      amplitude_of_singular_fct() << std::endl;

    // Output face elements used to compute amplitude of singularity
    sprintf(filename,"%s/integrand_elements%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());
    some_file.open(filename);
    nel=Face_mesh_for_singularity_integral_pt->nelement();
    for (unsigned e=0;e<nel;e++)
    {
      NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt=
	dynamic_cast<NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>
	(Face_mesh_for_singularity_integral_pt->element_pt(e));
      el_pt->get_contribution_integral(some_file); //output(some_file,npts);
    }
    some_file.close();
  }
  
  
  // Output flux elements
  sprintf(filename,"%s/flux_elements%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  
  nel = Traction_boundary_condition_mesh_pt->nelement();
  for (unsigned e=0; e<nel; e++)
  {
    npts = 5;
    Traction_boundary_condition_mesh_pt->finite_element_pt(e)->output(some_file,npts);
  }
  some_file.close();

  // ============================================================
  // plot solution along theta=0 and theta=pi/2 from the corner
  // ============================================================
  
  unsigned N = 1001;
  Vector<Vector<double> > coords_theta_0(N);
  Vector<Vector<double> > coords_theta_pi_2(N);

  // assign coordinates
  for(unsigned i=0; i<N; i++)
  {
    // make space for 2 coordinates
    coords_theta_0[i].resize(2);
    coords_theta_pi_2[i].resize(2);
    
    coords_theta_0[i][0] = Global_Physical_Variables::L_up +
      (double)(i)*Global_Physical_Variables::L_down/(double)(N);
    coords_theta_0[i][1] = 0.0;

    coords_theta_pi_2[i][0] = Global_Physical_Variables::L_up;
    coords_theta_pi_2[i][1] = (double)(i)*Global_Physical_Variables::H_up/(double)(N);
  }

  // create the visualisers for the two angles
  LineVisualiser line_visualiser_theta_0(   Bulk_mesh_pt, coords_theta_0);
  LineVisualiser line_visualiser_theta_pi_2(Bulk_mesh_pt, coords_theta_pi_2);

  // open file
  sprintf( filename, "%s/soln%i_theta=0.dat", Doc_info.directory().c_str(), Doc_info.number() );

  some_file.open(filename);

  // print out the points
  line_visualiser_theta_0.output(some_file);

  // done with this file
  some_file.close();

  // open file
  sprintf( filename, "%s/soln%i_theta=pi_2.dat", Doc_info.directory().c_str(), Doc_info.number() );

  some_file.open(filename);

  // print out the points
  line_visualiser_theta_pi_2.output(some_file);

  // done with this file
  some_file.close();
  
} // end_of_doc_solution




//=======================================================================
// hierher
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::check_condition_number()
{
 
  // Fake r_c equation?
  if (Global_Physical_Variables::Problem_type_for_check_condition_number == 2)
  {
    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>
      (Singular_fct_element_mesh_pt->element_pt(0));
   
    // Change r_C so that C is assigned directly
    double imposed_amplitude = Global_Physical_Variables::singular_amplitude_for_debug;
    el_pt->impose_singular_fct_amplitude(imposed_amplitude);
  }
  else // QUEHACERES probably remove this at some point
  {


    // Storage for the Jacobian
    CRDoubleMatrix matrix;
 
    // Storage for the (to be unused) residual vector
    DoubleVector residual;
 
    // Storage for the eigenvector associated with the largest eigenvalue
    DoubleVector eigenvector_max;
 
    // Storage for the eigenvector associated with the smallest eigenvalue
    DoubleVector eigenvector_min;
 
    // Get the Jacobian matrix and residual vector
    get_jacobian(residual,matrix);
 
    // Tell the user
    oomph_info << "\n=================================================="
	       << "\nBeginning eigenvalue/singular value calculation..."
	       << "\n=================================================="
	       << std::endl;
 
    // Storage for the largest eigenvalue/singular value
    double lambda_max = power_method(&matrix,eigenvector_max);
 
    // Storage for the smallest eigenvalue/singular value
    double lambda_min = inverse_power_method(&matrix,eigenvector_min);
 
    // Set the output precision
    oomph_info.stream_pt()->precision(10);
 
    // If we're just computing eigenvalues
    if (!PowerIterationHelperNamespace::Compute_singular_values)
    {
      // Output the eigenvalues
      oomph_info << "\nLargest eigenvalue: " << lambda_max
		 << "\nSmallest eigenvalue: " << lambda_min
		 << std::endl;
    }
    // If we're computing singular values instead
    else
    {
      // Output the singular values
      oomph_info << "\nLargest singular value: " << lambda_max
		 << "\nSmallest singular value: " << lambda_min
		 << "\nMatrix size and condition number: " 
		 << residual.nrow() << " " << lambda_max/lambda_min
		 << std::endl;
    }
 
    // Create a method of outputting the Jacobian
    // std::ofstream matrix_file("jacobian.dat");
 
    // Dump the matrix so we can check that we're doing the right thing...
    //matrix.sparse_indexed_output(matrix_file);
 
    // Close the file
    // matrix_file.close();
 
    // Output the eigenvector associated with the largest eigenvalue
    eigenvector_max.output("eigenvector_max.dat");
 
    // Also output the eigenvector associated with the smallest eigenvalue
    eigenvector_min.output("eigenvector_min.dat"); 
  } 

}


//=======================================================================
// hierher
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::check_residual_for_exact_non_singular_fe_soln()
{ 
  oomph_info << std::endl << std::endl;

  // Assign exact non-singular solution to nodal values
  unsigned nnod=Bulk_mesh_pt->nnode();
  
  for (unsigned j=0; j<nnod; j++)
  {
    Node* nod_pt = Bulk_mesh_pt->node_pt(j);
    Vector<double> x(2);
    x[0]=nod_pt->x(0);
    x[1]=nod_pt->x(1);

    Vector<double> u(3);
    Global_Physical_Variables::u_non_singular(x, u); // hierher break deliberately...

    for(unsigned i=0; i<2; i++)
    {
      nod_pt->set_value(i, u[i]);
    }
  }

  // STAGE 1: Evaluate the integral with exact non-singular part
  //------------------------------------------------------------

  // Check the computation of the integral that determines C.
  // If we set u_fe and its gradient to the non-singular part of
  // the solution this should be close to zero!
  unsigned nel = Face_mesh_for_singularity_integral_pt->nelement();
  
  for (unsigned e=0; e<nel; e++)
  {
    NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt =
      dynamic_cast<NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>
      (Face_mesh_for_singularity_integral_pt->element_pt(e));
    
    el_pt->exact_non_singular_fct_pt() = 
      &Global_Physical_Variables::u_and_gradient_non_singular;
  }
 
 
  ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
    dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>
    (Singular_fct_element_mesh_pt->element_pt(0));
  
  unsigned n_dof = el_pt->ndof();
  Vector<double> residuals(n_dof, 0.0);
  
  el_pt->fill_in_contribution_to_residuals(residuals);
 
  oomph_info << "Residual if we assign the exact non-singular part of the "
	     << "solution for area "
	     << Global_Physical_Variables::Uniform_element_area 
	     << " and number of elements: " 
	     << Bulk_mesh_pt->nelement() << " : " 
	     << fabs(residuals[0]) << endl;
 
  doc_solution();
  Doc_info.number()++;

  // Reset
  for (unsigned e=0; e<nel; e++)
  {
    NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt =
      dynamic_cast<NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>      
      (Face_mesh_for_singularity_integral_pt->element_pt(e));
    
    el_pt->exact_non_singular_fct_pt() = 0;
  }


  // STAGE 2: Compute integral from (previously assigned) exact nodal values
  //------------------------------------------------------------------------
  residuals.clear();
  residuals.resize(n_dof);
  el_pt->fill_in_contribution_to_residuals(residuals);
 
  oomph_info << "Residual if we assign the exact non-singular part of the "
	     << "solution to nodal values for area "
	     << Global_Physical_Variables::Uniform_element_area 
	     << " and number of elements: " 
	     << Bulk_mesh_pt->nelement() << " : " 
	     << fabs(residuals[0]) << endl;
 
  doc_solution();
  Doc_info.number()++;
 
}



//=======================================================================
// Impose amplitude of singular fct via fake eqn
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::impose_amplitude_runs()
{

  ofstream some_file;
  ofstream some_file_error;
  ofstream some_file_integral_of_error;
  char filename[100];

  sprintf(filename,"%s/rc_vs_amplitude.dat",
          Doc_info.directory().c_str());
  some_file.open(filename);

  sprintf(filename,"%s/integral_of_Z2_error.dat",
          Doc_info.directory().c_str());
  some_file_integral_of_error.open(filename);


  // Loop over imposed amplitudes
  double a_min=-1.0;
  unsigned nstep_aux=4;
  double d_ampl=(1.0-a_min)/double(nstep_aux);
  unsigned nstep=2*nstep_aux;
  double imposed_amplitude=a_min;
  for(unsigned incr = 0; incr < nstep; incr++)
  {
    double integral_of_Z2error = 0.0;
    
    sprintf(filename,"%s/elementwise_Z2error%i.dat",
	    Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file_error.open(filename);
   

    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>
      (Singular_fct_element_mesh_pt->element_pt(0));

    // Change r_C so that C is assigned directly
    el_pt->impose_singular_fct_amplitude(imposed_amplitude);
    some_file << imposed_amplitude << " ";
    some_file_integral_of_error << imposed_amplitude << " ";
   
    // Solve the eqns
    this->newton_solve();

    oomph_info << "have solved for imposed_amplitude = "
	       << imposed_amplitude << " output soln in "
	       << Doc_info.number() << std::endl;

    doc_solution();
    Doc_info.number()++;

    // Check residual (computed by actual integral!) in singular element  
    el_pt->dont_impose_singular_fct_amplitude();
  
    unsigned n_dof=el_pt->ndof();
    Vector<double> residuals(n_dof,0.0);  
    el_pt->fill_in_contribution_to_residuals(residuals);

    some_file <<  residuals[0] << endl;

    unsigned n = Bulk_mesh_pt->nelement();
    Vector<double> elementwise_Z2error(n);
    Mesh* mesh_copy_pt = Bulk_mesh_pt;
  
    // Use actual value without normalisation!
    Z2ErrorEstimator* z2_pt=dynamic_cast<Z2ErrorEstimator*>(
      Bulk_mesh_pt->spatial_error_estimator_pt());
    double backup=z2_pt->reference_flux_norm();
    z2_pt->reference_flux_norm() = 1.0;
  
    //Bulk_mesh_pt->spatial_error_estimator_pt()->
    z2_pt->get_element_errors(mesh_copy_pt,elementwise_Z2error);

    // Reset
    z2_pt->reference_flux_norm()=backup;
    
    // Plot constant pressure th'out element
    unsigned nplot = 3;
    for(unsigned i=0; i<n;i++)
    {
      FiniteElement* el_pt=
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
      unsigned num_plot_points=el_pt->nplot_points(nplot);
      Vector<double> s(2);
      some_file_error << el_pt->tecplot_zone_string(nplot);
      for(unsigned j=0; j<num_plot_points;j++)
      {
	el_pt->get_s_plot(j,nplot,s);
	Vector<double> x(2);
	x[0]=el_pt->interpolated_x(s,0);
	x[1]=el_pt->interpolated_x(s,1);
	some_file_error << x[0] << " " << x[1] << " " 
			<< elementwise_Z2error[i] << endl;
      }
      el_pt->write_tecplot_zone_footer(some_file_error,nplot);
      integral_of_Z2error += elementwise_Z2error[i]*el_pt->size();
    }

    some_file_integral_of_error << integral_of_Z2error << endl;
    some_file_error.close();
    imposed_amplitude += d_ampl;
  }
  
  some_file_integral_of_error.close();
  some_file.close();
}

// QUEHACERES function to validate the maths - we'll assume u is correct,
// and check dudx via finite-diff of u
void validate_dudx_via_finite_diff(DocInfo doc_info)
{
  // number of grid points in each direction
  const int N = 200;

  // step size
  double hx = 2.0/(double)(N);
  double hy = 1.0/(double)(N);
  
  char filename[100];
  ofstream dudx_validation;
  sprintf(filename, "%s/dudx_validation.dat", doc_info.directory().c_str() );
  dudx_validation.open(filename);
  
  // loop in x
  for(int i=0; i<N; i++)
  {
    // loop in y
    for(int j=0; j<N; j++)
    {
      // get grid coordinates
      Vector<double> x(2, 0.0);
      x[0] =  1.0 + 2.0*(double)(i)/(double)(N);
      x[1] = -0.5 + (double)(j)/(double)(N);

      dudx_validation << x[0] << " "
		      << x[1] << " ";
      
      // get x-coordinates of left and right neighbouring points
      Vector<double> x_left(2);
      x_left[0] = 1.0 + 2.0*(double)(i-1)/(double)(N);
      x_left[1] = x[1];

      Vector<double> x_right(2);
      x_right[0] = 1.0 + 2.0*(double)(i+1)/(double)(N);
      x_right[1] = x[1];

      Vector<double> x_up(2);
      x_up[0] = x[0];
      x_up[1] = -0.5 + (double)(j+1)/(double)(N);
      
      Vector<double> x_down(2);
      x_down[0] = x[0];
      x_down[1] = -0.5 + (double)(j-1)/(double)(N);
      
      // analytic solution and gradient
      Vector<double> u;
      Vector<double> u_left;
      Vector<double> u_right;
      Vector<double> u_down;
      Vector<double> u_up;
      
      DenseMatrix<double> dudx;
      DenseMatrix<double> dudx_dummy;
      
      // finite-diff gradient
      DenseMatrix<double> dudx_fd(2,2);
            
      // get analytic solutions
      Global_Physical_Variables::singular_fct_and_gradient(x,       u,       dudx);
      Global_Physical_Variables::singular_fct_and_gradient(x_left,  u_left,  dudx_dummy);
      Global_Physical_Variables::singular_fct_and_gradient(x_right, u_right, dudx_dummy);
      Global_Physical_Variables::singular_fct_and_gradient(x_up,    u_up,    dudx_dummy);
      Global_Physical_Variables::singular_fct_and_gradient(x_down,  u_down,  dudx_dummy);

      // output shit
      for(unsigned k=0; k<2; k++)
      {
	for(unsigned l=0; l<2; l++)
	{
	  dudx_validation << dudx(k,l) << " ";
	}
      }
      for(unsigned k=0; k<2; k++)
      {
	dudx_validation << u[k] << " ";
      }
      
      // do the central finite-diff
      // --------------------

      for(unsigned k=0; k<2; k++)
      {
	// duk_dx
	dudx_fd(k,0) = (u_right[k] - u_left[k]) / (2*hx);
	// duk_dy
	dudx_fd(k,1) = (u_up[k]    - u_down[k]) / (2*hy);

	dudx_validation << dudx_fd(k,0) << " "
			<< dudx_fd(k,1) << " ";
      }
      for(unsigned k=0; k<2; k++)
      {
	// error
	dudx_validation << dudx_fd(k,0) - dudx(k,0) << " "
			<< dudx_fd(k,1) - dudx(k,1) << " ";
      }
      dudx_validation << std::endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for backward step with impedance outflow bc
//=====================================================================
int main(int argc, char **argv)
{

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // output directory for results
  CommandLineArgs::specify_command_line_flag("--dir", &Global_Physical_Variables::dir);
  
  // finite-diff validation of dudx maths
  CommandLineArgs::specify_command_line_flag("--validate_dudx");

  // validation of singular stress
  CommandLineArgs::specify_command_line_flag("--validate_singular_stress");
  
  // Don't subtract off singularity
  CommandLineArgs::specify_command_line_flag("--dont_subtract_singularity");

  // don't apply any boundary conditions at the outflow
  CommandLineArgs::specify_command_line_flag("--do_nothing_outflow_bc");
  
  // Suppress singular term in exact solution
  CommandLineArgs::specify_command_line_flag("--suppress_sing_in_exact_soln");

  // Do runs where we impose the amplitude via the "fake" r_c equation
  CommandLineArgs::specify_command_line_flag(
    "--check_overall_behaviour_when_setting_amplitude_via_fake_rc_equation");

  // Check condition number and specify number that identifies which
  // case we're doing: 0: real problem; 1: FE-only; 2: fake r_c equation
  CommandLineArgs::specify_command_line_flag(
    "--check_condition_number",
    &Global_Physical_Variables::Problem_type_for_check_condition_number);

  // Check the r_c equation
  CommandLineArgs::specify_command_line_flag("--check_rc_equation");

  // Scaling factor for domain (defaults to 1.0)
  CommandLineArgs::specify_command_line_flag("--scaling_factor_for_domain",
    &Global_Physical_Variables::Scaling_factor_for_domain);
 
  // Use Dirichlet boundary conditions (allegedly makes condition
  // number (even) worse)
  CommandLineArgs::specify_command_line_flag(
    "--enforce_dirichlet_bcs_by_lagrange_multipliers");

  // Use Dirichlet boundary conditions (allegedly makes condition
  // number (even) worse)
  CommandLineArgs::specify_command_line_flag("--use_dirichlet_bcs");
  
  // uniform target element area
  CommandLineArgs::specify_command_line_flag("--element_area",
    &Global_Physical_Variables::Uniform_element_area);

  // high resolution element area in region around the corner
  CommandLineArgs::specify_command_line_flag("--high_res_element_area",
    &Global_Physical_Variables::High_res_element_area);
  
  // Constant in singular solution
  CommandLineArgs::specify_command_line_flag("--constant_in_singular_solution",
    &Global_Physical_Variables::Constant_on_singular_boundary);

  // fake singular amplitude for debug
  CommandLineArgs::specify_command_line_flag("--sing_amplitude",
    &Global_Physical_Variables::singular_amplitude_for_debug);
  
  // // Max. number of adaptations in code
  // unsigned max_adapt = 0;
  // CommandLineArgs::specify_command_line_flag(
  //  "--max_adapt",&max_adapt);

  // Parse command line
  CommandLineArgs::parse_and_assign(); 
 
  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();
  
  // hierher tidy
  if (CommandLineArgs::command_line_flag_has_been_set("--use_dirichlet_bcs"))
  {
    Global_Physical_Variables::Do_traction_problem = false;
  }
    
  // Build the problem with 
  StepProblem<ProjectableTaylorHoodElement<MyTNavierStokesElement<2,3> > > problem;

  // Set up doc info
  // ---------------
  
  // Set output directory
  problem.doc_info()->set_directory(Global_Physical_Variables::dir.c_str());
 
  // Step number
  problem.doc_info()->number() = 0;

  // die if no folder to output results (saves doing a tedious calc only to find there's no output!)
  problem.doc_info()->enable_error_if_directory_does_not_exist();

  
  if (CommandLineArgs::command_line_flag_has_been_set("--validate_dudx"))
  {
    validate_dudx_via_finite_diff(*problem.doc_info());
    exit(0); 
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--validate_singular_stress"))
  {
    problem.validate_stress();
    exit(0);
  }
  
  // Check condition number
  //=======================
  if (CommandLineArgs::command_line_flag_has_been_set("--check_condition_number"))
  {
    
    // Sorry this is a mess!
    if ((Global_Physical_Variables::Problem_type_for_check_condition_number == 1)
	&& (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity")))
    {
      oomph_info << "Please specify --dont_subtract_singularity\n";
      abort();
    }

    // Check condition number
    problem.check_condition_number();
    
    // Solve the bloody thing
    problem.newton_solve();
    problem.doc_solution();
    exit(0);
  }

  // Impose amplitude via fake eqn?
  //================================
  if (CommandLineArgs::command_line_flag_has_been_set
      ("--check_overall_behaviour_when_setting_amplitude_via_fake_rc_equation"))
  {
    problem.impose_amplitude_runs();
    exit(0);
  }
  

  // Check residual r_c
  //===================
  if (CommandLineArgs::command_line_flag_has_been_set
      ("--check_rc_equation"))
  {
    problem.check_residual_for_exact_non_singular_fe_soln();
    exit(0);
  }
  


  // Doing normal run: Just solve the bloody thing
  problem.newton_solve();
  problem.doc_solution();
  problem.doc_info()->number()++;
  

  // cout << "not done yet\n";
  // abort();

  

  //##################



  // bool doc_resid=false;

  // Global_Physical_Variables::shrink_the_mesh=false;
  // Global_Physical_Variables::scale=0.5;

  // // Solve, refine uniformly and keep going
  //  for (unsigned i=0;i<max_adapt;i++)
  //   {

  //    if (doc_resid)
  //     {
  //      if(i==max_adapt-1)
  //      {
  //       problem.impose_amplitude_runs();
  //       problem.check_residual_for_exact_non_singular_fe_soln();        
  //      }

  //     }
  //    else
  //     {
  //      // Solve the bloody thing
  //      problem.newton_solve();
  //      problem.doc_solution();
  //      problem.doc_info()->number()++;
  //     }


  //    if (i!=(max_adapt-1)) problem.refine_uniformly();
  //   }


}
