//  Copyright (c) 2014 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define HPX_LIMIT 8

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>
//#include <hpx/hpx_fwd.hpp>
//#include <hpx/lcos/when_all.hpp>
//#include <hpx/include/local_lcos.hpp>
//#include <hpx/util/unwrapped.hpp>

#define UW hpx::util::unwrapped
#define UW2 hpx::util::unwrapped2

// Simulates gas flow inside of a tube that is surrounded by an environment
// that maintains a constant temperature everywhere. Entropy is also assumed to
// be constant everywhere. Our governing equations are a simplified form of
// Euler laws for gas dynamics:
//
//      d(rho)/dt + d(rho*v)/dx = 0         (Continuity Law)
//      d(rho*v)/dt + d(rho*v^2+P)/dx = 0   (Conservation of Momentum
//
//      P = gamma^2*rho                     (Ideal Gas Law for Isothermal Flow)
//
// where rho is density, v is velocity, rho*v is momentum, P is pressure and
// gammas is the pressure constant. Density and momentum are conserved; the rest
// are primitive variables.
//
// This problem is representative of a class of non-linear hyperbolic PDEs (encounted
// frequently in the study of fluid dynamics) known as nonlinear conservation laws.
//
// We solve this system using explicit finite differencing methods and an
// operator-splitting approach; first advection is performed, and then source
// terms. We use a simple donor-cell upwind scheme for advection; it is 1st
// order in time and space. 
//
// Currently, the initial conditions are a gaussian density distribution at rest.
// Reflecting boundary conditions are used (e.g. closed tube). According to my
// textbook, the pressure force is supposed to cleave the gaussian blob in half,
// sending one dense concentration of gas towards each boundary. The blobs reflect
// off the boundaries and smash into each other, reforming a smaller density
// peak at the location of the initial gaussian peak. Lather, rinse, repeat.
//
// References: "Numerical Methods for Conservation Laws" (LeVeque), Lecture Notes
// of C.P. Dullemond/R. Kuiper.
 
///////////////////////////////////////////////////////////////////////////////

typedef std::size_t coord;

enum cell_type
{
    INTERIOR = 0,  // Interior grid point. 

    LEFT     = (1<<0),
    RIGHT    = (1<<1),
    BOUNDARY = (1<<2),
    OUTSIDE  = (1<<3),

    LEFT_BOUNDARY   = (LEFT | BOUNDARY),  // Leftmost grid point.
    LEFT_OUTSIDE    = (LEFT | OUTSIDE),   // Outside the domain, to the left.
    RIGHT_BOUNDARY  = (RIGHT | BOUNDARY), // Rightmost grid point.
    RIGHT_OUTSIDE   = (RIGHT | OUTSIDE)   // Outside the domain, to the right.
};

struct state 
{
    double rho;     // Density.
    double mom;     // Momentum.
    cell_type type;

    state(double r, double m, cell_type t = INTERIOR)
      : rho(r), mom(m), type(t) {}
};

// Create a future for a state that is currently ready.
hpx::future<state> ready(cell_type t) {
    return hpx::make_ready_future(state(0.0, 0.0, t));
}

struct solver
{
    typedef std::vector<hpx::shared_future<state> > space; 
    typedef std::vector<space> spacetime; 
    typedef std::vector<hpx::shared_future<double> > timestep_sizes;

    // Number of substeps. We have two here (due to operator splitting);
    // Advection, then source terms (just pressure for us).
    static constexpr coord substeps = 2;

    // CFL parameters.
    static constexpr double init_dt           = 0.1;  // Initial timestep size.
    static constexpr double dt_growth_limiter = 1.25; // Max timestep size growth per step.
    static constexpr double C                 = 0.4;  // Courant number.

    static constexpr double dx = 1.0; // Grid spacing.

    static constexpr double gamma = 7.0/5.0; // Pressure constant.

  private:
    coord nx; // Number of grid points; 0-indexed.
    coord nT; // Number of steps; 1-indexed (T=0 is initial state).

    spacetime U; // U[t][i] is the state of position i at substep t.

    timestep_sizes CFL_dt; // CFL_dt[T] is the timestep size (dt) for step T. 

  public:
    solver(coord nx_, coord nT_) : nx(nx_), nT(nT_ + 1), U(nT*substeps), CFL_dt(nT) {
        for (space& s : U) s.resize(nx);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Operators

    // Compute velocity at left cell interfaces by arithmetic mean. 
    static double velocity(state left, state middle) {
        double v_left   = left.mom/left.rho; 
        double v_middle = middle.mom/middle.rho;
        return 0.5*(v_middle+v_left);
    };

    // Compute the numerical flux at the left interface using the donor-cell scheme.
    static state flux(state left, state middle) {
        if ((middle.type & LEFT) || (middle.type == RIGHT_OUTSIDE)) {
            return state(0.0, 0.0, middle.type);
        }

        double v = velocity(left, middle);

        if (v > 0.0) {
            double fluxrho = left.rho*v;
            double fluxmom = (std::pow(left.mom, 2))/left.rho;
            return state(fluxrho, fluxmom, middle.type); 
        } else {
            double fluxrho = middle.rho*v;
            double fluxmom = (std::pow(middle.mom, 2))/middle.rho;
            return state(fluxrho, fluxmom, middle.type); 
        }
    }

    // Compute the pressure force (dP/dx) by central difference approximation.
    // We assume that the pressure outside the tube is equal to the pressure at
    // the boundary (e.g. reflecting boundary conditions); we're not using
    // ghost zones so we must add a 1/2 factor for the boundaries.
    static double pressure(state left, state middle, state right) {
        if      (middle.type == LEFT_BOUNDARY)
            return 0.25 * (gamma-1.0)*(right.rho-middle.rho);
        else if (middle.type == RIGHT_BOUNDARY)
            return 0.25 * (gamma-1.0)*(middle.rho-left.rho);
        else 
            return 0.5 * (gamma-1.0)*(right.rho-left.rho);
    };

    // Update for donor-cell advection.
    static state advection(double dt, state left, state middle, state right) {
        state fl_middle = flux(left, middle);  // F[t-1][i]
        state fl_right  = flux(middle, right); // F[t-1][i+1]

        double rho = middle.rho - (dt/dx)*(fl_right.rho-fl_middle.rho);
        double mom = middle.mom - (dt/dx)*(fl_right.mom-fl_middle.mom);
        return state(rho, mom, middle.type); // U[t+1][i] 
    };

    // Update for sources.
    static state sources(double dt, state left, state middle, state right) {
        double mom = middle.mom - (dt/dx)*pressure(left, middle, right);
        return state(middle.rho, mom, middle.type); // U[t+2][i] 
    };

    ///////////////////////////////////////////////////////////////////////////

    template <typename F>
    void initial_conditions(F f) {
        for (coord i = 0; i < nx; ++i) {
            state ui = f(i, nx);
            if (i == 0)    ui.type = LEFT_BOUNDARY;
            if (i == nx-1) ui.type = RIGHT_BOUNDARY;
            U[0][i] = U[1][i] = hpx::make_ready_future(ui);
        }

        CFL_dt[0] = hpx::make_ready_future(C*init_dt);
    }

    ///////////////////////////////////////////////////////////////////////////

    // Select a timestep size that enforces the CFL condition.
    double enforce_cfl(double T, std::vector<state> const& u) {
        double min_dt = 10000.0;

        for (coord i = 1; i < nx; ++i) {
            double v = std::fabs(velocity(u[i-1], u[i]));

            if (std::fabs(v-0.0) < 1e-16) // v == 0
                continue;

            double cs = (gamma-1)*gamma; // Speed of sound.

            min_dt = std::fmin(dx/(std::fabs(cs)+v), min_dt);
        }

        // We must divide out the courant number from last_dt. The .get here
        // will not wait; it's a dependency of our dependency.
        double last_dt = CFL_dt[T-1].get()/C; 

        double dt = C*std::fmin(min_dt, last_dt*dt_growth_limiter);

        return dt; 
    }

    ///////////////////////////////////////////////////////////////////////////

    // Retrieve a future referring to the state at time T.
    hpx::future<space> state_solution(coord T) {
        HPX_ASSERT((T >= 0)); HPX_ASSERT(T < nT);
        coord t = T*2-1; // t = substep, T = step; remember T is 1-indexed.
        return hpx::when_all(U[t]);
    } 

    // Retrieve a future referring to the timestep size at time T.
    hpx::shared_future<double> cfl_solution(coord T) {
        HPX_ASSERT((T >= 0)); HPX_ASSERT(T < nT);
        return CFL_dt[T]; 
    } 

    // Dataflow-style solver.
    hpx::future<space> solve() { 
        for (coord t = 2; t < nT*substeps; t += 2) {
            coord T = t/2; // t = substep, T = step

            using hpx::lcos::local::dataflow;   
 
            ///////////////////////////////////////////////////////////////////
            // Compute next timestep size

            // The CFL condition imposes an implicit global barrier; in a code
            // like this, it is the only timestep-wide dependency preventing us
            // from overlapping the computation of multiple timesteps. 
            CFL_dt[T] =
                hpx::when_all(U[t-1])
                    .then(UW2(
                        boost::bind(&solver::enforce_cfl, this, T, _1)
                     ));
    
            ///////////////////////////////////////////////////////////////////
            // Advection step: t-1 -> t
    
            for (coord i = 1; i < nx - 1; ++i)
                U[t][i] =
                    hpx::when_all(CFL_dt[T], U[t-1][i-1], U[t-1][i], U[t-1][i+1])
                        .then(UW2(
                            &advection
                         ));
    
            // Boundary conditions
    
            U[t][0]    = 
                hpx::when_all(CFL_dt[T], ready(LEFT_OUTSIDE), U[t-1][0], U[t-1][1])
                    .then(UW2(
                        &advection
                     ));
    
            U[t][nx-1] =
                hpx::when_all(CFL_dt[T], U[t-1][nx-2], U[t-1][nx-1], ready(RIGHT_OUTSIDE))
                    .then(UW2(
                        &advection
                     ));
    
            ///////////////////////////////////////////////////////////////////
            // Sources Step : t -> t+1
 
            for (coord i = 1; i < nx - 1; ++i)
                U[t+1][i] = dataflow(UW(&sources), CFL_dt[T], U[t][i-1], U[t][i], U[t][i+1]);
//                    hpx::when_all(CFL_dt[T], U[t][i-1], U[t][i], U[t][i+1])
//                        .then(UW(
//                            &sources
//                        )); 

            // Boundary conditions
    
            U[t+1][0]    = dataflow(UW(&sources), CFL_dt[T], ready(LEFT_OUTSIDE), U[t][0], U[t][1]);
//                hpx::when_all(CFL_dt[T], ready(LEFT_OUTSIDE), U[t][0], U[t][1])
//                    .then(UW(
//                        &sources
//                    )); 
    
            U[t+1][nx-1] = dataflow(UW(&sources), CFL_dt[T], U[t][nx-2], U[t][nx-1], ready(RIGHT_OUTSIDE));
//                hpx::when_all(CFL_dt[T], U[t][nx-2], U[t][nx-1], ready(RIGHT_OUTSIDE))
//                    .then(UW(
//                        &sources
//                    )); 
        }

        return state_solution(nT-1);
    }
};


int main()
{
    coord nx = 100;
    coord nT = 400;

    solver s(nx, nT);

    // Gaussian density distribution at rest (e.g. no initial momentum).
    s.initial_conditions(
        [] (coord i, coord max) {
            double const xmid = double(max-1)/2.0;
            double rho = 1.0 + 0.3*std::exp(-std::pow(i-xmid, 2)/std::pow(0.1*(max-1), 2));
            return state(rho, 0.0); 
        }
    ); 

    s.solve();

    std::ofstream U_out("U_hpx.dat");
    std::ofstream dt_out("dt_hpx.dat");

    for (coord T = 1; T < nT; ++T) {
        ///////////////////////////////////////////////////////////////////////
        // Write a block of state records

        solver::space u = s.state_solution(T).get();

        for (coord i = 0; i < u.size(); ++i) {
            state ui = u[i].get(); // Won't block; always ready.
            double v = (i != 0) ? solver::velocity(u[i-1].get(), ui) : 0; 
            U_out << (boost::format("%i %i %.12g %.12g %.12g\n") % T % i % ui.rho % ui.mom % v);
        }

        U_out << "\n";

        ///////////////////////////////////////////////////////////////////////
        // Write a dt record

        dt_out << (boost::format("%i %.12g\n") % T % s.cfl_solution(T).get());
        std::cout << (boost::format("STEP %i DT %.12g\n") % T % s.cfl_solution(T).get());
    }

    return 0;
}

