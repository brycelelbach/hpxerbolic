//  Copyright (c) 2014 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define HPX_LIMIT 8

#include <hpx/hpx.hpp>
#include <hpx/hpx_main.hpp>

#define UW hpx::util::unwrapped

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
    typedef hpx::shared_future<double> timestep_size;

    struct solution {
        space U;
        timestep_size dt;
    };

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
    const coord nx; // Number of grid points; 0-indexed.
    space U;        // U[i] is the state of position i

    timestep_size cfl_dt; // The timestep size (dt).

  public:
    // Generate the initial state, using a user-supplied initialization function.
    template <typename IC>
    solver(coord nx_, IC ic) : nx(nx_), U(nx) {
        for (coord i = 0; i < nx; ++i) {
            state ui = ic(i, nx);
            if (i == 0)    ui.type = LEFT_BOUNDARY;
            if (i == nx-1) ui.type = RIGHT_BOUNDARY;
            U[i] = hpx::make_ready_future(ui);
        }

        cfl_dt = hpx::make_ready_future(C*init_dt);
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

    // Select a timestep size that enforces the CFL condition.
    static double enforce_cfl(double last_dt, space const& u) {
        double min_dt = 10000.0;

        for (coord i = 1; i < u.size(); ++i) {
            // These calls to .get() won't block; the futures should be ready.
            double v = std::fabs(velocity(u[i-1].get(), u[i].get()));

            if (std::fabs(v-0.0) < 1e-16) // v == 0
                continue;

            double cs = (gamma-1)*gamma; // Speed of sound.

            min_dt = std::fmin(dx/(std::fabs(cs)+v), min_dt);
        }

        // We must divide out the Courant number from last_dt.
        double dt = C*std::fmin(min_dt, (last_dt/C)*dt_growth_limiter);

        return dt;
    }

    ///////////////////////////////////////////////////////////////////////////

    // Dataflow-style solver.
    solution step() {
        using hpx::lcos::local::dataflow;
        using hpx::util::unwrapped;
        using hpx::util::unwrapped2;

        // Operators
        auto EnforceCFL = unwrapped(&enforce_cfl);
        auto Advection  = unwrapped(&advection);
        auto Sources    = unwrapped(&sources);

        ///////////////////////////////////////////////////////////////////
        // Compute next timestep size

        // The CFL condition imposes an implicit global barrier; in a code
        // like this, it is the only timestep-wide dependency preventing us
        // from overlapping the computation of multiple timesteps.
        cfl_dt = dataflow(EnforceCFL, cfl_dt, hpx::when_all(U));

        ///////////////////////////////////////////////////////////////////
        // Advection step: T -> T+(1/2)

        for (coord i = 1; i < nx - 1; ++i)
            U[i] = dataflow(Advection, cfl_dt, U[i-1], U[i], U[i+1]);

        // Boundary conditions
        U[0]    = dataflow(Advection, cfl_dt, ready(LEFT_OUTSIDE), U[0], U[1]);
        U[nx-1] = dataflow(Advection, cfl_dt, U[nx-2], U[nx-1], ready(RIGHT_OUTSIDE));

        ///////////////////////////////////////////////////////////////////
        // Sources Step : T+(1/2) -> T+1

        for (coord i = 1; i < nx - 1; ++i)
            U[i] = dataflow(Sources, cfl_dt, U[i-1], U[i], U[i+1]);

        // Boundary conditions
        U[0]    = dataflow(Sources, cfl_dt, ready(LEFT_OUTSIDE), U[0], U[1]);
        U[nx-1] = dataflow(Sources, cfl_dt, U[nx-2], U[nx-1], ready(RIGHT_OUTSIDE));

        return solution{U, cfl_dt};
    }
};

int main()
{
    coord nx = 100;
    coord nT = 1000;

    solver sim(nx,
        // Gaussian density distribution at rest (e.g. no initial momentum).
        [] (coord i, coord nx) {
            double xmid = double(nx-1)/2.0;
            double dg = std::pow(0.1*(nx-1), 2);
            double rho = 1.0 + 0.3*std::exp(-std::pow(i-xmid, 2)/dg);
            return state(rho, 0.0);
        }
    );

    // This future is used to serialize I/O; any operation that writes to a file
    // will attach itself to the end of the chain.
    hpx::future<void> io_chain = hpx::make_ready_future();

    std::ofstream U_out("U_hpx.dat");
    std::ofstream dt_out("dt_hpx.dat");

    for (coord T = 1; T < nT; ++T) {
        solver::solution S = sim.step();

        using hpx::lcos::local::dataflow;
        using hpx::util::unwrapped;
        using hpx::util::unwrapped2;

        // Attach asynchronous output operations to the computational tree.
        io_chain = dataflow(unwrapped2(
            [T, &U_out] (std::vector<state> const& u) {
                // Write a state record.
                for (coord i = 0; i < u.size(); ++i) {
                    state ui = u[i];
                    double v = (i != 0) ? solver::velocity(u[i-1], ui) : 0;
                    U_out << ( boost::format("%i %i %.12g %.12g %.12g\n")
                             % T % i % ui.rho % ui.mom % v);
                }

                U_out << "\n";
            }),
            hpx::when_all(S.U),
            std::move(io_chain)
        );

        io_chain = dataflow(unwrapped(
            [T, &dt_out] (double dt) {
                // Write a timestep size record.
                dt_out << (boost::format("%i %.12g\n") % T % dt);

                std::cout << (boost::format("STEP %i DT %.12g\n") % T % dt);

                if ((dt <= 1e-8) || (dt >= 1e8))
                    std::cout << "ERROR: Timestep size is outside of tolerance, "
                                 "numeric instability suspected\n";
            }),
            S.dt,
            std::move(io_chain)
        );
    }

    // Now, we wait.
    io_chain.get();

    return 0;
}

