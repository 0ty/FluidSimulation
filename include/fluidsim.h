#ifndef FLUIDSIM_FLUIDSIM_H
#define FLUIDSIM_FLUIDSIM_H

#include <iostream>
#include <vector>
#include <array>
#include <tuple>

#include <Eigen/Dense>
#include <omp.h>

class FluidSim{
public:

    //////////////////// Typedef: ////////////////////

    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    enum class BC{
        NONE = 0,
        WALL = 1,
        INLET,
        OUTLET
    };
    using BCMatrix = Eigen::Matrix<BC, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    template<int N> using LGFieldPairs = std::array<std::pair<double(*)[3], const Matrix*>,N>; // array of pairs of local and global fields

    //////////////////// Dimensions: ////////////////////

    const int height;
    const int width;
    double dx;

    //////////////////// Fields: ////////////////////

    // primitive variables:
    Matrix rho;
    Matrix u;
    Matrix v;
    Matrix e;
    Matrix p;
    Matrix T;

    Matrix mu; // viscocity (not a constant value but a field here, so it can change within fluid (e.g. with temperature) )
    Matrix k; // thermal conductivity
    double g{0};

    // conserved variables [rho, rho*u, rho*v, E_tot] (contained within U):
    std::array<Matrix, 4> U;

    //////////////////// Fields for Simulation Step Calculation: ////////////////////

    std::array<Matrix, 4> flux_east; // change of conserved variables U through the east side of each cell in one time step
    std::array<Matrix, 4> flux_north;
    std::array<Matrix, 4> U_pre;
    std::array<Matrix, 4> U2;

    //////////////////// Standard Values: ////////////////////

    double std_rho{1.};
    double std_e{1.};
    double std_u{0.};
    double std_v{0.};

    double std_mu{0.3};
    double std_k{20};
    double std_g{0};

    //////////////////// Thermodynamic Parameters: ////////////////////

    double gamma{1.4}; // ratio of specific heats
    double R{287}; // gas constant

    //////////////////// Utility: ////////////////////

private:
    int alternate{0b01};
    int dxf{}; /// use forward differences for the outer x-derivative
    int dyf{}; // use forward differences for the outer y-derivative
    BCMatrix bcMatrix;


public:
    FluidSim(int height, int width, double dx);

    /// Calculate the remaining thermodynamic variables: Given {rho,e}, calculate {p,T}.
    ///
    /// State principle of thermodynamics: Local thermodynamic state is fixed by
    /// any two independent thermodynamic variables {rho,p,e,T,...}.
    /// I choose {rho,e} as the independent variables, so {p,T} can be determined from them.
    void calculate_p_T();

    /// component-wise version of calculate_p_T()
    void calculate_p_T(int i, int j);

    void resetFields();

    void resetParameters();

    void resetBoundaryConditions();

    /// Calculate the primitive variables [rho,u,v,e,p,T] (member variables)
    /// from the provided conserved variables U_in=[rho,rho*i,rho*v,E_tot].
    void conservativeToPrimative(const std::array<Matrix, 4>& U_in);

    /// Calculate the conserved variables U_in=[rho,rho*i,rho*v,E_tot] (passed by reference)
    /// from the primitive variables [rho,u,v,e,p,T] (member variables).
    void primitiveToConservative(std::array<Matrix, 4>& U_in) const;

    static void addBlobToField(Matrix& field, int i0, int j0, double value, int blobRadius, double clampValue=999.);

    void addBounaryCondition(BC BOUNDARYCONDITION,int i0, int j0, int radius);

    /// Calculates the average rho,e,u,v of all cells which are not walls
    /// If all cells are walls (Ncells==0), it returns the standard values.
    [[nodiscard]] std::tuple<double,double,double,double> getCellAverages() const;

    /// Converts the boundary condition matrix into an int-vector
    void extractBCMatrix(std::vector<int>& intBCMatrix) const;

    /// Applies Boundary conditions (to primitive fields only!)
    void applyBCtoCells();

    /// Advances the fields one time step according to the MacCormack scheme.
    FluidSim& step(double dt);

    [[nodiscard]] inline std::tuple<int,int,int,int> getPeriodicIndices(int i,int j) const{

        int i_next = (i+1+ height)%height; // adds +height so it also works with negative numbers
        int j_next = (j+1+ width)%width;
        int i_prev = (i-1+ height)%height;
        int j_prev = (j-1+ width)%width;

        return std::make_tuple(i_next, j_next, i_prev, j_prev);
    }

    inline void global_to_local(const Matrix& field_g, double field_l[3][3], int i, int j, int m, int n) const{
        field_l[m][n] = field_g((i + height)%height,(j + width)%width); //periodic neighbors
    }

    [[nodiscard]] std::tuple<double,double,double,double> getFieldTotals() const;

    /// Using the primitive variables [rho,u,v,e,p,T] (member variables),
    /// calculates the east and north fluxes (member variables),
    /// Boundary conditions are used without modifying the primitive variables.
    ///
    ///
    /// The terms are discretized according to the MacCormack method.
    /// It alternates forward and backward differences for the derivatives in certain terms.
    /// The bools dxf and dyf indicate the use of forward/backward differences:
    /// dxf == true/false -> outer x-derivative is forward/backward
    /// dyf == true/false -> outer y-derivative is forward/backward
    ///
    /// Coordinate system:
    ///
    /// +-----> j,v
    /// |
    /// |
    /// V
    /// i,u
    ///
    template<bool dxf, bool dyf>
    void calculateFluxes(){

#pragma omp parallel for collapse(2)
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                auto [i_next, j_next, i_prev, j_prev] = getPeriodicIndices(i, j);


                /// Copy all global fields into corresponding local fields
                /// (in which boundary conditions can be safely applied without overwriting the global fields)
                ///
                ///        GLOBAL MATRIX:                LOCAL MATRIX
                /// (i-1,j-1) (i-1,  j) (i-1,j+1)      (0,0) (0,1) (0,2)
                /// (  i,j-1) (  i,  j) (  i,j+1)  =>  (1,0) (1,1) (1,2)
                /// (i+1,j-1) (i-1,  j) (i+1,j+1)      (2,0) (2,0) (2,2)
                ///
                /// using C-arrays here instead of std::arrays because it is twice as fast.
                /// (because 2D C-arrays have different memory layout than 2D std::arrays: It's guaranteed to be contiguous).

                // local fields:
                double rho_l[3][3];
                double u_l[3][3];
                double v_l[3][3];
                double e_l[3][3];
                double p_l[3][3];
                double T_l[3][3];
                double mu_l[3][3];
                double k_l[3][3];


                LGFieldPairs<8> local_global_field_pairs{
                        std::pair{rho_l, &rho},
                        std::pair{u_l, &u},
                        std::pair{v_l, &v},
                        std::pair{e_l, &e},
                        std::pair{p_l, &p},
                        std::pair{T_l, &T},
                        std::pair{mu_l, &mu},
                        std::pair{k_l, &k}
                };

                // Copy fields into corresponding local fields
                for(auto field_pair : local_global_field_pairs){

                    auto [local_field, global_field] = field_pair;

                    for(int m=0; m<=2; ++m) {
                        for (int n =0; n<=2; ++n) {
                            global_to_local(*global_field, local_field, i-1+m, j-1+n, m, n);
                        }
                    }

                }


                //////////////////// EAST FLUXES: ////////////////////

                // Apply Boundary conditions to local fields:

                //  .  .  .
                //  .  x  .
                //  .  BC .
                switch (bcMatrix(i_next, j)) {
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i, j, 2, 1);
                        }
                        break;
                }

                //  .  BC  .
                //  .  x  .
                //  .  .  .
                switch (bcMatrix(i_prev, j)) {
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i, j, 0, 1);
                        }
                        break;
                }

                //  .  .  .
                //  .  BC .
                //  .  .  .
                switch(bcMatrix(i, j)){
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i_next, j, 1, 1);
                        }
                        break;
                }

                ////////// Navier-Stokes I //////////

                double massflux_east;
                if constexpr(dxf) massflux_east = -rho_l[2][1] * u_l[2][1]; // euler forward diff
                else              massflux_east = -rho_l[1][1] * u_l[1][1]; // euler backward diff

                ////////// Navier-Stokes II //////////

                double conv_u_east;
                if constexpr(dxf) conv_u_east = u_l[2][1] * massflux_east; // negative sign is in mass flux
                else              conv_u_east = u_l[1][1] * massflux_east;

                double pressure_x_east;
                if constexpr(dxf) pressure_x_east = -p_l[2][1];
                else              pressure_x_east = -p_l[1][1];

                double tau_xx_east;
                if constexpr(dxf) tau_xx_east = 2. / dx * mu_l[2][1] * (u_l[2][1] - u_l[1][1] );
                else              tau_xx_east = 2. / dx * mu_l[1][1] * (u_l[2][1] - u_l[1][1] );

                ////////// Navier-Stokes III //////////

                double conv_v_east;
                if constexpr(dxf) conv_v_east = v_l[2][1] * massflux_east;
                else              conv_v_east = v_l[1][1] * massflux_east;

                double tau_xy_vx_east;
                if constexpr(dxf) tau_xy_vx_east = 1. / dx * mu_l[2][1] * (v_l[2][1] - v_l[1][1] );
                else              tau_xy_vx_east = 1. / dx * mu_l[1][1] * (v_l[2][1] - v_l[1][1] );
                double tau_xy_uy_east;
                if constexpr(dxf) tau_xy_uy_east = 1. / dx * mu_l[2][1] * (u_l[2][2] - u_l[2][0] );
                else              tau_xy_uy_east = 1. / dx * mu_l[1][1] * (u_l[1][2] - u_l[1][0] );
                double tau_xy_east = tau_xy_vx_east;// + tau_xy_uy_east;

                ////////// Navier-Stokes IV //////////

                double conv_E_east;
                if constexpr(dxf) conv_E_east = massflux_east * (e_l[2][1] + (u_l[2][1] * u_l[2][1] + v_l[2][1] * v_l[2][1]) / 2. );
                else              conv_E_east = massflux_east * (e_l[1][1] + (u_l[1][1] * u_l[1][1] + v_l[1][1] * v_l[1][1]) / 2. );

                double pressure_E_east;
                if constexpr(dxf) pressure_E_east = u_l[2][1] * pressure_x_east;
                else              pressure_E_east = u_l[1][1] * pressure_x_east;

                double tau_xx_E_east;
                if constexpr(dxf) tau_xx_E_east = u_l[2][1] * tau_xx_east;
                else              tau_xx_E_east = u_l[1][1] * tau_xx_east;

                double tau_xy_E_east;
                if constexpr(dxf) tau_xy_E_east = v_l[2][1] * tau_xy_east;
                else              tau_xy_E_east = v_l[1][1] * tau_xy_east;

                double T_diffus_east;
                if constexpr(dxf) T_diffus_east = k_l[2][1] * (T_l[2][1] - T_l[1][1]);
                else              T_diffus_east = k_l[1][1] * (T_l[2][1] - T_l[1][1]);



                //////////////////// NORTH FLUXES: ////////////////////

                // Apply Boundary conditions to local fields:

                //  .  .  .
                //  .  x  BC
                //  .  .  .
                switch (bcMatrix(i, j_next)) {
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i, j, 1, 2);
                        }
                        break;
                }

                //  .  .  .
                //  BC x  .
                //  .  .  .
                switch (bcMatrix(i, j_prev)) {
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i, j, 1, 0);
                        }
                        break;
                }

                //  .  .  .
                //  .  BC .
                //  .  .  .
                switch(bcMatrix(i, j)){
                    case BC::NONE:
                        break;
                    case BC::WALL:
                        for (auto field_pair: local_global_field_pairs) {
                            auto [local_field, global_field] = field_pair;
                            global_to_local(*global_field, local_field, i, j_next, 1, 1);
                        }
                        break;
                }

                ////////// Navier-Stokes I //////////

                double massflux_north;
                if constexpr(dyf) massflux_north = -rho_l[1][2] * v_l[1][2];
                else              massflux_north = -rho_l[1][1] * v_l[1][1];

                ////////// Navier-Stokes II //////////

                double conv_u_north;
                if constexpr(dyf) conv_u_north = u_l[1][2] * massflux_north;
                else              conv_u_north = u_l[1][1] * massflux_north;

                double tau_xy_uy_north;
                if constexpr(dyf) tau_xy_uy_north = 1. / dx * mu_l[1][2] * (u_l[1][2] - u_l[1][1] );
                else              tau_xy_uy_north = 1. / dx * mu_l[1][1] * (u_l[1][2] - u_l[1][1] );
                double tau_xy_vx_north;
                if constexpr(dyf) tau_xy_vx_north = 1. / dx * mu_l[1][2] * (v_l[2][2] - v_l[0][2] );
                else              tau_xy_vx_north = 1. / dx * mu_l[1][1] * (v_l[2][1] - v_l[0][1] );
                double tau_xy_north = tau_xy_uy_north;

                ////////// Navier-Stokes III //////////

                double conv_v_north;
                if constexpr(dyf) conv_v_north = v_l[1][2] * massflux_north;
                else              conv_v_north = v_l[1][1] * massflux_north;

                double pressure_y_north;
                if constexpr(dyf) pressure_y_north = -p_l[1][2];
                else              pressure_y_north = -p_l[1][1];

                double tau_yy_north;
                if constexpr(dyf) tau_yy_north = 2. / dx * mu_l[1][2] * (v_l[1][2] - v_l[1][1] );
                else              tau_yy_north = 2. / dx * mu_l[1][1] * (v_l[1][2] - v_l[1][1] );

                ////////// Navier-Stokes IV //////////

                double conv_E_north;
                if constexpr(dyf) conv_E_north = massflux_north * (e_l[1][2] + (u_l[1][2] * u_l[1][2] + v_l[1][2] * v_l[1][2]) / 2. );
                else              conv_E_north = massflux_north * (e_l[1][1] + (u_l[1][1] * u_l[1][1] + v_l[1][1] * v_l[1][1]) / 2. );

                double pressure_E_north;
                if constexpr(dyf) pressure_E_north = v_l[1][2] * pressure_y_north;
                else              pressure_E_north = v_l[1][1] * pressure_y_north;

                double tau_yy_E_north;
                if constexpr(dyf) tau_yy_E_north = v_l[1][2] * tau_yy_north;
                else              tau_yy_E_north = v_l[1][1] * tau_yy_north;

                double tau_xy_E_north;
                if constexpr(dyf) tau_xy_E_north = u_l[1][2] * tau_xy_north;
                else              tau_xy_E_north = u_l[1][1] * tau_xy_north;

                double T_diffus_north;
                if constexpr(dyf) T_diffus_north = k_l[1][2] * (T_l[1][2] - T_l[1][1]);
                else              T_diffus_north = k_l[1][1] * (T_l[1][2] - T_l[1][1]);


                //////////////////// FINAL FLUXES: ////////////////////

                flux_east[0](i,j) = massflux_east;
                flux_north[0](i,j) = massflux_north;

                flux_east[1](i,j) = conv_u_east + pressure_x_east + tau_xx_east;
                flux_north[1](i,j) = conv_u_north + tau_xy_north;

                flux_east[2](i,j) = conv_v_east + tau_xy_east;
                flux_north[2](i,j) = conv_v_north + pressure_y_north + tau_yy_north;

                flux_east[3](i,j) = conv_E_east + pressure_E_east + tau_xx_E_east + tau_xy_E_east + T_diffus_east;
                flux_north[3](i,j) = conv_E_north + pressure_E_north + tau_yy_E_north + tau_xy_E_north + T_diffus_north;


            }
        }
    }

};



#endif //FLUIDSIM_FLUIDSIM_H
