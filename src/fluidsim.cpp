#include "fluidsim.h"

FluidSim::FluidSim(int height, int width, double dx)
        : height{height},
          width{width},
          dx{dx},
          U{Matrix(height,width), Matrix(height,width), Matrix(height,width), Matrix(height,width)},
          flux_east{Matrix(height,width), Matrix(height,width), Matrix(height,width), Matrix(height,width)},
          flux_north{Matrix(height,width), Matrix(height,width), Matrix(height,width), Matrix(height,width)},
          U_pre{Matrix(height,width), Matrix(height,width), Matrix(height,width), Matrix(height,width)},
          U2{Matrix(height,width), Matrix(height,width), Matrix(height,width), Matrix(height,width)}
{
    resetFields();
    resetParameters();
    bcMatrix = BCMatrix::Constant(height, width, BC::NONE);

}

void FluidSim::calculate_p_T() {
    p = e.cwiseProduct(rho) * (gamma - 1.);
    T = (gamma - 1) * e / R;
}

void FluidSim::calculate_p_T(int i, int j) {
    p(i,j) = e(i,j)*rho(i,j) * (gamma - 1.);
    T(i,j) = (gamma - 1) * e(i,j) / R;
}

void FluidSim::resetFields() {

    rho = FluidSim::Matrix::Constant(height,width, std_rho);
    e = FluidSim::Matrix::Constant(height,width, std_e);
    u = FluidSim::Matrix::Constant(height,width, std_u);
    v = FluidSim::Matrix::Constant(height,width, std_v);
    calculate_p_T();

    //conservativeToPrimative()
}

void FluidSim::resetParameters() {

    mu = FluidSim::Matrix::Constant(height,width, std_mu);
    k = FluidSim::Matrix::Constant(height,width, std_k);
    g = std_g;

}

void FluidSim::resetBoundaryConditions() {

    auto [rho_avg, e_avg, u_avg, v_avg] = getCellAverages();

#pragma omp parallel for collapse(2)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            if (bcMatrix(i, j) == BC::WALL) {
                rho(i,j) = rho_avg;
                e(i,j) = e_avg;
                u(i,j) = u_avg;
                v(i,j) = v_avg;
                calculate_p_T(i,j);
            }

        }
    }

    bcMatrix.setConstant(BC::NONE);
}

void FluidSim::conservativeToPrimative(const std::array<Matrix, 4> &U_in) {

    rho = U_in[0];
    u = U_in[1].cwiseQuotient(U_in[0]);
    v = U_in[2].cwiseQuotient(U_in[0]);
    e = U_in[3].cwiseQuotient(rho) - (u.cwiseProduct(u) + v.cwiseProduct(v) ) / 2.;
    p = (e.cwiseProduct(rho))*(gamma-1);
    T = (gamma-1)*e/R;
}

void FluidSim::primitiveToConservative(std::array<Matrix, 4> &U_in) const {

    U_in[0] = rho;
    U_in[1] = rho.cwiseProduct(u);
    U_in[2] = rho.cwiseProduct(v);
    U_in[3] = (e + (u.cwiseProduct(u) + v.cwiseProduct(v)) / 2).cwiseProduct(rho);
}

void FluidSim::addBlobToField(FluidSim::Matrix &field, int i0, int j0, double value, int blobRadius, double clampValue) {

    int f_height{ static_cast<int>(field.rows()) };
    int f_width{ static_cast<int>(field.cols()) };

//    if(blobRadius>f_height || blobRadius>f_width){
//        throw std::exception("Radius of blob cannot be larger than the matrix.");
//    }

    for (int i = -blobRadius; i <= blobRadius; ++i) {
        for (int j = -blobRadius; j <= blobRadius; ++j) {

            // modulo:
            int i_temp = (i0 + i + f_height) % f_height;
            int j_temp = (j0 + j + f_width) % f_width;

            double strength = 1-std::sqrt(i*i+j*j)/blobRadius;
            if(strength<0) strength = 0; // strength = between [0,1] depending on how far away from middle (i0,j0).

            field(i_temp, j_temp) += value * strength;

            if(field(i_temp, j_temp) > clampValue) field(i_temp, j_temp) = clampValue;
            if(field(i_temp, j_temp) < -clampValue) field(i_temp, j_temp) = -clampValue;

        }
    }

}

void FluidSim::addBounaryCondition(FluidSim::BC BOUNDARYCONDITION, int i0, int j0, int radius) {

//        if(radius>height || radius>width){
//            throw std::exception("Radius cannot be larger than the matrix.");
//        }

    auto [rho_avg, e_avg, u_avg, v_avg] = getCellAverages();

    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {

            int i_temp = (i0 + i + height) % height;
            int j_temp = (j0 + j + width) % width;

            // First, restore field values when the boundary condition of the cell was a wall,
            // because values on wall cells are garbage.
            if(bcMatrix(i_temp, j_temp) == BC::WALL){
                rho(i_temp, j_temp) = rho_avg;
                e(i_temp, j_temp) = e_avg;
                u(i_temp, j_temp) = u_avg;
                v(i_temp, j_temp) = v_avg;
                calculate_p_T(i_temp, j_temp);
            }

            if(std::sqrt(i*i+j*j) < radius) {
                bcMatrix(i_temp, j_temp) = BOUNDARYCONDITION;
            }

        }
    }

}

std::tuple<double, double, double, double> FluidSim::getCellAverages() const {

    double rho_avg{0};
    double e_avg{0};
    double u_avg{0};
    double v_avg{0};

    int Ncells{0}; // Number of cells that are not walls

#pragma omp parallel for collapse(2) reduction(+:rho_avg, e_avg, u_avg, v_avg, Ncells)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            if (bcMatrix(i, j) != BC::WALL) {
                rho_avg += rho(i,j);
                e_avg += e(i,j);
                u_avg += u(i,j);
                v_avg += v(i,j);
                ++Ncells;
            }

        }
    }

    rho_avg = Ncells !=0 ? rho_avg/Ncells : std_rho;
    e_avg = Ncells !=0 ? e_avg/Ncells : std_e;
    u_avg = Ncells !=0 ? u_avg/Ncells : std_u;
    v_avg = Ncells !=0 ? v_avg/Ncells : std_v;

    return std::make_tuple(rho_avg, e_avg, u_avg, v_avg);

}

void FluidSim::extractBCMatrix(std::vector<int> &intBCMatrix) const {

//    if( intBCMatrix.size() != bcMatrix.size() ){
//        throw std::exception("Provided int-vector must be same size (i.e. height*width) as fluidsim BC-Matrix");
//    }


#pragma omp parallel for collapse(2)
    for(int i=0; i<height; ++i) {
        for (int j = 0; j < width; ++j) {

            int idx = (i*width + j);

            intBCMatrix[idx] = static_cast<int>(bcMatrix(i,j));
        }
    }

}

void FluidSim::applyBCtoCells() {

#pragma omp parallel for collapse(2)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            auto [i_next, j_next, i_prev, j_prev] = getPeriodicIndices(i, j);

            if (bcMatrix(i, j) == BC::WALL) {
                rho(i,j) = std_rho;
                e(i,j) = std_e;
                u(i,j) = std_u;
                v(i,j) = std_v;
                calculate_p_T(i, j);
            }

            if (bcMatrix(i_prev, j) == BC::WALL) {
                u(i, j) = 0;
            }

            if (bcMatrix(i_next, j) == BC::WALL) {
                u(i, j) = 0;
            }

            if (bcMatrix(i, j_next) == BC::WALL) {
                v(i, j) = 0;
            }

            if (bcMatrix(i, j_prev) == BC::WALL) {
                v(i, j) = 0;
            }
        }
    }

}

std::tuple<double, double, double, double> FluidSim::getFieldTotals() const {

    double rho_total{0.};
    double rhou_total{0.};
    double rhov_total{0.};
    double E_total{0.};


#pragma omp parallel for collapse(2) reduction(+:rho_total, E_total, rhou_total, rhov_total)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if(bcMatrix(i, j) != BC::WALL){
                rho_total += U[0](i,j); // missing *dx*dy for actual integral, but ignore here since cells all have same size
                rhou_total += U[1](i,j);
                rhov_total += U[2](i,j);
                E_total += U[3](i,j);
            }
        }
    }

    return std::make_tuple(rho_total, rhou_total, rhov_total, E_total);
}

FluidSim &FluidSim::step(double dt) {

    applyBCtoCells();

    // The MacCormack method alternates forward and backward differences.
    // For this, the bits of the utility variable "alternate" are used.
    // (alternate increments every loop with modulo 4: 0,1,2,3,0,1,2,3,...).
    dxf = ((alternate&0b01)>>0); // 1,0,1,0,1,0,1,0,... (alternates every loop)
    dyf = ((alternate&0b10)>>1); // 0,1,1,0,1,1,0,0,... (alternates every 2nd loop)

    /////////////// Predictor: ///////////////

    primitiveToConservative(U);

    if(!dxf && !dyf)      calculateFluxes<false,false>();
    else if (!dxf && dyf) calculateFluxes<false,true>();
    else if (dxf && !dyf) calculateFluxes<true,false>();
    else if (dxf && dyf)  calculateFluxes<true,true>();

    // Apply fluxes (changes through surfaces of cells)
    for (int l = 0; l < 4; ++l) {
#pragma omp parallel for collapse(2)
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {

                auto [i_next, j_next, i_prev, j_prev] = getPeriodicIndices(i, j);

                U_pre[l](i, j) = U[l](i, j) + dt/dx*
                                              ( flux_east[l](i,j) - flux_east[l](i_prev,j)
                                                + flux_north[l](i,j) - flux_north[l](i,j_prev) );
            }
        }
    }

    // Apply body forces
#pragma omp parallel for collapse(2)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if(bcMatrix(i,j) == BC::NONE){
                U_pre[1](i,j) += rho(i,j)*g;
                U_pre[3](i,j) += rho(i,j)*g*u(i,j);
            }
        }
    }

    conservativeToPrimative(U_pre);
    applyBCtoCells();
    primitiveToConservative(U_pre);



    /////////////// Corrector: ///////////////

    // In the corrector the derivative approximation has to be opposite to the ones in the predictor.
    // E.g. if the predictor used forward differences for the x-derivative, it must now use backward differences.
    if(!dxf && !dyf)      calculateFluxes<true,true>();
    else if (!dxf && dyf) calculateFluxes<true,false>();
    else if (dxf && !dyf) calculateFluxes<false,true>();
    else if (dxf && dyf)  calculateFluxes<false,false>();

    // Apply fluxes (changes through surface of cells)
    for(int i=0; i<height; ++i) {
#pragma omp parallel for collapse(2)
        for (int j = 0; j < width; ++j) {
            for (int l = 0; l < 4; ++l) {

                auto [i_next, j_next, i_prev, j_prev] = getPeriodicIndices(i, j);

                U2[l](i, j) = 0.5 * ( U[l](i, j) + U_pre[l](i, j) + dt/dx*
                                                                    ( flux_east[l](i,j) - flux_east[l](i_prev,j)
                                                                      + flux_north[l](i,j) - flux_north[l](i,j_prev) )
                );

            }
        }
    }

    // Apply Body forces
#pragma omp parallel for collapse(2)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if(bcMatrix(i,j) == BC::NONE) {
                U2[1](i, j) += rho(i, j) * g;
                U2[3](i, j) += rho(i, j) * g * u(i, j);
            }
        }
    }

    U = U2;
    conservativeToPrimative(U);
    applyBCtoCells();
    primitiveToConservative(U);


    alternate = (alternate+1)%4;

    return *this;
}
