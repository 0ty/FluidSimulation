#include "uidata.h"


void UIdata::applyUIdata(FluidSim &fluidsim) {

    if(resetParameters){
        resetParameters = false;
        input_mu = fluidsim.std_mu;
        input_k = fluidsim.std_k;
        input_g = fluidsim.std_g;
        parameterChanged = true;
    }

    if(resetBoundaries){
        resetBoundaries = false;
        fluidsim.resetBoundaryConditions();
    }

    if(resetFields){
        resetFields = false;
        fluidsim.resetFields();
    }

    if(parameterChanged){
        parameterChanged = false;
        fluidsim.mu = FluidSim::Matrix::Constant(fluidsim.height,fluidsim.width, input_mu);
        fluidsim.k = FluidSim::Matrix::Constant(fluidsim.height,fluidsim.width, input_k);
        fluidsim.g = input_g;
    }

}

void UIdata::extractSimData(const FluidSim &fluidsim) {

    if(calculateFieldTotals && timetypes::sclock::now()-lastTimelineUpdate > timelineUpdateInterval) {
        timelineIdx = (timelineIdx - 1 + timelineLength) % timelineLength;

        auto [rho_total, rhou_total, rhov_total, E_total] = fluidsim.getFieldTotals();
        rho_timeline[timelineIdx] = rho_total;
        rhou_timeline[timelineIdx] = rhou_total;
        rhov_timeline[timelineIdx] = rhov_total;
        E_timeline[timelineIdx] = E_total;

        lastTimelineUpdate = timetypes::sclock::now();

    }
}
