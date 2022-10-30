#ifndef FLUIDSIM_UIDATA_H
#define FLUIDSIM_UIDATA_H

#include <chrono>
#include "fluidsim.h"
#include "timetypes.h"


/// Architecture of the user interface:
///
///    ┌──────────┐         ┌───────────────┐               ┌─────────────┐
///    │ Fluidsim │         │    UIdata     │               │     UI      │
///    │          │         │               │               │             │
///    │          ◄─────────applyUIdata     │               │             │
///    │          ┼────────►extractSimData  │               │             │
///    │          │         │               ◄─────►drawUIandExtractUIdata │
///    └──────────┘         └───────────────┘               └─────────────┘
///                                                        (can be replaced)
///
///   Fluidsim does not know about UIdata nor UI. It is the standalone simulation.
///
///   UIdata contains all the settings fed by the UI (by the user), and can apply them to the fluidsim.
///   It also extracts data from the fluidsim that can later be displayed in the UI.
///   UIdata does not know about the UI. Thus, the UI can be replaced with a different UI at anytime.
///
///   UI links to the UIdata (member variable reference) and feeds it with user input data (e.g. the value of the parameter sliders),
///   as well as extracts data to finally display it (e.g. extract total energy to plot it).
///
class UIdata {
public:
    ////////// Reset Buttons: //////////

    bool resetFields{false};
    bool resetBoundaries{false};
    bool resetParameters{false};

    ////////// Parameters: //////////

    bool parameterChanged{false};
    float input_mu{0.3};
    float input_k{20};
    float input_g{0.0};

    ////////// Timelines: //////////

    bool calculateFieldTotals{false};
    timetypes::tp lastTimelineUpdate{};
    std::chrono::milliseconds timelineUpdateInterval{20};
    static constexpr int timelineLength{100};
    int timelineIdx{0};
    float rho_timeline[timelineLength]{};
    float rhou_timeline[timelineLength]{};
    float rhov_timeline[timelineLength]{};
    float E_timeline[timelineLength]{};

    ////////// Selected Interaction Mode: //////////

    enum class INTERACTION_MODE{
        MOVE_FLUID,
        EDIT_WALL
    };
    INTERACTION_MODE interactionMode{INTERACTION_MODE::MOVE_FLUID};
    float brushSize{1};

    ////////// Mouse Status: //////////

    bool mouseLeftPressed{false};
    bool mouseRightPressed{false};

    ////////// Main Render Loop Information: //////////

    timetypes::ns fpsCycletime{};
    int simStepCount{};

public:

    void applyUIdata(FluidSim& fluidsim);

    void extractSimData(const FluidSim& fluidsim);

};


#endif //FLUIDSIM_UIDATA_H
