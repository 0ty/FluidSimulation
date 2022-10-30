#include <iostream>
#include <chrono>
#include <cmath>

#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "timetypes.h"
#include "fluidsim.h"
#include "ui.h"
#include "uidata.h"
#include "fieldsvis.h"


int main(int argc, char *argv[]) {
    // Optional command line arguments: height(int), width(int), scale(double)

    int height{60};
    int width{80};

    double scale{4.};

    if(argc >= 2) height = std::stoi(argv[1]);
    if(argc >= 3) width = std::stoi(argv[2]);
    if(argc >= 4) scale = std::stod(argv[3]);


    //////////////////////////////////////////////////


    int menu_width{300};
    sf::RenderWindow window{sf::VideoMode{static_cast<unsigned int>(menu_width + 3 * scale * width),
                                          static_cast<unsigned int>(2 * scale * height)}, "Fluid Simulation"};

    UIdata uidata{};
    UI ui{window, uidata, menu_width, static_cast<int>(height*2*scale) };


    double dx{1.};
    double dt{0.05};
    FluidSim fluidsim{height, width, dx};
    FluidSim::addBlobToField(fluidsim.rho,width*0.4,height*0.2,0.05,std::min(width,height)/2 );

    FieldsVis fieldsVis{fluidsim, scale};
    fieldsVis.arrangeSprites({static_cast<float>(menu_width),0});


    //////////////////////////////////////////////////


    // "render_dt" = real physical time a simulation step (with the unitless time dt) should take:
    timetypes::ns render_dt{ std::chrono::milliseconds{2} };

    // to calculate FPS of every loop (starts at zero every loop):
    timetypes::tp fpsTimer{};

    // "elapsedTime" = accumulates remaining simulation time.
    // Indicates how much the simulation lags behind the render / how much time the simulation needs to progress
    // to catch up with the render.
    // (It does not start at zero every loop; it might have some remaining simulation time left
    // from previous loop which has to be processed)
    timetypes::ns elapsedTime{0};
    timetypes::tp t2{};
    timetypes::tp t1{ timetypes::sclock::now()+std::chrono::milliseconds{200} }; // give 200ms of head start for first loop, since it's usually a bit slower, to avoid stall.


    // Main render loop:
    while (window.isOpen())
    {
        ui.handleEvents();

        //////////////////// Rendering: ////////////////////

        // Calculate time between renders:
        uidata.fpsCycletime = timetypes::sclock::now() - fpsTimer;
        fpsTimer = timetypes::sclock::now(); // reset timer


        window.clear(sf::Color{ 62,71,81 });
        fieldsVis.drawFields(window);
        fieldsVis.drawBrushCircle(uidata, window);
        ui.drawUIandExtractUIdata();
        window.display();

        // Calculate accumulated time that simulation has to process to catch up again:
        t2 = timetypes::sclock::now();
        elapsedTime += t2 - t1;
        t1 = timetypes::sclock::now();


        //////////////////// Simulation: ////////////////////

        uidata.applyUIdata(fluidsim);
        uidata.simStepCount = 0; // reset simStepCount

        // Simulation while-loop: Do as many simulation steps as necessary to catch up to render:
        // Reduce "elapsedTime" (how much simulation lags behind render) with every simulation step until the simulation has caught up.
        while (elapsedTime >= render_dt) {
            elapsedTime -= render_dt;
            ++uidata.simStepCount;

            fieldsVis.applyInteraction(uidata,window);

            fluidsim.step(dt);

            uidata.extractSimData(fluidsim);

        }
    }

}
