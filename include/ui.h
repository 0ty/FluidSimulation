#ifndef FLUIDSIM_UI_H
#define FLUIDSIM_UI_H

#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <imgui.h>
#include <imgui-SFML.h>

#include "uidata.h"

class UI {
public:
    UI(sf::RenderWindow& window, UIdata& uidata, int menu_width, int menu_height);

    void drawUIandExtractUIdata();

    void handleEvents();

private:
    sf::RenderWindow& window;
    int menu_width;
    int menu_height;
    UIdata &uidata; // link underlying data to write to or get information from
    sf::Event event{};
    sf::Clock deltaClock{};

    void mainLoopInformation() const;

    void resetButtons();

    void parameterSliders();

    void interactionSelect();

    void totalFieldQuantitiesTimelines() const;

    void seperator() const;

    /// Wrapper for min/max clamped ImGui::SliderFloat
    bool sliderFloatClamp(const char* label, float* v, float v_min, float v_max, const char* format = "%.3f", float power = 1.0f) const;

    void helpMarker(const char* desc) const;

    void totalFieldsHeaderHint() const;


};

#endif //FLUIDSIM_UI_H
