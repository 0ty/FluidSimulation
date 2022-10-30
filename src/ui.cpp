#include "ui.h"



UI::UI(sf::RenderWindow &window, UIdata &uidata, int menu_width, int menu_height)
        : window{window},
          uidata{uidata},
          menu_width{menu_width},
          menu_height{menu_height}
{
    ImGui::SFML::Init(window);
}

void UI::drawUIandExtractUIdata() {

    ImGui::SFML::Update(window, deltaClock.restart());

    ImGuiWindowFlags window_flags{};
    window_flags |= ImGuiWindowFlags_NoTitleBar;
    window_flags |= ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoResize;
    window_flags |= ImGuiWindowFlags_NoCollapse;
    ImGui::GetStyle().WindowRounding = 0.0f;
    ImGui::GetStyle().WindowBorderSize = 0.0f;

    ImGui::SetNextWindowPos(ImVec2(0,0));
    ImGui::SetNextWindowSize(ImVec2(menu_width, menu_height));
    ImGui::Begin("Menu", nullptr, window_flags);


    mainLoopInformation();
    seperator();

    resetButtons();
    seperator();

    interactionSelect();
    seperator();

    parameterSliders();
    seperator();

    totalFieldQuantitiesTimelines();
    seperator();



    ImGui::End();
    ImGui::SFML::Render(window);

}

void UI::handleEvents() {

    while (window.pollEvent(event)) {
        ///// ImGui events: /////
        ImGui::SFML::ProcessEvent(event);

        ///// SFML events: /////
        if (event.type == sf::Event::Closed) {
            window.close();
        }
        if (event.type == sf::Event::Resized){
            sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);
            window.setView(sf::View(visibleArea));
        }
        if (event.type == sf::Event::MouseButtonPressed){
            if(event.mouseButton.button == sf::Mouse::Left) uidata.mouseLeftPressed = true;
            if(event.mouseButton.button == sf::Mouse::Right) uidata.mouseRightPressed = true;
        }
        if (event.type == sf::Event::MouseButtonReleased){
            if(event.mouseButton.button == sf::Mouse::Left) uidata.mouseLeftPressed = false;
            if(event.mouseButton.button == sf::Mouse::Right) uidata.mouseRightPressed = false;
        }
    }

}

void UI::mainLoopInformation() const {

    // convert fpsCycletime into seconds with decimal places (i.e. underlying type is now double instead of int)
    double fpsCycletimeSeconds = (std::chrono::duration<double, std::ratio<1, 1>>{ uidata.fpsCycletime }).count();

    ImGui::Text( std::format("FPS: {:.1f}", 1./fpsCycletimeSeconds).c_str() );

    ImGui::Text( std::format("Simulation steps per frame: {}", uidata.simStepCount).c_str() );

}

void UI::resetButtons() {

    ImGui::BeginGroup();

    ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(3.f, 0.8f, 0.6f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(3.f, 1.0f, 0.8f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(3.f, 0.7f, 0.5f));
    if(ImGui::Button("Reset\nAll")){
        uidata.resetFields = true;
        uidata.resetBoundaries = true;
        uidata.resetParameters = true;
    }
    ImGui::PopStyleColor(3);

    ImGui::SameLine();
    if(ImGui::Button("Reset\nFields")) uidata.resetFields = true;
    ImGui::SameLine();
    if(ImGui::Button("Reset\nBoundaries")) uidata.resetBoundaries = true;
    ImGui::SameLine();
    if(ImGui::Button("Reset\nParameters")) uidata.resetParameters = true;

    ImGui::EndGroup();

}

void UI::parameterSliders() {

    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.55f);

    if (sliderFloatClamp("Viscocity", &(uidata.input_mu), 0.1, 0.8, "mu = %.2f"))
        uidata.parameterChanged = true;
    if (sliderFloatClamp("Thermal Conduct.", &(uidata.input_k), 20, 500, "k = %.0f"))
        uidata.parameterChanged = true;
    if (sliderFloatClamp("Gravity", &(uidata.input_g), 0, 0.0002, "g = %.5f"))
        uidata.parameterChanged = true;

    ImGui::PopItemWidth();

}

void UI::interactionSelect() {

    if (ImGui::Selectable("Move Fluid", uidata.interactionMode == UIdata::INTERACTION_MODE::MOVE_FLUID))
        uidata.interactionMode = UIdata::INTERACTION_MODE::MOVE_FLUID;

    if (ImGui::Selectable("Add/Remove Wall", uidata.interactionMode == UIdata::INTERACTION_MODE::EDIT_WALL))
        uidata.interactionMode = UIdata::INTERACTION_MODE::EDIT_WALL;

    ImGui::SameLine();
    helpMarker("Left click to add wall; Right click to remove.");

    if (uidata.interactionMode != UIdata::INTERACTION_MODE::MOVE_FLUID){
        sliderFloatClamp("Brush Size", &(uidata.brushSize), 1, 6, "%.0f");
    }

}

void UI::totalFieldQuantitiesTimelines() const {

    if (ImGui::CollapsingHeader("Total Field Quantities")){

        totalFieldsHeaderHint(); // show in beginning when header is not collaped

        uidata.calculateFieldTotals = true;

        ImVec2 graphSize{ImGui::GetWindowWidth() * 0.45f, 60.0f};

        ImGui::PlotLines("", uidata.rho_timeline, UIdata::timelineLength, uidata.timelineIdx,
                         std::format("Total Density\nrho = {:.2f}", uidata.rho_timeline[uidata.timelineIdx]).c_str(),
                         FLT_MAX, FLT_MAX, ImVec2(graphSize));
        ImGui::SameLine();
        ImGui::PlotLines("", uidata.E_timeline, UIdata::timelineLength, uidata.timelineIdx,
                         std::format("Total Energy\nE = {:.2f}", uidata.E_timeline[uidata.timelineIdx]).c_str(),
                         FLT_MAX, FLT_MAX, ImVec2(graphSize));


        ImGui::PlotLines("", uidata.rhou_timeline, UIdata::timelineLength, uidata.timelineIdx,
                         std::format("Total x-Momentum\nrho*u = {:.2f}", uidata.rhou_timeline[uidata.timelineIdx]).c_str(),
                         FLT_MAX, FLT_MAX, ImVec2(graphSize));
        ImGui::SameLine();
        ImGui::PlotLines("", uidata.rhov_timeline, UIdata::timelineLength, uidata.timelineIdx,
                         std::format("Total y-Momentum\nrho*v = {:.2f}", uidata.rhov_timeline[uidata.timelineIdx]).c_str(),
                         FLT_MAX, FLT_MAX, ImVec2(graphSize));

    } else{
        uidata.calculateFieldTotals = false;
        totalFieldsHeaderHint(); // show when header is collapsed
    }

}

void UI::seperator() const {
    ImGui::Dummy(ImVec2(0.0f, 5.0f));
    ImGui::Separator();
    ImGui::Dummy(ImVec2(0.0f, 5.0f));
}

bool UI::sliderFloatClamp(const char *label, float *v, float v_min, float v_max, const char *format, float power) const {
    float v_backup = *v;
    if (!ImGui::SliderFloat(label, v, v_min, v_max, format, power))
        return false;
    if (*v < v_min) *v = v_min;
    else if (*v > v_max) *v = v_max;
    return v_backup != *v;
}

void UI::helpMarker(const char *desc) const {
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void UI::totalFieldsHeaderHint() const {
    ImGui::SameLine();
    helpMarker("Calculation of total field quantities can be expensive. Collapse this header to improve performance.");
}


