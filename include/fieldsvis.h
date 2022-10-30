#ifndef FLUIDSIM_FIELDSVIS_H
#define FLUIDSIM_FIELDSVIS_H

#include <array>
#include <chrono>

#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "timetypes.h"
#include "fluidsim.h"
#include "matrixtexture.h"
#include "uidata.h"



class FieldsVis {

public:
    explicit FieldsVis(FluidSim& fluidsim, double scale=1.);

    void arrangeSprites(const sf::Vector2f& offset = {0,0}, int spritesPerRow = 3);

    void drawFields(sf::RenderWindow& window);

    void drawBrushCircle(const UIdata& uidata, sf::RenderWindow& window);

    void applyInteraction(const UIdata& uidata, const sf::RenderWindow& window);

    void moveFluid(const UIdata& uidata, sf::Sprite &sprite);

    void editWall(const UIdata& uidata);


private:
    FluidSim& fluidsim; // linking to fluidsim that is to be visualized
    double scale;
    std::vector<int> intBCMatrix;

    static constexpr int Nfields{5};
    std::array<MatrixTexture,Nfields> matrixTextures;
    std::array<sf::Sprite,Nfields> sprites{};

    ////////// For Interaction: //////////

    timetypes::tp lastMouseUpdate{};
    std::chrono::milliseconds mouseUpdateInterval{8};
    const sf::Sprite* last_hit_sprite_ptr{nullptr};
    sf::Vector2f mouseFieldPos{};
    sf::Vector2f lastMouseFieldPos{};
    sf::Vector2f mouseFieldPosDelta{};
    sf::Vector2f mouseVelocity{};

    double moveStrength{0.003};
    int moveBlobRadius{6};
    double moveClampValue{1.};

    sf::CircleShape brushCircle{0};

    ////////// Text: //////////

    sf::Font font{};
    std::array<sf::Text,Nfields> labels;

};


#endif //FLUIDSIM_FIELDSVIS_H