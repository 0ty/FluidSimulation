#include "fieldsvis.h"

FieldsVis::FieldsVis(FluidSim &fluidsim, double scale)
        : fluidsim{fluidsim},
          scale{scale},
          intBCMatrix(fluidsim.width*fluidsim.height),
          matrixTextures{
                  MatrixTexture{fluidsim.width,fluidsim.height},
                  MatrixTexture{fluidsim.width,fluidsim.height},
                  MatrixTexture{fluidsim.width,fluidsim.height},
                  MatrixTexture{fluidsim.width,fluidsim.height},
                  MatrixTexture{fluidsim.width,fluidsim.height}
          }
{

    for(int i=0; i<Nfields; ++i){
        sprites[i].setTexture(matrixTextures[i].getTexture());
        sprites[i].setScale(scale, scale);
    }
    arrangeSprites();

    brushCircle.setFillColor(sf::Color::Transparent);
    brushCircle.setOutlineColor(sf::Color{255,255,255,100});
    brushCircle.setOutlineThickness(2);

    font.loadFromFile("../resources/Inconsolata-Bold.ttf");
    for(sf::Text& label : labels){
        label.setFont(font);
        label.setCharacterSize(14);
        label.setOutlineThickness(1);
    }
    labels[0].setString("Density rho");
    labels[1].setString("Temperature T");
    labels[2].setString("Pressure p");
    labels[3].setString("Velocity u");
    labels[4].setString("Velocity v");
}

void FieldsVis::arrangeSprites(const sf::Vector2f &offset, int spritesPerRow) {

    int x_idx{0};
    int y_idx{0};

    for(int i=0; i<Nfields; ++i){
        x_idx = i%spritesPerRow;
        y_idx = i/spritesPerRow;

        ///// sprites position: /////

        sf::Vector2f pos{static_cast<float>(offset.x + x_idx * fluidsim.width * scale),
                         static_cast<float>(offset.y + y_idx * fluidsim.height * scale)};
        sprites[i].setPosition(pos);

        ///// label position: /////

        sf::Vector2f labelCornerVec{0.f,static_cast<float>(fluidsim.height * scale - labels[i].getCharacterSize())};
        sf::Vector2f labelMargin{0.08f*fluidsim.width, -0.08f*fluidsim.height};

        labels[i].setPosition(pos + labelCornerVec + labelMargin);

    }
}

void FieldsVis::drawFields(sf::RenderWindow &window) {

    fluidsim.extractBCMatrix(intBCMatrix);

    matrixTextures[0].updateTexture(fluidsim.rho.data(), intBCMatrix.data(), colormaps::ColorPallette::Viridis);
    matrixTextures[1].updateTexture(fluidsim.T.data(), intBCMatrix.data(), colormaps::ColorPallette::Magma);
    matrixTextures[2].updateTexture(fluidsim.p.data(), intBCMatrix.data(), colormaps::ColorPallette::Parula);
    matrixTextures[3].updateTexture(fluidsim.u.data(), intBCMatrix.data(), colormaps::ColorPallette::RtoG, true, -1, 1);
    matrixTextures[4].updateTexture(fluidsim.v.data(), intBCMatrix.data(), colormaps::ColorPallette::RtoG, true, -1, 1);

    for(const sf::Sprite& sprite : sprites) window.draw(sprite);

    for(const sf::Text& label : labels) window.draw(label);

}

void FieldsVis::drawBrushCircle(const UIdata &uidata, sf::RenderWindow &window) {

    if(uidata.interactionMode == UIdata::INTERACTION_MODE::EDIT_WALL){
        brushCircle.setRadius(scale*uidata.brushSize*0.95);
        brushCircle.setOrigin(brushCircle.getRadius(),brushCircle.getRadius());
        sf::Vector2i mousePixelPos = sf::Mouse::getPosition(window);
        sf::Vector2f mouseWorldPos = window.mapPixelToCoords(mousePixelPos);
        brushCircle.setPosition(mouseWorldPos);
        window.draw(brushCircle);
    }

}

void FieldsVis::applyInteraction(const UIdata &uidata, const sf::RenderWindow &window) {

    if ((timetypes::sclock::now() - lastMouseUpdate > mouseUpdateInterval) &&
        (uidata.mouseLeftPressed || uidata.mouseRightPressed)) {

        sf::Vector2i mousePixelPos = sf::Mouse::getPosition(window);
        sf::Vector2f mouseWorldPos = window.mapPixelToCoords(mousePixelPos);

        // Check all sprites/fields if our mouse is inside it -> apply interaction to this sprite/field:
        for (sf::Sprite &sprite: sprites) {

            // Mouse position (in terms of index i,j) within the current field
            mouseFieldPos = sprite.getInverseTransform().transformPoint(mouseWorldPos);

            // If mouse is within the current field (not outside):
            if (mouseFieldPos.x >= 0 && mouseFieldPos.x <= fluidsim.width &&
                mouseFieldPos.y >= 0 && mouseFieldPos.y <= fluidsim.height) {

                switch (uidata.interactionMode) {
                    case UIdata::INTERACTION_MODE::MOVE_FLUID:
                        moveFluid(uidata, sprite);
                        break;

                    case UIdata::INTERACTION_MODE::EDIT_WALL:
                        editWall(uidata);
                        break;

                }

                last_hit_sprite_ptr = &sprite;
                lastMouseFieldPos = mouseFieldPos;
                break; // If we found that the mouse is in one sprite -> no need to check remaining ones
            }


        }

        lastMouseUpdate = timetypes::sclock::now();
    }

}

void FieldsVis::moveFluid(const UIdata &uidata, sf::Sprite &sprite) {

    // For moving the fluid we need the change of the mouse position (Delta) between the last and current step.
    // But if the mouse was outside the current field (it was in a different field) in the last step,
    // we cannot compute the delta. -> Only apply this interaction if the mouse is in the same sprite as the previous step.
    if (&sprite == last_hit_sprite_ptr) {
        mouseFieldPosDelta = mouseFieldPos - lastMouseFieldPos;

        auto update_dt = timetypes::sclock::now() - lastMouseUpdate;
        auto dt_seconds = (std::chrono::duration<float, std::ratio<1>>(update_dt)).count();

        mouseVelocity = mouseFieldPosDelta / dt_seconds;

        fluidsim.addBlobToField(fluidsim.u, mouseFieldPos.y, mouseFieldPos.x,
                                mouseVelocity.y * moveStrength, moveBlobRadius, moveClampValue);
        fluidsim.addBlobToField(fluidsim.v, mouseFieldPos.y, mouseFieldPos.x,
                                mouseVelocity.x * moveStrength, moveBlobRadius, moveClampValue);
    }

}

void FieldsVis::editWall(const UIdata &uidata) {

    // left mouse pressed -> add wall
    if (uidata.mouseLeftPressed) {
        fluidsim.addBounaryCondition(FluidSim::BC::WALL, mouseFieldPos.y, mouseFieldPos.x, uidata.brushSize);
    }

    // right mouse pressed -> remove wall
    if (uidata.mouseRightPressed) {
        fluidsim.addBounaryCondition(FluidSim::BC::NONE, mouseFieldPos.y, mouseFieldPos.x, uidata.brushSize);
    }

}
