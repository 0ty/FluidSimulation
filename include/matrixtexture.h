#ifndef FLUIDSIM_MATRIXTEXTURE_H
#define FLUIDSIM_MATRIXTEXTURE_H

#include <vector>
#include <cmath>

#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "colormaps.h"


class MatrixTexture{
public:
    MatrixTexture(int width, int height);

    void updateTexture(const double matrix[], const int intBCMatrix[],
                       colormaps::ColorPallette colorPallette = colormaps::ColorPallette::Viridis,
                       bool fixedRange = false, double fixmin = 0, double fixmax = 1);

    const sf::Texture& getTexture() const;

private:
    const int width;
    const int height;
    std::vector<sf::Uint8> pixels;
    sf::Texture texture;

};

#endif //FLUIDSIM_MATRIXTEXTURE_H
