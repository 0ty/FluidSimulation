#include "matrixtexture.h"

MatrixTexture::MatrixTexture(int width, int height) : width{width}, height{height}, pixels(width*height*4,0) {
    texture.create(width, height);
    // texture.setSmooth(true);
}

void MatrixTexture::updateTexture(const double *matrix, const int *intBCMatrix, colormaps::ColorPallette colorPallette,
                                  bool fixedRange, double fixmin, double fixmax) {

    double min{matrix[0]};
    double max{matrix[0]};

    if (fixedRange) {
        min = fixmin;
        max = fixmax;

    } else {

        // min/max values are not given -> find min/max
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                int idx = (i + j*width);
                if (intBCMatrix[idx] == 1) continue;
                if (matrix[idx] > max) max = matrix[idx];
                if (matrix[idx] < min) min = matrix[idx];
            }
        }
        if (max == min) min -= 1.0;
    }

#pragma omp parallel for collapse(2)
    for(int i{0}; i < width; ++i){
        for(int j{0}; j < height; ++j){
            int idx = (i + j*width);
            int pixelidx = idx*4;

            double value = (matrix[idx]-min)/(max-min); // between [0,1]

            // Paint Wall black
            if(intBCMatrix[idx] == 1){
                pixels[pixelidx + 0] = 20;
                pixels[pixelidx + 1] = 20;
                pixels[pixelidx + 2] = 20;
                pixels[pixelidx + 3] = 255;
                continue;
            }

            // Paint NaN-value pink (i.e. when sim becomes unstable)
            if(std::isnan(value)){
                pixels[pixelidx + 0] = 255;
                pixels[pixelidx + 1] = 79;
                pixels[pixelidx + 2] = 226;
                pixels[pixelidx + 3] = 255;
                continue;
            }

            auto [cm_red, cm_blue, cm_green] = colormaps::valueToRGB(value, colorPallette);

            pixels[pixelidx + 0] = 255 * cm_red;
            pixels[pixelidx + 1] = 255 * cm_blue;
            pixels[pixelidx + 2] = 255 * cm_green;
            pixels[pixelidx + 3] = 255 * 1;
        }
    }

    texture.update(pixels.data());
}

const sf::Texture &MatrixTexture::getTexture() const {
    return texture;
}
