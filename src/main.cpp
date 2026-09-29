#include "Particles.h"

int main() {
    Particles pr(1000);

    std::vector<unsigned char> pixels = getPixels();

    pr.startPos();
    pr.render(pixels);
    savePPM(pixels);

    const int   TOTAL_STEPS = 300000;   // сколько шагов физики
    const int   SAVE_EVERY  = 1000;    // раз в сколько шагов писать кадр

    for (int step_id = 0; step_id < TOTAL_STEPS; ++step_id) {
        pr.proccessPositions();

        // записать кадр
        if (step_id % SAVE_EVERY == 0) {
            pr.render(pixels);
            savePPM(pixels);
        }
    }

    return 0;
}