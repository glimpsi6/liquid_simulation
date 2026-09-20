#include "Particles.h"

/*------------------------------------------------------------------*/

void Particles::startPos() {
    int cols = (int)std::sqrt(N);
    float spacing = h * 0.4f;
    float x0 = W_px * 0.3f;
    float y0 = H_px * 0.3f;
    for (int i = 0; i < N; ++i) {
        x[i] = x0 + (i % cols) * spacing;
        y[i] = y0 + (i / cols) * spacing;
    }
}

/*------------------------------------------------------------------*/

void Particles::proccessPositions() {
    //плотность
    for(int i = 0; i < N; i++) {
        rho[i] = 0.0f;
        for(int j = 0; j < N; j++) {
            float dx = x[i] - x[j];
            float dy = y[i] - y[j];
            float r2 = dx * dx + dy * dy;
            rho[i] += mass * coreW(r2, h * h);
        } 
    }

    //давление
    for(int i = 0; i < N; i++) {
        press[i] = k * (rho[i] - rho0);
        press[i] = std::max(press[i], 0.0f);
    }

    //силы
    for(int i = 0; i < N; i++) {
        ax[i] = 0;
        ay[i] = 0;
        for(int j = 0; j < N; j++) {
            float dx = x[i] - x[j];
            float dy = y[i] - y[j];
            float r2 = dx * dx + dy * dy;
            float h2 = h * h;
            if(r2 < h2) {
                float dw = dW(r2, h2, C);
                float rho_ij = 2.0f * rho[i] * rho[j] / (rho[i] + rho[j]);
                float fp = -mass * (press[i] + press[j]) / (2.0f * rho_ij) * dw;
                ax[i] += fp * dx;
                ay[i] += fp * dy;
                //вязкость
                float fv = mu * mass / rho[j] * ddW(r2, h2, C); 
                ax[i] += fv * (vx[j] - vx[i]);
                ay[i] += fv * (vy[j] - vy[i]);
            }
        }
    }

    //интеграция
    for(int i = 0; i < N; i++) {
        ax[i] += gx;
        ay[i] += gy;

        // защита от деления на ~0
        float invrho = (rho[i] > 1e-6f) ? 1.0f / rho[i] : 0.0f;
        vx[i] += ax[i] * invrho * dt;
        vy[i] += ay[i] * invrho * dt;

        // защита от NaN/inf (на случай, если что-то всё же проскочило)
        if (!std::isfinite(vx[i])) vx[i] = 0.0f;
        if (!std::isfinite(vy[i])) vy[i] = 0.0f;

        // клэмп скорости, чтобы частица не улетела на тысячи пикселей за кадр
        const float VMAX = 2000.0f;
        vx[i] = std::clamp(vx[i], -VMAX, VMAX);
        vy[i] = std::clamp(vy[i], -VMAX, VMAX);

        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
    }   

    //границы
    for(int i = 0; i < N; i++) {
        if (x[i] < 0.0f)  { x[i] = 0.0f;   vx[i] = -vx[i] * 0.5f; }
        if (x[i] > W_px)  { x[i] = W_px;   vx[i] = -vx[i] * 0.5f; }
        if (y[i] < 0.0f)  { y[i] = 0.0f;   vy[i] = -vy[i] * 0.5f; }
        if (y[i] > H_px)  { y[i] = H_px;   vy[i] = -vy[i] * 0.5f; }
    }
}

/*------------------------------------------------------------------*/

void Particles::render(std::vector<unsigned char>& pixels) {
    std::fill(pixels.begin(), pixels.end(), 0);
    for (int i = 0; i < N; ++i) {
        int px = (int)x[i];
        int py = (int)y[i];
        if (px < 0 || px >= W_px || py < 0 || py >= H_px) continue;
        int idx = (py * W_px + px) * 3;
        pixels[idx + 0] = 255;
        pixels[idx + 1] = 0;
        pixels[idx + 2] = 0;
    }
}

/*------------------------------------------------------------------*/

void writePPM(const char* filename, const std::vector<unsigned char>& pixels) {
    FILE* f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", W_px, H_px);
    fwrite(pixels.data(), 1, pixels.size(), f);
    fclose(f);
}

/*------------------------------------------------------------------*/

std::vector<unsigned char>& getPixels() {
    static std::vector<unsigned char> pixels(W_px * H_px * 3, 0);
    return pixels;
}

/*------------------------------------------------------------------*/

void savePPM(std::vector<unsigned char>& pixels)
{
    static int frame_number = 0;

    char filename[64];
    snprintf(filename, sizeof(filename), "frames/frame%04d.ppm", frame_number);
    writePPM(filename, pixels);
    frame_number++;
}

/*------------------------------------------------------------------*/

float coreW(float r2, float h2) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return (4.0f / (Pi * h2 * h2)) * diff * diff * diff;
}

/*------------------------------------------------------------------*/

float dW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return C * diff * diff;
}

/*------------------------------------------------------------------*/

float ddW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    return C * (h2 - r2) * (3.0f * r2 - h2); //todo если вязкость будет странной то тут коэф
}

/*------------------------------------------------------------------*/