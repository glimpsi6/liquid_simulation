#include "Particles.h"

/*------------------------------------------------------------------*/

void Particles::buildGrid() {
    for (auto& cell : grid) cell.clear();
    
    for (int i = 0; i < N; ++i) {
        int cx = (int)(x[i] / cellSize);
        int cy = (int)(y[i] / cellSize);
        cx = std::clamp(cx, 0, gridW - 1);
        cy = std::clamp(cy, 0, gridH - 1);
        grid[cy * gridW + cx].push_back(i);
    }
}

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
    const float h2 = h * h;
    const float VMAX = 2000.0f;
    
    buildGrid();
    
    // плотность
    for (int i = 0; i < N; ++i) {
        rho[i] = mass * coreW(0.0f, h2);  // самовклад
        int cx = std::clamp((int)(x[i] / cellSize), 0, gridW - 1);
        int cy = std::clamp((int)(y[i] / cellSize), 0, gridH - 1);
        for (int dy = -1; dy <= 1; ++dy) {
            int ny = cy + dy; if (ny < 0 || ny >= gridH) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int nx = cx + dx; if (nx < 0 || nx >= gridW) continue;
                for (int j : grid[ny * gridW + nx]) {
                    if (i == j) continue;
                    float ddx = x[i] - x[j];
                    float ddy = y[i] - y[j];
                    float r2 = ddx*ddx + ddy*ddy;
                    if (r2 >= h2) continue;
                    rho[i] += mass * coreW(r2, h2);
                }
            }
        }
    }
    
    // давление
    for (int i = 0; i < N; ++i) {
        press[i] = std::max(k * (rho[i] - rho0), 0.0f);
    }
    
    // силы
    for (int i = 0; i < N; ++i) {
        ax[i] = ay[i] = 0.0f;
        int cx = std::clamp((int)(x[i] / cellSize), 0, gridW - 1);
        int cy = std::clamp((int)(y[i] / cellSize), 0, gridH - 1);
        for (int dy = -1; dy <= 1; ++dy) {
            int ny = cy + dy; if (ny < 0 || ny >= gridH) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int nx = cx + dx; if (nx < 0 || nx >= gridW) continue;
                for (int j : grid[ny * gridW + nx]) {
                    if (i == j) continue;
                    float ddx = x[i] - x[j];
                    float ddy = y[i] - y[j];
                    float r2 = ddx*ddx + ddy*ddy;
                    if (r2 >= h2) continue;
                    
                    float dw = dW(r2, h2, C);
                    float rho_ij = 0.5f * (rho[i] + rho[j]);
                    float fp = -mass * (press[i] + press[j]) / (2.0f * rho_ij) * dw;
                    ax[i] += fp * ddx;
                    ay[i] += fp * ddy;
                    
                    float fv = mu * mass / rho_ij * ddW(r2, h2, C);
                    ax[i] += fv * (vx[j] - vx[i]);
                    ay[i] += fv * (vy[j] - vy[i]);
                }
            }
        }
    }
    
    // интеграция + защита
    for (int i = 0; i < N; ++i) {
        ax[i] += gx;
        ay[i] += gy;
        float invrho = (rho[i] > 1e-6f) ? 1.0f / rho[i] : 0.0f;
        vx[i] += ax[i] * invrho * dt;
        vy[i] += ay[i] * invrho * dt;
        if (!std::isfinite(vx[i])) vx[i] = 0.0f;
        if (!std::isfinite(vy[i])) vy[i] = 0.0f;
        vx[i] = std::clamp(vx[i], -VMAX, VMAX);
        vy[i] = std::clamp(vy[i], -VMAX, VMAX);
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
    }
    
    // границы
    for (int i = 0; i < N; ++i) {
        if (x[i] < 0.0f)  { x[i] = 0.0f;  vx[i] = -vx[i] * 0.5f; }
        if (x[i] > W_px)  { x[i] = W_px;  vx[i] = -vx[i] * 0.5f; }
        if (y[i] < 0.0f)  { y[i] = 0.0f;  vy[i] = -vy[i] * 0.5f; }
        if (y[i] > H_px)  { y[i] = H_px;  vy[i] = -vy[i] * 0.5f; }
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