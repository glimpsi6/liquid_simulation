#include "Particles.h"

/*------------------------------------------------------------------*/

int Particles::cellIndex(float px, float py) const {
    float fx = std::max(0.0f, px * invCellSize);
    float fy = std::max(0.0f, py * invCellSize);
    fx = std::min(fx, maxCx);
    fy = std::min(fy, maxCy);
    return int(fy) * gridW + int(fx); 
}

/*------------------------------------------------------------------*/

void Particles::buildGrid() {
    const int nCells = gridW * gridH;

    // Гистограмма: сколько частиц в каждой ячейке
    std::fill(cellStart.begin(), cellStart.end(), 0); 
    for (int i = 0; i < N; ++i) {
        const int c = cellIndex(x[i], y[i]);
        cellOf[i] = c;                                  
        ++cellStart[c + 1];
    }

    // Префиксные суммы на месте
    for (int c = 0; c < nCells; ++c)
        cellStart[c + 1] += cellStart[c];
    assert(cellStart[nCells] == N); // каждая частица попала ровно в одну ячейку

    // Курсоры и новое место каждой частицы
    std::copy_n(cellStart.begin(), nCells, cursor.begin());
    for (int i = 0; i < N; ++i)
        dest[i] = cursor[cellOf[i]]++;  

    // Физическая перестановка данных
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        const int d = dest[i];
        x2[d]  = x[i];   y2[d]  = y[i];
        vx2[d] = vx[i];  vy2[d] = vy[i];
    }
    std::swap(x, x2);   std::swap(y, y2);
    std::swap(vx, vx2); std::swap(vy, vy2);
}

/*------------------------------------------------------------------*/

void Particles::startPos() {
    int cols = (int)std::sqrt(N);
    float spacing = h * 0.2f;
    float x0 = W_px * 0.3f;
    float y0 = H_px * 0.3f;
    for (int i = 0; i < N; ++i) {
        x[i] = x0 + (i % cols) * spacing;
        y[i] = y0 + (i / cols) * spacing;
    }
}

/*------------------------------------------------------------------*/

void Particles::proccessPositions() {
    const float h2   = h * h;
    const float VMAX = 2000.0f;

    buildGrid();

    // --- плотность и давление ---
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        const float xi = x[i], yi = y[i];
        float rhoi = mass * coreW(0.0f, h2); // самовклад
        int cx = std::clamp((int)(x[i] / cellSize), 0, gridW - 1);
        int cy = std::clamp((int)(y[i] / cellSize), 0, gridH - 1);
        for (int ry = -1; ry <= 1; ++ry) {
            const int ny = cy + ry;
            if (ny < 0 || ny >= gridH) continue;
            const int c0 = ny * gridW + std::max(cx - 1, 0);
            const int c1 = ny * gridW + std::min(cx + 1, gridW - 1);
            for (int j = cellStart[c0]; j < cellStart[c1 + 1]; ++j) {   // одна непрерывная полоса из 1–3 ячеек
                if (j == i) continue;
                float ddx = xi - x[j];
                float ddy = yi - y[j];
                float r2  = ddx*ddx + ddy*ddy;
                if (r2 >= h2) continue;
                rhoi += mass * coreW(r2, h2);
            }
        }
        rho[i]   = rhoi;
        press[i] = k * (rho[i] - rho0);
    }

    // --- силы ---
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        ax[i] = ay[i] = 0.0f;
        int cx = std::clamp((int)(x[i] / cellSize), 0, gridW - 1);
        int cy = std::clamp((int)(y[i] / cellSize), 0, gridH - 1);
        const float xi = x[i], yi = y[i], vxi = vx[i], vyi = vy[i];
        const float rhoi = rho[i];
        float pi_rho2 = press[i] / (rhoi * rhoi);
        
        for (int ry = -1; ry <= 1; ++ry) {
            const int ny = cy + ry;
            if (ny < 0 || ny >= gridH) continue;
            const int c0 = ny * gridW + std::max(cx - 1, 0);
            const int c1 = ny * gridW + std::min(cx + 1, gridW - 1);
            for (int j = cellStart[c0]; j < cellStart[c1 + 1]; ++j) {   // одна непрерывная полоса из 1–3 ячеек
                if (j == i) continue;
                float ddx = xi - x[j];
                float ddy = yi - y[j];
                float r2  = ddx*ddx + ddy*ddy;
                if (r2 >= h2) continue;

                float dw     = dW(r2, h2, C);
                float rho_ij = 0.5f * (rhoi + rho[j]);
                float pj_rho2 = press[j] / (rho[j] * rho[j]);
                float fp = -mass * (pi_rho2 + pj_rho2) * dw;
                ax[i] += fp * ddx;
                ay[i] += fp * ddy;

                float fv = - mu * mass / rho_ij * ddW(r2, h2, C);
                ax[i] += fv * (vx[j] - vxi);
                ay[i] += fv * (vy[j] - vyi);
            }
        }
    }

    // --- интеграция + защита ---
    #pragma omp parallel for
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

    // --- границы ---
    #pragma omp parallel for
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
    namespace fs = std::filesystem;
    fs::path p(filename);
    if (p.has_parent_path()) fs::create_directories(p.parent_path());

    FILE* f = std::fopen(filename, "wb");
    if (!f) {
        std::fprintf(stderr, "fopen('%s') failed: %s (cwd: %s)\n",
                     filename, std::strerror(errno), fs::current_path().c_str());
        return;
    }

    fprintf(f, "P6\n%d %d\n255\n", W_px, H_px);
    fwrite(pixels.data(), 1, pixels.size(), f);
    std::fclose(f);
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

float Particles::coreW(float r2, float h2) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return (4.0f / (Pi * h2 * h2)) * diff * diff * diff;
}

/*------------------------------------------------------------------*/

float Particles::dW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return C * diff * diff;
}

/*------------------------------------------------------------------*/

float Particles::ddW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    return C * (h2 - r2) * (3.0f * r2 - h2); //todo если вязкость будет странной то тут коэф
}

/*------------------------------------------------------------------*/