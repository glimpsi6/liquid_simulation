#include <cstdio>
#include <vector>
#include <algorithm>

// constexpr float W_world = 1.0f;
// constexpr float H_world = 1.0f;
constexpr int   W_px    = 800;
constexpr int   H_px    = 600;
#define Pi 3.14159f 

struct Particles {
    int N;
    std::vector<float> x, y;       // позиции
    std::vector<float> vx, vy;     // скорости
    std::vector<float> ax, ay;     // ускорения
    std::vector<float> rho;        // плотности
    std::vector<float> press;      // давления
    float h;                       // радиус сглаживания
    float rho0;                    // плотность покоя
    float mass;                    // масса частицы
    float k;                       // жёсткость (газовая постоянная)
    float mu;                      // вязкость
    float C;                       // нормировочный коэффициент
    float gx;                      // гравитация по x
    float gy;                      // гравитация по y
    float dt;                      // шаг времени

    Particles(int n) : N(n),
        x(n), y(n),
        vx(n), vy(n),
        ax(n), ay(n),
        rho(n), press(n) {
        h = 20.0f; 
        rho0 = 1.0f;
        mass = 1.0f;
        k = 50.0f;
        mu = 5.0f;
        C = -24.0f / (Pi * pow(h, 8));
        gx = 0.0f;
        gy = -9.8f * 100.0f;
        dt = 0.0005f;
    }
};

void writePPM(const char* filename, const std::vector<unsigned char>& pixels);
std::vector<unsigned char>& getPixels();
void savePPM(std::vector<unsigned char>& pixels);
void proccessPositions(Particles& pr);

void writePPM(const char* filename, const std::vector<unsigned char>& pixels) {
    FILE* f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", W_px, H_px);
    fwrite(pixels.data(), 1, pixels.size(), f);
    fclose(f);
}

void startPos(Particles& pr) {
    int cols = (int)std::sqrt(pr.N);
    int rows = (pr.N + cols - 1) / cols;
    for (int i = 0; i < pr.N; ++i) {
        pr.y[i] = i / cols;
        pr.x[i] = i % cols;
    }
}

std::vector<unsigned char>& getPixels() {
    static std::vector<unsigned char> pixels(W_px * W_px * 3, 0);
    return pixels;
}

void savePPM(std::vector<unsigned char>& pixels)
{
    static int frame_number = 0;

    char filename[64];
    snprintf(filename, sizeof(filename), "frames/frame%04d.ppm", frame_number);
    writePPM(filename, pixels);
    frame_number++;
}

float coreW(float r2, float h2) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return (4.0f / (Pi * h2 * h2)) * diff * diff * diff;
}

float dW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    float diff = h2 - r2;
    return C * diff * diff;
}

float ddW(float r2, float h2, float C) {
    if (r2 >= h2) return 0.0f;
    return C * (h2 - r2) * (3.0f * r2 - h2); //todo если вязкость будет странной то тут коэф
}

void proccessPositions(Particles& pr) {
    //плотность
    for(int i = 0; i < pr.N; i++) {
        pr.rho[i] = 0.0f;
        for(int j = 0; j < pr.N; j++) {
            if(i == j) continue;
            float dx = pr.x[i] - pr.x[j];
            float dy = pr.y[i] - pr.y[j];
            float r2 = dx * dx + dy * dy;
            pr.rho[i] += pr.mass * coreW(r2, pr.h * pr.h);
        } 
    }

    //давление
    for(int i = 0; i < pr.N; i++) {
        pr.press[i] = pr.k * (pr.rho[i] - pr.rho0);
        pr.press[i] = std::max(pr.press[i], 0.0f);
    }

    //силы
    for(int i = 0; i < pr.N; i++) {
        pr.ax[i] = 0;
        pr.ay[i] = 0;
        for(int j = 0; j < pr.N; j++) {
            float dx = pr.x[i] - pr.x[j];
            float dy = pr.y[i] - pr.y[j];
            float r2 = dx * dx + dy * dy;
            float h2 = pr.h * pr.h;
            if(r2 < h2) {
                float dw = dW(r2, h2, pr.C);
                float fp = -pr.mass * (pr.press[i] + pr.press[j]) / (2.0f * pr.rho[j]) * dw;
                pr.ax[i] += fp * dx;
                pr.ay[i] += fp * dy;
                //вязкость
                float fv = pr.mu * pr.mass / pr.rho[j] * ddW(r2, h2, pr.C); 
                //float fv = pr.mu * pr.mass / pr.rho[j] * coreW(r2, h2);
                pr.ax[i] += fv * (pr.vx[j] - pr.vx[i]);
                pr.ay[i] += fv * (pr.vy[j] - pr.vy[i]);
            }
        }
    }

    //интеграция
    for(int i = 0; i < pr.N; i++) {
        pr.ax[i] += pr.gx;
        pr.ay[i] += pr.gy;
        pr.vx[i] += pr.ax[i] / pr.rho[i] * pr.dt;
        pr.vy[i] += pr.ay[i] / pr.rho[i] * pr.dt;
        pr.x[i]  += pr.vx[i] * pr.dt;
        pr.y[i]  += pr.vy[i] * pr.dt;
    }
    
    //границы
    for(int i = 0; i < pr.N; i++) {
        if (pr.x[i] < 0.0f)  { pr.x[i] = 0.0f;   pr.vx[i] = -pr.vx[i] * 0.5f; }
        if (pr.x[i] > W_px)  { pr.x[i] = W_px; pr.vx[i] = -pr.vx[i] * 0.5f; }
        if (pr.y[i] < 0.0f)  { pr.y[i] = 0.0f;   pr.vy[i] = -pr.vy[i] * 0.5f; }
        if (pr.y[i] > H_px)  { pr.y[i] = H_px; pr.vy[i] = -pr.vy[i] * 0.5f; }
    }
}

void render(const Particles& pr, std::vector<unsigned char>& pixels) {
    std::fill(pixels.begin(), pixels.end(), 0);

    for (int i = 0; i < pr.N; ++i) {

        // int px = (int)(pr.x[i] / W_world * W_px);
        // int py = (int)(pr.y[i] / H_world * H_px);
        // py = H_px - 1 - py;

        // if (px < 0 || px >= W_px || py < 0 || py >= H_px) continue;

        int idx = (pr.y[i] * W_px + pr.x[i]) * 3;
        pixels[idx + 0] = 255;
        pixels[idx + 1] = 0;
        pixels[idx + 2] = 0;
    }
}

int main() {
    Particles particles(100);

    std::vector<unsigned char> pixels = getPixels();

    startPos(particles);
    render(particles, pixels);
    savePPM(pixels);

    const int   TOTAL_STEPS = 20000;   // сколько шагов физики
    const int   SAVE_EVERY  = 100;    // раз в сколько шагов писать кадр

    for (int step_id = 0; step_id < TOTAL_STEPS; ++step_id) {
        proccessPositions(particles);

        // время от времени — записать кадр
        if (step_id % SAVE_EVERY == 0) {
            render(particles, pixels);
            savePPM(pixels);
        }
    }

    return 0;
}