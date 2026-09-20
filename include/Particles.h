#ifndef PARTICLES_H
#define PARTICLES_H

#include <cstdio>
#include <vector>
#include <algorithm>

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
        gy = 9.8f * 100.0f;
        dt = 0.0005f;
    }

    void startPos();
    void proccessPositions();
    void render(std::vector<unsigned char>& pixels);
};

void writePPM(const char* filename, const std::vector<unsigned char>& pixels);
std::vector<unsigned char>& getPixels();
void savePPM(std::vector<unsigned char>& pixels);
float coreW(float r2, float h2);
float dW(float r2, float h2, float C);
float ddW(float r2, float h2, float C);

#endif // PARTICLES_H