#ifndef PARTICLES_H
#define PARTICLES_H

#include <cstdio>
#include <vector>
#include <algorithm>

constexpr int   W_px    = 800;
constexpr int   H_px    = 600;
#define Pi 3.14159f 

struct Particles {
public:
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

    int cellSize;
    int gridW, gridH;
    std::vector<int> cellCount;     // сколько частиц в каждой ячейке
    std::vector<int> cellStart;     // индекс начала ячейки в particleIds
    std::vector<int> particleIds;   // id частиц, уложенные по ячейкам
    std::vector<int> cursor;        // временный, для второго прохода

    Particles(int n) : N(n),
        x(n), y(n),
        vx(n), vy(n),
        ax(n), ay(n),
        rho(n), press(n) {
        h = 10.0f; 
        rho0 = 203.0f;
        mass = 1.0f;
        k = 15000.0f;
        mu = 0.0f;
        C = -24.0f / (Pi * pow(h, 8));
        gx = 0.0f;
        gy = 9.8f;
        dt = 0.01f;

        cellSize = (int)h;
        gridW = (W_px / cellSize) + 2;   // +2 = запас на границу
        gridH = (H_px / cellSize) + 2;
        cellCount.assign(gridW * gridH, 0);
        cellStart.assign(gridW * gridH, 0);
        cursor.assign(gridW * gridH, 0);
        particleIds.reserve(N);
    }

    void startPos();
    void proccessPositions();
    void render(std::vector<unsigned char>& pixels);

private:
    float coreW(float r2, float h2);
    float dW(float r2, float h2, float C);
    float ddW(float r2, float h2, float C);
    void buildGrid();
};

void writePPM(const char* filename, const std::vector<unsigned char>& pixels);
std::vector<unsigned char>& getPixels();
void savePPM(std::vector<unsigned char>& pixels);

#endif // PARTICLES_H