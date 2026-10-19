#ifndef PARTICLES_H
#define PARTICLES_H

#include <cstdio>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <cassert>

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
    float invCellSize, maxCx, maxCy;
    std::vector<float> x2, y2, vx2, vy2;   // размер N
    std::vector<int>   cellStart;          // размер nCells + 1
    std::vector<int>   cursor;             // размер nCells
    std::vector<int>   cellOf, dest;       // размер N

    Particles(int n) : N(n),
        x(n), y(n),
        vx(n), vy(n),
        ax(n), ay(n),
        rho(n), press(n) {
        h = 10.0f; 
        rho0 = 203.0f;
        mass = 0.5f;
        k = 15000.0f;
        mu = 0.0f;
        C = -24.0f / (Pi * pow(h, 8));
        gx = 0.0f;
        gy = 9.8f;
        dt = 0.01f;

        cellSize = (int)h;
        gridW = (W_px / cellSize) + 2;   // +2 = запас на границу
        gridH = (H_px / cellSize) + 2;
        const int nCells = gridW * gridH;
        cellStart.assign(nCells + 1, 0);
        cursor.assign(nCells, 0);
        cellOf.assign(n, 0);
        dest.assign(n, 0);
        x2.assign(n, 0.0f);  y2.assign(n, 0.0f);
        vx2.assign(n, 0.0f); vy2.assign(n, 0.0f);

        invCellSize = 1.0f / cellSize;  
        maxCx = float(gridW - 1);  
        maxCy = float(gridH - 1);
    }

    void startPos();
    void proccessPositions();
    void render(std::vector<unsigned char>& pixels);

private:
    float coreW(float r2, float h2);
    float dW(float r2, float h2, float C);
    float ddW(float r2, float h2, float C);
    int cellIndex(float px, float py) const;
    void buildGrid();
};

void writePPM(const char* filename, const std::vector<unsigned char>& pixels);
std::vector<unsigned char>& getPixels();
void savePPM(std::vector<unsigned char>& pixels);

#endif // PARTICLES_H