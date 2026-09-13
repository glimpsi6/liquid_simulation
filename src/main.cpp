#include <cstdio>
#include <vector>

#define W 800
#define H 600

struct Particles {
    int N;
    std::vector<float> x;
    std::vector<float> y;

    Particles(int n) : N(n), x(n, 0.0f), y(n, 0.0f) {}
};

void writePPM(const char* filename, const std::vector<unsigned char>& pixels);
void initParticlesInx(Particles& pr, std::vector<unsigned char>& pixels);
std::vector<unsigned char>& getPixels();
void savePPM(std::vector<unsigned char>& pixels);

void writePPM(const char* filename, const std::vector<unsigned char>& pixels) {
    FILE* f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    fwrite(pixels.data(), 1, pixels.size(), f);
    fclose(f);
}

void initParticlesInx(Particles& pr, std::vector<unsigned char>& pixels) {
    int cols = sqrt(pr.N);
    int rows = (pr.N + cols - 1) / cols;
    for(int i = 0; i < pr.N; i++) {
        int row = i / cols;
        int col = i % cols;
        pr.x[i] = col;
        pr.y[i] = row;

        int idx = (row * W + col) * 3;
        pixels[idx + 0] = 255;  // R
        pixels[idx + 1] = 0;    // G
        pixels[idx + 2] = 0;    // B
    }
}

std::vector<unsigned char>& getPixels() {
    static std::vector<unsigned char> pixels(W * H * 3, 0);
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

int main() {
    Particles particles(50);

    std::vector<unsigned char> pixels = getPixels();

    initParticlesInx(particles, pixels);
    savePPM(pixels);

    return 0;
}