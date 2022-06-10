#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <glm/glm.hpp>
#include <omp.h>
#include "../include/Camera.h"
#include "../include/Scenes.h"
#include "../include/Render.h"
//#include</usr/local/opt/libomp/include/omp.h>

static const int WIDTH = 1024 / 2;
static const int HEIGHT = 768 / 2;
static const int SUB_SAMPLING_CAUSTICS = 8;
static const int SUB_SAMPLING_MONTE_CARLO = 100;
static const int SUB_SAMPLING_DIRECT_SPECULAR = 20;
static const int NUMBER_OF_PHOTONS_EMISSION = 1000000;

int savePPM(const char* file_name, const int width, const int height, unsigned char* data) {
    FILE *fp = fopen(file_name, "wb"); // b - binary mode
    fprintf(fp, "P6\n%d %d\n255\n", width, height);
    fwrite(data, 1, width * height * 3, fp);
    fclose(fp);
    return EXIT_SUCCESS;
}

const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d-%X", &tstruct);
    return buf;
}

int main(int argc, char const *argv[]) {
    time_t time_start, time_now, rendertime_start;
    time(&time_start);

    // The camera is used to cast appropriate initial rays
    Camera c(glm::vec3(0, 0, 3.2), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0), M_PI / 3, WIDTH, HEIGHT);
    glm::vec3 camera_plane_normal = glm::normalize(c.center - c.eye);

    // 3D objects are contained in the Scene object
    Scene s(argv[1]);
    Render render(&s);

    SpectralDistribution* irradiance_values = new SpectralDistribution[c.WIDTH * c.HEIGHT];
    unsigned char* pixels = new unsigned char[c.WIDTH * c.HEIGHT * 3]; // w * h * rgb

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-0.5, 0.5);

    std::cerr << "Building photon map ...\n";
    render.buildPhotonMap(NUMBER_OF_PHOTONS_EMISSION);

    std::cerr << "\n------------------ Rendering started ------------------\n" << std::endl;

    time(&rendertime_start);

    double prerender_time = difftime(rendertime_start, time_start);

    for (int x = 0; x < c.WIDTH; ++x) {
        // Parallellize the for loop with openMP.
        #pragma omp parallel for
        for (int y = 0; y < c.HEIGHT; ++y) {
            int ind = (x + y * c.WIDTH);
            SpectralDistribution sd;
            if (SUB_SAMPLING_DIRECT_SPECULAR) {
                for (int i = 0; i < SUB_SAMPLING_DIRECT_SPECULAR; ++i) {
                    Ray r = c.castRay(x, (c.HEIGHT - y - 1), dis(gen), dis(gen));
                    sd += render.traceRay(r, Render::WHITTED_SPECULAR) * glm::dot(r.direction, camera_plane_normal);
                }
                irradiance_values[ind] += sd / SUB_SAMPLING_DIRECT_SPECULAR * (2.f * M_PI);
            }
            sd = SpectralDistribution();
            if (SUB_SAMPLING_CAUSTICS) {
                for (int i = 0; i < SUB_SAMPLING_CAUSTICS; ++i) {
                    Ray r = c.castRay(x, (c.HEIGHT - y - 1), dis(gen), dis(gen));
                    sd += render.traceRay(r, Render::CAUSTICS) * glm::dot(r.direction, camera_plane_normal);
                }
                irradiance_values[ind] += sd / SUB_SAMPLING_CAUSTICS * (2.f * M_PI);
            }
            sd = SpectralDistribution();
            if (SUB_SAMPLING_MONTE_CARLO) {
                for (int i = 0; i < SUB_SAMPLING_MONTE_CARLO; ++i) {
                    Ray r = c.castRay(x, (c.HEIGHT - y - 1), dis(gen), dis(gen));
                    sd += render.traceRay(r, Render::MONTE_CARLO) * glm::dot(r.direction, camera_plane_normal);
                }
                irradiance_values[ind] += sd / SUB_SAMPLING_MONTE_CARLO * (2.f * M_PI);
            }
        }

		if((x + 1) % 16 == 0) {
			float rendering_percent = (x + 1) * 100.f / c.WIDTH;
			time(&time_now);
			double tL = difftime(time_now, rendertime_start);
			double tR = (tL / rendering_percent) * (100 - rendering_percent);
			int h = tR / (60 * 60);
			int m = (int(tR) % (60 * 60)) / 60;
			int s = int(tR) % 60;

			std::cerr << rendering_percent << " \% of rendering finished. \n";
			std::cerr << "Estimated time left is  " << h << " : " << m << " : " << s << "\n";
		}
    }

    std::cerr << "\n------------------ Rendering ended ------------------\n" << std::endl;

    // Convert to byte data
    // Gamma correction
    float gamma = 1 / 2.2;
    for (int x = 0; x < c.WIDTH; ++x) {
        for (int y = 0; y < c.HEIGHT; ++y) {
            int ind = (x + y * c.WIDTH);
			for (int o = 0; o < 3; ++o) 
            	pixels[ind * 3 + o] = char(int(glm::clamp(glm::pow(irradiance_values[ind][o],gamma), 0.0f, 1.0f) * 255));
        }
    }

    std::string date_time = currentDateTime();
    std::string file_name = date_time + ".ppm";

    // Save the image data to file
    savePPM(file_name.c_str(), WIDTH, HEIGHT, pixels);

    // Save information in a text file
    std::ofstream F(date_time + ".txt");
    F << "Rendered file information:\n\n";
    F << "File name                    : " + date_time + ".ppm\n";
    F << "Resolution                   : " + std::to_string(WIDTH) + " x " + std::to_string(HEIGHT) + "\n";
    F << "Caustic sub sampling         : " + std::to_string(SUB_SAMPLING_CAUSTICS) + "\n";
    F << "Monte Carlo sub sampling     : " + std::to_string(SUB_SAMPLING_MONTE_CARLO) + "\n";
    F << "Direct specular sub sampling : " + std::to_string(SUB_SAMPLING_DIRECT_SPECULAR) + "\n";
    F << "Emitted photons              : " + std::to_string(NUMBER_OF_PHOTONS_EMISSION) + "\n";
    F << "Objects in scene             : " + std::to_string(s.getNumberOfObjects()) + "\n";
    F << "Spheres in scene             : " + std::to_string(s.getNumberOfSpheres()) + "\n";
    F << "Triangles in scene           : " + std::to_string(s.getNumberOfTriangles()) + "\n";
    F << "Gamma                        : " + std::to_string(gamma) + "\n";
    F.close();


    delete [] irradiance_values;
    delete [] pixels;

    std::cout << '\a';
    return EXIT_SUCCESS;
}