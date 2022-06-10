#include "../include/Camera.h"
#include <iostream>
const float D = 0.5f;

Camera::Camera (glm::vec3 eye, glm::vec3 center, glm::vec3 up, float fov, int width, int height) : eye(eye), center(center), up(up), fov(fov), WIDTH(width), HEIGHT(height) {
	auto V = glm::lookAt(eye, center, up);
	auto P = glm::perspective(fov, float(WIDTH) / HEIGHT, 0.1f, 100.f);
	VP_inv = glm::inverse(V * P);
}

Ray Camera::castRay(int pixel_x, int pixel_y, float parameter_x, float parameter_y) {
	Ray r;
	r.origin = r.direction = glm::vec3(0);
	if (pixel_x < 0 || pixel_x > WIDTH - 1 || pixel_y < 0 || pixel_y > HEIGHT - 1 || std::fabs(parameter_x) > 0.5 || std::fabs(parameter_y) > 0.5) {
		std::cout << "error in castRay\n";
		return r;
	}
	
	glm::vec4 fr = VP_inv * glm::vec4(((pixel_x + parameter_x) / WIDTH - 0.5) * 2, ((pixel_y + parameter_y) / HEIGHT - 0.5) * 2, 1, 1 );
	glm::vec4 to = VP_inv * glm::vec4(((pixel_x + parameter_x) / WIDTH - 0.5) * 2, ((pixel_y + parameter_y) / HEIGHT - 0.5) * 2, -1, 1 );
	glm::vec3 direction = glm::normalize(glm::vec3(to) * to.w - glm::vec3(fr) * fr.w);
	r.origin = eye;
	r.direction = direction;
	r.material = Material::air();
	r.radiance = SpectralDistribution();
	r.radiance[0] = r.radiance[1] = r.radiance[2] = 1;
	return r;
}