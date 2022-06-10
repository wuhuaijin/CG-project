#include "../include/Object3D.h"
#include "../external_libraries/common_include/objloader.h"
#include "../external_libraries/common_include/vboindexer.h"
#include <random>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

const float eps = 1e-5;

// --- Object3D class functions --- //

Object3D::Object3D(Material* material) : material_(material) {}

Material Object3D::material() const {
	return material_ ? *material_ : Material();
}

// --- Mesh class functions --- //

Mesh::Mesh(glm::mat4 transform, const char* file_path, Material * material) : Object3D(material) {
	transform_ = transform;
	std::vector<glm::vec3> tmp_positions;
	std::vector<glm::vec2> tmp_uvs;
	std::vector<glm::vec3> tmp_normals;
	if(!loadOBJ(file_path, tmp_positions, tmp_uvs, tmp_normals))
		exit (EXIT_FAILURE);
	for (int i = 0; i < tmp_positions.size(); i++) {
		tmp_positions[i] = glm::vec3(transform_ * glm::vec4(tmp_positions[i], 1));
		tmp_normals[i] = glm::vec3(transform_ * glm::vec4(tmp_normals[i], 0));
	}
	indexVBO(tmp_positions, tmp_uvs, tmp_normals, indices_, positions_, uvs_, normals_);

	std::cout << "Building octree for mesh." << std::endl;
	ot_aabb_ = new OctTreeAABB(this);
	std::cout << "Octree built." << std::endl;
}

bool Mesh::intersect(IntersectionData* id, Ray r) const {
	return ot_aabb_->intersect(id, r);
}

glm::mat4 Mesh::getTransform() const {
	return transform_;
}

glm::vec3 Mesh::getMinPosition() const {
	auto ans = positions_[0];
	for (auto it : positions_) {
		ans.x = std::min(it.x, ans.x);
		ans.y = std::min(it.y, ans.y);
		ans.z = std::min(it.z, ans.z);
	}
	return ans;
}

glm::vec3 Mesh::getMaxPosition() const {
	auto ans = positions_[0];
	for (auto it : positions_) {
		ans.x = std::max(it.x, ans.x);
		ans.y = std::max(it.y, ans.y);
		ans.z = std::max(it.z, ans.z);
	}
	return ans;
}

int Mesh::getNumberOfTriangles() const {
	return indices_.size() / 3;
}

// --- Sphere class functions --- //

Sphere::Sphere(glm::vec3 position, float radius, Material* material) : Object3D(material), POSITION_(position), RADIUS_(radius) {}

bool Sphere::intersect(IntersectionData* id, Ray r) const {
	float p_half = glm::dot((r.origin - POSITION_), r.direction);
	float to_square = pow(p_half, 2) + pow(RADIUS_, 2) - pow(glm::length(r.origin - POSITION_), 2);
	if (to_square < 0)
		return false;
	float t = (-p_half < sqrt(to_square)) ? -p_half + sqrt(to_square) : -p_half - sqrt(to_square);
	if (t < 0)
		return false;
	glm::vec3 n = r.origin + t * r.direction - POSITION_;
	id->t = t;
	id->normal = glm::normalize(n);
	id->material = material();
	return true;
}

glm::vec3 Sphere::getPointOnSurface(float u, float v) const {
	float inclination = glm::acos(1 - 2 * u);
	float azimuth = 2 * M_PI * v;
	auto rand_dir = glm::vec3(1, 0, 0);
	rand_dir = glm::normalize(glm::rotate(rand_dir, inclination, glm::vec3(0,1,0)));
	rand_dir = glm::normalize(glm::rotate(rand_dir, azimuth, glm::vec3(1, 0, 0)));
	return POSITION_ + rand_dir * RADIUS_;
}

// --- Plane class functions --- //

Plane::Plane(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, Material* material) : Object3D(material), P0_(p0), P1_(p1), P2_(p2), NORMAL_(glm::normalize(glm::cross(p0 - p1, p0 - p2))), AREA_(glm::length(glm::cross(p0 - p1, p0 - p2))) {}

bool Plane::intersect(IntersectionData* id, Ray r) const {
	auto e1 = P1_ - P0_, e2 = P2_ - P0_;
	auto P = glm::cross(r.direction, e2);
	auto Q = glm::cross(r.origin - P0_, e1);
	auto det = glm::dot(e1, P);
	if (std::fabs(det) < eps)
		return false;
	auto u = glm::dot(r.origin - P0_, P) / det;
	auto v = glm::dot(r.direction, Q) / det;
	if (u < 0.f || u > 1.f || v < 0.f || v > 1.f)
		return false;
	auto t = glm::dot(e2, Q) / det;
	if (t <= eps)
		return false;
	id->t = t;
	id->normal = glm::normalize(glm::cross(e1, e2));
	id->material = material();
	return true;
}

glm::vec3 Plane::getPointOnSurface(float u, float v) const {
	return P0_ + u * (P1_ - P0_) + v * (P2_ - P0_);
}

float Plane::getArea() const {
	return AREA_;
}

glm::vec3 Plane::getNormal() const {
	return NORMAL_;
}

glm::vec3 Plane::getFirstTangent() const {
	return glm::normalize(P1_ - P0_);
}

// --- LightSource class functions --- //

LightSource::LightSource(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, float flux, SpectralDistribution color) : emitter_(p0, p1, p2, NULL), radiosity(flux / emitter_.getArea() * color) {}

bool LightSource::intersect(LightSourceIntersectionData* light_id, Ray r)
{
	IntersectionData id;
	if (!emitter_.intersect(&id, r))
		return false;
	light_id->normal = id.normal;
	light_id->t = id.t;
	light_id->radiosity = radiosity;
	light_id->area = getArea();
	return true;
}

glm::vec3 LightSource::getPointOnSurface(float u, float v) {
	return emitter_.getPointOnSurface(u, v);
}

float LightSource::getArea() const {
	return emitter_.getArea();
}

glm::vec3 LightSource::getNormal() const {
	return emitter_.getNormal();
}

Ray LightSource::shootLightRay() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0, 1);

	Ray ret;
	ret.origin = getPointOnSurface(dis(gen), dis(gen));
	glm::vec3 normal = emitter_.getNormal();
	glm::vec3 tangent = emitter_.getFirstTangent();
	float inclination = acos(dis(gen));
	float azimuth = 2 * M_PI * dis(gen);
	auto rand_dir = normal;
	rand_dir = glm::normalize(glm::rotate(rand_dir, inclination, tangent));
	rand_dir = glm::normalize(glm::rotate(rand_dir, azimuth, normal));
	ret.direction = rand_dir;
	ret.material = Material::air();
	return ret;
}