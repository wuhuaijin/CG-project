#include <string>
#include <iostream>
#include "../include/utils.h"

const float Photon::RADIUS = 0.1f;

SpectralDistribution::SpectralDistribution() {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = 0;
}

float SpectralDistribution::norm() const {
	float s = 0;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		s += data[i];
	return s / N_WAVELENGTHS;
}

std::ostream& operator << (std::ostream& os, const SpectralDistribution& t) {
	os << "[ ";
	for (int i = 0; i < t.N_WAVELENGTHS - 1; i++)
		os << t.data[i] << ", ";
	os << t.data[t.N_WAVELENGTHS - 1] << "]";
	return os;
}

SpectralDistribution operator * (float f, const SpectralDistribution& t) {
	return t * f;
}

float& SpectralDistribution::operator [] (const int i) {
	return data[i];
}

SpectralDistribution SpectralDistribution::operator + (const SpectralDistribution& t) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = data[i] + t.data[i];
	return ret;
}

SpectralDistribution SpectralDistribution::operator - (const SpectralDistribution& t) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = data[i] - t.data[i];
	return ret;
}

SpectralDistribution SpectralDistribution::operator ^ (const float& f) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = pow(data[i], f);
	return ret;
}

SpectralDistribution SpectralDistribution::operator / (const float& f) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = data[i] / f;
	return ret;
}

SpectralDistribution SpectralDistribution::operator*(const float& f) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = data[i] * f;
	return ret;
}

SpectralDistribution SpectralDistribution::operator*(const SpectralDistribution& t) const {
	SpectralDistribution ret;
	for (int i = 0; i < N_WAVELENGTHS; i++)
		ret.data[i] = data[i] * t.data[i];
	return ret;
}


SpectralDistribution SpectralDistribution::operator += (const SpectralDistribution& t) {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = data[i] + t.data[i];
	return *this;
}

SpectralDistribution SpectralDistribution::operator -= (const SpectralDistribution& t) {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = data[i] - t.data[i];
	return *this;
}

SpectralDistribution SpectralDistribution::operator *= (const float& f) {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = data[i] * f;
	return *this;
}

SpectralDistribution SpectralDistribution::operator /= (const float& f) {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = data[i] / f;
	return *this;
}

SpectralDistribution SpectralDistribution::operator *= (const SpectralDistribution& t) {
	for (int i = 0; i < N_WAVELENGTHS; i++)
		data[i] = data[i] * t.data[i];
	return *this;
}

Material Material::air() {
	Material air;
	air.transmissivity = air.refraction_index = 1;
	return air;
}

SpectralDistribution evaluateLambertianBRDF(glm::vec3 d1, glm::vec3 d2,	glm::vec3 normal, SpectralDistribution albedo) {
	return albedo / M_PI;
}

SpectralDistribution evaluatePerfectBRDF(SpectralDistribution albedo) {
	return albedo;
}

SpectralDistribution evaluateOrenNayarBRDF(glm::vec3 d1, glm::vec3 d2, glm::vec3 normal, SpectralDistribution albedo, float roughness) {
	float sigma2 = roughness * roughness;
	float A = 1 - 0.5 * sigma2 / (sigma2 + 0.57);
	float B = 0.45 * sigma2 / (sigma2 + 0.09);
	float theta = glm::acos(glm::dot(d2, normal));
	float theta_d1 = glm::acos(glm::dot(d1, normal));
	float alpha = glm::max(theta, theta_d1);
	float beta = glm::min(theta, theta_d1);

	return albedo / M_PI * (A + (B * glm::max(0.f, glm::dot(d1, d2))) * glm::sin(alpha) * glm::tan(beta));
}