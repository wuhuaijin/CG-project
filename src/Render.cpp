#include "../include/Render.h"

#include "../external_libraries/common_include/pugixml.h"
#include "../include/xmlTraverser.h"
#include "../include/Scenes.h"

#include <iostream>
#include <string>
#include <random>

Render::Render(Scene* scene_) {
	scene = scene_;
    GENERATOR = new std::mt19937(rd_());
    DISTRIBUTION = new std::uniform_real_distribution<float>(0, 1);
}

Render::~Render() {
	delete GENERATOR;
	delete DISTRIBUTION;
}

bool Render::intersect(IntersectionData* id, Ray r) {

    //这里如果是球面或者三角形的话可以直接算，但是不确定是否只有球面和三角形
    IntersectionData res_id;
    res_id.t = 1e9;

    // 遍历找到最小交点
    Object3D* res_obj = NULL;
    IntersectionData id_it;
    for (int i = 0; i < scene->object_list.size(); ++i) {
        Object3D* obj_it = scene->object_list[i];
        if (obj_it->intersect(&id_it, r)) {
            if (id_it.t < res_id.t) {
                res_id = id_it;
                res_obj = obj_it;
            }
        }
    }
    if (res_obj != NULL) {
        *id = res_id;
        return true;
    }
    else return false;
}

bool Render::intersectLight(LightSourceIntersectionData* light_id, Ray r) {
	LightSourceIntersectionData light_res_id;
	light_res_id.t = 1e9; 

    //先找光
	LightSource* res_light = NULL;
    LightSourceIntersectionData id_it;
	for (int i = 0; i < scene->light_list.size(); ++i) {
        LightSource* light_it = scene->light_list[i];
		if (light_it->intersect(&id_it,r)) {
            if (id_it.t < light_res_id.t) {
                light_res_id = id_it;
                res_light = light_it;
            }
		}
	}
    if (!res_light) return false;
    IntersectionData res_id;
    res_id.t = 1e9;

    // 遍历找到最小交点
    Object3D* res_obj = NULL;
    IntersectionData id_it_;
    for (int i = 0; i < scene->object_list.size(); ++i) {
        Object3D* obj_it = scene->object_list[i];
        if (obj_it->intersect(&id_it_, r)) {
            if (id_it_.t < res_id.t) {
                res_id = id_it_;
                res_obj = obj_it;
            }
        }
    }

	if (res_obj && res_id.t < light_res_id.t) {
		return false;
	}

	*light_id = light_res_id;
	return true;	
}

SpectralDistribution Render::traceDiffuseRay(Ray r, int render_mode, IntersectionData id, int iteration) {
	r.has_intersected = true;
    return traceLocalDiffuseRay(r, render_mode, id) + traceIndirectDiffuseRay(r, render_mode, id, iteration);
}


SpectralDistribution Render::traceLocalDiffuseRay(Ray r, int render_mode, IntersectionData id) {
	SpectralDistribution L_local;
	static const int num_samples = 1;
	for (int i = 0; i < scene->light_list.size(); ++i) {
		for (int j = 0; j < num_samples; ++j) {
			Ray shadow_ray = r;
			glm::vec3 differance = scene->light_list[i]->getPointOnSurface((*DISTRIBUTION)(*GENERATOR),(*DISTRIBUTION)(*GENERATOR)) - shadow_ray.origin;
			shadow_ray.direction = glm::normalize(differance);

			SpectralDistribution brdf;
			float cos_theta = glm::dot(shadow_ray.direction, id.normal);

			LightSourceIntersectionData shadow_ray_id;

			if (id.material.diffuse_roughness) {
				brdf = evaluateOrenNayarBRDF(-r.direction, shadow_ray.direction, id.normal,
					id.material.color_diffuse * id.material.reflectance * (1 -id.material.specular_reflectance), id.material.diffuse_roughness);
			}
			else {
				brdf = evaluateLambertianBRDF(-r.direction, shadow_ray.direction, id.normal,
					id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
            }
			if (intersectLight(&shadow_ray_id, shadow_ray)) {
				float cos_light_angle = glm::dot(shadow_ray_id.normal, -shadow_ray.direction);
				float light_solid_angle = shadow_ray_id.area / num_samples * glm::clamp(cos_light_angle, 0.0f, 1.0f) / glm::pow(glm::length(differance), 2) / (M_PI * 2);
				L_local += brdf * shadow_ray_id.radiosity * cos_theta * light_solid_angle;
			}
		}
	}
	return L_local;
}

SpectralDistribution Render::traceIndirectDiffuseRay(Ray r, int render_mode, IntersectionData id, int iteration) {
	SpectralDistribution L_indirect;
	static const int num_samples = 1;
	for (int i = 0; i < num_samples; ++i) {
		glm::vec3 tan = glm::normalize(glm::cross(id.normal, id.normal + glm::vec3(1,1,1)));

		float rand1 = (*DISTRIBUTION)(*GENERATOR);
		float rand2 = (*DISTRIBUTION)(*GENERATOR);

		// Uniform distribution over a hemisphere
		float inclination = acos(sqrt(rand1));
		float azimuth = 2 * M_PI * rand2;
		// Change the actual vector
		glm::vec3 random_direction = id.normal;
		random_direction = glm::normalize(glm::rotate(random_direction, inclination, tan));
		random_direction = glm::normalize(glm::rotate(random_direction, azimuth, id.normal));

		float cos_angle = glm::dot(random_direction, id.normal);

		SpectralDistribution brdf;
		if (id.material.diffuse_roughness) {
			brdf = evaluateOrenNayarBRDF(-r.direction, random_direction, id.normal,
				id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance), id.material.diffuse_roughness);
		}
		else {
			brdf = evaluateLambertianBRDF(-r.direction, random_direction, id.normal,
				id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
		}

		r.direction = random_direction;
		r.radiance *= M_PI * brdf; 
		L_indirect += traceRay(r, render_mode, iteration + 1) * M_PI * brdf;
	}
	return L_indirect / num_samples;
}


SpectralDistribution Render::traceRefractedRay(Ray r, int render_mode, IntersectionData id, int iteration, glm::vec3 offset, bool inside) {
	Ray recursive_ray = r;
	recursive_ray.has_intersected = true;

	glm::vec3 normal = inside ? -id.normal : id.normal;
	glm::vec3 perfect_refraction = glm::refract(r.direction, normal, r.material.refraction_index / id.material.refraction_index);
    if (inside)
		offset = -offset;

	if (perfect_refraction != glm::vec3(0)) { 
		float n1 = r.material.refraction_index;
		float n2 = id.material.refraction_index;
		float R_0 = pow((n1 - n2)/(n1 + n2), 2);
		float R = R_0 + (1 - R_0) * pow(1 - glm::dot(normal, -r.direction),5);

		Ray recursive_ray_reflected = recursive_ray;
		Ray recursive_ray_refracted = recursive_ray;

		recursive_ray_reflected.material = Material::air();
		recursive_ray_reflected.origin = r.origin + id.t * r.direction + offset;
		recursive_ray_refracted.material = id.material;
		recursive_ray_refracted.origin = r.origin + id.t * r.direction - offset;
		
		SpectralDistribution to_return;
		recursive_ray_reflected.direction = glm::reflect(r.direction, id.normal);
		recursive_ray_refracted.direction = perfect_refraction;

		SpectralDistribution brdf_specular = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance * R);
		SpectralDistribution brdf_refractive = evaluatePerfectBRDF(id.material.color_diffuse * id.material.reflectance * id.material.specular_reflectance * (1 - R));

		recursive_ray_reflected.radiance *= brdf_specular;
		recursive_ray_refracted.radiance *= brdf_refractive;

		SpectralDistribution reflected_part = traceRay(recursive_ray_reflected, render_mode, iteration + 1) * brdf_specular;
		SpectralDistribution refracted_part = traceRay(recursive_ray_refracted, render_mode, iteration + 1) * brdf_refractive;
		return reflected_part + refracted_part;
	}
	else {
        recursive_ray.origin = r.origin + id.t * r.direction + offset;
		SpectralDistribution brdf_specular = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance);
		recursive_ray.direction = glm::reflect(r.direction, id.normal);
		recursive_ray.radiance *= brdf_specular;
		return traceRay(recursive_ray, render_mode, iteration + 1) * brdf_specular;
	}
}


SpectralDistribution Render::traceSpecularRay(Ray r, int render_mode, IntersectionData id, int iteration) {
	r.has_intersected = true;
	SpectralDistribution specular = SpectralDistribution();
	r.direction = glm::reflect(r.direction, id.normal);
	SpectralDistribution brdf = evaluatePerfectBRDF(id.material.color_specular * id.material.reflectance * id.material.specular_reflectance);
	r.radiance *= brdf;
	specular += traceRay(r, render_mode, iteration + 1) * brdf;
	return specular;
}


SpectralDistribution Render::traceRay(Ray r, int render_mode, int iteration) {
	IntersectionData id;
	LightSourceIntersectionData light_id;

	if (intersectLight(&light_id, r)) {
        if (render_mode == WHITTED_SPECULAR) {
            return light_id.radiosity / (M_PI * 2);
        }
        else {
            return SpectralDistribution();
        }
    }
	else if (intersect(&id, r))	{
		float random = (*DISTRIBUTION)(*GENERATOR);
		float non_termination_probability = iteration == 0 ? 1.0 : 0.7;
        /// 可调整
		if (random > non_termination_probability || iteration > 20)
			return SpectralDistribution();
		glm::vec3 offset = id.normal * 0.00001f;
		bool inside = false;
		if (glm::dot(id.normal, r.direction) > 0) 
			inside = true;
		
		float transmissivity = id.material.transmissivity;
		float reflectance = id.material.reflectance;
		float specularity = id.material.specular_reflectance;

		SpectralDistribution total;
		if (1 - transmissivity) { // Completely or partly reflected
			Ray recursive_ray = r;
			recursive_ray.origin = r.origin + id.t * r.direction + (inside ? -offset : offset);
			SpectralDistribution specular_part = specularity ? traceSpecularRay(recursive_ray, render_mode, id, iteration) : SpectralDistribution();
			SpectralDistribution diffuse_part;
			switch (render_mode) {
				case CAUSTICS : {
					KDTreeNode ref_node;
					ref_node.p.position = r.origin + r.direction * id.t + offset;

					std::vector<KDTreeNode> closest_photons;
					photon_map.find_within_range(ref_node, Photon::RADIUS, std::back_insert_iterator<std::vector<KDTreeNode> >(closest_photons));
					
					SpectralDistribution photon_radiance;
					for (int i = 0; i < closest_photons.size(); ++i) {
						SpectralDistribution brdf;
						if (id.material.diffuse_roughness) {
							brdf = evaluateOrenNayarBRDF(-r.direction, closest_photons[i].p.direction_in, id.normal,
								id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance), id.material.diffuse_roughness);
						}
						else {
							brdf = evaluateLambertianBRDF(-r.direction, closest_photons[i].p.direction_in, id.normal,
								id.material.color_diffuse * id.material.reflectance * (1 - id.material.specular_reflectance));
						}
						float distance = glm::length(closest_photons[i].p.position - ref_node.p.position);
						float photon_area = Photon::RADIUS * Photon::RADIUS * M_PI;
						photon_radiance += closest_photons[i].p.delta_flux * (glm::length(distance) < Photon::RADIUS ? 1 : 0) / (photon_area * 2 * M_PI) * brdf * (2 * M_PI);
					}

					diffuse_part = photon_radiance;
					break;
				}
				case WHITTED_SPECULAR : {
					diffuse_part = SpectralDistribution();
					break;
				}
				case MONTE_CARLO : {
					diffuse_part = (1 - specularity) ? traceDiffuseRay(recursive_ray, render_mode, id, iteration) : SpectralDistribution();
					break;
				}
                case PHOTON_MAPPING : {
					if (r.has_intersected) {
						Photon p;
                        p.direction_in = -r.direction;
						p.position = recursive_ray.origin;
						float photon_area = Photon::RADIUS * Photon::RADIUS * M_PI;
						p.delta_flux = recursive_ray.radiance / non_termination_probability * photon_area * M_PI;
						KDTreeNode insert_node;
						insert_node.p = p;
						photon_map.insert(insert_node);
					}
					break;
				}
				default : {
					diffuse_part = SpectralDistribution();
					break;
				}
			}

			total += (specular_part + diffuse_part) * (1 - transmissivity);
		}
		if (transmissivity) { 
			SpectralDistribution transmitted_part = traceRefractedRay(r, render_mode, id, iteration, offset, inside);
			total += transmitted_part * transmissivity;
		}
		return total / non_termination_probability;
	}
	return SpectralDistribution();
}


void Render::buildPhotonMap(const int n_photons) {
    if (scene->light_list.size() == 0) {
        std::cout << "No light " << photon_map.size() << std::endl;
        return;
    }

    SpectralDistribution total_flux = SpectralDistribution();
    float total_flux_res = 0;
    for (int i = 0; i < scene->light_list.size(); ++i) {
        total_flux_res += scene->light_list[i]->radiosity.norm() * scene->light_list[i]->getArea();
        total_flux += scene->light_list[i]->radiosity * scene->light_list[i]->getArea();
    }
    for (int k = 0; k < 50; ++k) {
        #pragma omp parallel for
        ///并行50 可调整
        for (int i = 0; i < n_photons / 50; ++i) {
            int picked_light_source = 0;
            float accumulating_chance = 0;
            float random = (*DISTRIBUTION)(*GENERATOR);
            for (int i = 0; i < scene->light_list.size(); ++i) {
                float interval = scene->light_list[i]->radiosity.norm() * scene->light_list[i]->getArea() / total_flux_res;
                if (random > accumulating_chance && random < accumulating_chance + interval) { 
                    picked_light_source = i;
                    break;
                }
                else {
                    accumulating_chance += interval; 
                }
            }
            Ray r = scene->light_list[picked_light_source]->shootLightRay();
            r.has_intersected = false;
            SpectralDistribution delta_flux = total_flux / n_photons;
            float photon_area = Photon::RADIUS * Photon::RADIUS * M_PI;
            r.radiance = delta_flux / (photon_area * M_PI * 2);
            traceRay(r, PHOTON_MAPPING);
        }
    }
    std::cout << "Number of photons in scene: " << photon_map.size() << std::endl;
    std::cout << "Optimizing kd tree" << std::endl;
    photon_map.optimize();

}

