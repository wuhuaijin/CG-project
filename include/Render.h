#include <vector>
#include <map>
#include <random>

#include <glm/glm.hpp>

#include "utils.h"
#include "Object3D.h"
#include "Scenes.h"
#include "../external_libraries/common_include/kdtree++/kdtree.hpp"


class Render {

private:
    friend struct scene_traverser;
    Scene* scene;
    KDTree::KDTree<3, KDTreeNode> photon_map;
    std::random_device rd_;

    std::mt19937* GENERATOR;
    std::uniform_real_distribution<float>* DISTRIBUTION;
    
    SpectralDistribution traceDiffuseRay(Ray r, int render_mode, IntersectionData id, int iteration);
	SpectralDistribution traceLocalDiffuseRay(Ray r, int render_mode, IntersectionData id);
	SpectralDistribution traceIndirectDiffuseRay(Ray r, int render_mode, IntersectionData id, int iteration);
	
	SpectralDistribution traceSpecularRay(Ray r, int render_mode, IntersectionData id, int iteration);
	SpectralDistribution traceRefractedRay(Ray r, int render_mode, IntersectionData id, int iteration, glm::vec3 offset, bool inside);


	bool intersect(IntersectionData* id, Ray r);
    bool intersectLight(LightSourceIntersectionData* light_id, Ray r);

public:
    Render(Scene* scene);
    ~Render();

    enum RenderMode{
	  MONTE_CARLO, PHOTON_MAPPING, CAUSTICS, WHITTED_SPECULAR,
	};

    SpectralDistribution traceRay(Ray r, int render_mode, int iteration = 0);
	void buildPhotonMap(const int n_photons);


};