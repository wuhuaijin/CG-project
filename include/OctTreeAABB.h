#ifndef OCT_TREE_AABB
#define OCT_TREE_AABB

#include <vector>

#include <glm/glm.hpp>
#include "utils.h"

class Mesh;

struct AABB {
	bool intersect(Ray ray, IntersectionData* id) const;
	bool intersectTriangle(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) const;

	glm::vec3 mn, mx;
};

// A node of an octree. Has eight child nodes
class OctNodeAABB {
public:
	OctNodeAABB(OctNodeAABB* par, int d, Mesh* _mesh, glm::vec3 _aabb_min, glm::vec3 _aabb_max);
	~OctNodeAABB();

	bool intersect(IntersectionData* id, Ray r) const;

protected:
	Mesh* mesh;
	AABB aabb;
	std::vector<unsigned int> triangle_indices;

	OctNodeAABB* children_[8]; 
};

// An octree containing axis aligned bounding boxes
class OctTreeAABB : public OctNodeAABB{
public:
	OctTreeAABB(Mesh* mesh);
	~OctTreeAABB();
};

#endif