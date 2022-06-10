#ifndef BVH_TREE
#define BVH_TREE

#include <vector>

#include <glm/glm.hpp>
#include "utils.h"

class Mesh;

struct AABB {
	glm::vec3 mn, mx;
	AABB();
	AABB(Mesh* mesh, unsigned int a, unsigned int b, unsigned int c);
	float intersect(Ray ray, float t_smallest) const;
	void enclose(glm::vec3 p);
	void enclose(AABB tmp);
	bool empty();
	float area();
	glm::vec3 center();
};

// A node of an octree. Has eight child nodes
class BVHNode {
public:
	void build(int d, Mesh* _mesh, glm::vec3 _aabb_min, glm::vec3 _aabb_max);
	void _build(int d, Mesh* _mesh);
	~BVHNode();

	bool intersect(IntersectionData* id, Ray r) const;

protected:
	Mesh* mesh;
	AABB aabb;
	std::vector<unsigned int> tr_indices;
	BVHNode* ch[2]; 
};

// An octree containing axis aligned bounding boxes
class BVHTree : public BVHNode{
public:
	BVHTree(Mesh* mesh);
	~BVHTree();
};

#endif