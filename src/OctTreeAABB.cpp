#include "../include/OctTreeAABB.h"
#include "../include/Object3D.h"

#include "../external_libraries/common_include/boxOverlap.h"

#include <iostream>

// --- AABB --- //

bool AABB::intersect(Ray ray, IntersectionData* id) const {
    glm::vec3 origin = ray.origin;
    glm::vec3 direction = ray.direction;

    float tl, tr;
    
    float l, r;
    float x = origin.x, y = origin.y, z = origin.z;
    float dx = direction.x, dy = direction.y, dz = direction.z;
    
    l = (mn.x - x) / dx; r = (mx.x - x) / dx;
    if(l > r) std::swap(l, r);
    tl = l; tr = r;
    
    l = (mn.y - y) / dy; r = (mx.y - y) / dy;
    if(l > r) std::swap(l, r);
    tl = std::max(tl, l); tr = std::min(tr, r);
    if(tl > tr) return false;
    
    l = (mn.z - z) / dz; r = (mx.z - z) / dz;
    if(l > r) std::swap(l, r);
    tl = std::max(tl, l); tr = std::min(tr, r);
    if(tl > tr) return false;
    
	if(tl > id->t) return false; //the triangle will be in other BB

    return true;
}

bool AABB::intersectTriangle(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) const {
    glm::vec3 center_point = (mn + mx) / 2.f;
    glm::vec3 scale = (mx - center_point) / 1.f;
    float boxcenter[3] = {center_point[0], center_point[1], center_point[2] };
    float boxhalfsize[3] = {scale[0], scale[1], scale[2] };
    float triverts[3][3] = { {p0.x, p0.y, p0.z }, {p1.x, p1.y, p1.z }, {p2.x, p2.y, p2.z } };
    return triBoxOverlap( boxcenter, boxhalfsize, triverts) == 1;
}

// --- OctNode --- //

OctNodeAABB::OctNodeAABB(OctNodeAABB* par, int d, Mesh* _mesh,
        glm::vec3 _aabb_min, glm::vec3 _aabb_max) {
    aabb.mn = _aabb_min;
    aabb.mx = _aabb_max;
    mesh = _mesh;
	
    std::vector<unsigned int>& indices = (!par) ? mesh->indices_ : par->triangle_indices;
    for (int i = 0; i < indices.size(); i += 3) {
        if (aabb.intersectTriangle(mesh->positions_[indices[i]], mesh->positions_[indices[i + 1]], mesh->positions_[indices[i + 2]])) {
            for (int j = 0; j < 3; ++j) triangle_indices.push_back(indices[i + j]);
        }
    }

    if (d == 0 || triangle_indices.size() <= 3 * 16) {
        for (int i = 0; i < 8; ++i) children_[i] = NULL;
		return;
    }
	
    for (int i = 0; i < 8; ++i) {
        glm::vec3 child_mn = glm::vec3(
            i & 1 == 0 ? _aabb_min.x : (_aabb_min.x + _aabb_max.x) / 2.f,
            i & 2 == 0 ? _aabb_min.y : (_aabb_min.y + _aabb_max.y) / 2.f,
            i & 4 == 0 ? _aabb_min.z : (_aabb_min.z + _aabb_max.z) / 2.f);
        glm::vec3 child_mx = glm::vec3(
            i & 1 == 0 ? (_aabb_min.x + _aabb_max.x) / 2.f : _aabb_max.x,
            i & 2 == 0 ? (_aabb_min.y + _aabb_max.y) / 2.f : _aabb_max.y,
            i & 4 == 0 ? (_aabb_min.z + _aabb_max.z) / 2.f : _aabb_max.z);
        children_[i] = new OctNodeAABB(this, d - 1, mesh, child_mn, child_mx);
    }
}

OctNodeAABB::~OctNodeAABB() {
    for (int i = 0; i < 8; ++i) if (children_[i]) delete children_[i];
}

bool OctNodeAABB::intersect(IntersectionData* id, Ray r) const {
    if (triangle_indices.size() == 0) return false;
    float t_smallest = id? id->t : 1e9;
    bool fl = false;
    if (children_[0] == NULL) {//leaf
        // Check intersection for all triangles in this node
        for (int i = 0; i < triangle_indices.size(); i += 3) {
            // Moller-Trumbore intersection algorithm
            glm::vec3 p0 = mesh->positions_[triangle_indices[i]];
            glm::vec3 p1 = mesh->positions_[triangle_indices[i + 1]];
            glm::vec3 p2 = mesh->positions_[triangle_indices[i + 2]];

            glm::vec3 e1 = p1 - p0;
			glm::vec3 e2 = p2 - p0;
            glm::vec3 P = glm::cross(r.direction, e2);
            float det = glm::dot(e1, P);
            if(det > -0.00001 && det < 0.00001) continue;

            glm::vec3 T = r.origin - p0;
            float u = glm::dot(T, P) / det;
            if(u < 0.f || u > 1.f) continue;

            glm::vec3 Q = glm::cross(T, e1);
            float v = glm::dot(r.direction, Q) / det;
            if(v < 0.f || u + v > 1.f) continue;

           float t = glm::dot(e2, Q) / det;
            if(t > 0.00001 && t < t_smallest) { //ray intersection
                t_smallest = t;
                glm::vec3 n0 = mesh->normals_[triangle_indices[i]];
                glm::vec3 n1 = mesh->normals_[triangle_indices[i + 1]];
                glm::vec3 n2 = mesh->normals_[triangle_indices[i + 2]];
				id->t = t;
				id->normal = glm::normalize((1 - u - v) * n0 + u * n1 + v * n2);
				id->material = mesh->material();
                fl = true;
            }
        }

        return fl;
    }

    for (int i = 0; i < 8; ++i) if(children_[i]->aabb.intersect(r, id))
        fl |= children_[i]->intersect(id, r);
    
	return fl;
}

// --- OctTreeAABB --- //

OctTreeAABB::OctTreeAABB(Mesh* mesh) : OctNodeAABB(NULL, 8, mesh,
        mesh->getMinPosition(), mesh->getMaxPosition()) {}

OctTreeAABB::~OctTreeAABB() {}