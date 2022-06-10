#include "../include/BVHTree.h"
#include "../include/Object3D.h"

#include "../external_libraries/common_include/boxOverlap.h"

#include <iostream>

// --- AABB --- //

AABB::AABB() {
	float INF = 1e9;
	mn = glm::vec3(INF, INF, INF); 
	mx = glm::vec3(-INF, -INF, -INF); 
}

AABB::AABB(Mesh* mesh, unsigned int a, unsigned int b, unsigned int c) {
	float INF = 1e9;
	mn = glm::vec3(INF, INF, INF); 
	mx = glm::vec3(-INF, -INF, -INF); 
	enclose(mesh->positions_[a]);
	enclose(mesh->positions_[b]);
	enclose(mesh->positions_[c]);
}

float AABB::intersect(Ray ray, float t_smallest) const {
    glm::vec3 origin = ray.origin;
    glm::vec3 direction = ray.direction;

    float tl = 0, tr = 1e9;
    
    float l, r;
    float x = origin.x, y = origin.y, z = origin.z;
    float dx = direction.x, dy = direction.y, dz = direction.z;
    
    l = (mn.x - x) / dx; r = (mx.x - x) / dx;
    if(l > r) std::swap(l, r);
    tl = std::max(tl, l); tr = std::min(tr, r);
    if(tl > tr) return 1e9;
    
    l = (mn.y - y) / dy; r = (mx.y - y) / dy;
    if(l > r) std::swap(l, r);
    tl = std::max(tl, l); tr = std::min(tr, r);
    if(tl > tr) return 1e9;
    
    l = (mn.z - z) / dz; r = (mx.z - z) / dz;
    if(l > r) std::swap(l, r);
    tl = std::max(tl, l); tr = std::min(tr, r);
    if(tl > tr) return 1e9;
	
    if(tl < t_smallest) return tl;
	return 1e9;
}

void AABB::enclose(glm::vec3 p) { //enclose a point
	mn.x = std::min(mn.x, p.x);
	mn.y = std::min(mn.y, p.y);
	mn.z = std::min(mn.z, p.z);
	
	mx.x = std::max(mx.x, p.x);
	mx.y = std::max(mx.y, p.y);
	mx.z = std::max(mx.z, p.z);
}

void AABB::enclose(AABB p) { //enclose a bounding box
	mn.x = std::min(mn.x, p.mn.x);
	mn.y = std::min(mn.y, p.mn.y);
	mn.z = std::min(mn.z, p.mn.z);
	
	mx.x = std::max(mx.x, p.mx.x);
	mx.y = std::max(mx.y, p.mx.y);
	mx.z = std::max(mx.z, p.mx.z);
}

bool AABB::empty() {
	return mn.x > mx.x || mn.y > mx.y || mn.z > mx.z;
}

float AABB::area() { // surface area
    if(empty()) return 0.0f;
    glm::vec3 extent = mx - mn;
    return 2.0f * (extent.x * extent.z + extent.x * extent.y + extent.y * extent.z);
}

glm::vec3 AABB::center() { return (mx + mn) / 2.f; }


// --- BVHNode --- //

void BVHNode::build(int d, Mesh* _mesh, glm::vec3 _aabb_min, glm::vec3 _aabb_max) {
	tr_indices.assign(_mesh->indices_.begin(), _mesh->indices_.end());
	//aabb.mn = _aabb_min; aabb.mx = _aabb_max;
	_build(d, _mesh);
}

void BVHNode::_build(int d, Mesh* _mesh) {
    mesh = _mesh;

	int n_tr = tr_indices.size() / 3;
	AABB* tr_bb = new AABB[n_tr];
	glm::vec3* tr_c = new glm::vec3 [n_tr];
	for (int i = 0; i < n_tr; ++i) {
		tr_bb[i] = AABB(_mesh, tr_indices[i * 3], tr_indices[i * 3 + 1], tr_indices[i * 3 + 2]);
		tr_c[i] = tr_bb[i].center();
		aabb.enclose(tr_bb[i]);
	}

    if (d == 0 || n_tr <= 16) {
        ch[0] = ch[1] = NULL;
		delete []tr_bb;
		delete []tr_c;
		return;
    }

	bool fl = 0;
	float cost = 1e9;
	const int N_SPLIT = 32;
	std::vector<unsigned int> lprs, rprs;
	
	for (int o = 0; o < 3; ++o) {
		std::vector<unsigned int>* p = new std::vector<unsigned int> [N_SPLIT];
		AABB* bucket_bb = new AABB[N_SPLIT];

		for (int i = 0; i < n_tr; ++i) {
			int tmp;
			if(o == 0) tmp = (tr_c[i].x - aabb.mn.x) * N_SPLIT / (aabb.mx.x - aabb.mn.x);
			if(o == 1) tmp = (tr_c[i].y - aabb.mn.y) * N_SPLIT / (aabb.mx.y - aabb.mn.y);
			if(o == 2) tmp = (tr_c[i].z - aabb.mn.z) * N_SPLIT / (aabb.mx.z - aabb.mn.z);
			tmp = std::max(tmp, 0);
			tmp = std::min(tmp, N_SPLIT - 1);
			p[tmp].push_back(i);
		}

		int nowans = -1;
		for (int i = 0; i + 1 < N_SPLIT; ++i) {
			AABB nowl, nowr;
			int cntl = 0, cntr = 0;
			for (int j = 0; j <= i; ++j) {
				nowl.enclose(bucket_bb[j]); cntl += p[j].size();
			}
			for (int j = i + 1; j < N_SPLIT; ++j) {
				nowr.enclose(bucket_bb[j]); cntr += p[j].size();
			}

			if(cntl && cntr) {
				float c = nowl.area() * cntl + nowr.area() * cntr;
				if(c < cost) { cost = c; nowans = i; }
			}
		}

		if(nowans > -1) {
			fl = true; lprs.clear(); rprs.clear(); 
			for (int i = 0; i <= nowans; ++i) lprs.insert(lprs.end(), p[i].begin(), p[i].end());
			for (int i = nowans + 1; i < N_SPLIT; ++i) rprs.insert(rprs.end(), p[i].begin(), p[i].end());
		}

		delete []p;
		delete []bucket_bb;
	}
	
	ch[0] = new BVHNode();
	ch[1] = new BVHNode();
	
	if(!fl) {
		int d = n_tr / 2 - 1;
		for (int i = 0; i <= 3 * d; ++i) ch[0]->tr_indices.push_back(tr_indices[i]);
		for (int i = 3 * d + 1; i < 3 * n_tr; ++i) ch[1]->tr_indices.push_back(tr_indices[i]);
	}
	else {
		for (int i = 0; i < lprs.size(); ++i) 
			for (int j = 0; j < 3; ++j) ch[0]->tr_indices.push_back(tr_indices[lprs[i] * 3 + j]);
		for (int i = 0; i < rprs.size(); ++i) 
			for (int j = 0; j < 3; ++j) ch[1]->tr_indices.push_back(tr_indices[rprs[i] * 3 + j]);
		lprs.clear(); rprs.clear();
	}
	
	delete []tr_bb;
	delete []tr_c;

	ch[0]->_build(d - 1, mesh);
	ch[1]->_build(d - 1, mesh);
}

BVHNode::~BVHNode() {
    for (int i = 0; i < 2; ++i) if (ch[i]) delete ch[i];
}

bool BVHNode::intersect(IntersectionData* id, Ray r) const {	
    if (tr_indices.size() == 0) return false;
    bool fl = false;
    if (ch[0] == NULL || ch[1] == NULL) {//leaf
        // Check intersection for all triangles in this node
        for (int i = 0; i < tr_indices.size(); i += 3) {
            // Moller-Trumbore intersection algorithm
            glm::vec3 p0 = mesh->positions_[tr_indices[i]];
            glm::vec3 p1 = mesh->positions_[tr_indices[i + 1]];
            glm::vec3 p2 = mesh->positions_[tr_indices[i + 2]];

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
            if(t > 0.00001 && t < id->t) { //ray intersection
                glm::vec3 n0 = mesh->normals_[tr_indices[i]];
                glm::vec3 n1 = mesh->normals_[tr_indices[i + 1]];
                glm::vec3 n2 = mesh->normals_[tr_indices[i + 2]];
				id->t = t;
				id->normal = glm::normalize((1 - u - v) * n0 + u * n1 + v * n2);
				id->material = mesh->material();
                fl = true;
            }
        }

        return fl;
    }

	float lt = ch[0]->aabb.intersect(r, id->t);
	float rt = ch[1]->aabb.intersect(r, id->t);

    if(lt < rt && lt < id->t) {
        fl |= ch[0]->intersect(id, r);
		if(rt < id->t) fl |= ch[1]->intersect(id, r);
	}
	else if(rt <= lt && rt < id->t) {
        fl |= ch[1]->intersect(id, r);
		if(lt < id->t) fl |= ch[0]->intersect(id, r);
	}
    
	return fl;
}

// --- BVHTree --- //

BVHTree::BVHTree(Mesh* mesh) {
	build(16, mesh, mesh->getMinPosition(), mesh->getMaxPosition());
}

BVHTree::~BVHTree() {}