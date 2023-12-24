#include "Hittable.h"

void Sphere::Sample(HitRecord *hit_record, float &pdf) const
{
    // TODO: Add your code here.
    float randomNumber1 = get_random_float();
    float randomNumber2 = get_random_float();
    float theta, phi;

    theta = 2 * M_PI * randomNumber1;
    phi = M_PI / 2 * randomNumber2;

    Vec dir(sin(phi) * cos(theta), cos(phi), sin(phi) * sin(theta));
    hit_record->position = this->o_ + this->r_ * dir;
    hit_record->normal = dir;
    hit_record->material = material_;
    hit_record->emission = emission;
    pdf = 1;
}

void Quadric::Sample(HitRecord *hit_record, float &pdf) const
{
    // No need to implement this function
}

void Triangle::Sample(HitRecord *hit_record, float &pdf) const
{
    // TODO: Add your code here.
    float randomNumber1 = get_random_float();
    float randomNumber2 = get_random_float();

    Point p = a_ + randomNumber1 * (b_ - a_) + randomNumber2 * (c_ - a_);

    hit_record->position = p;
    hit_record->normal = n_a_;
    pdf = 1;
}

void CompleteTriangle::Sample(HitRecord *hit_record, float &pdf) const
{
    triangle_.Sample(hit_record, pdf);
    hit_record->material = material_;
    hit_record->emission = emission;
}

void Mesh::Sample(HitRecord *hit_record, float &pdf) const
{
    // TODO: Add your code here.
}

void HittableList::Sample(HitRecord *hit_record, float &pdf) const
{
    // TODO: Add your code here.
}

// Sphere
bool Sphere::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    Point oc = ray.o - this->o_;
    double A, B, C;
    A = 1;
    B = 2.0f * glm::dot(oc, ray.d);
    C = glm::dot(oc, oc) - this->r_ * this->r_;

    double discrim = (B * B) - (4 * A * C);
    if (discrim >= 0)
    {
        double roots = sqrt(discrim);
        float t = fmin((-B - roots) / (2.0f * A), (-B + roots) / (2.0f * A));
        if (t > 0)
        {
            hit_record->position = ray.At(t);
            hit_record->in_direction = ray.d;
            hit_record->normal = (hit_record->position - this->o_) / this->r_;
            hit_record->material = this->material_;
            hit_record->distance = glm::length(hit_record->position - ray.o);
            hit_record->reflection = glm::normalize(ray.d - 2 * (glm::dot(hit_record->normal, ray.d)) * hit_record->normal);
            return true;
        }
    }

    return false;
}

// Quadric
bool Quadric::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    glm::vec4 O(ray.o.x, ray.o.y, ray.o.z, 1);
    glm::vec4 D(ray.d.x, ray.d.y, ray.d.z, 0);
    float A, B, C;
    A = glm::dot(D, this->A_ * D);
    B = 2.0f * glm::dot(O, this->A_ * D);
    C = glm::dot(O, this->A_ * O);

    float discrim = (B * B) - (4 * A * C);
    if (discrim >= 0)
    {
        float roots = sqrt(discrim);
        float t = fmin((-B - roots) / (2.0f * A), (-B + roots) / (2.0f * A));
        glm::mat4x4 gradient_mat = this->A_ + glm::transpose(this->A_);
        if (t > 0)
        {
            hit_record->position = ray.At(t);
            hit_record->in_direction = ray.d;
            hit_record->material = this->material_;
            hit_record->distance = glm::length(hit_record->position - ray.o);

            glm::vec4 rayPos(hit_record->position, 1);
            hit_record->normal = glm::normalize(Vec(gradient_mat * rayPos));
            hit_record->reflection = glm::normalize(ray.d - 2 * (glm::dot(hit_record->normal, ray.d)) * hit_record->normal);
            return true;
        }
    }

    return false;
}

// Triangle
bool Triangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    // TODO: Add your code here.
    Vec norm = glm::normalize(glm::cross((this->b_ - this->a_), (this->c_ - this->a_)));

    float nd = glm::dot(norm, ray.d);
    float np = glm::dot(norm, ray.o);
    float d = glm::dot(norm, this->a_);
    float t = (d - np) / nd;

    if (t < 0.0001)
    {
        return false;
    }

    Point O = ray.At(t);

    Vec AB_test = glm::cross((O - this->a_), (O - this->b_));
    Vec BC_test = glm::cross((O - this->b_), (O - this->c_));
    Vec CA_test = glm::cross((O - this->c_), (O - this->a_));

    if (glm::dot(AB_test, BC_test) > 0 &&
        glm::dot(BC_test, CA_test) > 0 &&
        glm::dot(CA_test, AB_test) > 0)
    {
        hit_record->position = O;
        hit_record->in_direction = ray.d;
        hit_record->distance = glm::length(O - ray.o);
        if (phong_interpolation_)
        {
            float denom = glm::dot(glm::cross((this->b_ - this->a_), (this->c_ - this->a_)), norm);
            float alpha = glm::dot(glm::cross((this->c_ - this->b_), (O - this->b_)), this->n_b_) / denom;
            float beta = glm::dot(glm::cross((this->a_ - this->c_), (O - this->c_)), this->n_c_) / denom;
            float gamma = glm::dot(glm::cross((this->b_ - this->a_), (O - this->a_)), this->n_a_) / denom;

            hit_record->normal = glm::normalize(alpha * this->n_a_ + beta * this->n_b_ + gamma * this->n_c_);
            hit_record->reflection = glm::normalize(ray.d - 2 * (glm::dot(hit_record->normal, ray.d)) * hit_record->normal);
        }
        else
        {
            hit_record->normal = this->n_a_;
            hit_record->reflection = glm::normalize(ray.d - 2 * (glm::dot(hit_record->normal, ray.d)) * hit_record->normal);
        }
        // no need to set material in this function
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray &ray, HitRecord *hit_record) const
{
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret)
    {
        hit_record->material = material_;
    }
    return ret;
}

// Mesh
Mesh::Mesh(const std::string &file_path,
           const Material &material,
           bool phong_interpolation) : ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation)
{
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++)
    {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto &face : f_ind_)
    {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto &vertex_normal : vertex_normals_)
    {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto &face : f_ind_)
    {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min(1e5f, 1e5f, 1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto &vertex : vertices_)
    {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++)
    {
        InsertFace(root_, i);
    }

    area = 0.0f;
    for (auto &triangle : triangles_)
    {
        area += triangle.getArea();
    }
    emission = material_.emission;
}

bool Mesh::Hit(const Ray &ray, HitRecord *hit_record) const
{
    const bool brute_force = false;
    if (brute_force)
    {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_)
        {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record))
            {
                if (curr_hit_record.distance < min_dist)
                {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f)
        {
            hit_record->material = material_;
            return true;
        }
        return false;
    }
    else
    {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret)
        {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t> &face, const Point &bbox_min, const Point &bbox_max) const
{
    for (size_t idx : face)
    {
        const auto &pt = vertices_[idx];
        for (int i = 0; i < 3; i++)
        {
            if (pt[i] < bbox_min[i] + 1e-6f)
                return false;
            if (pt[i] > bbox_max[i] - 1e-6f)
                return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray &ray, const Point &bbox_min, const Point &bbox_max) const
{
    float t_min = -1e5f;
    float t_max = 1e5f;

    for (int i = 0; i < 3; i++)
    {
        if (glm::abs(ray.d[i]) < 1e-6f)
        {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f)
            {
                t_min = 1e5f;
                t_max = -1e5f;
            }
        }
        else
        {
            if (ray.d[i] > 0.f)
            {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else
            {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode *u, size_t face_idx)
{
    const Point &bbox_min = u->bbox_min;
    const Point &bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++)
    {
        for (size_t b = 0; b < 2; b++)
        {
            for (size_t c = 0; c < 2; c++)
            {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max))
                {
                    if (u->childs[child_idx] == nullptr)
                    {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode *child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs)
    {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode *u, const Ray &ray, HitRecord *hit_record) const
{
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max))
    {
        return false;
    }
    float distance = 1e5f;
    for (const auto &face_idx : u->face_index)
    {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto &child : u->childs)
    {
        if (child == nullptr)
        {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < distance)
            {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}

// Hittable list
void HittableList::PushHittable(const Hittable &hittable)
{
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray &ray, HitRecord *hit_record) const
{
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_)
    {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record))
        {
            if (curr_hit_record.distance < min_dist)
            {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}
