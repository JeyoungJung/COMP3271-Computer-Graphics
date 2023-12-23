#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &scene,
               int trace_depth);

Color Shade(const std::vector<LightSource> &light_sources,
            const Hittable &hittable_collection,
            const HitRecord &hit_record,
            int trace_depth)
{
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    // Ambient term
    color = hit_record.material.k_a * hit_record.material.ambient;
    HitRecord shadow_record;

    // Diffuse and Specular term
    // for each light source
    for (LightSource ls : light_sources)
    {
        Ray shadow_ray = Ray(hit_record.position, glm::normalize(ls.position - hit_record.position));
        Vec li = glm::normalize(Vec(ls.position - hit_record.position));
        float nli = glm::dot(hit_record.normal, li);

        // the dot product is more than 0 (less than 90 deg), the light will hit
        if (nli >= 0)
        {
            // if the shadow ray did not hit anything, compute the second part
            if (!hittable_collection.Hit(shadow_ray, &shadow_record))
            {
                Vec shoot_in_reversed = -1.0f * Vec(hit_record.in_direction);
                // rv cannot be negative
                float rv = fmax(0.0f, glm::dot(hit_record.reflection, shoot_in_reversed));
                Color ds = 1 * ls.intensity *
                           (hit_record.material.k_d * hit_record.material.diffuse * nli +
                            hit_record.material.k_s * hit_record.material.specular * pow(rv, hit_record.material.sh));
                color += ds;
            }
        }
    }
    // Reflected ray contribution
    if (trace_depth < kMaxTraceDepth)
    {
        if (hit_record.material.k_s > 0)
        {
            Ray reflected_ray = Ray(hit_record.position, glm::normalize(hit_record.reflection));
            Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth + 1);
            color += r_color * hit_record.material.k_s;
        }
    }
    // clamp color to 1.0
    color = Color(color.x > 1.0 ? 1.0 : color.x,
                  color.y > 1.0 ? 1.0 : color.y,
                  color.z > 1.0 ? 1.0 : color.z);

    return color;
}

Color TraceRay(const Ray &ray,
               const std::vector<LightSource> &light_sources,
               const Hittable &hittable_collection,
               int trace_depth)
{
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);
    // If ray hits, shade the pixel
    if (hittable_collection.Hit(ray, &record))
    {
        color = Shade(light_sources, hittable_collection, record, trace_depth);
    }
    return color;
}

int main()
{
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("/Users/jayjung/Dev/COMP3271-Computer-Graphics/PA2/");

    // Construct scene
    Scene scene(work_dir, "scene/monkey.toml");
    const Camera &camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f)
            {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f)
                {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++)
            {
                image[idx + i] = (uint8_t)(glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f)
            {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}
