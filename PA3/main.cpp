#include <iostream>
#include <vector>
#include <thread>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray &ray, const std::vector<LightSource> &light_sources, const Hittable &scene, int trace_depth);
Color Shade(const std::vector<LightSource> &light_sources,
            const Hittable &hittable_collection,
            const HitRecord &hit_record,
            int trace_depth);

void RenderImagePart(int tid, float *imageBuffer, int width, int height, int startLine, int windowSize, const Camera *camera, const Scene *scene, int spp)
{
    int endLine = startLine + windowSize;
    for (int i = startLine; i < endLine; i++)
    {
        int x = i % width;
        int y = height - i / width - 1;

        Color color(0.f, 0.f, 0.f);
        for (int j = 0; j < spp; j++)
        {
            float biasX = get_random_float();
            float biasY = get_random_float();
            Ray ray = camera->RayAt(float(x) + biasX, float(y) + biasY);
            color += TraceRay(ray, scene->light_sources_, scene->hittable_collection_, 1);
        }
        color /= float(spp);
        color = Color(color.x > 1.0 ? 1.0 : color.x,
                      color.y > 1.0 ? 1.0 : color.y,
                      color.z > 1.0 ? 1.0 : color.z);

        int idx = 3 * i;
        imageBuffer[idx] += color.r;
        imageBuffer[idx + 1] += color.g;
        imageBuffer[idx + 2] += color.b;
    }
}

int main()
{
    // TODO: Set your workdir (absolute path) here. Don't forget the last slash
    const std::string work_dir("/Users/jayjung/Dev/COMP3271-Computer-Graphics/PA3/");

    // Construct scene
    Scene scene(work_dir, "scene/teapot_area_light.toml");
    const Camera &camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;
    int pixels = width * height;

    std::vector<unsigned char> image(width * height * 4, 0);
    float *image_buffer = new float[width * height * 3];
    for (int i = 0; i < width * height * 3; ++i)
    {
        image_buffer[i] = 0.f;
    }

    int spp = 100;
    int NUM_THREADS = 32;

    float progress = 0.f;

    int windowSize = pixels / NUM_THREADS;
    int wsRemainder = pixels % NUM_THREADS;
    int startLine = 0;

    std::thread threads[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++)
    {
        if (i == 0)
        {
            threads[i] = std::thread(&RenderImagePart, i, image_buffer, width, height, startLine, windowSize + wsRemainder, &camera, &scene, spp);
        }
        else
        {
            threads[i] = std::thread(&RenderImagePart, i, image_buffer, width, height, startLine, windowSize, &camera, &scene, spp);
        }
        startLine += windowSize;
    }

    float currProgress = 0.f;
    for (int i = 0; i < NUM_THREADS; i++)
    {
        threads[i].join();
    }

    // copy from image_buffer to image
    for (int i = 0; i < width * height; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            image[4 * i + j] = (uint8_t)(glm::min(image_buffer[3 * i + j], 1.f - 1e-5f) * 256.f);
        }
        image[4 * i + 3] = 255;
    }

    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);

    lodepng::save_file(png, work_dir + "outputs/teapot_area_light.png");
}

Vec ToWorld(const Vec &localRay, const Vec &N)
{
    // This function transforms a vector from local coordinate to world coordinate.
    Vec perpendicularVec;

    if (abs(N.x) > abs(N.y))
    {
        perpendicularVec = glm::normalize(Vec(N.z, 0, -N.x));
    }
    else
    {
        perpendicularVec = glm::normalize(Vec(0, -N.z, N.y));
    }
    Vec orthogonalVec = glm::normalize(glm::cross(perpendicularVec, N));
    glm::mat3 transformMatrix = glm::mat3(perpendicularVec, N, orthogonalVec);

    return transformMatrix * localRay;
}

Vec SampleHemisphere(const HitRecord &hit_record)
{
    // This function randomly samples a direction on the hemisphere.
    // It will calls ToWorld() to transform the local coordinate to world coordinate.
    float randomNumber1 = get_random_float();
    float randomNumber2 = get_random_float();

    float sphericalTheta, sphericalPhi;
    sphericalTheta = 2 * M_PI * randomNumber1;
    sphericalPhi = M_PI / 2 * randomNumber2;

    Vec localDirection(sin(sphericalPhi) * cos(sphericalTheta), cos(sphericalPhi), sin(sphericalPhi) * sin(sphericalTheta));
    return ToWorld(localDirection, hit_record.normal);
}

Color Shade(const std::vector<LightSource> &light_sources,
            const Hittable &hittable_collection,
            const HitRecord &hit_record,
            int trace_depth)
{
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    if (trace_depth >= kMaxTraceDepth)
    {
        return color;
    }

    float emissionRay = hit_record.material.emission;
    float pdf = 1 / M_PI;
    Ray outRay(hit_record.position, SampleHemisphere(hit_record));
    float adjointRay = fmax(glm::dot(outRay.d, hit_record.normal), 0.0f);
    Color brdf = hit_record.material.k_d / M_PI * hit_record.material.diffuse;
    Color incomingRay = TraceRay(outRay, light_sources, hittable_collection, trace_depth + 1);
    color = (emissionRay + incomingRay * adjointRay * brdf / pdf);
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

    if (hittable_collection.Hit(ray, &record))
    {
        Color shadedColor = Shade(light_sources, hittable_collection, record, trace_depth);
        return shadedColor;
    }

    return color;
}