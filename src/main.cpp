#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include <iostream>
#include "camera.h"
#include "color.h"
#include "vec3.h"
#include "ray.h"

double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto half_b = dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return -1.0;
    return (-half_b - sqrt(discriminant) ) / a;
}

color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    if (depth <= 0) 
        return color(0, 0, 0);
    if(world.hit(r, 0.001, infinity, rec)) {
        point3 target = rec.p + rec.normal + random_in_hemisphere(rec.normal);
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);

}

int main() {
    // Image
    const auto aspect_ratio = 16.0 / 9.0;

    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int spp = 100;
    const int max_depth = 50;

    // World
    hittable_list world;
    world.add(make_shared<sphere>(point3(0,0,-1), 0.5));
    world.add(make_shared<sphere>(point3(0,-100.5,-1), 100));


    
    // Camera
    
    camera cam;
    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < spp; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, spp);
        }
    }

    std::cerr << "\nDone.\n";
}

