#include "SDL.h"
#include "SDL_syswm.h"
#include "as-camera-sdl/as-camera-sdl.h"
#include "as-camera/as-camera-controller.hpp"
#include "as/as-math-ops.hpp"
#include "as/as-view.hpp"
#include "bgfx-imgui/imgui_impl_bgfx.h"
#include "bgfx/bgfx.h"
#include "bgfx/platform.h"
#include "bx/math.h"
#include "bx/timer.h"
#include "file-ops.h"
#include "hsv-rgb.h"
#include "imgui.h"
#include "marching-cubes/marching-cubes.h"
#include "sdl-imgui/imgui_impl_sdl.h"

#include <algorithm>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

struct PosColorVertex
{
    as::vec3_t position;
    uint32_t abgr;
};

struct PosNormalUv
{
    as::vec3_t position;
    as::vec3_t normal;
    // as::vec2_t uv;
};

static PosColorVertex cube_vertices_col[] = {
    {as::vec3_t{-1.0f, 1.0f, 1.0f}, 0xff000000},
    {as::vec3_t{1.0f, 1.0f, 1.0f}, 0xff0000ff},
    {as::vec3_t{-1.0f, -1.0f, 1.0f}, 0xff00ff00},
    {as::vec3_t{1.0f, -1.0f, 1.0f}, 0xff00ffff},
    {as::vec3_t{-1.0f, 1.0f, -1.0f}, 0xffff0000},
    {as::vec3_t{1.0f, 1.0f, -1.0f}, 0xffff00ff},
    {as::vec3_t{-1.0f, -1.0f, -1.0f}, 0xffffff00},
    {as::vec3_t{1.0f, -1.0f, -1.0f}, 0xffffffff},
};

static const uint16_t cube_tri_list_col[] = {
    0, 1, 2, 1, 3, 2, 4, 6, 5, 5, 6, 7, 0, 2, 4, 4, 2, 6,
    1, 5, 3, 5, 7, 3, 0, 4, 1, 4, 5, 1, 2, 3, 6, 6, 3, 7,
};

// clang-format off
static PosNormalUv cube_vertices_norm[] = {
    // far face
    /*0*/{as::vec3_t{-1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// ftl
    /*1*/{as::vec3_t{1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// ftr
    /*2*/{as::vec3_t{-1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// fbl
    /*1*/{as::vec3_t{1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// ftr
    /*3*/{as::vec3_t{1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// fbr
    /*2*/{as::vec3_t{-1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, 0.0f, 1.0f}},// fbl

    // near face
    /*4*/{as::vec3_t{-1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// ntl
    /*6*/{as::vec3_t{-1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// nbl
    /*5*/{as::vec3_t{1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// ntr
    /*5*/{as::vec3_t{1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// ntr
    /*6*/{as::vec3_t{-1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// nbl
    /*7*/{as::vec3_t{1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, 0.0f, -1.0f}},// nbr

    // left face
    /*0*/{as::vec3_t{-1.0f, 1.0f, 1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// ftl
    /*2*/{as::vec3_t{-1.0f, -1.0f, 1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// fbl
    /*4*/{as::vec3_t{-1.0f, 1.0f, -1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// ntl
    /*4*/{as::vec3_t{-1.0f, 1.0f, -1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// ntl
    /*2*/{as::vec3_t{-1.0f, -1.0f, 1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// fbl
    /*6*/{as::vec3_t{-1.0f, -1.0f, -1.0f}, as::vec3_t{-1.0f, 0.0f, 0.0f}},// nbl

    // right face
    /*1*/{as::vec3_t{1.0f, 1.0f, 1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// ftr
    /*5*/{as::vec3_t{1.0f, 1.0f, -1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// ntr
    /*3*/{as::vec3_t{1.0f, -1.0f, 1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// fbr
    /*5*/{as::vec3_t{1.0f, 1.0f, -1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// ntr
    /*7*/{as::vec3_t{1.0f, -1.0f, -1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// nbr
    /*3*/{as::vec3_t{1.0f, -1.0f, 1.0f}, as::vec3_t{1.0f, 0.0f, 0.0f}},// fbr

    // top face
    /*0*/{as::vec3_t{-1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ftl
    /*4*/{as::vec3_t{-1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ntl
    /*1*/{as::vec3_t{1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ftr
    /*4*/{as::vec3_t{-1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ntl
    /*5*/{as::vec3_t{1.0f, 1.0f, -1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ntr
    /*1*/{as::vec3_t{1.0f, 1.0f, 1.0f}, as::vec3_t{0.0f, 1.0f, 0.0f}},// ftr

    // bottom face
    /*2*/{as::vec3_t{-1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// fbl
    /*3*/{as::vec3_t{1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// fbr
    /*6*/{as::vec3_t{-1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// nbl
    /*6*/{as::vec3_t{-1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// nbl
    /*3*/{as::vec3_t{1.0f, -1.0f, 1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// fbr
    /*7*/{as::vec3_t{1.0f, -1.0f, -1.0f}, as::vec3_t{0.0f, -1.0f, 0.0f}},// nbr
};
// clang-format on

static bgfx::ShaderHandle createShader(
    const std::string& shader, const char* name)
{
    const bgfx::Memory* mem = bgfx::copy(shader.data(), shader.size());
    const bgfx::ShaderHandle handle = bgfx::createShader(mem);
    bgfx::setName(handle, name);
    return handle;
}

struct Fps
{
    enum
    {
        MaxSamples = 20
    };

    int64_t samples[MaxSamples] = {};
    int head = 0;
    int tail = MaxSamples - 1;
    bool initialized = false;
};

namespace fps
{

int64_t calculateWindow(Fps& fps, const int64_t now)
{
    if (!fps.initialized && fps.head == fps.tail) {
        fps.initialized = true;
    }

    fps.samples[fps.head] = now;
    fps.head = (fps.head + 1) % fps.MaxSamples;

    int64_t result = -1;
    if (fps.initialized) {
        result = fps.samples[fps.tail] - fps.samples[fps.head];
        fps.tail = (fps.tail + 1) % fps.MaxSamples;
    }

    return result;
}

} // namespace fps

std::optional<bgfx::ProgramHandle> createShaderProgram(
    const char* vert_shader_path, const char* frag_shader_path)
{
    std::string vshader;
    if (!fileops::read_file(vert_shader_path, vshader)) {
        return {};
    }

    std::string fshader;
    if (!fileops::read_file(frag_shader_path, fshader)) {
        return {};
    }

    bgfx::ShaderHandle vsh = createShader(vshader, "vshader");
    bgfx::ShaderHandle fsh = createShader(fshader, "fshader");

    return bgfx::createProgram(vsh, fsh, true);
}

template<class T>
inline void hashCombine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{

template<>
struct hash<as::vec3_t>
{
    std::size_t operator()(const as::vec3_t& vec) const
    {
        size_t seed = 0;
        hashCombine(seed, vec.x);
        hashCombine(seed, vec.y);
        hashCombine(seed, vec.z);
        return seed;
    }
};

} // namespace std

struct Vec3EqualFn
{
    bool operator()(const as::vec3_t& lhs, const as::vec3_t& rhs) const
    {
        return as::vec::equal(lhs, rhs);
    }
};

int main(int argc, char** argv)
{
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    } else {
        const int width = 800;
        const int height = 600;
        const float aspect = float(width) / float(height);
        SDL_Window* window = SDL_CreateWindow(
            argv[0], SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width,
            height, SDL_WINDOW_SHOWN);

        if (window == nullptr) {
            printf(
                "Window could not be created! SDL_Error: %s\n", SDL_GetError());
            return 1;
        }

        SDL_SysWMinfo wmi;
        SDL_VERSION(&wmi.version);
        if (!SDL_GetWindowWMInfo(window, &wmi)) {
            return 1;
        }

        bgfx::renderFrame(); // single threaded mode

        bgfx::PlatformData pd{};
#if BX_PLATFORM_WINDOWS
        pd.nwh = wmi.info.win.window;
#elif BX_PLATFORM_OSX
        pd.nwh = wmi.info.cocoa.window;
#endif // BX_PLATFORM_WINDOWS ? BX_PLATFORM_OSX

        bgfx::Init bgfx_init;
        bgfx_init.type = bgfx::RendererType::Count; // auto choose renderer
        bgfx_init.resolution.width = width;
        bgfx_init.resolution.height = height;
        bgfx_init.resolution.reset = BGFX_RESET_VSYNC;
        bgfx_init.platformData = pd;
        bgfx::init(bgfx_init);

        bgfx::setViewClear(
            0, BGFX_CLEAR_COLOR | BGFX_CLEAR_DEPTH, 0x6495EDFF, 1.0f, 0);
        bgfx::setViewRect(0, 0, 0, width, height);

        const float gizmo_offset_percent = 0.85f;
        const float gizmo_size_percent = 0.1f;
        bgfx::setViewClear(1, BGFX_CLEAR_DEPTH);
        bgfx::setViewRect(
            1, width * gizmo_offset_percent, height * gizmo_offset_percent,
            width * gizmo_size_percent, height * gizmo_size_percent);

        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();

        ImGui_Implbgfx_Init(255);
#if BX_PLATFORM_WINDOWS
        ImGui_ImplSDL2_InitForD3D(window);
#elif BX_PLATFORM_OSX
        ImGui_ImplSDL2_InitForMetal(window);
#endif // BX_PLATFORM_WINDOWS ? BX_PLATFORM_OSX

        bgfx::VertexLayout pos_col_vert_layout;
        pos_col_vert_layout.begin()
            .add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
            .add(bgfx::Attrib::Color0, 4, bgfx::AttribType::Uint8, true)
            .end();

        bgfx::VertexLayout pos_norm_vert_layout;
        pos_norm_vert_layout.begin()
            .add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
            .add(bgfx::Attrib::Normal, 3, bgfx::AttribType::Float, true)
            // .add(bgfx::Attrib::TexCoord0, 2, bgfx::AttribType::Float)
            .end();

        bgfx::VertexBufferHandle cube_col_vbh = bgfx::createVertexBuffer(
            bgfx::makeRef(cube_vertices_col, sizeof(cube_vertices_col)),
            pos_col_vert_layout);
        bgfx::IndexBufferHandle cube_col_ibh = bgfx::createIndexBuffer(
            bgfx::makeRef(cube_tri_list_col, sizeof(cube_tri_list_col)));

        bgfx::VertexBufferHandle cube_norm_vbh = bgfx::createVertexBuffer(
            bgfx::makeRef(cube_vertices_norm, sizeof(cube_vertices_norm)),
            pos_norm_vert_layout);

        const bgfx::ProgramHandle program_norm =
            createShaderProgram(
                "shader/next/v_next.bin", "shader/next/f_next.bin")
                .value_or(bgfx::ProgramHandle(BGFX_INVALID_HANDLE));

        const bgfx::ProgramHandle program_col =
            createShaderProgram(
                "shader/simple/v_simple.bin", "shader/simple/f_simple.bin")
                .value_or(bgfx::ProgramHandle(BGFX_INVALID_HANDLE));

        const bgfx::UniformHandle u_light_dir =
            bgfx::createUniform("u_lightDir", bgfx::UniformType::Vec4, 1);
        const bgfx::UniformHandle u_camera_pos =
            bgfx::createUniform("u_cameraPos", bgfx::UniformType::Vec4, 1);

        as::vec3_t light_dir{0.0f, 1.0f, -1.0f};

        asc::Camera camera{};
        // initial camera position and orientation
        auto cam_start = as::point3_t{0.0f};
        camera.look_at = cam_start;
        // camera.focal_dist = -10.0f;

        // initial mouse state
        MouseState mouse_state = mouseState();

        // camera control structure
        asc::CameraControl camera_control{};
        camera_control.pitch = camera.pitch;
        camera_control.yaw = camera.yaw;
        // camera_control.dolly = -10.0f;

        // camera properties
        asc::CameraProperties camera_props{};
        camera_props.rotate_speed = 0.005f;
        camera_props.translate_speed = 10.0f;
        camera_props.look_smoothness = 5.0f;

        auto prev = bx::getHPCounter();

        const int dimension = 30;
        auto points = mc::createPointVolume(dimension);
        auto cellValues = mc::createCellValues(dimension);
        auto cellPositions = mc::createCellPositions(dimension);

        std::vector<as::vec3_t> filtered_verts;
        std::vector<as::vec3_t> filtered_norms;
        std::unordered_map<
            as::vec3_t, as::index_t, std::hash<as::vec3_t>, Vec3EqualFn>
            unique_verts;

        Fps fps;
        for (bool quit = false; !quit;) {
            SDL_Event current_event;
            while (SDL_PollEvent(&current_event) != 0) {
                updateCameraControlKeyboardSdl(
                    current_event, camera_control, camera_props);
                ImGui_ImplSDL2_ProcessEvent(&current_event);
                if (current_event.type == SDL_QUIT) {
                    quit = true;
                    break;
                }
            }

            ImGui_Implbgfx_NewFrame();
            ImGui_ImplSDL2_NewFrame(window);

            ImGui::NewFrame();

            auto freq = double(bx::getHPFrequency());
            int64_t time_window = fps::calculateWindow(fps, bx::getHPCounter());
            double framerate = time_window > -1
                                 ? (double)(fps.MaxSamples - 1)
                                       / (double(time_window) / freq)
                                 : 0.0;

            updateCameraControlMouseSdl(
                camera_control, camera_props, mouse_state);

            // frame dt
            auto now = bx::getHPCounter();
            auto delta = now - prev;
            prev = now;

            float dt = delta / freq;

            asc::update_camera(
                camera, camera_control, camera_props, dt,
                asc::Handedness::Left);

            // marching cube scene
            {
                float view[16];
                as::mat::to_arr(as::mat4::from_affine(camera.view()), view);
                const as::mat4_t persp = as::view::perspective_d3d_lh(
                    as::deg_to_rad(35.0f), float(width) / float(height), 0.01f,
                    100.0f);

                float proj[16];
                as::mat::to_arr(persp, proj);

                bgfx::setViewTransform(0, view, proj);

                auto marching_cube_begin = bx::getHPCounter();

                static bool analytical_normals = false;
                static bool draw_normals = false;

                /* static */ as::mat3_t cam_orientation = camera.transform().rotation;

                static float camera_adjust = (float(dimension) * 0.5f) + 1.0f;

                /* static */ const as::point3_t lookat = camera.look_at;
                const as::point3_t offset =
                    lookat
                    + cam_orientation * as::vec3_t::axis_z(camera_adjust);

                static float tesselation = 1.0f;
                static float scale = 14.0f;
                static float threshold = 4.0f; // initial

                generatePointData(
                    points, dimension, scale, tesselation, offset.as_vec());
                generateCellData(cellPositions, cellValues, points, dimension);

                auto triangles =
                    mc::march(cellPositions, cellValues, dimension, threshold);

                std::vector<as::index_t> indices;
                indices.resize(triangles.size() * 3);

                as::index_t index = 0;
                as::index_t unique = 0;
                for (const auto& tri : triangles) {
                    for (int64_t i = 0; i < 3; ++i) {
                        const auto vert = tri.verts_[i];
                        const auto norm = tri.norms_[i];
                        const auto exists = unique_verts.find(vert);
                        if (exists == std::end(unique_verts)) {
                            filtered_verts.push_back(vert);
                            filtered_norms.push_back(norm);
                            unique_verts.insert({vert, unique});
                            indices[index] = unique;
                            unique++;
                        } else {
                            indices[index] = exists->second;
                        }
                        index++;
                    }
                }

                uint32_t max_vertices = 32 << 10;
                bgfx::TransientVertexBuffer mc_triangle_tvb;
                bgfx::allocTransientVertexBuffer(
                    &mc_triangle_tvb, max_vertices, pos_norm_vert_layout);

                bgfx::TransientIndexBuffer tib;
                bgfx::allocTransientIndexBuffer(&tib, max_vertices);

                PosNormalUv* vertex = (PosNormalUv*)mc_triangle_tvb.data;
                int16_t* index_data = (int16_t*)tib.data;

                for (as::index_t i = 0; i < filtered_verts.size(); i++) {
                    vertex[i].normal = analytical_normals
                                         ? as::vec::normalize(filtered_norms[i])
                                         : as::vec3_t::zero();
                    vertex[i].position = filtered_verts[i];
                }

                for (as::index_t indice = 0; indice < indices.size();
                     indice++) {
                    index_data[indice] = indices[indice];
                }

                if (!analytical_normals) {
                    for (as::index_t indice = 0; indice < indices.size();
                         indice += 3) {
                        const as::vec3_t e1 =
                            filtered_verts[indices[indice]]
                            - filtered_verts[indices[indice + 1]];
                        const as::vec3_t e2 =
                            filtered_verts[indices[indice + 2]]
                            - filtered_verts[indices[indice + 1]];
                        const as::vec3_t normal = as::vec3::cross(e1, e2);

                        vertex[indices[indice]].normal += normal;
                        vertex[indices[indice + 1]].normal += normal;
                        vertex[indices[indice + 2]].normal += normal;
                    }

                    for (as::index_t i = 0; i < filtered_verts.size(); i++) {
                        vertex[i].normal = as::vec::normalize(vertex[i].normal);
                    }
                }

                float model[16];
                as::mat::to_arr(as::mat4_t::identity(), model);
                bgfx::setTransform(model);

                bgfx::setIndexBuffer(&tib, 0, indices.size());
                bgfx::setVertexBuffer(
                    0, &mc_triangle_tvb, 0, filtered_verts.size());
                bgfx::setState(BGFX_STATE_DEFAULT);
                bgfx::submit(0, program_norm);

                bgfx::TransientVertexBuffer mc_line_tvb;
                bgfx::allocTransientVertexBuffer(
                    &mc_line_tvb, max_vertices, pos_col_vert_layout);

                PosColorVertex* vertex_normals =
                    (PosColorVertex*)mc_line_tvb.data;

                for (as::index_t i = 0; i < filtered_verts.size(); i++) {
                    vertex_normals->position = vertex[i].position;
                    vertex_normals->abgr = 0xff000000;
                    vertex_normals++;
                    vertex_normals->position =
                        vertex[i].position + vertex[i].normal;
                    vertex_normals->abgr = 0xff000000;
                    vertex_normals++;
                }

                if (draw_normals) {
                    float identity[16];
                    as::mat::to_arr(as::mat4_t::identity(), identity);
                    bgfx::setTransform(identity);

                    bgfx::setState(BGFX_STATE_DEFAULT | BGFX_STATE_PT_LINES);

                    bgfx::setVertexBuffer(
                        0, &mc_line_tvb, 0, filtered_verts.size() * 2);
                    bgfx::submit(0, program_col);
                }

                static float rot = 0.0f;
                static float shear1 = 0.0f;
                static float shear2 = 0.0f;
                static float spin_speed = 0.0f;

                const auto a = as::mat4::shear_y(shear1, shear2);
                const auto m = as::mat4::from_mat3(as::mat3::rotation_y(rot));

                float cmodel[16];
                as::mat::to_arr(a, cmodel);
                rot += dt * spin_speed;

                bgfx::setState(BGFX_STATE_DEFAULT);

                bgfx::setUniform(u_light_dir, (void*)&light_dir, 1);
                bgfx::setUniform(u_camera_pos, (void*)&camera.look_at, 1);

                bgfx::setTransform(cmodel);

                bgfx::setVertexBuffer(0, cube_norm_vbh);
                bgfx::submit(0, program_norm);

                bgfx::TransientVertexBuffer cube_line_tvb;
                bgfx::allocTransientVertexBuffer(
                    &cube_line_tvb, max_vertices, pos_col_vert_layout);

                PosColorVertex* cube_normal =
                    (PosColorVertex*)cube_line_tvb.data;

                for (as::index_t i = 0; i < 36; i++) {
                    cube_normal->position = cube_vertices_norm[i].position;
                    cube_normal->abgr = 0xff000000;
                    cube_normal++;
                    cube_normal->position = cube_vertices_norm[i].position
                                          + cube_vertices_norm[i].normal;
                    cube_normal->abgr = 0xff000000;
                    cube_normal++;
                }

                bgfx::setTransform(cmodel);

                bgfx::setState(BGFX_STATE_DEFAULT | BGFX_STATE_PT_LINES);

                bgfx::setVertexBuffer(0, &cube_line_tvb, 0, 72);
                bgfx::submit(0, program_col);

                const double to_ms = 1000.0 / freq;
                auto marching_cube_time =
                    double(bx::getHPCounter() - marching_cube_begin);

                // project from world to screen
                // const as::vec4_t clip =
                //     persp * as::mat4::from_affine(camera.view()) * as::mat4::translation(a);
                // const as::vec3_t ndc = as::vec3::from_vec4(clip / clip.w);
                // const as::vec2_t window =
                //     (as::vec2::from_vec3(ndc) + as::vec2_t{1.0f}) / 2.0f;

                // ImGui::SetWindowPos(ImVec2(
                //     window.x * float(width),
                //     height - (window.y * float(height))));

                ImGui::Text("Framerate: ");
                ImGui::SameLine(100);
                ImGui::Text("%f", framerate);

                ImGui::Text("Marching Cube update: ");
                ImGui::SameLine(160);
                ImGui::Text("%f", marching_cube_time * to_ms);

                float light_dir_arr[3];
                as::vec::to_arr(light_dir, light_dir_arr);
                ImGui::InputFloat3("Light Dir", light_dir_arr, 3);
                light_dir = as::vec::from_arr(light_dir_arr);

                ImGui::SliderFloat("Spin Speed", &spin_speed, 0.0f, 10.0f);
                ImGui::SliderFloat("sheer1", &shear1, -10.0f, 10.0f);
                ImGui::SliderFloat("sheer2", &shear2, -10.0f, 10.0f);
                ImGui::SliderFloat("Threshold", &threshold, 0.0f, 10.0f);
                ImGui::SliderFloat("Back", &camera_adjust, 0.0f, 100.0f);
                ImGui::SliderFloat("Scale", &scale, 0.0f, 100.0f);
                ImGui::SliderFloat("Tesselation", &tesselation, 0.001f, 10.0f);
                ImGui::Checkbox("Draw Normals", &draw_normals);
                ImGui::Checkbox("Analytical Normals", &analytical_normals);

                ImGui::InputFloat("Camera Smoothness", &camera_props.look_smoothness, 0.0f, 100.0f);
            }

            // gizmo cube
            {
                float view[16];
                as::mat::to_arr(
                    as::mat4::from_mat3_vec3(
                        camera.view().rotation,
                        as::vec3_t::axis_z(10.0f)),
                    view);

                const float extent = 10.0f * aspect;

                float proj[16];
                as::mat::to_arr(
                    as::view::ortho_d3d_lh(
                        -extent, extent, -10.0f, 10.0f, 0.001f, 100.0f),
                    proj);

                bgfx::setViewTransform(1, view, proj);

                as::mat4_t rot = as::mat4::from_mat3(as::mat3::scale(4.0f));

                float model[16];
                as::mat::to_arr(rot, model);

                bgfx::setTransform(model);

                bgfx::setVertexBuffer(0, cube_col_vbh);
                bgfx::setIndexBuffer(cube_col_ibh);

                bgfx::submit(1, program_col);
            }

            filtered_verts.clear();
            filtered_norms.clear();
            unique_verts.clear();

            ImGui::Render();

            bgfx::frame();
        }

        mc::destroyCellValues(cellValues, dimension);
        mc::destroyCellPositions(cellPositions, dimension);
        mc::destroyPointVolume(points, dimension);

        bgfx::destroy(u_camera_pos);
        bgfx::destroy(u_light_dir);
        bgfx::destroy(cube_col_vbh);
        bgfx::destroy(cube_col_ibh);
        bgfx::destroy(cube_norm_vbh);
        bgfx::destroy(program_norm);
        bgfx::destroy(program_col);

        ImGui_ImplSDL2_Shutdown();
        ImGui_Implbgfx_Shutdown();

        ImGui::DestroyContext();
        bgfx::shutdown();

        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    return 0;
}
