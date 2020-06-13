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

#include <optional>

struct PosColorVertex
{
    float x;
    float y;
    float z;
    uint32_t abgr;
};

struct PosNormalUv
{
    as::vec3_t position;
    as::vec3_t normal;
    // as::vec2_t uv;
};

static PosColorVertex cube_vertices[] = {
    {-1.0f, 1.0f, 1.0f, 0xff000000},   {1.0f, 1.0f, 1.0f, 0xff0000ff},
    {-1.0f, -1.0f, 1.0f, 0xff00ff00},  {1.0f, -1.0f, 1.0f, 0xff00ffff},
    {-1.0f, 1.0f, -1.0f, 0xffff0000},  {1.0f, 1.0f, -1.0f, 0xffff00ff},
    {-1.0f, -1.0f, -1.0f, 0xffffff00}, {1.0f, -1.0f, -1.0f, 0xffffffff},
};

static const uint16_t cube_tri_list[] = {
    0, 1, 2, 1, 3, 2, 4, 6, 5, 5, 6, 7, 0, 2, 4, 4, 2, 6,
    1, 5, 3, 5, 7, 3, 0, 4, 1, 4, 5, 1, 2, 3, 6, 6, 3, 7,
};

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

        bgfx::VertexLayout pos_norm_uv_layout;
        pos_norm_uv_layout.begin()
            .add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
            .add(bgfx::Attrib::Normal, 3, bgfx::AttribType::Float, true)
            // .add(bgfx::Attrib::TexCoord0, 2, bgfx::AttribType::Float)
            .end();

        bgfx::VertexBufferHandle vbh = bgfx::createVertexBuffer(
            bgfx::makeRef(cube_vertices, sizeof(cube_vertices)),
            pos_col_vert_layout);
        bgfx::IndexBufferHandle ibh = bgfx::createIndexBuffer(
            bgfx::makeRef(cube_tri_list, sizeof(cube_tri_list)));

        const bgfx::ProgramHandle program_norm =
            createShaderProgram(
                "shader/next/v_next.bin", "shader/next/f_next.bin")
                .value_or(bgfx::ProgramHandle(BGFX_INVALID_HANDLE));

        const bgfx::ProgramHandle program_col =
            createShaderProgram(
                "shader/simple/v_simple.bin", "shader/simple/f_simple.bin")
                .value_or(bgfx::ProgramHandle(BGFX_INVALID_HANDLE));

        asc::Camera camera{};
        // initial camera position and orientation
        auto cam_start = as::vec3_t{0.0f};
        camera.look_at = cam_start;

        // initial mouse state
        MouseState mouse_state = mouseState();

        // camera control structure
        asc::CameraControl camera_control{};
        camera_control.pitch = camera.pitch;
        camera_control.yaw = camera.yaw;

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

            // main cube
            {
                float view[16];
                as::mat::to_arr(camera.view(), view);

                float proj[16];
                as::mat::to_arr(
                    as::view::perspective_d3d_lh(
                        as::deg_to_rad(35.0f), float(width) / float(height),
                        0.01f, 100.0f),
                    proj);

                bgfx::setViewTransform(0, view, proj);

                auto marching_cube_begin = bx::getHPCounter();

                /* static */ as::mat3_t cam_orientation =
                    as::mat3::from_mat4(camera.transform());

                static float camera_adjust = (float(dimension) * 0.5f) + 1.0f;

                /* static */ const as::vec3_t lookat = camera.look_at;
                const as::vec3_t offset =
                    lookat
                    + cam_orientation * as::vec3_t::axis_z(camera_adjust);

                static float scale = 14.0f;
                generatePointData(points, dimension, scale, offset);
                generateCellData(cellPositions, cellValues, points, dimension);

                static float threshold = 4.0f;
                auto triangles =
                    mc::march(cellPositions, cellValues, dimension, threshold);

                uint32_t max_vertices = 32 << 10;
                bgfx::TransientVertexBuffer tvb;
                bgfx::allocTransientVertexBuffer(
                    &tvb, max_vertices, pos_norm_uv_layout);

                PosNormalUv* vertex = (PosNormalUv*)tvb.data;

                int vertCount = 0;
                for (const auto& tri : triangles) {
                    const as::vec3_t e1 = tri.verts_[1] - tri.verts_[0];
                    const as::vec3_t e2 = tri.verts_[2] - tri.verts_[0];
                    const as::vec3_t normal = as::vec::abs(
                        as::vec::normalize(as::vec3::cross(e1, e2)));
                    for (const auto& vert : tri.verts_) {
                        vertex->position = vert;
                        vertex->normal = normal;
                        vertex++;
                        vertCount++;
                    }
                }

                float model[16];
                as::mat::to_arr(as::mat4_t::identity(), model);
                bgfx::setTransform(model);

                bgfx::setVertexBuffer(0, &tvb, 0, vertCount);
                bgfx::setState(BGFX_STATE_DEFAULT);
                bgfx::submit(0, program_norm);

                float cmodel[16];
                as::mat::to_arr(
                    as::mat4::from_mat3_vec3(cam_orientation, lookat), cmodel);
                bgfx::setTransform(cmodel);

                bgfx::setVertexBuffer(0, vbh);
                bgfx::setIndexBuffer(ibh);
                bgfx::submit(0, program_col);

                const double toMs = 1000.0 / freq;
                auto marching_cube_time =
                    double(bx::getHPCounter() - marching_cube_begin);

                ImGui::Text("Framerate: ");
                ImGui::SameLine(100);
                ImGui::Text("%f", framerate);

                ImGui::Text("Marching Cube update: ");
                ImGui::SameLine(160);
                ImGui::Text("%f", marching_cube_time * toMs);

                ImGui::SliderFloat("Threshold", &threshold, 0.0f, 10.0f);
                ImGui::SliderFloat("Back", &camera_adjust, 0.0f, 100.0f);
                ImGui::SliderFloat("Scale", &scale, 0.0f, 100.0f);
            }

            // gizmo cube
            {
                float view[16];
                as::mat::to_arr(
                    as::mat4::from_mat3_vec3(
                        as::mat3::from_mat4(camera.view()),
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

                bgfx::setVertexBuffer(0, vbh);
                bgfx::setIndexBuffer(ibh);

                bgfx::submit(1, program_col);
            }

            ImGui::Render();

            bgfx::frame();
        }

        mc::destroyCellValues(cellValues, dimension);
        mc::destroyCellPositions(cellPositions, dimension);
        mc::destroyPointVolume(points, dimension);

        bgfx::destroy(vbh);
        bgfx::destroy(ibh);
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
