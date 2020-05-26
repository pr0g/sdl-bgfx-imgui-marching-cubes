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

struct PosColorVertex
{
    float x;
    float y;
    float z;
    uint32_t abgr;
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

        bgfx::VertexBufferHandle vbh = bgfx::createVertexBuffer(
            bgfx::makeRef(cube_vertices, sizeof(cube_vertices)),
            pos_col_vert_layout);
        bgfx::IndexBufferHandle ibh = bgfx::createIndexBuffer(
            bgfx::makeRef(cube_tri_list, sizeof(cube_tri_list)));

        std::string vshader;
        if (!fileops::read_file("shader/v_simple.bin", vshader)) {
            return 1;
        }

        std::string fshader;
        if (!fileops::read_file("shader/f_simple.bin", fshader)) {
            return 1;
        }

        bgfx::ShaderHandle vsh = createShader(vshader, "vshader");
        bgfx::ShaderHandle fsh = createShader(fshader, "fshader");

        bgfx::ProgramHandle program = bgfx::createProgram(vsh, fsh, true);

        asc::Camera camera{};
        // initial camera position and orientation
        camera.look_at = as::vec3_t{22.84f, 16.43f, 37.43f};
        camera.pitch = 0.3f;
        camera.yaw = 3.6f;

        // initial mouse state
        MouseState mouse_state = mouseState();

        // camera control structure
        asc::CameraControl camera_control{};
        camera_control.pitch = camera.pitch;
        camera_control.yaw = camera.yaw;

        // camera properties
        asc::CameraProperties camera_props{};
        camera_props.rotate_speed = 0.01f;
        camera_props.translate_speed = 1.0f;
        camera_props.look_smoothness = 5.0f;

        auto prev = bx::getHPCounter();

        for (bool quit = false; !quit;) {
            SDL_Event currentEvent;
            while (SDL_PollEvent(&currentEvent) != 0) {
                updateCameraControlKeyboardSdl(
                    currentEvent, camera_control, camera_props);
                ImGui_ImplSDL2_ProcessEvent(&currentEvent);
                if (currentEvent.type == SDL_QUIT) {
                    quit = true;
                    break;
                }
            }

            ImGui_Implbgfx_NewFrame();
            ImGui_ImplSDL2_NewFrame(window);

            ImGui::NewFrame();
            ImGui::ShowDemoWindow(); // your drawing here
            ImGui::Render();

            updateCameraControlMouseSdl(
                camera_control, camera_props, mouse_state);

            auto now = bx::getHPCounter();
            auto delta = now - prev;
            prev = now;

            float dt = delta / (float)bx::getHPFrequency();

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
                        0.1f, 100.0f),
                    proj);

                bgfx::setViewTransform(0, view, proj);

                //

                const int dimension = 15;
                auto points = create_point_volume(dimension);
                auto cells = create_cell_volume(dimension);

                generate_point_data(
                    points, dimension, /*offset*/ as::vec3_t::zero(),
                    /*scale*/ 14.0f);
                generate_cell_data(cells, points, dimension);

                auto triangles = march(cells, dimension, 4.0f /*threshold*/);

                uint32_t maxVertices = 32 << 10;
                bgfx::TransientVertexBuffer tvb;
                bgfx::allocTransientVertexBuffer(
                    &tvb, maxVertices, pos_col_vert_layout);

                PosColorVertex* vertex = (PosColorVertex*)tvb.data;

                const float increment = 360.0f / (triangles.size() * 3);
                float h = 0.0f;

                int vertCount = 0;
                for (const auto& tri : triangles) {
                    for (const auto& vert : tri.verts_) {
                        vertex->x = vert.x;
                        vertex->y = vert.y;
                        vertex->z = vert.z;

                        float r, g, b;
                        HSVtoRGB(r, g, b, h, 1.0f, 1.0f);
                        uint8_t ri = r * 255;
                        uint8_t gi = g * 255;
                        uint8_t bi = b * 255;

                        vertex->abgr = ri | gi << 8 | bi << 16 | 0xff << 24;
                        vertex++;
                        vertCount++;
                        h += increment;
                    }
                }

                as::mat4_t rot = as::mat4_t::identity();

                float model[16];
                as::mat::to_arr(rot, model);

                bgfx::setTransform(model);

                bgfx::setVertexBuffer(0, &tvb, 0, vertCount);

                bgfx::submit(0, program);

                destroy_cell_volume(cells, dimension);
                destroy_point_volume(points, dimension);
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

                bgfx::submit(1, program);
            }

            bgfx::frame();
        }

        bgfx::destroy(vbh);
        bgfx::destroy(ibh);
        bgfx::destroy(program);

        ImGui_ImplSDL2_Shutdown();
        ImGui_Implbgfx_Shutdown();

        ImGui::DestroyContext();
        bgfx::shutdown();

        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    return 0;
}
