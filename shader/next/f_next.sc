$input v_normal

#include <../bgfx_shader.sh>

void main() {
    gl_FragColor = vec4(normalize(v_normal), 1.0);
}
