$input v_normal

uniform vec3 u_lightDir;

#include <../bgfx_shader.sh>

void main() {
    vec3 normal = normalize(v_normal);
    float light_amount = max(0.0, dot(normal, normalize(u_lightDir.xyz)));
    vec3 color = abs(normal);
    gl_FragColor = vec4(color * light_amount, 1.0);
}
