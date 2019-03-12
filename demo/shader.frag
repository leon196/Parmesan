uniform float time;
uniform vec2 synth_Resolution;

// begin

uniform sampler2D b0;
uniform sampler2D b1;
uniform int	PASSINDEX;

void main() {
  if (PASSINDEX == 0) {
    gl_FragColor = texture2D(b1, gl_FragCoord.xy / synth_Resolution);
  } else if (PASSINDEX == 1) {
    vec2 uv = (gl_FragCoord.xy-0.5*synth_Resolution)/synth_Resolution.y;
    vec3 color = vec3(1.0);
    uv.x += sin(time);
    color *= step(length(uv), 0.5);
    float blend = 0.9;
    gl_FragColor = texture2D(b0, gl_FragCoord.xy / synth_Resolution) * blend + (1.0-blend) * vec4(color, 1);
  }
}


//! <preset file="shader.preset" />
