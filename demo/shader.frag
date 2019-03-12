precision highp float;

uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;
#define synth_Resolution resolution
#define veda

// begin

uniform sampler2D b0;
uniform sampler2D b1;
uniform int	PASSINDEX;

#ifdef veda
#define PASSINDEX 1
#define b0 backbuffer
#endif


#define repeat(p,r) (mod(p,r)-r/2.)
const float PI = 3.14159;
mat2 rot (float a) { float c=cos(a), s=sin(a); return mat2(c,s,-s,c); }
float smoothmin (float a, float b, float r) { float h = clamp(.5+.5*(b-a)/r, 0., 1.); return mix(b, a, h)-r*h*(1.-h); }
float sdCylinderBox (vec2 p, vec2 r) { vec2 b = abs(p)-r; return min(0.0, max(b.x, b.y)) + length(max(b,0.0)); }
float random (in vec2 st) { return fract(sin(dot(st.xy,vec2(12.9898,78.233)))*43758.5453123); }
float hash(float n) { return fract(sin(n) * 1e4); }
float noise(vec3 x) {
    const vec3 step = vec3(110, 241, 171);
    vec3 i = floor(x);
    vec3 f = fract(x);
    float n = dot(i, step);
    vec3 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix( hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y),
               mix(mix( hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}
float fbm (vec3 p) {
  float amplitude = 0.5;
  float result = 0.0;
  for (float index = 0.0; index <= 3.0; ++index) {
    result += noise(p / amplitude) * amplitude;
    amplitude /= 2.;
  }
  return result;
}
vec3 look (vec3 eye, vec3 target, vec2 anchor) {
    vec3 forward = normalize(target-eye);
    vec3 right = normalize(cross(forward, vec3(0,1,0)));
    vec3 up = normalize(cross(right, forward));
    return normalize(forward + right * anchor.x + up * anchor.y);
}
void moda(inout vec2 p, float repetitions) {
	float angle = 2.*PI/repetitions;
	float a = atan(p.y, p.x) + angle/2.;
	a = mod(a,angle) - angle/2.;
	p = vec2(cos(a), sin(a))*length(p);
}

float map (vec3 pos) {
  float chilly = noise(pos * 2.);
  float salty = fbm(pos*20.);
  float spicy = chilly*.1 + salty*.01;

  float cell = .5;
  pos.z += time * .1;
  pos = repeat(pos, cell);
  float scene = min(
    sdCylinderBox(pos.xz, vec2(.05,.005)),
    sdCylinderBox(pos.yz, vec2(.02,.005)));
  cell = .01;
  scene = min(
    max(scene, (length(repeat(pos+cell/2., cell/2.))-cell/8.)),
    max(scene, (length(repeat(pos+cell/2., cell))-cell/4.))
  );

  return scene;
}

vec3 getNormal (vec3 pos) {
  vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
  return normalize( e.xyy*map( pos + e.xyy ) + e.yyx*map( pos + e.yyx ) + e.yxy*map( pos + e.yxy ) + e.xxx*map( pos + e.xxx ) );
}

vec3 background (vec2 uv) {
  vec3 c = mix(vec3(1), vec3(1, 0.439, 0.545), uv.y+.85);
  // vec3 c = mix(vec3(1), vec3(0), uv.y+.85);
  float cell = .05;
  vec2 id = floor(uv/cell);
  float t = length(id)*.1;
  float r = .02-.0025*abs(sin(t*2.));
  uv = repeat(uv,cell);
  uv *= rot(t);
  c = mix(c, vec3(1, 0.949, 0.439), clamp(
  // c = mix(c, vec3(1), clamp(
    .0005/
    abs(
      max(abs(uv.x),abs(uv.y))
    -r),
  0., 1.));
  return c;
}

vec3 sketch (vec2 uv) {
  vec3 b = background(uv);
  float t0 = smoothstep(.0, 1.,time);
  float t1 = time*PI*8.;
  float t2 = time*PI*2.;
  float b1 = (sin(time*PI*4.));
  // uv.y *= 1.+t0*.02*cos(t1);
  // uv.x *= 1.+t0*.1*sin(t1)*sin(t2);
  float r = .2;
  vec2 offset = vec2(.01,0);
  vec3 c = vec3(0);
  c.r += smoothstep(r+.001,r,length(uv+offset));
  c.gb += smoothstep(r+.001,r,length(uv-offset));
  return c;
}

vec4 raymarch (vec3 eye, vec3 ray, float dither) {
  vec4 result = vec4(eye, 0);
  float total = dither * .2;
  for (float index = 30.0; index > 0.0; --index) {
    result.xyz = eye + ray * total;
    float dist = map(result.xyz);
    if (dist < 0.001 + total * 1.0/synth_Resolution.y) {
      result.a = pow(index / 30.0, 1./2.);
      break;
    }
    dist *= 0.9 + 0.1 * dither;
    total += dist;
  }
  return result;
}

void main() {
  if (PASSINDEX == 0) {
    gl_FragColor = texture2D(b1, gl_FragCoord.xy / synth_Resolution);
  } else if (PASSINDEX == 1) {
    vec2 uv = (gl_FragCoord.xy-0.5*synth_Resolution)/synth_Resolution.y;
    float dither = random(uv);
    vec3 offset = vec3(.02,0,0);
    vec3 eyeLeft = vec3(.01,.01,-4.)-offset;
    vec3 eyeRight = vec3(.01,.01,-4.)+offset;
    vec3 at = vec3(0);
    vec4 resultLeft = raymarch(eyeLeft, look(eyeLeft, at+offset, uv), dither);
    vec4 resultRight = raymarch(eyeRight, look(eyeRight, at-offset, uv), dither);
    // vec3 normal = getNormal(result.xyz);
    vec3 color = vec3(resultLeft.a, vec2(resultRight.a));
    // gl_FragColor = texture2D(b0, gl_FragCoord.xy / synth_Resolution) * .9 + .1 * vec4(color, 1);
    gl_FragColor = vec4(color, 1);
  }
}


//! <preset file="shader.preset" />
