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
#define PASSINDEX 0
#define b0 backbuffer
#endif

const float PI = 3.14159;
#define beat (time*140./60.)

#define repeat(p,r) (mod(p,r)-r/2.)
mat2 rot (float a) { float c=cos(a), s=sin(a); return mat2(c,s,-s,c); }
float smoothmin (float a, float b, float r) { float h = clamp(.5+.5*(b-a)/r, 0., 1.); return mix(b, a, h)-r*h*(1.-h); }
float sdCylinderBox (vec2 p, vec2 r) { vec2 b = abs(p)-r; return min(0.0, max(b.x, b.y)) + length(max(b,0.0)); }
float sdTorus( vec3 p, vec2 t ) {  vec2 q = vec2(length(p.xz)-t.x,p.y); return length(q)-t.y; }
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
  // float chilly = noise(pos * 2.);
  // float salty = fbm(pos*20.);
  // float spicy = chilly*.1 + salty*.01;
  float scene = 1.0;
  float s0 = floor(beat/2.);
  float b0 = mod(beat/2., 3.0)/3.;
  int stage = 4;
  int fly = 0;

  s0 = floor(time/10.);
  b0 = 0.;

  // stage = int(step(34., time)) + int(step(48., time)) + int(step(75., time)) + int(step(102., time));
  // fly = int(step(102., time));

  if (fly == 1) {
      pos.xz *= rot(b0+s0);
      pos.yz *= rot(b0+s0);
      pos.yx *= rot(b0+s0);
  }

  if (stage == 0) {

    scene = length(pos)-1.0;

  } else if (stage == 1) {

    vec3 cell = vec3(.4,.5,.3);
    // pos.z += time * .5;
    vec3 id = floor(pos/cell);
    pos.xy *= rot(id.z);
    pos.x = repeat(pos.x, cell.x);
    pos.y = repeat(pos.y, cell.y);
    pos.z = repeat(pos.z, cell.z);
    scene = min(
      sdCylinderBox(pos.xz, vec2(.05,.005)),
      sdCylinderBox(pos.yz, vec2(.01,.005)));

  } else if (stage == 2) {

    float cell = 4.;
    pos.z += time;
    float id = floor(pos.z/cell);
    float tunnel = length(pos.xy)-1.+.2*sin(id);
    pos.z = repeat(pos.z, cell);
    float amplitude = 1.0;
    for (int index = 0; index < 5; ++index) {
      pos = abs(pos)-.8*amplitude;
      pos.zx *= rot(-.5*amplitude+id*3.);
      scene = min(scene, abs(max(pos.x, max(pos.y, pos.z)))-.15*amplitude);
      amplitude /= 2.0;
    }
    scene = max(0.0, -scene);
    scene = max(scene, -tunnel);

  } else if (stage == 3) {

      float amplitude = 1.0;
      float range = .7;
      float ay = .4;
      float ax = -.4;
      float az = -.8;
      float blend = .1;
      float radius = .05;
      float thin = .002;
      float wave = 0.75+0.25*sin(length(pos)*2.-beat*PI);
      for (int index = 0; index < 4; ++index) {
        float w = 4.*time*smoothstep(.2, 1., float(index)/6.);
        pos = abs(pos)-(range*wave)*amplitude;
        pos.xz *= rot(ay/amplitude+w);
        pos.yz *= rot(ax/amplitude+w);
        pos.yx *= rot(az/amplitude+w);
        scene = smoothmin(scene, sdTorus(pos, vec2(radius, thin)), blend);
        amplitude /= 2.;
      }

  } else if (stage == 4) {

    float amplitude = 1.0;
    float range = .1+.4*b0;
    float ay = .4+.1*b0+s0*.2;
    float ax = -.2-.4*b0+s0*4.;
    float az = -.5-.2*b0+s0;
    for (int index = 0; index < 6; ++index) {
      pos = abs(pos)-range*amplitude;
      pos.xz *= rot(ay*amplitude);
      pos.yz *= rot(ax*amplitude);
      pos.yx *= rot(az*amplitude);
      pos = abs(pos)-range*amplitude*.1;
      scene = min(scene, sdCylinderBox(pos.xz, vec2(.001*amplitude)));
    }
  }

  return scene;
}

vec3 getNormal (vec3 pos) {
  vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
  return normalize( e.xyy*map( pos + e.xyy ) + e.yyx*map( pos + e.yyx ) + e.yxy*map( pos + e.yxy ) + e.xxx*map( pos + e.xxx ) );
}

vec4 raymarch (vec3 eye, vec3 ray) {
  vec4 result = vec4(eye, 0);
  float dither = random(ray.xy);
  float total = dither * .2;
  for (float index = 50.0; index > 0.0; --index) {
    result.xyz = eye + ray * total;
    float dist = map(result.xyz);
    if (dist < 0.001 + total * 1.0/synth_Resolution.y) {
      result.a = index/50.;
      break;
    }
    dist *= 0.9 + 0.1 * dither;
    total += dist;
  }
  return result;
}

vec3 anaglyph (vec3 eye, vec3 at, vec2 uv) {
  vec3 offset = vec3(.02,0,0);
  vec3 eyeLeft = eye-offset;
  vec3 eyeRight = eye+offset;
  vec4 resultLeft = raymarch(eyeLeft, look(eyeLeft, at+offset, uv));
  vec4 resultRight = raymarch(eyeRight, look(eyeRight, at-offset, uv));
  return vec3(resultLeft.a, vec2(resultRight.a));
}

vec3 dots (vec3 eye, vec3 at, vec2 uv) {
  float lod = synth_Resolution.y/8.;
  vec2 uvpixel = floor(uv * lod + .5) / lod;
  vec4 result = raymarch(eye, look(eye, at, uvpixel));
  vec3 color = vec3(result.w);
  float depth = smoothstep(10.0, 0., length(eye-result.xyz));
  float cell = 1./lod;
  uv = repeat(uv+cell/2., cell);
  float shape = smoothstep(cell/3.+cell/10., cell/3.,length(uv));
  return color * shape;
}

void main() {

  if (PASSINDEX == 0) {

    vec2 uv = (gl_FragCoord.xy-0.5*synth_Resolution)/synth_Resolution.y;
    vec3 eye = vec3(.01,.01,-4.);
    vec3 at = vec3(0);
    vec3 color = dots(eye, at, uv);
    vec4 frame = texture2D(b0, gl_FragCoord.xy / synth_Resolution);
    gl_FragColor = frame*.9 + .1*vec4(color, 1);

  } else if (PASSINDEX == 1) {

    vec2 uv = gl_FragCoord.xy / synth_Resolution;
    // vec4 frame = texture2D(b1, uv);
    gl_FragColor = texture2D(b0, uv);

  }

}


//! <preset file="shader.preset" />
