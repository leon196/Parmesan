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
const int STAGE_INTRO = 0;
const int STAGE_TUNNEL = 1;
const int STAGE_CITY = 2;
const int STAGE_RING = 3;
const int STAGE_KIF = 4;
#define beat (time*140./60./2.)
#define repeat(p,r) (mod(p,r)-r/2.)
int getStage () {
  return STAGE_CITY;
  #ifdef veda
  float t = mod(time, 100.);
  #else
  float t = time;
  #endif
  return int(step(34., t)) + int(step(48., t)) + int(step(75., t)) + int(step(102., t));
}

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
float sinc( float x, float k ) {
    float a = PI * (k*x-1.0);
    return sin(a)/a;
}
float sequence (float a, float b) {
  // #ifdef veda
  // float t = mod(time, 20.);
  // #else
  float t = time;
  // #endif
  return smoothstep(b+.1,b,t) * smoothstep(a,a+.1,t);
}

float map (vec3 pos) {
  // float chilly = noise(pos * 2.);
  // float salty = fbm(pos*20.);
  // float spicy = chilly*.1 + salty*.01;
  float scene = 1.0;
  float s0 = floor(beat);
  float b0 = mod(beat, 3.0)/3.;

  // s0 = floor(time/10.);
  // b0 = 0.;

  // stage = int(step(34., time)) + int(step(48., time)) + int(step(75., time)) + int(step(102., time));
  // fly = int(step(102., time));

  int stage = getStage();

  if (stage == STAGE_INTRO) {

    float bounce = sinc(mod(beat * 2., 1.), 5.);
    pos.y += bounce * .25 * sequence(0., 7.0);
    pos.z += sinc(mod(beat, 1.), 10.) * sequence(7., 35.);
    pos.xy += vec2(
      clamp(mod(mix(time, -time, step(1., mod(time, 2.))), 1.), 0., 1.) * 2. - 1.,
      -abs(sin(time*5.0))+.5) * sequence(17., 35.);
    // pos.y += mod(mix(time, -time, step(1., mod(time, 2.))), 1.);
    scene = length(pos)-1.0;

  } else if (stage == STAGE_TUNNEL) {

    vec3 cell = vec3(.4,.5,.9);
    // pos.z += time * .5;
    vec3 id = floor(pos/cell);
    pos.xy *= rot(id.z);
    pos.x = repeat(pos.x, cell.x);
    pos.y = repeat(pos.y, cell.y);
    pos.z = repeat(pos.z, cell.z);
    scene = min(
      sdCylinderBox(pos.xz, vec2(.05,.005)),
      sdCylinderBox(pos.yz, vec2(.01,.005)));

  } else if (stage == STAGE_CITY) {

    float cell = 2.;
    vec3 p = pos;
    pos.z += time;
    float id = floor(pos.z/cell);
    float sid = sin(id);
    float tunnel = length(pos.xy)-.5;//+.2*sin(id);
    // tunnel = min(tunnel, (length(pos.xy)-.4+.1*sin(id)));
    pos.z = repeat(pos.z, cell);
    float amplitude = 1.0;
    for (int index = 0; index < 7; ++index) {
      pos = abs(pos)-(.7+.2*sid)*amplitude;
      pos.zx *= rot(-.5*amplitude+id);
      scene = min(scene, abs(abs(max(pos.x, max(pos.y, pos.z)))-(.4+.2*sid)*amplitude)-.1*amplitude);
      amplitude /= 2.0;
    }
    scene = max(0.0, -scene);
    scene = max(scene, -tunnel);

  } else if (stage == STAGE_RING) {

    float amplitude = 1.0;
    float range = .7;
    float ay = .4;
    float ax = -.4;
    float az = -.8;
    float wave = 0.75+0.25*sin(length(pos)*2.-2.7-beat*PI*2.);
    float blend = .1;
    float radius = .2;
    float thin = .05;
    for (int index = 0; index < 3; ++index) {
      float w = 4.*time*smoothstep(.2, 1., float(index)/6.);
      pos = abs(pos)-(range*wave)*amplitude;
      pos.xz *= rot(ay/amplitude+w);
      pos.yz *= rot(ax/amplitude+w);
      pos.yx *= rot(az/amplitude+w);
      scene = smoothmin(scene, sdTorus(pos, vec2(radius, thin)*amplitude), blend);
      amplitude /= 2.;
    }

  } else if (stage == STAGE_KIF) {

    pos.xz *= rot(b0+s0);
    pos.yz *= rot(b0+s0);
    pos.yx *= rot(b0+s0);
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
  vec2 e = vec2(1.0,-1.0)*0.5773*0.00005;
  return normalize( e.xyy*map( pos + e.xyy ) + e.yyx*map( pos + e.yyx ) + e.yxy*map( pos + e.yxy ) + e.xxx*map( pos + e.xxx ) );
}

vec4 raymarch (vec3 eye, vec3 ray) {
  vec4 result = vec4(eye, 0);
  float dither = random(ray.xy+fract(time));
  float total = dither * .2;
  for (float index = 50.0; index > 0.0; --index) {
    result.xyz = eye + ray * total;
    float dist = map(result.xyz);
    if (dist < 0.001 + total * 1./synth_Resolution.y) {
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

vec3 dots (float levelOfDetails) {
  vec2 uv = gl_FragCoord.xy/synth_Resolution;
  vec2 lod = synth_Resolution.xy/levelOfDetails;
  vec2 uvpixel = floor(uv * lod + .5) / lod;
  vec4 frame = texture2D(b0, uvpixel);
  vec2 cell = 1./lod;
  uv = repeat(uv+cell/2., cell);
  float angle = random(uvpixel) * PI * 2.;
  uv += vec2(cos(angle),sin(angle)) * cell.y/6.;
  uv.x *= synth_Resolution.x/synth_Resolution.y;
  float shape = smoothstep(cell.y/200., 0., length(uv)-cell.y/3.);
  return clamp(frame.rgb * shape, 0., 1.);
}

float getShadow (vec3 pos, vec3 at, float k) {
    vec3 dir = normalize(at-pos);
    float maxt = 1.;
    float f = 1.;
    float t = .001;
    for (int i = 0; i < 60; ++i) {
        float dist = map(pos + dir * t);
        if (dist < .001) return 0.;
        f = min(f, k * dist / t);
        t += dist;
        if (t >= maxt) break;
    }
    return f;
}

vec3 shade (vec3 view, vec4 pos) {
  float ao = pos.a;
  vec3 color = vec3(0);
  vec3 normal = getNormal(pos.xyz);
  int stage = getStage();
  float s0 = floor(beat);
  float b0 = fract(beat);
  if (stage == STAGE_RING) {

    color = vec3(0.980, 0.729, 0.478);
    color += vec3(0.980, 0.533, 0.478) * clamp(dot(normal, normalize(vec3(1,1,-1)))*.5+.5, 0., 1.);
    color += vec3(0.980, 0.078, 0.149) * pow(clamp(dot(normal, normalize(vec3(-2,-2,-4))), 0., 1.), 4.);

  } else if (stage == STAGE_CITY) {

    vec3 light = normalize(vec3(0.,1.,0.));
    // intensity = clamp(intensity, 0., 1.);
    // intensity *= b0;
    color += vec3(1) * abs(dot(normal, light));
    color += vec3(1) * pow(clamp(dot(normal, normalize(vec3(-1,-.1,-2))), 0., 1.), 2.);
    color += vec3(1) * pow(clamp(dot(normal, normalize(vec3(1,.1,-2))), 0., 1.), 2.);
    // color = 1.-color;
    // color = smoothstep(.9, 1., color);
    // color.r += .01/(abs(cos(pos.y-beat)));
    // color.gb += .01/(1.-max(.0,sin((pos.z-beat)*2.)));

  } else {

    color = vec3(.1);
    color += vec3(0.760, 0.925, 1) * clamp(dot(normal, normalize(vec3(1,2,-1))), 0., 1.);
    color += vec3(1, 0.917, 0.760) * pow(clamp(dot(normal, normalize(vec3(-2,-4,-1))), 0., 1.), 2.);
    color += vec3(0.925, 1, 0.760) * pow(clamp(dot(normal, normalize(vec3(-4,0,0))), 0., 1.), 4.);

  }
  color *= ao;
  return color;
}

void main() {

  if (PASSINDEX == 0) {

    vec2 uv = (gl_FragCoord.xy-0.5*synth_Resolution)/synth_Resolution.y;
    vec3 eye = vec3(-.01,-.01,-4.);
    // float angle = floor(beat) * 2.5465 * PI;
    // eye.xy += vec2(cos(angle), sin(angle)) * .2;
    vec3 at = vec3(0);
    // at.xy -= vec2(cos(angle), sin(angle)) * .1;
    vec3 view = look(eye, at, uv);
    vec4 result = raymarch(eye, view);
    vec3 color = shade(view, result);
    // vec3 color = anaglyph(eye, at, uv);

    // gl_FragColor = vec4(color, 1.);
    vec4 frame = texture2D(b0, gl_FragCoord.xy/synth_Resolution);
    float blend = .0;//.7 * (1.-mod(beat, 1.));
    gl_FragColor = frame*blend + (1.-blend)*vec4(color, 1);
    // gl_FragColor = vec4(color, 1);

  } else if (PASSINDEX == 1) {

    vec2 uv = gl_FragCoord.xy / synth_Resolution;
    // vec4 frame = texture2D(b1, uv);
    gl_FragColor = texture2D(b0, uv);

  }

}


//! <preset file="shader.preset" />
