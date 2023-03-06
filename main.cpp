#include "stdlib.h" // for calloc

#define STB_IMAGE_WRITE_IMPLEMENTATION
#pragma warning(disable:4996)
#include "stb_image_write.h"

typedef signed char        i8;
typedef short              i16;
typedef int                i32;
typedef long long          i64;
typedef unsigned char      u8;
typedef unsigned short     u16;
typedef unsigned int       u32;
typedef unsigned long long u64;

typedef float f32;
typedef double f64;

struct v2 { f32 x, y; };
struct v3 { f32 x, y, z; };
struct v4 { f32 x, y, z, w; };
struct tri { v3 vert[3]; };

v3 operator+(v3 a, v3 b) { return v3{ a.x + b.x, a.y + b.y, a.z + b.z }; }
v3 operator-(v3 a, v3 b) { return v3{ a.x - b.x, a.y - b.y, a.z - b.z }; }
v3 operator*(v3 a, v3 b) { return v3{ a.x * b.x, a.y * b.y, a.z * b.z }; }
v3 operator/(v3 a, v3 b) { return v3{ a.x / b.x, a.y / b.y, a.z / b.z }; }

v3 operator+(v3 a, f32 b) { return v3{ a.x + b, a.y + b, a.z + b }; }
v3 operator-(v3 a, f32 b) { return v3{ a.x - b, a.y - b, a.z - b }; }
v3 operator*(v3 a, f32 b) { return v3{ a.x * b, a.y * b, a.z * b }; }
v3 operator/(v3 a, f32 b) { return v3{ a.x / b, a.y / b, a.z / b }; }

v3 operator+(f32 a, v3 b) { return b + a; }
v3 operator-(f32 a, v3 b) { return b - a; }
v3 operator*(f32 a, v3 b) { return b * a; }
v3 operator/(f32 a, v3 b) { return b - a; }

v3 cross(v3 a, v3 b)
{
	return
	{
		a.y * b.z - a.z * b.y,
		-(a.x * b.z - a.z * b.x),
		a.x * b.y - a.y * b.x
	};
}

f32 dot(v3 a, v3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

f32 absf(f32 x) { return x > -x ? x : -x; }

f32 rcp_sqrt(f32 number) // infamous 'fast reciprocal square root' function
{
	i32 i = 0x5f3759df - ((*(i32*)&number) >> 1);
	f32 y = *(f32*)&i;
	y = y * (1.5f - (number * 0.5f * y * y)); // 1st iteration
	y = y * (1.5f - (number * 0.5f * y * y)); // 2nd iteration
	return y;
}

v3 normalize(v3 V) { return V * rcp_sqrt(dot(V, V)); }

struct visbuf
{
	i32 width, height;
	f32* z_buf;
	u32* prim_buf;
};

visbuf make_visbuf(i32 width, i32 height)
{
	visbuf result;
	result.width = width;
	result.height = height;
	result.prim_buf = (u32*)calloc(width * height, sizeof(u32));
	result.z_buf = (f32*)calloc(width * height, sizeof(f32));
	return result;
}

// "<=" seems to work better than "<" for whatever reason
// We get overlap with "<=" but get jagged edge artifacts with "<" along shared edges crossing pixel centers
// rearranges the expression:   z <= decode(vb.z_buf[i * vb.width + j])
f32 compare_z(f32 z, i32 i, i32 j, const visbuf& vb) { return (z + 1.f) * vb.z_buf[i * vb.width + j] <= 1.f; } 

f32 encode_z(f32 z) { return 1.f/(1.f + z);  }
f32 decode_z(f32 z) { return (1.f/z) - 1.f; }

v3 get_viewdir(i32 i, i32 j, i32 width, i32 height)
{
	return {
		(f32)(1 - width + 2 * j) / (f32)width,
		-(f32)(1 - height + 2 * i) / (f32)height,
		-1.f
	};
}

void rasterize(i32 i, i32 j, u32 k, tri* prims, visbuf& vb)
{
	v3 D = get_viewdir(i, j, vb.width, vb.height);
	tri prim = prims[k];
	v3 V0 = prim.vert[0];
	v3 V1 = prim.vert[1];
	v3 V2 = prim.vert[2];

	v3 N = cross(V1 - V0, V2 - V0); // Must calculate the normal before modifying the vertex positions!
	f32 NdotV0 = dot(N, V0);

	// case 0: if(NdotV0 < 0 && DdotN < 0) // in front of camera pointing towards it ( frontfacing ) 
	// case 1: if(NdotV0 < 0 && DdotN > 0) // behind camera pointed towards camera 
	// case 2: if(NdotV0 > 0 && DdotN > 0) // in front of camera pointed away (backfacing)
	// case 3: if(NdotV0 > 0 && DdotN < 0) // behind camera pointed away from camera

	if (NdotV0 > 0.f) return; // excludes case 2 and case 3

	//V0 = V0 / absf(prim.vert[0].z);	// Assume z is negative, so -z == abs(z), for now!
	//V1 = V1 / absf(prim.vert[1].z);	// Assume z is negative, so -z == abs(z), for now!
	//V2 = V2 / absf(prim.vert[2].z);	// Assume z is negative, so -z == abs(z), for now!

	v3 V01 = cross(V0, V1);
	v3 V12 = cross(V1, V2);
	v3 V20 = cross(V2, V0);

	f32 e01 = -dot(V01, D);
	f32 e12 = -dot(V12, D);
	f32 e20 = -dot(V20, D);

	if (!(e01 >= 0.f && e12 >= 0.f && e20 >= 0.f)) return;

	f32 DdotN = dot(D, N); 

	// For camerapos C = 0, triangle normal N, view dir D
	// Solve dot(N, C-v0) + t * dot(D, N) == 0, for t

	// t < 0 implies case 1, therefore, require that t > 0
	// The expression for t can be flipped to immediately get reciprocal z.
	// To immediately get 1/(1+z) from z = z_nom/z_den: z_den/(z_den+z_nom)
	f32 t = NdotV0 / DdotN;
	f32 z = t * -D.z;
	if (t > 0.f 
		&& 
		//z < decode_z(vb.z_buf[i*vb.width + j])
		compare_z(z, i, j, vb)
		)
	{
		vb.z_buf[i * vb.width + j] = encode_z(z);
		vb.prim_buf[i * vb.width + j] = k+1;
	}
}

struct ShadingData
{
	v3 barycentric;
	v3 pos;
	v3 viewdir;
	v3 normaldir;
};
ShadingData get_shading_data(i32 i, i32 j, tri* prims, const visbuf& vb)
{
	u32 primidx = vb.prim_buf[i * vb.width + j] - 1;
	tri prim = prims[primidx];
	v3 D = get_viewdir(i, j, vb.width, vb.height);
	f32 Z = decode_z(vb.z_buf[i * vb.width + j]);
	v3 pos = D * Z;
	v3 N = normalize(cross(prim.vert[1] - prim.vert[0], prim.vert[2] - prim.vert[0]));
	ShadingData result;
	result.barycentric = v3{};
	result.pos = pos;
	result.viewdir = normalize(D) * -1.f;
	result.normaldir = N;
	return result;
}

f32 gamma_encode(f32 Clinear)
{
	Clinear >= 0.f ? Clinear : 0.f;
	Clinear <= 1.f ? Clinear : 1.f;
	return 1. / rcp_sqrt(Clinear);
}

u32 HDRtoLDR(v3 RGB)
{
	u32 R = (u8)(gamma_encode(RGB.x) * 255.f);
	u32 G = (u8)(gamma_encode(RGB.y) * 255.f);
	u32 B = (u8)(gamma_encode(RGB.z) * 255.f);
	return 0xFF000000 | B << 16 | G << 8 | R;
}

u32 shade(i32 i, i32 j, tri* prims, const visbuf& vb)
{
	u32 primbufval = vb.prim_buf[i * vb.width + j]; // note: unsigned value, be careful with subtraction!

	v3 rgb_palette[10] =
	{
		v3{1.f, 1.f, 1.f},
		v3{1.f, 1.f, 1.f},
		v3{0.f, 1.f, 0.f},
		v3{0.f, 1.f, 0.f},
		v3{1.f, 1.f, 1.f},
		v3{1.f, 1.f, 1.f},
		v3{1.f, 1.f, 1.f},
		v3{1.f, 1.f, 1.f},
		v3{1.f, 0.f, 0.f},
		v3{1.f, 0.f, 0.f}
	};

	if (primbufval)
	{
		u32 primidx = primbufval - 1;
		tri prim = prims[primidx];
		ShadingData s = get_shading_data(i, j, prims, vb);
		f32 NdotV = dot(s.normaldir, s.viewdir);
		f32 diffuse = NdotV >= 0.f ? NdotV : 0.f;
		v3 baseColor = rgb_palette[primidx];
		v3 color = v3{ diffuse,diffuse,diffuse } *baseColor;
		u32 LDR = HDRtoLDR(color);
		return LDR;
	}
	else
		return 0xFF330011;// 0xFF000000;
}

tri transform(u32 k, tri* prim)
{
	tri result = tri{{
		prim[k].vert[0] * 0.5 - v3{0.f, 0.f, 0.5f},
		prim[k].vert[1] * 0.5 - v3{0.f, 0.f, 0.5f},
		prim[k].vert[2] * 0.5 - v3{0.f, 0.f, 0.5f}
	}};
	return result;
}

i32 main(i32 argc, char* argv[])
{
	const i32 width = 600, height = 600;
	const i32 num_prims = 10;
	v3 v[8]
	{
		{-1.f, -1.f, -1.f},
		{ 1.f, -1.f, -1.f},
		{ 1.f,  1.f, -1.f},
		{-1.f,  1.f, -1.f},
		{-1.f, -1.f, -2.f},
		{ 1.f, -1.f, -2.f},
		{ 1.f,  1.f, -2.f},
		{-1.f,  1.f, -2.f},
	};
	tri prims[num_prims] =
	{
		{ { v[0], v[1], v[5] }},
		{ { v[0], v[5], v[4] }},
		{ { v[1], v[2], v[5] }},
		{ { v[2], v[6], v[5] }},
		{ { v[2], v[7], v[6] }},
		{ { v[2], v[3], v[7] }},
		{ { v[4], v[5], v[6] }},
		{ { v[4], v[6], v[7] }},
		{ { v[0], v[4], v[3] }},
		{ { v[3], v[4], v[7] }}
	};

	tri transformed_prims[num_prims];

	u32* img_out = (u32*)calloc(width * height, sizeof(u32));
	if (!img_out) return 0;

	visbuf vb = make_visbuf(width, height);

	for (u32 k = 0; k < num_prims; k++)
		transformed_prims[k] = transform(k, prims); // transform per triangle, should transform per vertex instead...

	for (i32 i = 0; i < height; ++i)
		for (i32 j = 0; j < width; ++j)
			for (u32 k = 0; k < num_prims; k++)				
				rasterize(i, j, k, transformed_prims, vb);

	for (i32 i = 0; i < height; ++i)
		for (i32 j = 0; j < width; ++j)
			img_out[i*width + j] = shade(i, j, transformed_prims, vb);

	stbi_write_bmp("imgout.bmp", width, height, 4, img_out);
}