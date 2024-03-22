#pragma once

// high level settings
// #define TWOLEVEL
#define WORLDSIZE 128 // power of 2. Warning: max 512 for a 512x512x512x4 bytes = 512MB world!
// #define USE_SIMD
// #define USE_FMA3
// #define SKYDOME
// #define WHITTED
// #define DOF

// low-level / derived
#define WORLDSIZE2	(WORLDSIZE*WORLDSIZE)
#ifdef TWOLEVEL
// feel free to replace with whatever suits your two-level implementation,
// should you chose this challenge.
#define BRICKSIZE	8
#define BRICKSIZE2	(BRICKSIZE*BRICKSIZE)
#define BRICKSIZE3	(BRICKSIZE*BRICKSIZE*BRICKSIZE)
#define GRIDSIZE	(WORLDSIZE/BRICKSIZE)
#define VOXELSIZE	(1.0f/WORLDSIZE)
#else
#define GRIDSIZE	WORLDSIZE
#endif
#define GRIDSIZE2	(GRIDSIZE*GRIDSIZE)
#define GRIDSIZE3	(GRIDSIZE*GRIDSIZE*GRIDSIZE)

// define maximum number of things in the scene
#define MAX_POINT_LIGHTS 10
#define MAX_SPOT_LIGHTS 10
#define MAX_AREA_LIGHTS 10
#define MAX_DIRECTIONAL_LIGHTS 10
#define MAX_SPHERES 10
#define MAX_PLANES 4

namespace Tmpl8 {

	//Structs --- new
	struct Material
	{
		//Color components
		float3 albedo;
		float3 diffuse;
		float3 specular;

		//Other material components
		float roughness;
		float reflectionCoefficient;
		float refractionCoefficient;

		//More information
		bool isTransparent;

		Material() = default;
		Material(float3 alb, float3 diff, float3 spec, float rough, float refl, float refr, bool isTransp)
			:albedo(alb), diffuse(diff), specular(spec), roughness(rough), reflectionCoefficient(refl), refractionCoefficient(refr), isTransparent(isTransp) {}
	};

	

class Ray
{
public:
	Ray() = default;
	Ray( const float3 origin, const float3 direction, const float rayLength = 1e34f, const int rgb = 0 )
		: O( origin ), D( direction ), t( rayLength ), voxel( rgb )
	{
		// calculate reciprocal ray direction for triangles and AABBs
		// TODO: prevent NaNs - or don't
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
		Dsign = (float3( -copysign( 1.0f, D.x ), -copysign( 1.0f, D.y ), -copysign( 1.0f, D.z ) ) + 1) * 0.5f;
	}

	// reflection, made with the help of chatgpt
	float reflectionCoefficient;  // Reflectivity of the material (e.g., 0.0 for non-reflective, 1.0 for fully reflective)
	int recursionDepth = 0;           // Current recursion depth of the ray

	float3 refractedDirection; // Store the refracted ray direction

	float3 IntersectionPoint() const { return O + t * D; }
	float3 GetNormal() const;
	float3 GetAlbedo() const;
	float GetReflectivity() const;
	uint GetMaterial() const;
	float GetRefractivity() const;
	float3 GetAbsorption( const float3& I ) const; // TODO: implement

	// ray data
	float3 O;					// ray origin
	float3 rD;					// reciprocal ray direction
	float3 D = float3( 0 );		// ray direction
	float t = 1e34f;			// ray length
	float3 Dsign = float3( 1 );	// inverted ray direction signs, -1 or 1
	uint voxel = 0;				// 32-bit ARGB color of a voxelhit object index; 0 = NONE
	float3 normalAtIntersection = float3(0);

	Material rayMaterial;


private:
	// min3 is used in normal reconstruction.
	__inline static float3 min3( const float3& a, const float3& b )
	{
		return float3( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
	}
};

class Sphere
{
public:
	float3 position;
	float radius;
	Material sphereMaterial;

	Sphere() = default;
	Sphere(const float3& pos, const float rad, Material mat)
		: position(pos), radius(rad), sphereMaterial(mat) {}

	float Intersect(const Ray& ray) const; // Intersection point t returns
	float3 GetNormal(float3& I );
};

class Plane
{
public:
	float3 position;
	float3 normal;
	float3 size;
	Material planeMaterial;

	Plane() = default;

	float Intersect(const Ray& ray) const;
};

class Cube
{
public:
	Cube() = default;
	Cube( const float3 pos, const float3 size );
	float Intersect( const Ray& ray ) const;
	bool Contains( const float3& pos ) const;
	float3 b[2];
};

struct PointLight
{
public:
	float3 pos;
	float3 color;
};

struct SpotLight
{
public:
	float3 pos;
	float3 color;
	float3 direction;
	float cosCutoff;

};

class AreaLight
{
public:
	float3 position;
	float3 normal;
	float3 color;
	float area;
};

class DirectionalLight
{
	float3 normal;
	float3 color;
};

class Scene
{
public:
	struct DDAState
	{
		int3 step;				// 16 bytes
		uint X, Y, Z;			// 12 bytes
		float t;				// 4 bytes
		float3 tdelta;
		float dummy1 = 0;		// 16 bytes
		float3 tmax;
		float dummy2 = 0;		// 16 bytes, 64 bytes in total
	};

	Scene();
	void FindNearest( Ray& ray ) const;
	bool IsOccluded(  Ray& ray, float3 lightPos ) const;
	void Set( const uint x, const uint y, const uint z, const uint v );
	unsigned int* grid;
	Cube cube;

	float3 globalIllumination = float3(0.1f, 0.1f, 0.1f);

	Material standartPlaneMaterial = Material(float3(1.0f, 1.0f, 0.5f), float3(1.0f, 1.0f, 1.0f), float3(1.0f, 1.0f, 1.0f), 1.0f, 0.0f, 0.0f, false);
	int currentPlanes = 0;
	Plane* planeArray[MAX_PLANES];

	Material standartSphereMaterial = Material(float3(0.0f, 0.0f, 0.0f), float3(1.0f, 1.0f, 1.0f), float3(1.0f,1.0f,1.0f), 1.0f,  0.0f,  1.2f, true);
	int currentSpheres = 0;
	Sphere* sphereArray[MAX_SPHERES];


	//implement for better DOD later
	int currentPointLights = 0;
	PointLight* pointLightArray[MAX_POINT_LIGHTS];

	int currentSpotLights = 0;
	SpotLight* spotLightArray[MAX_SPOT_LIGHTS];

	int currentAreaLights = 0;
	AreaLight* areaLightArray[MAX_AREA_LIGHTS];

	int currentDirectionalLights = 0;
	DirectionalLight* directionLightArray[MAX_DIRECTIONAL_LIGHTS];

	//vectors for now
	std::vector<PointLight*> pointLightVector;
	std::vector<SpotLight*> spotLightVector;
	std::vector<AreaLight*> areaLightVector;
	std::vector<DirectionalLight*> directionalLightVector;

private:
	bool Setup3DDDA( const Ray& ray, DDAState& state ) const;
};

}