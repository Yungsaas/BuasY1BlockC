#pragma once
#include "stb_image.h"

constexpr int MAX_RECURSION_DEPTH = 4;

#define EPSILON 0.001f
#define AIR_FRESNEL 1.0f

namespace Tmpl8
{

	class Renderer : public TheApp
	{
	public:
		// game flow methods
		void Init();
		float3 Trace(Ray& ray, int _recursionDepth);
		// anti aliasing
		float3 TracePixel(float x, float y, int _recursionDepth, int numSamples);

		void Tick(float deltaTime);
		void UI();
		void Shutdown();

		// light calc
		void PointLightHandling(Ray& ray, PointLight* light, float3& I, const float3 L, float3& albedo, float3& N, float3& shadingColor);
		void SpotLightHandling(Ray& ray, SpotLight* light, float3& I, const float3 L, float3& albedo, float3& N, float3& shadingColor);
		void AreaLightHandling(Ray& ray, AreaLight* light, const float3& I, const float3& albedo, const float3& N, float3& shadingColor);

		// input handling
		void MouseUp(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseDown(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
		void MouseWheel(float y) { /* implement if you want to handle the mouse wheel */ }
		void KeyUp(int key) { /* implement if you want to handle keys */ }
		void KeyDown(int key) { /* implement if you want to handle keys */ }

		// reflections, made with the help of chatgpt
		Ray CreateReflectionRay(const float3& reflectionPoint, const float3& reflectionDirection, const float3& surfaceNormal, const float reflectionCoefficient, const int recursionDepth);

		// refractions, made with the help of chatgpt
		float3 Refract(const float3& incidentDir, const float3& normal, float n2)
		{
			// Refractive index of air
			const float n1 = 1.0f;

			float cosTheta1 = dot(-incidentDir, normal);
			float eta = n1 / n2;
			float sinTheta2Sq = eta * eta * (1.0f - cosTheta1 * cosTheta1);

			// Check for total internal reflection
			if (sinTheta2Sq >= 1.0f)
			{
				// Total internal reflection, no refracted ray
				return float3(0.0f);
			}

			float cosTheta2 = sqrt(1.0f - sinTheta2Sq);
			return eta * incidentDir + (eta * cosTheta1 - cosTheta2) * normal;
		}

	// data members
	int2 mousePos;
	float4* accumulator;
	Scene scene;
	Camera camera;
	int elapsedFrameCounter = 0;
	int numSamples = 1;

	//HDR sky
	const char* HDR_IMAGE_PATH = "assets/neon_photostudio_4k.hdr";
	int skyBpp = 0;
	int skyWidth, skyHeight;
	float* skyPixels = stbi_loadf(HDR_IMAGE_PATH, &skyWidth, &skyHeight, &skyBpp, 0);
};

} // namespace Tmpl8