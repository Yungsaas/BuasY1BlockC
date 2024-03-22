#include "precomp.h"

// YOU GET:
// 1. A fast voxel renderer in plain C/C++
// 2. Normals and voxel colors
// FROM HERE, TASKS COULD BE:							FOR SUFFICIENT
// * Materials:
//   - Reflections and diffuse reflections				<=== Done
//   - Transmission with Snell, Fresnel					<=== TODO
//   - Textures, Minecraft-style						<=== TODO
//   - Beer's Law
//   - Normal maps
//   - Emissive materials with postproc bloom
//   - Glossy reflections (BASIC)
//   - Glossy reflections (microfacet)
// * Light transport:
//   - Point lights										<=== Done
//   - Spot lights										<=== Done
//   - Area lights										<=== Done
//	 - Sampling multiple lights with 1 ray
//   - Importance-sampling
//   - Image based lighting: sky
// * Camera:
//   - Depth of field									<=== TODO
//   - Anti-aliasing									<=== Done
//   - Panini, fish-eye etc.
//   - Post-processing: now also chromatic				<=== TODO
//   - Spline cam, follow cam, fixed look-at cam
//   - Low-res cam with CRT shader
// * Scene:
//   - HDR skydome										<=== Done
//   - Spheres											<=== Done
//   - Smoke & trilinear interpolation
//   - Signed Distance Fields
//   - Voxel instances with transform
//   - Triangle meshes (with a BVH)
//   - High-res: nested grid
//   - Procedural art: shapes & colors
//   - Multi-threaded Perlin / Voronoi
// * Various:
//   - Object picking
//   - Ray-traced physics
//   - Profiling & optimization
// * GPU:
//   - GPU-side Perlin / Voronoi
//   - GPU rendering *not* allowed!
// * Advanced:
//   - Ambient occlusion
//   - Denoising for soft shadows
//   - Reprojection for AO / soft shadows
//   - Line lights, tube lights, ...
//   - Bilinear interpolation and MIP-mapping
// * Simple game:										
//   - 3D Arkanoid										<=== TODO
//   - 3D Snake?
//   - 3D Tank Wars for two players
//   - Chess
// REFERENCE IMAGES:
// https://www.rockpapershotgun.com/minecraft-ray-tracing
// https://assetsio.reedpopcdn.com/javaw_2019_04_20_23_52_16_879.png
// https://www.pcworld.com/wp-content/uploads/2023/04/618525e8fa47b149230.56951356-imagination-island-1-on-100838323-orig.jpg

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( accumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );
	// try to load a camera
	FILE* f = fopen( "camera.bin", "rb" );
	if (f)
	{
		fread( &camera, 1, sizeof( Camera ), f );
		fclose( f );
	}
	for (int i = 0; i < skyWidth * skyHeight * 3; i++) skyPixels[i] = sqrtf(skyPixels[i]);
}

// -----------------------------------------------------------
// Evaluate light transport -> made with help from chatGPT
// -----------------------------------------------------------
float3 Renderer::Trace( Ray& ray, int _recursionDepth )
{
	ray.recursionDepth = _recursionDepth;
	if (ray.recursionDepth >= MAX_RECURSION_DEPTH) return float3(0);

	scene.FindNearest( ray );
	float3 normalizedDirection = normalize(ray.D);

	if (ray.t == 1e34f) // if nothing is hit, return the hdr
	{

		float uSky = skyWidth * atan2f(normalizedDirection.z, normalizedDirection.x) * INV2PI - 0.5f;
		float vSky = skyHeight * acosf(normalizedDirection.y) * INVPI - 0.5f;

		uint u = static_cast<int>(uSky);
		uint v = static_cast<int>(vSky);

		uint skyIdx = u + v * skyWidth;
		if (skyIdx < 0 || skyIdx > uint(skyHeight*skyWidth) ) {
			skyIdx = 0;
		}

		return float3(skyPixels[skyIdx * 3], skyPixels[skyIdx * 3 + 1], skyPixels[skyIdx * 3 + 2]);
	}

	uint material = ray.GetMaterial();
	float3 I = ray.O + ray.t * ray.D; // <- this is the Intersection Point
	static const float3 L = normalize( float3( 1, 4, 0.5f ) );
	float3 albedo = ray.GetAlbedo();
	float3 N = ray.GetNormal();
	float3 shadingColor = float3(0.0f, 0.0f, 0.0f);

	// If the material is transparent
	if (ray.rayMaterial.isTransparent)
	{
		//// Compute reflection ray
		float3 reflectionDirection = normalize(reflect(ray.D, N));

		Ray reflectionRay(I + EPSILON * reflectionDirection, reflectionDirection);

		// Compute refraction ray
		float3 refractedDirection = Refract(ray.D, N, ray.GetRefractivity());
		Ray refractionRay(I + EPSILON * refractedDirection, refractedDirection);

		// Trace reflection and refraction rays recursively
		float3 reflectionColor = Trace(reflectionRay, _recursionDepth + 1);
		float3 refractionColor = Trace(refractionRay, _recursionDepth + 1);

		// Compute Fresnel term
		float R0 = pow((1.0f - ray.GetRefractivity()) / (1.0f + ray.GetRefractivity()), 2);
		float fresnel = R0 + (1.0f - R0) * pow(1.0f - dot(-ray.D, N), 5);

		// Combine reflection and refraction colors using Fresnel term
		return fresnel * reflectionColor + (1.0f - fresnel) * refractionColor;
	}


	// Reflection handling
	if (_recursionDepth < MAX_RECURSION_DEPTH)
	{
		float reflectionCoefficient = ray.rayMaterial.reflectionCoefficient;

		// Direct Reflection
		if (reflectionCoefficient > 0.0f)
		{
			Ray reflectionRay = CreateReflectionRay(I, ray.D, N, reflectionCoefficient, _recursionDepth);
			shadingColor += reflectionCoefficient * Trace(reflectionRay, _recursionDepth + 1);
		}
	}

	if (ray.GetReflectivity() == 0 && ray.GetRefractivity() == 0) // <---- only calculate light if...
	{
		// Loop through the point light vector
		for each (PointLight * light in scene.pointLightVector)
		{
			PointLightHandling(ray, light, I, L, albedo, N, shadingColor); // <-- this was made with the help of ChatGPT
		}

		// Normalize the shading color
		shadingColor /= std::max(1u, static_cast<unsigned int>(scene.pointLightVector.size()));

		// Loop through the spotlight vector
		for each (SpotLight * light in scene.spotLightVector)
		{
			SpotLightHandling(ray, light, I, L, albedo, N, shadingColor); // <-- this was made with the help of ChatGPT
		}

		// Normalize the shading color
		shadingColor /= std::max(1u, static_cast<unsigned int>(scene.spotLightVector.size()));

		// Loop through the area light vector
		for each (AreaLight * light in scene.areaLightVector)
		{
			AreaLightHandling(ray, light, I, albedo, N, shadingColor);
		}

		// Normalize the shading color
		shadingColor /= std::max(1u, static_cast<unsigned int>(scene.areaLightVector.size()));
	}

	float3 ambientColor = float3(0);
	if(ray.GetReflectivity() < 1.0f)
	{ 
	ambientColor = scene.globalIllumination * albedo; // ambient Lighting
	}

	float3 finalColor = shadingColor + ambientColor;

	if (finalColor.x >= 1.0f)
	{
		finalColor.x = 1.0f;
	}
	if (finalColor.y >= 1.0f)
	{
		finalColor.y = 1.0f;
	}
	if (finalColor.z >= 1.0f)
	{
		finalColor.z = 1.0f;
	}

	return finalColor;
}


// Shooting multiple rays to the same pixel for anti aliasing
float3 Renderer::TracePixel(float x, float y, int _recursionDepth, int numSamples)
{
	float3 totalColor = float3(0.0f, 0.0f, 0.0f);

	for (int s = 0; s < numSamples; ++s)
	{
		// Calculate jittered coordinates within the pixel

		float jitterX = Rand((0.9f / float(numSamples)) * (float(s) + 1.0f));// calculate jittered x-coordinate within the pixel
		float jitterY = Rand((0.9f / float(numSamples)) * (float(s) + 1.0f));// calculate jittered y-coordinate within the pixel

		// Trace the ray for the jittered coordinates
		float3 pixel = Trace(camera.GetPrimaryRay(x + jitterX, y + jitterY), 0);

		// Accumulate the color
		totalColor += pixel;
	}

	// Average the color
	return totalColor / numSamples;
}



// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick( float deltaTime )
{
	// pixel loop
	Timer t;
	elapsedFrameCounter++;

	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)

	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			float4 pixel = float4(TracePixel((float)x, (float)y, 0, numSamples), 0);

			// translate accumulator contents to rgb32 pixels
			screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8( &pixel );
			accumulator[x + y * SCRWIDTH] += pixel;
			const float4 outPixel = accumulator[x + y * SCRWIDTH] / float(elapsedFrameCounter);
			screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&outPixel);
		}
	}

	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf( "%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000 );
	// handle user input
	if (camera.HandleInput(deltaTime))
	{
		memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * 16);
		elapsedFrameCounter = 0;
	}
}

// -----------------------------------------------------------
// Update user interface (imgui) -> this code was written with the help of chatGPT
// -----------------------------------------------------------
void Renderer::UI()
{
	static int selectedPointLightIndex = -1; // Variable to store the selected point light index
	static int selectedSpotLightIndex = -1; // Variable to store the selected spot light index
	static int selectedAreaLightIndex = -1; // Variable to store the selected are light index

	// Menu
	ImGui::Begin("Menu");

	// Dropdown for selecting point lights in the scene
	if (ImGui::BeginCombo("Point Light List", (selectedPointLightIndex != -1) ? std::to_string(selectedPointLightIndex + 1).c_str() : "Select Point Light"))
	{
		// Add Point Light button inside the combo box
		if (ImGui::Selectable("Add Point Light"))
		{
			// Add a new PointLight to the vector
			PointLight* newLight = new PointLight();
			scene.pointLightVector.push_back(newLight);
			selectedPointLightIndex = scene.pointLightVector.size() - 1; // Select the newly added point light
			selectedSpotLightIndex = -1; // Reset selected spot light index
		}

		for (int i = 0; i < scene.pointLightVector.size(); i++)
		{
			bool isSelected = (selectedPointLightIndex == i);
			if (ImGui::Selectable(("Point Light " + std::to_string(i + 1)).c_str(), isSelected))
			{
				selectedPointLightIndex = i;
				selectedSpotLightIndex = -1; // Reset selected spot light index
			}
			if (isSelected)
				ImGui::SetItemDefaultFocus();
		}

		ImGui::EndCombo();
	}

	// Dropdown for selecting spotlights in the scene
	if (ImGui::BeginCombo("Spot Light List", (selectedSpotLightIndex != -1) ? std::to_string(selectedSpotLightIndex + 1).c_str() : "Select Spot Light"))
	{
		// Add Spot Light button inside the combo box
		if (ImGui::Selectable("Add Spot Light"))
		{
			// Add a new SpotLight to the vector
			SpotLight* newLight = new SpotLight();
			scene.spotLightVector.push_back(newLight);
			selectedSpotLightIndex = scene.spotLightVector.size() - 1; // Select the newly added spot light
			selectedPointLightIndex = -1; // Reset selected point light index
		}

		for (int i = 0; i < scene.spotLightVector.size(); i++)
		{
			bool isSelected = (selectedSpotLightIndex == i);
			if (ImGui::Selectable(("Spot Light " + std::to_string(i + 1)).c_str(), isSelected))
			{
				selectedSpotLightIndex = i;
				selectedPointLightIndex = -1; // Reset selected point light index
			}
			if (isSelected)
				ImGui::SetItemDefaultFocus();
		}

		ImGui::EndCombo();
	}

	// Dropdown for selecting area lights in the scene
	if (ImGui::BeginCombo("Area Light List", (selectedAreaLightIndex != -1) ? std::to_string(selectedAreaLightIndex + 1).c_str() : "Select Area Light"))
	{
		// Add Area Light button inside the combo box
		if (ImGui::Selectable("Add Area Light"))
		{
			// Add a new AreaLight to the vector
			AreaLight* newLight = new AreaLight();
			scene.areaLightVector.push_back(newLight);
			selectedAreaLightIndex = scene.areaLightVector.size() - 1; // Select the newly added area light
		}

		for (int i = 0; i < scene.areaLightVector.size(); i++)
		{
			bool isSelected = (selectedAreaLightIndex == i);
			if (ImGui::Selectable(("Area Light " + std::to_string(i + 1)).c_str(), isSelected))
			{
				selectedAreaLightIndex = i;
			}
			if (isSelected)
				ImGui::SetItemDefaultFocus();
		}

		ImGui::EndCombo();
	}

	// Display information or edit the selected area light in a new window
	if (selectedAreaLightIndex != -1)
	{
		ImGui::Begin("Area Light Editor");

		// Delete Area Light button
		if (ImGui::Button("Delete Area Light"))
		{
			// Delete the selected area light and update the index
			delete scene.areaLightVector[selectedAreaLightIndex];
			scene.areaLightVector.erase(scene.areaLightVector.begin() + selectedAreaLightIndex);

			// Adjust selected index
			selectedAreaLightIndex = -1; // No area light selected after deletion
		}

		// Close button 
		ImGui::SameLine();
		if (ImGui::Button("Close"))
		{
			selectedAreaLightIndex = -1; // No area light selected
		}

		if (selectedAreaLightIndex != -1) // Avoid trying to access lights that don't exist
		{
			// Editing the selected area light
			ImGui::Text("Editing Area Light: %s", ("Area Light " + std::to_string(selectedAreaLightIndex + 1)).c_str());
			ImGui::SliderFloat3("Position", &scene.areaLightVector[selectedAreaLightIndex]->position.x, -3.0f, 3.0f);
			ImGui::SliderFloat3("Normal", &scene.areaLightVector[selectedAreaLightIndex]->normal.x, -1.0f, 1.0f);
			ImGui::ColorEdit3("Color", &scene.areaLightVector[selectedAreaLightIndex]->color.x);
			ImGui::SliderFloat("Area", &scene.areaLightVector[selectedAreaLightIndex]->area, 0.1f, 10.0f);
		}

		ImGui::End(); // Close "Area Light Editor" window
	}

	ImGui::End(); // Close "Menu" window

	// Display information or edit the selected point light in a new window
	if (selectedPointLightIndex != -1)
	{
		ImGui::Begin("Point Light Editor");

		// Delete Point Light button
		if (ImGui::Button("Delete Point Light"))
		{
			// Delete the selected point light and update the index
			delete scene.pointLightVector[selectedPointLightIndex];
			scene.pointLightVector.erase(scene.pointLightVector.begin() + selectedPointLightIndex);

			// Adjust selected indices
			selectedPointLightIndex = -1; // No point light selected after deletion
		}

		// Close button 
		ImGui::SameLine();
		if (ImGui::Button("Close"))
		{
			selectedSpotLightIndex = -1;
			selectedPointLightIndex = -1; // No point light selected
		}

		if(selectedPointLightIndex != -1) // Avoid trying to access lights that dont exist
		{
		// Editing the selected point light
		ImGui::Text("Editing Point Light: %s", ("Point Light " + std::to_string(selectedPointLightIndex + 1)).c_str());
		ImGui::SliderFloat3("Position", &scene.pointLightVector[selectedPointLightIndex]->pos.x, -3.0f, 3.0f);
		ImGui::ColorEdit3("Color", &scene.pointLightVector[selectedPointLightIndex]->color.x);
		}

		ImGui::End(); // Close "Point Light Editor" window
	}

	// Display information or edit the selected spot light in a new window
	if (selectedSpotLightIndex != -1)
	{
		ImGui::Begin("Spot Light Editor");

		// Delete Spot Light button
		if (ImGui::Button("Delete Spot Light"))
		{
			// Delete the selected spot light and update the index
			delete scene.spotLightVector[selectedSpotLightIndex];
			scene.spotLightVector.erase(scene.spotLightVector.begin() + selectedSpotLightIndex);

			// Adjust selected indices
			selectedSpotLightIndex = -1; // No spot light selected after deletion
		}

		// Close button 
		ImGui::SameLine();
		if (ImGui::Button("Close"))
		{
			selectedSpotLightIndex = -1; // No spot light selected
		}

		if (selectedSpotLightIndex != -1) // Avoid trying to access lights that dont exist
		{ 
		// Editing the selected spot light
		ImGui::Text("Editing Spot Light: %s", ("Spot Light " + std::to_string(selectedSpotLightIndex + 1)).c_str());
		ImGui::SliderFloat3("Position", &scene.spotLightVector[selectedSpotLightIndex]->pos.x, -3.0f, 3.0f);
		ImGui::ColorEdit3("Color", &scene.spotLightVector[selectedSpotLightIndex]->color.x);
		ImGui::SliderFloat3("Direction", &scene.spotLightVector[selectedSpotLightIndex]->direction.x, -1.0f, 1.0f);
		ImGui::SliderFloat("Cone Cutoff", &scene.spotLightVector[selectedSpotLightIndex]->cosCutoff, 0.0f, 3.0f);
		}
		ImGui::End(); // Close "Spot Light Editor" window
	}

	ImGui::SliderFloat("Sphere Refraction", &scene.sphereArray[0]->sphereMaterial.refractionCoefficient, 1.0f, 4.0f);

	// ray query on mouse
	Ray r = camera.GetPrimaryRay((float)mousePos.x, (float)mousePos.y);
	scene.FindNearest(r);
	ImGui::Text("voxel: %i", r.voxel);

	if (ImGui::IsAnyItemActive())
	{
		// if imgui changes, empty accumulator
		memset(accumulator, 0, SCRWIDTH* SCRHEIGHT * 16);
		elapsedFrameCounter = 0;
	}
}

// -----------------------------------------------------------
// User wants to close down
// -----------------------------------------------------------
void Renderer::Shutdown()
{
	// save current camera
	//FILE* f;
	//f = fopen("camera.bin", "wb");
	//if (f)
	//{
	//	// File opened successfully, proceed with writing and closing
	//	fwrite(&camera, 1, sizeof(Camera), f);
	//	fclose(f);
	//}
	//else
	//{
	//	printf("camera.bin could not be opened\n");
	//}
}

void Tmpl8::Renderer::PointLightHandling(Ray& ray, PointLight* light, float3& I, const float3 L, float3& albedo, float3& N, float3& shadingColor)
{
	// Calculate light direction and distance
	float3 lightDirection = normalize(light->pos - I);
	float lightDistance = length(light->pos - I);

	if (ray.GetReflectivity() == 1.0f)
	{
		return;
	}


	// Calculate diffuse shading
	float diffuseIntensity = max(0.0f, dot(N, lightDirection));
	if (diffuseIntensity == 0) 
	{
		return;
	}

	// Check if the light source is behind the object
	if (dot(N, L) < 0.0f)
	{
		return; // Light source is behind the object, skip shading
	}

	// Shadow ray
	Ray shadowRay(I + 0.001f * N, lightDirection); // Offset starting point slightly to avoid self-intersection
	shadowRay.t = lightDistance;

	// If the shadow ray intersects an object, the point is in shadow
	if (scene.IsOccluded(shadowRay, light->pos)) {
		return; // Skip shading for this light, as the point is in shadow
	}

	float3 diffuseColor = albedo * light->color * diffuseIntensity;

	// Calculate specular shading (might need to modify this based on your material properties)
	if(diffuseIntensity > 0.0f){

		// Attenuation (inverse square law)
		float attenuation = 1.0f / (lightDistance * lightDistance);

		// Accumulate lighting contributions from each light
		shadingColor += (diffuseColor) * attenuation;
	}
	else 
	{
		// Only add diffuse color if there is no specular reflection
		shadingColor += diffuseColor;
	}
	return;
}

void Tmpl8::Renderer::SpotLightHandling(Ray& ray, SpotLight* light, float3& I, const float3 L, float3& albedo, float3& N, float3& shadingColor)
{
	// Calculate light direction and distance
	float lightDistance = length(light->pos - I);

	if (ray.GetReflectivity() == 1.0f)
	{
		return;
	}

	// Check if the point is within the spot cone
	float cosTheta = dot(-light->direction, light->direction);
	if (cosTheta < light->cosCutoff)
	{
		return; // Point is outside the spot cone, skip shading
	}

	// Shadow ray
	Ray shadowRay(I + 0.001f * N, light->direction); // Offset starting point slightly to avoid self-intersection
	shadowRay.t = lightDistance; 

	// If the shadow ray intersects an object, the point is in shadow
	if (scene.IsOccluded(shadowRay, light->pos)) {
		return; // Skip shading for this light, as the point is in shadow
	}

	// Calculate diffuse shading
	float diffuseIntensity = max(0.0f, dot(N, light->direction));
	float3 diffuseColor = albedo * light->color * diffuseIntensity;

	// Calculate specular shading (modify based on your material properties)
	float3 viewDirection = normalize(ray.O - I);
	float3 reflectDirection = reflect(-light->direction, N);
	float specularIntensity = pow(max(0.0f, dot(viewDirection, reflectDirection)), ray.GetReflectivity());
	float3 specularColor = light->color * specularIntensity;

	// Attenuation (inverse square law)
	float attenuation = 1.0f / (lightDistance * lightDistance);

	// Accumulate lighting contributions from each light
	shadingColor += (diffuseColor + specularColor) * attenuation;
	return;
}

void Renderer::AreaLightHandling(Ray& ray, AreaLight* light, const float3& I, const float3& albedo, const float3& N, float3& shadingColor)
{

	float3 lightDirection = normalize(light->position - I);
	float lightDistance = length(light->position - I);

	//Check if it is a perfect mirror
	if (ray.GetReflectivity() == 1.0f)
	{
		return;
	}

	// Shadow ray
	Ray shadowRay(I + 0.001f * N, lightDirection);
	shadowRay.t = lightDistance;

	// If the shadow ray intersects an object, the point is in shadow
	if (scene.IsOccluded(shadowRay, light->position))
	{
		return; // Skip shading for this light, as the point is in shadow
	}

	// Calculate diffuse shading
	float diffuseIntensity = max(0.0f, dot(N, lightDirection));
	float3 diffuseColor = albedo * light->color * diffuseIntensity;

	// Calculate specular shading
	float3 viewDirection = normalize(ray.O - I);
	float3 reflectDirection = reflect(-lightDirection, N);

	// Incorporate the area light's normal in the specular calculation
	float specularIntensity = pow(max(0.0f, dot(viewDirection, normalize(reflectDirection + light->normal))), ray.GetReflectivity());
	float3 specularColor = light->color * specularIntensity;

	// Attenuation (inverse square law)
	float attenuation = 1.0f / (lightDistance * lightDistance);

	// Accumulate lighting contributions from each area light
	shadingColor += (diffuseColor + specularColor) * attenuation * light->area;
}

Ray Tmpl8::Renderer::CreateReflectionRay(const float3& reflectionPoint, const float3& reflectionDirection, const float3& surfaceNormal, const float reflectionCoefficient, const int recursionDepth)
{
	// Calculate the reflected direction using the surface normal and incident direction
	float3 reflectedDirection = reflect(reflectionDirection, surfaceNormal);

	// Create a new ray for reflection
	Ray reflectionRay;
	reflectionRay.O = reflectionPoint + 0.001f * surfaceNormal;  // Offset starting point slightly to avoid self-intersection
	reflectionRay.D = reflectedDirection;
	reflectionRay.reflectionCoefficient = reflectionCoefficient;
	reflectionRay.recursionDepth ++;  // Increment recursion depth

	return reflectionRay;
}

