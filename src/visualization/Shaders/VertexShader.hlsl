// Basic vertex shader for rendering 3D meshes
// Input structure from the vertex buffer
struct VSInput
{
    float3 Position : POSITION;
};

// Output structure sent to the pixel shader
struct VSOutput
{
    float4 Position : SV_POSITION;
    float3 Color : COLOR;
};

// Constant buffer for transformation matrices
cbuffer TransformCB : register(b0)
{
    matrix World;
    matrix View;
    matrix Projection;
};

// Entry point of the vertex shader
VSOutput main(VSInput input)
{
    VSOutput output;
    
    // Transform the vertex position through world, view, and projection matrices
    input.Position += float3(0, 0, -50); // TODO: remove this
    float4 worldPosition = mul(float4(input.Position, 1.0f), World);
    float4 viewPosition = mul(worldPosition, View);
    output.Position = mul(viewPosition, Projection);
    
    // Generate a simple color based on position
    output.Color = normalize(input.Position) * 0.5f + 0.5f;
    
    return output;
} 