using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;
using Random = UnityEngine.Random;

public class GravitySimulator : MonoBehaviour
{
    public int numParticlesPerGalaxy;    
    public float G;              // Gravitational constant, adjust for simulation scale
    public float softeningFactor; // Softening factor to prevent extreme forces at close distances
    public int localRadius;
    public float GalaxySize;
    public float centralMass;
    public float theta;
    public float thresholdDistance;    // Distance within which repulsive force acts
    public float softeningRate;    // Strength of the repulsive force
    public float largeScalesSoftening;
    public float largeScaleG;
    public Vector2 minMaxMassValues;
    public Vector3 simulationBounds;
    public Vector3 Galaxy1Position;
    public Vector3 Galaxy2Position;
    public Vector3 Galaxy1Velocity;
    public Vector3 Galaxy2Velocity;
    public GravityObject gravityObject;
    private PointOctree<GravityObject> gravityObjects;
    private List<PointOctreeNode<GravityObject>> nodes;
    private List<GravityObject> gravList;
    private int frameCount;
    public void CreateGalaxyGauss(int numParticles, Vector3 galaxyPosition, Vector3 initialVelocity) {
    for (int i = 0; i < numParticles; i++) {
        // Generate a Gaussian distance from the center using the Box-Muller transform for central clustering
        float gaussianDistance = Mathf.Abs(BoxMullerGaussian() * GalaxySize / 3); // Scale to GalaxySize with clustering near center
        float angle = Random.Range(0, 2 * Mathf.PI);

        // Position in a circular disk around the galaxy center with slight vertical variation
        Vector3 position = galaxyPosition + new Vector3(
            gaussianDistance * Mathf.Cos(angle),
            Random.Range(-GalaxySize / 40, GalaxySize / 40), // Reduced thickness for stability
            gaussianDistance * Mathf.Sin(angle)
        );

        // Instantiate particle
        GravityObject particle = Instantiate(gravityObject, position, Quaternion.identity);
        particle.mass = 1;

        // Calculate initial tangential velocity, reduced based on distance
        Vector3 toCenter = galaxyPosition - particle.transform.position;
        Vector3 tangentialVelocity = Vector3.Cross(toCenter.normalized, Vector3.up) * Mathf.Sqrt(G * centralMass / toCenter.magnitude);
        tangentialVelocity *= 0.02f; // Scale down tangential velocity to reduce outward scattering

        // Apply tangential velocity and initial velocity of the galaxy
        particle.velocity = tangentialVelocity + initialVelocity;
        gravList.Add(particle);
    }
}
    private float BoxMullerGaussian() {
    float u1 = Random.Range(0.0001f, 1f); // Uniform(0,1] random values
    float u2 = Random.Range(0.0001f, 1f);
    float randStdNormal = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Sin(2.0f * Mathf.PI * u2); // Standard normal
    return randStdNormal;
    }
    public bool outOfBounds(GravityObject obj) {
    
    Vector3 pos = obj.transform.position;
    
    return pos.x < -simulationBounds.x || pos.x > simulationBounds.x ||
           pos.y < -simulationBounds.y || pos.y > simulationBounds.y ||
           pos.z < -simulationBounds.z || pos.z > simulationBounds.z;
    }
    public void CreateEqualDistribution(int numParticle){
        for (int i = 0; i < numParticle; i++)
        {
            GravityObject particle = Instantiate(gravityObject, new Vector3(Random.Range(-simulationBounds.x/2,simulationBounds.x/2),
                                                                            Random.Range(-simulationBounds.y/2,simulationBounds.y/2),
                                                                            Random.Range(-simulationBounds.z/2,simulationBounds.z/2)),
                                                                            Quaternion.identity);
            particle.mass = Random.Range(minMaxMassValues.x, minMaxMassValues.y);
            gravList.Add(particle);

        }
    }
    void Start()
    {
        gravList = new List<GravityObject>();
        gravityObjects = new PointOctree<GravityObject>(Math.Max(Math.Max(simulationBounds.x,simulationBounds.y),simulationBounds.z), transform.position, 1);
        CreateGalaxyGauss(numParticlesPerGalaxy,Galaxy1Position,Galaxy1Velocity);
        CreateGalaxyGauss(numParticlesPerGalaxy,Galaxy2Position,Galaxy2Velocity);
        // CreateEqualDistribution(numParticlesPerGalaxy);
        nodes = gravityObjects.GetAllNodes();
    }

    void Update(){
        List<GravityObject> toRemove = new List<GravityObject>();
        frameCount++;
        float timeStep = Time.fixedDeltaTime;
        int numChecks = 0;


        gravityObjects.rebuild(gravList);
        gravityObjects.getRoot().CalculateCenterOfMass();
        List<PointOctreeNode<GravityObject>> nodes = gravityObjects.GetAllNodes();
        NativeArray<NodeData> nodeData  = new NativeArray<NodeData>(nodes.Count,Allocator.TempJob);

        for (int i = 0; i < nodes.Count; i++)
        {
            var node = nodes[i];
            nodeData[i] = new NodeData{
                CenterOfMass=node.CenterOfMass,
                TotalMass=node.TotalMass,
                position=node.Center,
                Size=node.SideLength,
                isLeaf=!node.HasChildren,
            };
        }

        NativeArray<GravityObjectData> particleArray  = new NativeArray<GravityObjectData>(gravList.Count,Allocator.TempJob);

        int totalNearbyCount = 0;
        foreach (var particle in gravList) {
            totalNearbyCount += gravityObjects.GetNearby(particle.Position, localRadius).Length;
        }
        NativeArray<GravityObjectData> nearbyParticlesFlatArray = new NativeArray<GravityObjectData>(totalNearbyCount, Allocator.TempJob);
        NativeArray<int> neighborhoodIndex = new NativeArray<int>(gravList.Count + 1, Allocator.TempJob); // +1 for easy range handling

        int flatIndex = 0;
        for(int i = 0; i < gravList.Count; i++){
            var particle = gravList[i];
            particleArray[i] = new GravityObjectData{
                position = particle.getPosition(),
                velocity = particle.velocity,
                mass = particle.mass
            };
        
            neighborhoodIndex[i] = flatIndex;
            var neighbors = gravityObjects.GetNearby(particle.Position, localRadius);
            foreach (var neighbor in neighbors)
            {
                if(neighbor!=particle){
                    nearbyParticlesFlatArray[flatIndex++] = new GravityObjectData{
                        position = neighbor.Position,
                        velocity = neighbor.velocity,
                        mass = neighbor.mass
                    };
                }
            }
            numChecks+=neighbors.Length;
        }
        neighborhoodIndex[gravList.Count] = flatIndex;
        NativeArray<float3> forces = new NativeArray<float3>(gravList.Count, Allocator.TempJob);
        LocalForcesJob localForcesJob = new LocalForcesJob {
            particles = particleArray,
            nearbyParticlesFlatArray = nearbyParticlesFlatArray,
            neighborhoodIndex = neighborhoodIndex,
            forces = forces,
            G = G,
            baseSofteningFactor = softeningFactor,
            softeningRate = softeningRate,
            thresholdDistance = thresholdDistance
        };
        JobHandle localForcesJobHandle = localForcesJob.Schedule(particleArray.Length, 64);

        LargeScaleForceJob largeScaleForceJob = new LargeScaleForceJob{
            particles = particleArray,
            nodes = nodeData,
            forces = forces,
            G = largeScaleG,
            softeningFactor = largeScalesSoftening,
            theta = theta

        };
        JobHandle largeScaleForceJobHandle = largeScaleForceJob.Schedule(particleArray.Length,64,localForcesJobHandle);
        UpdatePositionJob updatePositionJob = new UpdatePositionJob{
            particles = particleArray,
            forces = forces,
            deltaTime = timeStep
        };
        JobHandle updatePositionJobHandle = updatePositionJob.Schedule(particleArray.Length, 64, largeScaleForceJobHandle);
        updatePositionJobHandle.Complete();
        for (int i = 0; i < gravList.Count; i++) {
            gravList[i].UpdatePosition(particleArray[i].position);
            if (outOfBounds(gravList[i])) {
                toRemove.Add(gravList[i]);
            }
        }
        neighborhoodIndex.Dispose();
        nearbyParticlesFlatArray.Dispose();
        particleArray.Dispose();
        forces.Dispose();
        nodeData.Dispose();

        gravList.RemoveAll(obj => toRemove.Contains(obj));
        // foreach (var obj in toRemove) {
        //         Destroy(obj.gameObject);
        // }
        toRemove.Clear();
        Debug.Log(numChecks+ " " + gravityObjects.getRoot().CenterOfMass);
    }

    void OnDrawGizmos() 
    {
        if(gravityObjects!=null){
        gravityObjects.DrawAllBounds();
	    gravityObjects.DrawAllObjects();
        }
    }
}
[BurstCompile]
public struct LocalForcesJob : IJobParallelFor {
    [Unity.Collections.ReadOnly] public NativeArray<GravityObjectData> particles;
    [Unity.Collections.ReadOnly] public NativeArray<GravityObjectData> nearbyParticlesFlatArray;
    [Unity.Collections.ReadOnly] public NativeArray<int> neighborhoodIndex;
    public NativeArray<float3> forces;
    public float G;
    public float baseSofteningFactor;
    public float softeningRate;     
    public float thresholdDistance;

    public void Execute(int index) {
        GravityObjectData particle = particles[index];
        float3 totalForce = float3.zero;

        int startIndex = neighborhoodIndex[index];
        int endIndex = neighborhoodIndex[index + 1];

        for (int i = startIndex; i < endIndex; i++) {
            GravityObjectData other = nearbyParticlesFlatArray[i];
            totalForce += ForceCalculator.CalculateSoftenedForce(
                particle,
                other,
                G,
                baseSofteningFactor,
                softeningRate,
                thresholdDistance
            );
        }

        forces[index] = totalForce;
    }
} 
[BurstCompile]
public struct LargeScaleForceJob : IJobParallelFor{
    [Unity.Collections.ReadOnly] public NativeArray<GravityObjectData> particles;
    [Unity.Collections.ReadOnly] public NativeArray<NodeData> nodes;
    public NativeArray<float3> forces;
    public float theta;
    public float G;
    public float softeningFactor;
    public void Execute(int index){
        GravityObjectData particle = particles[index];
        float3 totalForce = ForceCalculator.CalculateLargeScaleForce(particle,nodes,G,theta,softeningFactor);
        forces[index]+=totalForce;
    }
    
}


[BurstCompile]
public struct UpdatePositionJob : IJobParallelFor{
    public NativeArray<GravityObjectData> particles;
    [Unity.Collections.ReadOnly] public NativeArray<float3> forces;
    public float deltaTime;

    public void Execute(int index) {
        GravityObjectData particle = particles[index];

        // Step 1: Half-step velocity update
        float3 acceleration = forces[index] / particle.mass;
        particle.velocity += 0.5f * acceleration * deltaTime;

        // Step 2: Full position update
        particle.position += particle.velocity * deltaTime;

        // Step 3: Half-step velocity update with new acceleration
        particle.velocity += 0.5f * acceleration * deltaTime;

        // Store the updated particle data back
        particles[index] = particle;
    }
}
public static class ForceCalculator{
    [BurstCompile]
    public static float3 CalculateSoftenedForce(GravityObjectData particle, GravityObjectData other, float G, float baseSofteningFactor, float softeningRate, float thresholdDistance) {
    float3 direction = other.position - particle.position;
    float distance = math.length(direction);
    
    // Calculate softening factor based on distance
    float softeningFactor;
    if (distance < thresholdDistance) {
        softeningFactor = baseSofteningFactor * math.exp(softeningRate * (thresholdDistance - distance));
    } else {
        softeningFactor = baseSofteningFactor;
    }

    // Calculate gravitational force with distance-dependent softening
    float distanceSquared = distance * distance + softeningFactor * softeningFactor;
    float forceMagnitude = G * particle.mass * other.mass / distanceSquared;
    return math.normalize(direction) * forceMagnitude;
    }
    [BurstCompile]
    public static float3 CalculateLargeScaleForce(GravityObjectData particle, NativeArray<NodeData> nodes, float G, float theta, float softeningFactor){
    float3 totalForce = float3.zero;

    // Iterate over all nodes to compute forces
    for (int i = 0; i < nodes.Length; i++)
    {
        NodeData node = nodes[i];

        // Skip empty nodes
        if (node.TotalMass == 0) continue;

        // Check if the node contains the current particle's position
        bool isSameNode = math.any(math.abs(node.position - particle.position) < node.Size * 0.5f);

            // For inter-node interactions, use the Barnes-Hut approximation
        float3 direction = node.CenterOfMass - particle.position;
        float distance = math.length(direction);
        float s = node.Size;

            // Apply theta criterion
        if (s / distance < theta || node.isLeaf){
            float softenedDistance = distance + softeningFactor;
            float forceMagnitude = (G * particle.mass * node.TotalMass) / (softenedDistance * softenedDistance * softenedDistance);
            totalForce += direction * forceMagnitude;
        }
    }
    return totalForce;
    }
}
