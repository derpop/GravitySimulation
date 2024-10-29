using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;
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
        particle.mass = Random.Range(minMaxMassValues.x, minMaxMassValues.y);

        // Calculate initial tangential velocity, reduced based on distance
        Vector3 toCenter = galaxyPosition - particle.transform.position;
        Vector3 tangentialVelocity = Vector3.Cross(toCenter.normalized, Vector3.up) * Mathf.Sqrt(G * centralMass / toCenter.magnitude);
        tangentialVelocity *= 0.5f; // Scale down tangential velocity to reduce outward scattering

        // Apply tangential velocity and initial velocity of the galaxy
        particle.velocity = tangentialVelocity + initialVelocity;
        gravList.Add(particle);
    }
}

// Box-Muller transform to generate Gaussian-distributed values
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
        // CreateGalaxy(numParticlesPerGalaxy,Galaxy2Position,Galaxy2Velocity);
        // CreateEqualDistribution(numParticlesPerGalaxy/10);
        nodes = gravityObjects.GetAllNodes();
    }

    // Update is called once per frame
//     void Update()
// {
//     List<GravityObject> toRemove = new List<GravityObject>();
//     frameCount++;
//     List<PointOctreeNode<GravityObject>> nodes = gravityObjects.GetAllNodes();
//     float timeStep = Time.fixedDeltaTime;
//     int numChecks =0;
//     if(frameCount%refresh==0){
//         gravityObjects.rebuild(gravList);
//         nodes = gravityObjects.GetAllNodes();
//         gravityObjects.getRoot().CalculateCenterOfMass();
//         foreach (GravityObject obj in gravList)
//         {
//         GravityObject[] nearbyObjects = gravityObjects.GetNearby(obj.transform.position, localRadius);

//         numChecks += nearbyObjects.Length;
//         obj.UpdatePosition(timeStep);
//         if(outOfBounds(obj)){
//             toRemove.Add(obj);
//         }
//         }
//     }
//     gravList.RemoveAll(obj => toRemove.Contains(obj));
//     foreach (var obj in toRemove) {
//         UnityEngine.Object.Destroy(obj.gameObject); // Assuming each GravityObject has an associated GameObject
//     }
    
//     // Clear toRemove list after destruction
//     toRemove.Clear();
//     Debug.Log(numChecks);

// }
    void Update(){
        List<GravityObject> toRemove = new List<GravityObject>();
        frameCount++;
        float timeStep = Time.fixedDeltaTime;
        int numChecks = 0;

        List<float3> centersOfMass = new List<float3>();
        List<float> totalMasses = new List<float>();
        List<float> sideLengths = new List<float>();

        gravityObjects.rebuild(gravList);
        gravityObjects.getRoot().CalculateCenterOfMass();
        gravityObjects.getAllNodes(gravityObjects.getRoot(),centersOfMass,totalMasses,sideLengths);

        NativeArray<float3> nodeCentersOfMass = new NativeArray<float3>(centersOfMass.ToArray(), Allocator.TempJob);
        NativeArray<float> nodeTotalMasses = new NativeArray<float>(totalMasses.ToArray(), Allocator.TempJob);
        NativeArray<float> nodeSideLengths = new NativeArray<float>(sideLengths.ToArray(), Allocator.TempJob);

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
            softeningFactor = softeningFactor
        };
        JobHandle localForcesJobHandle = localForcesJob.Schedule(particleArray.Length, 64);
        LargeScaleForcesJob largeScaleForcesJob = new LargeScaleForcesJob{
            particles = particleArray,
            nodeCentersOfMass = nodeCentersOfMass,
            nodeTotalMasses = nodeTotalMasses,
            nodeSideLengths = nodeSideLengths,
            forces = forces,
            G = G,
            softeningFactor = softeningFactor,
            theta = theta   
        };
        JobHandle largeScaleForcesJobHandle = largeScaleForcesJob.Schedule(particleArray.Length, 64, localForcesJobHandle);
        UpdatePositionJob updatePositionJob = new UpdatePositionJob{
            particles = particleArray,
            forces = forces,
            deltaTime = timeStep
        };
        JobHandle updatePositionJobHandle = updatePositionJob.Schedule(particleArray.Length, 64, largeScaleForcesJobHandle);
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
        nodeCentersOfMass.Dispose();
        nodeTotalMasses.Dispose();
        nodeSideLengths.Dispose();

        gravList.RemoveAll(obj => toRemove.Contains(obj));
        // foreach (var obj in toRemove) {
        //         Destroy(obj.gameObject);
        // }
        toRemove.Clear();
        Debug.Log(numChecks);
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
    public float softeningFactor;

    public void Execute(int index){
        GravityObjectData particle = particles[index];
        float3 totalForce = float3.zero;
        
        int startIndex = neighborhoodIndex[index];
        int endIndex = neighborhoodIndex[index + 1];

        for (int i = startIndex; i < endIndex; i++) {
            GravityObjectData other = nearbyParticlesFlatArray[i];
            totalForce += GravityCalculations.CalculateLocalForce(particle, other, G, softeningFactor);
        }

        forces[index] = totalForce;
    }

}
[BurstCompile]
public struct LargeScaleForcesJob : IJobParallelFor{
    [Unity.Collections.ReadOnly] public NativeArray<GravityObjectData> particles;
    [Unity.Collections.ReadOnly] public NativeArray<float3> nodeCentersOfMass; 
    [Unity.Collections.ReadOnly]  public NativeArray<float> nodeTotalMasses;
    [Unity.Collections.ReadOnly]  public NativeArray<float> nodeSideLengths;
    public NativeArray<float3> forces;
    public float G; 
    public float theta;
    public float softeningFactor;
    public void Execute(int index) {
        GravityObjectData particle = particles[index];
        float3 totalForce = GravityCalculations.CalculateLargeScaleForce(
            particle,
            nodeCentersOfMass,
            nodeTotalMasses,
            nodeSideLengths,
            G,
            softeningFactor,
            theta
        );

        // Accumulate the large-scale force into the forces array
        forces[index] += totalForce;
    }
}
[BurstCompile]
public struct UpdatePositionJob : IJobParallelFor{
    public NativeArray<GravityObjectData> particles; 
    [Unity.Collections.ReadOnly] public NativeArray<float3> forces; 
    public float deltaTime;

    public void Execute(int index){

        GravityObjectData particle = particles[index];

        float3 acceleration = forces[index] / particle.mass;
        particle.velocity += acceleration * deltaTime;

        particle.position += particle.velocity * deltaTime;

        particles[index] = particle;
    }
}
