using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;
using Unity.Collections;
public interface IObject {
    Vector3 Position { get; set; }
    float mass {get; set;}
}
public struct GravityObjectData{
    public float3 position;
    public float3 velocity;
    public float mass;
}
public static class GravityCalculations {
    [BurstCompile]
    public static float3 CalculateLocalForce(GravityObjectData particle, GravityObjectData other, float G, float softeningFactor) {
        if (particle.position.Equals(other.position)) return float3.zero;

        float3 direction = other.position - particle.position;
        float distance = math.sqrt(math.lengthsq(direction));
        float forceMagnitude = G * particle.mass * other.mass / (distance * distance + softeningFactor*softeningFactor);

        return math.normalize(direction) * forceMagnitude;
    }

    [BurstCompile]
    public static float3 CalculateLargeScaleForce(GravityObjectData particle, NativeArray<float3> nodeCentersOfMass, NativeArray<float> nodeTotalMasses, NativeArray<float> nodeSideLengths, float G, float softeningFactor, float theta){
        float3 totalForce = float3.zero;
        for(int i = 0; i < nodeCentersOfMass.Length; i++){
            if(nodeTotalMasses[i]==0) continue;

            float3 direction = nodeCentersOfMass[i] - particle.position;
            float distance = math.length(direction);

            if(distance==0) continue;

            if ((nodeSideLengths[i] / distance) < theta) {
                float forceMagnitude = G * particle.mass * nodeTotalMasses[i] / (distance * distance + softeningFactor*softeningFactor);
                totalForce += math.normalize(direction) * forceMagnitude;
            }
        }
        return totalForce;
    }
}
public class GravityObject : MonoBehaviour, IObject
{
    public float mass {get; set;}       
    public Vector3 velocity = Vector3.zero;   

    public Vector3 Position { get; set; }

    /// <summary>
    /// Calculates and applies gravitational forces from other particles to update this particle's velocity.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <param name="timeStep">Time step for the Euler integration</param>
    public void ApplyForces(Vector3 force,float timeStep)
    {
        // Calculate total gravitational acceleration acting on this particle
        Vector3 totalAcceleration = force/mass;
        velocity += totalAcceleration * timeStep;
    }

    /// <summary>
    /// Updates the position of the particle based on its current velocity, using Euler's method.
    /// </summary>
    /// <param name="timeStep">Time step for position update</param>
    public void UpdatePosition(float3 newPos)
    {
        transform.position = newPos;
        Position = transform.position;
    }
    public Vector3 getPosition(){
        Position = transform.position;
        return Position;
    }
    public Vector3 getVelocity(){
        return velocity;
    }
}
